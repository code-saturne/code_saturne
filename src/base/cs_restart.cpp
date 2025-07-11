/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_FCNTL_H)
# include <fcntl.h>
#endif

#if defined(HAVE_UNISTD_H)
# include <unistd.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_io_num.h"

#include "base/cs_all_to_all.h"
#include "base/cs_base.h"
#include "base/cs_block_dist.h"
#include "base/cs_block_to_part.h"
#include "base/cs_file.h"
#include "base/cs_io.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_save.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_part_to_block.h"
#include "base/cs_parall.h"
#include "base/cs_timer.h"
#include "base/cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_restart.cpp
        Manage checkpoint / restart files.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

typedef struct _location_t {

  char             *name;             /* Location name */
  size_t            id;               /* Associated id in file */
  cs_lnum_t         n_ents;           /* Local number of entities */
  cs_gnum_t         n_glob_ents_f;    /* Global number of entities by file */
  cs_gnum_t         n_glob_ents;      /* Global number of entities */
  const cs_gnum_t  *ent_global_num;   /* Possibly shared global entity numbers,
                                         or nullptr */
  cs_gnum_t        *_ent_global_num;  /* Private global entity numbers,
                                         or nullptr */

} _location_t;

struct _cs_restart_t {

  char              *name;           /* Name of restart file */

  cs_io_t           *fh;             /* Pointer to associated file handle */
  int                rank_step;      /* Block rank step for parallel IO */
  int                min_block_size; /* Minimum block size for parallel IO */

  size_t             n_locations;    /* Number of locations */
  _location_t       *location;       /* Location definition array */

  cs_restart_mode_t  mode;           /* Read or write */

};

typedef struct {

  int    id;                 /* Id of the writer */

  char  *name;              /* Name of the checkpoint file */
  char  *path;              /* Full path to the checkpoint file */

  int    n_prev_files;      /* Number of times this file has already
                               been written */
  int    n_prev_files_tot;  /* Total number of times this file has already
                               been written */
  char **prev_files;        /* Names of the previous versions */

} _restart_multiwriter_t;

/*============================================================================
 * Prototypes for private functions
 *============================================================================*/

static int
_check_section(cs_restart_t           *restart,
               void                   *context,
               const char             *sec_name,
               int                     location_id,
               int                     n_location_vals,
               cs_restart_val_type_t   val_type);

static int
_read_section(cs_restart_t           *restart,
              void                   *context,
              const char             *sec_name,
              int                     location_id,
              int                     n_location_vals,
              cs_restart_val_type_t   val_type,
              void                   *val);

static void
_write_section(cs_restart_t           *restart,
               void                   *context,
               const char             *sec_name,
               int                     location_id,
               int                     n_location_vals,
               cs_restart_val_type_t   val_type,
               const void             *val);

/*============================================================================
 * Static global variables
 *============================================================================*/

#if defined(WIN32) || defined(_WIN32)
static const char _dir_separator = '\\';
#else
static const char _dir_separator = '/';
#endif

/* Monitoring info */

static int    _restart_n_opens[2] = {0, 0};
static double _restart_wtime[2] = {0.0, 0.0};

/* Do we have a restart directory ? */

static int _restart_present = -1;

/* Restart time steps and frequency */

static int    _checkpoint_mesh = 1;          /* checkpoint mesh if possible */
/* time step interval */
static int    _checkpoint_nt_interval = CS_RESTART_INTERVAL_DEFAULT;
static int    _checkpoint_nt_next = -1;      /* next forced time step */
static int    _checkpoint_nt_last = -1;      /* last checkpoint time step */
static double _checkpoint_t_interval = -1.;  /* physical time interval */
static double _checkpoint_t_next = -1.;      /* next forced time value */
static double _checkpoint_t_last = 0.;       /* last forced time value */
static double _checkpoint_wt_interval = -1.; /* wall-clock interval */
static double _checkpoint_wt_next = -1.;     /* next forced wall-clock value */
static double _checkpoint_wt_last = 0.;      /* wall-clock time of last
                                                checkpointing */
/* Are we restarting from a NCFD file ? */

static int    _restart_from_ncfd = 0;

/* Restart modification */

static void                        *_restart_context = nullptr;
static cs_restart_check_section_t  *_check_section_f = _check_section;
static cs_restart_read_section_t   *_read_section_f  = _read_section;
static cs_restart_write_section_t  *_write_section_f = _write_section;

static size_t        _n_locations_ref;    /* Number of locations */
static _location_t  *_location_ref;       /* Location definition array */

/* Restart multi writer */

static int                       _n_restart_directories_to_write = 1;
static int                       _n_restart_multiwriters          = 0;
static _restart_multiwriter_t  **_restart_multiwriter             = nullptr;

/* Special case for in-memory buffer, which supercedes file I/O, and
   concatenates all pseudo-file checkpoint/restart data to a single file */

static bool      _need_finalize = false;
static cs_io_t  *_checkpoint_serialized_memory = nullptr;
static cs_io_t  *_restart_serialized_memory = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Call cleanup operations for checkpoint/restart subsystem.
 *----------------------------------------------------------------------------*/

static void
_restart_finalize(void)
{
  if (_checkpoint_serialized_memory != nullptr)
    cs_io_finalize(&_checkpoint_serialized_memory);
  if (_restart_serialized_memory != nullptr)
    cs_io_finalize(&_restart_serialized_memory);
}

/*----------------------------------------------------------------------------
 * Compute number of values in a record
 *
 * parameters:
 *   r               <-- associated restart file pointer
 *   location_id     <-- location id
 *   n_location_vals <-- number of values per location
 *----------------------------------------------------------------------------*/

static cs_gnum_t
_compute_n_ents(const cs_restart_t  *r,
                size_t               location_id,
                size_t               n_location_vals)
{
  cs_gnum_t retval = 0;

  if (location_id == 0)
    retval = n_location_vals;

  else if (location_id > 0 && location_id <= r->n_locations)
    retval =   r->location[location_id-1].n_glob_ents_f
             * (cs_gnum_t)n_location_vals;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              (int)location_id, r->name);

  return retval;
}

/*----------------------------------------------------------------------------
 * Analyze the content of a restart file to build locations
 *
 * parameters:
 *   r <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

static void
_locations_from_index(cs_restart_t  *r)
{
  cs_io_sec_header_t h;

  size_t rec_id = 0;
  size_t index_size = 0;

  /* Initialization */

  index_size = cs_io_get_index_size(r->fh);

  /* Analyze records to determine locations */

  for (rec_id = 0; rec_id < index_size; rec_id++) {

    h = cs_io_get_indexed_sec_header(r->fh, rec_id);

    if (h.location_id > r->n_locations) {

      _location_t  *loc = nullptr;

      if (h.location_id != r->n_locations + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Restart file \"%s\" declares a location number %d\n"
                    "but no location %d has been declared."),
                  r->name, (int)(h.location_id),
                  (int)(r->n_locations + 1));

      CS_REALLOC(r->location, r->n_locations + 1, _location_t);

      loc = r->location + r->n_locations;
      CS_MALLOC(loc->name, strlen(h.sec_name) + 1, char);
      strcpy(loc->name, h.sec_name);

      loc->id = h.location_id;
      loc->n_ents = 0;
      loc->n_glob_ents = 0;

      cs_io_set_indexed_position(r->fh, &h, rec_id);
      cs_io_set_cs_gnum(&h, r->fh);
      cs_io_read_global(&h, &(loc->n_glob_ents_f), r->fh);

      loc->ent_global_num = nullptr;
      loc->_ent_global_num = nullptr;

      r->n_locations += 1;
    }

  }
}

/*----------------------------------------------------------------------------
 * Initialize a checkpoint / restart file management structure;
 *
 * parameters:
 *   r <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

static void
_add_file(cs_restart_t  *r)
{
  double timing[2];
  cs_file_access_t method;

  const char magic_string[] = "Checkpoint / restart, R0";
  const long echo = CS_IO_ECHO_NONE;

  timing[0] = cs_timer_wtime();

  /* In read mode, open file to detect header first */

#if defined(HAVE_MPI)
  {
    MPI_Info           hints;
    MPI_Comm           block_comm, comm;

    cs_file_get_default_comm(nullptr, &block_comm, &comm);

    r->rank_step = 1;
    r->min_block_size = cs_parall_get_min_coll_buf_size();
    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);

    if (r->mode == CS_RESTART_MODE_READ) {
      cs_file_get_default_access(CS_FILE_MODE_READ, &method, &hints);
      r->fh = cs_io_initialize_with_index(r->name,
                                          magic_string,
                                          method,
                                          echo,
                                          hints,
                                          block_comm,
                                          comm);
      _locations_from_index(r);
    }
    else {
      cs_file_get_default_access(CS_FILE_MODE_WRITE, &method, &hints);
      r->fh = cs_io_initialize(r->name,
                               magic_string,
                               CS_IO_MODE_WRITE,
                               method,
                               echo,
                               hints,
                               block_comm,
                               comm);
    }
  }
#else
  {
    if (r->mode == CS_RESTART_MODE_READ) {
      cs_file_get_default_access(CS_FILE_MODE_READ, &method);
      r->fh = cs_io_initialize_with_index(r->name,
                                          magic_string,
                                          method,
                                          echo);
      _locations_from_index(r);
    }
    else {
      cs_file_get_default_access(CS_FILE_MODE_WRITE, &method);
      r->fh = cs_io_initialize(r->name,
                               magic_string,
                               CS_IO_MODE_WRITE,
                               method,
                               echo);
    }
  }
#endif

  timing[1] = cs_timer_wtime();
  _restart_wtime[r->mode] += timing[1] - timing[0];

  _restart_n_opens[r->mode] += 1;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Read variable values defined on a mesh location.
 *
 * parameters:
 *   r           <-> associated restart file pointer
 *   header          <-- header associated with current position in file
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers (1 to n numbering)
 *   n_location_vals <-- number of values par location
 *   val_type        <-- data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_read_ent_values(cs_restart_t           *r,
                 cs_io_sec_header_t     *header,
                 cs_gnum_t               n_glob_ents,
                 cs_lnum_t               n_ents,
                 const cs_gnum_t         ent_global_num[],
                 int                     n_location_vals,
                 cs_restart_val_type_t   val_type,
                 cs_byte_t               vals[])
{
  cs_byte_t  *buffer = nullptr;

  cs_lnum_t  block_buf_size = 0;

  size_t  nbr_byte_ent;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_char:
    nbr_byte_ent = n_location_vals;
    break;
  case CS_TYPE_int:
    nbr_byte_ent = n_location_vals * sizeof(int);
    cs_io_set_cs_lnum(header, r->fh);
    break;
  case CS_TYPE_cs_gnum_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_gnum_t);
    cs_io_set_cs_gnum(header, r->fh);
    break;
  case CS_TYPE_cs_real_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_real_t);
    break;
  default:
    nbr_byte_ent = 0;
    assert(0);
  }

  cs_block_dist_info_t bi
    = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                  cs_glob_n_ranks,
                                  r->rank_step,
                                  r->min_block_size / nbr_byte_ent,
                                  n_glob_ents);

  cs_all_to_all_t *d
    = cs_all_to_all_create_from_block(n_ents,
                                      CS_ALL_TO_ALL_USE_DEST_ID,
                                      ent_global_num,
                                      bi,
                                      cs_glob_mpi_comm);

  /* Read blocks */

  block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]) * nbr_byte_ent;

  if (block_buf_size > 0)
    CS_MALLOC(buffer, block_buf_size, cs_byte_t);

  cs_io_read_block(header,
                   bi.gnum_range[0],
                   bi.gnum_range[1],
                   buffer,
                   r->fh);

 /* Distribute blocks on ranks */

  cs_all_to_all_copy_array(d,
                           header->elt_type,
                           n_location_vals,
                           true,  /* reverse */
                           buffer,
                           vals);

  /* Free buffer */

  CS_FREE(buffer);

  cs_all_to_all_destroy(&d);
}

/*----------------------------------------------------------------------------
 * Write variable values defined on a mesh location.
 *
 * parameters:
 *   r               <-> associated restart file pointer
 *   sec_name        <-- section name
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers (1 to n numbering)
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values par location
 *   val_type        <-- data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_write_ent_values(const cs_restart_t     *r,
                  const char             *sec_name,
                  cs_gnum_t               n_glob_ents,
                  cs_lnum_t               n_ents,
                  const cs_gnum_t        *ent_global_num,
                  int                     location_id,
                  int                     n_location_vals,
                  cs_restart_val_type_t   val_type,
                  const cs_byte_t        *vals)
{
  cs_lnum_t  block_buf_size = 0;

  cs_datatype_t elt_type     = CS_DATATYPE_NULL;
  size_t      nbr_byte_ent = 0;
  cs_byte_t  *buffer = nullptr;

  cs_block_dist_info_t bi;

  cs_part_to_block_t *d = nullptr;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_char:
    nbr_byte_ent = n_location_vals;
    elt_type = CS_CHAR;
    break;
  case CS_TYPE_int:
    nbr_byte_ent = n_location_vals * sizeof(int);
    elt_type = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;
    break;
  case CS_TYPE_cs_gnum_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_gnum_t);
    elt_type = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
    break;
  case CS_TYPE_cs_real_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_real_t);
    elt_type =   (sizeof(cs_real_t) == cs_datatype_size[CS_DOUBLE])
               ? CS_DOUBLE : CS_FLOAT;
    break;
  default:
    assert(0);
  }

  bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                   cs_glob_n_ranks,
                                   r->rank_step,
                                   r->min_block_size / nbr_byte_ent,
                                   n_glob_ents);

  d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                      bi,
                                      n_ents,
                                      ent_global_num);

  /* Distribute to blocks */

  block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]) * nbr_byte_ent;

  if (block_buf_size > 0)
    CS_MALLOC(buffer, block_buf_size, cs_byte_t);

  /* Distribute blocks on ranks */

  cs_part_to_block_copy_array(d,
                              elt_type,
                              n_location_vals,
                              vals,
                              buffer);

  /* Write blocks */

  cs_io_write_block_buffer(sec_name,
                           n_glob_ents,
                           bi.gnum_range[0],
                           bi.gnum_range[1],
                           location_id,
                           0,
                           n_location_vals,
                           elt_type,
                           buffer,
                           r->fh);

  /* Free buffer */

  CS_FREE(buffer);

  cs_part_to_block_destroy(&d);
}

#endif /* #if defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Swap values of a renumbered array when reading
 *
 * parameters:
 *   n_ents          --> number of local entities
 *   ini_ent_num     --> initial entity numbers
 *   n_location_vals --> number of values per entity
 *   val_type        --> data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_restart_permute_read(cs_lnum_t               n_ents,
                      const cs_gnum_t        *ini_ent_num,
                      int                     n_location_vals,
                      cs_restart_val_type_t   val_type,
                      cs_byte_t              *vals)
{
  cs_lnum_t ent_id, jj;

  cs_lnum_t ii = 0;

  /* Instructions */

  if (ini_ent_num == nullptr)
    return;

  switch (val_type) {

  case CS_TYPE_char:
    {
      char  *val_ord;
      char  *val_cur = (char *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, char);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      CS_FREE(val_ord);
    }
    break;

  case CS_TYPE_int:
    {
      int  *val_ord;
      int  *val_cur = (int *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, int);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      CS_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_gnum_t:
    {
      cs_gnum_t  *val_ord;
      cs_gnum_t  *val_cur = (cs_gnum_t *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, cs_gnum_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      CS_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_real_t:
    {
      cs_real_t  *val_ord;
      cs_real_t  *val_cur = (cs_real_t *)vals;

      CS_MALLOC (val_ord, n_ents * n_location_vals, cs_real_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      CS_FREE(val_ord);
    }
    break;

  default:
    assert(0);

  }
}

/*----------------------------------------------------------------------------
 * Swap values of a renumbered array when writing
 *
 * parameters:
 *   n_ents          --> number of local entities
 *   ini_ent_num     --> initial entity numbers
 *   n_location_vals --> number of values per entity
 *   val_type        --> data type
 *   vals            --> array of values
 *
 * returns:
 *   pointer to array of values in initial entity order
 *----------------------------------------------------------------------------*/

static cs_byte_t *
_restart_permute_write(cs_lnum_t               n_ents,
                       const cs_gnum_t        *ini_ent_num,
                       int                     n_location_vals,
                       cs_restart_val_type_t   val_type,
                       const cs_byte_t        *vals)
{
  cs_lnum_t  ent_id, jj;

  cs_lnum_t  ii = 0;

  /* Instructions */

  if (ini_ent_num == nullptr)
    return nullptr;

  switch (val_type) {

  case CS_TYPE_char:
    {
      char  *val_ord;
      const char  *val_cur = (const char *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, char);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[(ini_ent_num[ent_id] - 1) * n_location_vals + jj]
            = val_cur[ii++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  case CS_TYPE_int:
    {
      int  *val_ord;
      const int  *val_cur = (const int *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, int);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[(ini_ent_num[ent_id] - 1) * n_location_vals + jj]
            = val_cur[ii++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  case CS_TYPE_cs_gnum_t:
    {
      cs_gnum_t  *val_ord;
      const cs_gnum_t  *val_cur = (const cs_gnum_t *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, cs_gnum_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[(ini_ent_num[ent_id] - 1) * n_location_vals + jj]
            = val_cur[ii++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  case CS_TYPE_cs_real_t:
    {
      cs_real_t  *val_ord;
      const cs_real_t  *val_cur = (const cs_real_t *)vals;

      CS_MALLOC(val_ord, n_ents * n_location_vals, cs_real_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[(ini_ent_num[ent_id] - 1) * n_location_vals + jj]
            = val_cur[ii++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  default:
    assert(0);
    return nullptr;

  }
}

/*----------------------------------------------------------------------------
 * Find a given record in an indexed restart file.
 *
 * parameters:
 *   restart   <-- associated restart file pointer
 *   prefix    <-- prefix to name of record
 *   name      <-- base name of record
 *   postfix   <-- postfix to name of record
 *
 * returns:
 *   the id assigned to the record, or -1 if not found
 *----------------------------------------------------------------------------*/

static int
_restart_section_id(cs_restart_t     *restart,
                    const char       *prefix,
                    const char       *name,
                    const char       *postfix)
{
  size_t index_size = cs_io_get_index_size(restart->fh);

  char *_sec_name = nullptr;
  const char *sec_name = name;

  int rec_id = -1;

  if (prefix != nullptr || postfix != nullptr) {


    size_t sec_name_l = strlen(name);

    if (prefix != nullptr)
      sec_name_l += strlen(prefix);
    if (postfix != nullptr)
      sec_name_l += strlen(postfix);

    CS_MALLOC(_sec_name, sec_name_l + 1, char);
    sec_name = _sec_name;

    if (prefix != nullptr) {
      strcpy(_sec_name, prefix);
      strcat(_sec_name, name);
    }
    else
      strcpy(_sec_name, name);

    if (postfix != nullptr)
      strcat(_sec_name, postfix);

  }

  /* Search for the record in the index */

  for (rec_id = 0; rec_id < (int)index_size; rec_id++) {
    const char * cmp_name = cs_io_get_indexed_sec_name(restart->fh, rec_id);
    if (strcmp(cmp_name, sec_name) == 0)
      break;
  }

  if (rec_id >= (int)index_size) {
    bft_printf(_("  %s: section \"%s\" not present.\n"),
               restart->name, sec_name);
    rec_id = -1;
  }

  CS_FREE(_sec_name);

  return rec_id;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute default particle destination rank array in case of untracked
 * particles.
 *
 * For simplicity, those particles are simply distributed among ranks
 * (as it is possible to define a global numbering based on a space-filling
 * curve when generating the restart file, this may be made "geometrically"
 * balanced also).
 *
 * parameters:
 *   p_bi         <-- pointer to particles bock distribution info
 *   p_cell_num   <-- global cell number associated with each particle
 *                    (0 for untracked particles)
 *   comm         <-- associated MPI communicator
 *
 * returns:
 *   default rank array for particles (>= 0 for untracked particles)
 *----------------------------------------------------------------------------*/

static int *
_default_p_rank(cs_block_dist_info_t  *p_bi,
                const cs_gnum_t       *p_cell_num,
                MPI_Comm               comm)
{
  cs_lnum_t i;
  cs_block_dist_info_t free_particle_bi;

  int n_ranks = 0, rank_id = -1;

  cs_lnum_t _n_particles = 0, n_free_particles = 0;
  cs_gnum_t _n_g_free_particles = 0, n_g_free_particles = 0;

  cs_lnum_t *free_particle_ids = nullptr;

  fvm_io_num_t *free_particle_io_num = nullptr;
  const cs_gnum_t *free_particle_num = nullptr;

  int *default_rank = nullptr;

  /* Initialization */

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  /* Count number of untracked particles */

  _n_particles = p_bi->gnum_range[1] - p_bi->gnum_range[0];
  n_free_particles = 0;

  for (i = 0; i < _n_particles; i++) {
    if (p_cell_num[i] == 0)
      n_free_particles += 1;
  }

  _n_g_free_particles = n_free_particles;
  MPI_Allreduce(&_n_g_free_particles, &n_g_free_particles, 1,
                CS_MPI_GNUM, MPI_SUM, comm);

  /* Return if we do not have untracked particles */

  if (n_g_free_particles == 0)
    return nullptr;

  /* Initialize rank info */

  MPI_Comm_size(comm, &n_ranks);
  MPI_Comm_size(comm, &rank_id);
  free_particle_bi = cs_block_dist_compute_sizes(rank_id,
                                                 n_ranks,
                                                 0,
                                                 0,
                                                 n_g_free_particles);

  /* Define distribution of untracked particles based on scan;
   *
   *  The main objective of this function
   *  is to ensure some measure of load balancing. */

  CS_MALLOC(default_rank, _n_particles, int);
  for (i = 0; i < _n_particles; i++)
    default_rank[i] = -1;

  CS_MALLOC(free_particle_ids, n_free_particles, cs_lnum_t);

  n_free_particles = 0;
  for (i = 0; i < _n_particles; i++) {
    if (p_cell_num[i] == 0)
      free_particle_ids[n_free_particles++] = i;
  }

  free_particle_io_num = fvm_io_num_create_from_scan(n_free_particles);
  free_particle_num = fvm_io_num_get_global_num(free_particle_io_num);

  /* Determine rank based on global numbering */
  for (i = 0; i < n_free_particles; i++) {
    default_rank[free_particle_ids[i]]
      =    ((free_particle_num[i] - 1) / free_particle_bi.block_size)
         * free_particle_bi.rank_step;
  }

  free_particle_io_num = fvm_io_num_destroy(free_particle_io_num);
  CS_FREE(free_particle_ids);

  return default_rank;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the presence of a given section in a restart file.
 *
 * \param[in]       restart          associated restart file pointer
 * \param[in, out]  context          associated context
 * \param[in]       sec_name         section name
 * \param[in]       location_id      id of corresponding location
 * \param[in]       n_location_vals  number of values per location (interlaced)
 * \param[in]       val_type         value type
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

static int
_check_section(cs_restart_t           *restart,
               void                   *context,
               const char             *sec_name,
               int                     location_id,
               int                     n_location_vals,
               cs_restart_val_type_t   val_type)
{
  CS_UNUSED(context);

  cs_lnum_t n_ents;
  cs_gnum_t n_glob_ents;

  size_t rec_id;
  cs_io_sec_header_t header;

  size_t index_size = 0;

  index_size = cs_io_get_index_size(restart->fh);

  assert(restart != nullptr);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
  }

  else {
    if (location_id < 0 || location_id > (int)(restart->n_locations))
      return CS_RESTART_ERR_LOCATION;
    n_glob_ents = (restart->location[location_id-1]).n_glob_ents;
    if ((restart->location[location_id-1]).n_glob_ents_f != n_glob_ents)
      return CS_RESTART_ERR_LOCATION;
    n_ents  = (restart->location[location_id-1]).n_ents;
  }

  /* Search for the corresponding record in the index */

  for (rec_id = 0; rec_id < index_size; rec_id++) {
    const char * cmp_name = cs_io_get_indexed_sec_name(restart->fh, rec_id);
    if (strcmp(cmp_name, sec_name) == 0)
      break;
  }

  /* If the record was not found */

  if (rec_id >= index_size)
    return CS_RESTART_ERR_EXISTS;

  /*
    If the location does not fit: we search for a location of same
    name with the correct location.
  */

  header = cs_io_get_indexed_sec_header(restart->fh, rec_id);

  if (header.location_id != (size_t)location_id) {

    rec_id++;

    while (rec_id < index_size) {
      header = cs_io_get_indexed_sec_header(restart->fh, rec_id);
      if (   (strcmp(header.sec_name, sec_name) == 0)
          && (header.location_id == (size_t)location_id))
        break;
      rec_id++;
    }

    if (rec_id >= index_size)
      return CS_RESTART_ERR_LOCATION;
  }

  /* If the number of values per location does not match */

  if (   header.location_id > 0
      && header.n_location_vals != (size_t)n_location_vals)
    return CS_RESTART_ERR_N_VALS;
  else if (header.location_id == 0 && header.n_vals != n_ents)
    return CS_RESTART_ERR_N_VALS;

  /* If the type of value does not match */

  if (header.elt_type == CS_CHAR) {
    if (val_type != CS_TYPE_char)
      return CS_RESTART_ERR_VAL_TYPE;
  }
  else if (header.elt_type == CS_INT32 || header.elt_type == CS_INT64) {
    cs_io_set_cs_lnum(&header, restart->fh);
    if (val_type != CS_TYPE_int)
      return CS_RESTART_ERR_VAL_TYPE;
  }
  else if (header.elt_type == CS_UINT32 || header.elt_type == CS_UINT64) {
    if (val_type != CS_TYPE_cs_gnum_t && val_type != CS_TYPE_int)
      return CS_RESTART_ERR_VAL_TYPE;
  }
  else if (header.elt_type == CS_FLOAT || header.elt_type == CS_DOUBLE) {
    if (val_type != CS_TYPE_cs_real_t)
      return CS_RESTART_ERR_VAL_TYPE;
  }

  /* Return */

  return CS_RESTART_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a section from a restart file.
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   location_id      id of corresponding location
 * \param[in]   n_location_vals  number of values per location (interlaced)
 * \param[in]   val_type         value type
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

static int
_read_section(cs_restart_t           *restart,
              void                   *context,
              const char             *sec_name,
              int                     location_id,
              int                     n_location_vals,
              cs_restart_val_type_t   val_type,
              void                   *val)
{
  CS_UNUSED(context);

  cs_lnum_t n_ents;
  cs_gnum_t n_glob_ents;

  const cs_gnum_t  *ent_global_num;

  size_t rec_id, rec_id_tmp;
  cs_io_sec_header_t header;

  cs_lnum_t _n_location_vals = n_location_vals;
  size_t index_size = 0;

  index_size = cs_io_get_index_size(restart->fh);

  assert(restart != nullptr);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = nullptr;
  }

  else {
    if (location_id < 0 || location_id > (int)(restart->n_locations)) {
      bft_printf(_("  %s: location id %d for \"%s\" does not exist.\n"),
                 restart->name, location_id, sec_name);
      return CS_RESTART_ERR_LOCATION;
    }
    n_glob_ents = (restart->location[location_id-1]).n_glob_ents;
    if ((restart->location[location_id-1]).n_glob_ents_f != n_glob_ents) {
      bft_printf
        (_("  %s: location id %d for \"%s\" has "
           "size %llu, but %llu is expected.\n"),
         restart->name, location_id, sec_name,
         (unsigned long long)(restart->location[location_id-1]).n_glob_ents_f,
         (unsigned long long)n_glob_ents);
      return CS_RESTART_ERR_LOCATION;
    }
    n_ents  = (restart->location[location_id-1]).n_ents;
    ent_global_num = (restart->location[location_id-1]).ent_global_num;
  }

  /* Search for the corresponding record in the index */

  for (rec_id = 0; rec_id < index_size; rec_id++) {
    const char * cmp_name = cs_io_get_indexed_sec_name(restart->fh, rec_id);
    if (strcmp(cmp_name, sec_name) == 0)
      break;
  }

  /* If the record was not found */

  if (rec_id >= index_size) {
    bft_printf(_("  %s: section \"%s\" not present.\n"),
               restart->name, sec_name);
    return CS_RESTART_ERR_EXISTS;
  }

  /*
    If the location does not fit: we search for a location of same
    name with the correct location.
  */

  header = cs_io_get_indexed_sec_header(restart->fh, rec_id);

  if (header.location_id != (size_t)location_id) {

    rec_id_tmp = rec_id;
    rec_id++;

    while (rec_id < index_size) {
      header = cs_io_get_indexed_sec_header(restart->fh, rec_id);
      if (   (strcmp(header.sec_name, sec_name) == 0)
          && (header.location_id == (size_t)location_id))
        break;
      rec_id++;
    }

    if (rec_id >= index_size) {
      header = cs_io_get_indexed_sec_header(restart->fh, rec_id_tmp);
      bft_printf(_("  %s: section \"%s\" at location id %d but not at %d.\n"),
                 restart->name, sec_name,
                 (int)(header.location_id), (int)location_id);
      return CS_RESTART_ERR_LOCATION;
    }
  }

  /* If the number of values per location does not match */

  if (   header.location_id > 0
      && header.n_location_vals != (size_t)n_location_vals) {
    bft_printf(_("  %s: section \"%s\" has %d values per location and "
                 " not %d.\n"),
               restart->name, sec_name,
               (int)header.n_location_vals, (int)n_location_vals);
    return CS_RESTART_ERR_N_VALS;
  }
  else if (header.location_id == 0 && header.n_vals != n_ents) {
    bft_printf(_("  %s: section \"%s\" has %d values and not %d.\n"),
               restart->name, sec_name, (int)header.n_vals, (int)n_ents);
    return CS_RESTART_ERR_N_VALS;
  }

  /* If the type of value does not match */

  if (header.elt_type == CS_CHAR) {
    if (val_type != CS_TYPE_char) {
      bft_printf(_("  %s: section \"%s\" is not of character type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }
  else if (header.elt_type == CS_INT32 || header.elt_type == CS_INT64) {
    cs_io_set_cs_lnum(&header, restart->fh);
    if (val_type != CS_TYPE_int) {
      bft_printf(_("  %s: section \"%s\" is not of integer type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }
  else if (header.elt_type == CS_UINT32 || header.elt_type == CS_UINT64) {
    if (val_type != CS_TYPE_cs_gnum_t && val_type != CS_TYPE_int) {
      bft_printf(_("  %s: section \"%s\" is not of global number type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }
  else if (header.elt_type == CS_FLOAT || header.elt_type == CS_DOUBLE) {
    if (val_type != CS_TYPE_cs_real_t) {
      bft_printf(_("  %s: section \"%s\" is not of floating-point type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }

  /* Now set position in file to read data */

  cs_io_set_indexed_position(restart->fh, &header, rec_id);

  /* Now define conversion info */

  if (header.elt_type == CS_UINT32 || header.elt_type == CS_UINT64) {
    if (val_type == CS_TYPE_cs_gnum_t)
      cs_io_set_cs_gnum(&header, restart->fh);
    else if (val_type == CS_TYPE_int)
      cs_io_set_cs_lnum(&header, restart->fh);
  }
  else if (header.elt_type == CS_FLOAT || header.elt_type == CS_DOUBLE) {
    if (sizeof(cs_real_t) != cs_datatype_size[header.elt_type]) {
      if (sizeof(cs_real_t) == cs_datatype_size[CS_FLOAT])
        header.elt_type = CS_FLOAT;
      else
        header.elt_type = CS_DOUBLE;
    }
  }

  /* Section contents */
  /*------------------*/

  /* In single processor mode or for global values */

  if (cs_glob_n_ranks == 1 || location_id == 0) {

    cs_io_read_global(&header, val, restart->fh);

    if (ent_global_num != nullptr)
      _restart_permute_read(n_ents,
                            ent_global_num,
                            _n_location_vals,
                            val_type,
                            (cs_byte_t *)val);
  }

#if defined(HAVE_MPI)

  /* In parallel mode for a distributed mesh location */

  else if (n_glob_ents > 0)
    _read_ent_values(restart,
                     &header,
                     n_glob_ents,
                     n_ents,
                     ent_global_num,
                     _n_location_vals,
                     val_type,
                     (cs_byte_t *)val);

#endif /* #if defined(HAVE_MPI) */

  /* Return */

  return CS_RESTART_SUCCESS;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a section to a restart file.
 *
 * \param[in]  restart          associated restart file pointer
 * \param[in]  sec_name         section name
 * \param[in]  location_id      id of corresponding location
 * \param[in]  n_location_vals  number of values per location (interlaced)
 * \param[in]  val_type         value type
 * \param[in]  val              array of values
 */
/*----------------------------------------------------------------------------*/

static void
_write_section(cs_restart_t           *restart,
               void                   *context,
               const char             *sec_name,
               int                     location_id,
               int                     n_location_vals,
               cs_restart_val_type_t   val_type,
               const void             *val)
{
  CS_UNUSED(context);

  cs_lnum_t        n_ents;
  cs_gnum_t        n_tot_vals, n_glob_ents;
  cs_datatype_t    elt_type = CS_DATATYPE_NULL;

  const cs_gnum_t  *ent_global_num;

  cs_lnum_t _n_location_vals = n_location_vals;

  assert(restart != nullptr);

  n_tot_vals = _compute_n_ents(restart, location_id, n_location_vals);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = nullptr;
  }

  else {
    assert(location_id >= 0 && location_id <= (int)(restart->n_locations));
    n_glob_ents = (restart->location[location_id-1]).n_glob_ents;
    n_ents  = (restart->location[location_id-1]).n_ents;
    ent_global_num = (restart->location[location_id-1]).ent_global_num;
  }

  /* Set val_type */

  switch (val_type) {
  case CS_TYPE_char:
    elt_type = CS_CHAR;
    break;
  case CS_TYPE_int:
    elt_type = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;
    break;
  case CS_TYPE_cs_gnum_t:
    elt_type = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
    break;
  case CS_TYPE_cs_real_t:
    elt_type =   (sizeof(cs_real_t) == cs_datatype_size[CS_DOUBLE])
               ? CS_DOUBLE : CS_FLOAT;
    break;
  default:
    assert(0);
  }

  /* Section contents */
  /*------------------*/

  /* In single processor mode of for global values */

  if (location_id == 0)
    cs_io_write_global(sec_name,
                       n_tot_vals,
                       location_id,
                       0,
                       1,
                       elt_type,
                       val,
                       restart->fh);


  else if (cs_glob_n_ranks == 1 || n_glob_ents == 0) {

    cs_byte_t  *val_tmp = nullptr;

    if (ent_global_num != nullptr)
      val_tmp = _restart_permute_write(n_ents,
                                       ent_global_num,
                                       _n_location_vals,
                                       val_type,
                                       (const cs_byte_t *)val);
    cs_io_write_global(sec_name,
                       n_tot_vals,
                       location_id,
                       0,
                       _n_location_vals,
                       elt_type,
                       (val_tmp != nullptr) ? val_tmp : val,
                       restart->fh);

    if (val_tmp != nullptr)
      CS_FREE (val_tmp);
  }

#if defined(HAVE_MPI)

  /* In parallel mode for a distributed mesh location */

  else
    _write_ent_values(restart,
                      sec_name,
                      n_glob_ents,
                      n_ents,
                      ent_global_num,
                      location_id,
                      _n_location_vals,
                      val_type,
                      (const cs_byte_t *)val);

#endif /* #if defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a location definition with a possible reference definition.
 *
 * \param[in]  restart         associated restart file pointer
 * \param[in]  location_name   name associated with the location
 * \param[in]  n_glob_ents     global number of entities
 * \param[in]  n_ents          local number of entities
 * \param[in]  ent_global_num  global entity numbers, or nullptr
 *
 * \return  the location id assigned, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

static int
_add_location_check_ref(cs_restart_t     *restart,
                        const char       *location_name,
                        cs_gnum_t         n_glob_ents,
                        cs_lnum_t         n_ents,
                        const cs_gnum_t  *ent_global_num)
{
  int loc_id;

  if (restart->mode == CS_RESTART_MODE_READ) {

    /* Search for a reference location with the same name */

    for (loc_id = 0; loc_id < (int)(_n_locations_ref); loc_id++) {

      if ((strcmp((restart->location[loc_id]).name, location_name) == 0)) {

        const _location_t  *_loc_ref = _location_ref + loc_id;

        n_glob_ents = _loc_ref->n_glob_ents;
        n_ents = _loc_ref->n_ents;
        ent_global_num =_loc_ref->ent_global_num;

      }

    }

  }

  return cs_restart_add_location(restart, location_name,
                                 n_glob_ents, n_ents,
                                 ent_global_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update checkpoint directory with mesh.
 *
 * The mesh_output file is moved to restart/mesh_input if present.
 * Otherwise, if mesh_input is present and a file (not a directory),
 * is linked to restart/mesh_input (using a hard link if possible)
 */
/*----------------------------------------------------------------------------*/

static void
_update_mesh_checkpoint(void)
{
  if (cs_glob_rank_id < 1 && _checkpoint_mesh > 0) {

    if (cs_glob_rank_id < 1) {
      const char _checkpoint[] = "checkpoint";
      if (cs_file_mkdir_default(_checkpoint) != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("The %s directory cannot be created"), _checkpoint);
    }

    /* Move mesh_output if present */

    const char *opath_i[2] = {"mesh_input.csm",
                              "restart/mesh_input.csm"};
    const char opath_o[] = "mesh_output.csm";
    const char npath[] = "checkpoint/mesh_input.csm";

    if (cs_file_isreg(opath_o)) {
      int retval = rename(opath_o, npath);
      if (retval != 0) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("Failure moving %s to %s:\n"
                     "%s\n"),
                   opath_o, npath, strerror(errno));
      }
    }

    /* Otherwise link mesh_input if it is a regular file
       (in case of a mesh_input directory, do nothing, since a
       directory should be created only in case of multiple meshes,
       and concatenation leads to a mesh_output being created,
       unless the (advanced) user has explicitely deactivated that
       output, in which case we abide by that choice) */

    else if (cs_glob_mesh->modified < 1) {

#if defined(HAVE_LINKAT) && defined(HAVE_FCNTL_H)

      for (int i = 0; i < 2; i++) {
        if (cs_file_isreg(opath_i[i])) {
          int retval = linkat(AT_FDCWD, opath_i[i],
                              AT_FDCWD, npath, AT_SYMLINK_FOLLOW);

          if (retval != 0) {
            cs_base_warn(__FILE__, __LINE__);
            bft_printf(_("Failure hard-linking %s to %s:\n"
                         "%s\n"),
                       opath_i[i], npath, strerror(errno));

          }
          break;
        }
      }

#endif

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create and allocate a new multiwriter structure.
 *
 * \return  pointer to the allocated object.
 */
/*----------------------------------------------------------------------------*/

static _restart_multiwriter_t *
_restart_multiwriter_create(void)
{
  _restart_multiwriter_t *new_writer = nullptr;
  CS_MALLOC(new_writer, 1, _restart_multiwriter_t);

  new_writer->id = -1;
  new_writer->name = nullptr;
  new_writer->path = nullptr;
  new_writer->n_prev_files = -1;  /* set at 0 after first (single) output */
  new_writer->n_prev_files_tot = 0;
  new_writer->prev_files = nullptr;

  return new_writer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get a multiwriter using its id.
 *
 * Returns a pointer to the writer object. If no writer has that id,
 * a null pointer is returned.
 *
 * \param[in] id  id of the multiwriter (int)
 *
 * \return  pointer to the multiwriter object.
 */
/*----------------------------------------------------------------------------*/

static _restart_multiwriter_t *
_restart_multiwriter_by_id(const int id)
{
  _restart_multiwriter_t *writer = nullptr;
  if (id >= 0 && id < _n_restart_multiwriters)
    writer = _restart_multiwriter[id];

  return writer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a multiwriter based on a file name and path.
 *
 * Returns the id of the multiwriter.
 *
 * \param[in] name  name of the restart file (char *)
 * \param[in] path  path to the restart file (char *)
 *
 * \return  id of the newly added multiwriter
 */
/*----------------------------------------------------------------------------*/

static int
_add_restart_multiwriter(const char  name[],
                         const char  path[])
{
  /* Check if the writer already exists */
  for (int i = 0; i < _n_restart_multiwriters; i++) {
    if (strcmp(_restart_multiwriter[i]->name, name) == 0)
      return i;
  }

  /* Check that the file name is neither null or empty ("") */
  if (name == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              _("Null pointer provided as file name.\n"));
  else if (strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty file name was provided.\n"));

  /* Allocate or reallocate the array */
  if (_n_restart_multiwriters == 0)
    CS_MALLOC(_restart_multiwriter,
              _n_restart_multiwriters + 1,
              _restart_multiwriter_t *);
  else
    CS_REALLOC(_restart_multiwriter,
               _n_restart_multiwriters + 1,
               _restart_multiwriter_t *);

  /* Create the empty structure */
  _restart_multiwriter_t *new_writer = _restart_multiwriter_create();

  /* Set id */
  new_writer->id = _n_restart_multiwriters;

  /* Set name */
  size_t lname = strlen(name) + 1;
  CS_MALLOC(new_writer->name, lname, char);
  strcpy(new_writer->name, name);

  /* Set path to subdir, which can be null */
  const char *_path = path;
  if (_path != nullptr) {
    if (strlen(_path) == 0)
      _path = nullptr;
  }

  if (_path != nullptr) {
    size_t lpath = strlen(_path) + 1;
    CS_MALLOC(new_writer->path, lpath, char);
    strcpy(new_writer->path, _path);
  }

  _restart_multiwriter[_n_restart_multiwriters] = new_writer;

  /* Increment the number of writers, and return the newly created
   * writer index */

  _n_restart_multiwriters++;

  return new_writer->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Increment the multiwriter previous files list.
 *
 * \param[in] mw     pointer to the multiwriter object.
 * \param[in] fname  name of the file to add to the old files list
 *
 */
/*----------------------------------------------------------------------------*/

static void
_restart_multiwriter_increment(_restart_multiwriter_t  *mw,
                               const char               fname[])
{
  mw->n_prev_files++;
  mw->n_prev_files_tot++;

  if (mw->prev_files == nullptr)
    CS_MALLOC(mw->prev_files, mw->n_prev_files, char *);
  else
    CS_REALLOC(mw->prev_files, mw->n_prev_files, char *);

  mw->prev_files[mw->n_prev_files - 1] = nullptr;
  size_t lenf = strlen(fname) + 1;
  CS_MALLOC(mw->prev_files[mw->n_prev_files - 1], lenf, char);
  strcpy(mw->prev_files[mw->n_prev_files - 1], fname);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query checkpoint intervals.
 *
 * \param[out]  nt_interval  if non-null, time-step interval for checkpoint
 *                             if > 0 time step interval for checkpoint
 *                             if 0, default of 4 checkpoints per run
 *                             if -1, checkpoint at end
 *                             if -2, no checkpointing
 * \param[out]  t_interval   if non-null, time value checkpoint interval;
 *                           has priority over nt_interval if > -1
 * \param[out]  wt_interval  if non-null, wall-clock interval for checkpoints;
 *                           has priority over nt_interval if > -1
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_get_intervals(int     *nt_interval,
                                    double  *t_interval,
                                    double  *wt_interval)
{
  if (nt_interval != nullptr)
     *nt_interval = _checkpoint_nt_interval;
  if (t_interval != nullptr)
     *t_interval = _checkpoint_t_interval;
  if (wt_interval != nullptr)
     *wt_interval = _checkpoint_wt_interval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define default checkpoint interval.
 *
 * \param[in]  nt_interval  if > 0 time step interval for checkpoint
 *                          if 0, default of 4 checkpoints per run
 *                          if -1, checkpoint at end
 *                          if -2, no checkpointing
 * \param[in]  t_interval   if > 0, time value interval for checkpoint
 * \param[in]  wt_interval  if > 0, wall-clock interval for checkpoints
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_interval(int     nt_interval,
                                   double  t_interval,
                                   double  wt_interval)
{
  _checkpoint_nt_interval = nt_interval;
  _checkpoint_t_interval = t_interval;
  _checkpoint_wt_interval = wt_interval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define checkpoint behavior for mesh.
 *
 * If mesh checkpointing is active (default), upon writing the first
 * checkpoint file of a computation, a mesh_output file is moved to
 * checkpoint/mesh_input if present. If not present but a mesh_input
 * file (or link to file) is present, a hard link to that file is
 * added as checkpoint/mesh_input.
 *
 * A mesh_input directory is ignored, as it is normally only created when
 * multiple input files are appended, which leads to the output and thus
 * presence of a mesh_output file, unless explicitely deactivated by the
 * user.
 *
 * \param[in]  mode  if 0, do not checkpoint mesh
 *                   if 1, checkpoint mesh_output or mesh_input file
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_mesh_mode(int  mode)
{
  _checkpoint_mesh = mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define last forced checkpoint time step.
 *
 * \param[in]  nt_last  last time step for forced checkpoint
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_last_ts(int  nt_last)
{
  _checkpoint_nt_last = nt_last;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define next forced checkpoint time step.
 *
 * \param[in]  nt_next  next time step for forced checkpoint
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_ts(int  nt_next)
{
  _checkpoint_nt_next = nt_next;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define next forced checkpoint time value.
 *
 * \param[in]  t_next  next time value for forced checkpoint
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_tv(double  t_next)
{
  _checkpoint_t_next = t_next;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define next forced checkpoint wall-clock time value.
 *
 * \param[in]  wt_next  next wall-clock time value for forced checkpoint
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_set_next_wt(double  wt_next)
{
  _checkpoint_wt_next = wt_next;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if checkpointing is recommended at a given time.
 *
 * \param[in]  ts  time step status structure
 *
 * \return  true if checkpointing is recommended, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_restart_checkpoint_required(const cs_time_step_t  *ts)
{
  assert(ts != nullptr);

  int nt = ts->nt_cur - ts->nt_prev;
  double t = ts->t_cur - ts->t_prev;

  bool retval = false;

  if (_checkpoint_nt_interval > CS_RESTART_INTERVAL_NONE) {

    if (ts->nt_cur == ts->nt_max)    /* Output at the last time step */
      retval = true;

    else if (_checkpoint_nt_interval == CS_RESTART_INTERVAL_DEFAULT) {
      /* default interval: current number of expected time_steps for this run,
         with a minimum of 10. */
      int nt_def = (ts->nt_max - ts->nt_prev)/4;
      if (nt_def < 10)
        nt_def = 10;
      if (nt % nt_def == 0)
        retval = true;
    }

    else if (_checkpoint_nt_interval > CS_RESTART_INTERVAL_DEFAULT &&
             nt % _checkpoint_nt_interval == 0)
      retval = true;

    else if (_checkpoint_nt_interval > CS_RESTART_INTERVAL_DEFAULT &&
             _checkpoint_nt_last > -1) {
      if (ts->nt_cur >= _checkpoint_nt_interval + _checkpoint_nt_last)
        retval = true;
    }

  }

  if (_checkpoint_t_interval > 0
      && _checkpoint_t_last + _checkpoint_t_interval <= t)
    retval = true;

  else if (_checkpoint_wt_next >= 0) {
    double wt = cs_timer_wtime();
    if (wt >= _checkpoint_wt_next)
      retval = true;
  }

  else if (   (_checkpoint_nt_next >= 0 && _checkpoint_nt_next <= ts->nt_cur)
           || (_checkpoint_t_next >= 0 && _checkpoint_t_next <= ts->t_cur))
    retval = true;

  else if (_checkpoint_wt_interval >= 0) {
    double wt = cs_timer_wtime();
    cs_parall_bcast(0, 1, CS_DOUBLE, &wt);
    if (wt - _checkpoint_wt_last >= _checkpoint_wt_interval)
      retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate checkpointing has been done at a given time.
 *
 * This updates the status for future checks to determine
 * if checkpointing is recommended at a given time.
 *
 * \param[in]  ts  time step status structure
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_checkpoint_done(const cs_time_step_t  *ts)
{
  assert(ts != nullptr);

  double t = ts->t_cur - ts->t_prev;

  if (_checkpoint_nt_next >= 0 && _checkpoint_nt_next <= ts->nt_cur)
    _checkpoint_nt_next = -1;

  if (_checkpoint_t_next >= 0 && _checkpoint_t_next <= ts->t_cur)
    _checkpoint_t_next = -1.;

  if (_checkpoint_wt_next >= 0) {
    double wt = cs_timer_wtime();
    if (wt >= _checkpoint_wt_next)
      _checkpoint_wt_next = -1.;
  }

  if (_checkpoint_t_interval > 0
      && _checkpoint_t_last + _checkpoint_t_interval <= t)
    _checkpoint_t_last = ts->t_cur;

  if (_checkpoint_wt_interval >= 0) {
    double wt = cs_timer_wtime();
    cs_parall_bcast(0, 1, CS_DOUBLE, &wt);
    if (wt - _checkpoint_wt_last >= _checkpoint_wt_interval)
      _checkpoint_wt_last = wt;
  }
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Access raw restart data serialized in memory.
 *
 * If called previously, reinitialize memory data structure.
 *
 * \param[out]  nb    size of data
 * \param[out]  data  pointer to data
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_get_from_memory_serialized(size_t   *nb,
                                      void    **data)
{
  *nb = 0;
  *data = nullptr;

  cs_io_t *r_io = _checkpoint_serialized_memory;
  if (r_io != nullptr)
    cs_io_get_data_in_mem(r_io, nb, data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate restart will be done based on a serialized data in memory.
 *
 * The restart subsystem takes ownership of the given data
 *
 * \param[in]  nb    number of matching bytes for data
 * \param[in]  data  data buffer (ownership is relinquished by caller)
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_set_from_memory_serialized(size_t   nb,
                                      void    *data)
{
  const char magic_string[] = "Checkpoint / restart, R0";
  const long echo = CS_IO_ECHO_NONE;

  if (_restart_serialized_memory != nullptr)
    cs_io_finalize(&_restart_serialized_memory);

#if defined(HAVE_MPI)
  {
    MPI_Comm  block_comm, comm;

    cs_file_get_default_comm(nullptr, &block_comm, &comm);
    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);

    _restart_serialized_memory
      = cs_io_initialize_with_index_from_mem("restart",
                                             magic_string,
                                             CS_FILE_IN_MEMORY_SERIAL,
                                             echo,
                                             nb, data,
                                             block_comm, comm);
  }
#else
  _restart_serialized_memory
    = cs_io_initialize_with_index_from_mem("restart",
                                           magic_string,
                                           CS_FILE_IN_MEMORY_SERIAL,
                                           echo,
                                           nb, data);
#endif

  if (_need_finalize == false) {
    _need_finalize = true;
    cs_base_at_finalize(_restart_finalize);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Indicate checkpoint will be done to serialized data in memory.
 *
 * If called previously, reinitialize memory data structure.
 *
 * \param[in]  status  checkpoint to memory if true, to file otherwise.
 */
/*----------------------------------------------------------------------------*/

void
cs_checkpoint_set_to_memory_serialized(bool  status)
{
  if (_checkpoint_serialized_memory != nullptr)
    cs_io_finalize(&_checkpoint_serialized_memory);

  if (status) {

    const char magic_string[] = "Checkpoint / restart, R0";
    const long echo = CS_IO_ECHO_NONE;

#if defined(HAVE_MPI)

    MPI_Comm           block_comm, comm;
    cs_file_get_default_comm(nullptr, &block_comm, &comm);

    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);

    _checkpoint_serialized_memory = cs_io_initialize("checkpoint",
                                                     magic_string,
                                                     CS_IO_MODE_WRITE,
                                                     CS_FILE_IN_MEMORY_SERIAL,
                                                     echo,
                                                     MPI_INFO_NULL,
                                                     block_comm,
                                                     comm);
#else

    _checkpoint_serialized_memory = cs_io_initialize("checkpoint",
                                                     magic_string,
                                                     CS_IO_MODE_WRITE,
                                                     CS_FILE_IN_MEMORY_SERIAL,
                                                     echo);

#endif

  }

  if (_need_finalize == false) {
    _need_finalize = true;
    cs_base_at_finalize(_restart_finalize);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if we have a restart directory.
 *
 * \return  1 if a restart directory is present, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_present(void)
{
  if (_restart_present < 0) {
    if (cs_glob_rank_id < 1) {
      if (  _restart_serialized_memory != nullptr
          || cs_file_isdir("restart"))
        _restart_present = 1;
      else
        _restart_present = 0;
    }
    cs_parall_bcast(0, 1, CS_INT_TYPE, &_restart_present);
  }

  return _restart_present;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a restart file.
 *
 * \param[in]  name  file name
 * \param[in]  path  optional directory name for output, or null for default
 *                   (directory automatically created if necessary)
 * \param[in]  mode  read or write
 *
 * \returns  pointer to initialized restart file structure
 */
/*----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_create(const char         *name,
                  const char         *path,
                  cs_restart_mode_t   mode)
{
  cs_restart_t  * restart;

  double timing[2];

  char *_name = nullptr;
  size_t  ldir, lname, lext;

  const char  *_path = path;
  const char _restart[] = "restart";
  const char _checkpoint[] = "checkpoint";
  const char _extension[]  = ".csc";

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Ensure mesh checkpoint is updated on first call */

  if (    mode == CS_RESTART_MODE_WRITE
      && _restart_n_opens[mode] == 0) {
    _update_mesh_checkpoint();
  }

  /* Initializations */

  timing[0] = cs_timer_wtime();

  if (_path != nullptr) {
    if (strlen(_path) == 0)
      _path = nullptr;
  }

  if (_path == nullptr) {
    if (mode == CS_RESTART_MODE_WRITE)
      _path = _checkpoint;
    else
      _path = _restart;
  }

  /* Create 'checkpoint' directory or read from 'restart' directory */

  if (cs_glob_rank_id < 1 && _restart_serialized_memory == nullptr) {
    if (mode == CS_RESTART_MODE_WRITE) {
      if (cs_file_mkdir_default(_path) != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("The %s directory cannot be created"), _path);
    }
    else if (mode == CS_RESTART_MODE_READ) {
      if (cs_file_isdir(_path) == 0) {
        bft_error(__FILE__, __LINE__, 0,
                  _("The %s directory cannot be found"), _path);
      }
    }
  }
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Barrier(cs_glob_mpi_comm);
#endif

  ldir = strlen(_path);
  lname = strlen(name);

  CS_MALLOC(_name, ldir + lname + 2, char);

  strcpy(_name, _path);
  _name[ldir] = _dir_separator;
  _name[ldir+1] = '\0';
  strcat(_name, name);
  _name[ldir+lname+1] = '\0';

  /* Following the addition of an extension, we check for READ mode
   * if a file exists without the extension */

  if (   mode == CS_RESTART_MODE_READ
      && _restart_serialized_memory == nullptr) {

    if (cs_file_isreg(_name) == 0 && cs_file_endswith(name, _extension)) {
      CS_FREE(_name);

      lext  = strlen(_extension);
      CS_MALLOC(_name, ldir + lname + 2 - lext, char);
      strcpy(_name, _path);
      _name[ldir] = _dir_separator;
      _name[ldir+1] = '\0';
      strncat(_name, name, lname-lext);
      _name[ldir+lname-lext+1] = '\0';
    }

  }
  else if (   mode == CS_RESTART_MODE_WRITE
           && _restart_serialized_memory == nullptr) {

    /* Check if file already exists, and if so rename and delete if needed */
    int writer_id = _add_restart_multiwriter(name, _name);
    _restart_multiwriter_t *mw = _restart_multiwriter_by_id(writer_id);

    /* Rename an already existing file */
    if (cs_file_isreg(_name) && mw->n_prev_files > -1) {

      char _subdir[19];
      sprintf(_subdir, "previous_dump_%04d", mw->n_prev_files_tot);
      size_t lsdir = strlen(_subdir);

      char *_re_name = nullptr;
      CS_MALLOC(_re_name, ldir + lsdir + lname + 3, char);

      strcpy(_re_name, _path);
      _re_name[ldir] = _dir_separator;
      _re_name[ldir+1] = '\0';

      strcat(_re_name, _subdir);

      /* Check that the sub-directory exists or can be created */
      if (cs_glob_rank_id < 1) {
        if (cs_file_mkdir_default(_re_name) != 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("The %s directory cannot be created"), _re_name);
      }
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1)
        MPI_Barrier(cs_glob_mpi_comm);
#endif

      _re_name[ldir+lsdir+1] = _dir_separator;
      _re_name[ldir+lsdir+2] = '\0';
      strcat(_re_name, name);
      _re_name[ldir+lsdir+lname+2] = '\0';

      rename(_name, _re_name);

      _restart_multiwriter_increment(mw, _re_name);

      CS_FREE(_re_name);
    }
    else
      mw->n_prev_files = 0;
  }

  /* Allocate and initialize base structure */

  CS_MALLOC(restart, 1, cs_restart_t);

  CS_MALLOC(restart->name, strlen(_name) + 1, char);

  strcpy(restart->name, _name);

  CS_FREE(_name);

  /* Initialize other fields */

  restart->mode = mode;

  restart->fh = nullptr;

  restart->rank_step = 1;
  restart->min_block_size = 0;

  /* Initialize location data */

  restart->n_locations = 0;
  restart->location = nullptr;

  /* Open associated file, and build an index of sections in read mode */

  if (mode == CS_RESTART_MODE_READ) {
    if (_restart_serialized_memory == nullptr)
      _add_file(restart);
    else {
      restart->fh = _restart_serialized_memory;
      restart->rank_step = 1;
      restart->min_block_size = cs_parall_get_min_coll_buf_size();
      _locations_from_index(restart);
    }
  }

  else {
    if (_checkpoint_serialized_memory == nullptr)
      _add_file(restart);
    else {
      restart->fh = _checkpoint_serialized_memory;
      restart->rank_step = 1;
      restart->min_block_size = cs_parall_get_min_coll_buf_size();
    }
  }

  /* Add basic location definitions */

  _add_location_check_ref(restart, "cells",
                          mesh->n_g_cells, mesh->n_cells,
                          mesh->global_cell_num);
  _add_location_check_ref(restart, "interior_faces",
                          mesh->n_g_i_faces, mesh->n_i_faces,
                          mesh->global_i_face_num);
  _add_location_check_ref(restart, "boundary_faces",
                          mesh->n_g_b_faces, mesh->n_b_faces,
                          mesh->global_b_face_num);
  _add_location_check_ref(restart, "vertices",
                          mesh->n_g_vertices, mesh->n_vertices,
                          mesh->global_vtx_num);

  timing[1] = cs_timer_wtime();
  _restart_wtime[mode] += timing[1] - timing[0];

  return restart;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy structure associated with a restart file (and close the file).
 *
 * \param[in, out]  restart  pointer to restart file structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_destroy(cs_restart_t  **restart)
{
  cs_restart_mode_t   mode;

  cs_restart_t *r = *restart;

  double timing[2];

  timing[0] = cs_timer_wtime();

  assert(restart != nullptr);

  mode = r->mode;

  if (   r->fh != nullptr
      && r->fh != _checkpoint_serialized_memory
      && r->fh != _restart_serialized_memory)
    cs_io_finalize(&(r->fh));

  /* Free locations array */

  if (r->n_locations > 0) {
    size_t loc_id;
    for (loc_id = 0; loc_id < r->n_locations; loc_id++) {
      CS_FREE((r->location[loc_id]).name);
      CS_FREE((r->location[loc_id])._ent_global_num);
    }
  }
  if (r->location != nullptr)
    CS_FREE(r->location);

  /* Free remaining memory */

  CS_FREE(r->name);

  CS_FREE(*restart);

  timing[1] = cs_timer_wtime();
  _restart_wtime[mode] += timing[1] - timing[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the locations associated with a restart file.
 *
 * For each type of entity, the corresponding flag is set to true if the
 * associated number of entities matches the current value (and so that we
 * consider the mesh locations are the same), false otherwise.
 *
 * \param[out]  restart       associated restart file pointer
 * \param[out]  match_cell    matching cells flag
 * \param[out]  match_i_face  matching interior faces flag
 * \param[out]  match_b_face  matching boundary faces flag
 * \param[out]  match_vertex  matching vertices flag
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_check_base_location(const cs_restart_t  *restart,
                               bool                *match_cell,
                               bool                *match_i_face,
                               bool                *match_b_face,
                               bool                *match_vertex)
{
  size_t location_id;

  *match_cell = false;
  *match_i_face = false;
  *match_b_face = false;
  *match_vertex = false;

  assert(restart != nullptr);

  for (location_id = 0; location_id < 4; location_id++) {

    const _location_t *loc = restart->location + location_id;

    if (loc->n_glob_ents_f == loc->n_glob_ents) {
      if (location_id == 0)
        *match_cell = true;
      else if (location_id == 1)
        *match_i_face = true;
      else if (location_id == 2)
        *match_b_face = true;
      else if (location_id == 3)
        *match_vertex = true;
    }

    else if (cs_glob_rank_id <= 0) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The size of location \"%s\" associated with\n"
                   "the restart file \"%s\" is %llu and does not\n"
                   "correspond to that of the current mesh (%llu).\n"),
                 loc->name, restart->name,
                 (unsigned long long)loc->n_glob_ents_f,
                 (unsigned long long)loc->n_glob_ents);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a location definition.
 *
 * \param[in]  restart         associated restart file pointer
 * \param[in]  location_name   name associated with the location
 * \param[in]  n_glob_ents     global number of entities
 * \param[in]  n_ents          local number of entities
 * \param[in]  ent_global_num  global entity numbers, or nullptr
 *
 * \return  the location id assigned, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_add_location(cs_restart_t     *restart,
                        const char       *location_name,
                        cs_gnum_t         n_glob_ents,
                        cs_lnum_t         n_ents,
                        const cs_gnum_t  *ent_global_num)
{
  double timing[2];

  int loc_id;

  timing[0] = cs_timer_wtime();

  if (restart->mode == CS_RESTART_MODE_READ) {

    /* Search for a location with the same name */

    for (loc_id = 0; loc_id < (int)(restart->n_locations); loc_id++) {

      if ((strcmp((restart->location[loc_id]).name, location_name) == 0)) {

        (restart->location[loc_id]).n_glob_ents = n_glob_ents;

        (restart->location[loc_id]).n_ents  = n_ents;
        (restart->location[loc_id]).ent_global_num = ent_global_num;
        (restart->location[loc_id])._ent_global_num = nullptr;

        timing[1] = cs_timer_wtime();
        _restart_wtime[restart->mode] += timing[1] - timing[0];

        return loc_id + 1;

      }
    }

    if (loc_id >= ((int)(restart->n_locations)))
      bft_error(__FILE__, __LINE__, 0,
                _("The restart file \"%s\" references no\n"
                  "location named \"%s\"."),
                restart->name, location_name);

  }

  else {

    cs_datatype_t gnum_type
      = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;

    /* Create a new location */

    restart->n_locations += 1;

    CS_REALLOC(restart->location, restart->n_locations, _location_t);
    CS_MALLOC((restart->location[restart->n_locations-1]).name,
              strlen(location_name)+1,
              char);

    strcpy((restart->location[restart->n_locations-1]).name, location_name);

    (restart->location[restart->n_locations-1]).id = restart->n_locations;
    (restart->location[restart->n_locations-1]).n_glob_ents    = n_glob_ents;
    (restart->location[restart->n_locations-1]).n_glob_ents_f  = n_glob_ents;
    (restart->location[restart->n_locations-1]).n_ents         = n_ents;
    (restart->location[restart->n_locations-1]).ent_global_num = ent_global_num;
    (restart->location[restart->n_locations-1])._ent_global_num = nullptr;

    cs_io_write_global(location_name, 1, restart->n_locations, 0, 0,
                       gnum_type, &n_glob_ents,
                       restart->fh);

    timing[1] = cs_timer_wtime();
    _restart_wtime[restart->mode] += timing[1] - timing[0];

    return restart->n_locations;
  }

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a reference location definition with a private copy.
 *
 * \param[in]  location_name   name associated with the location
 * \param[in]  n_glob_ents     global number of entities
 * \param[in]  n_ents          local number of entities
 * \param[in]  ent_global_num  global entity numbers, or nullptr
 *
 * \return  the location id assigned, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_add_location_ref(const char       *location_name,
                            cs_gnum_t         n_glob_ents,
                            cs_lnum_t         n_ents,
                            const cs_gnum_t  *ent_global_num)
{
  /* Create a new location */

  _n_locations_ref += 1;

  CS_REALLOC(_location_ref, _n_locations_ref, _location_t);
  CS_MALLOC((_location_ref[_n_locations_ref-1]).name,
            strlen(location_name)+1,
            char);

  strcpy((_location_ref[_n_locations_ref-1]).name, location_name);

  if (ent_global_num != nullptr) {
    CS_MALLOC((_location_ref[_n_locations_ref-1])._ent_global_num,
              n_ents, cs_gnum_t);
    for (cs_lnum_t i = 0; i < n_ents; i++) {
      (_location_ref[_n_locations_ref-1])._ent_global_num[i]
        = ent_global_num[i];
    }
  }
  else
    (_location_ref[_n_locations_ref-1])._ent_global_num = nullptr;

  (_location_ref[_n_locations_ref-1]).id = _n_locations_ref;
  (_location_ref[_n_locations_ref-1]).n_glob_ents    = n_glob_ents;
  (_location_ref[_n_locations_ref-1]).n_glob_ents_f  = n_glob_ents;
  (_location_ref[_n_locations_ref-1]).n_ents         = n_ents;
  (_location_ref[_n_locations_ref-1]).ent_global_num
    = (_location_ref[_n_locations_ref-1])._ent_global_num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Clear reference location definitions with a private copy.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_clear_locations_ref(void)
{
  for (size_t loc_id = 0; loc_id < _n_locations_ref; loc_id++) {
    CS_FREE((_location_ref[loc_id]).name);
    CS_FREE((_location_ref[loc_id])._ent_global_num);
  }
  CS_FREE(_location_ref);
  _n_locations_ref = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a context to restart section check operations.
 *
 * This context may be used by the \ref cs_restart_check_section_t,
 * \ref cs_restart_read_section_t, and \ref cs_restart_write_section_t
 * type functions.
 *
 * Note that the lifecycle of the data pointed to must be handled separately
 * (and the pointer must remain valid until the matching restart structure
 * is destroyed.
 *
 * \param[in]  context  pointer to associated data, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_set_context(void  *context)
{
  _restart_context = context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a function to restart section check operations.
 *
 * This allows defining alternate operations when checking restart sections.
 *
 * \param[in]  func   associated function
 *
 * \return   pointer to previous function
 */
/*----------------------------------------------------------------------------*/

cs_restart_check_section_t  *
cs_restart_set_check_section_func(cs_restart_check_section_t  *func)
{
  cs_restart_check_section_t  *p = _check_section_f;
  _check_section_f = func;

  return p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a function and its input to all restart section
 *         read operations.
 *
 * This allows defining alternate operations when reading restart sections.
 *
 * \param[in]  func   associated function
 *
 * \return   pointer to previous function
 */
/*----------------------------------------------------------------------------*/

cs_restart_read_section_t  *
cs_restart_set_read_section_func(cs_restart_read_section_t  *func)
{
  cs_restart_read_section_t  *p = _read_section_f;
  _read_section_f = func;

  return p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a function and its input to all restart section
 *         write operations.
 *
 * This allows defining alternate operations when writing restart sections.
 *
 * \param[in]  func     associated hook function
 *
 * \return   pointer to previous function
 */
/*----------------------------------------------------------------------------*/

cs_restart_write_section_t  *
cs_restart_set_write_section_func(cs_restart_write_section_t  *func)
{
  cs_restart_write_section_t  *p = _write_section_f;
  _write_section_f = func;

  return p;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return name of restart file
 *
 * \param[in]  restart  associated restart file pointer
 *
 * \return  base name of restart file
 */
/*----------------------------------------------------------------------------*/

const char *
cs_restart_get_name(const cs_restart_t  *restart)
{
  assert(restart != nullptr);

  return restart->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return local number of elements associated with a
 *         given restart location.
 *
 * \param[in]  restart      associated restart file pointer
 * \param[in]  location_id  id of corresponding location
 *
 * \return  number of elements associated with location.
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_restart_get_n_location_elts(const cs_restart_t  *restart,
                               int                  location_id)
{
  cs_lnum_t retval = 0;

  assert(restart != nullptr);

  if (location_id > 0 && location_id <= (int)(restart->n_locations))
    retval = restart->location[location_id-1].n_ents;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the index associated with a restart file in read mode
 *
 * \param[in]  restart  associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_dump_index(const cs_restart_t  *restart)
{
  size_t loc_id;

  assert(restart != nullptr);

  for (loc_id = 0; loc_id < restart->n_locations; loc_id++) {
    const _location_t *loc = &(restart->location[loc_id]);
    bft_printf(_("  Location: %s\n"
                 "    (number: %03d, n_glob_ents: %llu)\n"),
               loc->name, (int)(loc->id),
               (unsigned long long)(loc->n_glob_ents));
  }
  if (restart->n_locations > 0)
    bft_printf("\n");

  /* Dump general file info, including index */

  bft_printf(_("  General information associated with the restart file:\n"));

  cs_io_dump(restart->fh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the presence of a given section in a restart file.
 *
 * \param[in]  restart          associated restart file pointer
 * \param[in]  sec_name         section name
 * \param[in]  location_id      id of corresponding location
 * \param[in]  n_location_vals  number of values per location (interlaced)
 * \param[in]  val_type         value type
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_check_section(cs_restart_t           *restart,
                         const char             *sec_name,
                         int                     location_id,
                         int                     n_location_vals,
                         cs_restart_val_type_t   val_type)
{
  assert(restart != nullptr);

  return _check_section_f(restart,
                          _restart_context,
                          sec_name,
                          location_id,
                          n_location_vals,
                          val_type);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a section from a restart file.
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   location_id      id of corresponding location
 * \param[in]   n_location_vals  number of values per location (interlaced)
 * \param[in]   val_type         value type
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_section(cs_restart_t           *restart,
                        const char             *sec_name,
                        int                     location_id,
                        int                     n_location_vals,
                        cs_restart_val_type_t   val_type,
                        void                   *val)
{
  double timing[2];

  timing[0] = cs_timer_wtime();

  assert(restart != nullptr);

  int retval = _read_section_f(restart,
                               _restart_context,
                               sec_name,
                               location_id,
                               n_location_vals,
                               val_type,
                               val);

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  /* Return */

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a section to a restart file.
 *
 * \param[in]  restart          associated restart file pointer
 * \param[in]  sec_name         section name
 * \param[in]  location_id      id of corresponding location
 * \param[in]  n_location_vals  number of values per location (interlaced)
 * \param[in]  val_type         value type
 * \param[in]  val              array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_section(cs_restart_t           *restart,
                         const char             *sec_name,
                         int                     location_id,
                         int                     n_location_vals,
                         cs_restart_val_type_t   val_type,
                         const void             *val)
{
  double timing[2];

  timing[0] = cs_timer_wtime();

  assert(restart != nullptr);

  _write_section_f(restart,
                   _restart_context,
                   sec_name,
                   location_id,
                   n_location_vals,
                   val_type,
                   val);

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read basic particles information from a restart file.
 *
 * This includes building a matching location and associated global numbering.
 *
 * \param[in]   restart      associated restart file pointer
 * \param[in]   name         name of particles set
 * \param[out]  n_particles  number of particles, or nullptr
 *
 * \return  the location id assigned to the particles, or -1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_particles_info(cs_restart_t  *restart,
                               const char    *name,
                               cs_lnum_t     *n_particles)
{
  double timing[2];

  cs_lnum_t c_id;
  int rec_id;
  cs_io_sec_header_t  header;

  cs_lnum_t  block_buf_size = 0;
  size_t  nbr_byte_ent = sizeof(cs_gnum_t);
  cs_lnum_t n_cells = (restart->location[CS_MESH_LOCATION_CELLS-1]).n_ents;
  cs_gnum_t n_glob_cells
    = (restart->location[CS_MESH_LOCATION_CELLS-1]).n_glob_ents;
  cs_gnum_t n_glob_particles = 0;

  int  *default_p_rank = nullptr;
  const cs_gnum_t  *g_cell_num
    = restart->location[CS_MESH_LOCATION_CELLS-1].ent_global_num;
  const cs_datatype_t int_type
    = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  int loc_id = -1;

  timing[0] = cs_timer_wtime();

  if (n_particles != nullptr)
    *n_particles = 0;

  /* Search for location with the same name */

  for (loc_id = 0; loc_id < (int)(restart->n_locations); loc_id++) {
    if ((strcmp((restart->location[loc_id]).name, name) == 0))
      break;
  }

  if (loc_id >= (int)(restart->n_locations))
    return -1;

  n_glob_particles = (restart->location[loc_id]).n_glob_ents_f;

  /* Search for the corresponding cell_num record in the index */

  rec_id = _restart_section_id(restart, nullptr, name, "_cell_num");

  if (rec_id < 0)
    return -1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int  *b_cell_rank, *p_cell_rank;
    cs_gnum_t  *part_cell_num = nullptr;
    cs_part_to_block_t *pbd = nullptr;
    cs_all_to_all_t *d = nullptr;

    /* Now read matching cell_num to an arbitrary block distribution */

    cs_block_dist_info_t cell_bi
      = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                    cs_glob_n_ranks,
                                    restart->rank_step,
                                    restart->min_block_size / nbr_byte_ent,
                                    n_glob_cells);

    cs_block_dist_info_t part_bi
      = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                    cs_glob_n_ranks,
                                    restart->rank_step,
                                    restart->min_block_size / nbr_byte_ent,
                                    n_glob_particles);

    /* Read data to blocks */

    block_buf_size = (part_bi.gnum_range[1] - part_bi.gnum_range[0]);

    if (block_buf_size > 0)
      CS_MALLOC(part_cell_num, block_buf_size, cs_gnum_t);

    header = cs_io_get_indexed_sec_header(restart->fh, rec_id);

    cs_io_set_indexed_position(restart->fh, &header, rec_id);

    cs_io_read_block(&header,
                     part_bi.gnum_range[0],
                     part_bi.gnum_range[1],
                     part_cell_num,
                     restart->fh);

    /* Build block distribution cell rank info */

    CS_MALLOC(b_cell_rank,
              (cell_bi.gnum_range[1] - cell_bi.gnum_range[0]),
              int);

    CS_MALLOC(p_cell_rank, n_cells, int);

    pbd = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                          cell_bi,
                                          n_cells,
                                          g_cell_num);

    for (c_id = 0; c_id < n_cells; c_id++)
      p_cell_rank[c_id] = cs_glob_rank_id;

    cs_part_to_block_copy_array(pbd,
                                int_type,
                                1,
                                p_cell_rank,
                                b_cell_rank);

    cs_part_to_block_destroy(&pbd);

    CS_FREE(p_cell_rank);

    /* Now build distribution structure */

    default_p_rank = _default_p_rank(&part_bi,
                                     part_cell_num,
                                     cs_glob_mpi_comm);

    cs_lnum_t  n_part_ents = 0;
    cs_gnum_t  *ent_global_num = nullptr;

    d = cs_block_to_part_create_by_adj_s(cs_glob_mpi_comm,
                                         part_bi,
                                         cell_bi,
                                         1,
                                         part_cell_num,
                                         b_cell_rank,
                                         default_p_rank,
                                         &n_part_ents,
                                         &ent_global_num);

    if (default_p_rank != nullptr)
      CS_FREE(default_p_rank);

    CS_FREE(b_cell_rank);

    (restart->location[loc_id])._ent_global_num = ent_global_num;
    (restart->location[loc_id]).ent_global_num
      = (restart->location[loc_id])._ent_global_num;

    (restart->location[loc_id]).n_glob_ents = n_glob_particles;
    (restart->location[loc_id]).n_ents = n_part_ents;

    cs_all_to_all_destroy(&d);

    CS_FREE(part_cell_num);

  }

#endif /* #if defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {

    (restart->location[loc_id]).n_glob_ents = n_glob_particles;
    (restart->location[loc_id]).n_ents = n_glob_particles;

  }

  if (n_particles != nullptr)
    *n_particles = (restart->location[loc_id]).n_ents;

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  return loc_id + 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read basic particles information from a restart file.
 *
 * \param[in]   restart                associated restart file pointer
 * \param[in]   particles_location_id  location id of particles set
 * \param[out]  particle_cell_id       local cell id to which particles belong
 * \param[out]  particle_coords        local particle coordinates (interleaved)
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_particles(cs_restart_t  *restart,
                          int            particles_location_id,
                          cs_lnum_t     *particle_cell_id,
                          cs_real_t     *particle_coords)
{
  double timing[2];
  char *sec_name = nullptr;

  cs_lnum_t n_cells = (restart->location[CS_MESH_LOCATION_CELLS-1]).n_ents;
  const cs_gnum_t *g_cells_num
    = (restart->location[CS_MESH_LOCATION_CELLS-1]).ent_global_num;

  const char *name = (restart->location[particles_location_id - 1]).name;
  const char *cell_num_postfix = "_cell_num";
  const char *coords_postfix = "_coords";

  int retcode = CS_RESTART_SUCCESS;

  cs_lnum_t  n_particles = (restart->location[particles_location_id - 1]).n_ents;

  /* Read particle coordinates */

  CS_MALLOC(sec_name, strlen(name) + strlen(coords_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, coords_postfix);

  retcode = cs_restart_read_section(restart,
                                    sec_name,
                                    particles_location_id,
                                    3,
                                    CS_TYPE_cs_real_t,
                                    particle_coords);

  CS_FREE(sec_name);

  if (retcode != CS_RESTART_SUCCESS)
    return retcode;

  /* Read particle cell id */

  CS_MALLOC(sec_name, strlen(name) + strlen(cell_num_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, cell_num_postfix);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t *g_part_cell_num;

    CS_MALLOC(g_part_cell_num, n_particles, cs_gnum_t);

    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      particles_location_id,
                                      1,
                                      CS_TYPE_cs_gnum_t,
                                      g_part_cell_num);

    timing[0] = cs_timer_wtime();

    cs_block_to_part_global_to_local(n_particles,
                                     0,
                                     n_cells,
                                     false,
                                     g_cells_num,
                                     g_part_cell_num,
                                     particle_cell_id);

    CS_FREE(g_part_cell_num);

    timing[1] = cs_timer_wtime();
    _restart_wtime[restart->mode] += timing[1] - timing[0];

  }

#endif /* #if defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      particles_location_id,
                                      1,
                                      CS_TYPE_int,
                                      particle_cell_id);
    for (cs_lnum_t i = 0; i < n_particles; i++)
      particle_cell_id[i] -= 1;
  }

  CS_FREE(sec_name);

  return retcode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write basic particles information to a restart file.
 *
 * This includes defining a matching location and associated global numbering,
 * then writing particle coordinates and cell ids.
 *
 * \param[in]  restart            associated restart file pointer
 * \param[in]  name               name of particles set
 * \param[in]  number_by_coords   if true, numbering is based on current
 *                                coordinates; otherwise, it is simply based
 *                                on local numbers, plus the sum of particles
 *                                on lower MPI ranks
 * \param[in]  n_particles        local number of particles
 * \param[in]  particle_cell_id   local cell id (0 to n-1) to which particles
 *                                belong
 * \param[in]  particle_coords    local particle coordinates (interleaved)
 *
 * \return  the location id assigned to the particles
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_write_particles(cs_restart_t     *restart,
                           const char       *name,
                           bool              number_by_coords,
                           cs_lnum_t         n_particles,
                           const cs_lnum_t  *particle_cell_id,
                           const cs_real_t  *particle_coords)
{
  double timing[2];
  cs_lnum_t i;

  cs_gnum_t n_glob_particles = n_particles;
  cs_gnum_t  *global_particle_num = nullptr;
  cs_gnum_t  *global_part_cell_num = nullptr;
  fvm_io_num_t  *io_num = nullptr;
  char *sec_name = nullptr;

  const char *cell_num_postfix = "_cell_num";
  const char *coords_postfix = "_coords";

  int loc_id = -1;

  timing[0] = cs_timer_wtime();

  /* Build global numbering */

  cs_parall_counter(&n_glob_particles, 1);

  if (number_by_coords)
    io_num = fvm_io_num_create_from_sfc(particle_coords,
                                        3,
                                        n_particles,
                                        FVM_IO_NUM_SFC_MORTON_BOX);

  else
    io_num = fvm_io_num_create_from_scan(n_particles);

  global_particle_num = fvm_io_num_transfer_global_num(io_num);
  fvm_io_num_destroy(io_num);

  /* Create a new location, with ownership of global numbers */

  loc_id = cs_restart_add_location(restart,
                                   name,
                                   n_glob_particles,
                                   n_particles,
                                   global_particle_num);

  (restart->location[loc_id-1])._ent_global_num = global_particle_num;
  assert((restart->location[loc_id-1]).ent_global_num == global_particle_num);

  /* Write particle coordinates */

  CS_MALLOC(sec_name, strlen(name) + strlen(coords_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, coords_postfix);

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  cs_restart_write_section(restart,
                           sec_name,
                           loc_id,
                           3,
                           CS_TYPE_cs_real_t,
                           particle_coords);

  timing[0] = cs_timer_wtime();

  CS_FREE(sec_name);

  /* Write particle cell location information */

  CS_MALLOC(global_part_cell_num, n_particles, cs_gnum_t);

  if (restart->location[CS_MESH_LOCATION_CELLS-1].ent_global_num != nullptr) {
    const cs_gnum_t  *g_cell_num
      = restart->location[CS_MESH_LOCATION_CELLS-1].ent_global_num;
    for (i = 0; i < n_particles; i++) {
      if (particle_cell_id[i] > -1)
        global_part_cell_num[i] = g_cell_num[particle_cell_id[i]];
      else
        global_part_cell_num[i] = 0;
    }
  }
  else {
    for (i = 0; i < n_particles; i++)
      global_part_cell_num[i] = particle_cell_id[i] + 1;
  }

  CS_MALLOC(sec_name, strlen(name) + strlen(cell_num_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, cell_num_postfix);

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  cs_restart_write_section(restart,
                           sec_name,
                           loc_id,
                           1,
                           CS_TYPE_cs_gnum_t,
                           global_part_cell_num);

  CS_FREE(sec_name);

  CS_FREE(global_part_cell_num);

  return loc_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a referenced location id section from a restart file.
 *
 * The section read from file contains the global ids matching the local
 * element ids of a given location. Global id's are transformed to local
 * ids by this function.
 *
 * In case global referenced ids read do not match those of local elements,
 * id_base - 1 is assigned to the corresponding local ids.
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   location_id      id of location on which id_ref is defined
 * \param[in]   ref_location_id  id of referenced location
 * \param[in]   ref_id_base      base of location entity id numbers
 *                               (usually 0 or 1)
 * \param[out]  ref_id           array of location entity ids
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_ids(cs_restart_t     *restart,
                    const char       *sec_name,
                    int               location_id,
                    int               ref_location_id,
                    cs_lnum_t         ref_id_base,
                    cs_lnum_t        *ref_id)
{
  cs_lnum_t  n_ents = 0;
  cs_gnum_t  *g_num;

  _location_t  *ref_location = nullptr;

  int retcode = CS_RESTART_SUCCESS;

  assert(restart != nullptr);

  /* Local number of elements for location id */

  if (location_id == 0)
    n_ents = 1;

  else if (location_id > 0 && location_id <= (int)(restart->n_locations))
    n_ents = restart->location[location_id-1].n_ents;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              (int)location_id, restart->name);

  if (ref_location_id > 0 && ref_location_id <= (int)(restart->n_locations))
    ref_location = restart->location + ref_location_id-1;

  else if (ref_location_id != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              (int)location_id, restart->name);

  /* Transform local to global ids */

  CS_MALLOC(g_num, n_ents, cs_gnum_t);

  /* Read associated data */

  retcode = cs_restart_read_section(restart,
                                    sec_name,
                                    location_id,
                                    1,
                                    CS_TYPE_cs_gnum_t,
                                    g_num);

  if (retcode == CS_RESTART_SUCCESS) {

    double timing[2];

    timing[0] = cs_timer_wtime();

    if (ref_location_id == 0 || ref_location->ent_global_num == nullptr) {
      cs_lnum_t i;
      for (i = 0; i < n_ents; i++)
        ref_id[i] = g_num[i] + ref_id_base - 1;
    }

    else /* if location_id > 0 */
      cs_block_to_part_global_to_local(n_ents,
                                       ref_id_base,
                                       ref_location->n_ents,
                                       false,
                                       ref_location->ent_global_num,
                                       g_num,
                                       ref_id);

    timing[1] = cs_timer_wtime();
    _restart_wtime[restart->mode] += timing[1] - timing[0];

  }

  CS_FREE(g_num);

  /* Return */

  return retcode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write a referenced location id section to a restart file.
 *
 * The section written to file contains the global ids matching the
 * local element ids of a given location.
 *
 * \param[in]  restart          associated restart file pointer
 * \param[in]  sec_name         section name
 * \param[in]  location_id      id of location on which id_ref is defined
 * \param[in]  ref_location_id  id of referenced location
 * \param[in]  ref_id_base      base of location entity id numbers
 *                              (usually 0 or 1)
 * \param[in]  ref_id           array of location entity ids
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_write_ids(cs_restart_t           *restart,
                     const char             *sec_name,
                     int                     location_id,
                     int                     ref_location_id,
                     cs_lnum_t               ref_id_base,
                     const cs_lnum_t        *ref_id)
{
  double timing[2];
  cs_lnum_t i, n_ents = 0;
  cs_gnum_t  *g_num;

  _location_t   *ref_location = nullptr;

  assert(restart != nullptr);

  /* Local number of elements for location id */

  if (location_id == 0)
    n_ents = 1;

  else if (location_id > 0 && location_id <= (int)(restart->n_locations))
    n_ents = restart->location[location_id-1].n_ents;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              (int)location_id, restart->name);

  if (ref_location_id > 0 && ref_location_id <= (int)(restart->n_locations))
    ref_location = restart->location + ref_location_id-1;

  else if (ref_location_id != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              (int)location_id, restart->name);

  /* Transform local to global ids */

  timing[0] = cs_timer_wtime();

  CS_MALLOC(g_num, n_ents, cs_gnum_t);

  if (ref_location_id == 0) {
    for (i = 0; i < n_ents; i++)
      g_num[0] = ref_id[0] - ref_id_base + 1;
  }

  else { /* if location_id > 0 */
    if (ref_location->ent_global_num != nullptr) {
      for (i = 0; i < n_ents; i++) {
        if (ref_id[i] >= ref_id_base)
          g_num[i] = ref_location->ent_global_num[ref_id[i] - ref_id_base];
        else
          g_num[i] = 0;
      }
    }
    else {
      for (i = 0; i < n_ents; i++) {
        if (ref_id[i] >= ref_id_base)
          g_num[i] = ref_id[i] - ref_id_base + 1;
        else
          g_num[i] = 0;
      }
    }

  }

  timing[1] = cs_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  /* Write associated data */

  cs_restart_write_section(restart,
                           sec_name,
                           location_id,
                           1,
                           CS_TYPE_cs_gnum_t,
                           g_num);

  CS_FREE(g_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a section from a restart file, when that section may have used
 *         a different name in a previous version.
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   old_name         old name
 * \param[in]   location_id      id of corresponding location
 * \param[in]   n_location_vals  number of values per location (interlaced)
 * \param[in]   val_type         value type
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_section_compat(cs_restart_t           *restart,
                               const char             *sec_name,
                               const char             *old_name,
                               int                     location_id,
                               int                     n_location_vals,
                               cs_restart_val_type_t   val_type,
                               void                   *val)
{
  int retval = CS_RESTART_SUCCESS;

  assert(location_id >= 0);

  /* Check for section with current name */

  retval = cs_restart_check_section(restart,
                                    sec_name,
                                    location_id,
                                    n_location_vals,
                                    val_type);

  /* Check for older name, read and return if present */

  if (retval == CS_RESTART_ERR_N_VALS || retval == CS_RESTART_ERR_EXISTS) {

    retval =cs_restart_check_section(restart,
                                     old_name,
                                     location_id,
                                     n_location_vals,
                                     val_type);

    if (retval == CS_RESTART_SUCCESS) {
      retval = cs_restart_read_section(restart,
                                       old_name,
                                       location_id,
                                       n_location_vals,
                                       val_type,
                                       val);
      return retval;
    }

  }

  /* Read with current name (if the section is not found,
     logging will refer to the current name) */

  retval = cs_restart_read_section(restart,
                                   sec_name,
                                   location_id,
                                   n_location_vals,
                                   val_type,
                                   val);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a cs_real_3_t vector section from a restart file, when that
 *         section may have used a different name and been non-interleaved
 *         in a previous version.
 *
 * This function assumes a mesh-base location (i.e. location_id > 0)
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   old_name_x       old name, x component
 * \param[in]   old_name_y       old name, y component
 * \param[in]   old_name_z       old name, z component
 * \param[in]   location_id      id of corresponding location
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_real_3_t_compat(cs_restart_t  *restart,
                                const char    *sec_name,
                                const char    *old_name_x,
                                const char    *old_name_y,
                                const char    *old_name_z,
                                int            location_id,
                                cs_real_3_t   *val)
{
  int retval = CS_RESTART_SUCCESS;

  assert(location_id > 0);

  /* Check for section with current name */

  retval = cs_restart_check_section(restart,
                                    sec_name,
                                    location_id,
                                    3,
                                    CS_TYPE_cs_real_t);

  /* Check for older name series, read and return if present */

  if (retval == CS_RESTART_ERR_N_VALS || retval == CS_RESTART_ERR_EXISTS) {

    retval = cs_restart_check_section(restart,
                                      old_name_x,
                                      location_id,
                                      1,
                                      CS_TYPE_cs_real_t);

    if (retval == CS_RESTART_SUCCESS) {

      cs_real_t *buffer = nullptr;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      CS_MALLOC(buffer, n_ents*3, cs_real_t);

      retval = cs_restart_read_section(restart,
                                       old_name_x,
                                       location_id,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       buffer);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_y,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_z,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*2);

      if (retval == CS_RESTART_SUCCESS) {
        for (i = 0; i < n_ents; i++) {
          val[i][0] = buffer[i];
          val[i][1] = buffer[i + n_ents];
          val[i][2] = buffer[i + n_ents*2];
        }
      }

      CS_FREE(buffer);

      return retval;

    }
  }

  /* Read with current name (if the section is not found,
     logging will refer to the current name) */

  retval = cs_restart_read_section(restart,
                                   sec_name,
                                   location_id,
                                   3,
                                   CS_TYPE_cs_real_t,
                                   val);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a cs_real_6_t tensor section from a restart file, when that
 *         section may have used a different name and been non-interleaved
 *         in a previous version.
 *
 * This function assumes a mesh-base location (i.e. location_id > 0)
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   old_name_xx      old name, xx component
 * \param[in]   old_name_yy      old name, yy component
 * \param[in]   old_name_zz      old name, zz component
 * \param[in]   old_name_xy      old name, xy component
 * \param[in]   old_name_yz      old name, yz component
 * \param[in]   old_name_xz      old name, xz component
 * \param[in]   location_id      id of corresponding location
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_real_6_t_compat(cs_restart_t  *restart,
                                const char    *sec_name,
                                const char    *old_name_xx,
                                const char    *old_name_yy,
                                const char    *old_name_zz,
                                const char    *old_name_xy,
                                const char    *old_name_yz,
                                const char    *old_name_xz,
                                int            location_id,
                                cs_real_6_t   *val)
{
  int retval = CS_RESTART_SUCCESS;

  assert(location_id > 0);

  /* Check for section with current name */

  retval = cs_restart_check_section(restart,
                                    sec_name,
                                    location_id,
                                    6,
                                    CS_TYPE_cs_real_t);

  /* Check for older name series, read and return if present */

  if (retval == CS_RESTART_ERR_N_VALS || retval == CS_RESTART_ERR_EXISTS) {

    retval = cs_restart_check_section(restart,
                                      old_name_xx,
                                      location_id,
                                      1,
                                      CS_TYPE_cs_real_t);

    if (retval == CS_RESTART_SUCCESS) {

      cs_real_t *buffer = nullptr;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      CS_MALLOC(buffer, n_ents*6, cs_real_t);

      retval = cs_restart_read_section(restart,
                                       old_name_xx,
                                       location_id,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       buffer);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_yy,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_zz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*2);
      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_xy,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*3);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_yz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*4);
      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_xz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*5);

      if (retval == CS_RESTART_SUCCESS) {
        for (i = 0; i < n_ents; i++) {
          val[i][0] = buffer[i];
          val[i][1] = buffer[i + n_ents];
          val[i][2] = buffer[i + n_ents*2];
          val[i][3] = buffer[i + n_ents*3];
          val[i][4] = buffer[i + n_ents*4];
          val[i][5] = buffer[i + n_ents*5];
        }
      }

      CS_FREE(buffer);

      return retval;

    }
  }

  /* Read with current name (if the section is not found,
     logging will refer to the current name) */

  retval = cs_restart_read_section(restart,
                                   sec_name,
                                   location_id,
                                   6,
                                   CS_TYPE_cs_real_t,
                                   val);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a cs_real_66_t tensor section from a restart file, when that
 *         section may have used a different name and been non-interleaved
 *         in a previous version.
 *
 * This function assumes a mesh-base location (i.e. location_id > 0)
 *
 * \param[in]   restart          associated restart file pointer
 * \param[in]   sec_name         section name
 * \param[in]   old_name_xx      old name, xx component
 * \param[in]   old_name_yy      old name, yy component
 * \param[in]   old_name_zz      old name, zz component
 * \param[in]   old_name_xy      old name, xy component
 * \param[in]   old_name_yz      old name, yz component
 * \param[in]   old_name_xz      old name, xz component
 * \param[in]   location_id      id of corresponding location
 * \param[out]  val              array of values
 *
 * \return  0 (CS_RESTART_SUCCESS) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_read_real_66_t_compat(cs_restart_t  *restart,
                                 const char    *sec_name,
                                 const char    *old_name_xx,
                                 const char    *old_name_yy,
                                 const char    *old_name_zz,
                                 const char    *old_name_xy,
                                 const char    *old_name_yz,
                                 const char    *old_name_xz,
                                 int            location_id,
                                 cs_real_66_t  *val)
{
  int retval = CS_RESTART_SUCCESS;

  assert(location_id > 0);

  /* Check for section with current name */

  retval = cs_restart_check_section(restart,
                                    sec_name,
                                    location_id,
                                    6,
                                    CS_TYPE_cs_real_t);

  /* Check for older name series, read and return if present */

  if (retval == CS_RESTART_ERR_N_VALS || retval == CS_RESTART_ERR_EXISTS) {

    retval = cs_restart_check_section(restart,
                                      old_name_xx,
                                      location_id,
                                      1,
                                      CS_TYPE_cs_real_t);

    if (retval == CS_RESTART_SUCCESS) {

      cs_real_t *buffer = nullptr;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      CS_MALLOC(buffer, n_ents*6, cs_real_t);

      retval = cs_restart_read_section(restart,
                                       old_name_xx,
                                       location_id,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       buffer);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_yy,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_zz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*2);
      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_xy,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*3);

      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_yz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*4);
      if (retval == CS_RESTART_SUCCESS)
        retval = cs_restart_read_section(restart,
                                         old_name_xz,
                                         location_id,
                                         1,
                                         CS_TYPE_cs_real_t,
                                         buffer + n_ents*5);

      if (retval == CS_RESTART_SUCCESS) {
        for (i = 0; i < n_ents; i++) {
          val[i][0][0] = buffer[i];
          val[i][1][1] = buffer[i + n_ents*7];
          val[i][2][2] = buffer[i + n_ents*14];
          val[i][3][3] = buffer[i + n_ents*21];
          val[i][4][4] = buffer[i + n_ents*28];
          val[i][5][5] = buffer[i + n_ents*35];
        }
      }

      CS_FREE(buffer);

      return retval;

    }
  }

  /* Read with current name (if the section is not found,
     logging will refer to the current name) */
  retval = cs_restart_read_section(restart,
                                   sec_name,
                                   location_id,
                                   3,
                                   CS_TYPE_cs_real_t,
                                   val);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log checkpoint/restart setup.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_log_setup(void)
{
  char cf_str[80]; cf_str[79] = '\0';
  const char *cf_yn[] = {N_("no"), N_("yes")};

  if (_checkpoint_nt_interval <= CS_RESTART_INTERVAL_NONE)
    strncpy(cf_str, _("never"), 79);
  else if (_checkpoint_nt_interval == CS_RESTART_INTERVAL_ONLY_AT_END)
    strncpy(cf_str, _("at end of computation only"), 79);
  else if (_checkpoint_nt_interval == CS_RESTART_INTERVAL_DEFAULT) {
    strncpy(cf_str, _("default (4 checkpoints max)"), 79);
  }
  else
    snprintf(cf_str, 79, _("every %d time steps"), _checkpoint_nt_interval);

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Checkpoint / restart:\n"
                  "--------------------\n"
                  "\n"
                  "  Run is a restart:     %s\n\n"
                  "  Checkpoint frequency: %s\n"),
                cf_yn[cs_restart_present()], cf_str);

  if (_checkpoint_t_interval > 0)
    cs_log_printf(CS_LOG_SETUP,
                  _("                      : every %g s (simulated)\n"),
                  _checkpoint_t_interval);
  else if (_checkpoint_wt_interval > 0)
    cs_log_printf(CS_LOG_SETUP,
                  _("                      : every %g s (wall-clock time)\n"),
                  _checkpoint_wt_interval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print statistics associated with restart files
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_print_stats(void)
{
  bft_printf(_("\n"
               "Checkpoint / restart files summary:\n"
               "\n"
               "  Number of files read:             %3d\n"
               "  Number of files written:          %3d\n"
               "\n"
               "  Elapsed time for reading:         %12.3f\n"
               "  Elapsed time for writing:         %12.3f\n"),
             _restart_n_opens[0], _restart_n_opens[1],
             _restart_wtime[0], _restart_wtime[1]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Checks if restart is done from a NCFD checkpoint file
 *
 * \return 0 if no, 1 if yes
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_check_if_restart_from_ncfd(cs_restart_t  *r)
{
  int inttmp[1000];
  int ierror
    = cs_restart_read_section_compat(r,
                                     "neptune_cfd:checkpoint:main:version",
                                     "version_fichier_suite_principal",
                                     CS_MESH_LOCATION_NONE,
                                     1,
                                     CS_TYPE_int,
                                     inttmp);

  if (ierror == 0) {
    bft_printf(_("Remark: restarting based on a NEPTUNE_CFD computation.\n"));
    _restart_from_ncfd = 1;
  }

  return _restart_from_ncfd;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns if restart is done from a NCFD checkpoint file
 *
 * \return 0 if no, 1 if yes
 */
/*----------------------------------------------------------------------------*/

int
cs_restart_is_from_ncfd(void)
{
  return _restart_from_ncfd;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the number of checkpoints to keep.
 *
 * This function sets the maximum number of checkpoint files to keep.
 * Beyond this maximum, the oldest checkpoints are deleted to limit the
 * number of saved checkpoints. The default value is 1.
 *
 * If more than one file is kept, last one is always named "<prefix>.csc",
 * while others are names "<prefix_%04d.csc".
 *
 * %04 provides an id using 4 digits (0 padding), and value is the order
 * of writing, starting with 0. Hence: prefix_0000.csc, prefix_0001.csc, ...
 *
 * \param[in]   n_checkpoints   number of checkpoints to save.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_set_n_max_checkpoints(int  n_checkpoints)
{
  if (n_checkpoints <= 0) {

    /* deactivate checkpointing */
    _checkpoint_nt_interval = CS_RESTART_INTERVAL_NONE;
    _n_restart_directories_to_write = 0;

  }
  else
    _n_restart_directories_to_write = n_checkpoints;

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove all previous checkpoints which are not to be retained.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_clean_multiwriters_history(void)
{
  /* Check that the structure is allocated */
  if (   _restart_multiwriter == nullptr
      || _n_restart_directories_to_write < 0)
    return;

  for (int i = 0; i < _n_restart_multiwriters; i++) {
    _restart_multiwriter_t *mw = _restart_multiwriter_by_id(i);

    int n_files_to_remove
      = mw->n_prev_files - _n_restart_directories_to_write + 1;

    if (n_files_to_remove > 0) {
      for (int ii = 0; ii < n_files_to_remove; ii++) {

        if (cs_glob_rank_id <= 0) {
          char *path = mw->prev_files[ii];
          if (cs_glob_rank_id <= 0)
            cs_file_remove(path);

          /* Try to remove directory (if it is empty) */
          for (int j = strlen(path)-1; j > -1; j--) {
            if (path[j] == _dir_separator) {
              if (j > 0) {
                path[j] = '\0';
                cs_file_remove(path);
              }
              break;
            }
          }
        }

        CS_FREE(mw->prev_files[ii]);

      }

      /* Rotate available paths */

      int ii = 0;
      for (int jj = n_files_to_remove; jj < mw->n_prev_files; jj++) {
        mw->prev_files[ii] = mw->prev_files[jj];
        mw->prev_files[jj] = nullptr;
      }

      mw->n_prev_files -= n_files_to_remove;
      /* No need for extra reallocation of mw->prev_files */
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy the multiwriter structure at the end of the computation.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_multiwriters_destroy_all(void)
{
  if (_restart_multiwriter != nullptr) {
    for (int i = 0; i < _n_restart_multiwriters; i++) {
      _restart_multiwriter_t *w = _restart_multiwriter[i];

      CS_FREE(w->name);
      CS_FREE(w->path);

      for (int j = 0; j < w->n_prev_files; j++)
        CS_FREE(w->prev_files[j]);
      CS_FREE(w->prev_files);

      CS_FREE(w);

    }
    CS_FREE(_restart_multiwriter);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
