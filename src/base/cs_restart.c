/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_block_to_part.h"
#include "cs_file.h"
#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_save.h"
#include "cs_mesh_location.h"
#include "cs_part_to_block.h"
#include "cs_parall.h"
#include "cs_timer.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_restart.c
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
                                         or NULL */
  cs_gnum_t        *_ent_global_num;  /* Private global entity numbers,
                                         or NULL */

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

static int _restart_present = 0;

/* Restart time steps and frequency */

static int    _checkpoint_mesh = 1;          /* checkpoint mesh if possible */
static int    _checkpoint_nt_interval = -1;  /* time step interval */
static int    _checkpoint_nt_next = -1;      /* next forced time step */
static int    _checkpoint_nt_last = -1;      /* last checkpoint time step */
static double _checkpoint_t_interval = -1.;  /* physical time interval */
static double _checkpoint_t_next = -1.;      /* next forced time value */
static double _checkpoint_t_last = 0.;       /* last forced time value */
static double _checkpoint_wt_interval = -1.; /* wall-clock interval */
static double _checkpoint_wt_next = -1.;     /* next forced wall-clock value */
static double _checkpoint_wt_last = 0.;      /* wall-clock time of last
                                                checkpointing */

/* Restart modification */

static void                        *_restart_context = NULL;
static cs_restart_check_section_t  *_check_section_f = _check_section;
static cs_restart_read_section_t   *_read_section_f  = _read_section;
static cs_restart_write_section_t  *_write_section_f = _write_section;

static size_t        _n_locations_ref;    /* Number of locations */
static _location_t  *_location_ref;       /* Location definition array */

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

      _location_t  *loc = NULL;

      if (h.location_id != r->n_locations + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Restart file \"%s\" declares a location number %d\n"
                    "but no location %d has been declared."),
                  r->name, (int)(h.location_id),
                  (int)(r->n_locations + 1));

      BFT_REALLOC(r->location, r->n_locations + 1, _location_t);

      loc = r->location + r->n_locations;
      BFT_MALLOC(loc->name, strlen(h.sec_name) + 1, char);
      strcpy(loc->name, h.sec_name);

      loc->id = h.location_id;
      loc->n_ents = 0;
      loc->n_glob_ents = 0;

      cs_io_set_indexed_position(r->fh, &h, rec_id);
      cs_io_set_cs_gnum(&h, r->fh);
      cs_io_read_global(&h, &(loc->n_glob_ents_f), r->fh);

      loc->ent_global_num = NULL;
      loc->_ent_global_num = NULL;

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
    int                block_rank_step, min_block_size;
    MPI_Info           hints;
    MPI_Comm           block_comm, comm;

    cs_file_get_default_comm(&block_rank_step, &min_block_size,
                             &block_comm, &comm);

    r->rank_step = block_rank_step;
    r->min_block_size = min_block_size;
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
  cs_byte_t  *buffer = NULL;

  cs_lnum_t  block_buf_size = 0;

  size_t  nbr_byte_ent;

  cs_block_dist_info_t bi;

  cs_block_to_part_t *d = NULL;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_char:
    nbr_byte_ent = n_location_vals;
    break;
  case CS_TYPE_cs_int_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_int_t);
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
    assert(0);
  }

  bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                   cs_glob_n_ranks,
                                   r->rank_step,
                                   r->min_block_size / nbr_byte_ent,
                                   n_glob_ents);

  d = cs_block_to_part_create_by_gnum(cs_glob_mpi_comm,
                                      bi,
                                      n_ents,
                                      ent_global_num);

  /* Read blocks */

  block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]) * nbr_byte_ent;

  if (block_buf_size > 0)
    BFT_MALLOC(buffer, block_buf_size, cs_byte_t);

  cs_io_read_block(header,
                   bi.gnum_range[0],
                   bi.gnum_range[1],
                   buffer,
                   r->fh);

 /* Distribute blocks on ranks */

  cs_block_to_part_copy_array(d,
                              header->elt_type,
                              n_location_vals,
                              buffer,
                              vals);

  /* Free buffer */

  BFT_FREE(buffer);

  cs_block_to_part_destroy(&d);
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

  cs_datatype_t elt_type = CS_DATATYPE_NULL;
  size_t      nbr_byte_ent;
  cs_byte_t  *buffer = NULL;

  cs_block_dist_info_t bi;

  cs_part_to_block_t *d = NULL;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_char:
    nbr_byte_ent = n_location_vals;
    elt_type = CS_CHAR;
    break;
  case CS_TYPE_cs_int_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_int_t);
    elt_type = (sizeof(cs_int_t) == 8) ? CS_INT64 : CS_INT32;
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
    BFT_MALLOC(buffer, block_buf_size, cs_byte_t);

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

  BFT_FREE(buffer);

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

  if (ini_ent_num == NULL)
    return;

  switch (val_type) {

  case CS_TYPE_char:
    {
      char  *val_ord;
      char  *val_cur = (char *)vals;

      BFT_MALLOC(val_ord, n_ents * n_location_vals, char);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      BFT_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_int_t:
    {
      cs_int_t  *val_ord;
      cs_int_t  *val_cur = (cs_int_t *)vals;

      BFT_MALLOC(val_ord, n_ents * n_location_vals, cs_int_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      BFT_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_gnum_t:
    {
      cs_gnum_t  *val_ord;
      cs_gnum_t  *val_cur = (cs_gnum_t *)vals;

      BFT_MALLOC(val_ord, n_ents * n_location_vals, cs_gnum_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      BFT_FREE(val_ord);
    }
    break;

  case CS_TYPE_cs_real_t:
    {
      cs_real_t  *val_ord;
      cs_real_t  *val_cur = (cs_real_t *)vals;

      BFT_MALLOC (val_ord, n_ents * n_location_vals, cs_real_t);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[ii++]
            = val_cur[(ini_ent_num[ent_id] - 1) * n_location_vals + jj];
      }

      for (ii = 0; ii < n_ents * n_location_vals; ii++)
        val_cur[ii] = val_ord[ii];

      BFT_FREE(val_ord);
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

  if (ini_ent_num == NULL)
    return NULL;

  switch (val_type) {

  case CS_TYPE_char:
    {
      char  *val_ord;
      const char  *val_cur = (const char *)vals;

      BFT_MALLOC(val_ord, n_ents * n_location_vals, char);

      for (ent_id = 0; ent_id < n_ents; ent_id++) {
        for (jj = 0; jj < n_location_vals; jj++)
          val_ord[(ini_ent_num[ent_id] - 1) * n_location_vals + jj]
            = val_cur[ii++];
      }

      return (cs_byte_t *)val_ord;
    }
    break;

  case CS_TYPE_cs_int_t:
    {
      cs_int_t  *val_ord;
      const cs_int_t  *val_cur = (const cs_int_t *)vals;

      BFT_MALLOC(val_ord, n_ents * n_location_vals, cs_int_t);

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

      BFT_MALLOC(val_ord, n_ents * n_location_vals, cs_gnum_t);

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

      BFT_MALLOC(val_ord, n_ents * n_location_vals, cs_real_t);

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
    return NULL;

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

  char *_sec_name = NULL;
  const char *sec_name = name;

  int rec_id = -1;

  if (prefix != NULL || postfix != NULL) {


    size_t sec_name_l = strlen(name);

    if (prefix != NULL)
      sec_name_l += strlen(prefix);
    if (postfix != NULL)
      sec_name_l += strlen(postfix);

    BFT_MALLOC(_sec_name, sec_name_l + 1, char);
    sec_name = _sec_name;

    if (prefix != NULL) {
      strcpy(_sec_name, prefix);
      strcat(_sec_name, name);
    }
    else
      strcpy(_sec_name, name);

    if (postfix != NULL)
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

  BFT_FREE(_sec_name);

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

  cs_lnum_t *free_particle_ids = NULL;

  fvm_io_num_t *free_particle_io_num = NULL;
  const cs_gnum_t *free_particle_num = NULL;

  int *default_rank = NULL;

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
    return NULL;

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

  BFT_MALLOC(default_rank, _n_particles, int);
  for (i = 0; i < _n_particles; i++)
    default_rank[i] = -1;

  BFT_MALLOC(free_particle_ids, n_free_particles, cs_lnum_t);

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
  BFT_FREE(free_particle_ids);

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

  assert(restart != NULL);

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
    if (val_type != CS_TYPE_cs_int_t)
      return CS_RESTART_ERR_VAL_TYPE;
  }
  else if (header.elt_type == CS_UINT32 || header.elt_type == CS_UINT64) {
    if (val_type != CS_TYPE_cs_gnum_t && val_type != CS_TYPE_cs_int_t)
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

  cs_int_t   n_ents;
  cs_gnum_t n_glob_ents;

  const cs_gnum_t  *ent_global_num;

  size_t rec_id, rec_id_tmp;
  cs_io_sec_header_t header;

  cs_int_t _n_location_vals = n_location_vals;
  size_t index_size = 0;

  index_size = cs_io_get_index_size(restart->fh);

  assert(restart != NULL);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = NULL;
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
    if (val_type != CS_TYPE_cs_int_t) {
      bft_printf(_("  %s: section \"%s\" is not of integer type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }
  else if (header.elt_type == CS_UINT32 || header.elt_type == CS_UINT64) {
    if (val_type != CS_TYPE_cs_gnum_t && val_type != CS_TYPE_cs_int_t) {
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
    else if (val_type == CS_TYPE_cs_int_t)
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

    if (ent_global_num != NULL)
      _restart_permute_read(n_ents,
                            ent_global_num,
                            _n_location_vals,
                            val_type,
                            val);
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

  cs_int_t _n_location_vals = n_location_vals;

  assert(restart != NULL);

  n_tot_vals = _compute_n_ents(restart, location_id, n_location_vals);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = NULL;
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
  case CS_TYPE_cs_int_t:
    elt_type = (sizeof(cs_int_t) == 8) ? CS_INT64 : CS_INT32;
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

    cs_byte_t  *val_tmp = NULL;

    if (ent_global_num != NULL)
      val_tmp = _restart_permute_write(n_ents,
                                       ent_global_num,
                                       _n_location_vals,
                                       val_type,
                                       val);
    cs_io_write_global(sec_name,
                       n_tot_vals,
                       location_id,
                       0,
                       _n_location_vals,
                       elt_type,
                       (val_tmp != NULL) ? val_tmp : val,
                       restart->fh);

    if (val_tmp != NULL)
      BFT_FREE (val_tmp);
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
 * \param[in]  ent_global_num  global entity numbers, or NULL
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

    const char _checkpoint[] = "checkpoint";
    if (cs_file_mkdir_default(_checkpoint) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("The %s directory cannot be created"), _checkpoint);

    /* Move mesh_output if present */

    const char opath_i[] = "mesh_input";
    const char opath_o[] = "mesh_output";
    const char npath[] = "checkpoint/mesh_input";

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

    else if (cs_glob_mesh->modified < 1 && cs_file_isreg(opath_i)) {

#if defined(HAVE_LINKAT) && defined(HAVE_FCNTL_H)

      int retval = linkat(AT_FDCWD, opath_i,
                          AT_FDCWD, npath, AT_SYMLINK_FOLLOW);

      if (retval != 0) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(_("Failure hard-linking %s to %s:\n"
                     "%s\n"),
                   opath_i, npath, strerror(errno));

      }

#endif

    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Indicate if a restart directory is present.
 *
 * Fortran interface
 *
 * subroutine dflsui (ntsuit, ttsuit, wtsuit)
 * *****************
 *
 * integer          ntsuit      : <-- : > 0: checkpoint time step interval
 *                              :     : 0: default interval
 *                              :     : -1: checkpoint at end
 *                              :     : -2: no checkpoint
 * double precision ttsuit      : <-- : if> 0, checkpoint time interval
 * double precision wtsuit      : <-- : if> 0, checkpoint wall time interval
 *----------------------------------------------------------------------------*/

void CS_PROCF (dflsui, DFLSUI)
(
 cs_int_t   *ntsuit,
 cs_real_t  *ttsuit,
 cs_real_t  *wtsuit
)
{
  cs_restart_checkpoint_set_defaults(*ntsuit, *ttsuit, *wtsuit);
}

/*----------------------------------------------------------------------------
 * Check if checkpointing is recommended at a given time.
 *
 * Fortran interface
 *
 * subroutine reqsui (iisuit)
 * *****************
 *
 * integer          iisuit      : --> : 0 if no restart required, 1 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF (reqsui, REQSUI)
(
 cs_int_t   *iisuit
)
{
  if (cs_restart_checkpoint_required(cs_glob_time_step))
    *iisuit = 1;
  else
    *iisuit = 0;
}

/*----------------------------------------------------------------------------
 * Indicate checkpointing has been done at a given time.
 *
 * This updates the status for future checks to determine
 * if checkpointing is recommended at a given time.
 *
 * Fortran interface
 *
 * subroutine indsui
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (stusui, STUSUI)
(
 void
)
{
  cs_restart_checkpoint_done(cs_glob_time_step);
}

/*----------------------------------------------------------------------------
 * Save output mesh for turbomachinery if needed
 *
 * Fortran interface
 *
 * subroutine trbsui
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (trbsui, TRBSUI)
(
 void
)
{
  cs_mesh_save(cs_glob_mesh, NULL, "checkpoint", "mesh");
}

/*----------------------------------------------------------------------------
 * Indicate if a restart directory is present.
 *
 * Fortran interface
 *
 * subroutine indsui (isuite)
 * *****************
 *
 * integer          isuite      : --> : 1 for restart, 0 otherwise
 *----------------------------------------------------------------------------*/

void CS_PROCF (indsui, INDSUI)
(
 cs_int_t   *isuite
)
{
  *isuite = cs_restart_present();
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

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
cs_restart_checkpoint_set_defaults(int     nt_interval,
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
  assert(ts != NULL);

  int nt = ts->nt_cur - ts->nt_prev;
  double t = ts->t_cur - ts->t_prev;

  bool retval = false;

  if (_checkpoint_nt_interval > -2) {

    if (ts->nt_cur == ts->nt_max)
      retval = true;

    else if (_checkpoint_nt_interval == 0) {
      /* default interval: current number of expected time_steps for this run,
         with a minimum of 10. */
      int nt_def = (ts->nt_max - ts->nt_prev)/4;
      if (nt_def < 10)
        nt_def = 10;
      if (nt % nt_def == 0)
        retval = true;
    }

    else if (_checkpoint_nt_interval > 0 && nt % _checkpoint_nt_interval == 0)
      retval = true;

    else if (_checkpoint_nt_interval > 0 && _checkpoint_nt_last > -1) {
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
  assert(ts != NULL);

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
    if (wt - _checkpoint_wt_last >= _checkpoint_wt_interval)
      _checkpoint_wt_last = cs_timer_wtime();
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
  if (! _restart_present) {
     if (cs_file_isdir("restart"))
       _restart_present = 1;
  }

  return _restart_present;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a restart file.
 *
 * \param[in]  name  file name
 * \param[in]  path  optional directory name for output, or NULL for default
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

  char *_name = NULL;
  size_t  ldir, lname;

  const char  *_path = path;
  const char _restart[] = "restart";
  const char _checkpoint[] = "checkpoint";

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Ensure mesh checkpoint is updated on first call */

  if (    mode == CS_RESTART_MODE_WRITE
      && _restart_n_opens[mode] == 0) {
    _update_mesh_checkpoint();
  }

  /* Initializations */

  timing[0] = cs_timer_wtime();

  if (_path != NULL) {
    if (strlen(_path) == 0)
      _path = NULL;
  }

  if (_path == NULL) {
    if (mode == CS_RESTART_MODE_WRITE)
      _path = _checkpoint;
    else
      _path = _restart;
  }

  /* Create 'checkpoint' directory or read from 'restart' directory */

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

  ldir = strlen(_path);
  lname = strlen(name);

  BFT_MALLOC(_name, ldir + lname + 2, char);

  strcpy(_name, _path);
  _name[ldir] = _dir_separator;
  _name[ldir+1] = '\0';
  strcat(_name, name);
  _name[ldir+lname+1] = '\0';

  /* Allocate and initialize base structure */

  BFT_MALLOC(restart, 1, cs_restart_t);

  BFT_MALLOC(restart->name, strlen(_name) + 1, char);

  strcpy(restart->name, _name);

  BFT_FREE(_name);

  /* Initialize other fields */

  restart->mode = mode;

  restart->fh = NULL;

  restart->rank_step = 1;
  restart->min_block_size = 0;

  /* Initialize location data */

  restart->n_locations = 0;
  restart->location = NULL;

  /* Open associated file, and build an index of sections in read mode */

  _add_file(restart);

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

  assert(restart != NULL);

  mode = r->mode;

  if (r->fh != NULL)
    cs_io_finalize(&(r->fh));

  /* Free locations array */

  if (r->n_locations > 0) {
    size_t loc_id;
    for (loc_id = 0; loc_id < r->n_locations; loc_id++) {
      BFT_FREE((r->location[loc_id]).name);
      BFT_FREE((r->location[loc_id])._ent_global_num);
    }
  }
  if (r->location != NULL)
    BFT_FREE(r->location);

  /* Free remaining memory */

  BFT_FREE(r->name);

  BFT_FREE(*restart);

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

  assert(restart != NULL);

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
 * \param[in]  ent_global_num  global entity numbers, or NULL
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
        (restart->location[loc_id])._ent_global_num = NULL;

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

    BFT_REALLOC(restart->location, restart->n_locations, _location_t);
    BFT_MALLOC((restart->location[restart->n_locations-1]).name,
               strlen(location_name)+1,
               char);

    strcpy((restart->location[restart->n_locations-1]).name, location_name);

    (restart->location[restart->n_locations-1]).id = restart->n_locations;
    (restart->location[restart->n_locations-1]).n_glob_ents    = n_glob_ents;
    (restart->location[restart->n_locations-1]).n_glob_ents_f  = n_glob_ents;
    (restart->location[restart->n_locations-1]).n_ents         = n_ents;
    (restart->location[restart->n_locations-1]).ent_global_num = ent_global_num;
    (restart->location[restart->n_locations-1])._ent_global_num = NULL;

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
 * \param[in]  ent_global_num  global entity numbers, or NULL
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

  BFT_REALLOC(_location_ref, _n_locations_ref, _location_t);
  BFT_MALLOC((_location_ref[_n_locations_ref-1]).name,
             strlen(location_name)+1,
             char);

  strcpy((_location_ref[_n_locations_ref-1]).name, location_name);

  if (ent_global_num != NULL) {
    BFT_MALLOC((_location_ref[_n_locations_ref-1])._ent_global_num,
               n_ents, cs_gnum_t);
    for (cs_lnum_t i = 0; i < n_ents; i++) {
      (_location_ref[_n_locations_ref-1])._ent_global_num[i]
        = ent_global_num[i];
    }
  }
  else
    (_location_ref[_n_locations_ref-1])._ent_global_num = NULL;

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
    BFT_FREE((_location_ref[loc_id]).name);
    BFT_FREE((_location_ref[loc_id])._ent_global_num);
  }
  BFT_FREE(_location_ref);
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
 * \param[in]  context  pointer to associated data, or NULL
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
  assert(restart != NULL);

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

  assert(restart != NULL);

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

  assert(restart != NULL);

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
  assert(restart != NULL);

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

  assert(restart != NULL);

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

  assert(restart != NULL);

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
 * \param[out]  n_particles  number of particles, or NULL
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

  int  *default_p_rank = NULL;
  const cs_gnum_t  *g_cell_num
    = restart->location[CS_MESH_LOCATION_CELLS-1].ent_global_num;
  const cs_datatype_t int_type
    = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

  int loc_id = -1;

  timing[0] = cs_timer_wtime();

  if (n_particles != NULL)
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

  rec_id = _restart_section_id(restart, NULL, name, "_cell_num");

  if (rec_id < 0)
    return -1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int  *b_cell_rank, *p_cell_rank;
    cs_gnum_t  *part_cell_num = NULL;
    cs_part_to_block_t *pbd = NULL;
    cs_block_to_part_t *d = NULL;

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
      BFT_MALLOC(part_cell_num, block_buf_size, cs_gnum_t);

    header = cs_io_get_indexed_sec_header(restart->fh, rec_id);

    cs_io_set_indexed_position(restart->fh, &header, rec_id);

    cs_io_read_block(&header,
                     part_bi.gnum_range[0],
                     part_bi.gnum_range[1],
                     part_cell_num,
                     restart->fh);

    /* Build block distribution cell rank info */

    BFT_MALLOC(b_cell_rank,
               (cell_bi.gnum_range[1] - cell_bi.gnum_range[0]),
               int);

    BFT_MALLOC(p_cell_rank, n_cells, int);

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

    BFT_FREE(p_cell_rank);

    /* Now build distribution structure */

    default_p_rank = _default_p_rank(&part_bi,
                                     part_cell_num,
                                     cs_glob_mpi_comm);

    d = cs_block_to_part_create_by_adj_s(cs_glob_mpi_comm,
                                         part_bi,
                                         cell_bi,
                                         1,
                                         part_cell_num,
                                         b_cell_rank,
                                         default_p_rank);

    if (default_p_rank != NULL)
      BFT_FREE(default_p_rank);

    BFT_FREE(b_cell_rank);

    (restart->location[loc_id])._ent_global_num
      = cs_block_to_part_transfer_gnum(d);
    (restart->location[loc_id]).ent_global_num
      = (restart->location[loc_id])._ent_global_num;

    (restart->location[loc_id]).n_glob_ents = n_glob_particles;
    (restart->location[loc_id]).n_ents = cs_block_to_part_get_n_part_ents(d);

    cs_block_to_part_destroy(&d);

    BFT_FREE(part_cell_num);

  }

#endif /* #if defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {

    (restart->location[loc_id]).n_glob_ents = n_glob_particles;
    (restart->location[loc_id]).n_ents = n_glob_particles;

  }

  if (n_particles != NULL)
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
  char *sec_name = NULL;

  cs_lnum_t n_cells = (restart->location[CS_MESH_LOCATION_CELLS-1]).n_ents;
  const cs_gnum_t *g_cells_num
    = (restart->location[CS_MESH_LOCATION_CELLS-1]).ent_global_num;

  const char *name = (restart->location[particles_location_id - 1]).name;
  const char *cell_num_postfix = "_cell_num";
  const char *coords_postfix = "_coords";

  int retcode = CS_RESTART_SUCCESS;

  cs_lnum_t  n_particles = (restart->location[particles_location_id - 1]).n_ents;

  /* Read particle coordinates */

  BFT_MALLOC(sec_name, strlen(name) + strlen(coords_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, coords_postfix);

  retcode = cs_restart_read_section(restart,
                                    sec_name,
                                    particles_location_id,
                                    3,
                                    CS_TYPE_cs_real_t,
                                    particle_coords);

  BFT_FREE(sec_name);

  if (retcode != CS_RESTART_SUCCESS)
    return retcode;

  /* Read particle cell id */

  BFT_MALLOC(sec_name, strlen(name) + strlen(cell_num_postfix) + 1, char);
  strcpy(sec_name, name);
  strcat(sec_name, cell_num_postfix);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t *g_part_cell_num;

    BFT_MALLOC(g_part_cell_num, n_particles, cs_gnum_t);

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

    BFT_FREE(g_part_cell_num);

    timing[1] = cs_timer_wtime();
    _restart_wtime[restart->mode] += timing[1] - timing[0];

  }

#endif /* #if defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      particles_location_id,
                                      1,
                                      CS_TYPE_cs_int_t,
                                      particle_cell_id);
    for (cs_lnum_t i = 0; i < n_particles; i++)
      particle_cell_id[i] -= 1;
  }

  BFT_FREE(sec_name);

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
  cs_gnum_t  *global_particle_num = NULL;
  cs_gnum_t  *global_part_cell_num = NULL;
  fvm_io_num_t  *io_num = NULL;
  char *sec_name = NULL;

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

  BFT_MALLOC(sec_name, strlen(name) + strlen(coords_postfix) + 1, char);
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

  BFT_FREE(sec_name);

  /* Write particle cell location information */

  BFT_MALLOC(global_part_cell_num, n_particles, cs_gnum_t);

  if (restart->location[CS_MESH_LOCATION_CELLS-1].ent_global_num != NULL) {
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

  BFT_MALLOC(sec_name, strlen(name) + strlen(cell_num_postfix) + 1, char);
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

  BFT_FREE(sec_name);

  BFT_FREE(global_part_cell_num);

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

  _location_t  *ref_location = NULL;

  int retcode = CS_RESTART_SUCCESS;

  assert(restart != NULL);

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

  BFT_MALLOC(g_num, n_ents, cs_gnum_t);

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

    if (ref_location_id == 0 || ref_location->ent_global_num == NULL) {
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

  BFT_FREE(g_num);

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

  _location_t   *ref_location = NULL;

  assert(restart != NULL);

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

  BFT_MALLOC(g_num, n_ents, cs_gnum_t);

  if (ref_location_id == 0) {
    for (i = 0; i < n_ents; i++)
      g_num[0] = ref_id[0] - ref_id_base + 1;
  }

  else { /* if location_id > 0 */
    if (ref_location->ent_global_num != NULL) {
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

  BFT_FREE(g_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read a section from a restart file, when that section may have used
 *         a different name in a previous version.
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
 * \param[in]   old_name_y       old name, z component
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

      cs_real_t *buffer = NULL;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      BFT_MALLOC(buffer, n_ents*3, cs_real_t);

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

      BFT_FREE(buffer);

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

      cs_real_t *buffer = NULL;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      BFT_MALLOC(buffer, n_ents*6, cs_real_t);

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

      BFT_FREE(buffer);

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

      cs_real_t *buffer = NULL;
      cs_lnum_t i;
      cs_lnum_t n_ents = (restart->location[location_id-1]).n_ents;

      BFT_MALLOC(buffer, n_ents*6, cs_real_t);

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

      BFT_FREE(buffer);

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

END_C_DECLS
