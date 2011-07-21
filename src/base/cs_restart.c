/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Manage checkpoint / restart files
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_file.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_block_to_part.h>
#include <fvm_part_to_block.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_io.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* Fortran API */
/* ----------- */

/*
 * "Usual" max name length (a longer name is possible but will
 * incur dynamic memory allocation.
 */

#define CS_RESTART_NAME_LEN   64

/*============================================================================
 * Local type definitions
 *============================================================================*/

typedef struct _location_t {

  char              *name;            /* Location name */
  size_t             id;              /* Associated id in file */
  fvm_lnum_t         n_ents;          /* Local number of entities */
  fvm_gnum_t         n_glob_ents_f;   /* Global number of entities by file */
  fvm_gnum_t         n_glob_ents;     /* Global number of entities */
  const fvm_gnum_t  *ent_global_num;  /* Global entity numbers, or NULL */

} _location_t;

struct _cs_restart_t {

  char              *name;         /* Name of restart file */

  cs_io_t           *fh;           /* Pointer to associated file handle */

  size_t             n_locations;  /* Number of locations */
  _location_t       *location;     /* Location definition array */

  cs_restart_mode_t  mode;         /* Read or write */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

#if defined(WIN32) || defined(_WIN32)
static const char _dir_separator = '\\';
#else
static const char _dir_separator = '/';
#endif

/* Minimum buffer size on rank 0 (to limit number of blocks
   when there is a large number of processors) */

static int cs_restart_def_buf_size = 1024*1024*8;

/* Monitoring info */

static int    _restart_n_opens[2] = {0, 0};
static double _restart_wtime[2] = {0.0, 0.0};

/* Array for Fortran API */

static size_t         _restart_pointer_size = 2;
static cs_restart_t  *_restart_pointer_base[2] = {NULL, NULL};
static cs_restart_t **_restart_pointer = _restart_pointer_base;

/* Do we have a restart directory ? */

static int _restart_present = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return an available id in the array of restart file pointers.
 *
 * The array may be allocated or reallocated if necessary.
 *
 * returns:
 *   available id in the array of restart file pointers
 *----------------------------------------------------------------------------*/

static int
_new_restart_id(void)
{
  size_t i, j;

  for (i = 0;
       i < _restart_pointer_size && _restart_pointer[i] != NULL;
       i++);

  /* If no slot is available, we allow for more restart files */

  if (i == _restart_pointer_size) {

    if (_restart_pointer == _restart_pointer_base) {
      BFT_MALLOC(_restart_pointer, _restart_pointer_size*2, cs_restart_t *);
      for (j = 0; j < _restart_pointer_size; j++) {
        _restart_pointer[j] = _restart_pointer_base[j];
        _restart_pointer_base[j] = NULL;
      }
    }
    else
      BFT_REALLOC(_restart_pointer, _restart_pointer_size*2, cs_restart_t *);

    for (j = _restart_pointer_size; j < _restart_pointer_size * 2; j++)
      _restart_pointer[j] = NULL;

    _restart_pointer_size *= 2;

  }

  return i;
}

/*----------------------------------------------------------------------------
 * Free a slot from the array of restart file pointers.
 *
 * The array may be freed or reallocated if possible.
 *
 * parameters:
 *   id <-- id to free in the array of restart file pointers
 *----------------------------------------------------------------------------*/

static void
_free_restart_id(int id)
{
  size_t i, j;

  const size_t restart_pointer_base_size = 2;

  _restart_pointer[id] = NULL;

  /* Revert from dynamic to static array if applicable and possible */

  if ((size_t)id >= restart_pointer_base_size) {

    for (i = restart_pointer_base_size;
         i < _restart_pointer_size && _restart_pointer[i] == NULL;
         i++);

    /* If no slot above static array size is used, revert to static array  */

    if (i == _restart_pointer_size) {

      for (j = 0; j < restart_pointer_base_size; j++)
        _restart_pointer_base[j] = _restart_pointer[j];

      _restart_pointer_size = restart_pointer_base_size;

      BFT_FREE(_restart_pointer[j]);

      _restart_pointer = _restart_pointer_base;
    }

  }
}

/*----------------------------------------------------------------------------
 * Compute number of values in a record
 *
 * parameters:
 *   r               <-- associated restart file pointer
 *   location_id     <-- location id
 *   n_location_vals <-- number of values per location
 *----------------------------------------------------------------------------*/

static size_t
_compute_n_ents(const cs_restart_t  *r,
                size_t               location_id,
                size_t               n_location_vals)
{
  size_t retval = 0;

  if (location_id == 0)
    retval = n_location_vals;

  else if (location_id > 0 && location_id <= r->n_locations)
    retval = r->location[location_id-1].n_glob_ents_f * n_location_vals;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              location_id, r->name);

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
      cs_io_set_fvm_gnum(&h, r->fh);
      cs_io_read_global(&h, &(loc->n_glob_ents_f), r->fh);

      loc->ent_global_num = NULL;

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

  const char magic_string[] = "Checkpoint / restart, R0";
  const long echo = CS_IO_ECHO_NONE;

  timing[0] = bft_timer_wtime();

  /* In read mode, open file to detect header first */

  if (r->mode == CS_RESTART_MODE_READ) {

#if defined(HAVE_MPI)
    r->fh = cs_io_initialize_with_index(r->name,
                                        magic_string,
                                        cs_glob_io_hints,
                                        echo,
                                        cs_glob_mpi_comm);
#else
    r->fh = cs_io_initialize_with_index(r->name, magic_string, 0, echo);
#endif

    _locations_from_index(r);
  }

  else {

#if defined(HAVE_MPI)
    r->fh = cs_io_initialize(r->name,
                             magic_string,
                             CS_IO_MODE_WRITE,
                             cs_glob_io_hints,
                             echo,
                             cs_glob_mpi_comm);
#else
    r->fh = cs_io_initialize(r->name,
                             magic_string,
                             CS_IO_MODE_WRITE,
                             0,
                             echo);
#endif
  }

  timing[1] = bft_timer_wtime();
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
_read_ent_values(cs_restart_t        *r,
                 cs_io_sec_header_t  *header,
                 fvm_gnum_t           n_glob_ents,
                 fvm_lnum_t           n_ents,
                 const fvm_gnum_t     ent_global_num[],
                 int                  n_location_vals,
                 cs_type_t            val_type,
                 cs_byte_t            vals[])
{
  cs_byte_t  *buffer = NULL;

  fvm_lnum_t  block_buf_size = 0;

  size_t  nbr_byte_ent, nbr_byte_val;

  fvm_block_to_part_info_t bi;

  fvm_block_to_part_t *d = NULL;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_cs_int_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_int_t);
    nbr_byte_val = sizeof(cs_int_t);
    cs_io_set_fvm_lnum(header, r->fh);
    break;
  case CS_TYPE_cs_real_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_real_t);
    nbr_byte_val = sizeof(cs_real_t);
    break;
  default:
    assert(val_type == CS_TYPE_cs_int_t || val_type == CS_TYPE_cs_real_t);
  }

  bi = fvm_block_to_part_compute_sizes(cs_glob_rank_id,
                                       cs_glob_n_ranks,
                                       0,
                                       cs_restart_def_buf_size / nbr_byte_ent,
                                       n_glob_ents);

  d = fvm_block_to_part_create_by_gnum(cs_glob_mpi_comm,
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

  fvm_block_to_part_copy_array(d,
                               header->elt_type,
                               n_location_vals,
                               buffer,
                               vals);

  /* Free buffer */

  BFT_FREE(buffer);

  fvm_block_to_part_destroy(&d);
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
_write_ent_values(const cs_restart_t  *r,
                  const char          *sec_name,
                  fvm_gnum_t           n_glob_ents,
                  fvm_lnum_t           n_ents,
                  const fvm_gnum_t    *ent_global_num,
                  int                  location_id,
                  int                  n_location_vals,
                  cs_type_t            val_type,
                  const cs_byte_t     *vals)
{
  fvm_lnum_t  block_buf_size = 0;

  fvm_datatype_t elt_type = FVM_DATATYPE_NULL;
  size_t      nbr_byte_ent;
  cs_byte_t  *buffer = NULL;

  fvm_part_to_block_info_t bi;

  fvm_part_to_block_t *d = NULL;

  /* Initialization */

  switch (val_type) {
  case CS_TYPE_cs_int_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_int_t);
    elt_type = (sizeof(cs_int_t) == 8) ? FVM_INT64 : FVM_INT32;
    break;
  case CS_TYPE_cs_real_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_real_t);
    elt_type =   (sizeof(cs_real_t) == fvm_datatype_size[FVM_DOUBLE])
               ? FVM_DOUBLE : FVM_FLOAT;
    break;
  default:
    assert(val_type == CS_TYPE_cs_int_t || val_type == CS_TYPE_cs_real_t);
  }

  bi = fvm_part_to_block_compute_sizes(cs_glob_rank_id,
                                       cs_glob_n_ranks,
                                       0,
                                       cs_restart_def_buf_size / nbr_byte_ent,
                                       n_glob_ents);

  d = fvm_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                       bi,
                                       n_ents,
                                       ent_global_num);

  /* Distribute to blocks */

  block_buf_size = (bi.gnum_range[1] - bi.gnum_range[0]) * nbr_byte_ent;

  if (block_buf_size > 0)
    BFT_MALLOC(buffer, block_buf_size, cs_byte_t);

  /* Distribute blocks on ranks */

  fvm_part_to_block_copy_array(d,
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

  fvm_part_to_block_destroy(&d);
}

#endif /* #if defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Convert read/write arguments from the Fortran API to the C API.
 *
 * parameters:
 *   numsui   <-- restart file id
 *   itysup   <-- location type code
 *   irtype   <-- integer or real
 *   r        <-- pointer to restart file handle
 *   location  <-- location id
 *   val_type <-- integer of real
 *   ierror   <-- 0 = success, < 0 = error
 *----------------------------------------------------------------------------*/

static void
_section_f77_to_c(const cs_int_t   *numsui,
                  const cs_int_t   *itysup,
                  const cs_int_t   *irtype,
                  cs_restart_t    **r,
                  int              *location,
                  cs_type_t        *val_type,
                  cs_int_t         *ierror)
{
  cs_int_t r_id = *numsui - 1;

  *ierror = CS_RESTART_SUCCES;

  /* Pointer to associated restart file handle */

  if (   r_id < 0
      || r_id > (cs_int_t)_restart_pointer_size
      || _restart_pointer[r_id] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Restart file number <%d> can not be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_RESTART_ERR_FILE_NUM;
    return;
  }

  else
    *r = _restart_pointer[r_id];

  /* Location associated with section */

  switch (*itysup) {

  case 0:
    *location = CS_RESTART_LOCATION_NONE;
    break;

  case 1:
    *location = CS_RESTART_LOCATION_CELL;
    break;

  case 2:
    *location = CS_RESTART_LOCATION_I_FACE;
    break;

  case 3:
    *location = CS_RESTART_LOCATION_B_FACE;
    break;

  case 4:
    *location = CS_RESTART_LOCATION_VERTEX;
    break;

  default:
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Location type <%d> given for a restart file section\n"
                 "is invalid using the Fortran API."), (int)(*itysup));
    *ierror = CS_RESTART_ERR_LOCATION;
    return;

  }

  /* Val_Type associated with section */

  switch (*irtype) {

  case 1:
    *val_type = CS_TYPE_cs_int_t;
    break;

  case 2:
    *val_type = CS_TYPE_cs_real_t;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Value type <%d> given for a restart file section\n"
                "is invalid using the Fortran API."), (int)(*irtype));
    *ierror = CS_RESTART_ERR_VAL_TYPE;
    return;

  }
}

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
_restart_permute_read(cs_int_t           n_ents,
                      const fvm_gnum_t  *ini_ent_num,
                      cs_int_t           n_location_vals,
                      cs_type_t          val_type,
                      cs_byte_t         *vals)
{
  cs_int_t ent_id, jj;

  cs_int_t ii = 0;

  /* Instructions */

  if (ini_ent_num == NULL)
    return;

  switch (val_type) {

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
    assert(val_type == CS_TYPE_cs_int_t || val_type == CS_TYPE_cs_real_t);

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
_restart_permute_write(cs_int_t           n_ents,
                       const fvm_gnum_t  *ini_ent_num,
                       cs_int_t           n_location_vals,
                       cs_type_t          val_type,
                       const cs_byte_t   *vals)
{
  cs_int_t  ent_id, jj;

  cs_int_t  ii = 0;

  /* Instructions */

  if (ini_ent_num == NULL)
    return NULL;

  switch (val_type) {

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
    assert(val_type == CS_TYPE_cs_int_t || val_type == CS_TYPE_cs_real_t);
    return NULL;

  }
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Open a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE OPNSUI (NOMSUI, LNGNOM, IREAWR, NUMSUI, IERROR)
 * *****************
 *
 * CHARACTER*       NOMSUI      : --> : Restart file name
 * INTEGER          LNGNOM      : --> : Restart file name length
 * INTEGER          IREAWR      : --> : 1: read; 2: write
 * INTEGER          NUMSUI      : <-- : Number of opened restart file
 * INTEGER          IERROR      : <-- : 0: success; < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (opnsui, OPNSUI)
(
 const char       *nomsui,
 const cs_int_t   *lngnom,
 const cs_int_t   *ireawr,
       cs_int_t   *numsui,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  char    *bufname;

  size_t   id;

  cs_restart_mode_t restart_mode;


  /* Initialization */

  *numsui = 0;
  *ierror = CS_RESTART_SUCCES;

  /* Handle name for C API */

  bufname = cs_base_string_f_to_c_create(nomsui, *lngnom);

  /* File creation options */

  {
    switch(*ireawr) {
    case 1:
      restart_mode = CS_RESTART_MODE_READ;
      break;
    case 2:
      restart_mode = CS_RESTART_MODE_WRITE;
      break;
    default:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The access mode of the restart file <%s>\n"
                   "must be equal to 1 (read) or 2 (write) and not <%d>."),
                 bufname, (int)(*ireawr));

      *ierror = CS_RESTART_ERR_MODE;
    }

  }

  /* Search for an available slot and create file */

  if (*ierror == CS_RESTART_SUCCES) {

    id = _new_restart_id();
    _restart_pointer[id] = cs_restart_create(bufname, restart_mode);

    /* Return the position of the handle in the array
     * (id + 1 to have a 1 to n numbering, more conventional in Fortran) */

    *numsui = id + 1;
  }
  else
    *numsui = -1;

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&bufname);
}

/*----------------------------------------------------------------------------
 * Close a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE CLSSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : <-> : number of restart file to close
 * INTEGER          IERROR      : <-- : 0: success; < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (clssui, CLSSUI)
(
 const cs_int_t   *numsui,
       cs_int_t   *ierror
)
{
  cs_int_t r_id = *numsui - 1;

  *ierror = CS_RESTART_SUCCES;

  /* Check that the file is valid */

  if (   r_id < 0
      || r_id > (cs_int_t)_restart_pointer_size
      || _restart_pointer[r_id] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Restart file number <%d> can not be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_RESTART_ERR_FILE_NUM;
    return;
  }

  /* Close file */

  cs_restart_destroy(_restart_pointer[r_id]);

  _free_restart_id(r_id);
}


/*----------------------------------------------------------------------------
 * Check the locations associated with a restart file.
 *
 * For each type of entity, return 1 if the associated number of entities
 * matches the current value (and so that we consider the mesh locations are
 * the same), 0 otherwise.
 *
 * Fortran interface
 *
 * SUBROUTINE TSTSUI (NUMSUI, INDCEL, INDFAC, INDFBR, INDSOM)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * INTEGER          INDCEL      : <-- : Matching cells flag
 * INTEGER          INDFAC      : <-- : Matching interior faces flag
 * INTEGER          INDFBR      : <-- : Matching boundary faces flag
 * INTEGER          INDSOM      : <-- : Matching vertices flag
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstsui, TSTSUI)
(
 const cs_int_t  *numsui,
       cs_int_t  *indcel,
       cs_int_t  *indfac,
       cs_int_t  *indfbr,
       cs_int_t  *indsom
)
{
  cs_bool_t  match_cell, match_i_face, match_b_face, match_vertex;

  cs_int_t   r_id   = *numsui - 1;

  /* Associated structure pointer */

  if (   r_id < 0
      || r_id > (cs_int_t)_restart_pointer_size
      || _restart_pointer[r_id] == NULL) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Information on the restart file number <%d> unavailable\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *indcel = 0;
    *indfac = 0;
    *indfbr = 0;
    *indsom = 0;
    return;
  }

  else {

    cs_restart_check_base_location(_restart_pointer[r_id],
                                   &match_cell, &match_i_face,
                                   &match_b_face, &match_vertex);

    *indcel = (match_cell == true ? 1 : 0);
    *indfac = (match_i_face == true ? 1 : 0);
    *indfbr = (match_b_face == true ? 1 : 0);
    *indsom = (match_vertex == true ? 1 : 0);

  }

}


/*----------------------------------------------------------------------------
 * Print index associated with a restart file in read mode
 *
 * Fortran interface
 *
 * SUBROUTINE INFSUI (NUMSUI)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 *----------------------------------------------------------------------------*/

void CS_PROCF (infsui, INFSUI)
(
 const cs_int_t  *numsui
)
{
  cs_int_t   r_id   = *numsui - 1;

  /* Associated structure pointer */

  if (   r_id < 0
      || r_id > (cs_int_t)_restart_pointer_size
      || _restart_pointer[r_id] == NULL) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Information on the restart file number <%d> unavailable\n"
                 "(file already closed or invalid number)."), (int)(*numsui));
  }
  else {

    cs_restart_dump_index(_restart_pointer[r_id]);

  }
}


/*----------------------------------------------------------------------------
 * Read a section from a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE LECSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * CHARACTER*       NOMRUB      : --> : Section name
 * INTEGER          LNGNOM      : --> : Section name length
 * INTEGER          ITYSUP      : --> : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * INTEGER          NBVENT      : --> : N. values per location entity
 * INTEGER          IRTYPE      : --> : 1 for integers, 2 for double precision
 * (?)              TABVAR      : <-> : Array of values to read
 * INTEGER          IERROR      : <-- : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecsui, LECSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *nbvent,
 const cs_int_t   *irtype,
       void       *tabvar,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  char    *bufname;

  cs_type_t   val_type;

  cs_restart_t  *restart;
  int          location_id;


  *ierror = CS_RESTART_SUCCES;

  /* Handle name for C API */

  bufname = cs_base_string_f_to_c_create(nomrub, *lngnom);

  /* Handle other arguments for C API */

  _section_f77_to_c(numsui,
                    itysup,
                    irtype,
                    &restart,
                    &location_id,
                    &val_type,
                    ierror);

  if (*ierror < CS_RESTART_SUCCES)
    return;

  /* Read section */

  *ierror = cs_restart_read_section(restart,
                                    bufname,
                                    location_id,
                                    *nbvent,
                                    val_type,
                                    tabvar);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&bufname);
}


/*----------------------------------------------------------------------------
 * Write a section to a restart file
 *
 * Fortran interface
 *
 * SUBROUTINE ECRSUI (NUMSUI, NOMRUB, LNGNOM, ITYSUP, NBVENT, IRTYPE, TABVAR)
 * *****************
 *
 * INTEGER          NUMSUI      : --> : Restart file number
 * CHARACTER*       NOMRUB      : --> : Section name
 * INTEGER          LNGNOM      : --> : Section name length
 * INTEGER          ITYSUP      : --> : Location type:
 *                              :     :  0: scalar (no location)
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: boundary faces
 *                              :     :  4: vertices (if available)
 * INTEGER          NBVENT      : --> : N. values per location entity
 * INTEGER          IRTYPE      : --> : 1 for integers, 2 for double precision
 * (?)              TABVAR      : --> : Array of values to write
 * INTEGER          IERROR      : <-- : 0: success, < 0: error code
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrsui, ECRSUI)
(
 const cs_int_t   *numsui,
 const char       *nomrub,
 const cs_int_t   *lngnom,
 const cs_int_t   *itysup,
 const cs_int_t   *nbvent,
 const cs_int_t   *irtype,
 const void       *tabvar,
       cs_int_t   *ierror
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  char *bufname;

  cs_type_t val_type;

  cs_restart_t *restart;
  int location_id;


  *ierror = CS_RESTART_SUCCES;

  /* Handle name for C API */

  bufname = cs_base_string_f_to_c_create(nomrub, *lngnom);

  /* Handle other arguments for C API */

  _section_f77_to_c(numsui,
                    itysup,
                    irtype,
                    &restart,
                    &location_id,
                    &val_type,
                    ierror);

  if (*ierror < CS_RESTART_SUCCES)
    return;

  /* Write section */

  cs_restart_write_section(restart,
                           bufname,
                           location_id,
                           *nbvent,
                           val_type,
                           tabvar);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&bufname);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if we have a restart directory.
 *
 * returns:
 *   1 if a restart directory is present, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_restart_present(void)
{
  if (! _restart_present) {
     if (bft_file_isdir("restart"))
       _restart_present = 1;
  }

  return _restart_present;
}

/*----------------------------------------------------------------------------
 * Initialize a restart file
 *
 * parameters:
 *   name <-- file name
 *   mode <-- read or write
 *
 * returns:
 *   pointer to initialized restart file structure
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_create(const char         *name,
                  cs_restart_mode_t   mode)
{
  cs_restart_t  * restart;

  char  *path = NULL;
  double timing[2];

  const cs_mesh_t  *mesh = cs_glob_mesh;

  timing[0] = bft_timer_wtime();

  /* Create 'checkpoint' directory or read from 'restart' directory */

  if (mode == CS_RESTART_MODE_WRITE) {

    const char  dir[] = "checkpoint";
    size_t  ldir, lname;

    ldir = strlen(dir);
    lname = strlen(name);

    if (bft_file_mkdir_default(dir) == 0) {
      BFT_MALLOC(path, ldir + lname + 2, char);
      strcpy(path, dir);
      path[ldir] = _dir_separator;
      path[ldir+1] = '\0';
      strcat(path, name);
      path[ldir+lname+1] = '\0';
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("The checkpoint directory cannot be created"));

  }
  else if (mode == CS_RESTART_MODE_READ) {

    const char dir[] = "restart";
    size_t  ldir, lname;

    ldir = strlen(dir);
    lname = strlen(name);

    if (bft_file_isdir(dir) == 1) {
      BFT_MALLOC(path, ldir + lname + 2, char);
      strcpy(path, dir);
      path[ldir] = _dir_separator;
      path[ldir+1] = '\0';
      strcat(path, name);
      path[ldir+lname+1] = '\0';
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                _("The restart directory cannot be found"));

  }

  /* Allocate and initialize base structure */

  BFT_MALLOC(restart, 1, cs_restart_t);

  BFT_MALLOC(restart->name, strlen(path) + 1, char);

  strcpy(restart->name, path);

  BFT_FREE(path);

  /* Initialize other fields */

  restart->mode = mode;

  restart->fh = NULL;

  /* Initialize location data */

  restart->n_locations = 0;
  restart->location = NULL;

  /* Open associated file, and build an index of sections in read mode */

  _add_file(restart);

  /* Add basic location definitions */

  cs_restart_add_location(restart, "cells",
                          mesh->n_g_cells, mesh->n_cells,
                          mesh->global_cell_num);
  cs_restart_add_location(restart, "interior_faces",
                          mesh->n_g_i_faces, mesh->n_i_faces,
                          mesh->global_i_face_num);
  cs_restart_add_location(restart, "boundary_faces",
                          mesh->n_g_b_faces, mesh->n_b_faces,
                          mesh->global_b_face_num);
  cs_restart_add_location(restart, "vertices",
                          mesh->n_g_vertices, mesh->n_vertices,
                          mesh->global_vtx_num);

  timing[1] = bft_timer_wtime();
  _restart_wtime[mode] += timing[1] - timing[0];

  return restart;
}

/*----------------------------------------------------------------------------
 * Destroy structure associated with a restart file (and close the file).
 *
 * parameters:
 *   restart <-- pointer to restart file structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_restart_t *
cs_restart_destroy(cs_restart_t  *restart)
{
  cs_restart_mode_t   mode;

  double timing[2];

  timing[0] = bft_timer_wtime();

  assert(restart != NULL);

  mode = restart->mode;

  if (restart->fh != NULL)
    cs_io_finalize(&(restart->fh));

  /* Free locations array */

  if (restart->n_locations > 0) {
    size_t loc_id;
    for (loc_id = 0; loc_id < restart->n_locations; loc_id++)
      BFT_FREE((restart->location[loc_id]).name);
  }
  if (restart->location != NULL)
    BFT_FREE(restart->location);

  /* Free remaining memory */

  BFT_FREE(restart->name);
  BFT_FREE(restart);

  timing[1] = bft_timer_wtime();
  _restart_wtime[mode] += timing[1] - timing[0];

  return NULL;
}

/*----------------------------------------------------------------------------
 * Check the locations associated with a restart file.
 *
 * For each type of entity, the correspondinf flag is set to true if the
 * associated number of entities matches the current value (and so that we
 * consider the mesh locations are the same), false otherwise.
 *
 * parameters:
 *   restart      <-- associated restart file pointer
 *   match_cell   <-- matching cells flag
 *   match_i_face <-- matching interior faces flag
 *   match_b_face <-- matching boundary faces flag
 *   match_vertex <-- matching vertices flag
 *----------------------------------------------------------------------------*/

void
cs_restart_check_base_location(const cs_restart_t  *restart,
                               cs_bool_t           *match_cell,
                               cs_bool_t           *match_i_face,
                               cs_bool_t           *match_b_face,
                               cs_bool_t           *match_vertex)
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

/*----------------------------------------------------------------------------
 * Add a location definition.
 *
 * parameters:
 *   restart        <-- associated restart file pointer
 *   location_name  <-- name associated with the location
 *   n_glob_ents    <-- global number of entities
 *   n_ents         <-- local number of entities
 *   ent_global_num <-- global entity numbers, or NULL
 *
 * returns:
 *   the location id assigned, or -1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_add_location(cs_restart_t      *restart,
                        const char        *location_name,
                        fvm_gnum_t         n_glob_ents,
                        fvm_lnum_t         n_ents,
                        const fvm_gnum_t  *ent_global_num)
{
  double timing[2];

  int loc_id;

  timing[0] = bft_timer_wtime();

  if (restart->mode == CS_RESTART_MODE_READ) {

    /* Search for a location with the same name */

    for (loc_id = 0; loc_id < (int)(restart->n_locations); loc_id++) {

      if ((strcmp((restart->location[loc_id]).name, location_name) == 0)) {

        (restart->location[loc_id]).n_glob_ents = n_glob_ents;

        (restart->location[loc_id]).n_ents  = n_ents;
        (restart->location[loc_id]).ent_global_num = ent_global_num;

        timing[1] = bft_timer_wtime();
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

    fvm_datatype_t gnum_type
      = (sizeof(fvm_gnum_t) == 8) ? FVM_UINT64 : FVM_UINT32;

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

    cs_io_write_global(location_name, 1, restart->n_locations, 0, 0,
                       gnum_type, &n_glob_ents,
                       restart->fh);

    timing[1] = bft_timer_wtime();
    _restart_wtime[restart->mode] += timing[1] - timing[0];

    return restart->n_locations;
  }

  timing[1] = bft_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  return -1;
}

/*----------------------------------------------------------------------------
 * Print the index associated with a restart file in read mode
 *
 * parameters:
 *   restart <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

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


/*----------------------------------------------------------------------------
 * Read a section from a restart file.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *   val             --> array of values
 *
 * returns: 0 (CS_RESTART_SUCCES) in case of success,
 *          or error code (CS_RESTART_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_restart_read_section(cs_restart_t  *restart,
                        const char    *sec_name,
                        int            location_id,
                        cs_int_t       n_location_vals,
                        cs_type_t      val_type,
                        void          *val)
{
  double timing[2];

  cs_int_t   n_ents;
  fvm_gnum_t n_glob_ents;

  const fvm_gnum_t  *ent_global_num;

  size_t rec_id, rec_id_tmp;
  cs_io_sec_header_t header;

  cs_int_t _n_location_vals = n_location_vals;
  size_t index_size = 0;

  timing[0] = bft_timer_wtime();

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
                 restart->name, sec_name, header.location_id, location_id);
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

  if (header.elt_type == FVM_INT32 || header.elt_type == FVM_INT64) {
    cs_io_set_fvm_lnum(&header, restart->fh);
    if (val_type != CS_TYPE_cs_int_t) {
      bft_printf(_("  %s: section \"%s\" is not of integer type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }
  else if (header.elt_type == FVM_FLOAT || header.elt_type == FVM_DOUBLE) {
    if (sizeof(cs_real_t) != fvm_datatype_size[header.elt_type]) {
      if (sizeof(cs_real_t) == fvm_datatype_size[FVM_FLOAT])
        header.elt_type = FVM_FLOAT;
      else
        header.elt_type = FVM_DOUBLE;
    }
    if (val_type != CS_TYPE_cs_real_t) {
      bft_printf(_("  %s: section \"%s\" is not of floating-point type.\n"),
                 restart->name, sec_name);
      return CS_RESTART_ERR_VAL_TYPE;
    }
  }

  /* Now set position in file to read data */

  cs_io_set_indexed_position(restart->fh, &header, rec_id);

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

  else
    _read_ent_values(restart,
                     &header,
                     n_glob_ents,
                     n_ents,
                     ent_global_num,
                     _n_location_vals,
                     val_type,
                     (cs_byte_t *)val);

#endif /* #if defined(HAVE_MPI) */

  timing[1] = bft_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

  /* Return */

  return CS_RESTART_SUCCES;
}

/*----------------------------------------------------------------------------
 * Write a section to a restart file.
 *
 * parameters:
 *   restart         <-- associated restart file pointer
 *   sec_name        <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   val_type        <-- value type
 *   val             <-- array of values
 *----------------------------------------------------------------------------*/

void
cs_restart_write_section(cs_restart_t  *restart,
                         const char    *sec_name,
                         int            location_id,
                         cs_int_t       n_location_vals,
                         cs_type_t      val_type,
                         const void    *val)
{
  double timing[2];

  cs_int_t         n_tot_vals, n_glob_ents, n_ents;
  fvm_datatype_t   elt_type = FVM_DATATYPE_NULL;

  const fvm_gnum_t  *ent_global_num;

  cs_int_t _n_location_vals = n_location_vals;

  timing[0] = bft_timer_wtime();

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
  case CS_TYPE_cs_int_t:
    elt_type = (sizeof(cs_int_t) == 8) ? FVM_INT64 : FVM_INT32;
    break;
  case CS_TYPE_cs_real_t:
    elt_type =   (sizeof(cs_real_t) == fvm_datatype_size[FVM_DOUBLE])
               ? FVM_DOUBLE : FVM_FLOAT;
    break;
  default:
    assert(val_type == CS_TYPE_cs_int_t || val_type == CS_TYPE_cs_real_t);
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


  else if (cs_glob_n_ranks == 1) {

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

  timing[1] = bft_timer_wtime();
  _restart_wtime[restart->mode] += timing[1] - timing[0];

#endif /* #if defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------
 * Print statistics associated with restart files
 *----------------------------------------------------------------------------*/

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
