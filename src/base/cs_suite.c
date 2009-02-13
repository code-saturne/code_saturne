/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_file.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

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

#include "cs_suite.h"

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

#define CS_SUITE_NAME_LEN   64

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

struct _cs_suite_t {

  char            *name;         /* Name of restart file */

  cs_io_t         *fh;           /* Pointer to associated file handle */

  size_t           n_locations;  /* Number of locations */
  _location_t     *location;     /* Location definition array */

  cs_suite_mode_t  mode;         /* Read or write */
};


/*============================================================================
 * Static global variables
 *============================================================================*/

/* Minimum buffer size on rank 0 (to limit number of blocks
   when there is a large number of processors) */

static int cs_suite_taille_buf_def = 1024*1024*8;

/* Array for Fortran API */

static size_t       _restart_pointer_size = 0;
static cs_suite_t **_restart_pointer = NULL;


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute number of values in a record
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   location_id     <-- location id
 *   n_location_vals <-- number of values per location
 *----------------------------------------------------------------------------*/

static size_t
_compute_n_ents(const cs_suite_t  *suite,
                size_t             location_id,
                size_t             n_location_vals)
{
  size_t retval = 0;

  if (location_id == 0)
    retval = n_location_vals;

  else if (location_id > 0 && location_id <= suite->n_locations)
    retval = suite->location[location_id-1].n_glob_ents_f * n_location_vals;

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Location number %d given for restart file\n"
                "\"%s\" is not valid."),
              location_id, suite->name);

  return retval;
}

/*----------------------------------------------------------------------------
 * Analyze the content of a restart file to build locations
 *
 * parameters:
 *   suite    <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

static void
_locations_from_index(cs_suite_t  *suite)
{
  cs_io_sec_header_t h;

  size_t rec_id = 0;
  size_t index_size = 0;

  /* Initialization */

  index_size = cs_io_get_index_size(suite->fh);

  /* Analyze records to determine locations */

  for (rec_id = 0; rec_id < index_size; rec_id++) {

    h = cs_io_get_indexed_sec_header(suite->fh, rec_id);

    if (h.location_id > suite->n_locations) {

      _location_t  *loc = NULL;

      if (h.location_id != suite->n_locations + 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Restart file \"%s\" declares a location number %d\n"
                    "but no location %d has been declared."),
                  suite->name, (int)(h.location_id),
                  (int)(suite->n_locations + 1));

      BFT_REALLOC(suite->location, suite->n_locations + 1, _location_t);

      loc = suite->location + suite->n_locations;
      BFT_MALLOC(loc->name, strlen(h.sec_name) + 1, char);
      strcpy(loc->name, h.sec_name);

      loc->id = h.location_id;
      loc->n_ents = 0;
      loc->n_glob_ents = 0;

      cs_io_set_indexed_position(suite->fh, &h, rec_id);
      cs_io_set_fvm_gnum(&h, suite->fh);
      cs_io_read_global(&h, &(loc->n_glob_ents_f), suite->fh);

      loc->ent_global_num = NULL;

      suite->n_locations += 1;
    }

  }
}

/*----------------------------------------------------------------------------
 * Initialize a checkpoint / restart file management structure;
 *
 * parameters:
 *   suite <-> associated restart file pointer
 *----------------------------------------------------------------------------*/

static void
_add_file(cs_suite_t  *suite)
{
  const char magic_string[] = "Checkpoint / restart, R0";
  const long echo = CS_IO_ECHO_NONE;

  /* In read mode, open file to detect header first */

  if (suite->mode == CS_SUITE_MODE_LECTURE) {

#if defined(FVM_HAVE_MPI)
    suite->fh = cs_io_initialize_with_index(suite->name,
                                            magic_string,
                                            0,
                                            echo,
                                            cs_glob_base_mpi_comm);
#else
    suite->fh = cs_io_initialize_with_index(suite->name, magic_string, 0, echo);
#endif

    _locations_from_index(suite);
  }

  else {

#if defined(FVM_HAVE_MPI)
    suite->fh = cs_io_initialize(suite->name,
                                 magic_string,
                                 CS_IO_MODE_WRITE,
                                 0,
                                 echo,
                                 cs_glob_base_mpi_comm);
#else
    suite->fh = cs_io_initialize(suite->name,
                                 magic_string,
                                 CS_IO_MODE_WRITE,
                                 0,
                                 echo);
#endif
  }
}

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Read variable values defined on a mesh location.
 *
 * parameters:
 *   suite           <-> associated restart file pointer
 *   header          <-- header associated with current position in file
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers (1 to n numbering)
 *   n_location_vals <-- number of values par location
 *   datatype        <-- data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_read_ent_values(cs_suite_t                *suite,
                 cs_io_sec_header_t        *header,
                 fvm_gnum_t                 n_glob_ents,
                 fvm_lnum_t                 n_ents,
                 const fvm_gnum_t           ent_global_num[],
                 int                        n_location_vals,
                 cs_type_t                  datatype,
                 cs_byte_t                  vals[])
{
  cs_byte_t  *buffer = NULL;

  fvm_lnum_t  block_buf_size = 0;

  size_t  nbr_byte_ent, nbr_byte_val;

  fvm_block_to_part_info_t bi;

  fvm_block_to_part_t *d = NULL;

  /* Initialization */

  switch (datatype) {
  case CS_TYPE_cs_int_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_int_t);
    nbr_byte_val = sizeof(cs_int_t);
    cs_io_set_fvm_lnum(header, suite->fh);
    break;
  case CS_TYPE_cs_real_t:
    nbr_byte_ent = n_location_vals * sizeof(cs_real_t);
    nbr_byte_val = sizeof(cs_real_t);
    break;
  default:
    assert(datatype == CS_TYPE_cs_int_t || datatype == CS_TYPE_cs_real_t);
  }

  bi = fvm_block_to_part_compute_sizes(cs_glob_base_rang,
                                       cs_glob_base_nbr,
                                       cs_suite_taille_buf_def / nbr_byte_ent,
                                       n_glob_ents);

  d = fvm_block_to_part_create_strided(cs_glob_base_mpi_comm,
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
                   suite->fh);

 /* Distribute blocks on ranks */

  fvm_block_to_part_copy_array(d,
                               header->elt_type,
                               n_location_vals,
                               buffer,
                               vals);

  BFT_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Write variable values defined on a mesh location.
 *
 * parameters:
 *   suite           <-> associated restart file pointer
 *   sec_name        <-- section name
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers (1 to n numbering)
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values par location
 *   datatype        <-- data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_write_ent_values(const cs_suite_t  *suite,
                  const char        *sec_name,
                  fvm_gnum_t         n_glob_ents,
                  fvm_lnum_t         n_ents,
                  const fvm_gnum_t  *ent_global_num,
                  int                location_id,
                  int                n_location_vals,
                  cs_type_t          datatype,
                  const cs_byte_t   *vals)
{
  fvm_lnum_t  block_buf_size = 0;


  fvm_datatype_t elt_type = FVM_DATATYPE_NULL;
  size_t      nbr_byte_ent;
  cs_byte_t  *buffer = NULL;

  fvm_part_to_block_info_t bi;

  fvm_part_to_block_t *d = NULL;

  /* Initialization */

  switch (datatype) {
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
    assert(datatype == CS_TYPE_cs_int_t || datatype == CS_TYPE_cs_real_t);
  }

  bi = fvm_part_to_block_compute_sizes(cs_glob_base_rang,
                                       cs_glob_base_nbr,
                                       cs_suite_taille_buf_def / nbr_byte_ent,
                                       n_glob_ents);

  d = fvm_part_to_block_create_strided(cs_glob_base_mpi_comm,
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
                           suite->fh);

  /* Free buffer */

  BFT_FREE(buffer);

  fvm_part_to_block_destroy(&d);
}

#endif /* #if defined(_CS_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Convert read/write arguments from the Fortran API to the C API.
 *
 * parameters:
 *   numsui   <-- restart file id
 *   itysup   <-- location type code
 *   irtype   <-- integer or real
 *   suite    <-- pointer to restart file handle
 *   support  <-- location id
 *   datatype <-- integer of real
 *   ierror   <-- 0 = success, < 0 = error
 *----------------------------------------------------------------------------*/

static void
_section_f77_to_c(const cs_int_t   *numsui,
                  const cs_int_t   *itysup,
                  const cs_int_t   *irtype,
                  cs_suite_t      **suite,
                  int              *support,
                  cs_type_t        *datatype,
                  cs_int_t         *ierror)
{
  cs_int_t indsui = *numsui - 1;

  *ierror = CS_SUITE_SUCCES;

  /* Pointer to associated restart file handle */

  if (   indsui < 0
      || indsui > (cs_int_t)_restart_pointer_size
      || _restart_pointer[indsui] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Restart file number <%d> can not be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_SUITE_ERR_NUM_FIC;
    return;
  }

  else
    *suite = _restart_pointer[indsui];

  /* Location associated with section */

  switch (*itysup) {

  case 0:
    *support = CS_SUITE_SUPPORT_SCAL;
    break;

  case 1:
    *support = CS_SUITE_SUPPORT_CEL;
    break;

  case 2:
    *support = CS_SUITE_SUPPORT_FAC_INT;
    break;

  case 3:
    *support = CS_SUITE_SUPPORT_FAC_BRD;
    break;

  case 4:
    *support = CS_SUITE_SUPPORT_SOM;
    break;

  default:
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Location type <%d> given for a restart file section\n"
                 "is invalid using the Fortran API."), (int)(*itysup));
    *ierror = CS_SUITE_ERR_SUPPORT;
    return;

  }

  /* Datatype associated with section */

  switch (*irtype) {

  case 1:
    *datatype = CS_TYPE_cs_int_t;
    break;

  case 2:
    *datatype = CS_TYPE_cs_real_t;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Value type <%d> given for a restart file section\n"
                "is invalid using the Fortran API."), (int)(*irtype));
    *ierror = CS_SUITE_ERR_TYPE_VAL;
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
 *   datatype        --> data type
 *   vals            --> array of values
 *----------------------------------------------------------------------------*/

static void
_restart_permute_read(cs_int_t           n_ents,
                      const fvm_gnum_t  *ini_ent_num,
                      cs_int_t           n_location_vals,
                      cs_type_t          datatype,
                      cs_byte_t         *vals)
{
  cs_int_t ent_id, jj;

  cs_int_t ii = 0;

  /* Instructions */

  if (ini_ent_num == NULL)
    return;

  switch (datatype) {

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
    assert(datatype == CS_TYPE_cs_int_t || datatype == CS_TYPE_cs_real_t);

  }
}

/*----------------------------------------------------------------------------
 * Swap values of a renumbered array when writing
 *
 * parameters:
 *   n_ents          --> number of local entities
 *   ini_ent_num     --> initial entity numbers
 *   n_location_vals --> number of values per entity
 *   datatype        --> data type
 *   vals            --> array of values
 *
 * returns:
 *   pointer to array of values in initial entity order
 *----------------------------------------------------------------------------*/

static cs_byte_t *
_restart_permute_write(cs_int_t           n_ents,
                       const fvm_gnum_t  *ini_ent_num,
                       cs_int_t           n_location_vals,
                       cs_type_t          datatype,
                       const cs_byte_t   *vals)
{
  cs_int_t  ent_id, jj;

  cs_int_t  ii = 0;

  /* Instructions */

  if (ini_ent_num == NULL)
    return NULL;

  switch (datatype) {

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
    assert(datatype == CS_TYPE_cs_int_t || datatype == CS_TYPE_cs_real_t);
    return NULL;

  }
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

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
  char    *nombuf;

  size_t   id;

  cs_suite_mode_t suite_mode;


  /* Initialization */

  *numsui = 0;
  *ierror = CS_SUITE_SUCCES;

  /* Handle name for C API */

  nombuf = cs_base_string_f_to_c_create(nomsui, *lngnom);

  /* File creation options */

  {
    switch(*ireawr) {
    case 1:
      suite_mode = CS_SUITE_MODE_LECTURE;
      break;
    case 2:
      suite_mode = CS_SUITE_MODE_ECRITURE;
      break;
    default:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The access mode of the restart file <%s>\n"
                   "must be equal to 1 (read) or 2 (write) and not <%d>."),
                 nombuf, (int)(*ireawr));

      *ierror = CS_SUITE_ERR_MODE;
    }

  }

  /* Search for an available slot */

  if (*ierror == CS_SUITE_SUCCES) {

    for (id = 0;
         id < _restart_pointer_size && _restart_pointer[id] != NULL;
         id++);

    /* If no slot is available, we allow for more restart files */

    if (id == _restart_pointer_size) {

      BFT_REALLOC(_restart_pointer, _restart_pointer_size * 2,
                  cs_suite_t *);
      for (id = _restart_pointer_size;
           id < _restart_pointer_size * 2;
           id++)
        _restart_pointer[id] = NULL;
      _restart_pointer_size *= 2;

    }

  }

  /* Create file */

  if (*ierror == CS_SUITE_SUCCES)
    _restart_pointer[id] = cs_suite_cree(nombuf, suite_mode);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&nombuf);

  /*
   * Return the position of the handle in the array
   * (id + 1 to have a 1 to n numbering, more conventional in F77)
  */

  if (*ierror == CS_SUITE_SUCCES)
    *numsui = id + 1;
  else
    *numsui = -1;
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
  cs_int_t indsui = *numsui - 1;

  *ierror = CS_SUITE_SUCCES;

  /* Check that the file is valid */

  if (   indsui < 0
      || indsui > (cs_int_t)_restart_pointer_size
      || _restart_pointer[indsui] == NULL) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Restart file number <%d> can not be closed\n"
                 "(file already closed or invalid number)."), (int)(*numsui));

    *ierror = CS_SUITE_ERR_NUM_FIC;
    return;
  }

  /* Close file */

  cs_suite_detruit(_restart_pointer[indsui]);
  _restart_pointer[indsui] = NULL;
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
  cs_bool_t  corresp_cel, corresp_fac, corresp_fbr, corresp_som;

  cs_int_t   indsui   = *numsui - 1;

  /* Associated structure pointer */

  if (   indsui < 0
      || indsui > (cs_int_t)_restart_pointer_size
      || _restart_pointer[indsui] == NULL) {

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

    cs_suite_verif_support_base(_restart_pointer[indsui],
                                &corresp_cel, &corresp_fac,
                                &corresp_fbr, &corresp_som);

    *indcel = (corresp_cel == true ? 1 : 0);
    *indfac = (corresp_fac == true ? 1 : 0);
    *indfbr = (corresp_fbr == true ? 1 : 0);
    *indsom = (corresp_som == true ? 1 : 0);

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
  cs_int_t   indsui   = *numsui - 1;

  /* Associated structure pointer */

  if (   indsui < 0
      || indsui > (cs_int_t)_restart_pointer_size
      || _restart_pointer[indsui] == NULL) {

    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Information on the restart file number <%d> unavailable\n"
                 "(file already closed or invalid number)."), (int)(*numsui));
  }
  else {

    cs_suite_affiche_index(_restart_pointer[indsui]);

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
  char    *nombuf;

  cs_type_t   datatype;

  cs_suite_t  *suite;
  int          location_id;


  *ierror = CS_SUITE_SUCCES;

  /* Handle name for C API */

  nombuf = cs_base_string_f_to_c_create(nomrub, *lngnom);

  /* Handle other arguments for C API */

  _section_f77_to_c(numsui,
                    itysup,
                    irtype,
                    &suite,
                    &location_id,
                    &datatype,
                    ierror);

  if (*ierror < CS_SUITE_SUCCES)
    return;

  /* Read section */

  *ierror = cs_suite_lit_rub(suite,
                             nombuf,
                             location_id,
                             *nbvent,
                             datatype,
                             tabvar);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&nombuf);
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
  char    *nombuf;

  cs_type_t     datatype;

  cs_suite_t   *suite;
  int           location_id;


  *ierror = CS_SUITE_SUCCES;

  /* Handle name for C API */

  nombuf = cs_base_string_f_to_c_create(nomrub, *lngnom);

  /* Handle other arguments for C API */

  _section_f77_to_c(numsui,
                    itysup,
                    irtype,
                    &suite,
                    &location_id,
                    &datatype,
                    ierror);

  if (*ierror < CS_SUITE_SUCCES)
    return;

  /* Write section */

  cs_suite_ecr_rub(suite,
                   nombuf,
                   location_id,
                   *nbvent,
                   datatype,
                   tabvar);

  /* Free memory if necessary */

  cs_base_string_f_to_c_free(&nombuf);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a restart file
 *
 * parameters:
 *   nom  <-- file name
 *   mode <-- read or write
 *
 * returns:
 *   pointer to initialized restart file structure
 *----------------------------------------------------------------------------*/

cs_suite_t *
cs_suite_cree(const char             *nom,
              cs_suite_mode_t         mode)
{
  cs_suite_t  * suite;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Allocate and initialize base structure */

  BFT_MALLOC(suite, 1, cs_suite_t);

  BFT_MALLOC(suite->name, strlen(nom) + 1, char);

  strcpy(suite->name, nom);

  /* Initialize other fields */

  suite->mode = mode;

  suite->fh = NULL;

  /* Initialize location data */

  suite->n_locations = 0;
  suite->location = NULL;

  /* Open associated file, and build an index of sections in read mode */

  _add_file(suite);

  /* Add basic location definitions */

  cs_suite_ajoute_support(suite, "cells",
                          mesh->n_g_cells, mesh->n_cells,
                          mesh->global_cell_num);
  cs_suite_ajoute_support(suite, "interior_faces",
                          mesh->n_g_i_faces, mesh->n_i_faces,
                          mesh->global_i_face_num);
  cs_suite_ajoute_support(suite, "boundary_faces",
                          mesh->n_g_b_faces, mesh->n_b_faces,
                          mesh->global_b_face_num);
  cs_suite_ajoute_support(suite, "vertices",
                          mesh->n_g_vertices, mesh->n_vertices,
                          mesh->global_vtx_num);

  return suite;
}

/*----------------------------------------------------------------------------
 * Destroy structure associated with a restart file (and close the file).
 *
 * parameters:
 *   suite <-- pointer to restart file structure
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

cs_suite_t *
cs_suite_detruit(cs_suite_t  *suite)
{
  assert(suite != NULL);

  if (suite->fh != NULL)
    cs_io_finalize(&(suite->fh));

  /* Free locations array */

  if (suite->n_locations > 0) {
    size_t loc_id;
    for (loc_id = 0; loc_id < suite->n_locations; loc_id++)
      BFT_FREE((suite->location[loc_id]).name);
  }
  if (suite->location != NULL)
    BFT_FREE(suite->location);

  /* Free remaining memory */

  BFT_FREE(suite->name);
  BFT_FREE(suite);

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
 *   suite        <-- associated restart file pointer
 *   corresp_cell <-- matching cells flag
 *   corresp_fac  <-- matching interior faces flag
 *   corresp_fbr  <-- matching boundary faces flag
 *   corresp_som  <-- matching vertices flag
 *----------------------------------------------------------------------------*/

void
cs_suite_verif_support_base(const cs_suite_t  *suite,
                            cs_bool_t         *corresp_cel,
                            cs_bool_t         *corresp_fac,
                            cs_bool_t         *corresp_fbr,
                            cs_bool_t         *corresp_som)
{
  size_t location_id;

  *corresp_cel = false;
  *corresp_fac = false;
  *corresp_fbr = false;
  *corresp_som = false;

  assert(suite != NULL);

  for (location_id = 0; location_id < 4; location_id++) {

    const _location_t *loc = suite->location + location_id;

    if (loc->n_glob_ents_f == loc->n_glob_ents) {
      if (location_id == 0)
        *corresp_cel = true;
      else if (location_id == 1)
        *corresp_fac = true;
      else if (location_id == 2)
        *corresp_fbr = true;
      else if (location_id == 3)
        *corresp_som = true;
    }

    else if (cs_glob_base_rang <= 0) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The size of location \"%s\" associated with\n"
                   "the restart file \"%s\" is %lu and does not\n"
                   "correspond to that of the current mesh (%lu).\n"),
                 loc->name, suite->name,
                 (unsigned long)loc->n_glob_ents_f,
                 (unsigned long)loc->n_glob_ents);
    }

  }
}

/*----------------------------------------------------------------------------
 * Add a location definition.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   location_name   <-- name associated with the location
 *   n_glob_ents     <-- global number of entities
 *   n_ents          <-- local number of entities
 *   ent_global_num  <-- global entity numbers, or NULL
 *
 * returns:
 *   the location id assigned to the location, or -1 in case of error
 *----------------------------------------------------------------------------*/

int
cs_suite_ajoute_support(cs_suite_t        *suite,
                        const char        *location_name,
                        fvm_gnum_t         n_glob_ents,
                        fvm_lnum_t         n_ents,
                        const fvm_gnum_t  *ent_global_num)
{
  int loc_id;

  if (suite->mode == CS_SUITE_MODE_LECTURE) {

    /* Search for a location with the same name */

    for (loc_id = 0; loc_id < (int)(suite->n_locations); loc_id++) {

      if ((strcmp((suite->location[loc_id]).name, location_name) == 0)) {

        (suite->location[loc_id]).n_glob_ents = n_glob_ents;

        (suite->location[loc_id]).n_ents  = n_ents;
        (suite->location[loc_id]).ent_global_num = ent_global_num;

        return loc_id + 1;

      }
    }

    if (loc_id >= ((int)(suite->n_locations)))
      bft_error(__FILE__, __LINE__, 0,
                _("The restart file \"%s\" references no\n"
                  "location named \"%s\"."),
                suite->name, location_name);

  }

  else {

    fvm_datatype_t gnum_type
      = (sizeof(fvm_gnum_t) == 8) ? FVM_UINT64 : FVM_UINT32;

    /* Create a new location */

    suite->n_locations += 1;

    BFT_REALLOC(suite->location, suite->n_locations, _location_t);
    BFT_MALLOC((suite->location[suite->n_locations-1]).name,
               strlen(location_name)+1,
               char);

    strcpy((suite->location[suite->n_locations-1]).name, location_name);

    (suite->location[suite->n_locations-1]).id             = suite->n_locations;
    (suite->location[suite->n_locations-1]).n_glob_ents    = n_glob_ents;
    (suite->location[suite->n_locations-1]).n_glob_ents_f  = n_glob_ents;
    (suite->location[suite->n_locations-1]).n_ents         = n_ents;
    (suite->location[suite->n_locations-1]).ent_global_num = ent_global_num;

    cs_io_write_global(location_name, 1, suite->n_locations, 0, 0,
                       gnum_type, &n_glob_ents,
                       suite->fh);

    return suite->n_locations;
  }

  return -1;
}

/*----------------------------------------------------------------------------
 * Print the index associated with a restart file in read mode
 *
 * parameters:
 *   suite <-- associated restart file pointer
 *----------------------------------------------------------------------------*/

void
cs_suite_affiche_index(const cs_suite_t  *suite)
{
  size_t loc_id;

  assert(suite != NULL);

  for (loc_id = 0; loc_id < suite->n_locations; loc_id++) {
    const _location_t *loc = &(suite->location[loc_id]);
    bft_printf(_("  Location: %s\n"
                 "    (number: %03d, n_glob_ents: %lu)\n"),
               loc->name, (int)(loc->id), (unsigned long)(loc->n_glob_ents));
  }
  if (suite->n_locations > 0)
    bft_printf("\n");

  /* Dump general file info, including index */

  bft_printf(_("  General information associated with the restart file:\n"));

  cs_io_dump(suite->fh);
}


/*----------------------------------------------------------------------------
 * Read a section from a restart file.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   nom_rub         <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   typ_val         <-- value type
 *   val             --> array of values
 *
 * returns: 0 (CS_SUITE_SUCCES) in case of success,
 *          or error code (CS_SUITE_ERR_xxx) in case of error
 *----------------------------------------------------------------------------*/

int
cs_suite_lit_rub(cs_suite_t  *suite,
                 const char  *nom_rub,
                 int          location_id,
                 cs_int_t     n_location_vals,
                 cs_type_t    typ_val,
                 void        *val)
{
  cs_int_t   n_glob_ents, n_ents;
  const fvm_gnum_t  *ent_global_num;

  size_t rec_id;
  cs_io_sec_header_t header;

  cs_int_t _n_location_vals = n_location_vals;
  size_t index_size = 0;

  index_size = cs_io_get_index_size(suite->fh);

  assert(suite != NULL);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = NULL;
  }

  else {
    if (location_id < 0 || location_id > (int)(suite->n_locations))
      return CS_SUITE_ERR_SUPPORT;
    if (   (suite->location[location_id-1]).n_glob_ents_f
        != (suite->location[location_id-1]).n_glob_ents)
      return CS_SUITE_ERR_SUPPORT;
    n_glob_ents = (suite->location[location_id-1]).n_glob_ents;
    n_ents  = (suite->location[location_id-1]).n_ents;
    ent_global_num = (suite->location[location_id-1]).ent_global_num;
  }

  /* Search for the corresponding record in the index */

  for (rec_id = 0; rec_id < index_size; rec_id++) {
    const char * cmp_name = cs_io_get_indexed_sec_name(suite->fh, rec_id);
    if (strcmp(cmp_name, nom_rub) == 0)
      break;
  }

  /* If the record was not found */

  if (rec_id >= index_size)
    return CS_SUITE_ERR_EXISTE;

  /*
    If the location does not fit: we search for a location of same
    name with the correct location.
  */

  header = cs_io_get_indexed_sec_header(suite->fh, rec_id);

  if (header.location_id != (size_t)location_id) {

    rec_id++;

    while (rec_id < index_size) {
      header = cs_io_get_indexed_sec_header(suite->fh, rec_id);
      if (   (strcmp(header.sec_name, nom_rub) == 0)
          && (header.location_id == (size_t)location_id))
        break;
      rec_id++;
    }

    if (rec_id >= index_size)
      return CS_SUITE_ERR_SUPPORT;
  }

  /* If the number of values per location does not match */

  if (   (   header.location_id > 0
          && header.n_location_vals != (size_t)n_location_vals)
      || (   header.location_id == 0 && header.n_vals != n_ents))
    return CS_SUITE_ERR_NBR_VAL;

  /* If the type of value does not match */

  if (header.elt_type == FVM_INT32 || header.elt_type == FVM_INT64) {
    cs_io_set_fvm_lnum(&header, suite->fh);
    if (typ_val != CS_TYPE_cs_int_t)
      return CS_SUITE_ERR_TYPE_VAL;
  }
  else if (header.elt_type == FVM_FLOAT || header.elt_type == FVM_DOUBLE) {
    if (sizeof(cs_real_t) != fvm_datatype_size[header.elt_type]) {
      if (sizeof(cs_real_t) == fvm_datatype_size[FVM_FLOAT])
        header.elt_type = FVM_FLOAT;
      else
        header.elt_type = FVM_DOUBLE;
    }
    if (typ_val != CS_TYPE_cs_real_t)
      return CS_SUITE_ERR_TYPE_VAL;
  }

  /* Now set position in file to read data */

  cs_io_set_indexed_position(suite->fh, &header, rec_id);

  /* Section contents */
  /*------------------*/

  /* In single processor mode or for global values */

  if (cs_glob_base_nbr == 1 || location_id == 0) {

    cs_io_read_global(&header, val, suite->fh);

    if (ent_global_num != NULL)
      _restart_permute_read(n_ents,
                            ent_global_num,
                            _n_location_vals,
                            typ_val,
                            val);
  }

#if defined(_CS_HAVE_MPI)

  /* In parallel mode for a distributed mesh location */

  else
    _read_ent_values(suite,
                     &header,
                     n_glob_ents,
                     n_ents,
                     ent_global_num,
                     _n_location_vals,
                     typ_val,
                     (cs_byte_t *)val);

#endif /* #if defined(_CS_HAVE_MPI) */

  /* Return */

  return CS_SUITE_SUCCES;
}

/*----------------------------------------------------------------------------
 * Write a section to a restart file.
 *
 * parameters:
 *   suite           <-- associated restart file pointer
 *   nom_rub         <-- section name
 *   location_id     <-- id of corresponding location
 *   n_location_vals <-- number of values per location (interlaced)
 *   typ_val         <-- value type
 *   val             <-- array of values
 *----------------------------------------------------------------------------*/

void
cs_suite_ecr_rub(cs_suite_t   *suite,
                 const char   *nom_rub,
                 int           location_id,
                 cs_int_t      n_location_vals,
                 cs_type_t     typ_val,
                 const void   *val)
{
  cs_int_t         n_tot_vals, n_glob_ents, n_ents;
  fvm_datatype_t   elt_type;

  const fvm_gnum_t  *ent_global_num;

  cs_int_t _n_location_vals = n_location_vals;

  assert(suite != NULL);

  n_tot_vals = _compute_n_ents(suite, location_id, n_location_vals);

  /* Check associated location */

  if (location_id == 0) {
    n_glob_ents = n_location_vals;
    n_ents  = n_location_vals;
    _n_location_vals = 1;
    ent_global_num = NULL;
  }

  else {
    assert(location_id >= 0 && location_id <= (int)(suite->n_locations));
    n_glob_ents = (suite->location[location_id-1]).n_glob_ents;
    n_ents  = (suite->location[location_id-1]).n_ents;
    ent_global_num = (suite->location[location_id-1]).ent_global_num;
  }

  /* Set datatype */

  switch (typ_val) {
  case CS_TYPE_cs_int_t:
    elt_type = (sizeof(cs_int_t) == 8) ? FVM_INT64 : FVM_INT32;
    break;
  case CS_TYPE_cs_real_t:
    elt_type =   (sizeof(cs_real_t) == fvm_datatype_size[FVM_DOUBLE])
               ? FVM_DOUBLE : FVM_FLOAT;
    break;
  default:
    assert(typ_val == CS_TYPE_cs_int_t || typ_val == CS_TYPE_cs_real_t);
  }

  /* Section contents */
  /*------------------*/

  /* In single processor mode of for global values */

  if (location_id == 0)
    cs_io_write_global(nom_rub,
                       n_tot_vals,
                       location_id,
                       0,
                       1,
                       elt_type,
                       val,
                       suite->fh);


  else if (cs_glob_base_nbr == 1) {

    cs_byte_t  *val_tmp = NULL;

    if (ent_global_num != NULL)
      val_tmp = _restart_permute_write(n_ents,
                                       ent_global_num,
                                       _n_location_vals,
                                       typ_val,
                                       val);

    cs_io_write_global(nom_rub,
                       n_tot_vals,
                       location_id,
                       0,
                       _n_location_vals,
                       elt_type,
                       (val_tmp != NULL) ? val_tmp : val,
                       suite->fh);

    if (val_tmp != NULL)
      BFT_FREE (val_tmp);
  }

#if defined(_CS_HAVE_MPI)

  /* In parallel mode for a distributed mesh location */

  else
    _write_ent_values(suite,
                      nom_rub,
                      n_glob_ents,
                      n_ents,
                      ent_global_num,
                      location_id,
                      _n_location_vals,
                      typ_val,
                      (const cs_byte_t *)val);

#endif /* #if defined(_CS_HAVE_MPI) */
}

/*----------------------------------------------------------------------------
 * Initialize the restart file Fortran API
 *----------------------------------------------------------------------------*/

void
cs_suite_f77_api_init(void)
{
  size_t ind;

  /* Allocation du tableau des pointeurs */

  _restart_pointer_size = 10;
  BFT_MALLOC(_restart_pointer, _restart_pointer_size, cs_suite_t *);

  /* Set pointers array to NULL */

  for (ind = 0; ind < _restart_pointer_size; ind++)
    _restart_pointer[ind] = NULL;
}


/*----------------------------------------------------------------------------
 * Finalize the restart file Fortran API
 *----------------------------------------------------------------------------*/

void
cs_suite_f77_api_finalize(void)
{
  size_t ind;

  /* Close files that are not closed yet */

  for (ind = 0; ind < _restart_pointer_size; ind++) {
    if (_restart_pointer[ind] != NULL)
      cs_suite_detruit (_restart_pointer[ind]);
  }

  /* Free array of pointers */

  _restart_pointer_size = 0;
  BFT_FREE(_restart_pointer);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
