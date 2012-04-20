/*============================================================================
 * \file Field management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/* Field key definitions */

typedef struct {

  unsigned char      def_val[8];   /* Default value container (int, double,
                                      or pointer to string) */
  int                type_flag;    /* Type flag */
  char               type_id;      /* i: int; d: double; s: str */
  bool               is_sub;       /* Indicate if the key is a sub-key (in
                                      which case def_val contains
                                      the parent key id */

} cs_field_key_def_t;

/* Field key value structures */

typedef struct {

  unsigned char      val[8];       /* Value container (int, double,
                                      or pointer) */
  bool               is_set;       /* Has this key been set for the
                                      present field ? */

} cs_field_key_val_t;

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Field definitions */

static int  _n_fields = 0;
static int  _n_fields_max = 0;
static cs_field_t  *_fields = NULL;
static cs_map_name_to_id_t  *_field_map = NULL;

/* Key definitions */

static int  _n_keys = 0;
static int  _n_keys_max = 0;
static cs_field_key_def_t  *_key_defs = NULL;
static cs_map_name_to_id_t  *_key_map = NULL;

/* Key values : _key_vals[field_id*_n_keys_max + key_id] */

static cs_field_key_val_t  *_key_vals = NULL;

/* Names for logging */

static const int _n_type_flags = 6;
static const int _type_flag_mask[] = {CS_FIELD_INTENSIVE,
                                      CS_FIELD_EXTENSIVE,
                                      CS_FIELD_VARIABLE,
                                      CS_FIELD_PROPERTY,
                                      CS_FIELD_POSTPROCESS,
                                      CS_FIELD_ACCUMULATOR,
                                      CS_FIELD_USER};
static const char *_type_flag_name[] = {N_("intensive"),
                                        N_("extensive"),
                                        N_("variable"),
                                        N_("property"),
                                        N_("postprocess"),
                                        N_("accumulator"),
                                        N_("user")};

/*============================================================================
 * Fortran function prototypes for subroutines from field.f90.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set global temporary scalar field pointer to null.
 *
 * function fldps2
 * ***************
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldps2, FLDPS2)
(
 void
);

/*----------------------------------------------------------------------------
 * Set global temporary scalar field pointer to a given array.
 *
 * function fldps3 (nval, val)
 * ***************
 *
 * integer          nval        : <-- : Number of field values
 * double precision val         : <-- : pointer to field values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldps3, FLDPS3)
(
 const cs_int_t   *nval,
 cs_real_t        *val
);

/*----------------------------------------------------------------------------
 * Set global temporary vector field pointer to null.
 *
 * function fldpv2
 * ***************
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldpv2, FLDPV2)
(
 void
);

/*----------------------------------------------------------------------------
 * Set global temporary vector field pointer to a given array.
 *
 * function fldpv3 (nval1, nval2, val)
 * ***************
 *
 * integer          nval1       : <-- : Number of values for first index
 * integer          nval2       : <-- : Number of values for second index
 * double precision val         : <-- : pointer to field values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldpv3, FLDPV3)
(
 const cs_int_t   *nval1,
 const cs_int_t   *nval2,
 cs_real_t        *val
);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Create a field descriptor.
 *
 * For fields with a dimension greater than 1, components are interleaved.
 *
 * parameters:
 *   name        <-- field name
 *   type_flag   <-- mask of field property and category values
 *   location_id <-- id of associated location
 *   dim         <-- field dimension (number of components)
 *   interleaved <-- if dim > 1, indicate if field is interleaved
 *
 * returns:
 *   pointer to new field.
 *
 *----------------------------------------------------------------------------*/

static cs_field_t *
_field_create(const char   *name,
              int           type_flag,
              int           location_id,
              int           dim,
              bool          interleaved)
{
  int key_id;
  int field_id = -1;
  const char *addr_0 = NULL, *addr_1 = NULL;

  cs_field_t *f =  NULL;

  /* Initialize if necessary */

  if (_field_map == NULL)
    _field_map = cs_map_name_to_id_create();

  else
    addr_0 = cs_map_name_to_id_reverse(_field_map, 0);

  if (strlen(name) == 0)
    bft_error(__FILE__, __LINE__, 0, _("Defining a field requires a name."));

  /* Find or insert entry in map */

  field_id = cs_map_name_to_id(_field_map, name);

  /* Move name pointers of previous fields if necessary
     (i.e. reallocation of map names array) */

  addr_1 = cs_map_name_to_id_reverse(_field_map, 0);

  if (addr_1 != addr_0) {
    int i;
    ptrdiff_t addr_shift = addr_1 - addr_0;
    for (i = 0; i < field_id; i++)
      (_fields + i)->name += addr_shift;
  }

  if (field_id == _n_fields)
    _n_fields = field_id + 1;

  /* Reallocate fields if necessary */

  if (_n_fields > _n_fields_max) {
    if (_n_fields_max == 0)
      _n_fields_max = 8;
    else
      _n_fields_max *= 2;
    BFT_REALLOC(_fields, _n_fields_max, cs_field_t);
    BFT_REALLOC(_key_vals, _n_keys_max*_n_fields_max, cs_field_key_val_t);
  }

  /* Check type flags and location id */

  if (   (type_flag & CS_FIELD_INTENSIVE)
      && (type_flag & CS_FIELD_EXTENSIVE))
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\"\n"
                " may not be defined as both intensive and extensive."),
              name);
  else if (location_id < 0 || location_id >= cs_mesh_location_n_locations())
    bft_error(__FILE__, __LINE__, 0,
              _("Mesh location %d associated with field \"%s\"\n"
                " has not been defined yet."),
              location_id, name);

  /* Assign field */

  f = _fields + field_id;

  f->name = cs_map_name_to_id_reverse(_field_map, field_id);

  f->id = field_id;
  f->type = type_flag;
  f->dim = dim;
  if (f->dim > 1)
    f->interleaved = interleaved;
  else
    f->interleaved = true;
  f->location_id = location_id;
  f->n_time_vals = 1;

  f->val = NULL;
  f->val_pre = NULL;

  f->bc_coeffs = NULL;

  f->is_owner = true;

  /* Mark key values as not set */

  for (key_id = 0; key_id < _n_keys_max; key_id++) {
    memset((_key_vals + (f->id*_n_keys_max + key_id))->val, 0, 8);
    (_key_vals + (f->id*_n_keys_max + key_id))->is_set = false;
  }

  return f;
}

/*----------------------------------------------------------------------------*
 * allocate and initialize a field values array.
 *
 * parameters:
 *   n_elts  <-- number of associated elements
 *   dim     <-- associated dimension
 *   val_old <-- pointer to previous array in case of reallocation
 *               (usually NULL)
 *
 * returns  pointer to new field values.
 *----------------------------------------------------------------------------*/

static cs_real_t *
_add_val(cs_lnum_t   n_elts,
         int         dim,
         cs_real_t  *val_old)
{
  cs_lnum_t  ii;
  cs_real_t  *val = val_old;

  BFT_REALLOC(val, n_elts*dim, cs_real_t);

  /* Initialize field. This should not be necessary, but when using
     threads with Open MP, this should help ensure that the memory will
     first be touched by the same core that will later operate on
     this memory, usually leading to better core/memory affinity. */

  if (dim == 1) {
#   pragma omp parallel for
    for (ii = 0; ii < n_elts; ii++)
      val[ii] = 0.0;
  }
  else {
    cs_lnum_t jj;
#   pragma omp parallel for private(jj)
    for (ii = 0; ii < n_elts; ii++) {
      for (jj = 0; jj < dim; jj++)
        val[ii*dim + jj] = 0.0;
    }
  }

  return val;
}

/*----------------------------------------------------------------------------
 * Find an id matching a key or define a new key and associated id.
 *
 * parameters:
 *   name <-- key name
 *
 * returns:
 *   id of associated key definition structure
 *----------------------------------------------------------------------------*/

static int
_find_or_add_key(const char  *name)
{
  int key_id;

  /* Initialize if necessary */

  if (_key_map == NULL)
    _key_map = cs_map_name_to_id_create();

  /* Find or insert entry in map */

  key_id = cs_map_name_to_id(_key_map, name);

  if (key_id == _n_keys)
    _n_keys = key_id + 1;

  /* Reallocate key definitions if necessary */

  if (_n_keys > _n_keys_max) {
    int field_id, _key_id;
    int _n_keys_max_prev = _n_keys_max;
    if (_n_keys_max == 0)
      _n_keys_max = 8;
    else
      _n_keys_max *= 2;
    BFT_REALLOC(_key_defs, _n_keys_max, cs_field_key_def_t);
    BFT_REALLOC(_key_vals, _n_keys_max*_n_fields_max, cs_field_key_val_t);
    for (field_id = _n_fields - 1; field_id >= 0; field_id--) {
      for (_key_id = _n_keys - 2; _key_id >= 0; _key_id--)
        _key_vals[field_id*_n_keys_max + _key_id]
          = _key_vals[field_id*_n_keys_max_prev + _key_id];
    }
    for (field_id = 0; field_id < _n_fields; field_id++) {
      memset((_key_vals + (field_id*_n_keys_max + key_id))->val, 0, 8);
      (_key_vals + (field_id*_n_keys_max + key_id))->is_set = false;
    }
  }

  return key_id;
}

/*----------------------------------------------------------------------------
 * Add type flag info to the current position in the setup log.
 *
 * parameters:
 *   type <-- type flag
 *----------------------------------------------------------------------------*/

static inline void
_log_add_type_flag(int type)
{
  int i;
  int n_loc_flags = 0;

  for (i = 0; i < _n_type_flags; i++) {
    if (type & _type_flag_mask[i]) {
      if (n_loc_flags == 0)
        cs_log_printf(CS_LOG_SETUP, " (%s", _(_type_flag_name[i]));
      else
        cs_log_printf(CS_LOG_SETUP, ", %s", _(_type_flag_name[i]));
      n_loc_flags++;
    }
  }

  if (n_loc_flags > 0)
    cs_log_printf(CS_LOG_SETUP, ")");
}

/*----------------------------------------------------------------------------
 * Free strings associated to a key.
 *----------------------------------------------------------------------------*/

static void
_cs_field_free_str(void)
{
  int key_id, f_id;

  for (key_id = 0; key_id < _n_keys; key_id++) {

    cs_field_key_def_t *kd = _key_defs + key_id;

    if (kd->type_id == 's') {
      for (f_id = 0; f_id < _n_fields; f_id++) {
        cs_field_key_val_t *kv = _key_vals + (f_id*_n_keys_max + key_id);
        char **s = (char **)(kv->val);
        BFT_FREE(*s);
      }
    }

  }
}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define a field.
 *
 * Fortran interface; use flddef (see cs_fieldt_f2c.f90)
 *
 * subroutine fldde1 (name, lname, iexten, itycat, ityloc, idim, ilved,
 * *****************
 *                    iprev, idfld)
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          iexten      : <-- : 1: intensive; 2: extensive
 * integer          itycat      : <-- : Field category (may be added)
 *                              :     :   4: variable
 *                              :     :   8: property
 *                              :     :  16: postprocess
 *                              :     :  32: accumulator
 *                              :     :  64: user
 * integer          ityloc      : <-- : Location type
 *                              :     :  0: none
 *                              :     :  1: cells
 *                              :     :  2: interior faces
 *                              :     :  3: interior faces
 *                              :     :  4: vertices
 * integer          idim        : <-- : Field dimension
 * integer          ilved       : <-- : 0: not intereaved; 1: interleaved
 * integer          iprev       : <-- : 0: no previous values, 1: previous
 * integer          ifield      : --> : id of defined field
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldde1, FLDDE1)
(
 const char       *name,
 const cs_int_t   *lname,
 const cs_int_t   *iexten,
 const cs_int_t   *itycat,
 const cs_int_t   *ityloc,
 const cs_int_t   *idim,
 const cs_int_t   *ilved,
 const cs_int_t   *iprev,
 cs_int_t         *ifield
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufname;
  int type_flag = 0;
  bool interleaved = (*ilved == 0) ? false : true;
  bool has_prev = (*iprev == 0) ? false : true;

  cs_field_t *f = NULL;

  bufname = cs_base_string_f_to_c_create(name, *lname);

  if (*iexten & 1)
    type_flag = CS_FIELD_EXTENSIVE;
  else if (*iexten & 2)
    type_flag = CS_FIELD_INTENSIVE;

  if (*itycat & 4)
    type_flag = type_flag | CS_FIELD_VARIABLE;
  if (*itycat & 8)
    type_flag = type_flag | CS_FIELD_PROPERTY;
  if (*itycat & 16)
    type_flag = type_flag | CS_FIELD_POSTPROCESS;
  if (*itycat & 32)
    type_flag = type_flag | CS_FIELD_ACCUMULATOR;
  if (*itycat == 64)
    type_flag = type_flag | CS_FIELD_USER;

  f = cs_field_create(bufname,
                      type_flag,
                      *ityloc,
                      *idim,
                      interleaved,
                      has_prev);

  cs_base_string_f_to_c_free(&bufname);

  *ifield = f->id;
}

/*----------------------------------------------------------------------------
 * Allocate field values
 *
 * Fortran interface
 *
 * subroutine fldalo (ifield)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldalo, FLDALO)
(
 const cs_int_t   *ifield
)
{
  cs_field_t *f = cs_field_by_id(*ifield);

  cs_field_allocate_values(f);
}

/*----------------------------------------------------------------------------
 * Map values to a field.
 *
 * Fortran interface
 *
 * subroutine fldmap (ifield, val, valp)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * cs_real_t*       val         : <-- : Pointer to field values array
 * cs_real_t*       valp        : <-- : Pointer to values at previous
 *                              :     : time step if field was defined
 *                              :     : with iprev = 1
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldmap, FLDMAP)
(
 const cs_int_t   *ifield,
 cs_real_t        *val,
 cs_real_t        *valp
)
{
  cs_field_t *f = cs_field_by_id(*ifield);

  cs_field_map_values(f, val, valp);
}

/*----------------------------------------------------------------------------
 * Map field boundary coefficient arrays.
 *
 * Fortran interface
 *
 * subroutine fldbcm (ifield, icpled, a, b, af, bf)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * cs_real_t*       a           : <-- : explicit BC coefficients array
 * cs_real_t*       b           : <-- : implicit BC coefficients array
 * cs_real_t*       af          : <-- : explicit flux BC coefficients array,
 *                              :     : or a (or NULL)
 * cs_real_t*       bf          : <-- : implicit flux BC coefficients array,
 *                              :     : or a (or NULL)
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldbcm, FLDBCM)
(
 const cs_int_t   *ifield,
 cs_real_t        *a,
 cs_real_t        *b,
 cs_real_t        *af,
 cs_real_t        *bf
)
{
  cs_real_t *_af = (af != a) ? af : NULL;
  cs_real_t *_bf = (bf != b) ? bf : NULL;

  cs_field_t *f = cs_field_by_id(*ifield);

  cs_field_map_bc_coeffs(f, a, b, _af, _bf);
}

/*----------------------------------------------------------------------------
 * Allocate arrays for all defined fields based on their location.
 *
 * Location sized must thus be known.
 *
 * Fields that do not own their data should all have been mapped at this
 * stage, and are checked.
 *
 * Fortran interface
 *
 * subroutine fldama
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldama, FLDAMA)
(
 void
)
{
  cs_field_allocate_or_map_all();
}

/*----------------------------------------------------------------------------
 * Retrieve field value pointer for a scalar field.
 *
 * Fortran interface; use fldpts
 *
 * function fldps1 (ifield, iprev)
 * ***************
 *
 * integer          ifield      : <-- : Field id
 * integer          iprev       : <-- : if 1, pointer to previous values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldps1, FLDPS1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *iprev
)
{
  cs_field_t *f = cs_field_by_id(*ifield);
  cs_real_t *val = NULL;

  if (*iprev == 0)
    val = (f->val_pre);
  else
    val = (f->val);

  if (val == NULL)
    CS_PROCF(fldps2, FLDPS2)();
  else {
    cs_int_t nval;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(*ifield);
    nval = n_elts[2];
    CS_PROCF(fldps3, FLDPS3)(&nval, val);
  }
}

/*----------------------------------------------------------------------------
 * Retrieve field value pointer for a vector field.
 *
 * Fortran interface; use fldpts
 *
 * function fldpv1 (ifield, iprev)
 * ***************
 *
 * integer          ifield      : <-- : Field id
 * integer          iprev       : <-- : if 1, pointer to previous values
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldpv1, FLDPV1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *iprev
)
{
  cs_field_t *f = cs_field_by_id(*ifield);
  cs_real_t *val = NULL;

  if (*iprev == 0)
    val = (f->val_pre);
  else
    val = (f->val);

  if (val == NULL)
    CS_PROCF(fldpv2, FLDPV2)();
  else {
    cs_int_t nval1, nval2;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(*ifield);
    if (f->interleaved) {
      nval1 = f->dim;
      nval2 = n_elts[2];
    }
    else {
      nval1 = n_elts[2];
      nval2 = f->dim;
    }
    CS_PROCF(fldpv3, FLDPV3)(&nval1, &nval2, val);
  }
}

/*----------------------------------------------------------------------------
 * Return an id associated with a given field name if present.
 *
 * If the field has not been defined previously, -1 is returned.
 *
 * Fortran interface
 *
 * subroutine fldfid (name,   lname,  ifield)
 * *****************
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          ifield      : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldfi1, FLDFI1)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ifield
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufname;

  bufname = cs_base_string_f_to_c_create(name, *lname);

  *ifield = cs_map_name_to_id_try(_field_map, bufname);

  cs_base_string_f_to_c_free(&bufname);
}

/*----------------------------------------------------------------------------
 * Return an id associated with a given key name if present.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * Fortran interface; use fldkid (see cs_fieldt_f2c.f90)
 *
 * subroutine fldki1 (name,   lname,  ikeyid)
 * *****************
 *
 * character*       name        : <-- : Key name
 * integer          lname       : <-- : Key name length
 * integer          ikey        : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldki1, FLDKI1)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ikey
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufname;

  bufname = cs_base_string_f_to_c_create(name, *lname);

  *ikey = cs_field_key_id_try(bufname);

  cs_base_string_f_to_c_free(&bufname);
}

/*----------------------------------------------------------------------------
 * Assign an integer value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldski (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * integer          value       : <-- : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldski, FLDSKI)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_int_t         *value
)
{
  int retval = 0;

  cs_field_t *f = cs_field_by_id(*ifield);

  retval = cs_field_set_key_int(f, *ikey, *value);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, *ikey);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning integer value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, *ikey, key);
  }
}

/*----------------------------------------------------------------------------
 * Return a integer value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldgki (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * integer          value       : --> : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgki, FLDGKI)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_int_t         *value
)
{
  const cs_field_t *f = cs_field_by_id(*ifield);
  *value = cs_field_get_key_int(f, *ikey);
}

/*----------------------------------------------------------------------------
 * Assign a floating point value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldskd (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * double precision value       : <-- : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldskd, FLDSKD)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_real_t        *value
)
{
  int retval = 0;

  cs_field_t *f = cs_field_by_id(*ifield);

  retval = cs_field_set_key_double(f, *ikey, *value);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, *ikey);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning real value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, *ikey, key);
  }
}

/*----------------------------------------------------------------------------
 * Return a floating point value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * subroutine fldgkd (ifield, ikey, value)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * double precision value       : --> : Associated value
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgkd, FLDGKD)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 cs_real_t        *value
)
{
  const cs_field_t *f = cs_field_by_id(*ifield);
  *value = cs_field_get_key_double(f, *ikey);
}

/*----------------------------------------------------------------------------
 * Assign a character string for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * Fortran interface; use fldsk1 (see cs_fieldt_f2c.f90)
 *
 * subroutine fldsk1 (ifield, ikey, str, lstr)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * character*       str         : <-- : Associated string
 * integer          lstr        : <-- : Associated string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldsk1, FLDSK1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 const char       *str,
 const cs_int_t   *lstr
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufstr;

  int retval = 0;

  cs_field_t *f = cs_field_by_id(*ifield);

  bufstr = cs_base_string_f_to_c_create(str, *lstr);

  retval = cs_field_set_key_str(f, *ikey, bufstr);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, *ikey);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning real value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, *ikey, key);
  }

  cs_base_string_f_to_c_free(&bufstr);
}

/*----------------------------------------------------------------------------
 * Return a character string for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * Fortran interface; use fldgk1 (see cs_fieldt_f2c.f90)
 *
 * subroutine fldgk1 (ifield, ikey, str, lstr)
 * *****************
 *
 * integer          ifield      : <-- : Field id
 * integer          ikey        : <-- : Key id
 * character*       str         : --> : Associated string
 * integer          lstr        : <-- : Associated string length
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldgk1, FLDGK1)
(
 const cs_int_t   *ifield,
 const cs_int_t   *ikey,
 char             *str,
 const cs_int_t   *lstr
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  cs_int_t  i, l;

  const char *s = NULL;
  const cs_field_t *f = cs_field_by_id(*ifield);
  s = cs_field_get_key_str(f, *ikey);

  l = strlen(s);

  for (i = 0; i < l && i < *lstr; i++)
    str[i] = s[i];
  for ( ; i < *lstr; i++)
    str[i] = ' ';
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of defined fields.
 *
 * \return  number of defined fields.
 */
/*----------------------------------------------------------------------------*/

int
cs_field_n_fields(void)
{
  return _n_fields;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a field descriptor.
 *
 * For fields with a dimension greater than 1, components are interleaved.
 *
 * \param[in]  name          field name
 * \param[in]  type_flag     mask of field property and category values
 * \param[in]  location_id   id of associated location
 * \param[in]  dim           field dimension (number of components)
 * \param[in]  interleaved   indicate if values ar interleaved
 *                           (ignored if number of components < 2)
 * \param[in]  has_previous  maintain values at the previous time step ?
 *
 * \return  pointer to new field.
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_field_create(const char   *name,
                int           type_flag,
                int           location_id,
                int           dim,
                bool          interleaved,
                bool          has_previous)
{
  cs_field_t  *f =  _field_create(name,
                                  type_flag,
                                  location_id,
                                  dim,
                                  interleaved);

  f->n_time_vals = has_previous ? 2 : 1;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate arrays for field values.
 *
 * \param[in, out]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_values(cs_field_t  *f)
{
  assert(f != NULL);

  if (f->is_owner) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);

    f->val = _add_val(n_elts[2], f->dim, f->val);

    /* Add previous time step values if necessary */
    if (f->n_time_vals > 1)
      f->val_pre = _add_val(n_elts[2], f->dim, f->val_pre);

  };
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map existing values to field descriptor.
 *
 * \param[in, out]  f            pointer to field structure
 * \param[in]       val          pointer to array of values
 * \param[in]       val_pre      pointer to array of previous values, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_field_map_values(cs_field_t   *f,
                    cs_real_t    *val,
                    cs_real_t    *val_pre)
{
  assert(f != NULL);

  if (f->is_owner) {
    BFT_FREE(f->val);
    BFT_FREE(f->val_pre);
    f->is_owner = false;
  }

  f->val = val;

  /* Add previous time step values if necessary */

  if (f->n_time_vals > 1)
    f->val_pre = val_pre;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate boundary condition coefficients arrays.
 *
 * For fields on location CS_MESH_LOCATION_CELLS, boundary conditions
 * are located on CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * Boundary condition coefficients are not currently supported for other
 * locations (though support could be added by mapping a boundary->location
 * indirection array in the cs_mesh_location_t structure).
 *
 * For multidimensional fields, arrays are assumed to have the same
 * interleaving behavior as the field, unless components are coupled.
 *
 * For multidimensional fields with coupled components, interleaving
 * is the norm, and implicit coefficients arrays are arrays of block matrices,
 * not vectors, so the number of entries for each boundary face is
 * dim*dim instead of dim.
 *
 * \param[in, out]  f             pointer to field structure
 * \param[in]       have_flux_bc  if true, flux bc coefficients (af and bf)
 *                                are added
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_bc_coeffs(cs_field_t  *f,
                            bool         have_flux_bc)
{
  /* Add boundary condition coefficients if required */

  cs_lnum_t a_mult = f->dim;
  cs_lnum_t b_mult = f->dim;

  if (f->type & CS_FIELD_VARIABLE) {
    int coupled = 0;
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > -1)
      coupled = cs_field_get_key_int(f, coupled_key_id);
    if (coupled)
      b_mult *= f->dim;
  }

  if (f->location_id == CS_MESH_LOCATION_CELLS) {

    const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(location_id);

    if (f->bc_coeffs == NULL) {

      BFT_MALLOC(f->bc_coeffs, 1, cs_field_bc_coeffs_t);

      f->bc_coeffs->location_id = location_id;

      BFT_MALLOC(f->bc_coeffs->a, n_elts[0]*a_mult, cs_real_t);
      BFT_MALLOC(f->bc_coeffs->b, n_elts[0]*b_mult, cs_real_t);

      if (have_flux_bc) {
        BFT_MALLOC(f->bc_coeffs->af,  n_elts[0]*a_mult, cs_real_t);
        BFT_MALLOC(f->bc_coeffs->bf,  n_elts[0]*b_mult, cs_real_t);
      }
      else {
        f->bc_coeffs->af = NULL;
        f->bc_coeffs->bf = NULL;
      }

    }

    else {

      BFT_REALLOC(f->bc_coeffs->a, n_elts[0]*a_mult, cs_real_t);
      BFT_REALLOC(f->bc_coeffs->b, n_elts[0]*b_mult, cs_real_t);

      if (have_flux_bc) {
        BFT_REALLOC(f->bc_coeffs->af,  n_elts[0]*a_mult, cs_real_t);
        BFT_REALLOC(f->bc_coeffs->bf,  n_elts[0]*b_mult, cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->af);
        BFT_FREE(f->bc_coeffs->bf);
      }

    }

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\"\n"
                " has location %d, which does not support BC coefficients."),
              f->name, f->location_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map existing field boundary condition coefficient arrays.
 *
 * For fields on location CS_MESH_LOCATION_CELLS, boundary conditions
 * are located on CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * Boundary condition coefficients are not currently supported for other
 * locations (though support could be added by mapping a boundary->location
 * indirection array in the cs_mesh_location_t structure).
 *
 * For multidimensional fields, arrays are assumed to have the same
 * interleaving behavior as the field, unless components are coupled.
 *
 * For multidimensional fields with coupled components, interleaving
 * is the norm, and implicit coefficients arrays are arrays of block matrices,
 * not vectors, so the number of entris for each boundary face is
 * dim*dim instead of dim.
 *
 * \param[in, out]  f   pointer to field structure
 * \param[in]       a   explicit BC coefficients array
 * \param[in]       b   implicit BC coefficients array
 * \param[in]       af  explicit flux BC coefficients array, or NULL
 * \param[in]       bf  implicit flux BC coefficients array, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_field_map_bc_coeffs(cs_field_t  *f,
                       cs_real_t   *a,
                       cs_real_t   *b,
                       cs_real_t   *af,
                       cs_real_t   *bf)
{
  /* Add boundary condition coefficients if required */

  if (f->location_id == CS_MESH_LOCATION_CELLS) {

    const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

    if (f->bc_coeffs == NULL) {
      BFT_MALLOC(f->bc_coeffs, 1, cs_field_bc_coeffs_t);
      f->bc_coeffs->location_id = location_id;
    }
    else {
      BFT_FREE(f->bc_coeffs->a);
      BFT_FREE(f->bc_coeffs->b);
      BFT_FREE(f->bc_coeffs->af);
      BFT_FREE(f->bc_coeffs->bf);
    }

    f->bc_coeffs->a = a;
    f->bc_coeffs->b = b;
    f->bc_coeffs->af = af;
    f->bc_coeffs->bf = bf;

  }

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\"\n"
                " has location %d, which does not support BC coefficients."),
              f->name, f->location_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all defined fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_destroy_all(void)
{
  int i;

  for (i = 0; i < _n_fields; i++) {
    cs_field_t  *f = _fields + i;
    if (f->is_owner) {
      BFT_FREE(f->val);
      BFT_FREE(f->val_pre);
    }
    if (f->bc_coeffs != NULL) {
      if (f->is_owner == true) {
        BFT_FREE(f->bc_coeffs->a);
        BFT_FREE(f->bc_coeffs->b);
        BFT_FREE(f->bc_coeffs->af);
        BFT_FREE(f->bc_coeffs->bf);
      }
      BFT_FREE(f->bc_coeffs);
    }
  }

  BFT_FREE(_fields);

  cs_map_name_to_id_destroy(&_field_map);

  _cs_field_free_str();

  BFT_FREE(_key_vals);

  _n_fields = 0;
  _n_fields_max = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate arrays for all defined fields based on their location.
 *
 * Location sized must thus be known.
 *
 * Fields that do not own their data should all have been mapped at this
 * stage, and are checked.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_or_map_all(void)
{
  int i;

  for (i = 0; i < _n_fields; i++) {
    cs_field_t  *f = _fields + i;
    if (f->is_owner)
      cs_field_allocate_values(f);
    else
      if (f->val == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Field \"%s\"\n"
                    " requires mapped values which have not been set."),
                  f->name);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a field based on its id.
 *
 * This function requires that a field of the given id is defined.
 *
 * \param[in]  id   field id
 *
 * \return  pointer to the field structure
 */
/*----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_id(int  id)
{
  if (id > -1 && id < _n_fields)
    return _fields + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Field with id %d is not defined."), id);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a field based on its name.
 *
 * This function requires that a field of the given name is defined.
 *
 * \param[in]  name  field name
 *
 * \return  pointer to the field structure
 */
/*----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_name(const char  *name)
{
  int id = cs_map_name_to_id_try(_field_map, name);

  if (id > -1)
    return _fields + id;
  else {
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\" is not defined."), name);
    return NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a field based on its name if present.
 *
 * If no field of the given name is defined, NULL is returned.
 *
 * \param[in]  name  field name
 *
 * \return  pointer to the field structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_name_try(const char  *name)
{
  int id = cs_map_name_to_id_try(_field_map, name);

  if (id > -1)
    return _fields + id;
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an id associated with a given key name.
 *
 * The key must have been defined previously.
 *
 * \param[in]  name   key name
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_key_id(const char  *name)
{
  int id = -1;

  if (_key_map != NULL)
    id = cs_map_name_to_id_try(_key_map, name);

  if (id < 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\" is not defined."), name);

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return an id associated with a given key name if present.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * \param[in]  name   key name
 *
 * \return  id associated with key, or -1
 */
/*----------------------------------------------------------------------------*/

int
cs_field_key_id_try(const char  *name)
{
  int id = -1;

  if (_key_map != NULL)
    id = cs_map_name_to_id_try(_key_map, name);

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a key for an integer value by its name and return an
 * associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * \param[in]  name            key name
 * \param[in]  default_value   default value associated with key
 * \param[in]  type_flag       mask associated with field types with which the
 *                             key may be associated, or 0
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_key_int(const char  *name,
                        int          default_value,
                        int          type_flag)
{
  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;
  int *def_val = (int *)(kd->def_val);

  *def_val = default_value;
  kd->type_flag = type_flag;
  kd->type_id = 'i';
  kd->is_sub = false;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a key for an floating point value by its name and return an
 * associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * \param[in]  name            key name
 * \param[in]  default_value   default value associated with key
 * \param[in]  type_flag       mask associated with field types with which
 *                             the key may be associated, or 0
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_key_double(const char  *name,
                           double       default_value,
                           int          type_flag)
{
  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;
  double *def_val = (double *)(kd->def_val);

  *def_val = default_value;
  kd->type_flag = type_flag;
  kd->type_id = 'd';
  kd->is_sub = false;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a key for an floating point value by its name and return an
 * associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * \param[in]  name            key name
 * \param[in]  default_value   default value associated with key
 * \param[in]  type_flag       mask associated with field types with which
 *                             the key may be associated, or 0
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_key_str(const char  *name,
                        const char  *default_value,
                        int          type_flag)
{
  int n_keys_init = _n_keys;

  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;

  char **def_val = (char **)(kd->def_val);

  /* Free possible previous allocation */
  if (n_keys_init == _n_keys)
    BFT_FREE(*def_val);

  if (default_value != NULL) {
    BFT_MALLOC(*def_val, strlen(default_value) + 1, char);
    strcpy(*def_val, default_value);
  }
  else
    *def_val = NULL;
  kd->type_flag = type_flag;
  kd->type_id = 's';
  kd->is_sub = false;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a sub key.
 *
 * The sub key is the same type as the parent key.
 *
 * For a given field, when querying a sub key's value and that value has not
 * been set, the query will return the value of the parent key.
 *
 * \param[in]  name            key name
 * \param[in]  parent_id       parent key id
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_sub_key(const char  *name,
                        int          parent_id)
{
  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;
  cs_field_key_def_t *pkd = _key_defs + parent_id;
  int *def_val = (int *)(kd->def_val);

  assert(parent_id > -1 && parent_id < _n_keys);

  *def_val = parent_id;
  kd->type_flag = pkd->type_flag;
  kd->type_id = pkd->type_id;
  kd->is_sub = true;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all defined field keys and associated values.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_destroy_all_keys(void)
{
  _n_keys = 0;
  _n_keys_max = 0;
  BFT_FREE(_key_defs);

  cs_map_name_to_id_destroy(&_key_map);

  BFT_FREE(_key_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the type flag associated with a given key id.
 *
 * If the key has not been defined previously, -1 is returned.
 *
 * \param[in]  key_id  id of associated key
 *
 * \return  type flag associated with key, or -1
 */
/*----------------------------------------------------------------------------*/

int
cs_field_key_flag(int key_id)
{
  int retval = -1;

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    retval = kd->type_flag;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a integer value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  value   value associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_set_key_int(cs_field_t  *f,
                     int          key_id,
                     int          value)
{
  int retval = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      int *_val = (int *)(kv->val);
      *_val = value;
      kv->is_set = true;
    }
  }
  else
    retval = CS_FIELD_INVALID_KEY_ID;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a integer value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  integer value associated with the key id for this field
 */
/*----------------------------------------------------------------------------*/

int
cs_field_get_key_int(const cs_field_t  *f,
                     int                key_id)
{
  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 'i')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      int retval = 0;
      if (kv->is_set)
        retval = *((int *)(kv->val));
      else if (kd->is_sub)
        retval = cs_field_get_key_int(f, *((int *)(kd->def_val)));
      else
        retval = *((int *)(kd->def_val));
      return retval;
    }
  }
  else
    errcode = CS_FIELD_INVALID_KEY_ID;

  if (errcode != CS_FIELD_OK) {
    const char *key = cs_map_name_to_id_reverse(_key_map, key_id);
    if (errcode == CS_FIELD_INVALID_CATEGORY)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" with type flag %d\n"
                  "has no value associated with key %d (\"%s\")."),
                f->name, f->type, key_id, key);
    else if (errcode == CS_FIELD_INVALID_TYPE)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" has keyword %d (\"%s\")\n"
                  "of type \"%c\" and not \"%c\"."),
                f->name, key_id, key, (_key_defs + key_id)->type_id, 'i');
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a floating point value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  value   value associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_set_key_double(cs_field_t  *f,
                        int          key_id,
                        double       value)
{
  int retval = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      double *_val = (double *)(kv->val);
      *_val = value;
      kv->is_set = true;
    }
  }
  else
    retval = CS_FIELD_INVALID_KEY_ID;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a floating point value for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  floating point value associated with the key id for this field
 */
/*----------------------------------------------------------------------------*/

double
cs_field_get_key_double(const cs_field_t  *f,
                        int                key_id)
{
  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 'd')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      double retval = 0.;
      if (kv->is_set)
        retval = *((double *)(kv->val));
      else if (kd->is_sub)
        retval = cs_field_get_key_double(f, *((int *)(kd->def_val)));
      else
        retval = *((double *)(kd->def_val));
      return retval;
    }
  }
  else
    errcode = CS_FIELD_INVALID_KEY_ID;

  if (errcode != CS_FIELD_OK) {
    const char *key = cs_map_name_to_id_reverse(_key_map, key_id);
    if (errcode == CS_FIELD_INVALID_CATEGORY)
      bft_error(__FILE__, __LINE__, 0,
                _("Field %s with type flag %d\n"
                  "has no value associated with key %d (%s)."),
                f->name, f->type, key_id, key);
    else if (errcode == CS_FIELD_INVALID_TYPE)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" has keyword %d (\"%s\")\n"
                  "of type \"%c\" and not \"%c\"."),
                f->name, key_id, key, (_key_defs + key_id)->type_id, 'd');
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return 0.;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a character string for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  str     string associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_set_key_str(cs_field_t  *f,
                     int          key_id,
                     const char  *str)
{
  int retval = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      char **_val = (char **)(kv->val);
      if (kv->is_set == false)
        *_val = NULL;
      BFT_REALLOC(*_val, strlen(str) + 1, char);
      strcpy(*_val, str);
      kv->is_set = true;
    }
  }
  else
    retval = CS_FIELD_INVALID_KEY_ID;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a string for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  pointer to character string associated with
 *          the key id for this field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_field_get_key_str(const cs_field_t  *f,
                     int                key_id)
{
  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 's')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const char *str = NULL;
      if (kv->is_set)
        str = *((const char **)(kv->val));
      else if (kd->is_sub)
        str = cs_field_get_key_str(f, *((int *)(kd->def_val)));
      else
        str = *((const char **)(kd->def_val));
      return str;
    }
  }
  else
    errcode = CS_FIELD_INVALID_KEY_ID;

  if (errcode != CS_FIELD_OK) {
    const char *key = cs_map_name_to_id_reverse(_key_map, key_id);
    if (errcode == CS_FIELD_INVALID_CATEGORY)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" with type flag %d\n"
                  "has no value associated with key %d (\"%s\")."),
                f->name, f->type, key_id, key);
    else if (errcode == CS_FIELD_INVALID_TYPE)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" has keyword %d (\"%s\")\n"
                  "of type \"%c\" and not \"%c\"."),
                f->name, key_id, key, (_key_defs + key_id)->type_id, 'i');
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all field definitions to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_defs(void)
{
  int i, j, cat_id;

  int n_cat_fields = 0;

  int mask_id_start = 2; /* _type_flag_*[CS_FIELD_VARIABLE] */
  int mask_id_end = 6;   /* _type_flag_*[CS_FIELD_USER] */
  int mask_prev = 0;

  if (_n_fields == 0)
    return;

  /* Fields by category */

  for (cat_id = mask_id_start; cat_id < mask_id_end + 1; cat_id++) {

    n_cat_fields = 0;

    for (i = 0; i < _n_fields; i++) {

      char ilv_c = ' ';

      const cs_field_t *f = _fields + i;

      if (f->type & mask_prev)
        continue;

      if (cat_id == mask_id_end || f->type & _type_flag_mask[cat_id]) {

        /* Print header for first field of each category */

        if (n_cat_fields == 0) {

          char tmp_s[4][64] =  {"", "", "", ""};

          cs_log_strpad(tmp_s[0], _("Field"), 24, 64);
          cs_log_strpad(tmp_s[1], _("Dim."), 4, 64);
          cs_log_strpad(tmp_s[2], _("Location"), 20, 64);
          cs_log_strpad(tmp_s[3], _("Id"), 4, 64);

          /* Print logging header */

          if (cat_id < mask_id_end)
            cs_log_printf(CS_LOG_SETUP,
                          _("\n"
                            "Fields of type: %s\n"
                            "---------------\n"), _(_type_flag_name[cat_id]));
          else
            cs_log_printf(CS_LOG_SETUP,
                          _("\n"
                            "Other fields:\n"
                            "-------------\n"));
          cs_log_printf(CS_LOG_SETUP, "\n");

          cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s Type flag\n"),
                        tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

          for (j = 0; j < 4; j++)
            memset(tmp_s[j], '-', 64);

          tmp_s[0][24] = '\0';
          tmp_s[1][4] = '\0';
          tmp_s[2][20] = '\0';
          tmp_s[3][4] = '\0';

          cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s ---------\n"),
                        tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

        }

        /* Print field info */

        if (f->interleaved == false)
          ilv_c = 'n';

        cs_log_printf(CS_LOG_SETUP,
                      "  %-24s %d %c  %-20s %-4d ",
                      f->name, f->dim, ilv_c,
                      _(cs_mesh_location_get_name(f->location_id)),
                      f->id);

        if (f->type != 0) {
          cs_log_printf(CS_LOG_SETUP, "%-4d", f->type);
          _log_add_type_flag(f->type);
          cs_log_printf(CS_LOG_SETUP, "\n");
        }
        else
          cs_log_printf(CS_LOG_SETUP, "0\n");

        n_cat_fields++;

      }

    } /* End of loop on fields */

    if (cat_id < mask_id_end)
      mask_prev += _type_flag_mask[cat_id];

  } /* End fo loop on categories */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given field to log file.
 *
 * \param[in]  f             pointer to field structure
 * \param[in]  log_keywords  log level for keywords (0: do not log,
 *                           1: log non-default values, 2: log all)
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_info(const cs_field_t  *f,
                  int                log_keywords)
{
  if (f == NULL)
    return;

  /* Global indicators */
  /*-------------------*/

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "  Field: \"%s\"\n"), f->name);

  if (log_keywords > 0)
    cs_log_printf(CS_LOG_SETUP, "\n");

  cs_log_printf(CS_LOG_SETUP,
                _("    Id:                         %d\n"
                  "    Type:                       %d"),
                f->id, f->type);

  if (f->type != 0) {
    int i;
    int n_loc_flags = 0;
    for (i = 0; i < _n_type_flags; i++) {
      if (f->type & _type_flag_mask[i]) {
        if (n_loc_flags == 0)
          cs_log_printf(CS_LOG_SETUP, " (%s", _(_type_flag_name[i]));
        else
          cs_log_printf(CS_LOG_SETUP, ", %s", _(_type_flag_name[i]));
        n_loc_flags++;
      }
    }
    if (n_loc_flags > 0)
      cs_log_printf(CS_LOG_SETUP, ")");
    cs_log_printf(CS_LOG_SETUP, "\n");
  }

  cs_log_printf(CS_LOG_SETUP, _("    Location:                   %s\n"),
                cs_mesh_location_get_name(f->location_id));

  if (f->dim == 1)
    cs_log_printf(CS_LOG_SETUP, _("    Dimension:                  1\n"));
  else if (f->interleaved == false)
    cs_log_printf(CS_LOG_SETUP,
                  _("    Dimension:                  %d (non-interleaved)\n"),
                  f->dim);
  else
    cs_log_printf(CS_LOG_SETUP,
                  _("    Dimension:                  %d (interleaved)\n"),
                  f->dim);

  if (f->is_owner == false)
    cs_log_printf(CS_LOG_SETUP,
                  _("    Values mapped from external definition\n"));

  if (_n_keys > 0 && log_keywords > 0) {
    int i;
    const char null_str[] = "(null)";
    cs_log_printf(CS_LOG_SETUP, _("\n    Associated key values:\n"));
    for (i = 0; i < _n_keys; i++) {
      int key_id = cs_map_name_to_id_try(_key_map,
                                         cs_map_name_to_id_key(_key_map, i));
      cs_field_key_def_t *kd = _key_defs + key_id;
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const char *key = cs_map_name_to_id_key(_key_map, i);
      if (kd->type_flag == 0 || (kd->type_flag & f->type)) {
        if (kd->type_id == 'i') {
          if (kv->is_set == true)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10d\n"),
                          key, *((int *)kv->val));
          else if (log_keywords > 1)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10d (default)\n"),
                          key, *((int *)kd->def_val));
        }
        else if (kd->type_id == 'd') {
          if (kv->is_set == true)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10.3g\n"),
                          key, *((double *)kv->val));
          else if (log_keywords > 1)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10.3g (default)\n"),
                          key, *((double *)kd->def_val));
        }
        else if (kd->type_id == 's') {
          const char *s;
          if (kv->is_set == true) {
            s = *((const char **)(kv->val));
            if (s == NULL)
              s = null_str;
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %s\n"), key, s);
          }
          else if (log_keywords > 1) {
            s = *(const char **)(kd->def_val);
            if (s == NULL)
              s = null_str;
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10s (default)\n"),
                          key, s);
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all defined fields to log file.
 *
 * \param[in]  log_keywords  log level for keywords (0: do not log,
 *                           1: log non-default values, 2: log all)
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_fields(int  log_keywords)
{
  int i, cat_id;
  const cs_field_t  *f;

  int n_cat_fields = 0;

  int mask_id_start = 2; /* _type_flag_*[CS_FIELD_VARIABLE] */
  int mask_id_end = 6;   /* _type_flag_*[CS_FIELD_USER] */
  int mask_prev = 0;

  if (_n_fields == 0)
    return;

  /* Fields by category */

  for (cat_id = mask_id_start; cat_id < mask_id_end + 1; cat_id++) {

    n_cat_fields = 0;

    for (i = 0; i < _n_fields; i++) {

      f = _fields + i;

      if (f->type & mask_prev)
        continue;

      if (cat_id == mask_id_end || f->type & _type_flag_mask[cat_id]) {

        if (n_cat_fields == 0) {
          if (cat_id < mask_id_end)
            cs_log_printf(CS_LOG_SETUP,
                          _("\n"
                            "Fields of type: %s\n"
                            "---------------\n"), _(_type_flag_name[cat_id]));
          else
            cs_log_printf(CS_LOG_SETUP,
                          _("\n"
                            "Other fields:\n"
                            "-------------\n"));
        }
        cs_field_log_info(f, log_keywords);
        n_cat_fields++;

      }

    } /* End of loop on fields */

    if (cat_id < mask_id_end)
      mask_prev += _type_flag_mask[cat_id];

  } /* End of loop on categories */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all key definitions to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_key_defs(void)
{
  int i;
  char tmp_s[4][64] =  {"", "", "", ""};

  if (_n_keys == 0)
    return;

  /* Print logging header */

  cs_log_strpad(tmp_s[0], _("Field"), 24, 64);
  cs_log_strpad(tmp_s[1], _("Default"), 12, 64);
  cs_log_strpad(tmp_s[2], _("Type"), 7, 64);
  cs_log_strpad(tmp_s[3], _("Id"), 4, 64);

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Defined field keys:\n"
                  "-------------------\n\n"));
  cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s Type flag\n"),
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  for (i = 0; i < 24; i++)
    tmp_s[0][i] = '-';
  tmp_s[0][24] = '\0';
  for (i = 0; i < 12; i++)
    tmp_s[1][i] = '-';
  tmp_s[1][12] = '\0';
  for (i = 0; i < 7; i++)
    tmp_s[2][i] = '-';
  tmp_s[2][7] = '\0';
  for (i = 0; i < 4; i++)
    tmp_s[3][i] = '-';
  tmp_s[3][4] = '\0';

  cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s ---------\n"),
                tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  for (i = 0; i < _n_keys; i++) {

    int key_id = cs_map_name_to_id_try(_key_map,
                                       cs_map_name_to_id_key(_key_map, i));
    cs_field_key_def_t *kd = _key_defs + key_id;
    const char *key = cs_map_name_to_id_key(_key_map, i);

    if (kd->type_id == 'i') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12d integer %-4d "),
                    key, *((int *)kd->def_val), key_id);
    }
    else if (kd->type_id == 'd') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12.3g real    %-4d "),
                    key, *((double *)kd->def_val), key_id);
    }
    else if (kd->type_id == 's') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12s string  %-4d "),
                    key, (*((const char **)kd->def_val)), key_id);
    }
    if (kd->type_flag == 0)
      cs_log_printf(CS_LOG_SETUP, "0\n");
    else {
      cs_log_printf(CS_LOG_SETUP, "%-4d", kd->type_flag);
      _log_add_type_flag(kd->type_flag);
      cs_log_printf(CS_LOG_SETUP, "\n");
    }

  } /* End of loop on keys */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given field key to log file.
 *
 * \param[in]  key_id        id of associated key
 * \param[in]  log_defaults  if true, log default field values in addition to
 *                           defined field values
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_key_vals(int   key_id,
                      bool  log_defaults)
{
  int i, cat_id;
  cs_field_key_def_t *kd;

  int mask_id_start = 2; /* _type_flag_*[CS_FIELD_VARIABLE] */
  int mask_id_end = 6;   /* _type_flag_*[CS_FIELD_USER] */
  int mask_prev = 0;
  const char null_str[] = "(null)";


  if (key_id < 0 || key_id >= _n_keys)
    return;

  kd = _key_defs + key_id;

  /* Global indicators */
  /*-------------------*/

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "  Key: \"%s\", values per field\n"
                  "  ----\n"),
                cs_map_name_to_id_reverse(_key_map, key_id));

  /* Loop on categories, building a mask for previous categories
     so as not to output data twice */

  for (cat_id = mask_id_start; cat_id < mask_id_end + 1; cat_id++) {

    for (i = 0; i < _n_fields; i++) {

      const cs_field_t *f = _fields + i;

      if (f->type & mask_prev)
        continue;

      if (cat_id == mask_id_end || f->type & _type_flag_mask[cat_id]) {

        cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);

        if (kd->type_flag == 0 || (kd->type_flag & f->type)) {
          if (kd->type_id == 'i') {
            if (kv->is_set == true)
              cs_log_printf(CS_LOG_SETUP, "    %-24s %d\n",
                            f->name, *((int *)kv->val));
            else if (log_defaults)
              cs_log_printf(CS_LOG_SETUP, _("    %-24s %-10d (default)\n"),
                            f->name, *((int *)kd->def_val));
          }
          else if (kd->type_id == 'd') {
            if (kv->is_set == true)
              cs_log_printf(CS_LOG_SETUP, _("    %-24s %-10.3g\n"),
                          f->name, *((double *)kv->val));
            else if (log_defaults)
              cs_log_printf(CS_LOG_SETUP, _("    %-24s %-10.3g (default)\n"),
                            f->name, *((double *)kd->def_val));
          }
          else if (kd->type_id == 's') {
            const char *s;
            if (kv->is_set == true) {
              s = *((const char **)(kv->val));
              if (s == NULL)
                s = null_str;
              cs_log_printf(CS_LOG_SETUP, _("    %-24s %s\n"), f->name, s);
            }
            else if (log_defaults) {
              s = *(const char **)(kd->def_val);
              if (s == NULL)
                s = null_str;
              cs_log_printf(CS_LOG_SETUP, _("    %-24s %-10s (default)\n"),
                            f->name, s);
            }
          }
        }
      }

    }

    if (cat_id < mask_id_end)
      mask_prev += _type_flag_mask[cat_id];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all given field keys to log file.
 *
 * \param[in]  log_defaults  if true, log default field values in addition to
 *                           defined field values
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_all_key_vals(bool  log_defaults)
{
  int i;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Defined key values per field:\n"
                  "-----------------------------\n\n"));

  for (i = 0; i < _n_keys; i++)
    cs_field_log_key_vals(i, log_defaults);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define base keys.
 *
 * Keys defined by this function are:
 *   "label"     (string)
 *   "post_vis"  (integer)
 *   "coupled"   (integer, restricted to CS_FIELD_VARIABLE)
 *   "moment_dt" (integer, restricted to CS_FIELD_PROPERTY);
 *
 * A recommened practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_define_keys_base(void)
{
  cs_field_define_key_str("label", NULL, 0);

  cs_field_define_key_int("post_vis", 0, 0);
  cs_field_define_key_int("coupled", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("moment_dt", -1, CS_FIELD_PROPERTY);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
