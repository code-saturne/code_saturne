/*============================================================================
 * Field management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

  unsigned char      def_val[8];   /* Default value container (int or double) */
  int                type_flag;    /* Type flag */
  char               type_id;      /* i: int; d: double */

} cs_field_key_def_t;

/* Field key value structures */

typedef struct {

  unsigned char      val[8];       /* Value container (int or double) */
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
                                      CS_FIELD_USER};
static const char *_type_flag_name[] = {N_("intensive"),
                                        N_("extensive"),
                                        N_("variable"),
                                        N_("property"),
                                        N_("postprocess"),
                                        N_("user")};

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
 *
 * returns:
 *   pointer to new field.
 *
 *----------------------------------------------------------------------------*/

static cs_field_t *
_field_create(const char   *name,
              int           type_flag,
              int           location_id,
              int           dim)
{
  int key_id;
  int field_id = -1;
  cs_field_t *f =  NULL;

  /* Initialize if necessary */

  if (_field_map == NULL)
    _field_map = cs_map_name_to_id_create();

  /* Find or insert entry in map */

  field_id = cs_map_name_to_id(_field_map, name);

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

  strncpy(f->label, f->name, 64);
  f->id = field_id;
  f->type = type_flag;
  f->dim = dim;
  f->interleaved = true;
  f->location_id = location_id;

  f->val = NULL;
  f->val_pre = NULL;
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
 *   n_elts <-- number of associated elements
 *   dim    <-- associated dimension
 *
 * returns  pointer to new field values.
 *----------------------------------------------------------------------------*/

static cs_real_t *
_add_val(cs_lnum_t  n_elts,
         int        dim)
{
  cs_lnum_t  ii;
  cs_real_t  *val;

  BFT_MALLOC(val, n_elts*dim, cs_real_t);

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
 */
/*----------------------------------------------------------------------------*/
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

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Map an existing field.
 *
 * Fortran interface
 *
 * subroutine fldmap (name,   lname,  ireawr, numsui, ierror)
 * *****************
 *
 * character*       name        : <-- : Field name
 * integer          lname       : <-- : Field name length
 * integer          iexten      : <-- : 0: intensive; 1: extensive
 * integer          itycat      : <-- : Field category
 *                              :     :  0: variable
 *                              :     :  1: property
 * integer          ityloc      : <-- : Location type
 *                              :     :  0: cells
 *                              :     :  1: interior faces
 *                              :     :  2: interior faces
 *                              :     :  3: vertices
 * integer          idim        : <-- : Field dimension
 * integer          ilved       : <-- : 0: not intereaved; 1: interleaved
 * integer          iprev       : <-- : 0: no previous values, 1: previous
 * cs_real_t*       val         : <-- : Pointer to field values array
 * cs_real_t*       valp        : <-- : Pointer to values at previous
 *                              :     : time step if iprev = 1
 * integer          idfld       : --> : id of mapped field
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldmap, FLDMAP)
(
 const char       *name,
 const cs_int_t   *lname,
 const cs_int_t   *iexten,
 const cs_int_t   *itycat,
 const cs_int_t   *ityloc,
 const cs_int_t   *idim,
 const cs_int_t   *ilved,
 const cs_int_t   *iprev,
 cs_real_t        *val,
 cs_real_t        *valp,
 cs_int_t         *idfld
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufname;
  int type_flag = 0;
  bool interleaved = (*ilved == 0) ? false : true;
  cs_real_t *_valp = NULL;

  cs_field_t *f = NULL;

  bufname = cs_base_string_f_to_c_create(name, *lname);

  if (*iexten)
    type_flag = CS_FIELD_EXTENSIVE;
  else
    type_flag = CS_FIELD_INTENSIVE;

  if (*itycat == 0) {
    type_flag = type_flag | CS_FIELD_VARIABLE;
    _valp = valp;
  }
  else if (*itycat == 1) {
    type_flag = type_flag | CS_FIELD_PROPERTY;
    _valp = valp;
  }
  else if (*itycat == 2)
    type_flag = type_flag | CS_FIELD_POSTPROCESS;
  else if (*itycat == 3)
    type_flag = type_flag | CS_FIELD_USER;

  if (*iprev)
    _valp = valp;

  f = cs_field_create_mapped(bufname,
                             type_flag,
                             *ityloc,
                             *idim,
                             interleaved,
                             val,
                             _valp);

  cs_base_string_f_to_c_free(&bufname);

  *idfld = f->id;
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

void CS_PROCF (fldfid, FLDFID)
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
 * Fortran interface
 *
 * subroutine fldkid (name,   lname,  ikeyid)
 * *****************
 *
 * character*       name        : <-- : key
 * integer          lname       : <-- : Key name length
 * integer          ikeyid      : --> : id of given key
 *----------------------------------------------------------------------------*/

void CS_PROCF (fldkid, FLDKID)
(
 const char       *name,
 const cs_int_t   *lname,
 cs_int_t         *ikeyid
 CS_ARGF_SUPP_CHAINE              /*   (possible 'length' arguments added
                                        by many Fortran compilers) */
)
{
  char *bufname;

  bufname = cs_base_string_f_to_c_create(name, *lname);

  *ikeyid = cs_field_key_id_try(bufname);

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
                bool          has_previous)
{
  cs_field_t  *f =  _field_create(name,
                                  type_flag,
                                  location_id,
                                  dim);

  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(location_id);

  f->val = _add_val(n_elts[2], dim);

  /* Add previous time step values if necessary */

  if (has_previous)
    f->val_pre = _add_val(n_elts[2], dim);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a field descriptor for an existing array of values.
 *
 * \param[in]  name         field name
 * \param[in]  type_flag    mask of field property and category values
 * \param[in]  location_id  id of associated location
 * \param[in]  dim          field dimension (number of components)
 * \param[in]  interleaved  true if field components are interleaved
 * \param[in]  val          pointer to array of values
 * \param[in]  val_pre      pointer to array of previous values, or NULL
 *
 * \return  pointer to new field.
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_field_create_mapped(const char   *name,
                       int           type_flag,
                       int           location_id,
                       int           dim,
                       bool          interleaved,
                       double       *val,
                       double       *val_pre)
{
  cs_field_t  *f =  _field_create(name,
                                  type_flag,
                                  location_id,
                                  dim);

  f->interleaved = interleaved;
  f->val = val;
  f->is_owner = false;

  if (val_pre != NULL)
    f->val_pre = val_pre;

  return f;
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
  }

  _n_fields = 0;
  _n_fields_max = 0;
  BFT_FREE(_fields);

  cs_map_name_to_id_destroy(&_field_map);

  BFT_FREE(_key_vals);
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
 * \brief Define a key for an integer value by its name and return an associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * \param[in]  name            key name
 * \param[in]  default_value   default value associated with key
 * \param[in]  type flag       mask associated with field types with which the
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
 * \param[in]  type flag       mask associated with field types with which
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
 * \param[in]  value   value associated with key
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
      int *_val = (int *)(kv->val);
      return *_val;
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
 * \param[in]  value   value associated with key
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
      double *_val = (double *)(kv->val);
      return *_val;
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

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to a given field to log file.
 *
 * \param[in]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_info(const cs_field_t  *f)
{
  if (f == NULL)
    return;

  /* Global indicators */
  /*-------------------*/

  bft_printf(_("\n"
               "  Field:\"%s\"\n"), f->name);

  bft_printf(_("\n"
               "    Label:                      %s\n"
               "    Id:                         %d\n"
               "    Type:                       %d"),
             f->label, f->id, f->type);

  if (f->type != 0) {
    int i;
    int n_loc_flags = 0;
    for (i = 0; i < _n_type_flags; i++) {
      if (f->type & _type_flag_mask[i]) {
        if (n_loc_flags == 0)
          bft_printf(_(" (%s"), _(_type_flag_name[i]));
        else
          bft_printf(_(", %s"), _(_type_flag_name[i]));
        n_loc_flags++;
      }
    }
    if (n_loc_flags > 0)
      bft_printf(_(")"));
    bft_printf(_("\n"));
  }

  if (f->dim == 1)
    bft_printf(_("    Dimension:                  1\n"));
  else if (f->interleaved == false)
    bft_printf(_("    Dimension:                  %d (non-interleaved)\n"),
               f->dim);
  else
    bft_printf(_("    Dimension:                  %d (interleaved)\n"),
               f->dim);

  if (f->is_owner == false)
    bft_printf(_("    Values mapped from external definition\n"));

  if (_n_keys > 0) {
    int i;
    bft_printf(_("\n    Associated key values:\n"));
    for (i = 0; i < _n_keys; i++) {
      int key_id = cs_map_name_to_id_try(_key_map,
                                         cs_map_name_to_id_key(_key_map, i));
      cs_field_key_def_t *kd = _key_defs + key_id;
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const char *key = cs_map_name_to_id_key(_key_map, i);
      if (kd->type_flag == 0 || (kd->type_flag & f->type)) {
        if (kd->type_id == 'i') {
          if (kv->is_set == true)
            bft_printf(_("      %-24s %10d\n"), key, *((int *)kv->val));
          else
            bft_printf(_("      %-24s %10d (default)\n"),
                       key, *((int *)kd->def_val));
        }
        else if (kd->type_id == 'd') {
          if (kv->is_set == true)
            bft_printf(_("      %-24s %10.3g\n"), key, *((double *)kv->val));
          else
            bft_printf(_("      %-24s %10.3g (default)\n"),
                       key, *((double *)kd->def_val));
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all defined fields to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_fields(void)
{
  int i, cat_id;
  const cs_field_t  *f;

  int n_cat_fields = 0;
  int *field_id = NULL;

  int mask_id_start = 2; /* _type_flag_*[CS_FIELD_VARIABLE] */
  int mask_id_end = 6;   /* _type_flag_*[CS_FIELD_USER] */

  if (_n_fields == 0)
    return;

  BFT_MALLOC(field_id, _n_fields, int);

  for (i = 0; i < _n_fields; i++)
    field_id[i] = cs_map_name_to_id_try(_field_map,
                                        cs_map_name_to_id_key(_field_map, i));

  /* Fields by category */

  for (cat_id = mask_id_start; cat_id < mask_id_end; cat_id++) {
    n_cat_fields = 0;
    for (i = 0; i < _n_fields; i++) {
      if (field_id[i] < 0)
        continue;
      f = _fields + field_id[i];
      if (f->type & _type_flag_mask[cat_id]) {
        if (n_cat_fields == 0)
          bft_printf(_("\n"
                       "Fields of type: %s\n"
                       "---------------\n"), _type_flag_name[cat_id]);
        cs_field_log_info(f);
        n_cat_fields++;
        field_id[i] = -1;
      }
    }
  }

  /* Other fields if present */

  n_cat_fields = 0;
  for (i = 0; i < _n_fields; i++) {
    if (field_id[i] < 0)
      continue;
    f = _fields + field_id[i];
    if (f->type & _type_flag_mask[cat_id]) {
      if (n_cat_fields == 0)
        bft_printf(_("\n"
                     "Other fields:\n"
                     "-------------\n"));
      cs_field_log_info(f);
      n_cat_fields++;
      field_id[i] = -1;
    }
  }

  BFT_FREE(field_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print info relative to all defined field keys to log file.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_log_keys(void)
{
  int i;
  char tmp_s[4][64] =  {"", "", "", ""};

  if (_n_keys == 0)
    return;

  cs_log_strpad(tmp_s[0], _("Key word"), 24, 64);
  cs_log_strpad(tmp_s[1], _("Default"), 12, 64);
  cs_log_strpad(tmp_s[2], _("Type"), 4, 64);
  cs_log_strpad(tmp_s[3], _("Id"), 4, 64);

  bft_printf(_("\n"
               "Defined field keys:\n"
               "-------------------\n\n"));
  bft_printf(_("  %s %s %s    %s Type flag\n"),
             tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  for (i = 0; i < 24; i++)
    tmp_s[0][i] = '-';
  for (i = 0; i < 12; i++)
    tmp_s[1][i] = '-';
  for (i = 0; i < 4; i++)
    tmp_s[2][i] = '-';
  for (i = 0; i < 4; i++)
    tmp_s[3][i] = '-';

  bft_printf(_("  %s %s %s    %s ---------\n"),
             tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

  for (i = 0; i < _n_keys; i++) {
    int key_id = cs_map_name_to_id_try(_key_map,
                                       cs_map_name_to_id_key(_key_map, i));
    cs_field_key_def_t *kd = _key_defs + key_id;
    const char *key = cs_map_name_to_id_key(_key_map, i);
    if (kd->type_id == 'i') {
      bft_printf
        (_("  %-24s %-12d integer %-4d %-4d"),
         key, *((int *)kd->def_val), key_id, kd->type_flag);
    }
    else if (kd->type_id == 'd') {
      bft_printf
        (_("  %-24s %-12.3g real    %-4d %-4d"),
         key, *((double *)kd->def_val), key_id, kd->type_flag);
    }
    if (kd->type_flag != 0) {
      int j;
      int n_loc_flags = 0;
      for (j = 0; j < _n_type_flags; j++) {
        if (kd->type_flag & _type_flag_mask[j]) {
          if (n_loc_flags == 0)
            bft_printf(_(" (%s"), _(_type_flag_name[j]));
          else
            bft_printf(_(", %s"), _(_type_flag_name[j]));
          n_loc_flags++;
        }
      }
      if (n_loc_flags > 0)
        bft_printf(_(")"));
    }
    bft_printf(_("\n"));
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
