/*============================================================================
 * Field management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field.c
        Field management.

  \struct cs_field_bc_coeffs_t

  \brief Field boundary condition descriptor (for variables)

  \var cs_field_bc_coeffs_t::location_id
       Id of matching location

  \var cs_field_bc_coeffs_t::icodcl
       Low-level BC type code
  \var cs_field_bc_coeffs_t::rcodcl1
       1st part of low-level BC values definition
  \var cs_field_bc_coeffs_t::rcodcl2
       2nd part of low-level BC values definition
  \var cs_field_bc_coeffs_t::rcodcl3
       3rd part of low-level BC values definition

  \var cs_field_bc_coeffs_t::a
       Explicit coefficient
  \var cs_field_bc_coeffs_t::b
       Implicit coefficient
  \var cs_field_bc_coeffs_t::af
       Explicit coefficient for flux
  \var cs_field_bc_coeffs_t::bf
       Implicit coefficient for flux
  \var cs_field_bc_coeffs_t::ad
       Explicit coefficient for divergence
  \var cs_field_bc_coeffs_t::bd
       Implicit coefficient for divergence
  \var cs_field_bc_coeffs_t::ac
       Explicit coefficient for convection
  \var cs_field_bc_coeffs_t::bc
       Implicit coefficient for convection

  \struct cs_field_t

  \brief Field descriptor

  Members of this field are publicly accessible, to allow for concise
  syntax, as it is expected to be used in many places.

  \var  cs_field_t::name
        Canonical name
  \var  cs_field_t::id
        Field id (based on order of field declaration, starting at 0)
  \var  cs_field_t::type
        Field type flag (sum of field mask constants, defining if a field
        is a variable, a property, ...)
  \var  cs_field_t::dim
        Field dimension (usually 1 for scalar, 3 for vector, or 6 for
        symmetric tensor)
  \var  cs_field_t::location_id
        Id of matching mesh location
  \var  cs_field_t::n_time_vals
        Number of time values
  \var  cs_field_t::vals
        vals[0][:] is a pointer to val
        vals[1][:] is a pointer to val_pre
        vals[p][:] is a pointer to p ith previous field values
  \var  cs_field_t::val
        For each active location, pointer to matching values array
  \var  cs_field_t::val_pre
        For each active location, pointer to matching previous values array
        (only if n_time_vals > 1)
  \var  cs_field_t::bc_coeffs
        Boundary condition coefficients, for variable type fields
  \var  cs_field_t::is_owner
        Ownership flag for values
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* Field descriptor allocation block size */

#define _CS_FIELD_S_ALLOC_SIZE       16

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Field key definition values */

union _key_val_t {
  int                         v_int;
  double                      v_double;
  void                       *v_p;
};

/* Field key definitions */

typedef struct {

  union _key_val_t              def_val;        /* Default value container (int,
                                                   double, or pointer to string
                                                   or structure) */
  cs_field_log_key_struct_t    *log_func;       /* print function for
                                                   structure */

  cs_field_log_key_struct_t    *log_func_default; /* default values log
                                                     function for structure */

  cs_field_clear_key_struct_t  *clear_func;     /* memory free function
                                                   for sub-structures */

  size_t                      type_size;        /* Type length for added types
                                                   (0 for 'i', 'd', or 's') */
  int                         type_flag;        /* Field type flag */
  char                        type_id;          /* i: int; d: double; s: str;
                                                   t: type */
  char                        log_id;           /* s: setup; n: none */

  bool                        is_sub;           /* Indicate if the key is a
                                                   sub-key (in which case
                                                   def_val contains the parent
                                                   key id */

} cs_field_key_def_t;

/* Field key value structures */

typedef struct {

  union _key_val_t   val;          /* Value container (int, double,
                                      or pointer) */
  char               is_set;       /* Has this key been set for the
                                      present field ? */
  char               is_locked;    /* Has this key been locked for the
                                      present field ? */

} cs_field_key_val_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Field definitions */

static int  _n_fields = 0;
static int  _n_fields_max = 0;
static cs_field_t  **_fields = NULL;
static cs_map_name_to_id_t  *_field_map = NULL;

/* Key definitions */

static int  _n_keys = 0;
static int  _n_keys_max = 0;
static cs_field_key_def_t  *_key_defs = NULL;
static cs_map_name_to_id_t  *_key_map = NULL;

static int _k_label = -1;

/* Key values : _key_vals[field_id*_n_keys_max + key_id] */

static cs_field_key_val_t  *_key_vals = NULL;

/* Names for logging */

static const int _n_type_flags = 8;
static const int _type_flag_mask[] = {CS_FIELD_INTENSIVE,
                                      CS_FIELD_EXTENSIVE,
                                      CS_FIELD_VARIABLE,
                                      CS_FIELD_PROPERTY,
                                      CS_FIELD_POSTPROCESS,
                                      CS_FIELD_ACCUMULATOR,
                                      CS_FIELD_USER,
                                      CS_FIELD_CDO};
static const char *_type_flag_name[] = {N_("intensive"),
                                        N_("extensive"),
                                        N_("variable"),
                                        N_("property"),
                                        N_("postprocess"),
                                        N_("accumulator"),
                                        N_("user"),
                                        N_("CDO")};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Names for components */

/*! \var field name extension for vector components */
const char *cs_glob_field_comp_name_3[] = {"[X]", "[Y]", "[Z]"};

/*! \var field name extension for symmetric tensor components */
const char *cs_glob_field_comp_name_6[] = {"[XX]", "[YY]", "[ZZ]",
                                           "[XY]", "[YZ]", "[XZ]"};

/*! \var field name extension for tensor components */
const char *cs_glob_field_comp_name_9[] = {"[XX]", "[XY]", "[XZ]",
                                           "[YX]", "[YY]", "[YZ]",
                                           "[ZX]", "[ZY]", "[ZZ]"};

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

int
cs_f_field_n_fields(void);

int
cs_f_field_id_by_name(const char *name);

int
cs_f_field_location(const cs_field_t *f);

int
cs_f_field_id_by_name_try(const char *name);

void
cs_f_field_get_name(int           id,
                    int           name_max,
                    const char  **name,
                    int          *name_len);

void
cs_f_field_get_dimension(int           id,
                         int           dim[1]);

void
cs_f_field_get_type(int           id,
                    int          *type);

int
cs_f_field_have_previous(int   id);

void
cs_f_field_set_n_previous(int  id,
                          int  n_previous);

void
cs_f_field_get_n_previous(int  id,
                          int  n_previous[1]);

void
cs_f_field_var_ptr_by_id(int          id,
                         int          pointer_type,
                         int          pointer_rank,
                         int          dim[2],
                         cs_real_t  **p);

void
cs_f_field_var_ptr_by_id_try(int          id,
                             int          pointer_type,
                             int          pointer_rank,
                             int          dim[2],
                             cs_real_t  **p);

void
cs_f_field_bc_coeffs_ptr_by_id(int          id,
                               int          pointer_type,
                               int          pointer_rank,
                               int          dim[3],
                               cs_real_t  **p);

void
cs_f_field_set_key_int(int  f_id,
                       int  k_id,
                       int  value);

void
cs_f_field_set_key_int_bits(int  f_id,
                            int  k_id,
                            int  mask);

void
cs_f_field_clear_key_int_bits(int  f_id,
                              int  k_id,
                              int  mask);

void
cs_f_field_set_key_double(int     f_id,
                          int     k_id,
                          double  value);

void
cs_f_field_set_key_str(int          f_id,
                       int          k_id,
                       const char  *str);

void
cs_f_field_get_key_str(int           f_id,
                       int           key_id,
                       int           str_max,
                       const char  **str,
                       int          *str_len);

void
cs_f_field_set_key_struct(int    f_id,
                          int    k_id,
                          void  *k_value);

void
cs_f_field_get_key_struct(int    f_id,
                          int    k_id,
                          void  *k_value);

void
cs_f_field_get_label(int           f_id,
                     int           str_max,
                     const char  **str,
                     int          *str_len);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a field descriptor.
 *
 * parameters:
 *   name        <-- field name
 *   type_flag   <-- mask of field property and category values
 *   location_id <-- id of associated location
 *   dim         <-- field dimension (number of components)
 *
 * returns:
 *   pointer to new field.
 *----------------------------------------------------------------------------*/

static cs_field_t *
_field_create(const char   *name,
              int           type_flag,
              int           location_id,
              int           dim)
{
  int key_id;
  int field_id = -1;
  size_t l = strlen(name);

  cs_field_t *f = cs_field_by_name_try(name);

  /* Check this name was not already used */

  if (f != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error creating field:\n"
                "  name:        \"%s\"\n"
                "  location_id: %d\n"
                "  dimension:   %d\n\n"
                "A field with that name has already been defined:\n"
                "  id:          %d\n"
                "  location_id: %d\n"
                "  dimension:   %d"),
              name, location_id, dim, f->id, f->location_id, f->dim);

  /* Initialize if necessary */

  if (_field_map == NULL)
    _field_map = cs_map_name_to_id_create();

  if (l == 0)
    bft_error(__FILE__, __LINE__, 0, _("Defining a field requires a name."));

  for (size_t i = 0; i < l; i++) {
    if (name[i] == '[' || name[i] == ']')
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" is not allowed,\n"
                  "as \'[\' and \']\' are reserved for component access."),
                name);
  }

  /* Insert entry in map */

  field_id = cs_map_name_to_id(_field_map, name);

  if (field_id == _n_fields)
    _n_fields = field_id + 1;

  /* Reallocate fields pointer if necessary */

  if (_n_fields > _n_fields_max) {
    if (_n_fields_max == 0)
      _n_fields_max = 8;
    else
      _n_fields_max *= 2;
    BFT_REALLOC(_fields, _n_fields_max, cs_field_t *);
    BFT_REALLOC(_key_vals, _n_keys_max*_n_fields_max, cs_field_key_val_t);
  }

  /* Allocate fields descriptor block if necessary
     (to reduce fragmentation and improve locality of field
     descriptors, they are allocated in blocks) */

  int shift_in_alloc_block = field_id % _CS_FIELD_S_ALLOC_SIZE;
  if (shift_in_alloc_block == 0)
    BFT_MALLOC(_fields[field_id], _CS_FIELD_S_ALLOC_SIZE, cs_field_t);
  else
    _fields[field_id] = _fields[field_id - shift_in_alloc_block]
                        + shift_in_alloc_block;

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

  f = _fields[field_id];

  f->name = cs_map_name_to_id_reverse(_field_map, field_id);

  f->id = field_id;
  f->type = type_flag;
  f->dim = dim;
  f->location_id = location_id;
  f->n_time_vals = 1;

  f->vals = NULL;
  f->val = NULL;
  f->val_pre = NULL;

  f->grad = NULL;

  f->bc_coeffs = NULL;

  f->is_owner = true;

  /* Mark key values as not set */

  for (key_id = 0; key_id < _n_keys_max; key_id++) {
    memset(&((_key_vals + (f->id*_n_keys_max + key_id))->val),
           0,
           sizeof(union _key_val_t));
    (_key_vals + (f->id*_n_keys_max + key_id))->is_set = 0;
    (_key_vals + (f->id*_n_keys_max + key_id))->is_locked = 0;
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
  cs_real_t  *val = val_old;

  BFT_REALLOC(val, n_elts*dim, cs_real_t);

  /* Initialize field. This should not be necessary, but when using
     threads with Open MP, this should help ensure that the memory will
     first be touched by the same core that will later operate on
     this memory, usually leading to better core/memory affinity. */

  const cs_lnum_t _n_elts = dim * n_elts;
# pragma omp parallel for if (_n_elts > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
    val[ii] = 0.0;

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
    int _n_keys_max_prev = _n_keys_max;
    if (_n_keys_max == 0)
      _n_keys_max = 8;
    else
      _n_keys_max *= 2;
    BFT_REALLOC(_key_defs, _n_keys_max, cs_field_key_def_t);
    BFT_REALLOC(_key_vals, _n_keys_max*_n_fields_max, cs_field_key_val_t);
    for (int field_id = _n_fields - 1; field_id >= 0; field_id--) {
      for (int _key_id = _n_keys - 2; _key_id >= 0; _key_id--)
        _key_vals[field_id*_n_keys_max + _key_id]
          = _key_vals[field_id*_n_keys_max_prev + _key_id];
    }
  }

  for (int field_id = 0; field_id < _n_fields; field_id++) {
    memset((&(_key_vals + (field_id*_n_keys_max + key_id))->val),
           0,
           sizeof(union _key_val_t));
    (_key_vals + (field_id*_n_keys_max + key_id))->is_set = 0;
    (_key_vals + (field_id*_n_keys_max + key_id))->is_locked = 0;
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
  for (int key_id = 0; key_id < _n_keys; key_id++) {

    cs_field_key_def_t *kd = _key_defs + key_id;

    if (kd->type_id == 's') {

      for (int f_id = 0; f_id < _n_fields; f_id++) {
        cs_field_key_val_t *kv = _key_vals + (f_id*_n_keys_max + key_id);
        BFT_FREE(kv->val.v_p);
      }

      if (kd->def_val.v_p != NULL)
        BFT_FREE(kd->def_val.v_p);

    } /* If the key is a "string" key */

  } /* Loop on keys */
}

/*----------------------------------------------------------------------------
 * Free structure associated to a key.
 *----------------------------------------------------------------------------*/

static void
_cs_field_free_struct(void)
{
  for (int key_id = 0; key_id < _n_keys; key_id++) {

    cs_field_key_def_t *kd = _key_defs + key_id;

    if (kd->type_id == 't') {

      for (int f_id = 0; f_id < _n_fields; f_id++) {
        cs_field_key_val_t *kv = _key_vals + (f_id*_n_keys_max + key_id);
        if (kd->clear_func != NULL)
          kd->clear_func(kv->val.v_p);
        BFT_FREE(kv->val.v_p);
      }

      if (kd->def_val.v_p != NULL)
        BFT_FREE(kd->def_val.v_p);

    } /* If the key is a "structure" key */

  } /* Loop on keys */
}

/*----------------------------------------------------------------------------
 * Check if a key may be used with a given field.
 *
 * If the key id is not valid, or the field category is not
 * compatible, a fatal error is provoked.
 *
 * parameters:
 *   f      <-- pointer to field structure
 *   key_id <-- id of associated key
 *
 * returns:
 *   associated error code
 *----------------------------------------------------------------------------*/

static int
_check_key(const cs_field_t  *f,
           int                key_id)
{
  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;

  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
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
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return errcode;
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return the number of defined fields.
 *
 * return:
 *   number of defined fields.
 *----------------------------------------------------------------------------*/

int
cs_f_field_n_fields(void)
{
  return cs_field_n_fields();
}

/*----------------------------------------------------------------------------
 * Return the id of a defined field based on its name.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name <-- field name
 *
 * returns:
 *   id the field structure
 *----------------------------------------------------------------------------*/

int
cs_f_field_id_by_name(const char *name)
{
  int retval;
  cs_field_t  *f = cs_field_by_name(name);

  retval = f->id;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the location of a field.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f <-- field
 *
 * returns:
 *   location
 *----------------------------------------------------------------------------*/

int
cs_f_field_location(const cs_field_t *f)
{
  int retval;

  retval = f->location_id;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the id of a defined field based on its name.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   name <-- field name
 *
 * returns:
 *   id the field structure
 *----------------------------------------------------------------------------*/

int
cs_f_field_id_by_name_try(const char *name)
{
  int retval;
  cs_field_t  *f = cs_field_by_name_try(name);

  if (f != NULL)
    retval = f->id;
  else
    retval = -1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the name of a field defined by its id.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id       <-- field id
 *   name_max <-- maximum name length
 *   name     --> pointer to associated length
 *   name_len --> length of associated length
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_name(int           id,
                    int           name_max,
                    const char  **name,
                    int          *name_len)
{
  const cs_field_t *f = cs_field_by_id(id);
  *name = f->name;
  *name_len = strlen(*name);

  if (*name_len > name_max) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving name from Field %d (\"%s\"):\n"
         "Fortran caller name length (%d) is too small for name \"%s\"\n"
         "(of length %d)."),
       f->id, f->name, name_max, *name, *name_len);
  }
}

/*----------------------------------------------------------------------------
 * Return the dimension of a field defined by its id.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id  <-- field id
 *   dim --> field dimension
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_dimension(int  id,
                         int  dim[1])
{
  const cs_field_t *f = cs_field_by_id(id);

  dim[0] = f->dim;
}

/*----------------------------------------------------------------------------
 * Return the number of previous values of a field.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id  <-- field id
 *   dim --> field dimension
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_n_previous(int  id,
                          int  n_previous[1])
{
  const cs_field_t *f = cs_field_by_id(id);

  n_previous[0] = f->n_time_vals - 1;
}

/*----------------------------------------------------------------------------
 * Return the type flag of a field defined by its id.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id   <-- field id
 *   type <-- field type flag
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_type(int           id,
                    int          *type)
{
  const cs_field_t *f = cs_field_by_id(id);

  *type = f->type;
}

/*----------------------------------------------------------------------------
 * Indicate if a field maintains values at previous time steps
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id  <-- field id
 *
 * returns:
 *   1 if previous values are available, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_f_field_have_previous(int  id)
{
  int retval = 0;
  const cs_field_t *f = cs_field_by_id(id);

  if (f->n_time_vals > 1)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Change a field's handling of values at previous time steps.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id         <-- field id
 *   n_previous <-- number of previous values to save
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_n_previous(int  id,
                          int  n_previous)
{
  cs_field_t *f = cs_field_by_id(id);

  cs_field_set_n_time_vals(f, n_previous + 1);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a field's variable values
 * (current var if previous does not exist)
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id           <-- field id
 *   pointer_type <-- 1: var; 2: var_p or var if var_p does not exists;
 *   pointer_rank <-- expected rank (1 for scalar, 2 for vector)
 *   dim          --> dimensions (indexes in Fortran order,
 *                    dim[i] = 0 if i unused)
 *   p            --> returned pointer
 *
 * returns:
 *   pointer to the field structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_f_field_var_ptr_by_id_try(int          id,
                             int          pointer_type,
                             int          pointer_rank,
                             int          dim[2],
                             cs_real_t  **p)
{
  cs_field_t *f = cs_field_by_id(id);
  int cur_p_rank = 1;

  dim[0] = 0;
  dim[1] = 0;
  *p = NULL;

  if (pointer_type == 1 || pointer_type == 2) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    if (pointer_type == 1 || f->val_pre == NULL)
      *p = f->val;
    else
      *p = f->val_pre;

    if (*p == NULL) /* Adjust dimensions to assist Fortran bounds-checking */
      _n_elts = 0;

    if (f->dim == 1)
      dim[0] = _n_elts;
    else {
      dim[0] = f->dim;
      dim[1] = _n_elts;
      cur_p_rank = 2;
    }

  }

  if (cur_p_rank != pointer_rank)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Fortran pointer of rank %d requested for values of field \"%s\",\n"
         "which have rank %d."),
       pointer_rank, f->name, cur_p_rank);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a field's variable values
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id           <-- field id
 *   pointer_type <-- 1: var; 2: var_prev; 3: var_prev2
 *   pointer_rank <-- expected rank (1 for scalar, 2 for vector)
 *   dim          --> dimensions (indexes in Fortran order,
 *                    dim[i] = 0 if i unused)
 *   p            --> returned pointer
 *
 * returns:
 *   pointer to the field structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_f_field_var_ptr_by_id(int          id,
                         int          pointer_type,
                         int          pointer_rank,
                         int          dim[2],
                         cs_real_t  **p)
{
  cs_field_t *f = cs_field_by_id(id);
  int cur_p_rank = 1;

  dim[0] = 0;
  dim[1] = 0;
  *p = NULL;

  if (pointer_type > f->n_time_vals)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Fortran pointer with %d previous values of field \"%s\",\n"
         "requests the %d previous values."),
       f->n_time_vals-1, f->name, pointer_type-1);

  if (pointer_type == 1 || pointer_type == 2 || pointer_type == 3) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    assert(pointer_type <= f->n_time_vals);

    *p = f->vals[pointer_type - 1];

    if (*p == NULL) /* Adjust dimensions to assist Fortran bounds-checking */
      _n_elts = 0;

    /* If dimension 1 is asked and field is of dimension one */
    if (f->dim == 1 && pointer_rank == 1)
      dim[0] = _n_elts;
    else {
      dim[0] = f->dim;
      dim[1] = _n_elts;
      cur_p_rank = 2;
    }

  }

  if (cur_p_rank != pointer_rank)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Fortran pointer of rank %d requested for values of field \"%s\",\n"
         "which have rank %d."),
       pointer_rank, f->name, cur_p_rank);
}

/*----------------------------------------------------------------------------
 * Return a pointer to a field's boundary condition coefficient values
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id           <-- field id
 *   pointer_type <-- 1: bc_coeffs->a;     2: bc_coeffs->b
 *                    3: bc_coeffs->af;    4: bc_coeffs->bf
 *                    5: bc_coeffs->ad;    6: bc_coeffs->bd
 *                    7: bc_coeffs->ac;    8: bc_coeffs->bc
 *   pointer_rank <-- expected rank (1 for scalar, 2 for vector)
 *   dim          <-- dimensions (indexes in Fortran order,
 *                    dim[i] = 0 if i unused)
 *   p            <-- returned pointer
 *
 * returns:
 *   pointer to the field structure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_f_field_bc_coeffs_ptr_by_id(int          id,
                               int          pointer_type,
                               int          pointer_rank,
                               int          dim[3],
                               cs_real_t  **p)
{
  cs_field_t *f = cs_field_by_id(id);
  int cur_p_rank = 1;

  dim[0] = 0;
  dim[1] = 0;
  dim[2] = 0;
  *p = NULL;

  const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(location_id);
  cs_lnum_t _n_elts = n_elts[2];

  assert(f->location_id == CS_MESH_LOCATION_CELLS);

  if (f->bc_coeffs == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s\"\n"
                " does not have associated BC coefficients."),
              f->name);

  if (f->bc_coeffs != NULL) {

    if (pointer_type == 1)
      *p = f->bc_coeffs->a;
    else if (pointer_type == 2)
      *p = f->bc_coeffs->b;
    else if (pointer_type == 3)
      *p = f->bc_coeffs->af;
    else if (pointer_type == 4)
      *p = f->bc_coeffs->bf;
    else if (pointer_type == 5)
      *p = f->bc_coeffs->ad;
    else if (pointer_type == 6)
      *p = f->bc_coeffs->bd;
    else if (pointer_type == 7)
      *p = f->bc_coeffs->ac;
    else if (pointer_type == 8)
      *p = f->bc_coeffs->bc;

    if (*p == NULL) /* Adjust dimensions to assist Fortran bounds-checking */
      _n_elts = 0;

    if (f->dim == 1 || pointer_type == 9 || pointer_type == 10)
      dim[0] = _n_elts;

    else {

      int coupled = 0;

      if (f->type & CS_FIELD_VARIABLE) {
        int coupled_key_id = cs_field_key_id_try("coupled");
        if (coupled_key_id > -1)
          coupled = cs_field_get_key_int(f, coupled_key_id);
      }

      if (coupled) {

        if (pointer_type == 1 || pointer_type == 3 || pointer_type == 5
            || pointer_type == 7) {
          dim[0] = f->dim;
          dim[1] = _n_elts;
          cur_p_rank = 2;
        }
        else { /* if (pointer_type == 2 || pointer_type == 4 || pointer_type == 6
                      || pointer_type == 8) */
          dim[0] = f->dim;
          dim[1] = f->dim;
          dim[2] = _n_elts;
          cur_p_rank = 3;
        }

      }
      else { /* uncoupled */

        dim[0] = f->dim;
        dim[1] = _n_elts;
        cur_p_rank = 2;

      }

    }

  }

  if (cur_p_rank != pointer_rank)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Fortran pointer of rank %d requested for BC coefficients of field\n"
         " \"%s\", which have rank %d."),
       pointer_rank, f->name, cur_p_rank);
}

/*----------------------------------------------------------------------------
 * Assign an integer value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id  <-- field id
 *   k_id  <-- key id
 *   value <-- associated value
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_int(int  f_id,
                       int  k_id,
                       int  value)
{
  int retval = 0;

  cs_field_t *f = cs_field_by_id(f_id);

  retval = cs_field_set_key_int(f, k_id, value);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, k_id);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning integer value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, k_id, key);
  }
}

/*----------------------------------------------------------------------------
 * Set integer bits matching a mask to 1 for a given key for a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id <-- field id
 *   k_id <-- key id
 *   mask <-- associated mask
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_int_bits(int  f_id,
                            int  k_id,
                            int  mask)
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_key_int_bits(f, k_id, mask);
}

/*----------------------------------------------------------------------------
 * Set integer bits matching a mask to 0 for a given key for a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id <-- field id
 *   k_id <-- key id
 *   mask <-- associated mask
 *----------------------------------------------------------------------------*/

void
cs_f_field_clear_key_int_bits(int  f_id,
                              int  k_id,
                              int  mask)
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_clear_key_int_bits(f, k_id, mask);
}

/*----------------------------------------------------------------------------
 * Assign a floating point value for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id  <-- field id
 *   k_id  <-- key id
 *   value <-- associated value
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_double(int     f_id,
                          int     k_id,
                          double  value)
{
  int retval = 0;

  cs_field_t *f = cs_field_by_id(f_id);

  retval = cs_field_set_key_double(f, k_id, value);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, k_id);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning real value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, k_id, key);
  }
}

/*----------------------------------------------------------------------------
 * Assign a character string for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id <-- field id
 *   k_id <-- key id
 *   str  <-- associated string
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_str(int          f_id,
                       int          k_id,
                       const char  *str)
{
  cs_field_t *f = cs_field_by_id(f_id);
  int retval = cs_field_set_key_str(f, k_id, str);

  if (retval != 0) {
    const char *key = cs_map_name_to_id_reverse(_key_map, k_id);
    bft_error(__FILE__, __LINE__, 0,
              _("Error %d assigning string value to Field \"%s\" with\n"
                "type flag %d with key %d (\"%s\")."),
              retval, f->name, f->type, k_id, key);
  }
}

/*----------------------------------------------------------------------------
 * Return a character string for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id    <-- field id
 *   k_id    <-- id of associated key
 *   str_max <-- maximum string length
 *   str     --> pointer to associated string
 *   str_len --> length of associated string
 *
 * returns:
 *   pointer to character string
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_key_str(int           f_id,
                       int           key_id,
                       int           str_max,
                       const char  **str,
                       int          *str_len)
{
  const cs_field_t *f = cs_field_by_id(f_id);
  *str = cs_field_get_key_str(f, key_id);

  if (str != NULL)
    *str_len = strlen(*str);
  else
    *str_len = 0;

  if (*str_len > str_max) {
    const char *key = cs_map_name_to_id_reverse(_key_map, key_id);
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving string from Field %d (\"%s\") and key %d (\"%s\"):\n"
         "Fortran caller string length (%d) is too small for string \"%s\"\n"
         "(of length %d)."),
       f->id, f->name, key_id, key, str_max, *str, *str_len);
  }
}

/*----------------------------------------------------------------------------
 * Assign a simple structure for a given key to a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id    <-- field id
 *   k_id    <-- id of associated key
 *   k_value --> pointer to structure
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_key_struct(int    f_id,
                          int    k_id,
                          void  *k_value)
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_key_struct(f, k_id, k_value);
}

/*----------------------------------------------------------------------------
 * Copy a structure for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id    <-- field id
 *   k_id    <-- id of associated key
 *   k_value --> pointer to structure
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_key_struct(int    f_id,
                          int    k_id,
                          void  *k_value)
{
  const cs_field_t *f = cs_field_by_id(f_id);

  cs_field_get_key_struct(f, k_id, k_value);
}

/*----------------------------------------------------------------------------
 * Return a label associated with a field.
 *
 * If the "label" key has been set for this field, its associated string
 * is returned. Otherwise, the field's name is returned.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   f_id    <-- field id
 *   str_max <-- maximum string length
 *   str     --> pointer to associated string
 *   str_len --> length of associated string
 *
 * returns:
 *   pointer to character string
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_label(int           f_id,
                     int           str_max,
                     const char  **str,
                     int          *str_len)
{
  const cs_field_t *f = cs_field_by_id(f_id);
  *str = cs_field_get_label(f);

  *str_len = strlen(*str);

  if (*str_len > str_max) {
    const char *key = cs_map_name_to_id_reverse(_key_map, _k_label);
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error retrieving string from Field %d (\"%s\") and key %d (\"%s\"):\n"
         "Fortran caller string length (%d) is too small for string \"%s\"\n"
         "(of length %d)."),
       f->id, f->name, _k_label, key, str_max, *str, *str_len);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
  cs_field_t  *f = _field_create(name,
                                 type_flag,
                                 location_id,
                                 dim);

  cs_base_check_bool(&has_previous);

  f->n_time_vals = has_previous ? 2 : 1;

  BFT_MALLOC(f->vals, f->n_time_vals, cs_real_t *);
  for (int i = 0; i < f->n_time_vals; i++)
    f->vals[i] = NULL;

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a field matching a given name and attributes,
 *        creating it if necessary.
 *
 * If a field with the same name but different attributes is present,
 * this is considered an error.
 *
 * The default number of time values associated with a field created through
 * this function is 1. To modify it, use \ref cs_field_set_n_time_vals.
 *
 * \param[in]  name          field name
 * \param[in]  type_flag     mask of field property and category values
 * \param[in]  location_id   id of associated location
 * \param[in]  dim           field dimension (number of components)
 * \param[in]  has_previous  maintain values at the previous time step ?
 *
 * \return  pointer to field
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_field_find_or_create(const char   *name,
                        int           type_flag,
                        int           location_id,
                        int           dim,
                        bool          has_previous)
{
  cs_field_t *f = cs_field_by_name_try(name);

  if (f != NULL) {

    if (   type_flag != f->type || location_id != f->location_id
        || dim != f->dim) {
      bft_error(__FILE__, __LINE__, 0,
                _("Mismatch in field definitions:\n"
                  "  name:        \"%s\"\n"
                  "  type_flag:   %d\n"
                  "  location_id: %d\n"
                  "  dimension:   %d\n\n"
                  "A previous definition for that has attributes:\n"
                  "  id:          %d\n"
                  "  type_flag:   %d\n"
                  "  location_id: %d\n"
                  "  dimension:   %d\n\n"),
                name, type_flag, location_id, dim,
                f->id, f->type, f->location_id, f->dim);
    }

  }
  else {

    f =  _field_create(name,
                       type_flag,
                       location_id,
                       dim);

    cs_base_check_bool(&has_previous);

    f->n_time_vals = has_previous ? 2 : 1;

    BFT_MALLOC(f->vals, f->n_time_vals, cs_real_t *);
    for (int i = 0; i < f->n_time_vals; i++)
      f->vals[i] = NULL;

  }

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Change the number of time values managed by a field.
 *
 * The minimum will never be below 1, as the current time is always handled.
 *
 * \param[in, out]  f            pointer to field structure
 * \param[in]       n_time_vals  number of time values to maintain
 */
/*----------------------------------------------------------------------------*/

void
cs_field_set_n_time_vals(cs_field_t  *f,
                         int          n_time_vals)
{
  assert(f != NULL);
  if (f == NULL)
    return;

  int _n_time_vals = n_time_vals;

  const int n_time_vals_ini = f->n_time_vals;

  if (_n_time_vals < 1)
    _n_time_vals = 1;

  else if (_n_time_vals > 3)
    bft_error(__FILE__, __LINE__, 0,
              "%s called for field \"%s\" with n_time_vals = %d\n"
              " but only values 1, 2 and 3 are currently supported.",
              __func__, f->name, n_time_vals);
  else
    _n_time_vals = n_time_vals;

  if (_n_time_vals == n_time_vals_ini)
    return;

  /* Update number of time values */

  f->n_time_vals = _n_time_vals;

  BFT_REALLOC(f->vals, f->n_time_vals, cs_real_t *);
  for (int i = n_time_vals_ini; i < f->n_time_vals; i++)
    f->vals[i] = NULL;

  /* If allocation or mapping has already been done */

  if (f->val != NULL) {
    if (n_time_vals_ini > _n_time_vals) {
      assert(n_time_vals_ini == 2 && _n_time_vals == 1);
      if (f->is_owner)
        BFT_FREE(f->val_pre);
      else
        f->val_pre = NULL;
    }
    else { /* if (n_time_vals_ini < _n_time_vals) */
      if (f->is_owner) {
        const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
        f->val_pre = _add_val(n_elts[2], f->dim, f->val_pre);
      }
    }
  }
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
    int ii;

    /* Initialization */

    for (ii = 0; ii < f->n_time_vals; ii++)
      f->vals[ii] = _add_val(n_elts[2], f->dim, f->vals[ii]);

    f->val = f->vals[0];
    if (f->n_time_vals > 1)
      f->val_pre = f->vals[1];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Map existing value arrays to field descriptor.
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
  if (f == NULL)
    return;

  if (f->is_owner) {
    BFT_FREE(f->val);
    BFT_FREE(f->val_pre);
    f->is_owner = false;
  }

  f->val = val;
  f->vals[0] = val;

  /* Add previous time step values if necessary */

  if (f->n_time_vals > 1) {
    f->val_pre = val_pre;
    f->vals[1] = val_pre;
  }
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
 * For multidimensional fields with coupled components, implicit b and bf
 * coefficient arrays are arrays of block matrices, not vectors, so the
 * number of entries for each boundary face is dim*dim instead of dim.
 *
 * \param[in, out]  f             pointer to field structure
 * \param[in]       have_flux_bc  if true, flux bc coefficients (af and bf)
 *                                are added
 * \param[in]       have_mom_bc   if true, div BC coefficients (ad and bd)
 *                                are added
 * \param[in]       have_conv_bc  if true, convection BC coefficients (ac and bc)
 *                                are added
 * \param[in]       have_exch_bc  if true, exchange boundary coefficients (hint
 *                                and hext) are added
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_bc_coeffs(cs_field_t  *f,
                            bool         have_flux_bc,
                            bool         have_mom_bc,
                            bool         have_conv_bc,
                            bool         have_exch_bc)
{
  /* Add boundary condition coefficients if required */

  cs_lnum_t a_mult = f->dim;
  cs_lnum_t b_mult = f->dim;

  cs_base_check_bool(&have_flux_bc);
  cs_base_check_bool(&have_mom_bc);
  cs_base_check_bool(&have_conv_bc);

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

      f->bc_coeffs->icodcl = NULL;
      f->bc_coeffs->rcodcl1 = NULL;
      f->bc_coeffs->rcodcl2 = NULL;
      f->bc_coeffs->rcodcl3 = NULL;

      BFT_MALLOC(f->bc_coeffs->a, n_elts[0]*a_mult, cs_real_t);
      BFT_MALLOC(f->bc_coeffs->b, n_elts[0]*b_mult, cs_real_t);

      if (have_flux_bc) {
        BFT_MALLOC(f->bc_coeffs->af, n_elts[0]*a_mult, cs_real_t);
        BFT_MALLOC(f->bc_coeffs->bf, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        f->bc_coeffs->af = NULL;
        f->bc_coeffs->bf = NULL;
      }

      if (have_mom_bc) {
        BFT_MALLOC(f->bc_coeffs->ad, n_elts[0]*a_mult, cs_real_t);
        BFT_MALLOC(f->bc_coeffs->bd, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        f->bc_coeffs->ad = NULL;
        f->bc_coeffs->bd = NULL;
      }

      if (have_conv_bc) {
        BFT_MALLOC(f->bc_coeffs->ac, n_elts[0]*a_mult, cs_real_t);
        BFT_MALLOC(f->bc_coeffs->bc, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        f->bc_coeffs->ac = NULL;
        f->bc_coeffs->bc = NULL;
      }

      if (have_exch_bc) {
        BFT_MALLOC(f->bc_coeffs->hint, n_elts[0], cs_real_t);
        BFT_MALLOC(f->bc_coeffs->_hext, n_elts[0], cs_real_t);
      }
      else {
        f->bc_coeffs->hint = NULL;
        f->bc_coeffs->_hext = NULL;
      }

    }

    else {

      BFT_REALLOC(f->bc_coeffs->a, n_elts[0]*a_mult, cs_real_t);
      BFT_REALLOC(f->bc_coeffs->b, n_elts[0]*b_mult, cs_real_t);

      if (have_flux_bc) {
        BFT_REALLOC(f->bc_coeffs->af, n_elts[0]*a_mult, cs_real_t);
        BFT_REALLOC(f->bc_coeffs->bf, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->af);
        BFT_FREE(f->bc_coeffs->bf);
      }

      if (have_mom_bc) {
        BFT_REALLOC(f->bc_coeffs->ad, n_elts[0]*a_mult, cs_real_t);
        BFT_REALLOC(f->bc_coeffs->bd, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->ad);
        BFT_FREE(f->bc_coeffs->bd);
      }

      if (have_conv_bc) {
        BFT_REALLOC(f->bc_coeffs->ac, n_elts[0]*a_mult, cs_real_t);
        BFT_REALLOC(f->bc_coeffs->bc, n_elts[0]*b_mult, cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->ac);
        BFT_FREE(f->bc_coeffs->bc);
      }

      if (have_exch_bc) {
        BFT_MALLOC(f->bc_coeffs->hint, n_elts[0], cs_real_t);
        BFT_MALLOC(f->bc_coeffs->_hext, n_elts[0], cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->hint);
        BFT_FREE(f->bc_coeffs->_hext);
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
 * \brief  Initialize boundary condition coefficients arrays.
 *
 * For fields on location CS_MESH_LOCATION_CELLS, boundary conditions
 * are located on CS_MESH_LOCATION_BOUNDARY_FACES.
 *
 * Boundary condition coefficients are not currently supported for other
 * locations (though support could be added by mapping a boundary->location
 * indirection array in the cs_mesh_location_t structure).
 *
 * For multidimensional fields with coupled components, implicit b and bf
 * coefficient arrays are arrays of block matrices, not vectors, so the
 * number of entries for each boundary face is dim*dim instead of dim.
 *
 * \param[in, out]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_init_bc_coeffs(cs_field_t  *f)
{
  /* Add boundary condition coefficients if required */

  cs_lnum_t dim = f->dim;

  int ifac;
  int coupled = 0;

  if (f->type & CS_FIELD_VARIABLE) {
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > -1)
      coupled = cs_field_get_key_int(f, coupled_key_id);
  }

  if (f->location_id == CS_MESH_LOCATION_CELLS) {

    const int location_id = CS_MESH_LOCATION_BOUNDARY_FACES;
    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(location_id);

    if (coupled == 0 && dim == 1) {

      for (ifac = 0; ifac < n_elts[0]; ifac++) {
        f->bc_coeffs->a[ifac] = 0.;
        f->bc_coeffs->b[ifac] = 1.;
      }

      if (f->bc_coeffs->af != NULL)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->af[ifac] = 0.;
          f->bc_coeffs->bf[ifac] = 0.;
        }

      if (f->bc_coeffs->ad != NULL) {
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ad[ifac] = 0.;
          f->bc_coeffs->bd[ifac] = 1.;
        }
      }

      if (f->bc_coeffs->ac != NULL) {
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ac[ifac] = 0.;
          f->bc_coeffs->bc[ifac] = 0.;
        }
      }

    }

    /* Coupled vectorial BCs */
    else if (coupled && dim == 3) {

      for (ifac = 0; ifac < n_elts[0]; ifac++) {
        f->bc_coeffs->a[ifac*dim] = 0.;
        f->bc_coeffs->a[ifac*dim + 1] = 0.;
        f->bc_coeffs->a[ifac*dim + 2] = 0.;
        f->bc_coeffs->b[ifac*dim*dim] = 1.;
        f->bc_coeffs->b[ifac*dim*dim + 1] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 2] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 3] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 4] = 1.;
        f->bc_coeffs->b[ifac*dim*dim + 5] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 6] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 7] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 8] = 1.;
      }

      if (f->bc_coeffs->af != NULL)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->af[ifac*dim] = 0.;
          f->bc_coeffs->af[ifac*dim + 1] = 0.;
          f->bc_coeffs->af[ifac*dim + 2] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 1] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 2] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 3] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 4] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 5] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 6] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 7] = 0.;
          f->bc_coeffs->bf[ifac*dim*dim + 8] = 0.;
        }

      if (f->bc_coeffs->ad != NULL)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ad[ifac*dim] = 0.;
          f->bc_coeffs->ad[ifac*dim + 1] = 0.;
          f->bc_coeffs->ad[ifac*dim + 2] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim] = 1.;
          f->bc_coeffs->bd[ifac*dim*dim + 1] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 2] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 3] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 4] = 1.;
          f->bc_coeffs->bd[ifac*dim*dim + 5] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 6] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 7] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 8] = 1.;
        }

      if (f->bc_coeffs->ac != NULL)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ac[ifac*dim] = 0.;
          f->bc_coeffs->ac[ifac*dim + 1] = 0.;
          f->bc_coeffs->ac[ifac*dim + 2] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 1] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 2] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 3] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 4] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 5] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 6] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 7] = 0.;
          f->bc_coeffs->bc[ifac*dim*dim + 8] = 0.;
        }

    }
    else {
      for (ifac = 0; ifac < n_elts[0]; ifac++) {
        for (int isou = 0; isou < dim ; isou++) {
          for (int jsou = 0; jsou < dim; jsou ++) {
            f->bc_coeffs->b[ifac*dim*dim + isou*dim +jsou] = 0.;
            if (isou == jsou) {
              f->bc_coeffs->b[ifac*dim*dim + isou*dim +jsou] = 1.;
            }
          }
        }
      }

      if (f->bc_coeffs->af != NULL) {
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          for (int isou = 0; isou < dim ; isou++) {
            f->bc_coeffs->af[ifac*dim + isou] = 0.;
            for (int jsou = 0; jsou < dim; jsou ++) {
              f->bc_coeffs->bf[ifac*dim*dim + isou*dim +jsou] = 0.;
            }
          }
        }
      }

      if (f->bc_coeffs->ad != NULL) {
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          for (int isou = 0; isou < dim ; isou++) {
            f->bc_coeffs->ad[ifac*dim + isou] = 0.;
            for (int jsou = 0; jsou < dim; jsou ++) {
              f->bc_coeffs->bd[ifac*dim*dim + isou*dim +jsou] = 0.;
              if (isou == jsou) {
                f->bc_coeffs->bd[ifac*dim*dim + isou*dim +jsou] = 1.;
              }
            }
          }
        }
      }

      if (f->bc_coeffs->ac != NULL) {
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          for (int isou = 0; isou < dim ; isou++) {
            f->bc_coeffs->ac[ifac*dim + isou] = 0.;
            for (int jsou = 0; jsou < dim; jsou ++) {
              f->bc_coeffs->bc[ifac*dim*dim + isou*dim +jsou] = 0.;
            }
          }
        }
      }
    }

    if (f->bc_coeffs->hint != NULL) {
      for (ifac = 0; ifac < n_elts[0]; ifac++) {
        f->bc_coeffs->hint[ifac] = 0.;
      }
    }

    if (f->bc_coeffs->_hext != NULL) {
      for (ifac = 0; ifac < n_elts[0]; ifac++) {
        f->bc_coeffs->_hext[ifac] = 0.;
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
 * \brief  Allocate arrays for field gradient.
 *
 * \param[in, out]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_gradient(cs_field_t  *f)
{
  assert(f != NULL);

  if (f->is_owner) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);

    /* Initialization */

    f->grad = _add_val(n_elts[2], 3*f->dim, f->grad);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set current field values to the given constant.
 *
 * \param[in, out]  f  pointer to field structure
 * \param[in]       c  assigned value
 */
/*----------------------------------------------------------------------------*/

void
cs_field_set_values(cs_field_t  *f,
                    cs_real_t    c)
{
  assert(f != NULL);
  if (f == NULL)
    return;

  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
  const cs_lnum_t _n_vals = n_elts[2]*f->dim;

# pragma omp parallel for if (_n_vals > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < _n_vals; ii++)
    f->val[ii] = c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current field values to previous values if applicable.
 *
 * For fields with only one time value, or values not allocated yet,
 * this is a no-op.
 *
 * \param[in, out]  f  pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_current_to_previous(cs_field_t  *f)
{
  assert(f != NULL);

  if (f->n_time_vals > 1) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    const cs_lnum_t _n_elts = n_elts[2];

#   pragma omp parallel if (_n_elts > CS_THR_MIN)
    {
      const int dim = f->dim;

      if (f->is_owner) {
        if (dim == 1) {
          for (int kk = f->n_time_vals - 1; kk > 0; kk--) {
#           pragma omp for
            for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
              f->vals[kk][ii] = f->vals[kk-1][ii];
          }
        }
        else {
          for (int kk = f->n_time_vals - 1; kk > 0; kk--) {
#           pragma omp for
            for (cs_lnum_t ii = 0; ii < _n_elts; ii++) {
              for (cs_lnum_t jj = 0; jj < dim; jj++)
                f->vals[kk][ii*dim + jj] = f->vals[kk-1][ii*dim + jj];
            }
          }
        }
      }
      else {
        if (dim == 1) {
#         pragma omp for
          for (cs_lnum_t ii = 0; ii < _n_elts; ii++)
            f->val_pre[ii] = f->val[ii];
        }
        else {
#         pragma omp for
          for (cs_lnum_t ii = 0; ii < _n_elts; ii++) {
            for (cs_lnum_t jj = 0; jj < dim; jj++)
              f->val_pre[ii*dim + jj] = f->val[ii*dim + jj];
          }
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all defined fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_destroy_all(void)
{
  for (int i = 0; i < _n_fields; i++) {
    cs_field_t  *f = _fields[i];
    if (f->is_owner) {
      if (f->vals != NULL) {
        int ii;
        for (ii = 0; ii < f->n_time_vals; ii++)
          BFT_FREE(f->vals[ii]);
      }
    }
    BFT_FREE(f->vals);

    if (f->grad != NULL)
      BFT_FREE(f->grad);

    if (f->bc_coeffs != NULL) {
      BFT_FREE(f->bc_coeffs->a);
      BFT_FREE(f->bc_coeffs->b);
      BFT_FREE(f->bc_coeffs->af);
      BFT_FREE(f->bc_coeffs->bf);
      BFT_FREE(f->bc_coeffs->ad);
      BFT_FREE(f->bc_coeffs->bd);
      BFT_FREE(f->bc_coeffs->ac);
      BFT_FREE(f->bc_coeffs->bc);
      BFT_FREE(f->bc_coeffs->hint);
      BFT_FREE(f->bc_coeffs->_hext);
      BFT_FREE(f->bc_coeffs);
    }
  }

  for (int i = 0; i < _n_fields; i++) {
    if (i % _CS_FIELD_S_ALLOC_SIZE == 0)
      BFT_FREE(_fields[i]);
  }

  BFT_FREE(_fields);

  cs_map_name_to_id_destroy(&_field_map);

  _cs_field_free_str();
  _cs_field_free_struct();

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
    cs_field_t  *f = _fields[i];
    if (f->is_owner)
      cs_field_allocate_values(f);
    else {
      if (f->val == NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("Field \"%s\"\n"
                    " requires mapped values which have not been set."),
                  f->name);
    }
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
    return _fields[id];
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
    return _fields[id];
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
    return _fields[id];
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a field based on a composite name.
 *
 * The name is expected to be of the form <name_prefix>_<name_suffix>.
 *
 * \param[in]  name_prefix  first part of field name
 * \param[in]  name_suffix  second part of field name
 *
 * \return  pointer to the field structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_composite_name(const char  *name_prefix,
                           const char  *name_suffix)
{
  cs_field_t *f = cs_field_by_composite_name_try(name_prefix,
                                                 name_suffix);

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Field \"%s_%s\" is not defined."), name_prefix, name_suffix);

  return f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to a field based on a composite name if present.
 *
 * The name is expected to be of the form <name_prefix>_<name_suffix>.
 * If no field of the given name is defined, NULL is returned.
 *
 * \param[in]  name_prefix  first part of field name
 * \param[in]  name_suffix  second part of field name
 *
 * \return  pointer to the field structure, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_field_t  *
cs_field_by_composite_name_try(const char  *name_prefix,
                               const char  *name_suffix)
{
  size_t lp = strlen(name_prefix);
  size_t ls = strlen(name_suffix);
  size_t lt = lp + ls + 1;

  char _buffer[196];
  char *buffer = _buffer;

  if (lt + 1> 196)
    BFT_MALLOC(buffer, lt+1, char);

  memcpy(buffer, name_prefix, lp);
  buffer[lp] = '_';
  memcpy(buffer + lp + 1, name_suffix, ls);
  buffer[lt] = '\0';

  int id = cs_map_name_to_id_try(_field_map, buffer);

  if (buffer != _buffer)
    BFT_FREE(buffer);

  if (id > -1)
    return _fields[id];
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of a defined field based on its name.
 *
 * If no field with the given name exists, -1 is returned.
 *
 * \param[in]  name   key name
 *
 * \return  id of the field, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

int
cs_field_id_by_name(const char *name)
{
  int id = cs_map_name_to_id_try(_field_map, name);

  return id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of a defined field and an associated component
 *  based on a component name.
 *
 * If no field with the given name exists, -1 is returned.
 *
 * \param[in]   name   field or field+component name
 * \param[out]  f_id   field id, or -1 if no match was found
 * \param[out]  c_id   component id, or -1 for all components
 */
/*----------------------------------------------------------------------------*/

void
cs_field_component_id_by_name(const char  *name,
                              int         *f_id,
                              int         *c_id)
{
  size_t l = strlen(name);

  *f_id = -1;
  *c_id = -1;

  /* Case with an extension */

  if (l > 3) {
    if (name[l-1] == ']') {
      size_t l0 = -1;
      char _name0[128];
      char *name0 = _name0;
      if (l >= 128)
        BFT_MALLOC(name0, l + 1, char);
      strcpy(name0, name);
      for (l0 = l-2; l0 > 0; l0--) {
        if (name0[l0] == '[') {
          name0[l0] = '\0';
          *f_id = cs_map_name_to_id_try(_field_map, name0);
          break;
        }
        else
          name0[l0] = toupper(name0[l0]);
      }
      if (*f_id > -1) {
        cs_field_t *f = cs_field_by_id(*f_id);
        const char **c_name;
        switch (f->dim) {
        case 3:
          c_name = cs_glob_field_comp_name_3;
          break;
        case 6:
          c_name = cs_glob_field_comp_name_6;
          break;
        case 9:
          c_name = cs_glob_field_comp_name_9;
          break;
        default:
          c_name = NULL;
        }
        if (c_name != NULL) {
          for (int _c_id = 0; *c_id < 0 &&_c_id < f->dim; _c_id++) {
            if (strcmp(name0 + l0 + 1, c_name[_c_id]) == 0)
              *c_id = _c_id;
          }
        }
        if (*c_id < 0 && l-l0 < 63) {
          char c_str[64], c_ref[64];
          strncpy(c_str, name0 + l0 + 1, 63);
          c_str[l - l0 - 2] = '\0';
          for (int _c_id = 0; *c_id < 0 &&_c_id < f->dim; _c_id++) {
            sprintf(c_ref, "%d", _c_id);
            if (strcmp(c_str, c_ref) == 0)
              *c_id = _c_id;
          }
        }
        if (*c_id < 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("Field \"%s\" does not have a component \"%s\"."),
                    f->name, name + l0 - 1);
      }
      if (name0 != _name0)
        BFT_FREE(name0);
    }
  }

  /* Case with no extension */

  if (*f_id == -1)
    *f_id = cs_map_name_to_id_try(_field_map, name);
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

  kd->def_val.v_int = default_value;
  kd->log_func = NULL;
  kd->type_size = 0;
  kd->type_flag = type_flag;
  kd->type_id = 'i';
  kd->log_id = 's';
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

  kd->def_val.v_double = default_value;
  kd->log_func = NULL;
  kd->type_size = 0;
  kd->type_flag = type_flag;
  kd->type_id = 'd';
  kd->log_id = 's';
  kd->is_sub = false;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a key for a string value by its name and return an
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

  /* Free possible previous allocation */
  if (n_keys_init == _n_keys)
    BFT_FREE(kd->def_val.v_p);

  if (default_value != NULL) {
    BFT_MALLOC(kd->def_val.v_p, strlen(default_value) + 1, char);
    strcpy(kd->def_val.v_p, default_value);
  }
  else
    kd->def_val.v_p = NULL;
  kd->log_func = NULL;
  kd->type_size = 0;
  kd->type_flag = type_flag;
  kd->type_id = 's';
  kd->log_id = 's';
  kd->is_sub = false;

  return key_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a key for a structure value by its name and return an
 * associated id.
 *
 * If the key has already been defined, its previous default value is replaced
 * by the current value, and its id is returned.
 *
 * \param[in]  name              key name
 * \param[in]  default_value     pointer to default value associated with key
 * \param[in]  log_func          pointer to logging function
 * \param[in]  log_func_default  pointer to default logging function
 * \param[in]  clear_func        pointer to substructures free function
 * \param[in]  size              sizeof structure
 * \param[in]  type_flag         mask associated with field types with which
 *                               the key may be associated, or 0
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_key_struct(const char                   *name,
                           const void                   *default_value,
                           cs_field_log_key_struct_t    *log_func,
                           cs_field_log_key_struct_t    *log_func_default,
                           cs_field_clear_key_struct_t  *clear_func,
                           size_t                        size,
                           int                           type_flag)
{
  int n_keys_init = _n_keys;

  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;

  /* Free possible previous allocation */
  if (n_keys_init == _n_keys)
    BFT_FREE(kd->def_val.v_p);

  if (default_value != NULL) {
    BFT_MALLOC(kd->def_val.v_p, size, unsigned char);
    memcpy(kd->def_val.v_p, default_value, size);
  }
  else
    kd->def_val.v_p = NULL;
  kd->log_func = log_func;
  kd->log_func_default = log_func_default;
  kd->clear_func = clear_func;
  kd->type_size = size;
  kd->type_flag = type_flag;
  kd->type_id = 't';
  kd->log_id = 's';
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

  assert(parent_id > -1 && parent_id < _n_keys);

  kd->def_val.v_int = parent_id;
  kd->type_flag = pkd->type_flag;
  kd->type_id = pkd->type_id;
  kd->log_id = pkd->log_id;
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
  int key_id;
  for (key_id = 0; key_id < _n_keys; key_id++) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    if (kd->type_id == 's' || kd->type_id == 't') {
      BFT_FREE(kd->def_val.v_p);
    }
  }

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
 * \brief Disable logging setup values associated with a given key.
 *
 * This is useful when a key is used not for setup purposes, but to track
 * values associated with a field, such as convergence or performance data.
 *
 * \param[in]  key_id  id of associated key
 */
/*----------------------------------------------------------------------------*/

void
cs_field_key_disable_setup_log(int  key_id)
{
  assert(key_id >= 0 && key_id < _n_keys);
  cs_field_key_def_t *kd = _key_defs + key_id;
  kd->log_id = 'n';
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query if a given key has been set for a field.
 *
 * If the key id is not valid, or the field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  true if the key has been set for this field, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_field_is_key_set(const cs_field_t  *f,
                    int                key_id)
{
  int errcode = _check_key(f, key_id);

  if (errcode == CS_FIELD_OK) {
    cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
    bool retval = false;
    if (kv->is_set)
      retval = true;
    return retval;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query if a given key has been locked for a field.
 *
 * If the key id is not valid, or the field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  true if the key has been locked for this field, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_field_is_key_locked(const cs_field_t  *f,
                       int                key_id)
{
  int errcode = _check_key(f, key_id);

  if (errcode == CS_FIELD_OK) {
    cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
    bool retval = false;
    if (kv->is_locked)
      retval = true;
    return retval;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lock a field relative to a given key.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_lock_key(cs_field_t  *f,
                  int          key_id)
{
  int retval = CS_FIELD_OK;

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      kv->is_locked = 1;
    }
  }
  else
    retval = CS_FIELD_INVALID_KEY_ID;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a integer value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 * If the data type does not match, CS_FIELD_INVALID_TYPE is returned.
 * If the key value has been locked, CS_FIELD_LOCKED is returned.
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

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;
  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 'i')
      retval = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      if (kv->is_locked)
        retval = CS_FIELD_LOCKED;
      else {
        kv->val.v_int = value;
        kv->is_set = 1;
      }
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

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;
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
        retval = kv->val.v_int;
      else if (kd->is_sub)
        retval = cs_field_get_key_int(f, kd->def_val.v_int);
      else
        retval = kd->def_val.v_int;
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

  return CS_FIELD_OK;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set integer bits matching a mask to 1 for a given key for a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 * If the data type does not match, CS_FIELD_INVALID_TYPE is returned.
 * If the key value has been locked, CS_FIELD_LOCKED is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  mask    mask associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_set_key_int_bits(cs_field_t  *f,
                          int          key_id,
                          int          mask)
{
  int value = cs_field_get_key_int(f, key_id);

  value |= mask;

  int retval = cs_field_set_key_int(f, key_id, value);
  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set integer bits matching a mask to 0 for a given key for a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 * If the data type does not match, CS_FIELD_INVALID_TYPE is returned.
 * If the key value has been locked, CS_FIELD_LOCKED is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  mask    mask associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_clear_key_int_bits(cs_field_t  *f,
                            int          key_id,
                            int          mask)
{
  int value = cs_field_get_key_int(f, key_id);

  value |= mask;
  value -= mask;

  int retval = cs_field_set_key_int(f, key_id, value);
  return retval;
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

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;
  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 'd')
      retval = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      if (kv->is_locked)
        retval = CS_FIELD_LOCKED;
      else {
        kv->val.v_double = value;
        kv->is_set = 1;
      }
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

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Field is not defined.", __func__);
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
        retval = kv->val.v_double;
      else if (kd->is_sub)
        retval = cs_field_get_key_double(f, kd->def_val.v_int);
      else
        retval = kd->def_val.v_double;
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

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;
  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 's')
      retval = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      if (kv->is_locked)
        retval = CS_FIELD_LOCKED;
      else {
        if (kv->is_set == 0)
          kv->val.v_p = NULL;
        BFT_REALLOC(kv->val.v_p, strlen(str) + 1, char);
        strcpy(kv->val.v_p, str);
        kv->is_set = 1;
      }
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

  if (f == NULL)
    return NULL;
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
        str = kv->val.v_p;
      else if (kd->is_sub)
        str = cs_field_get_key_str(f, kd->def_val.v_int);
      else
        str = kd->def_val.v_p;
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
 * \brief Assign a simple structure for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 * \param[in]  s       structure associated with key
 *
 * \return  0 in case of success, > 1 in case of error
 */
/*----------------------------------------------------------------------------*/

int
cs_field_set_key_struct(cs_field_t  *f,
                        int          key_id,
                        void        *s)
{
  int retval = CS_FIELD_OK;

  if (f == NULL)
    return CS_FIELD_INVALID_FIELD;
  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      retval = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 't')
      retval = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      if (kv->is_locked)
        retval = CS_FIELD_LOCKED;
      else {
        if (kv->is_set == 0)
          BFT_MALLOC(kv->val.v_p, kd->type_size, unsigned char);
        memcpy(kv->val.v_p, s, kd->type_size);
        kv->is_set = 1;
      }
    }
  }
  else
    retval = CS_FIELD_INVALID_KEY_ID;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a structure for a given key associated with a field.
 *
 * If the key id is not valid, or the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]   f       pointer to field structure
 * \param[in]   key_id  id of associated key
 * \param[out]  s       structure associated with key
 *
 * \return  pointer to structure associated with
 *          the key id for this field (same as s)
 */
/*----------------------------------------------------------------------------*/

const void *
cs_field_get_key_struct(const cs_field_t  *f,
                        const int          key_id,
                        void              *s)
{
  if (f == NULL)
    return NULL;
  assert(f->id >= 0 && f->id < _n_fields);

  int errcode = CS_FIELD_OK;

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 't')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const unsigned char *p = NULL;
      if (kv->is_set)
        p = kv->val.v_p;
      else if (kd->is_sub)
        p = cs_field_get_key_struct(f, kd->def_val.v_int, s);
      else
        p = kd->def_val.v_p;
      memcpy(s, p, kd->type_size);
      return s;
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
 * \brief Return a pointer to a simple structure for a given key to a field.
 *
 * If the key id is not valid, the value type or field category is not
 * compatible, or the structure has been locked, a fatal error is provoked.
 *
 * Note that using this function marks the field's value for this structure
 * as set, and if no values have been set yet, the structure is set to
 * default values.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  pointer to key structure in case of success, NULL in case of error
 */
/*----------------------------------------------------------------------------*/

void *
cs_field_get_key_struct_ptr(cs_field_t  *f,
                            int          key_id)
{
  if (f == NULL)
    return NULL;
  assert(f->id >= 0 && f->id < _n_fields);

  int errcode = CS_FIELD_OK;

  if (key_id > -1) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 't')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      void *p = NULL;
      if (kv->is_locked)
        errcode = CS_FIELD_LOCKED;
      else {
        if (kv->is_set == 0) {
          BFT_MALLOC(kv->val.v_p, kd->type_size, unsigned char);
          cs_field_get_key_struct(f, key_id, kv->val.v_p);
        }
        p = kv->val.v_p;
        kv->is_set = 1;
        return p;
      }
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
    else if (errcode == CS_FIELD_LOCKED)
      bft_error(__FILE__, __LINE__, 0,
                _("Field \"%s\" structure indicated by keyword %d (\"%s\")\n"
                  "has been locked.\n"
                  "use %s to access instead."),
                f->name, key_id, key, "cs_field_get_key_struct_const_ptr");
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a read-only pointer to a simple structure for a given key
 *        to a field.
 *
 * If the key id is not valid, the value type or field category is not
 * compatible, a fatal error is provoked.
 *
 * \param[in]  f       pointer to field structure
 * \param[in]  key_id  id of associated key
 *
 * \return  pointer to key structure in case of success, NULL in case of error
 */
/*----------------------------------------------------------------------------*/

const void *
cs_field_get_key_struct_const_ptr(const cs_field_t  *f,
                                  int                key_id)
{
  if (f == NULL)
    return NULL;
  assert(f->id >= 0 && f->id < _n_fields);

  int errcode = CS_FIELD_OK;

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else if (kd->type_id != 't')
      errcode = CS_FIELD_INVALID_TYPE;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const unsigned char *p = NULL;
      if (kv->is_set)
        p = kv->val.v_p;
      else if (kd->is_sub)
        p = cs_field_get_key_struct_const_ptr(f, kd->def_val.v_int);
      else
        p = kd->def_val.v_p;
      return p;
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
  int mask_id_end = 7;   /* _type_flag_*[CS_FIELD_CDO] */
  int mask_prev = 0;

  if (_n_fields == 0)
    return;

  /* Fields by category */

  for (cat_id = mask_id_start; cat_id < mask_id_end + 1; cat_id++) {

    size_t name_width = 24;

    /* First loop to determine name width */

    n_cat_fields = 0;

    for (i = 0; i < _n_fields; i++) {

      const cs_field_t *f = _fields[i];

      if (f->type & mask_prev)
        continue;

      size_t l = strlen(f->name);
      if (l > name_width)
        name_width = l;
    }

    if (name_width > 63)
      name_width = 63;

    /* Main loop */

    for (i = 0; i < _n_fields; i++) {

      char ilv_c = ' ';

      const cs_field_t *f = _fields[i];

      if (f->type & mask_prev)
        continue;

      if (cat_id == mask_id_end || f->type & _type_flag_mask[cat_id]) {

        char tmp_s[4][64] =  {"", "", "", ""};

        /* Print header for first field of each category */

        if (n_cat_fields == 0) {

          cs_log_strpad(tmp_s[0], _("Field"), name_width, 64);
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

          tmp_s[0][name_width] = '\0';
          tmp_s[1][4] = '\0';
          tmp_s[2][20] = '\0';
          tmp_s[3][4] = '\0';

          cs_log_printf(CS_LOG_SETUP, _("  %s %s %s %s ---------\n"),
                        tmp_s[0], tmp_s[1], tmp_s[2], tmp_s[3]);

        }

        /* Print field info */

        cs_log_strpad(tmp_s[0], f->name, name_width, 64);

        cs_log_strpad(tmp_s[1],
                      _(cs_mesh_location_get_name(f->location_id)),
                      20,
                      64);

        cs_log_printf(CS_LOG_SETUP,
                      "  %s %d %c  %s %-4d ",
                      tmp_s[0], f->dim, ilv_c,
                      tmp_s[1],
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

  } /* End of loop on categories */
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
                  "    Type:                       %d"), f->id, f->type);

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
  else
    cs_log_printf(CS_LOG_SETUP,
                  _("    Dimension:                  %d\n"),
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
      if (kd->log_id != 's')
        continue;
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      const char *key = cs_map_name_to_id_key(_key_map, i);
      if (kd->type_flag == 0 || (kd->type_flag & f->type)) {
        if (kd->type_id == 'i') {
          if (kv->is_set)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10d\n"),
                          key, kv->val.v_int);
          else if (log_keywords > 1)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10d (default)\n"),
                          key, kd->def_val.v_int);
        }
        else if (kd->type_id == 'd') {
          if (kv->is_set)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10.3g\n"),
                          key, kv->val.v_double);
          else if (log_keywords > 1)
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10.3g (default)\n"),
                          key, kd->def_val.v_double);
        }
        else if (kd->type_id == 's') {
          const char *s;
          if (kv->is_set) {
            s = kv->val.v_p;
            if (s == NULL)
              s = null_str;
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10s\n"), key, s);
          }
          else if (log_keywords > 1) {
            s = kd->def_val.v_p;
            if (s == NULL)
              s = null_str;
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10s (default)\n"),
                          key, s);
          }
        }
        else if (kd->type_id == 't') {
          const void *t;
          if (kv->is_set) {
            t = kv->val.v_p;
            if (kd->log_func != NULL) {
              cs_log_printf(CS_LOG_SETUP, _("      %-24s:\n"), key);
              kd->log_func(t);
            }
            else {
              cs_log_printf(CS_LOG_SETUP, _("      %-24s %-24p\n"), key, t);
            }
          }
          else if (log_keywords > 1) {
            t = kd->def_val.v_p;
            if (kd->log_func != NULL) {
              cs_log_printf(CS_LOG_SETUP, _("      %-24s: (default)\n"), key);
              kd->log_func(t);
            }
            else
              cs_log_printf(CS_LOG_SETUP, _("      %-24s %-24p (default)\n"),
                            key, t);
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

      f = _fields[i];

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

  cs_log_strpad(tmp_s[0], _("Key"), 24, 64);
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

  /* First loop on keys except structures */

  for (i = 0; i < _n_keys; i++) {

    int key_id = cs_map_name_to_id_try(_key_map,
                                       cs_map_name_to_id_key(_key_map, i));
    cs_field_key_def_t *kd = _key_defs + key_id;
    const char *key = cs_map_name_to_id_key(_key_map, i);

    if (kd->type_id == 'i') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12d integer %-4d "),
                    key, kd->def_val.v_int, key_id);
    }
    else if (kd->type_id == 'd') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12.3g real    %-4d "),
                    key, kd->def_val.v_double, key_id);
    }
    else if (kd->type_id == 's') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12s string  %-4d "),
                    key, (char *)(kd->def_val.v_p), key_id);
    }
    if (kd->type_id != 't') {
      if (kd->type_flag == 0)
        cs_log_printf(CS_LOG_SETUP, "0\n");
      else {
        cs_log_printf(CS_LOG_SETUP, "%-4d", kd->type_flag);
        _log_add_type_flag(kd->type_flag);
        cs_log_printf(CS_LOG_SETUP, "\n");
      }
    }

  } /* End of loop on keys */

  /* Second loop on keys structures */

  for (i = 0; i < _n_keys; i++) {

    int key_id = cs_map_name_to_id_try(_key_map,
                                       cs_map_name_to_id_key(_key_map, i));
    cs_field_key_def_t *kd = _key_defs + key_id;
    const char *key = cs_map_name_to_id_key(_key_map, i);

    if (kd->type_id == 't') {
      cs_log_printf(CS_LOG_SETUP,
                    _("  %-24s %-12s struct  %-4d "),
                    key, " ", key_id);

      if (kd->type_flag == 0)
        cs_log_printf(CS_LOG_SETUP, "0\n");
      else {
        cs_log_printf(CS_LOG_SETUP, "%-4d", kd->type_flag);
        _log_add_type_flag(kd->type_flag);
        cs_log_printf(CS_LOG_SETUP, "\n");
      }
    }
  } /* End of loop on keys */

  /* Third loop on keys structures for default values printing */

  char tmp_str[2][64] =  {"", ""};

  /* Print logging header */

  cs_log_strpad(tmp_str[0], _("Key"), 24, 64);
  cs_log_strpad(tmp_str[1], _("Default"), 12, 64);

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Default values for structure keys:\n"
                  "----------------------------------\n\n"));
  cs_log_printf(CS_LOG_SETUP, _("  %s %s Description\n"),
                tmp_str[0], tmp_str[1]);

  for (i = 0; i < 24; i++)
    tmp_str[0][i] = '-';
  tmp_str[0][24] = '\0';
  for (i = 0; i < 12; i++)
    tmp_str[1][i] = '-';
  tmp_str[1][12] = '\0';

  cs_log_printf(CS_LOG_SETUP,
                _("  %s %s -----------------------------------------\n"),
                tmp_str[0], tmp_str[1]);

  for (i = 0; i < _n_keys; i++) {

    int key_id = cs_map_name_to_id_try(_key_map,
                                       cs_map_name_to_id_key(_key_map, i));
    cs_field_key_def_t *kd = _key_defs + key_id;
    const void *t;

    if (kd->type_id == 't') {
      t = kd->def_val.v_p;

      if (kd->log_func_default != NULL)
        kd->log_func_default(t);
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

  /* First loop to determine field width */

  size_t name_width = 24;

  for (i = 0; i < _n_fields; i++) {
    const cs_field_t *f = _fields[i];
    size_t l = strlen(f->name);
    if (l > name_width)
      name_width = l;
  }
  if (name_width > 63)
    name_width = 63;

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

    /* Main loop on fields */

    for (i = 0; i < _n_fields; i++) {

      const cs_field_t *f = _fields[i];

      if (f->type & mask_prev)
        continue;

      if (cat_id == mask_id_end || f->type & _type_flag_mask[cat_id]) {

        char name_s[64] =  "";
        cs_log_strpad(name_s, f->name, name_width, 64);

        cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);

        if (kd->type_flag == 0 || (kd->type_flag & f->type)) {
          if (kd->type_id == 'i') {
            if (kv->is_set)
              cs_log_printf(CS_LOG_SETUP, "    %s %d\n",
                            name_s, kv->val.v_int);
            else if (log_defaults)
              cs_log_printf(CS_LOG_SETUP, _("    %s %-10d (default)\n"),
                            name_s, kd->def_val.v_int);
          }
          else if (kd->type_id == 'd') {
            if (kv->is_set)
              cs_log_printf(CS_LOG_SETUP, _("    %s %-10.3g\n"),
                          name_s, kv->val.v_double);
            else if (log_defaults)
              cs_log_printf(CS_LOG_SETUP, _("    %s %-10.3g (default)\n"),
                            name_s, kd->def_val.v_double);
          }
          else if (kd->type_id == 's') {
            const char *s;
            if (kv->is_set) {
              s = kv->val.v_p;
              if (s == NULL)
                s = null_str;
              cs_log_printf(CS_LOG_SETUP, _("    %s %s\n"), name_s, s);
            }
            else if (log_defaults) {
              s = kd->def_val.v_p;
              if (s == NULL)
                s = null_str;
              cs_log_printf(CS_LOG_SETUP, _("    %s %-10s (default)\n"),
                            name_s, s);
            }
          }
          else if (kd->type_id == 't') {
            if (kv->is_set) {
              cs_log_printf(CS_LOG_SETUP, _("\n    %s\n"), name_s);
              if (kd->log_func != NULL)
                kd->log_func(kv->val.v_p);
            }
            else if (log_defaults) {
              cs_log_printf(CS_LOG_SETUP, _("\n    %s (default)\n"), name_s);
              if (kd->log_func != NULL)
                kd->log_func(kd->def_val.v_p);
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
 *   "label"        (string)
 *   "log"          (integer)
 *   "post_vis"     (integer)
 *   "coupled"      (integer, restricted to CS_FIELD_VARIABLE)
 *   "moment_id"    (integer, restricted to
 *                   CS_FIELD_ACCUMULATOR | CS_FIELD_POSTPROCESS);
 *
 * A recommended practice for different submodules would be to use
 * "cs_<module>_key_init() functions to define keys specific to those modules.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_define_keys_base(void)
{
  cs_field_define_key_str("label", NULL, 0);
  _k_label = cs_field_key_id("label");

  cs_field_define_key_int("log", 0, 0);
  cs_field_define_key_int("post_vis", 0, 0);
  cs_field_define_key_int("coupled", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("moment_id", -1,
                          CS_FIELD_ACCUMULATOR | CS_FIELD_POSTPROCESS);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a label associated with a field.
 *
 * If the "label" key has been set for this field, its associated string
 * is returned. Otherwise, the field's name is returned.
 *
 * \param[in]  f       pointer to field structure
 *
 * \return  pointer to character string associated with label for this field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_field_get_label(const cs_field_t  *f)
{
  const char *label = cs_field_get_key_str(f, _k_label);

  if (label == NULL)
    label = f->name;

  return label;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
