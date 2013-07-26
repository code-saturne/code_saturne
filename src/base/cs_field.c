/*============================================================================
 * Field management.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

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
 * Additional doxygen documentation
 *============================================================================*/
/*!
  \file cs_field.c
        Field management.

  \struct cs_field_bc_coeffs_t

  \brief Field boundary condition descriptor (for variables)

  \var cs_field_bc_coeffs_t::location_id
       Id of matching location

  \var cs_field_bc_coeffs_t::a
       Explicit coefficient
  \var cs_field_bc_coeffs_t::b
       Implicit coefficient
  \var cs_field_bc_coeffs_t::af
       Explicit coefficient for flux
  \var cs_field_bc_coeffs_t::bf
       Implicit coefficient for flux

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
  \var  cs_field_t::interleaved
        are field value arrays interleaved ? (recommended for new developments,
        but mapped legacy fields may be non-interleaved)
  \var  cs_field_t::location_id
        Id of matching mesh location
  \var  cs_field_t::n_time_vals
        Number of time values (1 or 2)
  \var  cs_field_t::val
        For each active location, pointer to matching values array
  \var  cs_field_t::val_pre
        For each active location, pointer to matching previous values array
        (only if n_time_vals > 1)
  \var  cs_field_t::bc_coeffs
        Boundary condition coefficients, for variable type fields
  \var  cs_field_t::is_owner
        Ownership flag for values and boundary coefficients
*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Field key definitions */

typedef struct {

  unsigned char               def_val[8];   /* Default value container (int,
                                               double, or pointer to string
                                               or structure) */
  cs_field_log_key_struct_t  *log_func;     /* print function for structure */
  size_t                      type_size;    /* Type length for added types
                                               (0 for 'i', 'd', or 's') */
  int                         type_flag;    /* Field type flag */
  char                        type_id;      /* i: int; d: double; s: str;
                                               t: type */

  bool                        is_sub;       /* Indicate if the key is a sub-key
                                               (in which case def_val contains
                                               the parent key id */

} cs_field_key_def_t;

/* Field key value structures */

typedef struct {

  unsigned char      val[8];       /* Value container (int, double,
                                      or pointer) */
  bool               is_set;       /* Has this key been set for the
                                      present field ? */

} cs_field_key_val_t;

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
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

int
cs_f_field_n_fields(void);

int
cs_f_field_id_by_name(const char *name);

int
cs_f_field_id_by_name_try(const char *name);

void
cs_f_field_get_name(int           id,
                    int           name_max,
                    const char  **name,
                    int          *name_len);

void
cs_f_field_get_dimension(int           id,
                         int           dim[2]);

void
cs_f_field_var_ptr_by_id(int          id,
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


/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Create a field descriptor.
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
  f->interleaved = true;
  if (f->dim > 1 && interleaved == false)
    f->interleaved = false;
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

/*----------------------------------------------------------------------------
 * Free structure associated to a key.
 *----------------------------------------------------------------------------*/

static void
_cs_field_free_struct(void)
{
  int key_id, f_id;

  for (key_id = 0; key_id < _n_keys; key_id++) {

    cs_field_key_def_t *kd = _key_defs + key_id;

    if (kd->type_id == 't') {
      for (f_id = 0; f_id < _n_fields; f_id++) {
        cs_field_key_val_t *kv = _key_vals + (f_id*_n_keys_max + key_id);
        unsigned char **p = (unsigned char **)(kv->val);
        BFT_FREE(*p);
      }
    }

  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
 * Return the dimension of a field defined by its id.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id  <-- field id
 *   dim <-- field dimension and interleave flag
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_dimension(int  id,
                         int  dim[2])
{
  const cs_field_t *f = cs_field_by_id(id);

  dim[0] = f->dim;
  dim[1] = (f->interleaved) ? 1 : 0;
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
 * Return a pointer to a field's variable values
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   id           <-- field id
 *   pointer_type <-- 1: var; 2: var_p;
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

  dim[1] = 0;
  dim[2] = 0;
  *p = NULL;

  if (pointer_type == 1 || pointer_type == 2) {

    const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(f->location_id);
    cs_lnum_t _n_elts = n_elts[2];

    if (pointer_type == 1)
      *p = f->val;
    else
      *p = f->val_pre;

    if (*p == NULL) /* Adjust dimensions to assist Fortran bounds-checking */
      _n_elts = 0;

    if (f->dim == 1)
      dim[0] = _n_elts;
    else if (f->interleaved) {
      dim[0] = f->dim;
      dim[1] = _n_elts;
      cur_p_rank = 2;
    }
    else {
      dim[0] = _n_elts;
      dim[1] = f->dim;
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
 *   pointer_type <-- 1: bc_coeffs->a;   2: bc_coeffs->b
 *                    3: bc_coeffs->af;  4: bc_coeffs->bf
 *                    5: bc_coeffs->ad;  6: bc_coeffs->bd
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

  dim[1] = 0;
  dim[2] = 0;
  dim[3] = 0;
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

  if (f->type & CS_FIELD_VARIABLE) {

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

    if (*p == NULL) /* Adjust dimensions to assist Fortran bounds-checking */
      _n_elts = 0;

    if (f->dim == 1)
      dim[0] = _n_elts;

    else {

      int coupled = 0;
      int coupled_key_id = cs_field_key_id_try("coupled");

      if (coupled_key_id > -1)
        coupled = cs_field_get_key_int(f, coupled_key_id);

      if (coupled) {

        if (pointer_type == 1 || pointer_type == 3 || pointer_type == 5) {
          dim[0] = f->dim;
          dim[1] = _n_elts;
          cur_p_rank = 2;
        }
        else { /* if (pointer_type == 2 || pointer_type == 4 || pointer_type == 6) */
          dim[0] = f->dim;
          dim[1] = f->dim;
          dim[2] = _n_elts;
          cur_p_rank = 3;
        }

      }
      else { /* uncoupled */

        if (f->interleaved) {
          dim[0] = f->dim;
          dim[1] = _n_elts;
        }
        else {
          dim[0] = _n_elts;
          dim[1] = f->dim;
        }
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

  *str_len = strlen(*str);

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
 *
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
 *
 *----------------------------------------------------------------------------*/

void
cs_f_field_get_key_struct(int    f_id,
                          int    k_id,
                          void  *k_value)
{
  const cs_field_t *f = cs_field_by_id(f_id);

  cs_field_get_key_struct(f, k_id, k_value);
}

/*! \endcond (end ignore by Doxygen) */

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
 * is the norm, and implicit b and bf coefficient arrays are arrays of
 * block matrices, not vectors, so the number of entries for each boundary
 * face is dim*dim instead of dim.
 *
 * \param[in, out]  f             pointer to field structure
 * \param[in]       have_flux_bc  if true, flux bc coefficients (af and bf)
 *                                are added
 * \param[in]       have_mom_bc   if true, div BC coefficients (ad and bd)
 *                                are added
 */
/*----------------------------------------------------------------------------*/

void
cs_field_allocate_bc_coeffs(cs_field_t  *f,
                            bool         have_flux_bc,
                            bool         have_mom_bc)
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

      if (have_mom_bc) {
        BFT_MALLOC(f->bc_coeffs->ad,  n_elts[0]*a_mult, cs_real_t);
        BFT_MALLOC(f->bc_coeffs->bd,  n_elts[0]*b_mult, cs_real_t);
      }
      else {
        f->bc_coeffs->ad = NULL;
        f->bc_coeffs->bd = NULL;
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

      if (have_mom_bc) {
        BFT_REALLOC(f->bc_coeffs->ad,  n_elts[0]*a_mult, cs_real_t);
        BFT_REALLOC(f->bc_coeffs->bd,  n_elts[0]*b_mult, cs_real_t);
      }
      else {
        BFT_FREE(f->bc_coeffs->ad);
        BFT_FREE(f->bc_coeffs->bd);
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
 * For multidimensional fields, arrays are assumed to have the same
 * interleaving behavior as the field, unless components are coupled.
 *
 * For multidimensional fields with coupled components, interleaving
 * is the norm, and implicit b and bf coefficient arrays are arrays of
 * block matrices, not vectors, so the number of entries for each boundary
 * face is dim*dim instead of dim.
 *
 * \param[in, out]  f             pointer to field structure
 * \param[in]       have_flux_bc  if true, flux bc coefficients (af and bf)
 *                                are initialized
 * \param[in]       have_mom_bc   if true, div BC coefficients (ad and bd)
 *                                are initialized
 */
/*----------------------------------------------------------------------------*/

void
cs_field_init_bc_coeffs(cs_field_t  *f,
                        bool         have_flux_bc,
                        bool         have_mom_bc)
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

      if (have_flux_bc)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->af[ifac] = 0.;
          f->bc_coeffs->bf[ifac] = 0.;
        }

      if (have_mom_bc)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ad[ifac] = 0.;
          f->bc_coeffs->bd[ifac] = 1.;
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
        f->bc_coeffs->b[ifac*dim*dim + 3] = 1.;
        f->bc_coeffs->b[ifac*dim*dim + 4] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 5] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 6] = 1.;
        f->bc_coeffs->b[ifac*dim*dim + 7] = 0.;
        f->bc_coeffs->b[ifac*dim*dim + 8] = 0.;
      }

      if (have_flux_bc)
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

      if (have_mom_bc)
        for (ifac = 0; ifac < n_elts[0]; ifac++) {
          f->bc_coeffs->ad[ifac*dim] = 0.;
          f->bc_coeffs->ad[ifac*dim + 1] = 0.;
          f->bc_coeffs->ad[ifac*dim + 2] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim] = 1.;
          f->bc_coeffs->bd[ifac*dim*dim + 1] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 2] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 3] = 1.;
          f->bc_coeffs->bd[ifac*dim*dim + 4] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 5] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 6] = 1.;
          f->bc_coeffs->bd[ifac*dim*dim + 7] = 0.;
          f->bc_coeffs->bd[ifac*dim*dim + 8] = 0.;
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
        BFT_FREE(f->bc_coeffs->ad);
        BFT_FREE(f->bc_coeffs->bd);
      }
      BFT_FREE(f->bc_coeffs);
    }
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
  kd->log_func = NULL;
  kd->type_size = 0;
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
  kd->log_func = NULL;
  kd->type_size = 0;
  kd->type_flag = type_flag;
  kd->type_id = 'd';
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
  kd->log_func = NULL;
  kd->type_size = 0;
  kd->type_flag = type_flag;
  kd->type_id = 's';
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
 * \param[in]  name            key name
 * \param[in]  default_value   pointer to default value associated with key
 * \param[in]  log_func        pointer to logging function
 * \param[in]  size            sizeof structure
 * \param[in]  type_flag       mask associated with field types with which
 *                             the key may be associated, or 0
 *
 * \return  id associated with key
 */
/*----------------------------------------------------------------------------*/

int
cs_field_define_key_struct(const char                 *name,
                           const void                 *default_value,
                           cs_field_log_key_struct_t  *log_func,
                           size_t                      size,
                           int                         type_flag)
{
  int n_keys_init = _n_keys;

  int key_id = _find_or_add_key(name);

  cs_field_key_def_t *kd = _key_defs + key_id;

  unsigned char **def_val = (unsigned char **)(kd->def_val);

  /* Free possible previous allocation */
  if (n_keys_init == _n_keys)
    BFT_FREE(*def_val);

  if (default_value != NULL) {
    BFT_MALLOC(*def_val, size, unsigned char);
    memcpy(*def_val, default_value, size);
  }
  else
    *def_val = NULL;
  kd->log_func = log_func;
  kd->type_size = size;
  kd->type_flag = type_flag;
  kd->type_id = 't';
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
  int key_id;
  for (key_id = 0; key_id < _n_keys; key_id++) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    if (kd->type_id == 't') {
      unsigned char **def_val = (unsigned char **)(kd->def_val);
      BFT_FREE(*def_val);
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
  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

  if (key_id > -1 && key_id < _n_keys) {
    cs_field_key_def_t *kd = _key_defs + key_id;
    assert(key_id < _n_keys);
    if (kd->type_flag != 0 && !(f->type & kd->type_flag))
      errcode = CS_FIELD_INVALID_CATEGORY;
    else {
      cs_field_key_val_t *kv = _key_vals + (f->id*_n_keys_max + key_id);
      bool retval = false;
      if (kv->is_set)
        retval = true;
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
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Field keyword with id %d is not defined."),
                key_id);
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a integer value for a given key to a field.
 *
 * If the key id is not valid, CS_FIELD_INVALID_KEY_ID is returned.
 * If the field category is not compatible with the key (as defined
 * by its type flag), CS_FIELD_INVALID_CATEGORY is returned.
 * If the data type does not match, CS_FIELD_INVALID_TYPE is returned.
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
    else if (kd->type_id != 'i')
      retval = CS_FIELD_INVALID_TYPE;
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
    else if (kd->type_id != 'd')
      retval = CS_FIELD_INVALID_TYPE;
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
    else if (kd->type_id != 's')
      retval = CS_FIELD_INVALID_TYPE;
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
      unsigned char **_val = (unsigned char **)(kv->val);
      if (kv->is_set == false) {
                    BFT_MALLOC(*_val, kd->type_size, unsigned char);
                        }
      memcpy(*_val, s, kd->type_size);
      kv->is_set = true;
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
  int errcode = CS_FIELD_OK;

  assert(f->id >= 0 && f->id < _n_fields);

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
        p = *((const unsigned char **)(kv->val));
      else if (kd->is_sub)
        p = cs_field_get_key_struct(f, *((int *)(kd->def_val)), s);
      else
        p = *((const unsigned char **)(kd->def_val));
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

    size_t name_width = 24;

    /* First loop to determine name width */

    n_cat_fields = 0;

    for (i = 0; i < _n_fields; i++) {

      const cs_field_t *f = _fields + i;

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

      const cs_field_t *f = _fields + i;

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

        if (f->interleaved == false)
          ilv_c = 'n';

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
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10s\n"), key, s);
          }
          else if (log_keywords > 1) {
            s = *(const char **)(kd->def_val);
            if (s == NULL)
              s = null_str;
            cs_log_printf(CS_LOG_SETUP, _("      %-24s %-10s (default)\n"),
                          key, s);
          }
        }
        else if (kd->type_id == 't') {
          const void *t;
          if (kv->is_set == true) {
            t = (const void *)(kv->val);
            if (kd->log_func != NULL) {
              cs_log_printf(CS_LOG_SETUP, _("      %-24s:\n"), key);
              kd->log_func(t);
            }
            else {
              cs_log_printf(CS_LOG_SETUP, _("      %-24s %-24p\n"), key, t);
            }
          }
          else if (log_keywords > 1) {
            t = (const void *)(kd->def_val);
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

  /* First loop on keys execpt structures */

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
    const void *t;

    if (kd->type_id == 't') {
      t = *(const void **)(kd->def_val);

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

      if (kd->log_func != NULL)
        kd->log_func(t);
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
          else if (kd->type_id == 't') {
            const void *t;
            if (kv->is_set == true) {
              t = *(const void **)(kv->val);
              cs_log_printf(CS_LOG_SETUP, _("    %-24s\n"),
                            f->name);
              if (kd->log_func != NULL)
                kd->log_func(t);
            }
            else if (log_defaults) {
              t = *(const void **)(kd->def_val);
              cs_log_printf(CS_LOG_SETUP, _("    %-24s (default)\n"),
                            f->name);
              if (kd->log_func != NULL)
                kd->log_func(t);
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
 *   "log"       (integer)
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
  cs_field_define_key_int("log", 0, 0);
  cs_field_define_key_int("coupled", 0, CS_FIELD_VARIABLE);
  cs_field_define_key_int("moment_dt", -1, CS_FIELD_PROPERTY);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
