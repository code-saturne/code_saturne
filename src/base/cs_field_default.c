/*============================================================================
 * Field utility functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_equation_param.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field.h"
#include "cs_field_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field_default.c
        Field utility functions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Dimensions and pointers (mappable to legacy Fortran array,
   which is why this is not simply allocated locally for each field) */

static cs_lnum_t    _n_vars_bc = 0;
static cs_lnum_t    _n_b_faces = 0;
static int         *_icodcl = NULL;
static cs_real_t   *_rcodcl = NULL;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_field_update_bcs_ptr(int          dim_icodcl[2],
                          int          dim_rcodcl[3],
                          int        **p_icodcl,
                          cs_real_t  **p_rcodcl);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update field icodcl/rcodcl map and pointers to matching base arrays.
 *
 * This function is intended for use by Fortran wrappers.
 *
 * parameters:
 *   dim_icodcl --> dimensions of icodcl array
 *   dim_rcodcl --> dimensions of rcodcl array
 *   p_icodcl   --> pointer to icodcl base array
 *   p_rcodcl   --> pointer to rcodcl base array
 *----------------------------------------------------------------------------*/

void
cs_f_field_update_bcs_ptr(int          dim_icodcl[2],
                          int          dim_rcodcl[3],
                          int        **p_icodcl,
                          cs_real_t  **p_rcodcl)
{
  cs_field_build_bc_codes_all();

  dim_icodcl[0] = _n_b_faces;
  dim_icodcl[1] = _n_vars_bc;

  dim_rcodcl[0] = _n_b_faces;
  dim_rcodcl[1] = _n_vars_bc;
  dim_rcodcl[2] = 3;

  *p_icodcl = _icodcl;
  *p_rcodcl = _rcodcl;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a field shared between CDO and legacy schemes. This field is
 *         related to a general solved variable, with default options.
 *
 * \param[in]  name          field name
 * \param[in]  label         field default label, or empty
 * \param[in]  location_id   id of associated location
 * \param[in]  dim           field dimension
 * \param[in]  has_previous  no if lower than 1
 *
 * \return  newly defined field id
 */
/*----------------------------------------------------------------------------*/

int
cs_variable_cdo_field_create(const char  *name,
                             const char  *label,
                             int          location_id,
                             int          dim,
                             int          has_previous)
{
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE | CS_FIELD_CDO;

  /* If cmp_id > -1 then this is an existing field. This situation may happen
     with CDO field if a previous creation was made in the Fortran part */

  int cmp_id = cs_field_id_by_name(name);

  /* Conversion from int to bool (done in C to avoid spurious behavior with
     a boolean variable defined in the FORTRAN part */

  bool  previous = (has_previous < 1) ? false : true;
  cs_field_t *f = cs_field_find_or_create(name,
                                          field_type,
                                          location_id,
                                          dim,
                                          previous);  /* has_previous */

  if (cmp_id == -1) { /* Do not modify post_flag if the field was previously
                         created */

    const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
    cs_field_set_key_int(f, cs_field_key_id("log"), 1);
    cs_field_set_key_int(f, cs_field_key_id("post_vis"), post_flag);

    if (label != NULL) {
      if (strlen(label) > 0)
        cs_field_set_key_str(f, cs_field_key_id("label"), label);
    }

  }

  return f->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add field defining a general solved variable, with default options.
 *
 * \param[in]  name         field name
 * \param[in]  label        field default label, or empty
 * \param[in]  location_id  id of associated location
 * \param[in]  dim          field dimension
 *
 * \return  newly defined field id
 */
/*----------------------------------------------------------------------------*/

int
cs_variable_field_create(const char  *name,
                         const char  *label,
                         int          location_id,
                         int          dim)
{
  int field_type = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

  int cmp_id = cs_field_id_by_name(name);

  if (cmp_id > -1)
    bft_error(__FILE__, __LINE__, 0,
              _("Error defining variable \"%s\";\n"
                "this name is already reserved for field with id %d."),
              name, cmp_id);

  cs_field_t *f = cs_field_create(name,
                                  field_type,
                                  location_id,
                                  dim,
                                  true);  /* has_previous */

  const int post_flag = CS_POST_ON_LOCATION | CS_POST_MONITOR;
  cs_field_set_key_int(f, cs_field_key_id("log"), 1);
  cs_field_set_key_int(f, cs_field_key_id("post_vis"), post_flag);

  if (label != NULL) {
    if (strlen(label) > 0)
      cs_field_set_key_str(f, cs_field_key_id("label"), label);
  }

  if (dim > 1)
    cs_field_set_key_int(f, cs_field_key_id("coupled"), 1);

  return f->id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Access a field's equation parameters.
 *
 * If the equation parameters were never initialized, they will be initialized
 * based on the current defaults.
 *
 * If the field does not have associated equation paremeters (i.e. is not
 * a variable field or is a CDO field (which is referenced by but does not
 * directly reference equations), NULL is returned.
 *
 * \param[in, out]  f  pointer to associated field
 *
 * \return  pointer to field's equation parameters, or NULL
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_field_get_equation_param(cs_field_t  *f)
{
  cs_equation_param_t *eqp = NULL;

  static int k_id = -1;
  if (k_id < 0)
    k_id = cs_field_key_id_try("var_cal_opt");

  if (k_id >= 0 && (f->type & CS_FIELD_VARIABLE))
    eqp = cs_field_get_key_struct_ptr(f, k_id);

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Access a field's equation parameters for read only.
 *
 * If the equation parameters were never initialized, the current default
 * parameters will be returned instead.
 *
 * If the field does not have associated equation parameters (i.e. is not
 * a variable field or is a CDO field (which is referenced by but does not
 * directly reference equations), NULL is returned.
 *
 * \param[in]  f  pointer to associated field
 *
 * \return  const-qualified pointer to field's equation parameters, or NULL
 */
/*----------------------------------------------------------------------------*/

const cs_equation_param_t *
cs_field_get_equation_param_const(const cs_field_t  *f)
{
  const cs_equation_param_t *eqp = NULL;

  static int k_id = -1;
  if (k_id < 0)
    k_id = cs_field_key_id_try("var_cal_opt");

  if (k_id >= 0 && (f->type & CS_FIELD_VARIABLE))
    eqp = cs_field_get_key_struct_const_ptr(f, k_id);

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and map boundary condition coefficients for all
 *        variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_build_bc_codes_all(void)
{
  const cs_lnum_t n_b_faces
    = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES)[0];

  /* Determine number of solved variables */

  cs_lnum_t n_vars = 0;

  const int kv = cs_field_key_id("variable_id");
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    int var_id_p1 = (f->type & CS_FIELD_VARIABLE) ?
      cs_field_get_key_int(f, kv) : 0;

    if (var_id_p1 > 0)
      var_id_p1 += f->dim - 1;

    if (var_id_p1 > n_vars)
      n_vars = var_id_p1;

  }

  /* Allocate or remap only if needed */

  if (n_vars == _n_vars_bc && n_b_faces == _n_b_faces) {
    if (   (_icodcl != NULL && n_b_faces > 0)
        || (_icodcl == NULL && n_b_faces == 0))
      return;
  }

  /* Update and map otherwise */

  _n_vars_bc = n_vars;
  _n_b_faces = n_b_faces;

  BFT_REALLOC(_icodcl, n_vars*n_b_faces, int);
  BFT_REALLOC(_rcodcl, 3*n_vars*n_b_faces, cs_real_t);

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    int var_id = (f->type & CS_FIELD_VARIABLE) ?
      cs_field_get_key_int(f, kv) -1 : -1;

    if (var_id > -1) {

      int * icodcl  = _icodcl + n_b_faces*var_id;
      cs_real_t *rcodcl1 = _rcodcl + n_b_faces*var_id;
      cs_real_t *rcodcl2 = _rcodcl + n_b_faces*(n_vars+var_id);
      cs_real_t *rcodcl3 = _rcodcl + n_b_faces*(2*n_vars+var_id);

      if (f->bc_coeffs != NULL) {

        f->bc_coeffs->icodcl  = icodcl;
        f->bc_coeffs->rcodcl1 = rcodcl1;
        f->bc_coeffs->rcodcl2 = rcodcl2;
        f->bc_coeffs->rcodcl3 = rcodcl3;

        /* For multi-dimensional arrays, access using

           f->bcs->icodcl[coo_id*n_b_faces + face_id]

           and equivalent for icodcl, rcodc1, rcodcl2, rcodcl3. */

      }

      /* Initialize icodcl and rcodcl values to defaults
         (even if they are not mapped to C, to avoid
         issues in some Fortran initialization loops
         which do not do the appropriate checks for
         variable type. */

      for (cs_lnum_t coo_id = 0; coo_id < f->dim; coo_id++) {
        cs_lnum_t s_id = n_b_faces * coo_id;
        cs_lnum_t e_id = n_b_faces * (coo_id+1);

        for (cs_lnum_t i = s_id; i < e_id; i++) {
          icodcl[i] = 0;
          rcodcl1[i] = cs_math_infinite_r;
          rcodcl2[i] = cs_math_infinite_r;
          rcodcl3[i] = 0;
        }
      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Deallocate and unmap boundary condition coefficients for all
 *        variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_free_bc_codes_all(void)
{
  _n_vars_bc = 0;
  _n_b_faces = 0;

  const int kv = cs_field_key_id("variable_id");
  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    int var_id = (f->type & CS_FIELD_VARIABLE) ?
      cs_field_get_key_int(f, kv) : -1;

    if (var_id > -1 && f->bc_coeffs != NULL) {
      f->bc_coeffs->icodcl = NULL;
      f->bc_coeffs->rcodcl1 = NULL;
      f->bc_coeffs->rcodcl2 = NULL;
      f->bc_coeffs->rcodcl3 = NULL;
    }

  }

  BFT_FREE(_icodcl);
  BFT_FREE(_rcodcl);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
