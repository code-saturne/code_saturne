/*============================================================================
 * Interpolation function handling.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field_default.h"
#include "cs_gradient.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_interpolate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_interpolate.c
        Interpolation function handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate values defined on a mesh location at a given set of
 *        points using a P0 interpolation.
 *
 * This function allows unlocated points (with point_location < 0),
 * to which the value 0 is assigned.
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 */
/*----------------------------------------------------------------------------*/

void
cs_interpolate_from_location_p0(void                *input,
                                cs_datatype_t        datatype,
                                int                  val_dim,
                                cs_lnum_t            n_points,
                                const cs_lnum_t      point_location[],
                                const cs_real_3_t    point_coords[],
                                const void          *location_vals,
                                void                *point_vals)
{
  CS_UNUSED(input);
  CS_UNUSED(point_coords);

  switch(datatype) {

  case CS_INT32:
    {
      const int32_t *l_vals = (const int32_t *)location_vals;
      int32_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  case CS_INT64:
    {
      const int64_t *l_vals = (const int64_t *)location_vals;
      int64_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  case CS_DOUBLE:
    {
      const cs_real_t *l_vals = (const cs_real_t *)location_vals;
      cs_real_t *p_vals = point_vals;
      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t e_id = point_location[i];
        if (e_id > -1) {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = l_vals[e_id*val_dim + j];
        }
        else {
          for (cs_lnum_t j = 0; j < val_dim; j++)
            p_vals[i*val_dim + j] = 0;
        }
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Function %s does not currently handle %s data type."),
              __func__, cs_datatype_name[datatype]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate values defined on a mesh location at a given set of
 *        points using a P1 interpolation.
 *
 * This function assumes the input is a field name. If no field matches
 * the name, a "generic" interpolation (assuming homogeneous Neumann boundary
 * conditions) is used).
 *
 * The P1 interpolation is based on a local least-squares gradient,
 * so it is assumed that ghost values of the localtion_vals array
 * are synchronized if present.
 *
 * If the field's boundary values (i.e. associated field) are known,
 * they are used in the interpolation. Otherwise, if boundary conditions
 * are defined, they are used. When neither boundary values nor boundary
 * conditions are known, homogeneous Neumann boundary conditions are assumed.
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 */
/*----------------------------------------------------------------------------*/

void
cs_interpolate_from_location_p1(void                *input,
                                cs_datatype_t        datatype,
                                int                  val_dim,
                                cs_lnum_t            n_points,
                                const cs_lnum_t      point_location[],
                                const cs_real_3_t    point_coords[],
                                const void          *location_vals,
                                void                *point_vals)
{
  /* If used with a non-real argument type, use P0 interpolation */

  if (   datatype != CS_REAL_TYPE
      || (val_dim != 1 && val_dim != 3 && val_dim != 6)) {
    cs_interpolate_from_location_p0(input,
                                    datatype,
                                    val_dim,
                                    n_points,
                                    point_location,
                                    point_coords,
                                    location_vals,
                                    point_vals);
    return;
  }

  /* Main usage */

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)fvq->cell_cen;

  cs_halo_type_t halo_type
    = (m->cell_cells_idx != NULL) ? CS_HALO_EXTENDED : CS_HALO_STANDARD;

  cs_field_t *f = NULL;
  const cs_real_t *bc_coeff_a = NULL;
  const cs_real_t *bc_coeff_b = NULL;

  if (input != NULL) {
    const char *name = input;
    f = cs_field_by_name_try(name);
    if (f != NULL) {
      int kbf = cs_field_key_id_try("boundary_value_id");
      int bf_id = cs_field_get_key_int(f, kbf);
      if (bf_id > -1) {
        const cs_field_t *bf = cs_field_by_id(bf_id);
        bc_coeff_a = bf->val;
      }
      else if (f->bc_coeffs != NULL) {
        bc_coeff_a = f->bc_coeffs->a;
        bc_coeff_b = f->bc_coeffs->b;
        if (f->dim > 1 && f->type & CS_FIELD_VARIABLE) {
          int coupled = 0;
          int coupled_key_id = cs_field_key_id_try("coupled");
          if (coupled_key_id > -1)
            coupled = cs_field_get_key_int(f, coupled_key_id);
          if (coupled == 0) {   /* not handled in this case */
            bc_coeff_a = NULL;
            bc_coeff_b = NULL;
          }
        }
      }
      if (f->type & CS_FIELD_VARIABLE) {
        const cs_equation_param_t *eqp
          = cs_field_get_equation_param_const(f);
        cs_gradient_type_t gradient_type = CS_GRADIENT_LSQ;
        cs_gradient_type_by_imrgra(eqp->imrgra,
                                   &gradient_type,
                                   &halo_type);
      }
    }
  }

  switch(val_dim) {
  case 1:
    {
      const cs_real_t *c_vals = (const cs_real_t *)location_vals;
      cs_real_t *p_vals = (cs_real_t *)point_vals;

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t grad[3];
          cs_gradient_scalar_cell(m,
                                  fvq,
                                  c_id,
                                  halo_type,
                                  bc_coeff_a,
                                  bc_coeff_b,
                                  c_vals,
                                  NULL,
                                  grad);

          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          p_vals[i] = c_vals[c_id] + grad[0]*d[0]
                                   + grad[1]*d[1]
                                   + grad[2]*d[2];
        }
        else
          p_vals[i] = 0;
      }

    }
    break;

  case 3:
    {
      const cs_real_3_t *c_vals = (const cs_real_3_t *)location_vals;
      cs_real_3_t *p_vals = (cs_real_3_t *)point_vals;

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t grad[3][3];
          cs_gradient_vector_cell(m,
                                  fvq,
                                  c_id,
                                  halo_type,
                                  (const cs_real_3_t *)bc_coeff_a,
                                  (const cs_real_33_t *)bc_coeff_b,
                                  c_vals,
                                  NULL,
                                  grad);

          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          for (cs_lnum_t j = 0; j < 3; j++) {
            p_vals[i][j] = c_vals[c_id][j] + grad[j][0]*d[0]
                                           + grad[j][1]*d[1]
                                           + grad[j][2]*d[2];
          }
        }
        else {
          for (cs_lnum_t j = 0; j < 6; j++)
            p_vals[i][j] = 0;
        }
      }
    }
    break;

  case 6:
    {
      const cs_real_6_t *c_vals = (const cs_real_6_t *)location_vals;
      cs_real_6_t *p_vals = (cs_real_6_t *)point_vals;

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t grad[6][3];
          cs_gradient_tensor_cell(m,
                                  fvq,
                                  c_id,
                                  halo_type,
                                  (const cs_real_6_t *)bc_coeff_a,
                                  (const cs_real_66_t *)bc_coeff_b,
                                  c_vals,
                                  NULL,
                                  grad);
          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          for (cs_lnum_t j = 0; j < 6; j++) {
            p_vals[i][j] = c_vals[c_id][j] + grad[j][0]*d[0]
                                           + grad[j][1]*d[1]
                                           + grad[j][2]*d[2];
          }
        }
        else {
          for (cs_lnum_t j = 0; j < 6; j++)
            p_vals[i][j] = 0;
        }
      }

    }
    break;

  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
