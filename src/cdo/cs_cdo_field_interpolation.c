/*============================================================================
 * Functions to handle field interpolation with CDO schemes
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_array.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_equation.h"
#include "cs_equation_priv.h"
#include "cs_timer_stats.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_field_interpolation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDO_FIELD_INTERPOLATION_DBG   0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

cs_flag_t       _field_interpolation_flag = 0;
cs_equation_t  *_field_interpolation_scalar_c2v_eq = NULL;
cs_equation_t  *_field_interpolation_scalar_c2f_eq = NULL;

/*============================================================================
 * Private variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate equation(s) used for the field interpolation
 *         Two choices are available which can be combined:
 *         CS_CDO_FIELD_INTERPOLATION_SCALAR_C2V
 *         CS_CDO_FIELD_INTERPOLATION_SCALAR_C2F
 *
 * \param[in]      mode            kind of interpolation to perform
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_activate(cs_flag_t     mode)
{
  /* Store which kind of interpolation will be called */

  _field_interpolation_flag = mode;

  cs_property_t  *pty = cs_property_by_name("unity");
  if (pty == NULL) {
    pty = cs_property_add("unity", CS_PROPERTY_ISO);
    cs_property_def_iso_by_value(pty, "cells", 1.0);
  }

  if (mode & CS_CDO_FIELD_INTERPOLATION_SCALAR_C2V) {

    /* Add a new equation to build a cell --> vertices interpolation */

    _field_interpolation_scalar_c2v_eq
      = cs_equation_add("scalar_c2v_field_interpolation",
                        "scalar_c2v_field_interpolation",
                        CS_EQUATION_TYPE_PREDEFINED,
                        1,
                        CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_t  *eqp
      = cs_equation_get_param(_field_interpolation_scalar_c2v_eq);

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vcb");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RTOL, "1e-4");

    /* Add a diffusion term (Poisson eq.) */

    cs_equation_add_diffusion(eqp, pty);

  }

  if (mode & CS_CDO_FIELD_INTERPOLATION_SCALAR_C2F) {

    /* Add a new equation to build a cell --> faces interpolation */

    _field_interpolation_scalar_c2f_eq
      = cs_equation_add("scalar_c2f_field_interpolation",
                        "scalar_c2f_field_interpolation",
                        CS_EQUATION_TYPE_PREDEFINED,
                        1,
                        CS_PARAM_BC_HMG_NEUMANN);

    cs_equation_param_t  *eqp
      = cs_equation_get_param(_field_interpolation_scalar_c2f_eq);

    cs_equation_param_set(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_fb");
    cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
    cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "k_cycle");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_param_set(eqp, CS_EQKEY_ITSOL_RTOL, "1e-4");

    /* Add a diffusion term (Poisson eq.) */

    cs_equation_add_diffusion(eqp, pty);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate an array defined at vertices from an array defined at
 *         cells
 *
 * \param[in]      mesh            pointer to a mesh structure
 * \param[in]      cell_values     values at cells
 * \param[in, out] vtx_values      interpolated values at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_cell_to_vertices(const cs_mesh_t    *mesh,
                                            const cs_real_t    *cell_values,
                                            cs_real_t          *vtx_values)
{
  if (vtx_values == NULL)
    return; /* Should be allocated */

  if (_field_interpolation_scalar_c2v_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Equation related to the interpolation of cell array to"
              " vertices is not allocated.", __func__);

  cs_equation_t  *eq = _field_interpolation_scalar_c2v_eq;

  assert(CS_SPACE_SCHEME_CDOVCB == cs_equation_get_space_scheme(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
   * The linear solver is called inside and the field value is updated inside
   */

  cs_cdovcb_scaleq_interpolate(mesh,
                               cell_values,
                               eq->field_id,
                               eq->param,
                               eq->builder,
                               eq->scheme_context);

  /* Copy the computed solution into the given array at vertices */

  cs_field_t  *f = cs_field_by_id(eq->field_id);
  cs_array_real_copy(mesh->n_vertices, f->val, vtx_values);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Interpolate an array defined at faces from an array defined at
 *         cells
 *
 * \param[in]      mesh            pointer to a mesh structure
 * \param[in]      cell_values     values at cells
 * \param[in, out] face_values     interpolated values at faces
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_field_interpolation_cell_to_faces(const cs_mesh_t    *mesh,
                                         const cs_real_t    *cell_values,
                                         cs_real_t          *face_values)
{
  if (face_values == NULL)
    return; /* Should be allocated */

  if (_field_interpolation_scalar_c2f_eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Equation related to the interpolation of cell array to"
              " faces is not allocated.", __func__);

  cs_equation_t  *eq = _field_interpolation_scalar_c2f_eq;

  assert(CS_SPACE_SCHEME_CDOFB == cs_equation_get_space_scheme(eq));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Allocate, build and solve the algebraic system:
   * The linear solver is called inside and the field value is updated inside
   */

  cs_cdofb_scaleq_interpolate(mesh,
                              cell_values,
                              eq->field_id,
                              eq->param,
                              eq->builder,
                              eq->scheme_context);

  /* Copy the computed solution into the given array on faces */

  cs_real_t *f_values = cs_cdofb_scaleq_get_face_values(eq->scheme_context,
                                                        false);

  cs_array_real_copy(mesh->n_i_faces + mesh->n_b_faces, f_values, face_values);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
