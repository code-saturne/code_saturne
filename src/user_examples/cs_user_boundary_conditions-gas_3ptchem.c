/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

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
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for boundary condition definitions.
 */
/*----------------------------------------------------------------------------*/

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set boundary conditions to be applied.
 *
 * This function is called just before \ref cs_user_finalize_setup, and
 * boundary conditions can be defined in either of those functions,
 * depending on whichever is considered more readable or practical for a
 * given use.
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions_setup(cs_domain_t  *domain)
{
  CS_UNUSED(domain);

  /*! [setup_inlet] */
  const cs_zone_t *z = cs_boundary_zone_by_name("inlet");

  cs_real_t q_burner = cs_notebook_parameter_value_by_name("q_burner");

  cs_real_t q_inlet = q_burner / 360.;                             /* W */
  cs_real_t m_inlet = q_inlet / cs_glob_combustion_model->pcigas;  /* kg/s */
  cs_real_t r_inlet = 0.15;                                        /* m */

  cs_boundary_conditions_open_set_mass_flow_rate_by_value(z, m_inlet);

  cs_boundary_conditions_inlet_set_turbulence_hyd_diam(z, 2.*r_inlet);
  /*! [setup_inlet] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", access
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]\n\n
 *
 * Only the icodcl code values from the first component are used in the case
 * of vector or tensor fields, so the icodcl values can be defined as for
 * a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  /*! [loc_var_dec] */

  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_real_3_t *cdgfbo
    = (const cs_real_3_t *)(domain->mesh_quantities->b_face_cog);
  const cs_real_t *gxyz = cs_glob_physical_constants->gravity;

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  cs_combustion_model_t  *cm = cs_glob_combustion_model;

  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const int iflmab =  cs_field_get_key_int(CS_F_(vel), kbmasf);
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  /* Oxydant temperature needed by the combustion model
     when there is no oxydant inlet */

  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int icg = 0; icg < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; icg++)
    coefg[icg] = 0;
  coefg[1] = 1;

  cm->tinoxy = fp->t0;
  cm->hinoxy = cs_gas_combustion_t_to_h(coefg, cm->tinoxy);

  /*! [loc_var_dec] */

  /* Fuel inlet */

  /*! [inlet] */
  {
    const cs_zone_t *z = cs_boundary_zone_by_name("inlet");

    /* Fuel flow inlet type */
    cs_glob_bc_pm_info->ientfu[z->id] = 1;

    /* Inlet Temperature in K */
    cm->gas->tinfue = fp->t0;
  }
  /*! [inlet] */

  /* Open boundary: free inlet/outlet with backflow conditions.

     outlet: zero flux for velocity and temperature, prescribed pressure. */

  /*! [open] */
  {
    const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;

    cs_real_t *brom = CS_F_(rho_b)->val;

    int *p_icodcl = CS_F_(p)->bc_coeffs->icodcl;
    cs_real_t *p_rcodcl1 = CS_F_(p)->bc_coeffs->rcodcl1;
    cs_real_t *fm_rcodcl1 = CS_F_(fm)->bc_coeffs->rcodcl1;
    cs_real_t *fp2m_rcodcl1 = CS_F_(fp2m)->bc_coeffs->rcodcl1;
    cs_real_t *h_rcodcl1 = NULL;
    if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == 1) {
      h_rcodcl1 = CS_F_(h)->bc_coeffs->rcodcl1;
    }

    const cs_real_t xdh = 3.0;     /* acos(-1.d0) * 1 ( domain circumference) */
    const cs_real_t xitur = 0.02;  /* turbulence intensity */

    /* Loop on zone faces */

    const cs_zone_t *z = cs_boundary_zone_by_name("open");

    for (cs_lnum_t idx = 0; idx < z->n_elts; idx++) {
      const cs_lnum_t face_id = z->elt_ids[idx];

      bc_type[face_id] = CS_FREE_INLET;

      /* Imposed pressure */

      const cs_real_t pimp
        = fp->p0 + fp->ro0*cs_math_3_distance_dot_product(fp->xyzp0,
                                                          cdgfbo[face_id],
                                                          gxyz);
      p_icodcl[face_id] = 1;
      p_rcodcl1[face_id] = pimp;

      /* Outlet: Neumann conditions by defauls */

      if (bmasfl[face_id] > 0)
        continue;

      /* Inlet (backflow): Dirichlet conditions */

      const cs_lnum_t c_id = b_face_cells[face_id];

      cs_real_t uref2 = cs_math_3_square_norm(cvar_vel[c_id]);
      uref2 = cs_math_fmax(uref2, cs_math_epzero);

      cs_turbulence_bc_inlet_turb_intensity(face_id, uref2, xitur, xdh);

      /* Mixture fraction */
      fm_rcodcl1[face_id] = 0;

      /* Variance */
      fp2m_rcodcl1[face_id] = 0;

      /* Enthalpy */
      if (h_rcodcl1 != NULL)
        h_rcodcl1[face_id] = cm->hinoxy;

      /* Soot */
      if (cm->isoot == 1) {
        CS_F_(fsm)->bc_coeffs->rcodcl1[face_id] = 0;
        CS_F_(npm)->bc_coeffs->rcodcl1[face_id] = 0;
      }

      /* Density */
      brom[face_id] = fp->p0/(  cs_physical_constants_r
                              * cm->tinoxy/cm->gas->wmolg[1]);
    }
  }
  /*! [open] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
