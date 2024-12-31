/*============================================================================
 * Gas combustion model boundary conditions.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_boundary.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_combustion_gas.h"
#include "cs_combustion_model.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"

/* Prototypes for Fortran functions */

extern "C" int
cs_f_flamelet_rho_idx(void);

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_combustion_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_boundary_conditions.cpp

  \brief Gas combustion model boundary conditions.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type and macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destructor for coal boundary conditions inlet structure.
 *
 * \param[ci_p]  pointer to associated structure
 *
 * \return: null pointer
 */
/*----------------------------------------------------------------------------*/

static void *
_destroy_inlet(void  *ci_p)
{
  if (ci_p != NULL) {

    cs_combustion_bc_inlet_t *ci = (cs_combustion_bc_inlet_t *)ci_p;

    CS_FREE(ci);

  }

  return NULL;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to coal boundary conditions inlet structure.
 *
 * If no such structure was previously present, it is created and linked
 * to the matching open boundary condition inlet.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_combustion_bc_inlet_t *
cs_combustion_boundary_conditions_get_inlet(const  cs_zone_t   *zone)
{
  auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
              (cs_boundary_conditions_get_model_inlet(zone));

  /* Add and initialize coal inlet if not present */

  if (ci == NULL) {
    CS_MALLOC(ci, 1, cs_combustion_bc_inlet_t);

    ci->zone = zone;

    ci->ientox = 0;
    ci->ientfu = 0;
    ci->ientgf = 0;
    ci->ientgb = 0;

    ci->fment = nan("");  /* Set to usable value if used */
    ci->tkent = nan("");
    ci->tgf = nan("");

    /* Now register in high-level boundary conditions mechanism */

    cs_boundary_conditions_assign_model_inlet(zone,
                                              (void *)ci,
                                              (void *)_destroy_inlet);

    /* Link with boundary mechanism;
       This will allow simplifying later loops */

    assert(cs_glob_boundaries != NULL);

    cs_boundary_t *bdy = cs_glob_boundaries;
    if (cs_boundary_id_by_zone_id(bdy, zone->id) < 0)
      cs_boundary_add(bdy, CS_BOUNDARY_INLET, zone->name);

  }

  return (ci);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions(int  bc_type[])
{
  /* Initialization
   * ============== */

  const int *pm_flag = cs_glob_physical_model_flag;
  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* Boundary conditions for H */
  cs_real_t *rcodcl1_h = nullptr;

  if (   pm_flag[CS_COMBUSTION_3PT] > -1
      || pm_flag[CS_COMBUSTION_SLFM] == 1
      || pm_flag[CS_COMBUSTION_SLFM] == 3) {
    if (CS_F_(h) != nullptr)
      rcodcl1_h = CS_F_(h)->bc_coeffs->rcodcl1;
  }

  /* Boundary conditions mixture fraction and its variance */
  cs_real_t *rcodcl1_ifm = nullptr, *rcodcl1_ifp2m_ifsqm = nullptr;
  cs_real_t *rcodcl1_ipvm = nullptr;

  /* Soot model */
  cs_real_t *rcodcl1_ifsm = nullptr, *rcodcl1_inpm = nullptr;

  cs_field_t *f;

  f = cs_field_by_name_try("mixture_fraction");
  if (f != nullptr)
    rcodcl1_ifm = f->bc_coeffs->rcodcl1;

  /* Remark of SLFM model: a priori, 2nd moment of mixture fraction and
     progress variable are unknown until mixture fraction is solved.
     Particular treatment needed in cs_solve_transported_variables, not here */

  int mode_fp2m = 0;
  if (pm_flag[CS_COMBUSTION_SLFM] >= 0) {
    if (cm->mode_fp2m == 1)
      mode_fp2m = 1;
  }

  if (mode_fp2m == 0) {
    f = cs_field_by_name_try("mixture_fraction_variance");
    if (f != nullptr)
      rcodcl1_ifp2m_ifsqm = f->bc_coeffs->rcodcl1;
  }
  else if (mode_fp2m == 1) {
    f = cs_field_by_name_try("mixture_fraction_2nd_moment");
    if (f != nullptr)
        rcodcl1_ifp2m_ifsqm = f->bc_coeffs->rcodcl1;
  }

  if (pm_flag[CS_COMBUSTION_SLFM] >= 2) {
    f = cs_field_by_name_try("progress_variable");
    if (f != nullptr)
      rcodcl1_ipvm = f->bc_coeffs->rcodcl1;
  }

  if (cm->isoot == 1) {
    rcodcl1_ifsm = cs_field_by_name("soot_mass_fraction")->bc_coeffs->rcodcl1;
    rcodcl1_inpm
      = cs_field_by_name("soot_precursor_number")->bc_coeffs->rcodcl1;
  }

  /* Loop on inlet boundaries
     ======================== */

  assert(cs_glob_boundaries != NULL);

  const cs_boundary_t *bdy = cs_glob_boundaries;

  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

    const cs_lnum_t n_elts = z->n_elts;
    const cs_lnum_t *elt_ids = z->elt_ids;

    auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    cs_real_t qimp = cs_boundary_conditions_open_get_mass_flow_rate(z);

    cs_real_t h_in = -HUGE_VALF;
    cs_real_t ifm_in = 0, ifp2m_ifsqm_in = 0, ipvm_in = 0;

    /* Fuel inlet at tinfue */
    if (ci->ientfu == 1) {
      h_in = cm->hinfue;
      ifm_in = 1.;
      if (mode_fp2m == 1)
        ifp2m_ifsqm_in = 1.;
      ipvm_in = 1.;
    }

    /* Oxydant inlet at tinoxy */
    if (ci->ientox == 1) {
      h_in = cm->hinoxy;
      ifm_in = 0.;
    }

    if (ci->ientfu == 1 || ci->ientox == 1) {

      if (qimp < 0)
        continue;

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
        cs_lnum_t elt_id = elt_ids[elt_idx];
        if (   bc_type[elt_id] == CS_INLET
            || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

          rcodcl1_ifm[elt_id] = ifm_in;      // mean mixture fraction

          // mean mixture fraction variance or 2nd moment of mixture fraction
          if (rcodcl1_ifp2m_ifsqm != nullptr)
            rcodcl1_ifp2m_ifsqm[elt_id] = ifp2m_ifsqm_in;

          if (rcodcl1_ipvm != nullptr)
            rcodcl1_ipvm[elt_id] = ipvm_in;  // progress variable

          if (rcodcl1_h != nullptr)
            rcodcl1_h[elt_id] = h_in;        // mixture enthalpy

          if (rcodcl1_ifsm != nullptr) {     // soot model
            rcodcl1_ifsm[elt_id] = 0;
            rcodcl1_inpm[elt_id] = 0;
          }

        }

      } /* loop on zone faces */

    } /* Test on zone type */

  } /* loop on zones */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density at boundary for pulverized coal combustion.
 *
 * This is based on boundary condition definitions, but is called at an
 * earlier stage in the time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_density(void)
{
  /* Initialization
   * ============== */

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_nreal_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;
  const cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const cs_real_t pther = cs_glob_fluid_properties->pther;

  const cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  /* Mass density on edges for all faces */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    brom[face_id] = crom[cell_id];
  }

  /* Recompute density at inlets */

  assert(cs_glob_boundaries != NULL);

  const cs_boundary_t *bdy = cs_glob_boundaries;

  /* loop on boundary zones, ignore non-inlet zones */

  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {

    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);

    const cs_lnum_t n_elts = z->n_elts;
    const cs_lnum_t *elt_ids = z->elt_ids;

    auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    if (ci->ientfu == 1 || ci->ientox == 1) {

      cs_real_t rho_b_in = 0;

      if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] > -1) {

        if (ci->ientfu == 1)
          rho_b_in = pther/(cs_physical_constants_r*cm->tinfue/cm->wmolg[0]);
        else if (ci->ientox ==1)
          rho_b_in = pther/(cs_physical_constants_r*cm->tinoxy/cm->wmolg[1]);

      }
      else if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] > -1) {

        /* Fortran index:
         * - if (ci->ientfu == 1)
         *   flamelet_library(flamelet_rho, 1, 1, 1, nzm)
         * - if (ci->ientox == 1)
         *   flamelet_library(flamelet_rho, 1, 1, 1, 1)
         *
         *  TODO: to complete C/C++ migration,
         *        replace call to cs_f_flamelet_rho_idx()
         *        with "cm->flamelet_rho_idx", or
         *        add cs_flamelet_library_idx function handling all indexes.
        */

        int flamelet_rho_idx = cs_f_flamelet_rho_idx();

        if (ci->ientfu == 1)
          flamelet_rho_idx += (cm->nzm - 1)*(cm->nzvar)*(cm->nki)*(cm->nxr);

        rho_b_in = cm->flamelet_library_p[flamelet_rho_idx];

      }

      for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
        cs_lnum_t face_id = elt_ids[elt_idx];
        cs_lnum_t cell_id = b_face_cells[face_id];
        const cs_real_t vs = cs_math_3_dot_product(cvar_vel[cell_id],
                                                   f_n[face_id]);

        if (vs < 0) // inflow
          brom[face_id] = rho_b_in;
      }

    } /* Test on zone type */

  } /* loop on boundary zones */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
