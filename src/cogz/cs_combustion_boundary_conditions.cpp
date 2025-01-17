/*============================================================================
 * Gas combustion model boundary conditions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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

#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "base/cs_boundary.h"
#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_parameters.h"
#include "base/cs_parameters_check.h"
#include "base/cs_physical_constants.h"
#include "base/cs_prototypes.h"
#include "base/cs_time_step.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"

#include "cogz/cs_combustion_gas.h"
#include "cogz/cs_combustion_ht_convert.h"
#include "pprt/cs_combustion_model.h"
#include "pprt/cs_physical_model.h"

/* Prototypes for Fortran functions */

extern "C" int
cs_f_flamelet_rho_idx(void);

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cogz/cs_combustion_boundary_conditions.h"

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
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* Boundary conditions for H */
  cs_real_t *rcodcl1_h = nullptr;

  if (   pm_flag[CS_COMBUSTION_3PT] > -1
      || pm_flag[CS_COMBUSTION_SLFM] == 1
      || pm_flag[CS_COMBUSTION_SLFM] == 3) {
    if (CS_F_(h) != nullptr)
      rcodcl1_h = CS_F_(h)->bc_coeffs->rcodcl1;
  }

  /* Boundary conditions mixture fraction and its variance */
  cs_real_t *rcodcl1_fm = nullptr, *rcodcl1_fp2m_fsqm = nullptr;
  cs_real_t *rcodcl1_pvm = nullptr;

  /* Soot model */
  cs_real_t *rcodcl1_ifsm = nullptr, *rcodcl1_inpm = nullptr;

  cs_field_t *f;

  f = cs_field_by_name_try("mixture_fraction");
  if (f != nullptr)
    rcodcl1_fm = f->bc_coeffs->rcodcl1;

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
      rcodcl1_fp2m_fsqm = f->bc_coeffs->rcodcl1;
  }
  else if (mode_fp2m == 1) {
    f = cs_field_by_name_try("mixture_fraction_2nd_moment");
    if (f != nullptr)
        rcodcl1_fp2m_fsqm = f->bc_coeffs->rcodcl1;
  }

  if (pm_flag[CS_COMBUSTION_SLFM] >= 2) {
    f = cs_field_by_name_try("progress_variable");
    if (f != nullptr)
      rcodcl1_pvm = f->bc_coeffs->rcodcl1;
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
    cs_real_t fm_in = 0, fp2m_fsqm_in = 0, pvm_in = 0;

    if ((ci->ientfu != 1 && ci->ientox != 1) || qimp < 0)
      continue;

    cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
    for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
      coefg[0] = 0;

    /* Fuel inlet at tinfue */
    if (ci->ientfu == 1) {
      if (   cm->type == CS_COMBUSTION_3PT_ADIABATIC
          || cm->type == CS_COMBUSTION_3PT_PERMEATIC) {
        coefg[0] = 1.; coefg[1] = 0.; coefg[2] = 0.;
        cm->hinfue = cs_gas_combustion_t_to_h(coefg, cm->tinfue);
      }
      h_in = cm->hinfue;
      fm_in = 1.;
      if (mode_fp2m == 1)
        fp2m_fsqm_in = 1.;
      pvm_in = 1.;
    }

    /* Oxydant inlet at tinoxy */
    if (ci->ientox == 1) {
      if (   cm->type == CS_COMBUSTION_3PT_ADIABATIC
          || cm->type == CS_COMBUSTION_3PT_PERMEATIC) {
        coefg[0] = 0.; coefg[1] = 1.; coefg[2] = 0.;
        cm->hinoxy = cs_gas_combustion_t_to_h(coefg, cm->tinoxy);
      }
      h_in = cm->hinoxy;
      fm_in = 0.;
    }

    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
      cs_lnum_t elt_id = elt_ids[elt_idx];
      if (   bc_type[elt_id] == CS_INLET
          || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

        rcodcl1_fm[elt_id] = fm_in;      // mean mixture fraction

        // mean mixture fraction variance or 2nd moment of mixture fraction
        if (rcodcl1_fp2m_fsqm != nullptr)
          rcodcl1_fp2m_fsqm[elt_id] = fp2m_fsqm_in;

        if (rcodcl1_pvm != nullptr)
          rcodcl1_pvm[elt_id] = pvm_in;  // progress variable

        if (rcodcl1_h != nullptr)
          rcodcl1_h[elt_id] = h_in;        // mixture enthalpy

        if (rcodcl1_ifsm != nullptr) {     // soot model
          rcodcl1_ifsm[elt_id] = 0;
          rcodcl1_inpm[elt_id] = 0;
        }

      }

    } /* loop on zone faces */

  } /* loop on zones */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion with EBU model.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_ebu(int  bc_type[])
{
  /* Initialization
   * ============== */

  const int ebu_model = cs_glob_physical_model_flag[CS_COMBUSTION_EBU];
  assert(ebu_model > -1);

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* Boundary conditions mass fraction of fresh gas */
  cs_real_t *rcodcl1_ygfm = nullptr;
  rcodcl1_ygfm = cs_field_by_name("fresh_gas_fraction")->bc_coeffs->rcodcl1;

  /* Boundary conditions for H */
  cs_real_t *rcodcl1_h = nullptr;

  if (ebu_model == 1 || ebu_model == 3) {
    if (CS_F_(h) != nullptr)
      rcodcl1_h = CS_F_(h)->bc_coeffs->rcodcl1;
  }

  /* Boundary conditions for mixture fraction  */
  cs_real_t *rcodcl1_fm = nullptr;

  if (ebu_model == 2 || ebu_model == 3)
    rcodcl1_fm = cs_field_by_name("mixture_fraction")->bc_coeffs->rcodcl1;

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

    cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
    for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
      coefg[0] = 0;

    cs_real_t ygfm_in = 0, h_in = 0, fm_in = ci->fment;

    /* Cold premix or dilution */
    if (ci->ientgf == 1) {
      coefg[0] = ci->fment;
      coefg[1] = 1. - ci->fment;
      cs_real_t tgasf = ci->tkent;
      cs_real_t hgasf = cs_gas_combustion_t_to_h(coefg, tgasf);
      h_in = hgasf;
      ygfm_in = 1.;
    }

    /* Burned gas inlet (pilot flame) */
    else if (ci->ientgb == 1) {
      cs_real_t fs_1 = cm->fs[0];
      coefg[0] = fmax(0, (ci->fment-fs_1)/(1.-fs_1));
      coefg[2] = (ci->fment - coefg[0])/fs_1;
      coefg[1] = 1.0 - coefg[0] - coefg[2];
      cs_real_t tgasb = ci->tkent;
      cs_real_t hgasb = cs_gas_combustion_t_to_h(coefg, tgasb);
      h_in = hgasb;
      ygfm_in = 0.;
    }

    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
      cs_lnum_t elt_id = elt_ids[elt_idx];
      if (   bc_type[elt_id] == CS_INLET
          || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

        rcodcl1_ygfm[elt_id] = ygfm_in;     // fresh gas fraction

        if (rcodcl1_fm != nullptr)
          rcodcl1_fm[elt_id] = fm_in;       // mixture fraction

        if (rcodcl1_h != nullptr)
          rcodcl1_h[elt_id] = h_in;         // mixture enthalpy

      }

    } /* loop on zone faces */

  } /* loop on zones */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic boundary condition for gas combustion with Libby-Williams
 *        model.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_lw(int  bc_type[])
{
  /* Initialization
   * ============== */

  const int lw_model = cs_glob_physical_model_flag[CS_COMBUSTION_LW];
  assert(lw_model > -1);

  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* Boundary conditions of mass fraction of fuel */
  cs_real_t *rcodcl1_yfm
    = cs_field_by_name("mass_fraction")->bc_coeffs->rcodcl1;

  /* Boundary conditions of mass fraction variance of fuel */
  cs_real_t *rcodcl1_yfp2m
    = cs_field_by_name("mass_fraction_variance")->bc_coeffs->rcodcl1;

  /* Boundary conditions of mixture fraction */
  cs_real_t *rcodcl1_fm
    = cs_field_by_name("mixture_fraction")->bc_coeffs->rcodcl1;

  /* Boundary conditions of mixture fraction variance */
  cs_real_t *rcodcl1_fp2m
    = cs_field_by_name("mixture_fraction_variance")->bc_coeffs->rcodcl1;

  /* Boundary conditions for H */
  cs_real_t *rcodcl1_coyfp = nullptr;
  if (lw_model >= 2) {
    rcodcl1_coyfp
      = cs_field_by_name("mass_fraction_covariance")->bc_coeffs->rcodcl1;
  }

  /* Boundary conditions for H */
  cs_real_t *rcodcl1_h = nullptr;

  if (lw_model == 1 || lw_model == 3 || lw_model == 5) {
    if (CS_F_(h) != nullptr)
      rcodcl1_h = CS_F_(h)->bc_coeffs->rcodcl1;
  }

  /* Min and max values at inlets */

  cs_real_t fm_min = HUGE_VALF, h_min = HUGE_VALF;
  cs_real_t fm_max = -HUGE_VALF, h_max = -HUGE_VALF;

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

    cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
    for (int i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; i++)
      coefg[0] = 0;

    cs_real_t yfm_in = 0, h_in = 0;

    /* Enthlapy of gas mix. Assumes inlet is either ientgf = 1 or ientgf = 1 */

    if (ci->ientgf == 1) {       /* Fresh gas inlet */
      coefg[0] = ci->fment;
      coefg[1] = 1. - ci->fment;
      cs_real_t tgasf = ci->tkent;
      cs_real_t hgasf = cs_gas_combustion_t_to_h(coefg, tgasf);
      h_in = hgasf;
      yfm_in = ci->fment;
    }
    else if (ci->ientgb == 1) {  /* Burned gas inlet */
      coefg[0] = ci->fment;
      coefg[1] = 1. - ci->fment;
      cs_real_t tgasb = ci->tkent;
      cs_real_t hgasb = cs_gas_combustion_t_to_h(coefg, tgasb);
      h_in = hgasb;
      yfm_in = 0;
    }

    cs_real_t fm_in = ci->fment;

    for (cs_lnum_t elt_idx = 0; elt_idx < n_elts; elt_idx++) {
      cs_lnum_t elt_id = elt_ids[elt_idx];
      if (   bc_type[elt_id] == CS_INLET
          || bc_type[elt_id] == CS_CONVECTIVE_INLET) {

        rcodcl1_yfm[elt_id]   = yfm_in;     // fuel mass fraction
        rcodcl1_yfp2m[elt_id] = 0.;         // mass fraction variance
        rcodcl1_fm[elt_id]    = fm_in;      // mixture fraction
        rcodcl1_fp2m[elt_id]  = 0.;         // mixture fraction variance

        if (rcodcl1_coyfp != nullptr)
          rcodcl1_coyfp[elt_id] = 0.;       // mass fraction covariance

        if (rcodcl1_h != nullptr)
          rcodcl1_h[elt_id] = h_in;         // mixture enthalpy

      }

    } /* loop on zone faces */

    /* Update min/max */

    if (fm_in < fm_min) {
      fm_min = fm_in;
      h_min = h_in;
    }
    if (fm_in > fm_max) {
      fm_max = fm_in;
      h_max = h_in;
    }

  } /* loop on zones */

  cs_parall_min_loc_vals(1, &fm_min, & h_min);
  cs_parall_max_loc_vals(1, &fm_max, & h_max);

  cm->lw.fmin = fm_min;
  cm->lw.fmax = fm_max;
  cm->lw.hmin = h_min;
  cm->lw.hmax = h_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute density at boundary for gas combustion.
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
/*!
 * \brief Compute density at boundary for gas combustion, using
 *        EBU or Libby-Williams models.
 *
 * This is based on boundary condition definitions, but is called at an
 * earlier stage in the time step.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_density_ebu_lw(void)
{
  /* Initialization
   * ============== */

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_nreal_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;
  const cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

  const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  const cs_real_t p0 = cs_glob_fluid_properties->p0;

  const cs_real_t *crom = CS_F_(rho)->val;
  cs_real_t *brom = CS_F_(rho_b)->val;

  /* Mass density on edges for all faces */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    brom[face_id] = crom[cell_id];
  }

  if (cs_glob_time_step->nt_cur <= 1)
    return;

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

    if (ci->ientgb == 1 || ci->ientgf == 1) {

      cs_real_t coefg[3] = {ci->fment, 1.-ci->fment, 0};
      if (ci->ientgb == 1) {
        cs_real_t fs_1 = cm->fs[0];
        coefg[0] = fmax(0, (ci->fment-fs_1) / (1.-fs_1));
        coefg[2] = (ci->fment-coefg[0]) / fs_1;
        coefg[1] = 1. - coefg[0] - coefg[2];
      }
      assert(cm->n_gas_species <= 3);

      cs_real_t nbmol = 0;
      for (int igg = 0; igg < cm->n_gas_species; igg++)
        nbmol += coefg[igg]/cm->wmolg[igg];
      cs_real_t masmg = 1./ nbmol;
      cs_real_t temsmm = ci->tkent / masmg;

      cs_real_t rho_b_in = p0 / (cs_physical_constants_r*temsmm);

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
/*!
 * \brief Compute mean inlet enthalpy at boundary for
 *        EBU and Libby-Williams models.
 *
 * \param[out]  fmm  mean inlet mixture fraction
 * \param[out]  tkm  mean inlet mixture temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_boundary_conditions_mean_inlet_ebu_lw(cs_real_t  *fmm,
                                                    cs_real_t  *tkm)
{
  assert(cs_glob_boundaries != NULL);

  cs_real_t zsum[3] = {0, 0, 0};

  /* loop on boundary zones, ignore non-inlet zones */

  const cs_boundary_t *bdy = cs_glob_boundaries;
  for (int bdy_idx = 0; bdy_idx < bdy->n_boundaries; bdy_idx += 1) {
    if (! (bdy->types[bdy_idx] & CS_BOUNDARY_INLET))
      continue;

    const cs_zone_t *z = cs_boundary_zone_by_id(bdy->zone_ids[bdy_idx]);
    auto ci = reinterpret_cast<cs_combustion_bc_inlet_t *>
                (cs_boundary_conditions_get_model_inlet(z));

    const cs_real_t qimp = cs_boundary_conditions_open_get_mass_flow_rate(z);

    zsum[0] += qimp * ci->fment;
    zsum[1] += qimp * ci->tkent;
    zsum[2] += qimp;
  }
  cs_parall_sum(3, CS_REAL_TYPE, zsum);

  if (zsum[2] > cs_math_epzero) {
    *fmm = zsum[0] / zsum[2];
    *tkm = zsum[1] / zsum[2];
  }
  else {
    *fmm = 0;
    *tkm = cs_glob_fluid_properties->t0;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
