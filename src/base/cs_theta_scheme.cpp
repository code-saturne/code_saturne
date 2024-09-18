/*============================================================================
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *===========================================================================*/

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
 *---------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_restart_default.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_theta_scheme.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 * parameters:
 *   iappel  <--  call id (before or after cs_physical_properties_update)
 *
 * returns:
 *---------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var(const cs_lnum_t  iappel)
{
  /* Initialization */

  const int key_t_ext_id = cs_field_key_id("time_extrapolated");
  const int kthetvs = cs_field_key_id("diffusivity_extrapolated");

  cs_real_t *cpro_rho_mass = nullptr;
  cs_real_t *bpro_rho_mass = nullptr;

  if (cs_field_by_name_try("density_mass") != nullptr)
    cpro_rho_mass = cs_field_by_name_try("density_mass")->val;
  if (cs_field_by_name_try("boundary_density_mass") != nullptr)
    bpro_rho_mass = cs_field_by_name_try("boundary_density_mass")->val;

  const cs_velocity_pressure_param_t  *vp_param
    = cs_glob_velocity_pressure_param;

  const int kimasf = cs_field_key_id("inner_mass_flux_id");
  const int kbmasf = cs_field_key_id("boundary_mass_flux_id");
  int iflmas = cs_field_get_key_int(CS_F_(vel), kimasf);
  int iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);

  cs_real_t *imasfl = cs_field_by_id(iflmas)->val;
  cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;
  cs_real_t *imasfl_pre = cs_field_by_id(iflmas)->val_pre;
  cs_real_t *bmasfl_pre = cs_field_by_id(iflmab)->val_pre;

  const cs_real_t *crom   = CS_F_(rho)->val;
  const cs_real_t *brom   = CS_F_(rho_b)->val;

  cs_real_t *cpro_viscl  = CS_F_(mu)->val;
  cs_real_t *cproa_viscl = CS_F_(mu)->val_pre;
  cs_real_t *cpro_visct  = CS_F_(mu_t)->val;
  cs_real_t *cproa_visct = CS_F_(mu_t)->val_pre;

  cs_field_t *f_cp = CS_F_(cp);
  cs_real_t *cpro_cp = nullptr;
  cs_real_t *cproa_cp = nullptr;
  if (f_cp != nullptr) {
    cpro_cp  = f_cp->val;
    cproa_cp = f_cp->val_pre;
  }

  const int time_order = cs_glob_time_scheme->time_order;

  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces   = cs_glob_mesh->n_i_faces;
  const cs_lnum_t n_b_faces   = cs_glob_mesh->n_b_faces;

  const int k_sca = cs_field_key_id("scalar_id");

  /* At the very start of the time step */
  if (iappel == 1) {

    /*----------------------------------------------------------------------
     * Store the previous mass flux (n-1->n) in *_mass_flux_prev
     * Note: if there is no previous values, nothing is done for explicit
     * schemes (istmpf=0) a specific treatment is done because previous value
     * is used as a work array.
     *----------------------------------------------------------------------*/

    if ((cs_glob_time_scheme->istmpf > 0) || (vp_param->staggered == 1)) {
      if (vp_param->itpcol == 0) {
        cs_field_current_to_previous(cs_field_by_id(iflmas));
        cs_field_current_to_previous(cs_field_by_id(iflmab));
      }
      else {
        for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
          cs_real_t flux      = imasfl[face_id];
          imasfl[face_id]     = 2.*imasfl[face_id] - imasfl_pre[face_id];
          imasfl_pre[face_id] = flux;
        }
        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
          cs_real_t flux      = bmasfl[face_id];
          bmasfl[face_id]     = 2.*bmasfl[face_id] - bmasfl_pre[face_id];
          bmasfl_pre[face_id] = flux;
        }
      }
    }

    /*----------------------------------------------------------------------
     * If required, the density at time step n-1 is updated
     * Note that for VOF and dilatable algorithms, density at time step n-2
     * is also updated
     * Note that at the begining of the calculation, previous values have
     * been initialized by inivar
     *----------------------------------------------------------------------*/

    if (cs_glob_fluid_properties->irovar > 0) {
      if (time_order == 2 && vp_param->itpcol == 1) {
        const cs_real_t *croma  = CS_F_(rho)->val_pre;
        const cs_real_t *cromaa = CS_F_(rho)->vals[2];
        const cs_real_t *broma  = CS_F_(rho_b)->val_pre;
        const cs_real_t *bromaa = CS_F_(rho_b)->vals[2];

        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++) {
          cpro_rho_mass[c_id] = 3.*crom[c_id] - 3.*croma[c_id]
                                + cromaa[c_id];
        }

        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
          bpro_rho_mass[face_id] = 3.*brom[face_id] - 3.*broma[face_id]
                                   + bromaa[face_id];
        }
      }

      cs_field_current_to_previous(CS_F_(rho));
      cs_field_current_to_previous(CS_F_(rho_b));

      if (((cs_glob_velocity_pressure_model->idilat > 1
            || cs_glob_vof_parameters->vof_model > 0) && vp_param->itpcol == 0)
          || cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 3) {
        for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
          cpro_rho_mass[c_id] = crom[c_id];
        for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
          bpro_rho_mass[face_id] = brom[face_id];
      }
    }

    /* Extrapolate viscosities and heat capacity if required */
    if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0)
      cs_field_current_to_previous(CS_F_(mu));
    if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0)
      cs_field_current_to_previous(CS_F_(mu_t));
    if (f_cp != nullptr) {
      if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0)
        cs_field_current_to_previous(f_cp);
    }

    /* Extrapolate scalars if required */
    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
        continue;
      int i_sca = cs_field_get_key_int(f, k_sca);
      if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0))
        cs_field_current_to_previous(f);
    }

  }

  /* Just after cs_physical_properties_update (and thus before
     cs_solve_navier_stokes) */
  else if (iappel == 2) {

    /*----------------------------------------------------------------------
     * Update previous values : only done if values are unavailable at the
     * first time step or when reading the restart file
     *----------------------------------------------------------------------*/

    if (cs_restart_present() == 1) {
      /* Update densities, viscosities and heat capacity if required */
      if (cs_glob_fluid_properties->irovar > 0) {
        if (cs_restart_get_field_read_status(CS_F_(rho)->id) == 0)
          cs_field_current_to_previous(CS_F_(rho));
        if (cs_restart_get_field_read_status(CS_F_(rho_b)->id) == 0)
          cs_field_current_to_previous(CS_F_(rho_b));
      }
      if (cs_restart_get_field_read_status(CS_F_(mu)->id) == 0)
        cs_field_current_to_previous(CS_F_(mu));
      if (cs_restart_get_field_read_status(CS_F_(mu_t)->id) == 0)
        cs_field_current_to_previous(CS_F_(mu_t));
      if (f_cp != nullptr) {
        if (cs_restart_get_field_read_status(f_cp->id) == 0)
          cs_field_current_to_previous(f_cp);
      }

      /* Update scalars diffusivity if required */
      for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
        cs_field_t *f = cs_field_by_id(f_id);
        if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
          continue;
        int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
        if ((k_f_id > 0) && (cs_restart_get_field_read_status(k_f_id) == 0))
          cs_field_current_to_previous(cs_field_by_id(k_f_id));
      }
    }

    /*----------------------------------------------------------------------
     * Extrapolation of viscosities, heat capacity and scalars diffusivity
     * for the theta scheme. Values in the previous time step are kept before
     * the end the time iteration
     *----------------------------------------------------------------------*/

    /* Extrapolate viscosities and heat capacity */
    if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0) {
      const cs_real_t theta = cs_glob_time_scheme->thetvi;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t _viscol = cpro_viscl[c_id];
        cpro_viscl[c_id] = (1. + theta) * cpro_viscl[c_id]
                           - theta * cproa_viscl[c_id];
        cproa_viscl[c_id] = _viscol;
      }
    }
    if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0) {
      const cs_real_t theta = cs_glob_time_scheme->thetvi;
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t _viscot = cpro_visct[c_id];
        cpro_visct[c_id] = (1. + theta) * cpro_visct[c_id]
                           - theta * cproa_visct[c_id];
        cproa_visct[c_id] = _viscot;
      }
    }
    if (f_cp != nullptr) {
      if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0) {
        const cs_real_t theta = cs_glob_time_scheme->thetcp;
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t _cp = cpro_cp[c_id];
          cpro_cp[c_id] = (1. + theta) * cpro_cp[c_id]
                              - theta * cproa_cp[c_id];
          cproa_cp[c_id] = _cp;
        }
      }
    }

    /* Extrapolate scalars diffusivity */
    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
        continue;
      int i_sca = cs_field_get_key_int(f, k_sca);

      /* Is a scalar and should be computed with theta scheme */
      if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0)) {

        const cs_real_t theta = cs_field_get_key_double(f, kthetvs);
        int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
        cs_real_t *cpro_visca = cs_field_by_id(k_f_id)->val;
        cs_real_t *cproa_visca = cs_field_by_id(k_f_id)->val_pre;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          cs_real_t _visca = cpro_visca[c_id];
          cpro_visca[c_id] = (1. + theta) * cpro_visca[c_id]
                             - theta * cproa_visca[c_id];
          cproa_visca[c_id] = _visca;
        }

      }
    }

    return;

  }

  /* Just after cs_solve_navier_stokes, in loops for U/P and ALE */
  else if (iappel == 3) {

    /*----------------------------------------------------------------------
     * Only update an unique mass flux.
     * For the second order scheme (istmpf=2), theta is set to 0.5.
     * For the explicit scheme (istmpf=0), the mass flux is set to the
     * previous one F(n) in the last iteration. Another modification will
     * be done in iappel = 4.
     * Nothing is done for standard scheme (istmpf=1).
     * In case of iterations in cs_solve_navier_stokes, the following process
     * is applied in all iterations except the last one if istmpf=0
     *----------------------------------------------------------------------*/

    if ((cs_glob_time_scheme->istmpf == 2) && (vp_param->itpcol == 1)) {
      cs_real_t theta = 0.5;
      cs_real_t aa = 1./(2.-theta);
      cs_real_t bb = (1.-theta)/(2.-theta);

      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
        imasfl[face_id] = aa * imasfl[face_id] + bb * imasfl[face_id];
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        bmasfl[face_id] = aa * bmasfl[face_id] + bb * bmasfl[face_id];
    }
    else if (cs_glob_time_scheme->istmpf == 0) {
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
        imasfl[face_id] = imasfl_pre[face_id];
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        bmasfl[face_id] = bmasfl_pre[face_id];
    }

    return;

  }

  /* After cs_solve_navier_stokes and loops for U/P and ALE */
  else if (iappel == 4) {

    /*----------------------------------------------------------------------
     * Only update an unique mass flux.
     * Nothing is done for standard and second order schemes (istmpf=1 and 2)
     * For the explicit scheme (istmpf=0), F(n+1) is saved in imasfl_pre but
     * computations are done with F(n) in imasfl
     * In case of iterations in cs_solve_navier_stokes, the following process
     * is applied in all iterations if istmpf differs from 0 and only in the
     * last iteration if istmpf=0. By doing so, from the second sub-iteration
     * the computation will be done with F(n+1) and no longer with F(n).
     * It is assumed here that the user chose to do sub-iterations also to
     * have an implicit mass flux
     *----------------------------------------------------------------------*/

    if (cs_glob_time_scheme->istmpf == 0) {
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {
        cs_real_t flux      = imasfl[face_id];
        imasfl[face_id]     = imasfl_pre[face_id];
        imasfl_pre[face_id] = flux;
      }
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
        cs_real_t flux      = bmasfl[face_id];
        bmasfl[face_id]     = bmasfl_pre[face_id];
        bmasfl_pre[face_id] = flux;
      }
    }

    return;

  }

  /* Just after scalai */
  else if (iappel == 5) {

    /*----------------------------------------------------------------------
     * Previous modifications are reverted in order to be ready for the next
     * iteration in time.
     * Nothing is done for standard and second order schemes (istmpf=1 and 2)
     * For the explicit scheme (istmpf=0), imasfl is set to F(n+1)
     *----------------------------------------------------------------------*/

    /* Restore the mass flux */
    if (cs_glob_time_scheme->istmpf == 0) {
      for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++)
        imasfl[face_id] = imasfl_pre[face_id];
      for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
        bmasfl[face_id] = bmasfl_pre[face_id];
    }

    /* Restore physical properties */
    if (cs_field_get_key_int(CS_F_(mu), key_t_ext_id) > 0)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_viscl[c_id] = cproa_viscl[c_id];
    if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cpro_visct[c_id] = cproa_visct[c_id];
    if (f_cp != nullptr) {
      if (cs_field_get_key_int(f_cp, key_t_ext_id) > 0) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_cp[c_id] = cproa_cp[c_id];
      }
    }

    /* Restore scalars diffusivity */
    for (int f_id = 0; f_id < cs_field_n_fields(); f_id ++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE) || (f->type & CS_FIELD_CDO))
        continue;
      int i_sca = cs_field_get_key_int(f, k_sca);

      /* Is a scalar and should be computed with theta scheme */
      if ((i_sca > 0) && (cs_field_get_key_int(f, key_t_ext_id) > 0)) {

        int k_f_id = cs_field_get_key_int(f, cs_field_key_id("diffusivity_id"));
        cs_real_t *cpro_visca = cs_field_by_id(k_f_id)->val;
        cs_real_t *cproa_visca = cs_field_by_id(k_f_id)->val_pre;

        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
          cpro_visca[c_id] = cproa_visca[c_id];

      }
    }

    return;

  }
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
