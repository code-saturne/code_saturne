/*============================================================================
 * Solve the Navier-Stokes equations, source term convection
 * diffusion equations for scalars ... .
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_1d_wall_thermal.h"
#include "cs_ale.h"
#include "cs_array.h"
#include "cs_ast_coupling.h"
#include "cs_atmo.h"
#include "cs_array.h"
#include "cs_assert.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_coupling.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_compute_thermo_pressure_density.h"
#include "cs_ctwr_source_terms.h"
#include "cs_dilatable_scalar_diff_st.h"
#include "cs_fan.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_head_losses.h"
#include "cs_lagr.h"
#include "cs_lagr_head_losses.h"
#include "cs_log.h"
#include "cs_mobile_structures.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties_default.h"
#include "cs_prototypes.h"
#include "cs_porous_model.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_solve.h"
#include "cs_sat_coupling.h"
#include "cs_solve_navier_stokes.h"
#include "cs_solve_transported_variables.h"
#include "cs_syr_coupling.h"
#include "cs_thermal_model.h"
#include "cs_theta_scheme.h"
#include "cs_time_step_compute.h"
#include "cs_turbulence_htles.h"
#include "cs_turbulence_ke.h"
#include "cs_turbulence_kw.h"
#include "cs_turbulence_model.h"
#include "cs_turbulence_rij.h"
#include "cs_turbulence_sa.h"
#include "cs_turbulence_v2f.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_mass_injection.h"
#include "cs_wall_distance.h"
#include "cs_wall_condensation_1d_thermal.h"
#include "cs_wall_condensation.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_solve_all.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_solve_all.c
        Solve the Navier-Stokes equations,
        source term convection diffusion equations for scalars ... .
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

extern cs_real_t *cs_glob_ckupdc;
extern int cs_glob_ast_coupling_n_couplings;

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_atleca(void);

void
cs_f_atr1vf(void);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update previous values for variables fields
 *
 * \param[in]  itrale       ALE iteration number
 * \param[in]  n_fields     number of fields
 * \param[in]  n_cells_ext  number of cells with ghosts
 */
/*----------------------------------------------------------------------------*/

static void
_update_previous_values(const int        itrale,
                        const int        n_fields,
                        const cs_lnum_t  n_cells_ext)
{
  const int kst = cs_field_key_id_try("source_term_id");
  const int kstprv = cs_field_key_id_try("source_term_prev_id");
  const int key_buoyant_id = cs_field_key_id_try("coupled_with_vel_p");

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *fld = cs_field_by_id(f_id);

    if (!(fld->type & CS_FIELD_VARIABLE))
      continue;

    if (fld->type & CS_FIELD_CDO)
      continue;

    if (   fld != CS_F_(p)
        || cs_glob_velocity_pressure_model->idilat != 2)
      cs_field_current_to_previous(fld);

    // For buoyant scalar with source termes, current to previous for them
    const int st_prv_id = cs_field_get_key_int(fld, kstprv);
    const int coupled_with_vel_p_fld = cs_field_get_key_int(fld, key_buoyant_id);
    if (itrale <= 1 || st_prv_id < 0 || coupled_with_vel_p_fld != 1)
      continue;
    const int st_id = cs_field_get_key_int(fld, kst);
    cs_array_real_copy(n_cells_ext,
                       cs_field_by_id(st_id)->val,
                       cs_field_by_id(st_prv_id)->val);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute the pseudo tensorial time step
 *        if needed for the pressure solving
 *
 * \param[in]  m        pointer to mesh structure
 * \param[in]  ncepdc   number of cellules with head lossses
 * \param[in]  icepdc   ids of cells on which to apply head losses
 * \param[in]  ckupdc   head loss coefficients at matchin cells
 */
/*----------------------------------------------------------------------------*/

static void
_compute_tensorial_time_step(const cs_mesh_t   *m,
                             cs_lnum_t          ncepdc,
                             const cs_lnum_t    icepdc[],
                             const cs_real_6_t  ckupdc[])
{
  const cs_lnum_t  n_cells = m->n_cells;

  const cs_real_t *dt = CS_F_(dt)->val;
  cs_real_6_t *dttens = (cs_real_6_t *)cs_field_by_name("dttens")->val;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (int ii = 0; ii < 3; ii++)
      dttens[c_id][ii] = dt[c_id];
    for (int ii = 3; ii < 6; ii++)
      dttens[c_id][ii] = 0.0;
  }

  // dttens = (1/dt + Kpdc)^-1
# pragma omp parallel for if (ncepdc > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < ncepdc; ii++) {
    const cs_lnum_t c_id = icepdc[ii];
    const cs_real_t hdls[6] = {ckupdc[ii][0] + 1.0/dt[c_id],
                               ckupdc[ii][1] + 1.0/dt[c_id],
                               ckupdc[ii][2] + 1.0/dt[c_id],
                               ckupdc[ii][3],
                               ckupdc[ii][4],
                               ckupdc[ii][5]};

    cs_math_sym_33_inv_cramer(hdls, dttens[c_id]);
  }

  if (m->halo != NULL) {
    cs_mesh_sync_var_sym_tens(dttens);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief solve energy and variables equations when
 *        scalar and momentum they are coupled in case of buoyancy
 *
 * \param[in]  n_scal      number of scalar
 * \param[in]  iterns      number of iteration
 * \param[in]  n_cells     number of cells
 * \param[in]  scalar_idx  idx of scalar
 */
/*----------------------------------------------------------------------------*/

static void
_solve_coupled_vel_p_variables_equation(const int        n_scal,
                                        const int        iterns,
                                        const cs_lnum_t  n_cells,
                                        const int        scalar_idx[])
{
  if (cs_glob_velocity_pressure_model->n_buoyant_scal < 1)
    return;

  const cs_equation_param_t *eqp_vel
    = cs_field_get_equation_param_const(CS_F_(vel));

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  SOLVING ENERGY AND SCALARS EQUATIONS\n"
         "  ====================================\n\n"));
  }

  // Update coupled with dynamic scalar(s)
  cs_solve_transported_variables(iterns);

  // Diffusion terms for weakly compressible algorithm
  if (cs_glob_velocity_pressure_model->idilat > 3)
    cs_dilatable_scalar_diff_st(iterns);

  /* Update the density and turbulent viscosity
     ------------------------------------------ */

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  COMPUTATION OF PHYSICAL QUANTITIES\n"
         "  ==================================\n\n"));
  }

  // Disable solid cells in fluid_solid mode
  if (cs_glob_velocity_pressure_model->fluid_solid)
    cs_porous_model_set_has_disable_flag(1);
  cs_physical_properties_update(iterns);

  // Correct the scalar to ensure scalar conservation
  const cs_real_t *crom = CS_F_(rho)->val;
  const int key_buoyant_id = cs_field_key_id_try("coupled_with_vel_p");
  // Correction only made for the collocated time-scheme (Li Ma phd)
  cs_field_t *rho_mass = cs_field_by_name_try("density_mass");
  if (   rho_mass != NULL
      && cs_glob_velocity_pressure_param->itpcol == 1) {
    for (int ii = 0; ii < n_scal; ii++) {
      const cs_field_t *f = cs_field_by_id(scalar_idx[ii]);
      if ((f->type & CS_FIELD_CDO))
        continue;
      const int coupled_with_vel_p_fld
        = cs_field_get_key_int(f, key_buoyant_id);
      if (coupled_with_vel_p_fld != 1)
        continue;
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        f->val[c_id] = f->val[c_id]*rho_mass->val[c_id]/crom[c_id];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Save of the pressure and temperature (if we have internal energy)
 *
 * \param[in]  n_cells  number of cells
 */
/*----------------------------------------------------------------------------*/

static void
_update_pressure_temperature(cs_lnum_t n_cells)
{
  if (cs_glob_velocity_pressure_model->idilat != 2)
    return;

  cs_real_t *cvar_pr = CS_F_(p)->val;
  const cs_real_t *cvara_pr = CS_F_(p)->val_pre;

  cs_real_t c = 1.0;
  if (cs_glob_time_scheme->time_order == 2)
    c = 2.0;

  // Saving pressure corrected in resopv as ancient pressure
  cs_field_current_to_previous(CS_F_(p));

  if (   cs_glob_thermal_model->thermal_variable
      == CS_THERMAL_MODEL_INTERNAL_ENERGY) {
    cs_field_t *temp = cs_field_by_name_try("temperature");
    if (temp != NULL)
      cs_array_real_copy(n_cells, temp->val, temp->val_pre);
  }

  /* Saving the thermodynamic pressure at time n+1
   * coherent with the equation of state */
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    cvar_pr[c_id] = c * cvar_pr[c_id] - (c - 1) * cvara_pr[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loop on all solver equation except turbulence.
 *
 * \param[out] vel_verbosity  verbosity for velocity
 * \param[out] italim         implicit coupling iteration number
 * \param[out] itrfin         indicator for last iteration of implicit couplin
 * \param[out] ineefl         for ALE
 * \param[out] itrfup         indication of iteration
 * \param[out] must_return    if it is done
 */
/*----------------------------------------------------------------------------*/

static void
_solve_most(int              n_var,
            int              n_scal,
            int              isvhb,
            int              itrale,
            int              vel_verbosity,
            int             *italim,
            int             *itrfin,
            int             *ineefl,
            int             *itrfup,
            bool            *must_return,
            const int        scalar_idx[],
            const cs_real_t  ckupdc[][6],
            cs_real_t        htot_cond[])
{
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_field_t *th_f = cs_thermal_model_field();
  bool _must_return = *must_return;

  int *isostd = NULL;
  cs_real_t *hbord = NULL;
  cs_real_t *theipb = NULL;
  cs_real_t *visvdr = NULL;

  CS_MALLOC_HD(isostd, n_b_faces+1, int, cs_alloc_mode);

  /* Loop on cs_solve_navier_stokes for speed/pressure coupling
   * we stop at ntrup  or when we have converged itrfup equal zero
   * indicates that we need to redo an iteration for Syrthes, T1D or radiation. */
  *itrfup = 1;

  if (   cs_glob_time_scheme->isno2t > 0
      || cs_glob_velocity_pressure_param->nterup > 1)
    if (   cs_syr_coupling_n_couplings() > 0
        || cs_get_glob_1d_wall_thermal()->nfpt1t > 0
        || cs_glob_rad_transfer_params->type > 0)
      *itrfup = 0;

  if (*italim == 1)
    cs_field_build_bc_codes_all();

  if (isvhb > -1)
    BFT_MALLOC(hbord, n_b_faces, cs_real_t);

  if (th_f != NULL)
    BFT_MALLOC(theipb, n_b_faces, cs_real_t);

  if (   cs_glob_turb_model->itytur == 4
      && cs_glob_turb_les_model->idries == 1)
    BFT_MALLOC(visvdr, n_cells_ext, cs_real_t);

  int icvrge = 0, inslst = 0, iterns = 1;

  const cs_wall_condensation_t *wall_cond = cs_glob_wall_condensation;
  // Total number of cells with condensation source term
  const cs_gnum_t nftcdt = (cs_gnum_t)wall_cond->nfbpcd;

  while (iterns <= cs_glob_velocity_pressure_param->nterup) {

    // Call user BCs and computes BC coefficients
    cs_boundary_conditions_set_coeffs(n_var,
                                      iterns,
                                      isvhb,
                                      itrale,
                                      *italim,
                                      *itrfin,
                                      *ineefl,
                                      *itrfup,
                                      isostd,
                                      visvdr,
                                      hbord,
                                      theipb,
                                      nftcdt);

    if (nftcdt > 0) {
      cs_real_t *coefap = th_f->bc_coeffs->a;
      cs_real_t *cofafp = th_f->bc_coeffs->af;
      cs_real_t *cofbfp = th_f->bc_coeffs->bf;

      /* Pass the heat transfer computed by the Empiric laws
       * of the COPAIN condensation to impose the heat transfer
       * at the wall due to condensation for the enthalpy scalar. */
      for (cs_lnum_t ii = 0 ; ii < wall_cond->nfbpcd; ii++) {
        const cs_lnum_t iz = wall_cond->izzftcd[ii];
        const cs_lnum_t face_id = wall_cond->ifbpcd[ii];

        /* Enthalpy Boundary condition associated
           to the heat transfer due to condensation. */
        cofafp[face_id] = -wall_cond->hpcond[ii]*coefap[face_id];
        cofbfp[face_id] =  wall_cond->hpcond[ii];
        if (wall_cond->iztag1d[iz] == 2)
          hbord[face_id] = htot_cond[ii];
      }
    }

    /* Ground-atmosphere interface
       --------------------------- */

    /* FIXME why only we have atmo humid ?
       Deardorff force-restore model */
    if (   cs_glob_atmo_option->soil_model == 1
        && cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
      cs_soil_model();

    /* After coefficient are computed, we can easily deduce the terms
       to send for boundaries coupling (such as with Syrthes) */
    if (*itrfin == 1 && *itrfup == 1) {
      if (isvhb > -1)
        cs_syr_coupling_send_boundary(hbord, theipb);

      if (th_f != NULL && cs_glob_1d_wall_thermal->nfpt1t > 0) {
        cs_boundary_conditions_coupling_t_out(hbord, theipb);

        if (   wall_cond->icondb == 0
            || cs_glob_rad_transfer_params->type > 0)
          cs_boundary_conditions_coupling_t_in();
      }

      // 1-D thermal model coupling with condensation
      if (nftcdt > 0 && wall_cond->nztag1d == 1)
        cs_wall_condensation_1d_thermal_compute_temperature();

      // 0-D thermal model coupling with condensation
      if (wall_cond->icondv == 0)
        cs_wall_condensation_0d_thermal_solve();
    }

    /* Compute wall distance
     * TODO it has to be moved before cs_physical_properties_update, for
     * that bc types have to be known (itypfb)
     *
     * (New algorithm. the old one is in cs_boundary_condition_set_coeffs)
     * In ALE, this computation is done only for the first step */

    if (*italim == 1) {
      /*  Wall distance is computed if:
       *   - it has to be updated
       *   - we need it
       * In case there is no wall, distance is a big value. */
      if (   cs_glob_wall_distance_options->need_compute == 1
          && cs_glob_wall_distance_options->is_up_to_date == 0) {

        if (cs_glob_wall_distance_options->method != 2)
          cs_wall_distance(iterns);
        else if (cs_glob_wall_distance_options->method == 2)
          cs_wall_distance_geometric();
        // Wall distance is not updated except if ALE is switched on
        if (cs_glob_ale == CS_ALE_NONE)
          cs_get_glob_wall_distance_options()->is_up_to_date = 1;
      }
    }

    /* Compute y+ if needed and Van Driest damping */
    if (visvdr != NULL)
      cs_wall_distance_yplus(visvdr);

    if (   cs_restart_present() == 0
        && cs_glob_time_step->nt_prev == cs_glob_time_step->nt_max)
      _must_return = true;

    if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_FROZEN_CONTINUOUS_PHASE)
      _must_return = true;

    if (cs_glob_ale > CS_ALE_NONE)
      if (itrale == 0) {
        _must_return = true;
        cs_ale_solve_mesh_velocity(iterns);
      }

    if (_must_return) {
      BFT_FREE(hbord);
      BFT_FREE(theipb);
      BFT_FREE(visvdr);
      CS_FREE_HD(isostd);
      cs_field_free_bc_codes_all();
      *must_return = _must_return;
      return;
    }

    /* Compute velocity when not frozen:
       - We solve velocity and turbulence
       - We assume that all phases are frozen, or non are.
       --------------------------------------------------- */

    cs_time_control_t *vp_tc
      = &(cs_get_glob_velocity_pressure_param()->time_control);
    const cs_time_step_t *ts = cs_glob_time_step;
    bool _active_dyn = cs_time_control_is_active(vp_tc, ts);
    if (_active_dyn) {

      /* Solve momentum and mass equation
         -------------------------------- */

      // In case of buoyancy, scalars and momentum are coupled
      _solve_coupled_vel_p_variables_equation(n_scal,
                                              iterns,
                                              n_cells,
                                              scalar_idx);

      if (vel_verbosity > 0) {
        bft_printf
          (_(" ------------------------------------------------------------\n\n"
             "  SOLVING NAVIER-STOKES EQUATIONS (sub iter: %d)\n"
             "  ===============================\n\n"),
           iterns);
      }

      cs_solve_navier_stokes(iterns,
                             &icvrge,
                             itrale,
                             isostd,
                             ckupdc);

      if (   cs_glob_time_scheme->istmpf == 2
          && cs_glob_velocity_pressure_param->itpcol == 1)
        cs_theta_scheme_update_var(3);

      // If is the last iteration : inslst = 1
      if (   icvrge == 1
          || iterns == cs_glob_velocity_pressure_param->nterup) {
        /* If we need to do a new iteration for SYRTHES,
         * radiation, 1D thermal wall...
         * and that we are at the last iteration in ALE!

         *...then, we reset the convergence indicators to zero */
        if (*itrfup == 0 && *itrfin == 1) {
          *itrfup = 1;
          icvrge = 0;
          iterns--;
        }
        else
          inslst = 1;

        // For explicit mass flux
        if (cs_glob_time_scheme->istmpf == 0 && inslst == 0)
          cs_theta_scheme_update_var(3);
      }

    } // End velocity computation
    if (inslst == 1)
      break;

    iterns ++;

  }  // end while

  _update_pressure_temperature(n_cells);

  const cs_equation_param_t *eqp_vel
    = cs_field_get_equation_param_const(CS_F_(vel));

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "\n"
         "  COMPUTATION OF CFL AND FOURIER\n"
         "  ==============================\n\n"));
  }

  cs_courant_fourier_compute();

  *must_return = _must_return;

  BFT_FREE(hbord);
  BFT_FREE(theipb);
  BFT_FREE(visvdr);
  CS_FREE_HD(isostd);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve all turbulence equations.
 *
 * \param[in]  n_cells      number of cells
 * \param[in]  n_cells_ext  number of cells with ghost
 * \param[in]  verbosity    verbosiy of velocity
 */
/*----------------------------------------------------------------------------*/

static void
_solve_turbulence(cs_lnum_t   n_cells,
                  cs_lnum_t   n_cells_ext,
                  int         verbosity)
{
  if (   verbosity > 0
      && (   cs_glob_turb_model->itytur == 2
          || cs_glob_turb_model->itytur == 3
          || cs_glob_turb_model->itytur == 5
          || cs_glob_turb_model->iturb == CS_TURB_K_OMEGA)) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  SOLVING TURBULENT VARIABLES EQUATIONS\n"
         "  =====================================\n\n"));
  }

  if (   cs_glob_turb_model->itytur == 2
      || cs_glob_turb_model->itytur == 5) {
    cs_real_t *prdv2f = NULL;
    if (cs_glob_turb_model->itytur == 5)
      BFT_MALLOC(prdv2f, n_cells_ext, cs_real_t);
    cs_turbulence_ke(-1, prdv2f);

    if (cs_glob_turb_model->itytur == 5)
      cs_turbulence_v2f(prdv2f);

    BFT_FREE(prdv2f);

    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_ep = CS_F_(eps)->val;
    const cs_real_t *cvara_k = CS_F_(k)->val_pre;
    const cs_real_t *cvara_ep = CS_F_(eps)->val_pre;

    if (   (cs_glob_turb_rans_model->ikecou == 0)
        && (cs_glob_time_step_options->idtvar > CS_TIME_STEP_STEADY)) {

      const cs_real_t relaxk
        = cs_field_get_equation_param_const(CS_F_(k))->relaxv;
      const cs_real_t relaxe
        = cs_field_get_equation_param_const(CS_F_(eps))->relaxv;
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cvar_k[c_id] = relaxk*cvar_k[c_id] + (1.0-relaxk)*cvara_k[c_id];
        cvar_ep[c_id] = relaxe*cvar_ep[c_id] + (1.0-relaxe)*cvara_ep[c_id];
      }
    }

    // HTLES
    if (cs_glob_turb_model->hybrid_turb == 4)
      cs_turbulence_htles();

  }
  else if (cs_glob_turb_model->itytur == 3) {
    if (cs_glob_turb_model->iturb == CS_TURB_RIJ_EPSILON_EBRSM)
      cs_turbulence_rij_solve_alpha(CS_F_(alp_bl)->id, -1, cs_turb_xcl);
    cs_turbulence_rij(-1);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {
    cs_turbulence_kw(-1);
    cs_real_t *cvar_k = CS_F_(k)->val;
    cs_real_t *cvar_omg = CS_F_(omg)->val;
    const cs_real_t *cvara_k = CS_F_(k)->val_pre;
    const cs_real_t *cvara_omg = CS_F_(omg)->val_pre;
    if (   cs_glob_turb_rans_model->ikecou == 0
        && cs_glob_time_step_options->idtvar > CS_TIME_STEP_STEADY) {
      const cs_real_t relaxk
        = cs_field_get_equation_param_const(CS_F_(k))->relaxv;
      const cs_real_t relaxw
        = cs_field_get_equation_param_const(CS_F_(omg))->relaxv;
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cvar_k[c_id]   = relaxk*cvar_k[c_id] + (1.0-relaxk)*cvara_k[c_id];
        cvar_omg[c_id] = relaxw*cvar_omg[c_id] + (1.0-relaxw)*cvara_omg[c_id];
      }
    }

    // HTLES
    if (cs_glob_turb_model->hybrid_turb == 4)
      cs_turbulence_htles();
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {
    cs_turbulence_sa();
    cs_real_t *cvar_nusa = CS_F_(nusa)->val;
    const cs_real_t *cvara_nusa = CS_F_(nusa)->val_pre;
    if (cs_glob_time_step_options->idtvar > CS_TIME_STEP_STEADY) {
      const cs_real_t relaxn
        = cs_field_get_equation_param_const(CS_F_(nusa))->relaxv;
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        cvar_nusa[c_id]
          = relaxn*cvar_nusa[c_id]+(1.0-relaxn)*cvara_nusa[c_id];
    }
  }
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Resolution of incompressible Navier Stokes, scalar transport
 *        equations... for a time step.
 *
 * \param[in]     itrale        ALE iteration number
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_all(int  itrale)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  const cs_equation_param_t *eqp_p
    = cs_field_get_equation_param_const(CS_F_(p));

  const cs_equation_param_t *eqp_vel
    = cs_field_get_equation_param_const(CS_F_(vel));

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;

  /* Storage indicator of a scalar and its associated exchange coefficiento.
   * For the moment, we only store in the SYRTHES coupling case, or with the
   * 1D wall thermal module, and assume a single thermal scalar is coupled. */
  int isvhb = -1;

  /* storage id of scalar fields */
  int n_var = 0;
  int n_scal = 0;
  const int n_fields = cs_field_n_fields();

  int scalar_idx[n_fields];
  const int keysca = cs_field_key_id("scalar_id");

  {
    const cs_field_t *f_th = cs_thermal_model_field();
    const int kcpsyr = cs_field_key_id("syrthes_coupling");
    const int n_syr_couplings = cs_syr_coupling_n_couplings();
    const int nfpt1d = cs_get_glob_1d_wall_thermal()->nfpt1t;

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE))
        continue;
      if (f->type & CS_FIELD_CDO)
        continue;
      n_var++;
      const int sc_id = cs_field_get_key_int(f, keysca) - 1;
      if (sc_id < 0)
        continue;
      scalar_idx[n_scal] = f_id;
      n_scal++;
      if (n_syr_couplings > 0) {
        if (cs_field_get_key_int(f, kcpsyr) == 1)
          isvhb = f_id;
      }
      else if (nfpt1d > 0) {
        if (f == f_th) {
          isvhb = f_id;
        }
      }
    }
  }

  cs_real_t *cvar_pr = CS_F_(p)->val;

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  INITIALIZATIONS\n"
         "  ===============\n\n"));
  }

  // Compute z ground
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] != CS_ATMO_OFF)
    cs_atmo_z_ground_compute();

  if (   (   cs_glob_thermal_model->thermal_variable
          == CS_THERMAL_MODEL_TEMPERATURE)
      || (   cs_glob_thermal_model->thermal_variable
          == CS_THERMAL_MODEL_INTERNAL_ENERGY))
    cs_thermal_model_init();

  /* At the beginning of computation we reset the pressure
     ----------------------------------------------------- */

  /* We do this over an infinite period of time, because often the mass flux field
   *   initial is not zero divergence (CL included) and obtaining
   *   of a flow with zero divergence consistent with the stationary constraint
   *   may take a few steps.
   * Note that the pressure is taken up in the Stokes stage.
   * We do not do this in the case of taking into account the pressure
   *   hydrostatic, nor in the case of compressible */

  if (   (cs_restart_present() == 0)
      && (cs_glob_time_step->nt_cur <= cs_glob_time_step->nt_ini)
      && (cs_glob_velocity_pressure_param->iphydr != 1)
      && (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
      && (cs_glob_velocity_pressure_model->idilat < 2)) {

    if (eqp_p->verbosity > 1)
      bft_printf("Reinitialization of pressure at iteration %d\n\n",
                 cs_glob_time_step->nt_cur);
    const cs_real_t *xyzp0 = fp->xyzp0;
    const cs_real_t *gxyz = cs_glob_physical_constants->gravity;
    const cs_real_t pred0 = fp->pred0;
    const cs_real_t p0 = fp->p0, ro0 = fp->ro0;
    cs_real_t *cpro_prtot = cs_field_by_name("total_pressure")->val;
#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      cvar_pr[c_id] = pred0;
      cpro_prtot[c_id] = p0 + ro0*cs_math_3_distance_dot_product(xyzp0,
                                                                 cell_cen[c_id],
                                                                 gxyz);
    }
  }

  /* Halo synchronization (only variables require this) */
  if (m->halo != NULL) {
    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (!(f->type & CS_FIELD_VARIABLE))
        continue;

      if ((f->type & CS_FIELD_CDO))
        continue;

      cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD, f->val, f->dim);
      if (m->have_rotation_perio) {
        if (f->dim == 3)
          cs_halo_perio_sync_var_vect(m->halo,  CS_HALO_STANDARD, f->val, 3);
        else if (f->dim == 6)
          cs_halo_perio_sync_var_vect(m->halo,  CS_HALO_STANDARD, f->val, 6);
        else if (f->dim > 1)
          bft_error(__FILE__, __LINE__, 0,
                    _("field %s of dimension %d\n"
                      "cannot handle rotational periodicity\n"),
                    f->name, f->dim);
      }
    }

    static bool first_pass = true;
    if (first_pass) {
      if (   cs_glob_velocity_pressure_param->iphydr == 1
          && cs_field_by_name_try("volume_forces") != NULL) {
        cs_real_t *frcxt = cs_field_by_name_try("volume_forces")->val;
        cs_halo_sync_var_strided(m->halo, CS_HALO_STANDARD, frcxt, 3);
        cs_halo_perio_sync_var_vect(m->halo,  CS_HALO_STANDARD, frcxt, 3);
      }

      if (   cs_glob_velocity_pressure_param->icalhy == 1
          || cs_glob_velocity_pressure_model->idilat == 3)
        cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, CS_F_(rho)->val);
      first_pass = false;
    }

  }

  /* Temporal update of previous values (mass flux, density, ...)
   *  ------------------------------------------------------------
   *  We exit before SCHTMP as otherwise for 2nd order in time the value of the
   *  mass flux at the previous time is overwritten by the value at the current
   *  time step.
   *  When ntmabs = 0, there is no issue since all mass fluxes are 0. */

  if (   cs_restart_present()
      && cs_glob_time_step->nt_prev == cs_glob_time_step->nt_max)
    return;

  // If itrale = 0, we are initializing ALE; do not touch the mass flux either.
  if (itrale > 0)
    cs_theta_scheme_update_var(1);

  /* Update location of code_saturne/code_saturne coupling interfaces
     ---------------------------------------------------------------- */

  if (cs_sat_coupling_n_couplings() > 0)
    cs_sat_coupling_locate_all();

  /* Compute physical quantities, when
     - They are time-varying
     - They may change upon a compuitation restart
     --------------------------------------------- */

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  COMPUTATION OF PHYSICAL QUANTITIES\n"
         "  ==================================\n\n"));
  }

  if (cs_glob_velocity_pressure_model->fluid_solid)
    cs_porous_model_set_has_disable_flag(1);

  cs_physical_properties_update(-1);

  if (itrale > 0)
    cs_theta_scheme_update_var(2);

  /* Compute head loss coeffs.
     we do it even if there is no head loss on the local rank in case
     a user-defiend function requires collective operations such as
     computing a min or max value of a variable. */

  const int ncpdct = cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_HEAD_LOSS);
  const cs_lnum_t ncepdc = cs_volume_zone_n_type_cells(CS_VOLUME_ZONE_HEAD_LOSS);

  cs_lnum_t *icepdc = NULL;
  cs_real_6_t *ckupdc = NULL;

  if (ncpdct > 0) {
    BFT_MALLOC(icepdc, ncepdc, cs_lnum_t);
    cs_volume_zone_select_type_cells(CS_VOLUME_ZONE_HEAD_LOSS, icepdc);

    BFT_MALLOC(ckupdc, ncepdc, cs_real_6_t);

    cs_head_losses_compute(ckupdc);

    if (cs_glob_lagr_reentrained_model->iflow == 1)
      cs_lagr_head_losses(n_cells,
                          icepdc,
                          cs_boundary_conditions_get_bc_type(),
                          ckupdc);
  }

  /* Current to previous for variables
     --------------------------------- */

  _update_previous_values(itrale, n_fields, n_cells);

  /* Compute time step if variable
     ----------------------------- */

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  COMPUTATION OF CFL, FOURIER AND VARIABLE DT\n"
         "  ===========================================\n\n"));
  }

  cs_local_time_step_compute(itrale);
  const int nalinf = cs_glob_ale_n_ini_f;
  const int nbaste = cs_glob_ast_coupling_n_couplings;
  if (nbaste > 0 && itrale > nalinf)
    cs_ast_coupling_exchange_time_step(CS_F_(dt)->val);

  if (eqp_p->idften & CS_ANISOTROPIC_DIFFUSION)
    _compute_tensorial_time_step(m,
                                 ncepdc,
                                 icepdc,
                                 ckupdc);

  BFT_FREE(icepdc);

  /* Adjustment of Pth pressure and rho
   * volume mass for the variable density algorithm */

  if (   cs_glob_fluid_properties->ipthrm == 1
      || cs_glob_velocity_pressure_model->idilat == 3)
    cs_compute_thermo_pressure_density();

  /* Setup boundary conditions
     ------------------------- */

  if (eqp_vel->verbosity > 0) {
    bft_printf
      (_(" ------------------------------------------------------------\n\n"
         "  SETTING UP THE BOUNDARY CONDITIONS\n"
         "  ==================================\n\n"));
  }

  /* ALE method: start of loop for implying the movement of structures.
                 itrfin=0 indicates that we need to redo an iteration
                 for Syrthes, T1D or radiation. */
  int itrfup = 1;
  int italim = 1;
  int itrfin = 1;
  int ineefl = 0;

  if (   cs_glob_ale >  CS_ALE_NONE
      && itrale > nalinf
      && cs_glob_mobile_structures_i_max > 1) {
    /* Indicate if we need to return to the initial state at the end
       of an  ALE iteration. */
    ineefl = 1;

    if (   cs_syr_coupling_n_couplings() > 0
        || cs_get_glob_1d_wall_thermal()->nfpt1t > 0
        || cs_glob_rad_transfer_params->type > 0)
      itrfin = 0;
  }

  /* Evaluate mass source term coefficients
     (called on all ranks in case user calls global operations). */

  if (cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_MASS_SOURCE_TERM) > 0)
    cs_volume_mass_injection_eval();

  /* Fill the condensation arrays spcond for the sink term of condensation
   * and hpcond the thermal exchange coefficient associated to the phase
   * change (gas phase to liquid phase)
   * ---------------------------------------------------------------------- */

  cs_real_t *htot_cond = NULL;
  cs_wall_condensation_t *wall_cond = cs_get_glob_wall_condensation();
  if ((wall_cond->icondb == 0) || (wall_cond->icondv == 0)) {
    // Condensation source terms arrays initialized
    for (int ii = 0; ii < wall_cond->nfbpcd; ii++) {
      wall_cond->hpcond[ii] = 0.0;
      for (int ivar = 0; ivar < n_var; ivar++) {
        wall_cond->itypcd[ivar*wall_cond->nfbpcd + ii] = 0;
        wall_cond->spcond[ivar*wall_cond->nfbpcd + ii] = 0.0;
      }
    }

    for (cs_lnum_t ii = 0; ii < wall_cond->ncmast; ii++) {
      wall_cond->flxmst[ii] = 0.0;
      for (int ivar = 0; ivar < n_var; ivar++) {
        wall_cond->itypst[ivar*wall_cond->ncmast + ii] = 0;
        wall_cond->svcond[ivar*wall_cond->ncmast + ii] = 0.0;
      }
    }

    cs_user_wall_condensation(3);

    // Use empiric correlations to compute heat
    // and mass transfer due to wall condensation
    BFT_MALLOC(htot_cond, wall_cond->nfbpcd, cs_real_t);
    cs_wall_condensation_compute(htot_cond);
  }

  bool must_return = false;
  bool need_new_solve = true;
  cs_time_control_t *vp_tc
    = &(cs_get_glob_velocity_pressure_param()->time_control);
  const cs_time_step_t *ts = cs_glob_time_step;
  bool _active_dyn = cs_time_control_is_active(vp_tc, ts);

  while (need_new_solve) {

    _solve_most(n_var,
                n_scal,
                isvhb,
                itrale,
                eqp_vel->verbosity,
                &italim,
                &itrfin,
                &ineefl,
                &itrfup,
                &must_return,
                scalar_idx,
                ckupdc,
                htot_cond);

    need_new_solve = false;

    if (must_return) {
      BFT_FREE(htot_cond);
      BFT_FREE(ckupdc);
      return;
    }

    /* Computation on non-frozen velocity field, continued */

    if (_active_dyn && cs_glob_ale > CS_ALE_NONE) {

      /* Movement of structures in ALE and test implicit loop */

      const int n_structs = cs_mobile_structures_get_n_structures();
      if (n_structs > 0 || nbaste > 0) {
        cs_mobile_structures_displacement(itrale, italim, &itrfin);
        if (itrfin != -1) {
          italim++;
          need_new_solve = true;
        }
      }

    }

  } // End loop on need_new_solve (_solve_most)

  /* Computation on non-frozen velocity field, continued */

  if (_active_dyn) {

    // We pass in cs_theta_scheme_update_var only in explicit
    if (cs_glob_time_scheme->istmpf == 0)
      cs_theta_scheme_update_var(4);

    /* Solve turbulence
       ---------------- */

    _solve_turbulence(n_cells, n_cells_ext, eqp_vel->verbosity);

  } // end if _active_dyn

  BFT_FREE(htot_cond);

  // Re Enable solid cells in fluid_solid mode
  if (cs_glob_velocity_pressure_model->fluid_solid)
    cs_porous_model_set_has_disable_flag(0);

  /* Solve scalars
     ------------- */

  if (n_scal > 0 && cs_glob_rad_transfer_params->type > 0) {
    if (eqp_vel->verbosity > 0) {
      bft_printf
        (_(" ------------------------------------------------------------\n\n"
           "  SOLVING THERMAL RADIATIVE TRANSFER\n"
           "  ==================================\n\n"));
    }

    if (   cs_glob_atmo_option->radiative_model_1d == 1
        && cs_glob_physical_model_flag[CS_ATMOSPHERIC] > CS_ATMO_CONSTANT_DENSITY)
      cs_f_atr1vf();

    cs_rad_transfer_solve(cs_boundary_conditions_get_bc_type());
  }

  if (n_scal > 0) {
    if (eqp_vel->verbosity > 0) {
      bft_printf
        (_(" ------------------------------------------------------------\n\n"
           "  SOLVING ENERGY AND SCALARS EQUATIONS\n"
           "  ==================================\n\n"));
    }

    // Update non-buoyant scalar(s)
    cs_solve_transported_variables(-1);

    // Diffusion terms for weakly compressible algorithm
    if (cs_glob_velocity_pressure_model->idilat > 3)
      cs_dilatable_scalar_diff_st(-1);
  }

  cs_field_free_bc_codes_all();
  BFT_FREE(ckupdc);

  /* Handle mass flux, viscosity, density, and specific heat for theta-scheme
     ------------------------------------------------------------------------ */

  cs_theta_scheme_update_var(5);

  /* Update flow through fans
     ------------------------ */

  if (cs_fan_n_fans() > 0) {
    const int kimasf = cs_field_key_id("inner_mass_flux_id");
    const int kbmasf  = cs_field_key_id("boundary_mass_flux_id");
    cs_real_t *b_massflux
      = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf))->val;
    cs_real_t *i_massflux
      = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf))->val;

    cs_fan_compute_flows(m,
                         mq,
                         i_massflux,
                         b_massflux,
                         CS_F_(rho)->val,
                         CS_F_(rho_b)->val);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
