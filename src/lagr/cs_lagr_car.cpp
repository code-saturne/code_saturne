/*============================================================================
 * Compute particle characteristics: Tp, TL and PI
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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_printf.h"
#include "bft/bft_error.h"
#include "bft/bft_mem.h"

#include "fvm/fvm_periodicity.h"

#include "atmo/cs_atmo.h"
#include "base/cs_base.h"
#include "base/cs_defs.h"
#include "base/cs_field_operator.h"
#include "base/cs_halo.h"
#include "base/cs_interface.h"
#include "base/cs_math.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_random.h"
#include "base/cs_rotation.h"
#include "base/cs_search.h"
#include "base/cs_timer_stats.h"
#include "base/cs_thermal_model.h"
#include "base/cs_velocity_pressure.h"
#include "turb/cs_turbulence_model.h"

#include "base/cs_field.h"

#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_new.h"
#include "lagr/cs_lagr_particle.h"
#include "lagr/cs_lagr_stat.h"
#include "lagr/cs_lagr_precipitation_model.h"
#include "lagr/cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_car.h"

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute particle characteristics (except force_p): Tp, TL and PI
 * as well as covariance and variance tensors for the stochastic model
 *
 * \param[in] iprev                time step indicator for fields
 *                                   0: use fields at current time step
 *                                   1: use fields at previous time step
 * \param[in]  phase_id            carrier phase id
 * \param[in]  ip                  particle index in set
 * \param[in]  nor                 current step id (for 2nd order scheme)
 * \param[in]  dt_part             time step associated to the particle
 * \param[out] taup                dynamic characteristic time
 * \param[out] tlag                fluid characteristic Lagrangian time scale
 * \param[out] piil                term in integration of up sde
 * \param[out] bx                  turbulence characteristics
 * \param[out] tempct              thermal characteristic time
 * \param[out] beta                for the extended scheme
 * \param[out] vagaus              gaussian random variables
 * \param[out] br_gaus             gaussian random variables
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_car(int                iprev,
            int                phase_id,
            cs_lnum_t          ip,
            int                nor,
            const cs_real_t    dt_part,
            cs_real_t         *taup,
            cs_real_3_t        tlag,
            cs_real_3_t        piil,
            cs_real_33_t       bx,
            cs_real_2_t        tempct,
            cs_real_3_t        beta,
            cs_real_3_t       *vagaus,
            cs_real_6_t        br_gaus)
{
  /* Particles management */

  int cell_wise_integ = cs_glob_lagr_time_scheme->cell_wise_integ;

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;
  unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

  /* FIXME we may still need to do computations here */
  if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
    return;

  cs_lagr_extra_module_t *extra_i = cs_get_lagr_extra_module();
  cs_lagr_extra_module_t *extra = extra_i;
  int n_phases = extra->n_phases;

  iprev = cs::min(iprev, extra->vel->n_time_vals -1);

  /* Initialization
     ---------------*/

  bool turb_disp_model = false;
  if (   cs_glob_lagr_model->modcpl > 0
      && cs_glob_time_step->nt_cur > cs_glob_lagr_model->modcpl
      && cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt)
    turb_disp_model = true;

  cs_real_t rec    = 1000.0;
  cs_real_t c0     = cs_turb_crij_c0;
  cs_real_t cl     = 1.0 / (0.5 + 0.75 * c0);
  cs_real_t cb     = 0.8;
  cs_real_t d6spi  = 6.0 / cs_math_pi;
  cs_real_t d1s3   = 1.0 / 3.0;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_real_t diftl0 = -1;
  if (   cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] == 2)
    diftl0 = cs_field_get_key_double(cs_field_by_name("enthalpy"),
                                     cs_field_key_id("diffusivity_ref"));

  /* Compute Tp and Tc in case of thermal model
     -------------------------------------------*/

  cs_real_t p_diam   = cs_lagr_particle_get_real(particle, p_am,
                                                 CS_LAGR_DIAMETER);
  cs_lnum_t cell_id  = cs_lagr_particle_get_lnum(particle, p_am,
                                                 CS_LAGR_CELL_ID);

  cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);
  cs_real_t p_rom  = p_mass * d6spi / pow(p_diam, 3.0);

  cs_real_t rel_vel_norm = 0.;
  cs_real_t romf           = extra_i[phase_id].cromf->val[cell_id];
  cs_real_t xnul           = extra_i[phase_id].viscl->val[cell_id] / romf;
  cs_real_t *part_vel_seen =
    cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, ip, CS_LAGR_VELOCITY_SEEN);
  cs_real_t *part_vel      =
    cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, ip, CS_LAGR_VELOCITY);

  for (int idim_ = 0; idim_ < 3; idim_++) {
    rel_vel_norm += (part_vel_seen[3*phase_id + idim_] - part_vel[idim_])
                  * (part_vel_seen[3*phase_id + idim_] - part_vel[idim_]);
  }

  rel_vel_norm = sqrt(rel_vel_norm);

  /* Compute the local Reynolds number */
  cs_real_t rep = rel_vel_norm * p_diam / xnul; /* local Reynolds number */

  cs_real_t d2 = cs_math_sq(p_diam); /* drag coefficient */

  cs_real_t fdr;
  if (rep <= rec)
    fdr = 18.0 * xnul * (1.0 + 0.15 * pow (rep, 0.687)) / d2;
  else
    fdr = 0.44 * 3.0 / 4.0 * rel_vel_norm / p_diam;

  /* Tp computation */
  *taup = p_rom / romf / fdr;

  /* Added-mass term? */
  if (cs_glob_lagr_time_scheme->iadded_mass == 1)
    *taup *= (  1.0
                 + 0.5 * cs_glob_lagr_time_scheme->added_mass_const
                       * romf / p_rom);

  /* Tp user computation */

  cs_user_lagr_rt(phase_id, ip, rep, rel_vel_norm,
                  romf, p_rom, xnul, taup, dt_part);

  /* Tc computation  computed only at the first phase*/

  if (   (   (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
              && cs_glob_lagr_specific_physics->solve_temperature == 1)
         || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
         || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR)
      && phase_id == 0) {

    /* Fluid Cp */
    cs_real_t xcp;
    if (extra->icp > 0)
      xcp = extra->cpro_cp->val[0];
    else
      xcp = cs_glob_fluid_properties->cp0;

    /* Local Nu computation */
    /* a priori in gas or pulverized coal combustion,
       diffusivity is always constant */
    cs_real_t xrkl;
    if (diftl0 >= 0)
      xrkl = diftl0 / romf;
    else if (extra->cpro_viscls != nullptr)
      xrkl = extra->cpro_viscls->val[cell_id] / (romf * xcp);
    else
      xrkl = extra->visls0 / (romf * xcp);

    cs_real_t prt  = xnul / xrkl;
    cs_real_t fnus = 2.0 + 0.55 * pow (rep, 0.5) * pow (prt, (d1s3));

    cs_real_t p_cp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP);

    /* Thermal characteristic time Tc computation */
    tempct[0] = d2 * p_rom * p_cp / (fnus * 6.0 * romf * xcp * xrkl);

    /* User computation for Tc */
    cs_user_lagr_rt_t(ip, rep, rel_vel_norm, romf, p_rom,
                      xnul, xcp, xrkl, tempct, dt_part);

    /* Implicit source term for return thermal coupling */
    tempct[1] = fnus * cs_math_pi * p_diam * xrkl * romf;
  }

  /* Compute TL
     ---------- */

  /* Compute turbulent relaxation time and diffusion term associated to the
   * particle and based on turbulence model */

  if (cs_glob_lagr_model->idistu == 1) {
    if (phase_id == 0) {
      /* -> Stochastic draws are made for all phase on the first loop */
      if(nor == 1) {
        cs_random_normal(3 * (n_phases + 2), (cs_real_t*)vagaus);
        if (cs_glob_lagr_time_scheme->t_order== 2) {
          /* save vagaus */
          cs_real_t *_v_gaus =
            cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                     CS_LAGR_V_GAUSS);
          for (cs_lnum_t i = 0; i < n_phases + 2; i++) {
            for (cs_lnum_t id = 0; id < 3; id++)
               _v_gaus[3 * i + id] = vagaus[i][id];
          }
        }
      }
      else {
        /* Get previously drawn _v_gaus */
        cs_real_t *_v_gaus =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
                                                   CS_LAGR_V_GAUSS);
        for (cs_lnum_t i = 0; i < n_phases + 2; i++) {
          for (cs_lnum_t id = 0; id < 3; id++)
            vagaus[i][id] = _v_gaus[3 * i + id];
        }
      }
    }

    /* Calculation of TL BX and potentially beta */
    cs_real_t energi = extra_i[phase_id].cvar_k->val[cell_id];
    if (   extra_i[phase_id].lagr_time->val[cell_id] > cs_math_epzero
        && energi > cs_math_epzero) {
      if (turb_disp_model) {

        for (cs_lnum_t id = 0; id < 3; id++) {
          tlag[id]  =
            cs::max(extra_i[phase_id].anisotropic_lagr_time[cell_id][id],
                    cs_math_epzero);
          cs_real_t bxi = extra_i[phase_id].anisotropic_bx[cell_id][id];
          if (bxi > 0.0)
            bx[id][nor-1] = sqrt(bxi);
          else
            bx[id][nor-1] = 0.0;
        }
        /* Compute beta_i in the local referential in the first phase
         *
         * Note: the extended time scheme should
         * compute grad(Tl*_i) = grad(Tl)/b_i - Tl grad(b_i)/b_i^2
         *
         * In the following, only the first term is kept
         * to avoid computing grad(b_i) where no particles are present.
         *
         * The other possibility would be to compute b_i everywhere
         * (taking <u_pi> from the statistic "correctly" initialized)
         *
         */
        if (cs_glob_lagr_time_scheme->extended_t_scheme !=0 && phase_id == 0) {
          for (cs_lnum_t id = 0; id < 3; id++) {
            if (cs::abs(tlag[id] - *taup)
                < cs_math_epzero * *taup)
              beta[id] = 0.;
            else
              beta[id] = extra_i[phase_id].grad_lagr_time_r_et[cell_id][id]
                           * cs_math_pow2(bx[id][nor-1])
                           / (tlag[id]-*taup);
          }
        }
      }
      else { // fluid particles
        tlag[0] = extra_i[phase_id].lagr_time->val[cell_id];
        tlag[0] = cs::max(tlag[0], cs_math_epzero);

        if (cs_glob_lagr_model->idiffl == 0) {
          cs_real_t uvwdif = 0.;
          for (int id =0; id < 3; id++)
            uvwdif += cs_math_sq(part_vel_seen[id] - part_vel[id]);
          uvwdif = sqrt((3.0 * uvwdif) / (2.0 * energi));
          tlag[0] /= (1.0 + cb * uvwdif);
        }

        tlag[1] = tlag[0];
        tlag[2] = tlag[0];

        bx[0][nor-1] = sqrt (c0 * cl * energi / tlag[0]);
        bx[1][nor-1] = bx[0][nor-1];
        bx[2][nor-1] = bx[0][nor-1];
        /* Compute beta_i in the global referential */
        if (cs_glob_lagr_time_scheme->extended_t_scheme !=0 && phase_id == 0) {
          for (cs_lnum_t id = 0; id < 3; id++) {
            if (cs::abs(tlag[id] - *taup)
                < cs_math_epzero * *taup)
              beta[id] = 0.;
            else
              beta[id] = extra->grad_lagr_time[cell_id][id]
                           * cs_math_pow2(bx[id][nor-1])
                           / (tlag[id]-*taup);
          }
        }
      }
    }
    else {

      /* FIXME we may still need to do computations here */
      if (cs_lagr_particles_get_flag(p_set, ip, CS_LAGR_PART_FIXED))
        return;

      for (int id = 0; id < 3; id++) {
        tlag[id] = cs_math_epzero;
        bx[id][nor-1] = 0.0;
      }
      if (extra->grad_lagr_time_r_et != nullptr && phase_id == 0) {
        for (int id = 0; id < 3; id++)
            beta[id] = 0.;
      }
    }
  } // end if idistu == 1
  else {
    if (phase_id == 0) {
      for (cs_lnum_t i = 0; i < n_phases + 2; i++) {
        for (cs_lnum_t id = 0; id < 3; id++)
          vagaus[i][id] = 0.0;
      }
      if (extra->grad_lagr_time_r_et != nullptr) {
        for (int id = 0; id < 3; id++)
          beta[id] = 0.;
      }
    }
    for (cs_lnum_t id = 0; id < 3; id++ ) {
      tlag[id] = cs_math_epzero;
      bx[id][nor-1] = 0.0;
    }
  }

  if ( cs_glob_lagr_brownian->lamvbr == 1) {

    /* -> Stochastic draws anly in the first phase_id*/
    if(nor == 1 && phase_id == 0) {
      cs_random_normal(6, (cs_real_t*)br_gaus);
      if (cs_glob_lagr_time_scheme->t_order== 2) {
        /* Save br_gaus */
        cs_real_t *_br_gaus =
          cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
              CS_LAGR_BR_GAUSS);
        for (cs_lnum_t id = 0; id < 6; id++)
          _br_gaus[id] = br_gaus[id];
      }
    }
    else if(phase_id == 0) {
      /* Get previously drawn _br_gaus */
      cs_real_t *_br_gaus =
        cs_lagr_particle_attr_get_ptr<cs_real_t>(particle, p_am,
            CS_LAGR_V_GAUSS);
      for(int id = 0; id < 6; id++)
        br_gaus[id] = _br_gaus[id];
    }
  }
  else if (phase_id == 0) {
    for(int id = 0; id < 6; id++)
      br_gaus[id] = 0.;
  }

  /* Compute Pii
     ----------- */

  if (turb_disp_model && cs_glob_lagr_model->cs_used == 0) {
    /* Initialize piil to 0. */
    for (cs_lnum_t id = 0; id < 3; id++)
      piil[id] = 0.;

    if (   extra_i[phase_id].iturb == CS_TURB_RIJ_EPSILON_SSG
        || extra_i[phase_id].iturb == CS_TURB_K_EPSILON) {

      /* Modèle Arcen & Tanière (2009) */
      for (cs_lnum_t id = 0; id < 3; id++) {
        for (int j = 0; j < 3; j++){
          piil[id] += extra_i[phase_id].grad_cov_skp[3 * id + j][cell_id][j];
        }
      }

    }
    else {
      bft_error
        (__FILE__, __LINE__, 0,
        _("Lagrangian turbulent dispersion is not compatible with\n"
        "the selected turbulence model.\n"
        "\n"
        "Turbulent dispersion is taken into account with idistu = %d\n"
        " Activated turbulence model is %d, when only k-eps, LES, Rij-eps,\n"
        " V2f or k-omega are handled."),
        (int)cs_glob_lagr_model->idistu,
        (int)extra_i[phase_id].iturb);
    }
    for (int id = 0; id < 3; id++) {
      for (int j = 0; j < 3; j++){
        /* Instantaneous velocity model */
        piil[id] -= (part_vel_seen[3 * phase_id + j]
            - extra_i[phase_id].vel->vals[iprev][3 * cell_id + j])
          * extra_i[phase_id].grad_vel[cell_id][id][j];
      }
    }
  }
  else if (cs_glob_lagr_model->cs_used == 0) {
    for (cs_lnum_t id = 0; id < 3; id++) {
      piil[id] = - extra_i[phase_id].grad_pr[cell_id][id] / romf + grav[id];
      for (int j = 0; j < 3; j++) {
        piil[id]
          -= (part_vel_seen[3 * phase_id + j]
              - extra_i[phase_id].vel->vals[iprev][3 * cell_id + j])
          * extra_i[phase_id].grad_vel[cell_id][id][j];
      }
    }
  }

  /* Compute Pii for code_saturne models
     ----------- */
  else {
    /* Compute: II = ( -grad(P)/Rom(f) + g) */

    for (int id = 0; id < 3; id++)
      piil[id] = - extra->grad_pr[cell_id][id] / romf + grav[id];

    if (turb_disp_model) {
      cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

      /* Compute: II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Us>) + g ) */
      if (stat_w->vals[cell_wise_integ][cell_id] >
            cs_glob_lagr_stat_options->threshold) {
        /* add grad(<Vf>)*(<Up>-<Us>) if there is enough particles */

        int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

        cs_field_t *stat_vel
          = cs_lagr_stat_get_moment(stat_type,
                                    CS_LAGR_STAT_GROUP_PARTICLE,
                                    CS_LAGR_MOMENT_MEAN,
                                    0,
                                    -1);

        stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN);

        cs_field_t *stat_vel_s
          = cs_lagr_stat_get_moment(stat_type,
                                    CS_LAGR_STAT_GROUP_PARTICLE,
                                    CS_LAGR_MOMENT_MEAN,
                                    0,
                                    -1);

        for (int id = 0; id < 3; id++) {

          for (cs_lnum_t i = 0; i < 3; i++) {
            cs_real_t vpm   = stat_vel->vals[cell_wise_integ][cell_id*3 + i];
            cs_real_t vsm   = stat_vel_s->vals[cell_wise_integ][cell_id*3 + i];
            piil[id] += extra->grad_vel[cell_id][id][i] * (vpm - vsm);
          }
        }
      }
    }

    /* Add buoyancy effects based on Boussinesq approximation */
    if (    (    cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_COAL
              || cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_CTWR
              || (   cs_glob_lagr_model->physical_model == CS_LAGR_PHYS_HEAT
                  && cs_glob_lagr_specific_physics->solve_temperature_seen == 1))
         && cs_field_by_name_try("thermal_expansion") != nullptr
         && cs_glob_velocity_pressure_model->idilat == 0) {

      const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
      cs_real_t temp_ref = phys_pro->t0;
      cs_real_t temp_s   =
        cs_lagr_particles_get_real(p_set, ip, CS_LAGR_TEMPERATURE_SEEN);

      cs_real_t expansion_coef
        = cs_field_by_name("thermal_expansion")->val[cell_id];

      cs_real_t  _tkelvi = cs_physical_constants_celsius_to_kelvin;
      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
        /* potential temp at ref */
        cs_real_t pref = cs_glob_atmo_constants->ps;
        cs_real_t rair = phys_pro->r_pg_cnst;
        cs_real_t cp0 = phys_pro->cp0;
        cs_real_t rscp = rair/cp0;
        temp_ref = cs_glob_atmo_option->meteo_t0 *
          pow(pref/ cs_glob_atmo_option->meteo_psea, rscp) - _tkelvi;
      }
      else if (cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_KELVIN)
        temp_ref -= _tkelvi;

      cs_real_t buoyancy_fac = - expansion_coef * (temp_s - temp_ref);

      for (int id = 0; id < 3; id++)
        piil[id] += buoyancy_fac * grav[id];

      if (   cs_glob_lagr_time_scheme->interpol_field > 0
          && extra->grad_tempf != nullptr) {
        /* Interpolate the local hydrostatic pressure gradient so its is in
         * equillibrium with the interpolated temperature at the position of the
         * particle and not in the center of the cell */
        cs_real_t *part_coord =
          cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, ip, CS_LAGR_COORDS);
        cs_real_t *cell_cen = cs_glob_mesh_quantities->cell_cen[cell_id];
        for (int id = 0; id < 3; id++) {
          for (int i = 0; i < 3; i++)
            piil[id] += expansion_coef * grav[id] * extra->grad_tempf[cell_id][i]
                      * (part_coord[i] - cell_cen[i]);
        }
      }
    }

    /* Add Coriolis effect */
    /* FIXME Is mean Coriolis effects partially stored in grad_pr ? */
    if (cs_glob_rotation->omega > cs_math_epzero)
      cs_rotation_add_coriolis_v(cs_glob_rotation, -2., part_vel_seen,
                                 piil);

    /* Add particle back effect seen by mean fluid velocity
     * (two way coupling terms)
     * ==================================================== */
    if (   cs_glob_lagr_source_terms->ltsdyn == 1
        && cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING
        && cs_glob_time_step->nt_cur > cs_glob_time_step->nt_prev) {

        /* Add: II =  -alpha_p rho_p <(U_s -U_p)/tau_p> / (alpha_f rho_f)  */
      const cs_real_3_t *lagr_st_vel = (const cs_real_3_t *)
        cs_field_by_name("lagr_st_velocity")->vals[cell_wise_integ];

      for (int id = 0; id < 3; id++)
        piil[id] += lagr_st_vel[cell_id][id] / romf;
    }
  }

  /* Note that force_p is computed in cs_lagr_get_force_p as user may modify
   * it based on quantities associated to any phase */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute external force impacting the particle
 *
 * \param[in]  dt_part             time step associated to the particle
 * \param[in]  ip                  particle index in set
 * \param[in]  taup                dynamic characteristic time
 * \param[in]  tlag                fluid characteristic Lagrangian time scale
 * \param[in]  piil                term in integration of up sde
 * \param[in]  bx                  turbulence characteristics
 * \param[in]  tempct              thermal characteristic time
 * \param[in]  beta                for the extended scheme
 * \param[in]  tsfext              info for return coupling source terms
 * \param[in]  vagaus              gaussian random variables
 * \param[in]  force_p             user external force field (m/s^2)$
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_get_force_p(const cs_real_t    dt_part,
                    cs_lnum_t          ip,
                    cs_real_t         *taup,
                    cs_real_3_t       *tlag,
                    cs_real_3_t       *piil,
                    cs_real_33_t      *bx,
                    cs_real_t          tsfext,
                    cs_real_3_t       *vagaus,
                    cs_real_3_t        force_p)
{
  /* Management of user external force field
     ---------------------------------------- */
  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;

  const cs_real_t *grav = cs_glob_physical_constants->gravity;

  cs_lnum_t cell_id  = cs_lagr_particles_get_lnum(p_set, ip,
                                                  CS_LAGR_CELL_ID);
  cs_real_t d6spi  = 6.0 / cs_math_pi;
  cs_real_t p_diam   = cs_lagr_particles_get_real(p_set, ip,
                                                  CS_LAGR_DIAMETER);
  cs_real_t p_mass = cs_lagr_particles_get_real(p_set, ip, CS_LAGR_MASS);
  cs_real_t p_rom  = p_mass * d6spi / pow(p_diam, 3.0);
  cs_real_t romf           = extra->cromf->val[cell_id];

  for (int id = 0; id < 3; id++)
    force_p[id] = 0.;

  cs_user_lagr_ef(dt_part,
                  ip,
                  taup,
                  tlag,
                  piil,
                  bx,
                  tsfext,
                  vagaus,
                  extra->grad_pr[cell_id],
                  extra->grad_vel[cell_id],
                  p_rom,
                  force_p);

  cs_real_t added_mass_const = cs_glob_lagr_time_scheme->added_mass_const;

  /* Finalize forces on particles:
   *
   *  (- pressure gradient /p_rom +ext forces + g) . taup
   *
   * */
  if (cs_glob_lagr_time_scheme->iadded_mass == 0) {
    for (int id = 0; id < 3; id++) {
      force_p[id] = (- extra->grad_pr[cell_id][id] / p_rom
                  + grav[id] + force_p[id]);
    }
  }
  /* Added-mass term ? */
  else {
    /* TODO make force_p influenced by each phase as they may have different
     * density for added mass effects*/
    for (int id = 0; id < 3; id++) {
      force_p[id] = (- extra->grad_pr[cell_id][id] / p_rom
        * (1.0 + 0.5 * added_mass_const)
        / (1.0 + 0.5 * added_mass_const * romf / p_rom)
        + grav[id] + force_p[id]);
    }
  }
  /* Add Coriolis effect */
  /* FIXME Is mean Coriolis effects partially stored in grad_pr ? */
  if (   cs_glob_lagr_model->cs_used == 1
      && cs_glob_rotation->omega > cs_math_epzero) {

    cs_real_t *part_vel
      = cs_lagr_particles_attr_get_ptr<cs_real_t>(p_set, ip,
                                                  CS_LAGR_VELOCITY);
    cs_rotation_add_coriolis_v(cs_glob_rotation, -2., part_vel, force_p);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
