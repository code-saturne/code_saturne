/*============================================================================
 * Methods for auxiliary mean fluid quantities (Lagrangian time and gradients)
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
 * Functions dealing with auxiliary mean fluid quantities
 *============================================================================*/

#include "base/cs_defs.h"

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

#include "bft/bft_mem.h"
#include "base/cs_array.h"
#include "alge/cs_balance.h"
#include "alge/cs_face_viscosity.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_math.h"
#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_tracking.h"
#include "lagr/cs_lagr_stat.h"
#include "mesh/cs_mesh.h"
#include "base/cs_parameters.h"
#include "base/cs_physical_constants.h"
#include "pprt/cs_physical_model.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "lagr/cs_lagr_aux_mean_fluid_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond (DOXYGEN_SHOULD_SKIP_THIS) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Private function definitions
 *============================================================================*/

 /*----------------------------------------------------------------------------*/
   /*!
   * \brief Create a field name based on a radical and a phase, scalar, or a
   *        non condensable index.
   *
   * \param[in]  field_radical  pointer to a string containing the radical
   * \param[in]  index  int containing an index value
   *
   * \return field_name  pointer to a string containing the constructed field name
   */
 /*----------------------------------------------------------------------------*/

/* This function is also defined in cs_lagr.c ... we can define it in cs_lagr.h and
call it here if really needed */

static char *
_field_name_aux(const char *field_radical, const int index)
{
  char *field_name;
  if (index > -1) {
    CS_MALLOC(field_name, strlen(field_radical) + 2 + 1, char);
    sprintf(field_name, "%s_%1d", field_radical, index + 1);
  } else {
    CS_MALLOC(field_name, strlen(field_radical) + 1, char);
    sprintf(field_name, "%s", field_radical);
  }
  return field_name;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute gradient of particle covariance.
 *
 *  - particle velocity and particle velocity seen covariance
 *  - particle velocity seen variance
 *
 * \param[in]  phase_id        carrier phase id
 * \param[out] grad_cov_skp    gradient of particle velocity and
 *                             particle velocity seen covariance
 *
 * \param[out] grad_cov_sk     gradient of particle velocity seen covariance
 */
/*----------------------------------------------------------------------------*/

void
compute_particle_covariance_gradient(int          phase_id,
                                     cs_real_3_t *grad_cov_skp[9],
                                     cs_real_3_t *grad_cov_sk[6])
{
  assert (cs_glob_lagr_model->cs_used == 0);

  cs_lagr_extra_module_t *extra_i = cs_glob_lagr_extra_module;

  const cs_mesh_t *mesh = cs_glob_mesh;
  /* Compute gradients of covariance */
  /* Various initializations  */
  cs_var_cal_opt_t var_cal_opt;
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_field_get_key_struct(extra_i[phase_id].pressure,
                          key_cal_opt_id,
                          &var_cal_opt);

  /* Now compute the gradients */
  cs_real_t *f_inter_cov;
  CS_MALLOC(f_inter_cov, cs_glob_mesh->n_cells_with_ghosts, cs_real_t);

  /* Get the variable we want to compute the gradients
   * from (covariance velocity seen/velocity) */
  int stat_type =
    cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN_VELOCITY_COV);

  cs_field_t *stat_cov_skp =
    cs_lagr_stat_get_moment(stat_type,
                            CS_LAGR_STAT_GROUP_PARTICLE,
                            CS_LAGR_MOMENT_MEAN,
                            0,
                            -phase_id-1);

  for (cs_lnum_t i = 0; i < 9; i++){
    for (int iel_ = 0; iel_ < mesh->n_cells; iel_++){
      f_inter_cov[iel_] = stat_cov_skp->val[9 * iel_ + i];
    }
    if (mesh->halo != nullptr)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_inter_cov);

    cs_gradient_scalar("intermediate cov skp gradient [Lagrangian module]",
                        CS_GRADIENT_GREEN_ITER,
                        CS_HALO_STANDARD,
                        1,
                        0,
                        0,
                        1,
                        var_cal_opt.verbosity,
                        static_cast<cs_gradient_limit_t>(var_cal_opt.imligr),
                        var_cal_opt.epsrgr,
                        var_cal_opt.climgr,
                        0,
                        nullptr,
                        f_inter_cov,
                        nullptr,
                        nullptr,
                        grad_cov_skp[i]);
  }
  /* Get the variable we want to compute the gradients from (velocity seen) */
  stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY_SEEN);

  cs_field_t *stat_cov_sk = cs_lagr_stat_get_moment(stat_type,
                                                    CS_LAGR_STAT_GROUP_PARTICLE,
                                                    CS_LAGR_MOMENT_VARIANCE,
                                                    0,
                                                    -phase_id-1);

  for (int i = 0; i < 6; i++) {
    for (int iel_ = 0; iel_ < mesh->n_cells; iel_++){
      f_inter_cov[iel_] = stat_cov_sk->val[6 * iel_ + i];
    }
    if (mesh->halo != nullptr)
      cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, f_inter_cov);

    cs_gradient_scalar("intermediate cov sk gradient [Lagrangian module]",
                        CS_GRADIENT_GREEN_ITER,
                        CS_HALO_STANDARD,
                        1,
                        0, // n_r_sweeps (number of reconstruction sweeps)
                        0,
                        1,
                        var_cal_opt.verbosity,
                        static_cast<cs_gradient_limit_t>(var_cal_opt.imligr),
                        var_cal_opt.epsrgr,
                        var_cal_opt.climgr,
                        0,
                        nullptr,
                        f_inter_cov,
                        nullptr,
                        nullptr,
                        grad_cov_sk[i]);
  }

  CS_FREE(f_inter_cov);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute anisotropic fluid quantities for complete model (modpl == 1).
 *
 *  - Anisotropic Lagragian time
 *  - Anisotropic Diffusion matrix
 *  - Anisotropic gradient of Lagrangian time in the relativ basis
 *
 * \param[in]   iprev                  time step indicator for fields
 *                                       0: use fields at current time step
 *                                       1: use fields at previous time step
 * \param[in]   phase_id               carrier phase id
 * \param[out]  anisotropic_lagr_time  anisotropic Lagragian time scale (modcpl)
 * \param[out]  anisotropic_bx         anisotropic diffusion term (if modcpl)
 * \param[out]  grad_lagr_time_r_et    anisotropic Lagrangian time gradient in
 *                                     relative basis
 * \param[in]   grad_lagr_time         Lagrangian time gradient
 *
 */
/*----------------------------------------------------------------------------*/

void
compute_anisotropic_prop(int            iprev,
                         int            phase_id,
                         cs_real_3_t   *anisotropic_lagr_time,
                         cs_real_3_t   *anisotropic_bx,
                         cs_real_3_t   *grad_lagr_time_r_et,
                         cs_real_3_t   *grad_lagr_time)

{
  int cell_wise_integ = cs_glob_lagr_time_scheme->cell_wise_integ;
  cs_lagr_extra_module_t *extra_i = cs_glob_lagr_extra_module;

  cs_real_t c0     = cs_turb_crij_c0;
  cs_real_t cl     = 1.0 / (0.5 + 0.75 * c0);
  cs_real_t cbcb   = 0.64;

  cs_real_t *energi = extra_i[phase_id].cvar_k->val;
  cs_real_t *dissip = extra_i[phase_id].cvar_ep->val;

  /* Compute anisotropic value of the lagr_time and diffusion matrix */
  /* complete model */
  cs_real_3_t dir;
  cs_real_3_t bbi;
  cs_real_t mean_uvwdif;

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


  cs_field_t *stat_w = cs_lagr_stat_get_stat_weight(0);

  for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
    if (dissip[cell_id] > 0.0 && energi[cell_id] > 0.0) {
      if (stat_w->vals[cell_wise_integ][cell_id] >
            cs_glob_lagr_stat_options->threshold) {
        /* compute mean relative velocity <Up> - <Us>*/
        for (int i = 0; i < 3; i++)
          dir[i] = stat_vel->vals[cell_wise_integ][cell_id * 3 + i]
                 - stat_vel_s->vals[cell_wise_integ][cell_id * 3 + i];

          /* Compute and store the mean relative velocity square
         * |<U_r>|^2 = |<Up>-<Us>|^2*/
        mean_uvwdif = cs_math_3_square_norm(dir);

        mean_uvwdif = (3.0 * mean_uvwdif) / (2.0 * energi[cell_id]);

        /* FIXME add proper isotropic behavior */

        /* Crossing trajectory in the direction of "<u_f>-<u_p>"
         * and in the span-wise direction */
        cs_real_t an, at;
        an = (1.0 + cbcb * mean_uvwdif);
        at = (1.0 + 4.0 * cbcb * mean_uvwdif);

        bbi[0] = sqrt(an); /* First direction, n,
                              in the local reference frame */
        bbi[1] = sqrt(at); /* Second and third direction,
                              orthogonal to n */
        bbi[2] = sqrt(at);

        /* Compute the timescale in parallel and transverse directions */
        for (int id = 0; id < 3; id++)
          anisotropic_lagr_time[cell_id][id] = extra_i[phase_id].lagr_time->val[cell_id]
                                    / bbi[id];

        /* Compute the main direction in the global reference
         * frame */
         cs_math_3_normalize(dir, dir);

        /* Compute k_tilde in each cell */
        cs_real_t ktil = 0;
        if (extra_i[phase_id].itytur == 3) {
          cs_real_t *rij = &(extra_i[phase_id].cvar_rij->vals[iprev][6*cell_id]);
          /* Note that n.R.n = R : n(x)n */
          cs_real_t rnn = cs_math_3_sym_33_3_dot_product(dir, rij, dir);
          cs_real_t tr_r = cs_math_6_trace(rij);
          // bbn * R : n(x)n + bbt * R : (1 - n(x)n)
          ktil = 3.0 * (rnn * bbi[0] + (tr_r -rnn) * bbi[1])
               / (2.0 * (bbi[0] + bbi[1] + bbi[2]));
          /* bbi[1] == bbi[2] is used */
        }
        else if (   extra_i[phase_id].itytur == 2
                 || extra_i[phase_id].itytur == 4
                 || extra_i[phase_id].itytur == 5
                 || extra_i[phase_id].iturb == CS_TURB_K_OMEGA) {
          ktil  = energi[cell_id];
        }

        for (int i = 0; i <3; i++) {
          anisotropic_bx[cell_id][i] = cl * (  (c0  * ktil )
              + ((ktil- energi[cell_id]/ bbi[i])
                * 2.0 / 3.0)) / anisotropic_lagr_time[cell_id][i];
          anisotropic_bx[cell_id][i] =
              cs::max(anisotropic_bx[cell_id][i], cs_math_epzero);
        }
        if (grad_lagr_time_r_et != nullptr) {
          cs_real_33_t trans_m;
          /* Rotate the frame of reference with respect to the
           * mean relative velocity direction.
           * This referential differs the referential associated directly
           * to the particle for non spheric particle
           * TODO extend extended scheme to non spheric particle*/

          // The rotation axis is the result of the cross product between
          // the new direction vector and the main axis.
          cs_real_t n_rot[3];
          /* the direction in the local reference frame "_r" is (1, 0, 0)
           * by convention */
          const cs_real_t dir_r[3] = {1.0, 0.0, 0.0};

          // Use quaternion (cos(theta/2), sin(theta/2) n_rot)
          // where n_rot = dir ^ dir_r normalised
          // so also       dir ^ (dir + dir_r)
          //
          // cos(theta/2) = || dir + dir_r|| / 2
          cs_real_t dir_p_dir_r[3] = {dir[0] + dir_r[0],
            dir[1] + dir_r[1],
            dir[2] + dir_r[2]};
          cs_real_t dir_p_dir_r_normed[3];
          cs_math_3_normalize(dir_p_dir_r, dir_p_dir_r_normed);

          /* dir ^(dir + dir_r) / || dir + dir_r|| = sin(theta/2) n_rot
           * for the quaternion */
          cs_math_3_cross_product(dir, dir_p_dir_r_normed, n_rot);

          /* quaternion, could be normalized afterwards
           *
           * Note that the division seems stupid but is not
           * in case of degenerated case where dir is null
           * */
          const cs_real_t euler[4] =
          {  cs_math_3_norm(dir_p_dir_r)
            / (cs_math_3_norm(dir) + cs_math_3_norm(dir_r)),
            n_rot[0],
            n_rot[1],
            n_rot[2]};

          trans_m[0][0] = 2.*(euler[0]*euler[0]+euler[1]*euler[1]-0.5);
          trans_m[0][1] = 2.*(euler[1]*euler[2]+euler[0]*euler[3]);
          trans_m[0][2] = 2.*(euler[1]*euler[3]-euler[0]*euler[2]);
          trans_m[1][0] = 2.*(euler[1]*euler[2]-euler[0]*euler[3]);
          trans_m[1][1] = 2.*(euler[0]*euler[0]+euler[2]*euler[2]-0.5);
          trans_m[1][2] = 2.*(euler[2]*euler[3]+euler[0]*euler[1]);
          trans_m[2][0] = 2.*(euler[1]*euler[3]+euler[0]*euler[2]);
          trans_m[2][1] = 2.*(euler[2]*euler[3]-euler[0]*euler[1]);
          trans_m[2][2] = 2.*(euler[0]*euler[0]+euler[3]*euler[3]-0.5);

          /* transform grad(Tl) in the local
           * reference frame */

          cs_math_33_3_product(trans_m, grad_lagr_time[cell_id],
                                        grad_lagr_time_r_et[cell_id]);
          for (int i = 0; i < 3; i++)
            grad_lagr_time_r_et[cell_id][i] /= bbi[i];
        }
      } /* end stat weight > treshold */
      else {
        if (grad_lagr_time_r_et != nullptr) {
          for (int i = 0; i < 3; i++)
            grad_lagr_time_r_et[cell_id][i] = grad_lagr_time[cell_id][i];
        }
        for (int i = 0; i < 3; i++) {
          anisotropic_bx[cell_id][i] = cl * c0  * energi[cell_id] /
                                     extra_i[phase_id].lagr_time->val[cell_id];
          anisotropic_lagr_time[cell_id][i] =
            extra_i[phase_id].lagr_time->val[cell_id];
        }
      }
    } /* end k > threshold et eps  > treshold */
    else {
      for (int i = 0; i < 3; i++) {
        anisotropic_bx[cell_id][i] = 0;
        anisotropic_lagr_time[cell_id][i] = cs_math_epzero;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute auxilary mean fluid quantities.
 *
 *  - Lagrangian time (if modcpl == 1 also anisotropic values)
 *  - gradient of total pressure
 *  - velocity gradient
 *  - temperature gradient
 *  - Lagrangian time gradient (also gradient of anisotropic values if needed)
 *  - diffusion matix
 *
 * \param[in]   iprev                  time step indicator for fields
 *                                      0: use fields at current time step
 *                                      1: use fields at previous time step
 * \param[in]   phase_id               carrier phase id
 * \param[out]  lagr_time              Lagragian time scale
 * \param[out]  grad_pr                pressure gradient
 * \param[out]  grad_vel               velocity gradient
 * \param[out]  grad_tempf             fluid temperature gradient
 * \param[out]  anisotropic_lagr_time  anisotropic Lagragian time scale (modcpl)
 * \param[out]  anisotropic_bx         anisotropic diffusion term (if modcpl)
 * \param[out]  grad_lagr_time_r_et    anisotropic Lagrangian time gradient in
 *                                     relative basis
 * \param[out]  grad_lagr_time         Lagrangian time gradient
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_aux_mean_fluid_quantities(int            iprev, // FIXME compute at current and previous ?
                                  int            phase_id,
                                  cs_field_t    *lagr_time,
                                  cs_real_3_t   *grad_pr,
                                  cs_real_33_t  *grad_vel,
                                  cs_real_3_t   *grad_tempf,
                                  cs_real_3_t   *anisotropic_lagr_time,
                                  cs_real_3_t   *anisotropic_bx,
                                  cs_real_3_t   *grad_lagr_time_r_et,
                                  cs_real_3_t   *grad_lagr_time)
{
  /* TODO compute cell properties in coherence with iprev */
  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_lagr_extra_module_t *extra_i = cs_glob_lagr_extra_module;
  cs_lagr_extra_module_t *extra = extra_i;

  const cs_real_t *grav  = cs_glob_physical_constants->gravity;

  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();

  bool turb_disp_model = false;
  if (   cs_glob_lagr_model->modcpl > 0
      && cs_glob_time_step->nt_cur > cs_glob_lagr_model->modcpl
      && cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt)
    turb_disp_model = true;

  /* Use pressure gradient of NEPTUNE_CFD if needed */
  if (cs_glob_lagr_model->cs_used == 0) {
    cs_real_t *cpro_pgradlagr = cs_field_by_name("lagr_pressure_gradient")->val;

    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
      for (cs_lnum_t id = 0; id < 3; id++)
        grad_pr[iel][id] = cpro_pgradlagr[3*iel + id];

    if (turb_disp_model || cs_glob_lagr_model->shape > 0
        || cs_glob_lagr_time_scheme->interpol_field > 0) {

      char *f_name = nullptr;
      f_name = _field_name_aux("lagr_velocity_gradient", phase_id);
      cs_real_33_t *cpro_vgradlagr
        = (cs_real_33_t *)(cs_field_by_name(f_name)->val);
      CS_FREE(f_name);

      if (cpro_vgradlagr != nullptr && grad_vel != nullptr) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t i = 0; i < 3; i++) {
            for (cs_lnum_t j = 0; j < 3; j++)
              grad_vel[c_id][i][j] = cpro_vgradlagr[c_id][i][j];
          }
        }
      }
    }
    compute_particle_covariance_gradient(phase_id,
                                         extra_i[phase_id].grad_cov_skp,
                                         extra_i[phase_id].grad_cov_sk);
  }

  if (cs_glob_lagr_model->cs_used == 1) {

    cs_real_t *wpres = nullptr;

    /* Hydrostatic pressure algorithm? */
    int hyd_p_flag = cs_glob_velocity_pressure_param->iphydr;

    cs_real_3_t *f_ext = nullptr;
    // Warning: with neptune_cfd we do not enter here
    if (hyd_p_flag == 1 && cs_glob_lagr_model->cs_used == 1)
      f_ext = (cs_real_3_t *)(cs_field_by_name("volume_forces")->val);

    cs_real_t *solved_pres = extra->pressure->val;

    /* retrieve 2/3 rho^{n} k^{n} from solved pressure field for EVM models */
    assert(turb_model != nullptr);
    if (turb_model->order <= CS_TURB_FIRST_ORDER
        && cs_glob_turb_rans_model->igrhok == 0) {
      CS_MALLOC(wpres, n_cells_with_ghosts, cs_real_t);
      int time_id = (extra->cvar_k->n_time_vals > 1) ? 1 : 0;
      const cs_real_t *cvar_k = extra->cvar_k->vals[time_id];
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        wpres[c_id] =  solved_pres[c_id]
                     -  2./3. * extra_i[phase_id].cromf->val[c_id]
                      * cvar_k[c_id];
      }
    }
    else {
      wpres = solved_pres;
    }

    /* Parameters for gradient computation
     * =================================== */

    cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;
    cs_halo_type_t halo_type = CS_HALO_STANDARD;

    /* Get the calculation option from the pressure field */

    const cs_equation_param_t *eqp =
      cs_field_get_equation_param_const(extra_i[phase_id].pressure);

    cs_gradient_type_by_imrgra(eqp->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_real_t *weight = nullptr;
    cs_internal_coupling_t  *cpl = nullptr;
    int w_stride = 1;

    if (eqp->iwgrec == 1) {
      /* Weighted gradient coefficients */
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(extra_i[phase_id].pressure, key_id);
      if (diff_id > -1) {
        cs_field_t *weight_f = cs_field_by_id(diff_id);
        weight = weight_f->val;
        w_stride = weight_f->dim;
      }
      /* Internal coupling structure */
      key_id = cs_field_key_id_try("coupling_entity");
      if (key_id > -1) {
        int coupl_id = cs_field_get_key_int(extra_i[phase_id].pressure, key_id);
        if (coupl_id > -1)
          cpl = cs_internal_coupling_by_id(coupl_id);
      }
    } else if (eqp->iwgrec == 0) {
      if (eqp->idiff > 0) {
        int key_id = cs_field_key_id_try("coupling_entity");
        if (key_id > -1) {
          int coupl_id = cs_field_get_key_int(extra_i[phase_id].pressure, key_id);
          if (coupl_id > -1)
            cpl = cs_internal_coupling_by_id(coupl_id);
        }
      }
    }

    /* Compute pressure gradient
     * ========================= */

    cs_gradient_scalar("pressure [Lagrangian module]",
                       gradient_type,
                       halo_type,
                       1,
                       eqp->nswrgr,
                       hyd_p_flag,
                       w_stride,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       f_ext,
                       extra_i[phase_id].pressure->bc_coeffs,
                       wpres,
                       weight,
                       cpl,
                       grad_pr);

    if (wpres != solved_pres)
      CS_FREE(wpres);

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
      if(cs_glob_velocity_pressure_model->idilat == 0) {
        cs_real_t *romf = extra_i[phase_id].cromf->val;
        for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++)
        {
          for (cs_lnum_t i = 0; i < 3; i++)
            grad_pr[cell_id][i] += romf[cell_id] * grav[i];
        }
      }
      else {
        cs_real_t ro0 = cs_glob_fluid_properties->ro0;
        for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++)
        {
          for (cs_lnum_t i = 0; i < 3; i++)
            grad_pr[cell_id][i] += ro0 * grav[i];
        }
      }
    }

    cs_field_t *f_grad_pr = cs_field_by_name_try("lagr_pressure_gradient");

    /* Store pressure gradient for postprocessing purpose */
    if (f_grad_pr != nullptr) {
      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++)
          f_grad_pr->val[3*cell_id + i] = grad_pr[cell_id][i];
      }
    }

    // TODO add user momentum source terms to grad_pr

    /* Compute viscous terms
     * ===================== */
    if (cs_glob_lagr_model->viscous_terms) {

      cs_field_t *f_vel = CS_F_(vel);
      cs_equation_param_t eqp_loc;
      /* local copy of equation parameters for velocity */
      cs_field_get_key_struct(f_vel, cs_field_key_id("var_cal_opt"), &eqp_loc);
      const cs_mesh_t  *m = cs_glob_mesh;
      cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;

      cs_real_3_t *cvar_vel = (cs_real_3_t *)f_vel->val;
      const cs_real_3_t *cvar_vela = (const cs_real_3_t *)f_vel->val_pre;

      const cs_lnum_t n_b_faces = m->n_b_faces;
      const cs_lnum_t n_i_faces = m->n_i_faces;

      cs_real_t *i_visc, *b_visc;
      CS_MALLOC(i_visc, n_i_faces, cs_real_t);
      CS_MALLOC(b_visc, n_b_faces, cs_real_t);

      cs_face_viscosity(m,
                        mq,
                        eqp_loc.imvisf,
                        CS_F_(mu)->val, // cell viscosity
                        i_visc,
                        b_visc);

      cs_real_t *i_massflux = nullptr;
      cs_real_t *b_massflux = nullptr;

      CS_MALLOC(i_massflux, n_i_faces, cs_real_t);
      CS_MALLOC(b_massflux, n_b_faces, cs_real_t);
      cs_array_real_fill_zero(n_i_faces, i_massflux);
      cs_array_real_fill_zero(n_b_faces, b_massflux);

      cs_velocity_pressure_model_t *vp_model = cs_get_glob_velocity_pressure_model();

      cs_real_3_t *_div_mu_gradvel = nullptr;
      cs_real_3_t *div_mu_gradvel = nullptr;
      cs_field_t *f_visc_forces
        = cs_field_by_name_try("algo:viscous_shear_divergence");

      if (f_visc_forces != nullptr)
        div_mu_gradvel = (cs_real_3_t *)f_visc_forces->val;
      else {
        CS_MALLOC(_div_mu_gradvel,m->n_cells_with_ghosts, cs_real_3_t);
        div_mu_gradvel = _div_mu_gradvel;
      }

      /* Do not consider convective terms */
      eqp_loc.iconv = 0;
      cs_real_t *i_secvis= nullptr;
      cs_real_t *b_secvis= nullptr;

      //TODO: compute it
      if (vp_model->ivisse == 1) {
        CS_MALLOC(i_secvis, n_i_faces, cs_real_t);
        CS_MALLOC(b_secvis, n_b_faces, cs_real_t);
      }

      cs_array_real_fill_zero(3*n_cells_with_ghosts, (cs_real_t *)div_mu_gradvel);

      /* Compute - div(mu_gradu) */
      cs_balance_vector(cs_glob_time_step_options->idtvar,
                        f_vel->id,
                        0,//imasac,
                        1,//inc,
                        0, //TODO vp_model->ivisse,
                        &eqp_loc,
                        cvar_vel,
                        cvar_vela,
                        f_vel->bc_coeffs,
                        nullptr, // bc_coeffs_solve
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        i_secvis,
                        b_secvis,
                        nullptr,
                        nullptr,
                        nullptr,
                        0,
                        nullptr,
                        nullptr,
                        nullptr,
                        div_mu_gradvel);

      const cs_real_t  *cell_f_vol = mq->cell_vol;

      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
        cs_real_t dvol = 1. / cell_f_vol[cell_id];
        for (cs_lnum_t i = 0; i < 3; i++) {
          /* mind that div_mu_gradvel contains -div(mu gradvel) */
          div_mu_gradvel[cell_id][i] *= -dvol;
          grad_pr[cell_id][i] -= div_mu_gradvel[cell_id][i];
        }
      }

      CS_FREE(i_massflux);
      CS_FREE(b_massflux);
      CS_FREE(_div_mu_gradvel);
      CS_FREE(i_visc);
      CS_FREE(b_visc);
      CS_FREE(i_secvis);
      CS_FREE(b_secvis);
    }

    /* Compute velocity gradient
       ========================= */

    if (turb_disp_model || cs_glob_lagr_model->shape > 0
        || cs_glob_lagr_time_scheme->interpol_field > 0) {
      cs_field_gradient_vector(extra_i[phase_id].vel,
                               0,
                               1, /* inc */
                               grad_vel);
    }
  }

  /* Compute temperature gradient
     ============================ */

  if (   cs_glob_lagr_model->physical_model != CS_LAGR_PHYS_OFF
      && extra->temperature != nullptr
      && cs_glob_lagr_time_scheme->interpol_field > 0
      && phase_id == 0)
    cs_field_gradient_scalar(extra->temperature,
                             0,
                             1, /* inc */
                             grad_tempf);

  /* Compute Lagrangian time gradient
     ================================ */

  if (cs_glob_lagr_model->idistu > 0) {

    cs_real_t cl     = 1.0 / (0.5 + 0.75 * cs_turb_crij_c0);

    cs_real_t *energi = extra_i[phase_id].cvar_k->val;
    cs_real_t *dissip = extra_i[phase_id].cvar_ep->val;

    if (extra_i[phase_id].itytur == 3) {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        energi[cell_id] = 0.5 * (  extra_i[phase_id].cvar_rij->val[6*cell_id]
                                 + extra_i[phase_id].cvar_rij->val[6*cell_id+1]
                                 + extra_i[phase_id].cvar_rij->val[6*cell_id+2]);

    }
    else if (extra_i[phase_id].iturb == CS_TURB_K_OMEGA) {

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        dissip[cell_id] = extra_i[phase_id].cmu * energi[cell_id]
                                     * extra_i[phase_id].cvar_omg->val[cell_id];

    }
    else if (extra_i[phase_id].itytur != 2
          && extra_i[phase_id].itytur != 4
          && extra_i[phase_id].itytur != 5) {

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

    /* Initialize cell Lagrangian time
     *
     * Note: the extended time scheme should
     * compute grad(Tl*_i) = grad(Tl)/b_i - Tl grad(b_i)/b_i^2
     *
     * In the following, only the first term is kept
     * to avoid computing grad(b_i) where no particles are present.
     *
     * The other possibility would be to compute b_i everywhere
     * (taking <u_pi> from the statistic "correctly" initialized)
     */

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      if (dissip[cell_id] > 0.0 && energi[cell_id] > 0.0) {

        cs_real_t tl  = cl * energi[cell_id] / dissip[cell_id];
        tl  = cs::max(tl, cs_math_epzero);

        lagr_time->val[cell_id] = tl;

      }
      else
        lagr_time->val[cell_id] = cs_math_epzero;

    }

    if (grad_lagr_time != nullptr)
      cs_field_gradient_scalar(lagr_time,
                               0, /* use_previous_t */
                               1, /* inc: not an increment */
                               grad_lagr_time);

    if (turb_disp_model)
      compute_anisotropic_prop(iprev,
                               phase_id,
                               anisotropic_lagr_time,
                               anisotropic_bx,
                               grad_lagr_time_r_et,
                               grad_lagr_time);


    else if (grad_lagr_time_r_et != nullptr) {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
        for (int i = 0; i < 3; i++)
          grad_lagr_time_r_et[cell_id][i] = grad_lagr_time[cell_id][i];
      }
    }
  }
  else { // idistu == 0
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      lagr_time->val[cell_id] = cs_math_epzero;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
