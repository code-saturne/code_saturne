/*============================================================================
 * Methods for auxiliary mean fluid quantities (Lagrangian time and gradients)
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

/*============================================================================
 * Functions dealing with auxiliary mean fluid quantities
 *============================================================================*/

#include "cs_defs.h"

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

#include "bft_mem.h"
#include "cs_array.h"
#include "cs_balance.h"
#include "cs_face_viscosity.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_stat.h"
#include "cs_mesh.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_turbulence_model.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_aux_mean_fluid_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond (DOXYGEN_SHOULD_SKIP_THIS) */

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute auxilary mean fluid quantities.
 *
 *  - Lagragian time
 *  - gradient of total pressure
 *  - velocity gradient
 *  - Lagragian time gradient
 *
 * \param[out]  lagr_time              Lagragian time scale
 * \param[out]  grad_pr                pressure gradient
 * \param[out]  grad_vel               velocity gradient
 * \param[out]  grad_lagr_time         Lagrangian time gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_aux_mean_fluid_quantities(cs_field_t    *lagr_time,
                                  cs_real_3_t   *grad_pr,
                                  cs_real_33_t  *grad_vel,
                                  cs_real_3_t   *grad_lagr_time)
{
  cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  const cs_real_t *grav  = cs_glob_physical_constants->gravity;

  const cs_turb_model_t  *turb_model = cs_get_glob_turb_model();

  bool turb_disp_model = false;
  if (   cs_glob_lagr_model->modcpl > 0
      && cs_glob_time_step->nt_cur > cs_glob_lagr_model->modcpl
      && cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->idstnt)
    turb_disp_model = true;

  /* Use pressure gradient of NEPTUNE_CFD if needed */
  if (cs_field_by_name_try("velocity_1") != NULL) {
    cs_real_t *cpro_pgradlagr = cs_field_by_name("lagr_pressure_gradient")->val;

    for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++)
      for (cs_lnum_t id = 0; id < 3; id++)
        grad_pr[iel][id] = cpro_pgradlagr[3*iel + id];

    if (turb_disp_model || cs_glob_lagr_model->shape > 0) {

      cs_real_33_t *cpro_vgradlagr
        = (cs_real_33_t *)(cs_field_by_name("lagr_velocity_gradient")->val);

      if (cpro_vgradlagr != NULL && grad_vel != NULL) {
        for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
          for (cs_lnum_t i = 0; i < 3; i++) {
            for (cs_lnum_t j = 0; j < 3; j++)
              grad_vel[c_id][i][j] = cpro_vgradlagr[c_id][i][j];
          }
        }
      }
    }
  }

  else {

    cs_real_t *wpres = NULL;

    /* Hydrostatic pressure algorithm? */
    int hyd_p_flag = cs_glob_velocity_pressure_param->iphydr;

    cs_real_3_t *f_ext = NULL;
    if (hyd_p_flag == 1)
      f_ext = (cs_real_3_t *)(cs_field_by_name("volume_forces")->val);

    cs_real_t *solved_pres = extra->pressure->val;

    /* retrieve 2/3 rho^{n} k^{n} from solved pressure field for EVM models */
    assert(turb_model != NULL);
    if (turb_model->order <= CS_TURB_FIRST_ORDER
        && cs_glob_turb_rans_model->igrhok == 0) {
      BFT_MALLOC(wpres, n_cells_with_ghosts, cs_real_t);

      int time_id = (extra->cvar_k->n_time_vals > 1) ? 1 : 0;
      const cs_real_t *cvar_k = extra->cvar_k->vals[time_id];
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        wpres[c_id] =  solved_pres[c_id]
                     -  2./3. * extra->cromf->val[c_id]
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
      cs_field_get_equation_param_const(extra->pressure);

    cs_gradient_type_by_imrgra(eqp->imrgra,
                               &gradient_type,
                               &halo_type);

    cs_real_t *weight = NULL;
    cs_internal_coupling_t  *cpl = NULL;
    int w_stride = 1;

    if (eqp->iwgrec == 1) {
      /* Weighted gradient coefficients */
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(extra->pressure, key_id);
      if (diff_id > -1) {
        cs_field_t *weight_f = cs_field_by_id(diff_id);
        weight = weight_f->val;
        w_stride = weight_f->dim;
      }
      /* Internal coupling structure */
      key_id = cs_field_key_id_try("coupling_entity");
      if (key_id > -1) {
        int coupl_id = cs_field_get_key_int(extra->pressure, key_id);
        if (coupl_id > -1)
          cpl = cs_internal_coupling_by_id(coupl_id);
      }
    } else if (eqp->iwgrec == 0) {
      if (eqp->idiff > 0) {
        int key_id = cs_field_key_id_try("coupling_entity");
        if (key_id > -1) {
          int coupl_id = cs_field_get_key_int(extra->pressure, key_id);
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
                       extra->pressure->bc_coeffs,
                       wpres,
                       weight,
                       cpl,
                       grad_pr);

    if (wpres != solved_pres)
      BFT_FREE(wpres);

    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++)
          grad_pr[cell_id][i] += ro0 * grav[i];
      }
    }

    cs_field_t *f_grad_pr = cs_field_by_name_try("lagr_pressure_gradient");

    /* Store pressure gradient for postprocessing purpose */
    if (f_grad_pr != NULL) {
      for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
        for (cs_lnum_t i = 0; i < 3; i++)
          f_grad_pr->val[3*cell_id + i] = grad_pr[cell_id][i];
      }
    }

    //TODO add user momentum source terms to grad_pr

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
      BFT_MALLOC(i_visc, n_i_faces, cs_real_t);
      BFT_MALLOC(b_visc, n_b_faces, cs_real_t);

      cs_face_viscosity(m,
                        mq,
                        eqp_loc.imvisf,
                        CS_F_(mu)->val, // cell viscosity
                        i_visc,
                        b_visc);

      cs_real_t *i_massflux = NULL;
      cs_real_t *b_massflux = NULL;

      BFT_MALLOC(i_massflux, n_i_faces, cs_real_t);
      BFT_MALLOC(b_massflux, n_b_faces, cs_real_t);
      cs_array_real_fill_zero(n_i_faces, i_massflux);
      cs_array_real_fill_zero(n_b_faces, b_massflux);

      cs_velocity_pressure_model_t *vp_model = cs_get_glob_velocity_pressure_model();

      cs_real_3_t *_div_mu_gradvel = NULL;
      cs_real_3_t *div_mu_gradvel = NULL;
      cs_field_t *f_visc_forces
        = cs_field_by_name_try("algo:viscous_shear_divergence");

      if (f_visc_forces != NULL)
        div_mu_gradvel = (cs_real_3_t *)f_visc_forces->val;
      else {
        BFT_MALLOC(_div_mu_gradvel,m->n_cells_with_ghosts, cs_real_3_t);
        div_mu_gradvel = _div_mu_gradvel;
      }

      /* Do not consider convective terms */
      eqp_loc.iconv = 0;
      cs_real_t *i_secvis= NULL;
      cs_real_t *b_secvis= NULL;

      //TODO: compute it
      if (vp_model->ivisse == 1) {
        BFT_MALLOC(i_secvis, n_i_faces, cs_real_t);
        BFT_MALLOC(b_secvis, n_b_faces, cs_real_t);
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
                        i_massflux,
                        b_massflux,
                        i_visc,
                        b_visc,
                        i_secvis,
                        b_secvis,
                        NULL,
                        NULL,
                        NULL,
                        0,
                        NULL,
                        NULL,
                        NULL,
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

      BFT_FREE(i_massflux);
      BFT_FREE(b_massflux);
      BFT_FREE(_div_mu_gradvel);
      BFT_FREE(i_visc);
      BFT_FREE(b_visc);
      BFT_FREE(i_secvis);
      BFT_FREE(b_secvis);
    }

    /* Compute velocity gradient
       ========================= */

    if (turb_disp_model || cs_glob_lagr_model->shape > 0
        || cs_glob_lagr_time_scheme->interpol_field) {
      cs_field_gradient_vector(extra->vel,
                               0,
                               1, /* inc */
                               grad_vel);
    }
  }

  /* Compute Lagrangian time gradient
     ================================ */
  if(cs_glob_lagr_model->idistu > 0) {

    cs_real_t cl     = 1.0 / (0.5 + 0.75 * cs_turb_crij_c0);

    cs_real_t *energi = extra->cvar_k->val;
    cs_real_t *dissip = extra->cvar_ep->val;

    if (extra->itytur == 3) {
      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        energi[cell_id] = 0.5 * (  extra->cvar_rij->val[6*cell_id]
                                 + extra->cvar_rij->val[6*cell_id+1]
                                 + extra->cvar_rij->val[6*cell_id+2]);

    }
    else if (extra->iturb == CS_TURB_K_OMEGA) {

      for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
        dissip[cell_id] = extra->cmu * energi[cell_id]
                                     * extra->cvar_omg->val[cell_id];

    }
    else if (extra->itytur != 2 && extra->itytur != 4 && extra->itytur != 5) {

      bft_error
        (__FILE__, __LINE__, 0,
         _("Lagrangian turbulent dispersion is not compatible with\n"
           "the selected turbulence model.\n"
           "\n"
           "Turbulent dispersion is taken into account with idistu = %d\n"
           " Activated turbulence model is %d, when only k-eps, LES, Rij-eps,\n"
           " V2f or k-omega are handled."),
         (int)cs_glob_lagr_model->idistu,
         (int)extra->iturb);

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
     *
     */

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      if (dissip[cell_id] > 0.0 && energi[cell_id] > 0.0) {

        cs_real_t tl  = cl * energi[cell_id] / dissip[cell_id];
        tl  = CS_MAX(tl, cs_math_epzero);

        lagr_time->val[cell_id] = tl;

      }
      else
        lagr_time->val[cell_id] = cs_math_epzero;

    }

    if (grad_lagr_time != NULL)
      cs_field_gradient_scalar(lagr_time,
                               0, /* use_previous_t */
                               1, /* inc: not an increment */
                               grad_lagr_time);

  }
  else { // idistu == 0
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      lagr_time->val[cell_id] = cs_math_epzero;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
