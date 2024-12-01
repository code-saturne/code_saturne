/*============================================================================
 * Initialize fields.
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
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_array.h"
#include "cs_array_reduce.h"
#include "cs_assert.h"
#include "cs_cf_model.h"
#include "cs_cf_thermo.h"
#include "cs_ctwr_initialize.h"
#include "cs_dispatch.h"
#include "cs_elec_model.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_gas_mix.h"
#include "cs_gui.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_ibm.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parameters_check.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_porosity_from_scan.h"
#include "cs_prototypes.h"
#include "cs_scalar_clipping.h"
#include "cs_solve_navier_stokes.h"
#include "cs_time_step.h"
#include "cs_turbulence_init.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_initialize_fields.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_ppiniv0(void);

void
cs_f_user_initialization_wrapper(cs_real_t  dt[]);

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_initialize_fields.cpp

  \brief Various field initializations.
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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computed variable and property initialization.
 *
 * The time step, the indicator of wall distance computation are also
 * initialized just before reading a restart file or use the user
 * initializations.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_0(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const int n_fields = cs_field_n_fields();

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kdflim = cs_field_key_id("diffusion_limiter_id");
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");

  /* Initialize temperature to reference temperature if present */
  if (CS_F_(t) != nullptr) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->t0,
                             CS_F_(t)->val);
  }

  /* Initialize boundary temperature to "marker" if present */
  if (CS_F_(t_b) != nullptr) {
    cs_array_real_set_scalar(m->n_b_faces,
                             -cs_math_big_r,
                             CS_F_(t_b)->val);
  }

  /* Time step
     --------- */

  /* dt might be used on the halo cells during the ALE initialization
     otherwise dt is synchronized in the pressure correction step. */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           cs_glob_time_step->dt_ref,
                           CS_F_(dt)->val);

  /* Initialize physical properties
     ------------------------------ */

  /* Density */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           fluid_props->ro0,
                           CS_F_(rho)->val);
  cs_array_real_set_scalar(m->n_b_faces,
                           fluid_props->ro0,
                           CS_F_(rho_b)->val);

  /* Boussinesq (set to impossible value at this stage) */
  if (cs_glob_velocity_pressure_model->idilat == 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             -1.,
                             cs_field_by_name("thermal_expansion")->val);
  }

  /* Molecular viscosity */
  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           fluid_props->viscl0,
                           CS_F_(mu)->val);

  /* Specific heat */
  if (fluid_props->icp >= 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->cp0,
                             CS_F_(cp)->val);
  }

  /* The total pressure will be initialized to P0 + rho.g.r in inivar
     if the user has not done a prior initialization;
     so mark to "uninit" value here */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             - cs_math_infinite_r,
                             cs_field_by_name("total_pressure")->val);
  }

  /* Initialization of mix_mol_mas with default values (air)
     (used in cs_cf_thermo_default_init) */
  if (cs_glob_physical_model_flag[CS_GAS_MIX] >= 0) {
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->xmasmr,
                             cs_field_by_name("mix_mol_mas")->val);
  }

  /* Default initialisations for the compressible model */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
    /* In compressible, for now, the temperature is not solved but is a field of
       type variable anyway. The reference value has to be taken into account. */
    cs_array_real_set_scalar(m->n_cells_with_ghosts,
                             fluid_props->t0,
                             CS_F_(t)->val);

    /* Default isochoric specific heat (cv0), total energy and density */
    cs_cf_thermo_default_init();

    /* Default diffusivity for total energy */
    double visls_0 = cs_field_get_key_double(CS_F_(t), kvisl0);
    visls_0 /= fluid_props->cv0;
    cs_field_set_key_double(CS_F_(e_tot), kvisl0, visls_0);
  }

  /* Scalars diffusivity */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (!(f->type & CS_FIELD_VARIABLE) || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = cs_field_get_key_int(f, keysca) - 1;
    if (scalar_id < 0)
      continue;

    int ifcvsl = cs_field_get_key_int(f, kivisl);
    double visls_0 = cs_field_get_key_double(f, kvisl0);
    /* Diffusivity at current time step (and previous one if second order) */
    if (ifcvsl >= 0) {
      cs_field_t *f_vsl = cs_field_by_id(ifcvsl);
      for (int i = 0; i < f_vsl->n_time_vals; i++) {
        cs_array_real_set_scalar(m->n_cells_with_ghosts,
                                 visls_0,
                                 f_vsl->vals[i]);
      }
    }
  }

  /* Mesh viscosity for ALE */

  if (cs_glob_ale != CS_ALE_NONE) {
    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(CS_F_(mesh_u));
    const int idftnp = eqp->idften;

    if (idftnp & CS_ANISOTROPIC_LEFT_DIFFUSION) {
      cs_real_t ref_val[6] = {1., 1., 1., 0., 0., 0.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              6,
                              ref_val,
                              cs_field_by_name("mesh_viscosity")->val);
    }
    else if (idftnp & CS_ISOTROPIC_DIFFUSION) {
      cs_array_real_set_scalar(m->n_cells_with_ghosts,
                               1,
                               cs_field_by_name("mesh_viscosity")->val);
    }

    cs_gui_mesh_viscosity();
  }

  /* Porosity */

  if (cs_glob_porous_model > 0) {
    cs_real_t *porosity = cs_field_by_name("porosity")->val;

    if (cs_glob_porosity_from_scan_opt->compute_porosity_from_scan)
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 0., porosity);

    else if (cs_glob_porosity_ibm_opt->porosity_mode > 0) {
      cs_field_t *f_ifp = cs_field_by_name_try("i_face_porosity");
      cs_field_t *f_bfp = cs_field_by_name_try("b_face_porosity");
      if (f_ifp != nullptr)
        cs_array_real_set_scalar(m->n_i_faces, 1., f_ifp->val);
      if (f_bfp != nullptr)
        cs_array_real_set_scalar(m->n_b_faces, 1., f_bfp->val);
    }

    else
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 1., porosity);

    /* Tensorial porosity */
    if (cs_glob_porous_model == 2) {
      cs_real_t ref_val[6] = {1., 1., 1., 0., 0., 0.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              6,
                              ref_val,
                              cs_field_by_name("tensorial_porosity")->val);
    }
  }

  /* Reference value for P*
     ---------------------- */

  cs_array_real_set_scalar(m->n_cells_with_ghosts,
                           cs_glob_fluid_properties->pred0,
                           CS_F_(p)->val);

  /* VoF
     --- */

  /* Note that the initialization here should not be strictly
     necessary, as cs_scalar_clipping should handle it, but
     doing so here will avoid activating the clipping
     counters; maybe a better manner of doing things would
     be for the counters to be deactivated at this step,
     which comes before user settings */

  if (cs_glob_vof_parameters->vof_model > 0) {
    const int kscmin = cs_field_key_id_try("min_scalar_clipping");

    cs_field_t *f_vf = cs_field_by_name("void_fraction");
    cs_real_t clvfmn = cs_field_get_key_double(f_vf, kscmin);

    cs_array_real_set_scalar(m->n_cells_with_ghosts, clvfmn, f_vf->val);
  }

  /* Turbulence
     ---------- */

  cs_turbulence_init_by_ref_quantities();

  /* Clip scalar variables (except k-epsilon, see above)
     --------------------------------------------------- */

  /* Clipping of non-variance scalars */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;
    if (cs_field_get_key_int(f, kscavr) >= 0)  /* is a variance */
      continue;

    cs_scalar_clipping(f);
  }

  /* Clipping of variances  */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (   (! (f->type & CS_FIELD_VARIABLE))
        || (f->type & CS_FIELD_CDO)) /* CDO variables handled elsewhere */
      continue;
    if (cs_field_get_key_int(f, keysca) < 0)
      continue;

    if (cs_field_get_key_int(f, kscavr) >= 0)
      cs_scalar_clipping(f);
  }

  /* ALE
     --- */

  if (cs_glob_ale > CS_ALE_NONE) {
    cs_field_t *f_coord0 = cs_field_by_name("vtx_coord0");
    cs_array_real_copy(m->n_vertices*3, m->vtx_coord, f_coord0->val);
  }

  /* Specific field initializations
     ------------------------------ */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    /* Diffusion limiters;
       Note that allowing a default value initialization for 1D-fields
       would allow handling this upon field definition */

    const int dl_id = cs_field_get_key_int(f, kdflim);
    if (dl_id > -1) {
      cs_field_t *f_dl = cs_field_by_id(dl_id);
      cs_array_real_set_scalar(m->n_cells_with_ghosts, 1.0, f_dl->val);
    }

  }

  /* Previous value initializations
     ------------------------------

     This is not restricted to variables
     (for example Cp, Cv, and mass fluxes may be handled here also,
     depending on the time stepping options). */

  for (int i = 0; i < n_fields; i++) {
    cs_field_t *f = cs_field_by_id(i);

    if (f->type & CS_FIELD_CDO) /* CDO variables handled elsewhere */
      continue;

    for (int n_p = f->n_time_vals; n_p > 1; n_p--)
      cs_field_current_to_previous(f);
  }
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variable, time step, and wall distance fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_initialize_fields_stage_1(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;

  const cs_lnum_t n_cells = m->n_cells;
  const int n_fields = cs_field_n_fields();
  const int *pm_flag = cs_glob_physical_model_flag;

  const int keysca = cs_field_key_id("scalar_id");
  const int kscavr = cs_field_key_id("first_moment_id");
  const int kscmin = cs_field_key_id("min_scalar_clipping");
  const int kscmax = cs_field_key_id("max_scalar_clipping");
  const int kclvfl = cs_field_key_id("variance_clipping");
  const int kwgrec = cs_field_key_id_try("gradient_weighting_id");

  const bool have_restart_aux
    = (    cs_restart_present()
        && cs_glob_restart_auxiliary->read_auxiliary > 0);

  /* Initialize gradient weighting fields
     ------------------------------------ */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (   !(f->type & CS_FIELD_VARIABLE)
        || f->type & CS_FIELD_CDO
        || f->dim != 1)
      continue;

    const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
    if (eqp->idiff < 1 || eqp->iwgrec != 1)
      continue;

    const cs_field_t *f_wg = cs_field_by_id(cs_field_get_key_int(f, kwgrec));
    if (f_wg->dim == 6) {
      const cs_real_t  t_val[6] = {1., 1., 1., 0., 0., 0.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              6,
                              t_val,
                              f_wg->val);
    }
    else if (f_wg->dim == 1) {
      const cs_real_t  s_val[1] = {1.};
      cs_array_real_set_value(m->n_cells_with_ghosts,
                              1,
                              s_val,
                              f_wg->val);
    }
  }

  /* Pre-initialization for specific physical models
     ----------------------------------------------- */

  if (pm_flag[CS_PHYSICAL_MODEL_FLAG] > 0) {
    cs_f_ppiniv0();

    if (pm_flag[CS_COOLING_TOWERS] >= 0)
      cs_ctwr_fields_init0();

    if (pm_flag[CS_JOULE_EFFECT] >= 1 || pm_flag[CS_ELECTRIC_ARCS] >= 1)
      cs_elec_fields_initialize(m);
  }

  /* GUI and user initialization
     --------------------------- */

  cs_gui_initial_conditions();

  {
    cs_real_t *dt_val = nullptr;
    if (CS_F_(dt) != nullptr)
      dt_val = CS_F_(dt)->val;
    cs_f_user_initialization_wrapper(dt_val);
  }
  cs_user_initialization(cs_glob_domain);

  /* Second stage of initialization for specific physical models
     -----------------------------------------------------------
     (after the user function call) */

  // Cooling towers
  if (pm_flag[CS_COOLING_TOWERS] >= 0)
    cs_ctwr_fields_init1();

  // Gas mixture modelling in presence of noncondensable gases and
  // condensable gas as stream.
  if (pm_flag[CS_GAS_MIX] >= 0)
    cs_gas_mix_initialization();

  // Compressible
  // Has to be called AFTER the gas mix initialization because the
  // mixture composition is taken into account in the thermodynamic
  // law, if gas mix specific physics is enabled.
  if (pm_flag[CS_COMPRESSIBLE] >= 0)
    cs_cf_initialize();

  /* VoF model
   ----------- */

  if (cs_glob_vof_parameters->vof_model > 0) {
    cs_vof_compute_linear_rho_mu(m);
    // density is stored at the two previous time steps
    for (int t_i = 0; t_i < 2; t_i++) {
      cs_field_current_to_previous(CS_F_(rho));
      cs_field_current_to_previous(CS_F_(rho_b));
    }
  }

  /* Compressible model
   -------------------- */

  if (   pm_flag[CS_COMPRESSIBLE] >= 0
      && have_restart_aux == false) {

    const cs_cf_model_t *cf_model = cs_glob_cf_model;
    const int ithvar = cf_model->ithvar;

    if (   ithvar !=  60000 && ithvar != 100000
        && ithvar != 140000 && ithvar != 150000 && ithvar != 210000) {
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("in compressible thermodynamics at initialization"),
         _("Unexpected value of the indicator ithvar (%d).\n\n"
           "Two and only two independant variables among\n"
           "P, rho, T and E (except T and E) should be imposed."),
         ithvar);
    }

    cs_real_t wrk0[1], wrk1[1], wrk2[1], wrk3[3];
    cs_cf_thermo(ithvar, -1, wrk0, wrk1, wrk2, &wrk3);

  }

  /* Pressure / Total pressure initialization
   * ----------------------------------------
   *
   * Standard case:
   * If the user has initialized the total pressure Ptot, P* is initialized
   * accordingly, only if the user has speficied the reference point.
   * (all values of the total pressure have to be initialized).
   * Otherwise, the total pressure is initialized using P*,
   * Ptot = P* + P0 + rho.g.r
   *
   * In case of restart without auxiliary, Ptot is recomputed with P*.
   * (For EVM models, the shift by 2/3*rho*k is missing)
   * In case of restart with auxiliary, nothing needs to be done.
   *
   * Compressible:
   * The total pressure field does not need to be defined. The solved pressure is
   * the total pressure. */

  const cs_field_t *f_pr_tot = cs_field_by_name_try("total_pressure");
  if (   f_pr_tot != nullptr && CS_F_(p) != nullptr
      && pm_flag[CS_COMPRESSIBLE] < 0) {

    const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
    const cs_real_t *gravity = cs_glob_physical_constants->gravity;
    const cs_real_t *cpro_prtot = f_pr_tot->val;
    cs_real_t *cvar_pr = CS_F_(p)->val;

    int uprtot = 0;

    if (   fluid_props->ixyzp0 > -1
        && have_restart_aux == false) {

      uprtot = 1;
      cs_real_t valmin = HUGE_VALF, valmax = - HUGE_VALF;
      cs_array_reduce_minmax(m->n_cells,
                             cpro_prtot,
                             valmin,
                             valmax);
      if (valmin <= -0.5*cs_math_infinite_r)
        uprtot = 0;

      cs_parall_min(1, CS_INT_TYPE, &uprtot);
    }

    if (uprtot > 0) {
      const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

      /* Copy global arrays to local ones to enable lambda capture for dispatch */
      const cs_real_t ro0 = fp->ro0;
      const cs_real_t pred0_m_p0 = fp->pred0 - fp->p0;
      const cs_real_t g[3] = {gravity[0], gravity[1], gravity[2]};
      const cs_real_t xyzp0[3] = {fp->xyzp0[0], fp->xyzp0[1], fp->xyzp0[1]};

      cs_dispatch_context ctx;
      ctx.parallel_for(n_cells, [=] CS_F_HOST_DEVICE (cs_lnum_t c_id) {
        cvar_pr[c_id] =    cpro_prtot[c_id]
                         - ro0 * cs_math_3_distance_dot_product(xyzp0,
                                                                cell_cen[c_id],
                                                                g)
                         + pred0_m_p0;
      });
    }

    else if (have_restart_aux == false)
      cs_solve_navier_stokes_update_total_pressure(m, mq, fp, gravity);

  }

  /* Clip turbulent quantities (initialized by user or restart)
   * ---------------------------------------------------------- */

  if (cs_turbulence_init_clip_and_verify() > 0) {
    cs_parameters_error(CS_ABORT_DELAYED,
                        _("in variables initialization"),
                        _("Errors in turbulent quantities initialization"));
  }

  /* Clip scalars (initialized by user or restart)
   * ---------------------------------------------
   *
   * If the user has modified values in `cs_user_initialization`:
   * - If values are "correct" (i.e. within prescribed bounds), the
   *   initialization is admissible, and is clipped to ensure it is
   *   coherent relative to the code's clipping mode.
   * - If values are outside of prescribed bounds, abort computation.
   *
   * The same logic is used in case of a computation restart, to ensure
   * the same behavior between a computation where a user modifies a
   * variable in `cs_user_initialization` an a computation where no
   * modification occurs.
   * Otherwise, values have already been clipped after default initialization.
   *
   * To summarize:
   * - Abort if values are out of admissible bounds.
   * - Clip when using admissible user-modified or restarted values.
   * - Calues have already been clipped in other cases.
   */

  /* Clip scalars first, as they may be needed to clip variances */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (   !(f->type & CS_FIELD_VARIABLE)
        || f->type & CS_FIELD_CDO
        || f->dim != 1)
      continue;
    int scalar_id = cs_field_get_key_int(f, keysca) - 1;
    if (scalar_id < 0)
      continue;

    if (cs_field_get_key_int(f, kscavr) >= 0)  /* is a variance */
      continue;

    /* Non-variance scalars */

    cs_real_t scminp = cs_field_get_key_double(f, kscmin);
    cs_real_t scmaxp = cs_field_get_key_double(f, kscmax);

    if (scminp <= scmaxp) {

      cs_real_t valmin = HUGE_VALF, valmax = - HUGE_VALF;
      cs_array_reduce_minmax(m->n_cells,
                             f->val,
                             valmin,
                             valmax);

      cs_parall_min_scalars(valmin);
      cs_parall_max_scalars(valmin);

      // Check coherence for clippings of non-variance scalars.
      if (valmin >= scminp && valmax <= scmaxp)
        cs_scalar_clipping(f);

      else {
        cs_parameters_error(CS_ABORT_DELAYED,
                            _("in variables initialization"),
                            _("Scalar quantities out of bounds for \"%s\":\n"
                              "  Minimum value       = %g\n"
                              "  'min_scal_clipping' = %g\n"
                              "  Maximum value       = %g\n"
                              "  'max_scal_clipping' = %g"),
                            f->name, valmin, scminp, valmax, scmaxp);

      }
    }

  } // loop on fields

  /* Now clip variances. */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);

    if (   !(f->type & CS_FIELD_VARIABLE)
        || f->type & CS_FIELD_CDO)
      continue;
    int scalar_id = cs_field_get_key_int(f, keysca) - 1;
    if (scalar_id < 0)
      continue;

    if (cs_field_get_key_int(f, kscavr) < 0)  /* not a variance */
      continue;

    /* Variance scalars */

    cs_real_t scminp = cs_field_get_key_double(f, kscmin);
    cs_real_t scmaxp = cs_field_get_key_double(f, kscmax);

    if (scminp <= scmaxp) {
      cs_real_t valmin = HUGE_VALF, valmax = - HUGE_VALF;
      cs_array_reduce_minmax(m->n_cells,
                             f->val,
                             valmin,
                             valmax);

      cs_parall_min_scalars(valmin);
      cs_parall_max_scalars(valmin);

      int iclvfl = cs_field_get_key_int(f, kclvfl);

      // Check coherence for variance clippings.
      // For iclvfl = 1, only check positivity, otherwise it well be
      // difficult do have a correct initialization.

      if (valmin >= scminp && valmax <= scmaxp) {

        if (iclvfl == 0) {
          // We could clip if valmin >= 0, but by definition,
          // this would not add anything.
          if (valmin < 0)
            cs_parameters_error(CS_ABORT_DELAYED,
                                _("in variables initialization"),
                                _("Negative variance for \"%s\":\n"
                                  "  Minimum value = %g"),
                                f->name, valmin);
        }
        else if (iclvfl == 1) {
          // Here we clip to be coherent with the scalar's value.
          if (valmin >= 0)
            cs_scalar_clipping(f);
          else
            cs_parameters_error(CS_ABORT_DELAYED,
                                _("in variables initialization"),
                                _("Negative variance for \"%s\":\n"
                                  "  Minimum value = %g"),
                                f->name, valmin);
        }
        else if (iclvfl == 2) {
          cs_real_t vfmin = fmax(scminp, 0);
          cs_real_t vfmax = scmaxp;
          // We could clip when valmin >= vfmin and valmax <= vfmax
          // but by definition, this would add nothing.
          if (valmin < vfmin || valmax > vfmax)
            cs_parameters_error
              (CS_ABORT_DELAYED,
               _("in variables initialization"),
               _("Variance out of bounds or negative for \"%s\":\n"
                 "  Minimum value       = %g\n"
                 "  'min_scal_clipping' = %g\n"
                 "  Maximum value       = %g\n"
                 "  'max_scal_clipping' = %g\n\n"
                 "  Clipping mode       = %d"),
               f->name, valmin, scminp, valmax, scmaxp, iclvfl);
        }
      }
    }

  } // loop on fields

  cs_parameters_error_barrier();

  cs_user_extra_operations_initialize(cs_glob_domain);

  cs_log_printf(CS_LOG_DEFAULT,
                _("\n"
                  " ** VARIABLES INITIALIZATION\n"
                  "    ------------------------\n\n"));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
