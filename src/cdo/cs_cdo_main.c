/*============================================================================
 * Functions for solving equations with CDO discretizations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"

#include "cs_ale.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_boundary_zone.h"
#include "cs_cdo_assembly.h"
#include "cs_cdo_toolbox.h"
#include "cs_cdo_system.h"
#include "cs_control.h"
#include "cs_defs.h"
#include "cs_domain.h"
#include "cs_domain_op.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_equation_system.h"
#include "cs_gwf.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_maxwell.h"
#include "cs_navsto_system.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_pressure_correction.h"
#include "cs_prototypes.h"
#include "cs_solid_selection.h"
#include "cs_solidification.h"
#include "cs_thermal_system.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_volume_zone.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_main.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local definitions
 *============================================================================*/

static int  _cdo_ts_id = -1;
static long long  cs_cdo_setup_time = 0.;

static cs_property_t  *cs_dt_pty = NULL;

static bool _initialized_setup = false;
static bool _initialized_structures = false;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cdo_solve_steady_state_domain(void);

void
cs_f_cdo_solve_unsteady_state_domain(void);

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if this is the last iteration
 *
 * \param[in]  ts    pointer to a cs_time_step_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_is_last_iter(const cs_time_step_t   *ts)
{
  if (ts->t_max > 0) /* t_max has been set */
    if (ts->t_cur + ts->dt[0] >= ts->t_max)
      return true;

  if (ts->nt_max > 0) /* nt_max has been set */
    if (ts->nt_cur + 1 >= ts->nt_max)
      return true;

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Check if one needs to solve the steady-state thermal equation
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_needs_solving_steady_state_thermal(void)
{
  if (cs_thermal_system_is_activated() == false)
    return false;

  cs_equation_param_t  *thm_eqp =
    cs_equation_get_param(cs_equation_by_name(CS_THERMAL_EQNAME));

  if (cs_equation_param_has_time(thm_eqp))
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Check if one needs to solve the thermal equation
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

static inline bool
_needs_solving_thermal(void)
{
  if (cs_thermal_system_is_activated() == false)
    return false;

  cs_flag_t  thm_model = cs_thermal_system_get_model();

  /* Is there an advection term arising from the Navier--Stokes ? */

  if (thm_model & CS_THERMAL_MODEL_NAVSTO_ADVECTION)
    return false; /* This is managed inside the function
                     cs_navsto_system_compute_steady_state() */
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute equations which user-defined and steady-state
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_steady_user_equations(cs_domain_t   *domain)
{
  int  n_equations = cs_equation_get_n_equations();

  for (int eq_id = 0; eq_id < n_equations; eq_id++) {

    cs_equation_t  *eq = cs_equation_by_id(eq_id);

    if (cs_equation_is_steady(eq)) {

      cs_equation_type_t  type = cs_equation_get_type(eq);

      if (type == CS_EQUATION_TYPE_USER) {

        cs_flag_t  eq_flag = cs_equation_get_flag(eq);

        if ((eq_flag & CS_EQUATION_USER_TRIGGERED) == 0) {

          if (cs_equation_uses_new_mechanism(eq))
            cs_equation_solve_steady_state(domain->mesh, eq);

          else { /* Deprecated */

            /* Define the algebraic system */

            cs_equation_build_system(domain->mesh, eq);

            /* Solve the algebraic system */

            cs_equation_solve_deprecated(eq);

          }

        } /* Not triggered by user */

      } /* User-defined equation */

    } /* Steady-state equation */

  } /* Loop on equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute user-defined equation which are time-dependent
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 * \param[in]       nt_cur     current number of iteration done
 */
/*----------------------------------------------------------------------------*/

static void
_compute_unsteady_user_equations(cs_domain_t   *domain,
                                 int            nt_cur)
{
  const int  n_equations = cs_equation_get_n_equations();

  if (nt_cur >= 0) {

    for (int eq_id = 0; eq_id < n_equations; eq_id++) {

      cs_equation_t  *eq = cs_equation_by_id(eq_id);

      if (!cs_equation_is_steady(eq)) {

        cs_equation_type_t  type = cs_equation_get_type(eq);

        if (type == CS_EQUATION_TYPE_USER) {

          cs_flag_t  eq_flag = cs_equation_get_flag(eq);

          if ((eq_flag & CS_EQUATION_USER_TRIGGERED) == 0) {

            if (cs_equation_uses_new_mechanism(eq))

              /* By default, a current to previous operation is
                 performed */

              cs_equation_solve(true, domain->mesh, eq);

            else { /* Deprecated */

              /* Define the algebraic system */

              cs_equation_build_system(domain->mesh, eq);

              /* Solve domain */

              cs_equation_solve_deprecated(eq);

            }

          } /* Not triggered by user */

        } /* User-defined equation */

      } /* Unsteady equations */

    } /* Loop on equations */

  } /* nt_cur > 0 */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve only steady-state equations
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_solve_steady_state_domain(cs_domain_t  *domain)
{
  if (domain->cdo_context->mode == CS_DOMAIN_CDO_MODE_ONLY) {

    /* Otherwise log is called from the FORTRAN part */

    if (!cs_equation_needs_steady_state_solve()) {
      cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
      cs_log_printf(CS_LOG_DEFAULT,
                    "# Iter: 0 >> Initial state");
      cs_log_printf(CS_LOG_DEFAULT, "\n%s\n", cs_sep_h1);

      /* Extra operations and post-processing of the computed solutions */

      cs_post_time_step_begin(domain->time_step);

      cs_domain_post(domain);

      cs_post_time_step_end();

      return;
    }
  }

  bool  do_output = cs_domain_needs_log(domain, false);

  /* Output information */

  if (domain->only_steady) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
    cs_log_printf(CS_LOG_DEFAULT, "#      Solve steady-state problem(s)\n");
    cs_log_printf(CS_LOG_DEFAULT, "%s", cs_sep_h1);
  }
  else if (do_output) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
    cs_log_printf(CS_LOG_DEFAULT,
                  "# Iter: 0 >> Solve only requested steady-state equations");
    cs_log_printf(CS_LOG_DEFAULT, "\n%s\n", cs_sep_h1);
  }

  /* User-defined update/settings of physical properties */

  cs_user_physical_properties(domain);

  /* Predefined equation for the computation of the wall distance */

  if (cs_walldistance_is_activated())
    cs_walldistance_compute(domain->mesh,
                            domain->time_step,
                            domain->connect,
                            domain->cdo_quantities);

  /* If the problem is globally unsteady, only steady-state equations are
     solved */

  /* 1. Thermal module */

  if (_needs_solving_steady_state_thermal())
    cs_thermal_system_compute_steady_state(domain->mesh,
                                           domain->connect,
                                           domain->cdo_quantities,
                                           domain->time_step);

  /* 2. Groundwater flow module */

  if (cs_gwf_is_activated())
    cs_gwf_compute_steady_state(domain->mesh,
                                domain->time_step,
                                domain->connect,
                                domain->cdo_quantities);

  /* 3. Maxwell module */

  if (cs_maxwell_is_activated())
    cs_maxwell_compute_steady_state(domain->mesh,
                                    domain->time_step,
                                    domain->connect,
                                    domain->cdo_quantities);

  /* 4. Navier-Stokes module */

  if (cs_navsto_system_is_activated())
    cs_navsto_system_compute_steady_state(domain->mesh,
                                          domain->connect,
                                          domain->cdo_quantities,
                                          domain->time_step);

  /* User-defined equations */

  _compute_steady_user_equations(domain);

  /* Extra operations and post-processing of the computed solutions */

  cs_post_time_step_begin(domain->time_step);

  if (domain->only_steady) { /* Force writer activation and mesh output */

    cs_post_activate_writer(CS_POST_WRITER_ALL_ASSOCIATED, true);
    cs_post_write_meshes(domain->time_step);

  }

  cs_domain_post(domain);

  cs_post_time_step_end();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]  ts      pointer to a cs_time_step_t structure
 * \param[in, out]  ts_opt  pointer to a cs_time_step_options_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_define_current_time_step(cs_time_step_t           *ts,
                          cs_time_step_options_t   *ts_opt)
{
  assert(cs_dt_pty != NULL);
  assert(cs_property_is_uniform(cs_dt_pty));

  const double  t_cur = ts->t_cur;

  ts->dt[2] = ts->dt[1];
  ts->dt[1] = ts->dt[0];
  ts->dt[0] = cs_property_get_cell_value(0, t_cur, cs_dt_pty);

  if (cs_property_is_steady(cs_dt_pty) == false) {  /* dt_cur may change */

    /* Update time_options */

    double  dtmin = CS_MIN(ts_opt->dtmin, ts->dt[0]);
    double  dtmax = CS_MAX(ts_opt->dtmax, ts->dt[0]);

    ts_opt->dtmin = dtmin;
    ts_opt->dtmax = dtmax;

    /* TODO: Check how the following value is set in FORTRAN
     * domain->time_options.dtref = 0.5*(dtmin + dtmax); */

    if (ts->dt_ref < 0) /* Should be the initial val. */
      ts->dt_ref = ts->dt[0];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_solve_domain(cs_domain_t  *domain)
{
  const cs_time_step_t  *ts = domain->time_step;
  const int  nt_cur = ts->nt_cur;

  bool  do_output = cs_domain_needs_log(domain, true);

  /* Output information */

  if (do_output) {

    const double  t_cur = ts->t_cur;
    const double  dt_cur = ts->dt[0];

    cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
    cs_log_printf(CS_LOG_DEFAULT, "# Iter: %d >>"
                  " Solve domain from time=%6.4e to %6.4e; dt=%5.3e",
                  nt_cur + 1, t_cur, t_cur + dt_cur, dt_cur);
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);

  }

  cs_domain_set_stage(domain, CS_DOMAIN_STAGE_TIME_STEP_BEGIN);

  /* User-defined update/settings of physical properties */

  cs_user_physical_properties(domain);

  /* Solve predefined systems */

  if (cs_solidification_is_activated()) {

    cs_solidification_compute(domain->mesh,
                              domain->connect,
                              domain->cdo_quantities,
                              domain->time_step);

  }
  else {

    /* 1. Thermal module */

    if (_needs_solving_thermal())
      cs_thermal_system_compute(true, /* current to previous */
                                domain->mesh,
                                domain->connect,
                                domain->cdo_quantities,
                                domain->time_step);

    /* 2. Groundwater flow module */

    if (cs_gwf_is_activated())
      cs_gwf_compute(domain->mesh,
                     domain->time_step,
                     domain->connect,
                     domain->cdo_quantities);

    /* 3. Maxwell module */

    if (cs_maxwell_is_activated())
      cs_maxwell_compute(domain->mesh,
                         domain->time_step,
                         domain->connect,
                         domain->cdo_quantities);

    /* 4. Navier-Stokes module */

    if (cs_navsto_system_is_activated())
      cs_navsto_system_compute(domain->mesh,
                               domain->connect,
                               domain->cdo_quantities,
                               domain->time_step);

  }

  cs_domain_set_stage(domain, CS_DOMAIN_STAGE_TIME_STEP_END);

  /* User-defined update/settings of physical properties */

  cs_user_physical_properties(domain);

  /* User-defined equations */

  _compute_unsteady_user_equations(domain, nt_cur);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary the setup of all major structures:
 *          cs_domain_t structure, all equations and all properties
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

static void
_log_setup(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return;

  /* Output domain settings */

  cs_domain_setup_log(domain);

  if (domain->verbosity > -1) {

    /* Advection fields */

    cs_advection_field_log_setup();

    /* Properties */

    cs_property_log_setup();

    /* Summary of the thermal module */

    cs_thermal_system_log_setup();

    /* Summary of the groundwater module */

    cs_gwf_log_setup();

    /* Summary of the Maxwell module */

    cs_maxwell_log_setup();

    /* Summary of the Navier-Stokes system */

    cs_navsto_system_log_setup();

    /* Summary of the solidification module */

    cs_solidification_log_setup();

  } /* domain->verbosity > 0 */

  /* Summary of the setup for each equation */

  cs_equation_log_setup();

  /* Summary of the setup for each system of equations */

  cs_equation_system_log_setup();

  /* SLES options:
   * ------------
   * Part common with FV schemes
   * Options attached to the cs_sles_t structure (not the cs_sles_param_t
   * structure which is only used in the CDO part)
   */

  cs_log_printf(CS_LOG_SETUP, " Additional settings for SLES\n%s", cs_sep_h1);
  cs_sles_log(CS_LOG_SETUP);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the performance monitoring for the major steps
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

static void
_log_monitoring(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return;

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\nSummary of the performance monitoring for the CDO part\n"
                "------------------------------------------------------\n\n");

  long long  connect_time = cs_cdo_connect_get_time_perfo();
  long long  quant_time = cs_cdo_quantities_get_time_perfo();

  cs_log_printf(CS_LOG_PERFORMANCE, " %-38s Connect  Quantities  Total\n", " ");
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f %9.3f %9.3f seconds\n",
                "<CDO/Setup> Runtime", connect_time*1e-9, quant_time*1e-9,
                cs_cdo_setup_time*1e-9);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f seconds\n",
                "<CDO/ExtraOperations> Runtime", domain->tcp.nsec*1e-9);

  cs_log_printf(CS_LOG_PERFORMANCE, "\n %-35s %9.3f seconds\n",
                "<CDO> Total runtime", domain->tca.nsec*1e-9);

  cs_equation_system_log_monitoring();

  cs_equation_log_monitoring();
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve only steady-state equations
 */
/*----------------------------------------------------------------------------*/

void
cs_f_cdo_solve_steady_state_domain(void)
{
  _solve_steady_state_domain(cs_glob_domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve all the equations of a computational domain for one time step
 */
/*----------------------------------------------------------------------------*/

void
cs_f_cdo_solve_unsteady_state_domain(void)
{
  _solve_domain(cs_glob_domain);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the computational domain when CDO/HHO schemes are
 *         activated and cs_user_model() has been called.
 *         At this stage of the settings, mesh quantities and adjacencies are
 *         not defined. Only the major moddeling options are set. The related
 *         equations and main properties have been added.
 *
 * \param[in, out]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_setup(cs_domain_t   *domain)
{
  if (cs_domain_get_cdo_mode(domain) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  /* Timer statistics */

  cs_timer_t t0 = cs_timer_time();
  _cdo_ts_id = cs_timer_stats_id_by_name("cdo");
  if (_cdo_ts_id < 0)
    _cdo_ts_id = cs_timer_stats_create("stages", "cdo", "cdo");

  /* Store the fact that the CDO/HHO module is activated */

  cs_domain_cdo_log(domain);

  /* A property can be called easily from everywhere:
   * cs_property_by_name("time_step")
   * This is useful when building linear system.
   */

  cs_dt_pty = cs_property_by_name("time_step");

  /* Check that predefined properties have been created */

  assert(cs_property_by_name("unity") != NULL);
  assert(cs_dt_pty != NULL);

  cs_timer_stats_start(_cdo_ts_id);

  /* Add an automatic boundary zone gathering all "wall" boundaries */

  cs_boundary_def_wall_zones(domain->boundaries);

  /* First setup stage of the cs_domain_t structure
   * - Define extra domain boundaries
   * - Setup predefined equations
   * - Create fields
   */

  cs_domain_initialize_setup(domain);

  _initialized_setup = true;

  /* Monitoring */

  cs_timer_stats_stop(_cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(domain->tca, domain->tca, time_count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build additional connectivities and quantities when CDO/HHO schemes
 *         are activated.
 *         Finalize the setup and from the settings, define the structures
 *         related to equations and modules
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  m        pointer to a cs_mesh_t struct.
 * \param[in]       mq       pointer to a cs_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_structures(cs_domain_t           *domain,
                             cs_mesh_t             *m,
                             cs_mesh_quantities_t  *mq)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: cs_domain_t structure is not allocated.\n", __func__);

  domain->mesh = m;
  domain->mesh_quantities = mq;

  if (cs_domain_get_cdo_mode(domain) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  cs_timer_t  t0 = cs_timer_time();

  /* Timer statistics */

  cs_timer_stats_start(_cdo_ts_id);

  cs_domain_init_cdo_structures(domain);

  /* Last user setup stage */

  cs_domain_finalize_user_setup(domain);

  /* Assign to a cs_equation_t structure a list of functions to manage this
   * structure during the computation.
   * The set of functions chosen for each equation depends on the parameters
   * specifying the cs_equation_t structure */

  domain->only_steady = cs_equation_set_functions();

  /* Initialize and set main members before building and solving systems of
     equations (should be done after cs_equation_set_functions()) */

  cs_equation_system_set_functions();

  if (domain->only_steady)
    domain->is_last_iter = true;

  else { /* Setup the time step if not already done */

    if (cs_dt_pty == NULL)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Please check your settings.\n"
                " Unsteady computation but no current time step defined.\n",
                __func__);

    if (cs_dt_pty->n_definitions == 0)
      /* No definition available yet. Try a definition by value */
      cs_domain_automatic_time_step_settings(domain);

    if (cs_property_is_uniform(cs_dt_pty) == false)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Please check your settings.\n"
                " Unsteady computation with a non-uniform time step.\n",
                __func__);

  }

  /* Last setup stage */

  cs_domain_finalize_module_setup(domain);

  /* Initialization of the default post-processing for the computational
     domain */

  cs_domain_post_init(domain);

  /* Define builder structures for equations */

  cs_equation_define_builders(m);

  /* Define the context structure for each equation. One relies on the given
   * settings (low-level function pointers
   * Do the same thing for systems of equations and the NavSto module if needed
   */

  cs_equation_define_context_structures();
  cs_equation_system_define();

  if (cs_navsto_system_is_activated())
    cs_navsto_system_define_context(m);

  /* SLES settings:
   * -------------
   *  Define associated cs_sles_t structure to solve the linear systems.
   *  This should be done after the creation of fields and the definition of
   *  the scheme context and system helper structures.
   */

  cs_equation_set_sles();
  cs_equation_system_set_sles();

  if (cs_navsto_system_is_activated())
    cs_navsto_system_set_sles();

  /* Summary of the settings */

  _log_setup(domain);

  /* Flush log files */

  cs_log_printf_flush(CS_LOG_DEFAULT);
  cs_log_printf_flush(CS_LOG_SETUP);
  cs_log_printf_flush(CS_LOG_PERFORMANCE);

  /* Monitoring */

  cs_timer_stats_stop(_cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  _initialized_structures = true;

  CS_TIMER_COUNTER_ADD(domain->tca, domain->tca, time_count);
  cs_cdo_setup_time += domain->tca.nsec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free all structures allocated during the resolution of CDO/HHO
 *          schemes
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_finalize(cs_domain_t    *domain)
{
  if (cs_domain_get_cdo_mode(domain) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  /* Timer statistics */

  cs_timer_stats_start(_cdo_ts_id);

  /* Write a restart file if needed */

  cs_domain_write_restart(domain);

  /* Clean up for restart multiwriters */

  cs_restart_clean_multiwriters_history();

  /* Print monitoring information */

  _log_monitoring(domain);

  /* Free common structures relatated to equations */

  cs_equation_finalize_sharing(domain->cdo_context->vb_scheme_flag,
                               domain->cdo_context->vcb_scheme_flag,
                               domain->cdo_context->eb_scheme_flag,
                               domain->cdo_context->fb_scheme_flag,
                               domain->cdo_context->hho_scheme_flag);

  /* Free the memory related to equations */

  cs_equation_destroy_all();

  /* Free the memory related to systems of equations */

  cs_equation_system_destroy_all();

  /* Free memory related to advection fields */

  cs_advection_field_destroy_all();

  /* Free the memory used inside modules using CDO schemes */
  /* ----------------------------------------------------- */

  /* Solid cells */

  cs_solid_selection_free();

  /* ALE */

  cs_ale_destroy_all();

  /* Free memory related to the thermal module */

  cs_thermal_system_destroy();

  /* Free memory related to the groundwater flow module */

  cs_gwf_destroy_all();

  /* Free memory related to the Maxwell module */

  cs_maxwell_destroy_all();

  /* Navier-Stokes system */

  cs_navsto_system_destroy();

  /* CDO resolved Pressure correction coupled with FV*/

  cs_pressure_correction_cdo_destroy_all();

  /* Solidification module */

  cs_solidification_destroy_all();

  /* Free the memory related to "low-level" structures shared among schemes */
  /* ---------------------------------------------------------------------- */

  cs_cdo_assembly_finalize();

  cs_cdo_system_destroy_all();

  cs_cdo_toolbox_finalize();

  /* Set flag to OFF */

  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_OFF);

  cs_log_printf(CS_LOG_DEFAULT,
                "\n  Finalize and free CDO-related structures.\n");

  _initialized_setup = false;
  _initialized_structures = false;

  /* Free CDO structures related to geometric quantities and connectivity */

  domain->cdo_quantities = cs_cdo_quantities_free(domain->cdo_quantities);
  domain->connect = cs_cdo_connect_free(domain->connect);

  cs_timer_stats_stop(_cdo_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if CDO has been initialized.
 *
 * \param[in, out]  setup       indicator if setup has been initialized,
 *                              or NULL if not queried
 * \param[in, out]  structures  indicator if structures have been initialized,
 *                              or NULL if not queried
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_is_initialized(bool  *setup,
                      bool  *structures)
{
  if (setup != NULL)
    *setup = _initialized_setup;

  if (structures != NULL)
    *structures = _initialized_structures;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main program for running a simulation with CDO kernel
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_main(cs_domain_t   *domain)
{
  if (cs_domain_get_cdo_mode(domain) == CS_DOMAIN_CDO_MODE_OFF)
    return;
  if (cs_equation_get_n_equations() < 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n  No equation to solve. Immediate exit\n");
    return;
  }

  cs_timer_t t0 = cs_timer_time();

  /* Timer statistics */

  cs_timer_stats_start(_cdo_ts_id);

  /* Read a restart file if needed */

  cs_domain_read_restart(domain);

  /* Force the activation of writers for postprocessing */

  cs_post_activate_writer(CS_POST_WRITER_ALL_ASSOCIATED, true);

  /*  Build high-level structures and create algebraic systems
      Set the initial values of the fields and properties. */

  cs_domain_initialize_systems(domain);

  /* Remark: cs_user_initialization is called at the end of the previous
     function */

  /* Initialization for user-defined extra operations. Should be done
     after the domain initialization if one wants to overwrite the field
     initialization for instance */

  cs_user_extra_operations_initialize(cs_glob_domain);

  /* Output information */

  cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
  cs_log_printf(CS_LOG_DEFAULT, "#      Start main loop\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", cs_sep_h1);

  /* ======================================== */
  /* Compute first the steady-state equations */
  /* ======================================== */

  _solve_steady_state_domain(domain);

  cs_domain_set_stage(domain, CS_DOMAIN_STAGE_BEFORE_TIME_LOOP);

  /* User-defined update/settings of physical properties after solving the
     steady-state equations */

  cs_user_physical_properties(domain);

  /* ============== */
  /* Main time loop */
  /* ============== */

  while (cs_domain_needs_iteration(domain)) {

    /* Define the current time step */

    _define_current_time_step(domain->time_step, &(domain->time_options));

    /* Check if this is the last iteration */

    domain->is_last_iter = _is_last_iter(domain->time_step);

    /* Build and solve equations related to the computational domain */

    _solve_domain(domain);

    /* Increment time */

    cs_domain_increment_time(domain);

    /* Increment time steps */

    cs_domain_increment_time_step(domain);

    /* Extra operations and post-processing of the computed solutions */

    cs_post_time_step_begin(domain->time_step);

    cs_domain_post(domain);

    cs_post_time_step_end();

    /* Read a control file if present */

    cs_control_check_file();

    /* Add a checkpoint if needed */

    cs_domain_write_restart(domain);

    /* Clean up for restart multiwriters */

    cs_restart_clean_multiwriters_history();

    cs_timer_stats_increment_time_step();

  }

  cs_domain_set_stage(domain, CS_DOMAIN_STAGE_AFTER_TIME_LOOP);

  /* User-defined update/settings of physical properties (finalization stage) */

  cs_user_physical_properties(domain);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(domain->tca, domain->tca, time_count);

  cs_timer_stats_stop(_cdo_ts_id);
  if (cs_glob_rank_id <= 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", cs_sep_h1);
    cs_log_printf(CS_LOG_DEFAULT, "#\tExit CDO core module\n");
    cs_log_printf(CS_LOG_DEFAULT, "%s", cs_sep_h1);
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
