/*============================================================================
 * Routines for solving equations with CDO discretizations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_control.h"
#include "cs_defs.h"
#include "cs_domain.h"
#include "cs_domain_boundary.h"
#include "cs_domain_op.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_gwf.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_navsto_system.h"
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

static int  cs_cdo_ts_id;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

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

        if (cs_equation_uses_new_mechanism(eq))
          cs_equation_solve_steady_state(domain->mesh, eq);

        else { /* Deprecated */

          double  dt_cur = 0.; /* Useless for steady-state equations */

          /* Define the algebraic system */
          cs_equation_build_system(domain->mesh,
                                   domain->time_step,
                                   dt_cur,
                                   eq);

          /* Solve the algebraic system */
          cs_equation_solve_deprecated(eq);

        }

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

  if (nt_cur > 0) {

    for (int eq_id = 0; eq_id < n_equations; eq_id++) {

      cs_equation_t  *eq = cs_equation_by_id(eq_id);

      if (!cs_equation_is_steady(eq)) {

        cs_equation_type_t  type = cs_equation_get_type(eq);

        if (type == CS_EQUATION_TYPE_USER) {

          if (cs_equation_uses_new_mechanism(eq))
            cs_equation_solve(domain->mesh, domain->dt_cur, eq);

          else { /* Deprecated */

            /* Define the algebraic system */
            cs_equation_build_system(domain->mesh,
                                     domain->time_step,
                                     domain->dt_cur,
                                     eq);

            /* Solve domain */
            cs_equation_solve_deprecated(eq);

          }

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
  bool  do_output = cs_domain_needs_log(domain);

  /* Output information */
  if (domain->only_steady) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
    cs_log_printf(CS_LOG_DEFAULT, "#      Solve steady-state problem(s)\n");
    cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
  }
  else if (do_output) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
    cs_log_printf(CS_LOG_DEFAULT,
                  "-ite- 0; >> Solve only steady-state equations if needed");
    cs_log_printf(CS_LOG_DEFAULT, "\n%s\n", lsepline);
  }

  /* Predefined equation for the computation of the wall distance */
  if (cs_walldistance_is_activated())
    cs_walldistance_compute(domain->mesh,
                            domain->time_step,
                            domain->connect,
                            domain->cdo_quantities);

  /* If the problem is globally unsteady, only steady-state equations are
     solved */
  if (cs_gwf_is_activated())
    cs_gwf_compute_steady_state(domain->mesh,
                                domain->time_step,
                                domain->connect,
                                domain->cdo_quantities);

  if (cs_navsto_system_is_activated())
    cs_navsto_system_compute_steady_state(domain->mesh);

  /* User-defined equations */
  _compute_steady_user_equations(domain);

  /* Extra operations and post-processing of the computed solutions */
  cs_domain_post(domain);

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
  int  nt_cur = domain->time_step->nt_cur;
  bool  do_output = cs_domain_needs_log(domain);

  /* Output information */
  if (do_output) {

    double  t_cur = domain->time_step->t_cur;

    cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
    cs_log_printf(CS_LOG_DEFAULT,
                  "-ite- %d >> Solve domain from time=%6.4e to %6.4e;"
                  " dt=%5.3e",
                  nt_cur, t_cur, t_cur + domain->dt_cur, domain->dt_cur);
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
  }

  if (cs_gwf_is_activated())
    cs_gwf_compute(domain->mesh,
                   domain->time_step,
                   domain->dt_cur,
                   domain->connect,
                   domain->cdo_quantities);

  if (cs_navsto_system_is_activated())
    cs_navsto_system_compute(domain->mesh, domain->dt_cur);

  /* User-defined equations */
  _compute_unsteady_user_equations(domain, nt_cur);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary the setup of all major structures:
 *          cs_domain_t structure,
 *          all equations and all properties
 *
 * \param[in]   domain    pointer to the cs_domain_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

static void
_log_setup(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return;

  cs_cdo_connect_summary(domain->connect);
  cs_cdo_quantities_summary(domain->cdo_quantities);

  /* Output information */
  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of domain settings\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);

  /* Boundaries of the domain */
  cs_domain_boundary_log_setup();

  /* Time step summary */
  cs_log_printf(CS_LOG_SETUP, "\n  Time step information\n");
  if (domain->only_steady)
    cs_log_printf(CS_LOG_SETUP, "  >> Steady-state computation");

  else { /* Time information */

    cs_log_printf(CS_LOG_SETUP, "  >> Time step status:");
    if (domain->time_options.idtvar == 0)
      cs_log_printf(CS_LOG_SETUP, "  constant\n");
    else if (domain->time_options.idtvar == 1)
      cs_log_printf(CS_LOG_SETUP, "  variable in time\n");
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid idtvar value for the CDO module.\n"));

    cs_xdef_log(domain->time_step_def);

    if (domain->time_step->t_max > 0.)
      cs_log_printf(CS_LOG_SETUP, "%-30s %5.3e\n",
                    "  >> Final simulation time:", domain->time_step->t_max);
    if (domain->time_step->nt_max > 0)
      cs_log_printf(CS_LOG_SETUP, "%-30s %9d\n",
                    "  >> Final time step:", domain->time_step->nt_max);

  }
  cs_log_printf(CS_LOG_SETUP, "\n");

  /* Summary for each equation */
  cs_equation_log_setup();

  if (domain->verbosity > 0) {

    /* Properties */
    cs_property_log_setup();

    /* Advection fields */
    cs_advection_field_log_setup();

    /* Summary of the groundwater module */
    cs_gwf_log_setup();

    /* Summary of the Navier-Stokes system */
    cs_navsto_system_log_setup();

  } /* Domain->verbosity > 0 */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the computational domain when CDO/HHO schemes are
 *         activated
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
  cs_cdo_ts_id = cs_timer_stats_create("stages", "cdo", "cdo");
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Store the fact that the CDO/HHO module is activated */
  cs_domain_cdo_log(domain);

  /* Add predefined properties */
  cs_property_t  *pty = cs_property_add("unity", CS_PROPERTY_ISO);

  cs_property_def_iso_by_value(pty, "cells", 1.0);

  cs_timer_t t0 = cs_timer_time();

  /* Add a boundary zone gathering all "wall" boundaries */
  if (cs_navsto_system_is_activated() ||
      cs_walldistance_is_activated())
    cs_domain_boundary_def_wall_zones();

  /* According to the settings, add or not predefined equations:
      >> Wall distance
      >> Groundwater flows
      >> Mesh deformation
      >> Navier-Stokes system
  */
  cs_domain_setup_predefined_equations(domain);

  /* Add variables related to user-defined and predefined equations */
  cs_equation_create_fields();
  cs_advection_field_create_fields();

  /* Set the scheme flag for the computational domain */
  cs_domain_set_scheme_flags(domain);

  /* Monitoring */
  cs_timer_stats_stop(cs_cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(domain->tcs, domain->tcs, time_count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the structures related to the computational domain when
 *         CDO/HHO schemes are activated
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
  if (cs_domain_get_cdo_mode(domain) == CS_DOMAIN_CDO_MODE_OFF)
    return;

  cs_timer_t  t0 = cs_timer_time();

  /* Timer statistics */
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Last setup stage */
  cs_domain_finalize_setup(domain, m, mq);

  /* Initialization default post-processing for the computational domain */
  cs_domain_post_init(domain);

  /* Summary of the settings */
  _log_setup(domain);

  /* Output information */
  cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "#      Start main loop\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);

  /* Flush listing and setup.log files */
  cs_log_printf_flush(CS_LOG_DEFAULT);
  cs_log_printf_flush(CS_LOG_SETUP);
  cs_log_printf_flush(CS_LOG_PERFORMANCE);

  /* Monitoring */
  cs_timer_stats_stop(cs_cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(domain->tcs, domain->tcs, time_count);

  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Setup> Runtime", domain->tcs.wall_nsec*1e-9);
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
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Finalize user-defined extra operations */
  cs_user_cdo_end_extra_op(domain);

  /* Write a restart file if needed */
  cs_domain_write_restart(domain);

  /* Print monitoring information */
  cs_equation_log_monitoring();

  /* Free memory related to equations */
  cs_equation_destroy_all();

  /* Free memory related to advection fields */
  cs_advection_field_destroy_all();

  /* Free memory related to properties */
  cs_property_destroy_all();

  /* Free memory related to the groundwater flow module */
  cs_gwf_destroy_all();

  /* Navier-Stokes system */
  cs_navsto_system_destroy();

  /* Free common structures relatated to equations */
  cs_equation_common_free(domain->cdo_context);

  /* Set flag to OFF */
  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_OFF);

  cs_log_printf(CS_LOG_DEFAULT,
                "\n  Finalize and free CDO-related structures.\n");

  cs_timer_stats_stop(cs_cdo_ts_id);
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

  cs_timer_t t0 = cs_timer_time();

  /* Timer statistics */
  cs_timer_stats_start(cs_cdo_ts_id);

  /*  Build high-level structures and create algebraic systems
      Set the initial values of the fields and properties */
  cs_domain_initialize_systems(domain);

  /* Read a restart file if needed */
  cs_domain_read_restart(domain);

  /* Initialization for user-defined extra operations. Should be done
     after the domain initialization if one wants to overwrite the field
     initialization for instance */
  cs_user_cdo_start_extra_op(cs_glob_domain);

  /* Build and solve equations related to the computational domain in case of
     steady-state equations */
  _solve_steady_state_domain(domain);

  /* Main time loop */
  if (domain->time_step->nt_cur == 0)
    domain->time_step->nt_cur = 1; /* Same numbering as Finite Volume part */

  while (cs_domain_needs_iteration(domain)) {

    /* Define the current time step */
    cs_domain_define_current_time_step(domain);

    /* Build and solve equations related to the computational domain */
    _solve_domain(domain);

    /* Increment time */
    cs_domain_increment_time(domain);

    /* Extra operations and post-processing of the computed solutions */
    cs_domain_post(domain);

    /* Increment time */
    cs_domain_increment_time_step(domain);

    /* Read a control file if present */
    cs_control_check_file();

    /* Add a checkpoint if needed */
    cs_domain_write_restart(domain);

    cs_timer_stats_increment_time_step();

  }

  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Post> Runtime", domain->tcp.wall_nsec*1e-9);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(time_count, cs_glob_domain->tcs, time_count);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO> Total runtime", time_count.wall_nsec*1e-9);

  cs_timer_stats_stop(cs_cdo_ts_id);
  if (cs_glob_rank_id <= 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
    cs_log_printf(CS_LOG_DEFAULT, "#\tExit CDO core module\n");
    cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
