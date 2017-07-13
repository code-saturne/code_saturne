/*============================================================================
 * Routines for solving equations with CDO discretizations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"

#include "cs_defs.h"
#include "cs_base.h"
#include "cs_math.h"
#include "cs_timer.h"
#include "cs_log.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_mesh_location.h"
#include "cs_timer_stats.h"
#include "cs_volume_zone.h"

/* CDO module */
#include "cs_cdo.h"
#include "cs_domain.h"
#include "cs_equation.h"
#include "cs_gwf.h"
#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_sla.h"
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

static const char cs_cdoversion[] = "0.9.2";
static int  cs_cdo_ts_id;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the computational domain when CDO/HHO schemes are
 *         activated
 *
 * \param[in]  activation_mode   integer for tagging a mode
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_setup(int  activation_mode)
{
  cs_cdo_activation_mode = activation_mode;

  if (activation_mode == CS_CDO_OFF)
    return;

  /* Timer statistics */
  cs_cdo_ts_id = cs_timer_stats_create("stages", "cdo", "cdo");
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Store the fact that the CDO/HHO module is activated */
  cs_log_printf(CS_LOG_DEFAULT,
                "\n -msg- CDO/HHO module is activated *** Experimental ***");
  cs_log_printf(CS_LOG_DEFAULT, "\n -msg- CDO.tag  %s\n", cs_cdoversion);

  cs_timer_t t0 = cs_timer_time();

  /* Initialization of several modules */
  cs_math_set_machine_epsilon(); /* Compute and set machine epsilon */
  cs_quadrature_setup();         /* Compute constant used in quadrature rules */

  /* User-defined settings and default initializations
     WARNING: Change the order of call to the following routines with care
              This may incur bugs */
  cs_user_cdo_add_mesh_locations();

  /* Create a new structure for the computational domain */
  cs_glob_domain = cs_domain_create();

  /* Initialize the new computational domain
     - Set the default boundary and potentially add new boundary to the
       computational domain
     - Add predefined and user-defined equations
     - Activate CDO-related submodules
     - Add user-defined properties and advection fields
     - Set the time step
  */
  cs_user_cdo_init_setup(cs_glob_domain);

  /* Update mesh locations */
  cs_domain_update_mesh_locations(cs_glob_domain);

  /* Advanced settings (numerical scheme, hodge operators, solvers...).
     This call must be done before the field creation since the support
     of variable field depends on the choice of the numerical scheme. */
  cs_user_cdo_numeric_settings();

  /* Add variables related to user-defined and predefined equations */
  cs_equation_create_fields();
  cs_advection_field_create_fields();

  /* According to the settings, add or not predefined equations:
      >> Wall distance
      >> Groundwater flows
  */
  cs_domain_setup_predefined_equations(cs_glob_domain);

  /* Overwrite predefined settings if this is what user wants */
  cs_user_cdo_numeric_settings();

  /* Set the scheme flag for the computational domain */
  cs_domain_set_scheme_flag(cs_glob_domain);

  /* Monitoring */
  cs_timer_stats_stop(cs_cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(cs_glob_domain->tcs, cs_glob_domain->tcs, time_count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the structures related to the computational domain when
 *         CDO/HHO schemes are activated
 *
 * \param[in, out]  m     pointer to a cs_mesh_t struct.
 * \param[in]       mq    pointer to a cs_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_initialize_structures(cs_mesh_t             *m,
                             cs_mesh_quantities_t  *mq)
{
  if (cs_cdo_activation_mode == CS_CDO_OFF)
    return;

  cs_timer_t  t0 = cs_timer_time();

  /* Timer statistics */
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Last setup stage */
  cs_domain_finalize_setup(cs_glob_domain, m, mq);

  /* Initialization for user-defined extra operations */
  cs_user_cdo_start_extra_op(cs_glob_domain);

  /* Sumary of the settings */
  cs_cdo_connect_summary(cs_glob_domain->connect);
  cs_cdo_quantities_summary(cs_glob_domain->cdo_quantities);
  cs_domain_summary(cs_glob_domain);

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

  CS_TIMER_COUNTER_ADD(cs_glob_domain->tcs, cs_glob_domain->tcs, time_count);

  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Setup> Runtime", cs_glob_domain->tcs.wall_nsec*1e-9);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free all structures allocated during the resolution of CDO/HHO
 *          schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_finalize(void)
{
  if (cs_cdo_activation_mode == CS_CDO_OFF)
    return;

  cs_domain_t  *domain = cs_glob_domain;

  /* Timer statistics */
  cs_timer_stats_start(cs_cdo_ts_id);

  /* Finalize user-defined extra operations */
  cs_user_cdo_end_extra_op(domain);

  /* Write a restart file */
  cs_domain_write_restart(domain);

  /* Free cs_domain_structure (imply several operations to free memory) */
  domain = cs_domain_free(domain);
  assert(domain == NULL);
  cs_glob_domain = NULL;

  cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "#\tExit CDO Module\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
  cs_log_printf_flush(CS_LOG_DEFAULT);

  cs_timer_stats_stop(cs_cdo_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main program for running a simulation with CDO kernel
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_main(void)
{
  assert(cs_cdo_activation_mode != CS_CDO_OFF);

  cs_timer_t t0 = cs_timer_time();

  /* Timer statistics */
  cs_timer_stats_start(cs_cdo_ts_id);

  /*  Build high-level structures and create algebraic systems */
  cs_domain_t  *domain = cs_glob_domain;

  if (cs_restart_present())
    cs_domain_read_restart(domain);

  cs_domain_initialize_systems(domain);

  while (cs_domain_needs_iteration(domain)) { // Main time loop

    /* Define the current time step */
    cs_domain_define_current_time_step(domain);

    /* Build and solve equations related to the computational domain */
    cs_domain_solve(domain);

    /* Extra operations and post-processing of the computed solutions */
    cs_domain_process_after_solve(domain);

    /* Increment time */
    cs_domain_increment_time(domain);
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
  if (cs_glob_rank_id <= 0)
    printf("\n  --> Exit CDO module: simulation completed.\n\n");

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
