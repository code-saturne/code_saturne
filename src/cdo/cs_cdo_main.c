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
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*=============================================================================
 * Local constant and enum definitions
 *============================================================================*/

static const char cs_cdoversion[] = "0.9";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup a computational domain within the CDO framework
 *
 * \param[in, out]  m     pointer to a cs_mesh_t struct.
 * \param[in]       mq    pointer to a cs_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_t *
_setup_domain(cs_mesh_t             *m,
              cs_mesh_quantities_t  *mq)
{
  cs_timer_t t0 = cs_timer_time();

  /* Initialization of several modules */
  cs_math_set_machine_epsilon(); /* Compute and set machine epsilon */
  cs_quadrature_setup();         /* Compute constant used in quadrature rules */

  /* User-defined settings and default initializations
     WARNING: Change the order of call to the following routines with care
              This may incur bugs */

  /* Determine which location are already built */
  int n_mesh_locations_ini = cs_mesh_location_n_locations();

  /* Potentially add new mesh locations */
  cs_user_cdo_add_mesh_locations();

  /* Build all new mesh locations which are not set yet */
  int n_mesh_locations = cs_mesh_location_n_locations();
  for (int  i = n_mesh_locations_ini; i < n_mesh_locations; i++)
    cs_mesh_location_build(m, i);
  n_mesh_locations_ini = cs_mesh_location_n_locations();

  /* - Create and initialize a new computational domain
     - Set the default boundary and potentially add new boundary to the
       computational domain
     - Add predefined and user-defined equations
     - Set the time step
  */
  cs_domain_t  *domain = cs_domain_init(m, mq);

  /* Build all new mesh locations which are not set yet */
  n_mesh_locations = cs_mesh_location_n_locations();
  for (int  i = n_mesh_locations_ini; i < n_mesh_locations; i++)
    cs_mesh_location_build(m, i);

  cs_volume_zone_build_all(true);

  /* Advanced settings (numerical scheme, hodge operators, solvers...).
     This call must be done before the field creation since the support
     of variable field depends on the choice of the numerical scheme. */
  cs_user_cdo_numeric_settings(domain);

  /* Add variables related to user-defined and predefined equations */
  cs_domain_create_fields(domain);

  /* According to the settings, add or not predefined equations:
      >> Wall distance
      >> Groundwater flows
  */
  cs_domain_setup_predefined_equations(domain);
  cs_user_cdo_numeric_settings(domain); /* Overwrite predefined settings
                                           if this is what user wants */

  /* Set the definition of user-defined properties and/or advection
     fields */
  cs_user_cdo_set_domain(domain);

  /* Initialize post-processing */
  cs_post_activate_writer(-1,     /* default writer (volume mesh)*/
                          true);  /* activate if 1 */
  cs_post_write_meshes(NULL);     /* time step management structure set to NULL
                                     => Time-independent output is considered */

  /* Last setup stage */
  cs_domain_last_setup(domain);

  /* Initialization for user-defined extra operations */
  cs_user_cdo_start_extra_op(domain);

  /* Sumary of the settings */
  cs_cdo_connect_summary(domain->connect);
  cs_cdo_quantities_summary(domain->cdo_quantities);
  cs_domain_summary(domain);

  /* Monitoring */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Setup> Runtime", time_count.wall_nsec*1e-9);

  return domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free all structure allocated during the resolution with CDO schemes
 *
 * \param[in, out]  domain  pointer to a cs_domain_t structure pointer
 */
/*----------------------------------------------------------------------------*/

static void
_finalize(cs_domain_t  *domain)
{
  /* Finalize user-defined extra operations */
  cs_user_cdo_end_extra_op(domain);

  /* Write a restart file */
  cs_domain_write_restart(domain);

  /* Free cs_domain_structure (imply severals operation to free memory) */
  domain = cs_domain_free(domain);
  assert(domain == NULL);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main program for running a simulation with CDO kernel
 *
 * \param[in, out]  m     pointer to a cs_mesh_t struct.
 * \param[in]       mq    pointer to a cs_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_main(cs_mesh_t             *m,
            cs_mesh_quantities_t  *mq)
{
  /* Timer statistics */
  const int  cdo_ts_id = cs_timer_stats_create("stages", "cdo", "cdo");

  /* Output information */
  cs_log_printf(CS_LOG_DEFAULT, "\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "\tStart CDO Module  *** Experimental ***\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "\n -msg- Version.Tag  %s\n", cs_cdoversion);
  cs_log_printf(CS_LOG_SETUP,"\n <cdo-settings>\n");

  cs_timer_t t0 = cs_timer_time();
  cs_timer_stats_start(cdo_ts_id);

  /*  Build high-level structures and create algebraic systems */
  cs_domain_t  *domain = _setup_domain(m, mq);

  cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "      Start main loop on time iteration\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);

  /* Flush listing and setup.log files */
  cs_log_printf_flush(CS_LOG_DEFAULT);
  cs_log_printf_flush(CS_LOG_SETUP);

  while (cs_domain_needs_iterate(domain)) { // Main time loop

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

  /* Free main CDO structures */
  _finalize(domain);

  cs_log_printf(CS_LOG_SETUP,"\n<cdo-settings>\n");
  cs_log_printf(CS_LOG_DEFAULT, "\n%s", lsepline);
  cs_log_printf(CS_LOG_DEFAULT, "\tExit CDO Module\n");
  cs_log_printf(CS_LOG_DEFAULT, "%s", lsepline);
  cs_log_printf_flush(CS_LOG_DEFAULT);

  cs_timer_stats_stop(cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO> Total runtime", time_count.wall_nsec*1e-9);

  if (cs_glob_rank_id <= 0)
    printf("\n  --> Exit CDO module: simulation completed.\n\n");

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
