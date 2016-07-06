/*============================================================================
 * Routines for solving equations with CDO discretizations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "bft_printf.h"

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

/* CDO module */
#include "cs_cdo.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_domain.h"
#include "cs_equation.h"
#include "cs_groundwater.h"
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

static const char cs_cdoversion[] = "0.6";

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

  /* Initialization of several modules */
  cs_math_set_machine_epsilon(); /* Compute and set machine epsilon */
  cs_quadrature_setup();         /* Compute constant used in quadrature rules */

  /* Output information */
  bft_printf("\n");
  bft_printf("%s", lsepline);
  bft_printf("\tStart CDO Module  *** Experimental ***\n");
  bft_printf("%s", lsepline);
  bft_printf("\n -msg- Version.Tag  %s\n", cs_cdoversion);

  cs_timer_t t0 = cs_timer_time();
  cs_timer_stats_start(cdo_ts_id);

  /*  Build high-level structures and create algebraic systems */
  cs_domain_t  *domain = _setup_domain(m, mq);

  bft_printf("\n%s", lsepline);
  bft_printf("      Start main loop on time iteration\n");
  bft_printf("%s", lsepline);

  while (cs_domain_needs_iterate(domain)) { // Main time loop

    /* Define the current time step */
    cs_domain_define_current_time_step(domain);

    /* Build and solve equations related to the computational domain */
    cs_domain_solve(domain);

    /* Extra operations and post-processing of the computed solutions */
    cs_domain_postprocess(domain);

    /* Increment time */
    cs_domain_increment_time(domain);
    cs_timer_stats_increment_time_step();

  }

  /* Free main CDO structures */
  _finalize(domain);

  bft_printf("\n%s", lsepline);
  bft_printf("\tExit CDO Module\n");
  bft_printf("%s", lsepline);

  cs_timer_stats_stop(cdo_ts_id);
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                "t--> CDO total runtime                 %12.3f s\n",
                time_count.wall_nsec*1e-9);

  printf("\n  --> Exit: simulation completed for the CDO module\n\n");

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
