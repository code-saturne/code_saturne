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

#include "cs_base.h"
#include "cs_timer.h"
#include "cs_log.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_mesh_location.h"

/* CDO module */
#include "cs_cdo.h"
#include "cs_quadrature.h"
#include "cs_sla.h"
#include "cs_param.h"
#include "cs_equation.h"
#include "cs_domain.h"
#include "cs_walldistance.h"
#include "cs_groundwater.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdofb_scaleq.h"

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

static const char cs_cdoversion[] = "0.4";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup a computation within the CDO framework
 *
 * \param[in, out]  m     pointer to a cs_mesh_t struct.
 * \param[in]       mq    pointer to a cs_quantities_t struct.
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_t *
_setup(cs_mesh_t             *m,
       cs_mesh_quantities_t  *mq)
{
  /* Output information */
  bft_printf("\n");
  bft_printf("%s", lsepline);
  bft_printf("\tStart CDO Module  *** Experimental ***\n");
  bft_printf("%s", lsepline);
  bft_printf("\n -msg- Version.Tag  %s\n", cs_cdoversion);

  /* Initialization of several modules */
  cs_set_eps_machine();      /* Compute and set epsilon machine */
  cs_quadrature_setup();     /* Compute constant used in quadrature rules */

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

  /* Advanced settings (numerical scheme, hodge operators, solvers...) */
  cs_user_cdo_numeric_settings(domain);

  /* Add variables related to user-defined and predefined equations */
  cs_domain_create_fields(domain);

  /* According to the settings, add or not predefined equations:
      >> Wall distance
      >> Groundwater flows
  */
  cs_domain_setup_predefined_equations(domain);

  /* Set the definition of user-defined properties and/or advection
     fields */
  cs_user_cdo_set_domain(domain);

  /* Initialize post-processing */
  cs_post_activate_writer(-1,     /* default writer (volume mesh)*/
                          true);  /* activate if 1 */
  cs_post_write_meshes(NULL);     /* time step management structure set to NULL
                                     => Time-idenpendent output is considered */

  /* Last setup stage */
  cs_domain_last_setup(domain);

  /* Sumary of the settings */
  cs_cdo_connect_summary(domain->connect);
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
_finalize(cs_domain_t  **domain)
{
  /* Write a restart file */
  cs_domain_write_restart(*domain);

  /* Free temporary buffers allocated for each kind of numerical used */
  cs_cdovb_scaleq_finalize();
  cs_cdofb_scaleq_finalize();

  *domain = cs_domain_free(*domain);
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
  cs_timer_t  t0, t1;
  cs_timer_counter_t  time_count;

  /* Build high-level structures */
  t0 = cs_timer_time();

  /* Create algebraic systems */
  cs_domain_t  *domain = _setup(m, mq);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                "  -t-    CDO setup runtime                    %12.3f s\n",
                time_count.wall_nsec*1e-9);

  while (cs_domain_needs_iterate(domain)) {

    t0 = cs_timer_time();

    /* Define the current time step */
    cs_domain_define_current_time_step(domain);

    /* Solve equations */
    cs_domain_solve(domain);

    t1 = cs_timer_time();
    time_count = cs_timer_diff(&t0, &t1);

    /* Output */
    if (domain->time_step->nt_cur % domain->output_freq == 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "  -t-    CDO solver runtime (iter: %d)        %12.3f s\n",
                    domain->time_step->nt_cur, time_count.wall_nsec*1e-9);

    /* Increment time */
    cs_domain_increment_time(domain);

  } // Time loop

  /* Free main CDO structures */
  t0 = cs_timer_time();

  _finalize(&domain);
  assert(domain == NULL);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  -t-    Free CDO structures                  %12.3f s\n"),
                time_count.wall_nsec*1e-9);

  bft_printf("\n%s", lsepline);
  bft_printf("\tExit CDO Module\n");
  bft_printf("%s", lsepline);
  printf("\n  --> Exit CDO module\n\n");

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
