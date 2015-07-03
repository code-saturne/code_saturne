/*============================================================================
 * Routines for solving equations with CDO discretizations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include <fvm_defs.h>

#include <cs_base.h>
#include <cs_timer.h>
#include <cs_log.h>
#include <cs_post.h>
#include <cs_restart.h>
#include <cs_time_plot.h>
#include <cs_ext_neighborhood.h>
#include <cs_prototypes.h>
#include <cs_mesh_location.h>

/* CDO module */
#include "cs_cdo.h"
#include "cs_quadrature.h"
#include "cs_param.h"
#include "cs_param_eq.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_sla.h"
#include "cs_cdovb_codits.h"
#include "cs_cdofb_codits.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_main.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOEQ_VB  0
#define CS_CDOEQ_FB  1
#define CS_N_TYPES_OF_CDOEQS 2

/*=============================================================================
 * Local constant and enum definitions
 *============================================================================*/

static const char cs_cdoversion[] = "0.1";

static  _Bool  cs_cdo_do_navsto = false;
static  int  cs_cdo_n_equations = 0;
static  int  n_cdo_equations_by_type[CS_N_TYPES_OF_CDOEQS];

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize algebraic for equation to solve
 *
 * \param[in]  m         pointer to a cs_mesh_t struct.
 * \param[in]  mq        pointer to a cs_quantities_t struct.
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_create_algebraic_systems(const cs_mesh_t             *m,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *cdoq)
{
  int  i, eq_id;
  int  cdo_eq_counters[CS_N_TYPES_OF_CDOEQS];

  for (i = 0; i < CS_N_TYPES_OF_CDOEQS; i++) {
    n_cdo_equations_by_type[i] = 0;
    cdo_eq_counters[i] = 0;
  }

  /* Retrieve general information */
  cs_param_eq_get_info(&cs_cdo_do_navsto, &cs_cdo_n_equations);

  /* First loop: count the number of equations of each case */
  for (eq_id = 0; eq_id < cs_cdo_n_equations; eq_id++) {

    cs_space_scheme_t  space_scheme = cs_param_eq_get_space_scheme(eq_id);

    /* Up to now only this type of equation is handled */
    assert(cs_param_eq_get_type(eq_id) == CS_PARAM_EQ_TYPE_SCAL);

    /* Build algebraic system */
    switch (space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      n_cdo_equations_by_type[CS_CDOEQ_VB] += 1;
      break;
    case CS_SPACE_SCHEME_CDOFB:
      n_cdo_equations_by_type[CS_CDOEQ_FB] += 1;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid space scheme. Stop creating algebraic systems \n"));

    } /* space_scheme */

  } /* Loop on equations */

  /* Allocate structures */
  cs_cdovb_codits_create_all(n_cdo_equations_by_type[CS_CDOEQ_VB]);
  cs_cdofb_codits_create_all(n_cdo_equations_by_type[CS_CDOEQ_FB]);

  /* Initialize algebraic system related to each equation */
  for (eq_id = 0; eq_id < cs_cdo_n_equations; eq_id++) {

    const cs_param_eq_t  *eq = cs_param_eq_get_by_id(eq_id);

    /* Up to now only this type of equation is handled */
    assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);
    assert(eq->algo_info.type == CS_PARAM_EQ_ALGO_ITSOL);

    /* Build algebraic system */
    switch (eq->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      bft_printf("\n -cdo- SpaceDiscretization >> %s >> CDO.VB\n",
                 eq->name);
      cs_cdovb_codits_init(eq, m, cdo_eq_counters[CS_CDOEQ_VB]);
      cs_sla_itsol_update_buffer_sizes(cdoq->n_vertices, eq->itsol_info);
      cdo_eq_counters[CS_CDOEQ_VB] += 1;
      break;

    case CS_SPACE_SCHEME_CDOFB:
      bft_printf("\n -cdo- SpaceDiscretization >> %s >> CDO.FB\n",
                 eq->name);
      cs_cdofb_codits_init(eq, m, cdo_eq_counters[CS_CDOEQ_FB]);
      cs_sla_itsol_update_buffer_sizes(cdoq->n_faces, eq->itsol_info);
      cdo_eq_counters[CS_CDOEQ_FB] += 1;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid space scheme. Stop creating algebraic systems.\n"));

    } /* space_scheme */

  } /* Loop on equations */

  cs_sla_itsol_alloc_buffers();

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve Navier-Stokes equations and/or additional equations
 *
 * \param[in]  m         pointer to a cs_mesh_t struct.
 * \param[in]  mq        pointer to a cs_quantities_t struct.
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]  cdoq      pointer to a cs_cdo_quantities_t struct.
 * \param[in]  tcur      current physical time of the simulation
 */
/*----------------------------------------------------------------------------*/

static void
_solve(const cs_mesh_t              *m,
       const cs_mesh_quantities_t   *mq,
       const cs_cdo_connect_t       *connect,
       const cs_cdo_quantities_t    *cdoq,
       double                        tcur)
{
  int  i, eq_id;
  int  cdo_eq_counters[CS_N_TYPES_OF_CDOEQS];

  for (i = 0; i < CS_N_TYPES_OF_CDOEQS; i++)
    cdo_eq_counters[i] = 0;

  /* Solve each equation */
  for (eq_id = 0; eq_id < cs_cdo_n_equations; eq_id++) {

    cs_space_scheme_t  space_scheme = cs_param_eq_get_space_scheme(eq_id);

    bft_printf("\n");
    bft_printf(lsepline);
    bft_printf("  Solve equation %s\n", cs_param_eq_get_name(eq_id));
    bft_printf(lsepline);

    /* Up to now only this type of equation is handled */
    assert(cs_param_eq_get_type(eq_id) == CS_PARAM_EQ_TYPE_SCAL);

    /* Build algebraic system */
    switch (space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      cs_cdovb_codits_solve(m, connect, cdoq,
                            tcur,
                            cdo_eq_counters[CS_CDOEQ_VB]);
      cs_cdovb_codits_post(connect, cdoq, cdo_eq_counters[CS_CDOEQ_VB]);
      cdo_eq_counters[CS_CDOEQ_VB] += 1;
      break;

    case CS_SPACE_SCHEME_CDOFB:
      cs_cdofb_codits_solve(m, connect, cdoq,
                            tcur,
                            cdo_eq_counters[CS_CDOEQ_FB]);
      cs_cdofb_codits_post(connect, cdoq, cdo_eq_counters[CS_CDOEQ_FB]);
      cdo_eq_counters[CS_CDOEQ_FB] += 1;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid space scheme. Stop solving algebraic systems.\n"));

    } /* space_scheme */

  } /* Loop on equations */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Main program for running a simulation with CDO kernel
 *
 * \param[inout]  m     pointer to a cs_mesh_t struct.
 * \param[in]     mq    pointer to a cs_quantities_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_main(cs_mesh_t             *m,
            cs_mesh_quantities_t  *mq)
{
  int  i, time_iter;
  cs_timer_t  t0, t1;
  cs_timer_counter_t  time_count;

  // TODO: add time managment
  int  n_time_steps = 1;
  double  dt = 0.;

  /* Build high-level structures */
  t0 = cs_timer_time();

  /* Output information */
  bft_printf("\n");
  bft_printf(lsepline);
  bft_printf("     Start CDO Module  *** Experimental ***\n");
  bft_printf(lsepline);
  bft_printf("\n -msg- Version.Tag  %s\n", cs_cdoversion);

  /* Build additional connectivity using DEC matrices */
  cs_cdo_connect_t  *connect = cs_cdo_connect_build(m);

  cs_cdo_connect_resume(connect);

  /* Build additional mesh quantities in a seperate structure */
  cs_cdo_quantities_t  *cdoq = cs_cdo_quantities_build(m, mq, connect);

  cs_param_pty_set_default();

  /* User-defined settings (Keep this order of call) */
  cs_user_cdo_setup();            // initial setup
  cs_user_cdo_numeric_settings(); // advanced setup
  cs_user_cdo_hodge_settings();   // advanced setup

  /* Add variables related to user-defined equations */
  cs_param_eq_add_fields();

  /* Add user-defined material properties */
  cs_param_pty_add_fields();

  /* Build all new mesh locations which are not set yet */
  cs_mesh_location_build(m, -1);

  /* Resume the settings */
  cs_param_pty_resume_all();
  cs_param_eq_resume_all();

  /* Initialize post-processing */
  cs_post_activate_writer(-1,     /* default writer (volume mesh)*/
                          true);  /* activate if 1 */
  cs_post_write_meshes(NULL);     /* time step management structure set to NULL
                                     => Time-idenpendent output is considered */

  /* Initialization of several modules */
  cs_set_eps_machine();      /* Compute and set epsilon machine */
  cs_quadrature_setup();     /* Compute constant used in quadrature rules */
  cs_toolbox_init(4*m->n_cells);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  -t-    CDO setup runtime                    %12.3f s\n"),
                time_count.wall_nsec*1e-9);

  /* Create algebraic systems */
  t0 = cs_timer_time();

  _create_algebraic_systems(m, connect, cdoq);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  -t-    Creation of CDO systems              %12.3f s\n"),
                time_count.wall_nsec*1e-9);

  /* Loop on time iterations */
  for (time_iter = 0; time_iter < n_time_steps; time_iter++) {

    double  tcur = time_iter*dt;

    /* Solve linear systems */
    t0 = cs_timer_time();

    _solve(m, mq, connect, cdoq, tcur);

    t1 = cs_timer_time();
    time_count = cs_timer_diff(&t0, &t1);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  -t-    CDO solver runtime (iter: %d)        %12.3f s\n"),
                  time_iter, time_count.wall_nsec*1e-9);

    /* Extra operations */
    t0 = cs_timer_time();

    cs_user_cdo_extra_op(m, mq, connect, cdoq, tcur);

    t1 = cs_timer_time();
    time_count = cs_timer_diff(&t0, &t1);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  -t-    CDO extra op. (iter: %d)             %12.3f s\n"),
                  time_iter, time_count.wall_nsec*1e-9);

  } /* Loop on time steps */

  /* Free main CDO structures */
  t0 = cs_timer_time();

  cs_cdo_n_equations = 0;
  for (i = 0; i < CS_N_TYPES_OF_CDOEQS; i++)
    n_cdo_equations_by_type[i] = 0;

  cs_cdovb_codits_free_all();
  cs_cdofb_codits_free_all();
  cs_param_pty_free_all();
  cs_param_eq_free_all();

  cs_toolbox_finalize();
  cs_sla_itsol_free_buffers();

  cdoq = cs_cdo_quantities_free(cdoq);
  connect = cs_cdo_connect_free(connect);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  -t-    Free CDO structures                  %12.3f s\n"),
                time_count.wall_nsec*1e-9);

  bft_printf("\n ============================================\n");
  bft_printf("                Exit CDO Module\n");
  bft_printf(" ============================================\n");
  printf("\n  --> Exit CDO module\n\n");

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
