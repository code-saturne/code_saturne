/*============================================================================
 * User solver
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_solver-heat-equation.c
 *
 * \brief Set user solver for heat equation.
 *
 * See \subpage user_solver for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set user solver
 */
/*----------------------------------------------------------------------------*/

/*! [set_solver] */
int
cs_user_solver_set(void)
{
  return 1;
}
/*! [set_solver] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Main call to user solver
 *
 * \param[in] mesh pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_user_solver(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mesh_quantities)
{
  /*! [local_variables] */

  cs_int_t    i, iter, n_iter;

  cs_real_t   x0, xL, t0, tL, L;
  cs_real_t   r;

  cs_real_t  *t = NULL, *t_old = NULL, *t_sol = NULL;

  cs_restart_t    *restart, *checkpoint;
  cs_time_plot_t  *time_plot;

  const cs_int_t   n = mesh->n_cells;
  const cs_real_t  pi = 4.*atan(1.);

  const char  var_name[] = "temperature";

  /*! [local_variables] */

  /*! [initialization] */

  /* Initialization */

  BFT_MALLOC(t, n, cs_real_t);
  BFT_MALLOC(t_old, n, cs_real_t);
  BFT_MALLOC(t_sol, n, cs_real_t);

  x0 =  1.e30;
  xL = -1.e30;

  for (i = 0; i < mesh->n_b_faces; i++) {
    cs_real_t  x_face = mesh_quantities->b_face_cog[3*i];
    if (x_face < x0) x0 = x_face;
    if (x_face > xL) xL = x_face;
  }

  L = xL - x0; /* it is assumed that dx is constant and x0 = 0, XL =1 */

  t0 = 0.;
  tL = 0.;

  r = 0.2; /* Fourier number */

  n_iter = 100000;

  for (i = 0; i < n; i++)
    t_old[i] = sin(pi*(0.5+i)/n);

  /*! [initialization] */

  /*! [restart] */

  /* ------- */
  /* Restart */
  /* ------- */

  restart = cs_restart_create("main",                 /* file name */
                              NULL,                   /* force directory */
                              CS_RESTART_MODE_READ);  /* read mode */

  cs_restart_read_section(restart,   /* restart file */
                          var_name,  /* buffer name */
                          1,         /* location id */
                          1,         /* number of values per location */
                          2,         /* value type */
                          t_old);    /* buffer */

  cs_restart_destroy(&restart);

  /*! [restart] */

  /*! [time_monitoring] */

  /* --------------- */
  /* Time monitoring */
  /* --------------- */

  time_plot = cs_time_plot_init_probe(var_name,          /* variable name */
                                      "probes_",         /* file prefix */
                                      CS_TIME_PLOT_DAT,  /* file format */
                                      true,              /* use iter. numbers */
                                      -1.,               /* force flush */
                                      0,                 /* buffer size */
                                      1,                 /* number of probes */
                                      NULL,              /* probes list */
                                      NULL,              /* probes coord. */
                                      NULL);             /* probes names */

  /*! [time_monitoring] */

  /*! [calculation] */

  /* ----------- */
  /* Calculation */
  /* ----------- */

  /* Heat equation resolution by Finite Volume method */

  for (iter = 0; iter < n_iter; iter++) {

    /* 1D Finite Volume scheme, with constant dx */

    t[0] = t_old[0] + r*(t_old[1] - 3.*t_old[0] + 2.*t0);

    for (i = 1; i < n-1;  i++)
      t[i] = t_old[i] + r*(t_old[i+1] - 2.*t_old[i] + t_old[i-1]);

    t[n-1] = t_old[n-1] + r*(2.*tL - 3.*t_old[n-1] + t_old[n-2]);

    /* Update previous value of the temperature */

    for (i = 0; i < n; i++)
      t_old[i] = t[i];

    /* Analytical solution at the current time */

    for (i = 0; i < n; i++)
      t_sol[i] = sin(pi*(0.5+i)/n)*exp(-r*pi*pi*(iter+1)/(n*n));

    /* Plot maximum temperature value (center-value) */

    cs_time_plot_vals_write(time_plot,   /* time plot structure */
                            iter,        /* current iteration number */
                            -1.,         /* current time */
                            1,           /* number of values */
                            &(t[n/2]));  /* values */

  }

  /*! [calculation] */

  /*! [checkpoint] */

  /* --------- */
  /* Chekpoint */
  /* --------- */

  checkpoint = cs_restart_create("main",                  /* file name */
                                 NULL,                    /* force directory */
                                 CS_RESTART_MODE_WRITE);  /* write mode */

  cs_restart_write_section(checkpoint,  /* restart file */
                           var_name,    /* buffer name */
                           1,           /* location id */
                           1,           /* number of values per location */
                           2,           /* value type */
                           t);          /* buffer */

  cs_restart_destroy(&checkpoint);

  /*! [checkpoint] */

  /*! [post_processing] */

  /* --------------- */
  /* Post-processing */
  /* --------------- */

  cs_post_activate_writer(-1,     /* default writer */
                          true);  /* activate if 1 */

  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    var_name,                /* variable name */
                    1,                       /* variable dimension */
                    false,                   /* interleave if true */
                    true,                    /* define on parents */
                    CS_POST_TYPE_cs_real_t,  /* type */
                    t,                       /* value on cells */
                    NULL,                    /* value on interior faces */
                    NULL,                    /* value on boundary faces */
                    NULL);                   /* time-independent output */

  cs_post_write_var(CS_POST_MESH_VOLUME,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    "solution",              /* variable name */
                    1,                       /* variable dimension */
                    false,                   /* interleave if true */
                    true,                    /* define on parents */
                    CS_POST_TYPE_cs_real_t,  /* type */
                    t_sol,                   /* value on cells */
                    NULL,                    /* value on interior faces */
                    NULL,                    /* value on boundary faces */
                    NULL);                   /* time-independent output */

  /*! [post_processing] */

  /*! [finalization] */

  /* Finalization */

  BFT_FREE(t);
  BFT_FREE(t_old);
  BFT_FREE(t_sol);

  cs_time_plot_finalize(&time_plot);

  /*! [finalization] */

  return;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
