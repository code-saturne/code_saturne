/*============================================================================
 * Gradient reconstruction quality tests
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_ext_neighborhood.h"
#include "cs_gradient.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gradient_quality.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print a separator line in a log file
 *----------------------------------------------------------------------------*/

static void
_print_separator(void)
{
  int i;
  char separator[81];

  for (i = 0; i < 80; i++)
    separator[i] = '-';
  separator[80] = '\0';

  bft_printf("%s\n", separator);
}

/*----------------------------------------------------------------------------
 * Define variable values for test using sin(x+2y+3z) function.
 *
 * parameters:
 *   m     <-- pointer to mesh structure
 *   mq    <-- pointer to mesh quantities structure
 *   var   --> pointer to cell values array
 *   coefa --> pointer to coefa boundary condition array
 *   coefb --> pointer to coefb boundary condition array
 *----------------------------------------------------------------------------*/

static void
_sine_x_2y_3z_test_values(const cs_mesh_t             *m,
                          const cs_mesh_quantities_t  *mq,
                          cs_real_t                   *var,
                          cs_real_t                   *coefa,
                          cs_real_t                   *coefb)
{
  cs_lnum_t cell_id, face_id;
  double xx, yy, zz;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;

  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    xx = cell_cen[cell_id][0];
    yy = cell_cen[cell_id][1];
    zz = cell_cen[cell_id][2];
    var[cell_id] = sin(xx + 2.*yy + 3.*zz);
  }

  for (face_id = 0; face_id < n_b_faces; face_id++) {
    xx = b_face_cog[face_id][0];
    yy = b_face_cog[face_id][1];
    zz = b_face_cog[face_id][2];
    coefa[face_id] = sin(xx + 2.*yy + 3.*zz);
    coefb[face_id] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Compute error for test using sin(x+2y+3z) function.
 *
 * parameters:
 *   m     <-- pointer to mesh structure
 *   mq    <-- pointer to mesh quantities structure
 *   grad  <-> reconstructed gradient in, absolute error out
 *----------------------------------------------------------------------------*/

static void
_sine_x_2y_3z_test_error(const cs_mesh_t             *m,
                         const cs_mesh_quantities_t  *mq,
                         cs_real_t                   *grad)
{
  cs_lnum_t cell_id;
  double xx, yy, zz;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_3_t  *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;

  for (cell_id = 0; cell_id < n_cells_ext; cell_id++) {
    xx = cell_cen[cell_id][0];
    yy = cell_cen[cell_id][1];
    zz = cell_cen[cell_id][2];
    grad[cell_id]                 -=        cos(xx + 2.*yy + 3.*zz);
    grad[n_cells_ext   + cell_id] -=  2.0 * cos(xx + 2.*yy + 3.*zz);
    grad[n_cells_ext*2 + cell_id] -=  3.0 * cos(xx + 2.*yy + 3.*zz);
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Run several quality tests for gradients
 *----------------------------------------------------------------------------*/

void
cs_gradient_quality(void)
{
  cs_lnum_t face_id;

  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_int_t  *isympa = NULL;
  cs_real_t  *var = NULL, *ktvar=NULL;
  cs_real_t *coefa = NULL, *coefb = NULL, *grad = NULL;

  assert(m != NULL);

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Initialization */

  BFT_MALLOC(isympa, n_b_faces, cs_int_t);
  BFT_MALLOC(var, n_cells_ext, cs_real_t);
  BFT_MALLOC(ktvar, n_cells_ext, cs_real_t);
  BFT_MALLOC(coefa, n_b_faces, cs_real_t);
  BFT_MALLOC(coefb, n_b_faces, cs_real_t);
  BFT_MALLOC(grad, n_cells_ext*3, cs_real_t);

  /* Symmetry type (value 0 avoids extrapolating the gradient on boundary faces). */

  for (face_id = 0; face_id < n_b_faces; face_id++)
    isympa[face_id] = 0;

  _print_separator();
  bft_printf(_("\n"
               " Checking gradient reconstruction quality\n"
               " ========================================\n\n"));
  _print_separator();

  /* Analytical function: sin(x+2y+3z) */

  _sine_x_2y_3z_test_values(m, mq, var, coefa, coefb);

  /* Activate default writer */

  cs_post_activate_writer(-1, true);

  /* default gradient options */

  const cs_int_t iale = 1; /* set ALE indicator to 1 so as to force recompute of
                              boundary cells contribution at each gradient call */

  const cs_int_t ivar = 0, inc = 1, idimtr = 0, iphydp = 0, ipond = 0, iccocg = 1;
  const cs_int_t imobil = 0, nswrgp = 100, iwarnp = 0;
  const cs_real_t epsrgp = 1.e-5, climgp = 1.5, extrap = 0.;

  /* Compute gradient of analytical function using the following options:
     imrgra = 0
     imrgra = 1 (standard neighborhood)
     imrgra = 2 (extended neighborhood)
     imrgra = 4 (extended neighborhood)
     imrgra = 3 (reduced extended neighborhood) */

  const cs_int_t imrgra[] = {0, 1, 2, 4, 3};
  const cs_int_t imligp[] = {-1, 1, 1, -1, 1};
  const char *grd_name[] = {N_("Grad_RC"),
                            N_("Grad_LSQ"),
                            N_("Grad_LSQ_Ext"),
                            N_("Grad_LSQ_RC"),
                            N_("Grad_LSQ_ExtRed")};
  const char *grd_err_name[] = {N_("Err_Grad_RC"),
                                N_("Err_Grad_LSQ"),
                                N_("Err_Grad_LSQ_Ext"),
                                N_("Err_Grad_LSQ_RC"),
                                N_("Err_Grad_LSQ_ExtRed")};

  for (int test_id = 0; test_id < 5; test_id++) {

    /* Reduce extended gradient if required */

    if (imrgra[test_id] == 3) {
      double anomax = 3.1415*0.25; /* standard default value */
      cs_ext_neighborhood_reduce(m, mq, anomax);
    }

    /* Recontruct gradient */

    CS_PROCF(cgdcel, CGDCEL) (&ivar,
                              &(imrgra[test_id]),
                              &inc,
                              &iccocg,
                              &imobil,
                              &iale,
                              &nswrgp,
                              &idimtr,
                              &iphydp,
                              &ipond,
                              &iwarnp,
                              &(imligp[test_id]),
                              &epsrgp,
                              &extrap,
                              &climgp,
                              isympa,
                              NULL,
                              NULL,
                              NULL,
                              coefa,
                              coefb,
                              var,
                              NULL,
                              grad);

    /* Postprocess gradient */

    cs_post_write_var(-1,                      /* mesh_id */
                      _(grd_name[test_id]),
                      3,                       /* var_dim */
                      false,                   /* interlace */
                      true,                    /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      -1,                      /* nt_cur_abs */
                      0.,                      /* t_cur_abs */
                      grad,                    /* cel_vals */
                      NULL,                    /* i_face_vals */
                      NULL);                   /* b_face_vals */

    /* Compute absolute error */

    _sine_x_2y_3z_test_error(m, mq, grad);

    /* Postprocess error */

    cs_post_write_var(-1,                      /* mesh_id */
                      _(grd_err_name[test_id]),
                      3,                       /* var_dim */
                      false,                   /* interlace */
                      true,                    /* use_parent */
                      CS_POST_TYPE_cs_real_t,
                      -1,                      /* nt_cur_abs */
                      0.,                      /* t_cur_abs */
                      grad,                    /* cel_vals */
                      NULL,                    /* i_face_vals */
                      NULL);                   /* b_face_vals */

  } /* End of loop on tests */

  BFT_FREE(isympa);
  BFT_FREE(var);
  BFT_FREE(ktvar);
  BFT_FREE(coefa);
  BFT_FREE(coefb);
  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
