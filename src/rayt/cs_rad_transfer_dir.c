/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_timer.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_dir.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_rad_transfer_dir.c */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Private variable
 *============================================================================*/

const int x = 0;
const int y = 1;
const int z = 2;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize quadratures.
 *----------------------------------------------------------------------------*/

static void
_init_quadrature(void)
{
  if (cs_glob_rad_transfer_params->sxyz == NULL)
    BFT_MALLOC(cs_glob_rad_transfer_params->sxyz,
               cs_glob_rad_transfer_params->ndirs,
               cs_real_3_t);

  if (cs_glob_rad_transfer_params->angsol == NULL)
    BFT_MALLOC(cs_glob_rad_transfer_params->angsol,
               cs_glob_rad_transfer_params->ndirs,
               cs_real_t);
}

/*----------------------------------------------------------------------------
 * Vector normalization
 *
 * parameters:
 *   vec    <-> Vector to be normalized
 *----------------------------------------------------------------------------*/

inline static void
_normve(cs_real_t vec[3])
{
  cs_real_t norm;

  norm = cs_math_3_norm(vec);

  for (int i = 0; i < 3; i++)
    vec[i] = vec[i] / norm;
}

/*----------------------------------------------------------------------------
 * Compute the solid angle associated to a given direction of the Tn quadrature
 *
 * parameters:
 *   posnod    <-- Node coordinate
 *----------------------------------------------------------------------------*/

inline static cs_real_t
_lhuilier(cs_real_t  posnod[3][3])
{
  cs_real_t a, b, c, p;
  cs_real_t sol_angle;

  /* Renormalisation for dot products    */

  _normve(posnod[0]);
  _normve(posnod[1]);
  _normve(posnod[2]);

   /* segment lengths of the curvilinear triangle (R=1)  */

  a = acos(cs_math_3_dot_product(posnod[0], posnod[1]));
  b = acos(cs_math_3_dot_product(posnod[1], posnod[2]));
  c = acos(cs_math_3_dot_product(posnod[2], posnod[0]));

   /* perimeter  */

   p = 0.5 * (a + b + c);

   /* Solid angle */

   sol_angle  = 4.0 * atan(sqrt(  tan(p / 2.0)
                                * tan((p - a) / 2.0)
                                * tan((p - b) / 2.0)
                                * tan((p - c) / 2.0)
                               ));

   return sol_angle;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a quadrature Sn or Tn
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_dir(void)
{
  cs_rad_transfer_params_t *cs_rp = cs_glob_rad_transfer_params;

  /* Initializations  */

  switch(cs_rp->i_quadrature){

  case CS_RAD_QUADRATURE_S4: /* quadrature S4 : 24 directions  */
    cs_rp->ndirs = 3;
    break;

  case CS_RAD_QUADRATURE_S6: /* quadrature S6 : 48 directions  */
    cs_rp->ndirs = 6;
    break;

  case CS_RAD_QUADRATURE_S8: /* quadrature S8 : 80 directions  */
    cs_rp->ndirs = 10;
    break;

  case 4: /* Quadrature T2 : 32 directions  */
    cs_rp->ndirs = 4;
    break;

  case 5: /* Quadrature T4 : 128 directions */
    cs_rp->ndirs = 16;
    break;

  case 6: /* Quadrature Tn : 8 n^2 directions    */
    cs_rp->ndirs = pow(cs_rp->ndirec, 2.0);
    break;

  case 7: /* Quadrature 120 directions (LC11)    */
    cs_rp->ndirs = 15;
    break;

  case 8: /* Quadrature 48 directions (DCT020-2468)   */
    cs_rp->ndirs = 6;
    break;

  default:
    assert(0);
    break;

  }

  _init_quadrature();

  cs_real_t vec[10]  = {0.0};
  cs_real_t weight[5] = {0.0};

  /* =================================================================
   * Quadrature Sn : n(n+2) directions
   * ================================================================= */

  /* quadrature S4 : 24 directions  */

  if (cs_rp->i_quadrature == 1) {
    vec[0]    = 0.2958759;
    vec[1]    = 0.9082483;
    weight[0]  = 0.5235987;
    cs_rp->sxyz[0][0] = vec[0];
    cs_rp->sxyz[1][0] = vec[0];
    cs_rp->sxyz[2][0] = vec[1];
    cs_rp->sxyz[0][1] = vec[0];
    cs_rp->sxyz[1][1] = vec[1];
    cs_rp->sxyz[2][1] = vec[0];
    cs_rp->sxyz[0][2] = vec[1];
    cs_rp->sxyz[1][2] = vec[0];
    cs_rp->sxyz[2][2] = vec[0];
    cs_rp->angsol[0]  = weight[0];
    cs_rp->angsol[1]  = weight[0];
    cs_rp->angsol[2]  = weight[0];
  }

  /* quadrature S6 : 48 directions  */

  else if (cs_rp->i_quadrature == 2) {
    vec[0]    = 0.183867;
    vec[1]    = 0.6950514;
    vec[2]    = 0.9656013;
    weight[0]  = 0.1609517;
    weight[1]  = 0.3626469;

    cs_rp->sxyz[0][0] = vec[0];
    cs_rp->sxyz[1][0] = vec[0];
    cs_rp->sxyz[2][0] = vec[0];
    cs_rp->sxyz[3][0] = vec[1];
    cs_rp->sxyz[4][0] = vec[1];
    cs_rp->sxyz[5][0] = vec[2];
    cs_rp->sxyz[0][1] = vec[0];
    cs_rp->sxyz[1][1] = vec[1];
    cs_rp->sxyz[2][1] = vec[2];
    cs_rp->sxyz[3][1] = vec[0];
    cs_rp->sxyz[4][1] = vec[1];
    cs_rp->sxyz[5][1] = vec[0];
    cs_rp->sxyz[0][2] = vec[2];
    cs_rp->sxyz[1][2] = vec[1];
    cs_rp->sxyz[2][2] = vec[0];
    cs_rp->sxyz[3][2] = vec[1];
    cs_rp->sxyz[4][2] = vec[0];
    cs_rp->sxyz[5][2] = vec[0];
    cs_rp->angsol[0]  = weight[0];
    cs_rp->angsol[1]  = weight[1];
    cs_rp->angsol[2]  = weight[0];
    cs_rp->angsol[3]  = weight[1];
    cs_rp->angsol[4]  = weight[1];
    cs_rp->angsol[5]  = weight[0];
  }

  /* quadrature S8 : 80 directions  */

  else if (cs_rp->i_quadrature == 3) {
    vec[0]     = 0.1422555;
    vec[1]     = 0.5773503;
    vec[2]     = 0.8040087;
    vec[3]     = 0.9795543;
    weight[0]   = 0.0992284;
    weight[1]   = 0.1712359;
    weight[2]   = 0.4617179;

    cs_rp->sxyz[0][0] = vec[0];
    cs_rp->sxyz[1][0] = vec[0];
    cs_rp->sxyz[2][0] = vec[0];
    cs_rp->sxyz[3][0] = vec[0];
    cs_rp->sxyz[4][0] = vec[1];
    cs_rp->sxyz[5][0] = vec[1];
    cs_rp->sxyz[6][0] = vec[1];
    cs_rp->sxyz[7][0] = vec[2];
    cs_rp->sxyz[8][0] = vec[2];
    cs_rp->sxyz[9][0] = vec[3];
    cs_rp->sxyz[0][1] = vec[0];
    cs_rp->sxyz[1][1] = vec[1];
    cs_rp->sxyz[2][1] = vec[2];
    cs_rp->sxyz[3][1] = vec[3];
    cs_rp->sxyz[4][1] = vec[0];
    cs_rp->sxyz[5][1] = vec[1];
    cs_rp->sxyz[6][1] = vec[2];
    cs_rp->sxyz[7][1] = vec[0];
    cs_rp->sxyz[8][1] = vec[1];
    cs_rp->sxyz[9][1] = vec[0];
    cs_rp->sxyz[0][2] = vec[3];
    cs_rp->sxyz[1][2] = vec[2];
    cs_rp->sxyz[2][2] = vec[1];
    cs_rp->sxyz[3][2] = vec[0];
    cs_rp->sxyz[4][2] = vec[2];
    cs_rp->sxyz[5][2] = vec[1];
    cs_rp->sxyz[6][2] = vec[0];
    cs_rp->sxyz[7][2] = vec[1];
    cs_rp->sxyz[8][2] = vec[0];
    cs_rp->sxyz[9][2] = vec[0];
    cs_rp->angsol[0] = weight[1];
    cs_rp->angsol[1] = weight[0];
    cs_rp->angsol[2] = weight[0];
    cs_rp->angsol[3] = weight[1];
    cs_rp->angsol[4] = weight[0];
    cs_rp->angsol[5] = weight[2];
    cs_rp->angsol[6] = weight[0];
    cs_rp->angsol[7] = weight[0];
    cs_rp->angsol[8] = weight[0];
    cs_rp->angsol[9] = weight[1];
  }

  /* =================================================================
   * Quadrature Tn : 8n^2 directions
   * ================================================================= */

  /* Quadrature T2 : 32 directions  */

  else if (cs_rp->i_quadrature == 4) {
    vec[0]   = 0.2357022604;
    vec[1]   = 0.9428090416;
    vec[2]   = 0.5773502692;
    weight[0] = 0.5512855984;
    weight[1] = 0.3398369095;

    cs_rp->sxyz[0][0] = vec[0];
    cs_rp->sxyz[1][0] = vec[1];
    cs_rp->sxyz[2][0] = vec[2];
    cs_rp->sxyz[3][0] = vec[0];
    cs_rp->sxyz[0][1] = vec[0];
    cs_rp->sxyz[1][1] = vec[0];
    cs_rp->sxyz[2][1] = vec[2];
    cs_rp->sxyz[3][1] = vec[1];
    cs_rp->sxyz[0][2] = vec[1];
    cs_rp->sxyz[1][2] = vec[0];
    cs_rp->sxyz[2][2] = vec[2];
    cs_rp->sxyz[3][2] = vec[0];
    cs_rp->angsol[0] = weight[1];
    cs_rp->angsol[1] = weight[1];
    cs_rp->angsol[2] = weight[0];
    cs_rp->angsol[3] = weight[1];
  }

  /* Quadrature T4 : 128 directions */

  else if (cs_rp->i_quadrature == 5) {
    vec[0]    = 0.0990147543;
    vec[1]    = 0.4923659639;
    vec[2]    = 0.2357022604;
    vec[3]    = 0.123091491;
    vec[4]    = 0.8616404369;
    vec[5]    = 0.6804138174;
    vec[6]    = 0.5773502692;
    vec[7]    = 0.272165527;
    vec[8]    = 0.990147543;
    vec[9]    = 0.9428090416;
    weight[0]  = 0.0526559083;
    weight[1]  = 0.0995720042;
    weight[2]  = 0.0880369928;
    weight[3]  = 0.1320249278;
    weight[4]  = 0.155210815;

    cs_rp->sxyz[0][0]  = vec[0];
    cs_rp->sxyz[1][0]  = vec[1];
    cs_rp->sxyz[2][0]  = vec[2];
    cs_rp->sxyz[3][0]  = vec[3];
    cs_rp->sxyz[4][0]  = vec[4];
    cs_rp->sxyz[5][0]  = vec[5];
    cs_rp->sxyz[6][0]  = vec[6];
    cs_rp->sxyz[7][0]  = vec[7];
    cs_rp->sxyz[8][0]  = vec[3];
    cs_rp->sxyz[9][0]  = vec[8];
    cs_rp->sxyz[10][0] = vec[9];
    cs_rp->sxyz[11][0] = vec[4];
    cs_rp->sxyz[12][0] = vec[5];
    cs_rp->sxyz[13][0] = vec[1];
    cs_rp->sxyz[14][0] = vec[2];
    cs_rp->sxyz[15][0] = vec[0];
    cs_rp->sxyz[0][1]  = vec[0];
    cs_rp->sxyz[1][1]  = vec[3];
    cs_rp->sxyz[2][1]  = vec[2];
    cs_rp->sxyz[3][1]  = vec[1];
    cs_rp->sxyz[4][1]  = vec[3];
    cs_rp->sxyz[5][1]  = vec[7];
    cs_rp->sxyz[6][1]  = vec[6];
    cs_rp->sxyz[7][1]  = vec[5];
    cs_rp->sxyz[8][1]  = vec[4];
    cs_rp->sxyz[9][1]  = vec[0];
    cs_rp->sxyz[10][1] = vec[2];
    cs_rp->sxyz[11][1] = vec[1];
    cs_rp->sxyz[12][1] = vec[5];
    cs_rp->sxyz[13][1] = vec[4];
    cs_rp->sxyz[14][1] = vec[9];
    cs_rp->sxyz[15][1] = vec[8];
    cs_rp->sxyz[0][2]  = vec[8];
    cs_rp->sxyz[1][2]  = vec[4];
    cs_rp->sxyz[2][2]  = vec[9];
    cs_rp->sxyz[3][2]  = vec[4];
    cs_rp->sxyz[4][2]  = vec[1];
    cs_rp->sxyz[5][2]  = vec[5];
    cs_rp->sxyz[6][2]  = vec[6];
    cs_rp->sxyz[7][2]  = vec[5];
    cs_rp->sxyz[8][2]  = vec[1];
    cs_rp->sxyz[9][2]  = vec[0];
    cs_rp->sxyz[10][2] = vec[2];
    cs_rp->sxyz[11][2] = vec[3];
    cs_rp->sxyz[12][2] = vec[7];
    cs_rp->sxyz[13][2] = vec[3];
    cs_rp->sxyz[14][2] = vec[2];
    cs_rp->sxyz[15][2] = vec[0];
    cs_rp->angsol[0]   = weight[0];
    cs_rp->angsol[1]   = weight[1];
    cs_rp->angsol[2]   = weight[2];
    cs_rp->angsol[3]   = weight[1];
    cs_rp->angsol[4]   = weight[1];
    cs_rp->angsol[5]   = weight[3];
    cs_rp->angsol[6]   = weight[4];
    cs_rp->angsol[7]   = weight[3];
    cs_rp->angsol[8]   = weight[1];
    cs_rp->angsol[9]   = weight[0];
    cs_rp->angsol[10]  = weight[2];
    cs_rp->angsol[11]  = weight[1];
    cs_rp->angsol[12]  = weight[3];
    cs_rp->angsol[13]  = weight[1];
    cs_rp->angsol[14]  = weight[2];
    cs_rp->angsol[15]  = weight[0];
  }

  /* Quadrature Tn: 8*n^2 directions     */

  else if (cs_rp->i_quadrature == 6) {

    int nquad = cs_rp->ndirec;

    /* Compute X position and Z position of the center of all the sub
     * triangles of the main 2D triangle

     * Here for T2:      z
     *                  /\
     *                 /  \
     *                /  1 \         ! level 1 (1 triangle)
     *               /______\
     *              /\      /\
     *             /  \ 2  /  \      ! level 2 (3 triangles)
     *            /  1 \  /  3 \
     *           /______\/______\
     *         x                 y
     */

    /* Max number of points of a level */

    int max_npt = 2 * nquad - 1;

    cs_real_3_t *xyz3d;
    BFT_MALLOC(xyz3d, nquad * max_npt, cs_real_3_t);

    /* z position */

    for (int ii = 0; ii < nquad; ii++) {

      int lev     = nquad - ii;
      int lev_npt = 2 * lev - 1;

      for (int jj = 0; jj <= lev_npt; jj += 2)
        xyz3d[max_npt * (lev - 1) + jj][2] = (1.0 + 3.0 * ii) / (3.0 * nquad);

      for (int jj = 1; jj <= lev_npt; jj += 2)
        xyz3d[max_npt * (lev - 1) + jj][2] = (2.0 + 3.0 * ii) / (3.0 * nquad);

    }

    /* y position
     * y position y for each point of a given level increases of 1/(3*(nbpoint+1)/2)
     * in the same column of a triangle. So we fill y3d by column (diagonal going from
     * the top to the bottom left) */

    int jj   = 0;
    int int1 = 0;

    for (int tri = 1; tri <= max_npt; tri++) {

      for (int lev = jj; lev < nquad; lev++)
        xyz3d[lev * max_npt + tri - 1][1] = (tri + jj) / (3.0 * (max_npt + 1.0) / 2.0);

      int1++;

      if (int1 == 2) {
        int1 = 0;
        jj++;

      }

    }

    /* x position */
    for (int ii = 0; ii < nquad * max_npt; ii++)
      xyz3d[ii][0]  = 1.0 - xyz3d[ii][2] - xyz3d[ii][1];

    /* write xyz3d in cs_rp->sxyz   */

    jj = 0;

    for (int lev = 1; lev <= nquad; lev++) {

      int lev_npt = 2 * lev - 1;

      for (int ii = 0; ii < lev_npt; ii++) {
        cs_rp->sxyz[jj][0] = xyz3d[(lev - 1) * max_npt + ii][0];
        cs_rp->sxyz[jj][1] = xyz3d[(lev - 1) * max_npt + ii][1];
        cs_rp->sxyz[jj][2] = xyz3d[(lev - 1) * max_npt + ii][2];
        jj++;
      }

    }

    /* Vectors normalisation     */

    for (int ii = 0; ii < cs_rp->ndirs; ii++)
      _normve(cs_rp->sxyz[ii]);

    BFT_FREE(xyz3d);
    xyz3d = NULL;

    /* All the directions are now known, weights will be computed */

    /* Compute nodes positions (the number at node are the number of the node of a
     * given level)
     *                 z
     *                 1           ! level 1 (1 point)
     *                 /\
     *                /  \
     *              1/____\2       ! level 2 (2 points)
     *              /\    /\
     *             /  \  /  \
     *           1/____\/____\3    ! level 3 (3 points)
     *          x      2       y */

    /* Max number of points for a level    */

    max_npt = nquad + 1;

    BFT_MALLOC(xyz3d, max_npt * (nquad + 1), cs_real_3_t);

    /* z position */

    int lev = 0;
    int ipt = 0;
    xyz3d[max_npt * lev + ipt][2] = 1.0;

    for (int iquad = 0; iquad < nquad; iquad++) {
      lev++;

      /* There are "lev + 1" points at level "lev". Beware, C language is 0-based */
      int lev_npt = lev + 1;
      for (ipt = 0; ipt < lev_npt; ipt++)
        xyz3d[max_npt * lev + ipt][2] = xyz3d[max_npt * (lev - 1)][2] - 1.0 / nquad;
    }

    /* y position */
    /* We fill points in y by column (diagonal going from the top to the bottom left) */
    lev = 0;
    ipt = 0;
    xyz3d[max_npt * lev + ipt][1] = 0.0;

    for (int iquad = 0; iquad < nquad; iquad++) {
      lev++;

      /* There are "lev + 1" points at level "lev". Beware, C language is 0-based */
      int lev_npt = lev + 1;
      for (ipt = 0; ipt < lev_npt; ipt++)
        xyz3d[max_npt * lev + ipt][1] = ipt * (1.0 - xyz3d[max_npt * lev][2]) / lev;

    }

    /* x position (every points exist on the plane x+y+z=1) */

    for (int ii = 0; ii < nquad * max_npt; ii++)
      xyz3d[ii][0]  = 1.0 - xyz3d[ii][2] - xyz3d[ii][1];

    /* Compute the surface and the weight of each triangle     */

    int ntri = 0;

    cs_real_33_t posnod;

    /* Number of triangle levels */

    for (lev = 0; lev < nquad; lev++) {

      /* First triangle of the level */
      posnod[0][0] = xyz3d[max_npt * lev][0];
      posnod[0][1] = xyz3d[max_npt * lev][1];
      posnod[0][2] = xyz3d[max_npt * lev][2];
      posnod[1][0] = xyz3d[max_npt * (lev + 1)][0];
      posnod[1][1] = xyz3d[max_npt * (lev + 1)][1];
      posnod[1][2] = xyz3d[max_npt * (lev + 1)][2];
      posnod[2][0] = xyz3d[max_npt * (lev + 1) + 1][0];
      posnod[2][1] = xyz3d[max_npt * (lev + 1) + 1][1];
      posnod[2][2] = xyz3d[max_npt * (lev + 1) + 1][2];
      cs_rp->angsol[ntri] = _lhuilier(posnod);
      ntri++;

      /* Number of double triangle for this level */

      for (int ii = 0; ii < lev; ii++) {
        posnod[0][0] = xyz3d[max_npt * lev + ii][0];
        posnod[0][1] = xyz3d[max_npt * lev + ii][1];
        posnod[0][2] = xyz3d[max_npt * lev + ii][2];
        posnod[1][0] = xyz3d[max_npt * lev + ii + 1][0];
        posnod[1][1] = xyz3d[max_npt * lev + ii + 1][1];
        posnod[1][2] = xyz3d[max_npt * lev + ii + 1][2];
        posnod[2][0] = xyz3d[max_npt * (lev + 1) + ii + 1][0];
        posnod[2][1] = xyz3d[max_npt * (lev + 1) + ii + 1][1];
        posnod[2][2] = xyz3d[max_npt * (lev + 1) + ii + 1][2];
        cs_rp->angsol[ntri] = _lhuilier(posnod);
        ntri++;

        posnod[0][0] = xyz3d[max_npt * lev       + ii + 1][0];
        posnod[0][1] = xyz3d[max_npt * lev       + ii + 1][1];
        posnod[0][2] = xyz3d[max_npt * lev       + ii + 1][2];
        posnod[1][0] = xyz3d[max_npt * (lev + 1) + ii + 1][0];
        posnod[1][1] = xyz3d[max_npt * (lev + 1) + ii + 1][1];
        posnod[1][2] = xyz3d[max_npt * (lev + 1) + ii + 1][2];
        posnod[2][0] = xyz3d[max_npt * (lev + 1) + ii + 2][0];
        posnod[2][1] = xyz3d[max_npt * (lev + 1) + ii + 2][1];
        posnod[2][2] = xyz3d[max_npt * (lev + 1) + ii + 2][2];
        cs_rp->angsol[ntri] = _lhuilier(posnod);
        ntri++;
      }

    }

    BFT_FREE(xyz3d);

  }

  else if (cs_rp->i_quadrature == 7) {
    vec[0]    = 0.963560905;
    vec[1]    = 0.189143308;
    vec[2]    = 0.772965714;
    vec[3]    = 0.448622338;
    vec[4]    = 0.686596043;
    vec[5]    = 0.239106143;
    vec[6]    = 0.879538138;
    vec[7]    = 0.475828397;
    vec[8]    = 0.0;
    weight[0]  = 16 * atan(1.0) / 96.0;
    weight[1]  = 8  * atan(1.0) / 96.0;
    cs_rp->sxyz[0][0]  = vec[0];
    cs_rp->sxyz[1][0]  = vec[1];
    cs_rp->sxyz[2][0]  = vec[1];
    cs_rp->sxyz[3][0]  = vec[2];
    cs_rp->sxyz[4][0]  = vec[3];
    cs_rp->sxyz[5][0]  = vec[3];
    cs_rp->sxyz[6][0]  = vec[4];
    cs_rp->sxyz[7][0]  = vec[4];
    cs_rp->sxyz[8][0]  = vec[5];
    cs_rp->sxyz[9][0]  = vec[6];
    cs_rp->sxyz[10][0] = vec[7];
    cs_rp->sxyz[11][0] = vec[8];
    cs_rp->sxyz[12][0] = vec[8];
    cs_rp->sxyz[13][0] = vec[7];
    cs_rp->sxyz[14][0] = vec[6];
    cs_rp->sxyz[0][1]  = vec[1];
    cs_rp->sxyz[1][1]  = vec[0];
    cs_rp->sxyz[2][1]  = vec[1];
    cs_rp->sxyz[3][1]  = vec[3];
    cs_rp->sxyz[4][1]  = vec[2];
    cs_rp->sxyz[5][1]  = vec[2];
    cs_rp->sxyz[6][1]  = vec[5];
    cs_rp->sxyz[7][1]  = vec[4];
    cs_rp->sxyz[8][1]  = vec[4];
    cs_rp->sxyz[9][1]  = vec[7];
    cs_rp->sxyz[10][1] = vec[6];
    cs_rp->sxyz[11][1] = vec[6];
    cs_rp->sxyz[12][1] = vec[7];
    cs_rp->sxyz[13][1] = vec[8];
    cs_rp->sxyz[14][1] = vec[8];
    cs_rp->sxyz[0][2]  = vec[1];
    cs_rp->sxyz[1][2]  = vec[1];
    cs_rp->sxyz[2][2]  = vec[0];
    cs_rp->sxyz[3][2]  = vec[3];
    cs_rp->sxyz[4][2]  = vec[3];
    cs_rp->sxyz[5][2]  = vec[2];
    cs_rp->sxyz[6][2]  = vec[4];
    cs_rp->sxyz[7][2]  = vec[5];
    cs_rp->sxyz[8][2]  = vec[4];
    cs_rp->sxyz[9][2]  = vec[8];
    cs_rp->sxyz[10][2] = vec[8];
    cs_rp->sxyz[11][2] = vec[7];
    cs_rp->sxyz[12][2] = vec[6];
    cs_rp->sxyz[13][2] = vec[6];
    cs_rp->sxyz[14][2] = vec[7];
    for (int ii = 0; ii < 9; ii++)
      cs_rp->angsol[ii]   = weight[0];
    for (int ii = 9; ii < 15; ii++)
      cs_rp->angsol[ii]   = weight[1];
  }

  else if (cs_rp->i_quadrature == 8) {
    vec[0]   = 0.939848342;
    vec[1]   = 0.241542019;
    vec[2]   = 0.6817799;
    vec[3]   = 0.265240148;
    weight[0] = 0.2437531319;
    weight[1] = 0.2798456437;

    cs_rp->sxyz[0][0] = vec[0];
    cs_rp->sxyz[1][0] = vec[2];
    cs_rp->sxyz[2][0] = vec[2];
    cs_rp->sxyz[3][0] = vec[1];
    cs_rp->sxyz[4][0] = vec[3];
    cs_rp->sxyz[5][0] = vec[1];
    cs_rp->sxyz[0][1] = vec[1];
    cs_rp->sxyz[1][1] = vec[3];
    cs_rp->sxyz[2][1] = vec[2];
    cs_rp->sxyz[3][1] = vec[1];
    cs_rp->sxyz[4][1] = vec[2];
    cs_rp->sxyz[5][1] = vec[0];
    cs_rp->sxyz[0][2] = vec[1];
    cs_rp->sxyz[1][2] = vec[2];
    cs_rp->sxyz[2][2] = vec[3];
    cs_rp->sxyz[3][2] = vec[0];
    cs_rp->sxyz[4][2] = vec[2];
    cs_rp->sxyz[5][2] = vec[1];
    cs_rp->angsol[0] = weight[0];
    cs_rp->angsol[1] = weight[1];
    cs_rp->angsol[2] = weight[1];
    cs_rp->angsol[3] = weight[0];
    cs_rp->angsol[4] = weight[1];
    cs_rp->angsol[5] = weight[0];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
