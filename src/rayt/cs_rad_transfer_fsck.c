/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

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
#include "cs_time_step.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_fsck.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/* \file cs_rad_transfer_fsck.c.f90 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local static variables
 *============================================================================*/

static int ipass = 0;

static cs_real_t  *gi;
static cs_real_t  *tt;
static cs_real_t  *kpco2;
static cs_real_t  *kph2o;
static cs_real_t  *wv;
static cs_real_t  *dwv;
static cs_real_t  *kmfs;
static cs_real_t  *gq;

/*=============================================================================
 * Local const variables
 *============================================================================*/

const int ng = 100;
const int nt     = 38;
const int nconc  = 5;
const int nband  = 450;
const int maxit  = 1000;
const int imlinear    = 0;
const int im3dspline  = 7;
const int im1dspline  = 2;
const int im2dspline  = 6;
const cs_real_t eps    = 3e-14;
const cs_real_t x_kg[5] = {0.0001, 0.25, 0.5, 0.75, 1.0};

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Convert a line of values in a file into an array of cs_real_t and
 *        returns the number of values read
 *
 * \param[in]  radfile     pointer to file for reading values
 * \param[in]  values      array for values storage
 * \param[in]  nvalues     number of values read
 */
/*----------------------------------------------------------------------------*/

static inline void
_line_to_array(FILE      *radfile,
               cs_real_t  values[],
               cs_int_t  *nvalues)
{
  char line[256];
  fgets(line, 256, radfile);
  int index = 0;

  /* now read all info along the line */
  while (strlen(line) > 1) {
    char temp[256];
    /* store next value in string format */
    sscanf(line, "%s", temp);
    /* read and convert to double precision */
    sscanf(temp, "%lf", &(values[index]));
    /* increase index in array */
    index++;

    int l = strlen(temp);
    int i = 0;
    while (line[i] == ' ')
      i++;
    snprintf(temp, 256, "%s", &line[l+i]);
    strcpy(line, temp);
  }

  *nvalues = index;

  return;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief number of segents
 *
 * \param[in]     trad          Reference temperature
 * \param[in]     t             Local Temperature
 * \param[in]     xco2          CO2 volume fraction
 * \param[in]     xh2o          H2O volume fraction
 * \param[in]     interp_method Interpolation method
 * \param[inout]  itx[4][4]     itx[0]: TPlanck; iTx[1]: Tloc;
 *                              iTx[2]: xCO2; iTx[3] = xH2O
 */
/*----------------------------------------------------------------------------*/

inline static void
_gridposnbsg1(cs_real_t  trad,
              cs_real_t  t,
              cs_real_t  xco2,
              cs_real_t  xh2o,
              int        interp_method,
              int        itx[4][4])
{

  int itrad[4] = {0};
  int ita[4]   = {0};
  int ico2a[4] = {0};
  int ih2oa[4] = {0};

  /* Mole fraction interpolation: determine xCO2 */

  int i = 0;
  int j = nconc - 1;
  int ip = (i + j) / 2;
  while (j - i > 1) {
    if (xco2 < x_kg[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  /*  spline on x */

  if (interp_method == im2dspline) {

    if ((i > 0) && (i < nconc - 2)) {
      ico2a[0] = i - 1;
      ico2a[1] = i;
      ico2a[2] = i + 1;
      ico2a[3] = i + 2;
    }
    else if (i == 0) {
      ico2a[0] = 0;
      ico2a[1] = 1;
      ico2a[2] = 2;
      ico2a[3] = 3;
    }
    else if (i == nconc - 2) {
      ico2a[0] = nconc - 4;
      ico2a[1] = nconc - 3;
      ico2a[2] = nconc - 2;
      ico2a[3] = nconc - 1;
    }
    else
      cs_log_printf(CS_LOG_DEFAULT, _("x grid failure, i = %d"), i);
  }
  /* linear in x */
  else {
    if (i < 0) //not possible...
      i = 0;
    else if (i > nconc - 2)
      i = nconc - 2;
    ico2a[0] = i;
    ico2a[1] = i + 1;
  }

  /* Mole fraction interpolation: determine xH2O */

  i  = 0;
  j  = nconc - 1;
  ip = (i + j) / 2;
  while (j - i > 1) {
    if (xh2o < x_kg[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  /*  spline on x */
  if (interp_method == im2dspline) {

    if ((i > 0) && (i < nconc - 2)) {
      ih2oa[0] = i - 1;
      ih2oa[1] = i;
      ih2oa[2] = i + 1;
      ih2oa[3] = i + 2;
    }
    else if (i == 0) {
      ih2oa[0] = 0;
      ih2oa[1] = 1;
      ih2oa[2] = 2;
      ih2oa[3] = 3;
    }
    else if (i == nconc - 2) {
      ih2oa[0] = nconc - 4;
      ih2oa[1] = nconc - 3;
      ih2oa[2] = nconc - 2;
      ih2oa[3] = nconc - 1;
    }
    else
      cs_log_printf(CS_LOG_DEFAULT, _("x grid failure, i = %d"), i);
  }
  /* linear in x */
  else {
    if (i < 0) //not possible...
      i = 0;
    else if (i > nconc - 2)
      i = nconc - 2;
    ih2oa[0] = i;
    ih2oa[1] = i + 1;
  }

  /* Temperature interpolation */

  i  = 0;
  j  = nt - 1;
  ip = (i + j) / 2;
  while (j - i > 1) {
    if (t < tt[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  /* spline in t */
  if (interp_method != 0) {

    if ((i > 0) && (i < nt - 2)) {
      ita[0] = i - 1;
      ita[1] = i;
      ita[2] = i + 1;
      ita[3] = i + 2;
    }
    else if (i == 0) {
      ita[0]  = 0;
      ita[1]  = 1;
      ita[2]  = 2;
      ita[3]  = 3;
    }
    else if (i == nt - 1) {
      ita[0] = nt - 4;
      ita[1] = nt - 3;
      ita[2] = nt - 2;
      ita[3] = nt - 1;
    }
    else
      cs_log_printf(CS_LOG_DEFAULT, _("t grid failure, i = %d"), i);
  }
  /* linear in t */
  else {
    if (i < 0)
      i  = 0;
    else if (i > nt - 2)
      i = nt - 2;
    ita[0] = i;
    ita[1] = i + 1;
  }

  /* Determine closest Planck temperature in database
     lower than local temperature */

  i  = 0;
  j  = nt - 1;
  ip = (i + j) / 2;
  while (j - i > 1) {
    if (trad < tt[ip])
      j = ip;
    else
      i = ip;
    ip = (i + j) / 2;
  }

  /* spline in t */
  if (interp_method != 0) {
    if ((i > 0) && (i < nt - 2)) {
      itrad[0] = i - 1;
      itrad[1] = i;
      itrad[2] = i + 1;
      itrad[3] = i + 2;
    }
    else if (i == 0) {
      itrad[0] = 0;
      itrad[1] = 1;
      itrad[2] = 2;
      itrad[3] = 3;
    }
    else if (i == nt - 2) {
      itrad[0] = nt - 4;
      itrad[1] = nt - 3;
      itrad[2] = nt - 2;
      itrad[3] = nt - 1;
    }
    else
      cs_log_printf(CS_LOG_DEFAULT, _("t grid failure, i = %d"), i);
  }
  /* linear in t */
  else {
    if (i < 0)
        i = 0;
    else if (i > nt - 1) {
      i = nt - 1;
    }
    itrad[0] = i;
    itrad[1] = i + 1;
  }

  /* attribution itx */

  for (int k = 0; k < 4; k++) {
    itx[0][k] = itrad[k];
    itx[1][k] = ita[k];
    itx[2][k] = ico2a[k];
    itx[3][k] = ih2oa[k];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief This subroutine evaluates the cubic spline function
 *
 *    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
 *
 *    where  x(i) .lt. u .lt. x(i+1), using horner's rule
 *
 *  if  u .lt. x(1) then  i = 1  is used.
 *  if  u .ge. x(n) then  i = n  is used.
 *
 *  if  u  is not in the same interval as the previous call, then a
 *  binary search is performed to determine the proper interval.
 *
 * \param[in]     n            Number of data points
 * \param[in]     u            Abscissa at which the spline is to be evaluated
 * \param[in]     x,y          Arrays of data abscissas and ordinates
 * \param[in]     b,c,d        Arrays of spline coefficients computed by spline
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
_seval(int       n,
       cs_real_t u,
       cs_real_t x[],
       cs_real_t y[],
       cs_real_t b[],
       cs_real_t c[],
       cs_real_t d[])
{
  int i = 0;
  if (i > n)
    i = 0;

  /*  Binary search */
  if (   u < x[i]
      || u > x[i + 1]) {
    i = 0;
    int j = n;
    int k = (j + i)/2;
    while (j > i + 1) {
      if (u < x[k])
        j = k;
      else
        i = k;
    }
  }
  /*  Evaluate spline */
  cs_real_t dx = u - x[i];

  return y[i] + dx * (b[i] + dx * (c[i] + dx * d[i]));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  1-d monotonic spline interpolation: coefficients
 *  The coefficients b(i), c(i), and d(i), i = 1,2,...,n are computed
 *  for a monotonically varying cubic interpolating spline
 *
 *    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
 *    for  x(i) .le. x .le. x(i+1)
 *    with y(i+1).ge.y(i) (all i) or y(i+1).lt.y(i) (all i)
 *
 ********************************************************************************
 *
 * THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
 * BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
 *
 ********************************************************************************
 *
 *    Y(I) = S(X(I))
 *    B(I) = SP(X(I))
 *    C(I) = SPP(X(I))/2
 *    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
 *
 *  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
 *  TO EVALUATE THE SPLINE.
 *
 * \param[in]     n           Number of data points or knots (n.ge.2)
 * \param[in]     x           Abscissa of the knots in strictly increasing order
 * \param[in]     y           Ordinates of the knots
 * \param[out]    b,c,d       Arrays of spline coefficients as defined above
 */
/*----------------------------------------------------------------------------*/

inline static void
_splmi(int        n,
       cs_real_t  x[],
       cs_real_t  y[],
       cs_real_t  b[],
       cs_real_t  c[],
       cs_real_t  d[])
{
  cs_real_t delta[n], h[n], al[n], be[n];

  int nm1 = n - 1;

  if (n < 2)
    return;

  if (n >= 3) {

    /*  Calculate the h(i) and delta(i) */
    for (int i = 0; i < nm1; i++) {
      al[i]    = 0.0;
      be[i]    = 0.0;
      h[i]     = x[i + 1] - x[i];
      delta[i] = (y[i + 1] - y[i]) / h[i];
    }

    /* calculate first values for al and be by 3-point difference */
    if (CS_ABS(delta[0]) < eps)
      al[0] =  (  pow((h[0] + h[1]), 2.0) * y[1]
                - pow (h[0], 2.0) * y[2]
                - h[1] * (2.0 * h[0] + h[1]) * y[0])
             / (h[1] * (h[0] + h[1]) * (y[1] - y[0]));

    for (int i = 1; i < nm1; i++) {
      if (CS_ABS(delta[0]) < eps)
        al[i] =  (  pow(h[i - 1], 2.0) * y[i + 1]
                  + (pow (h[i], 2.0) - pow (h[i - 1], 2.0)) * y[i]
                  - pow (h[i], 2.0) * y[i - 1])
               / (h[i - 1] * (h[i] + h[i - 1]) * (y[i + 1] - y[i]));

    }

    int nm2 = n - 2;
    for (int i = 0; i < nm2; i++) {
      if (CS_ABS(delta[0]) < eps)
        be[i] =  (  pow(h[i], 2.0) * y[i + 2]
                  + (pow (h[i + 1], 2.0) - pow (h[i], 2.0)) * y[i + 1]
                  - pow (h[i + 1], 2.0) * y[i])
               / (h[i + 1] * (h[i] + h[i + 1]) * (y[i + 1] - y[i]));

    }

    if (CS_ABS(delta[0]) < eps) {
      be[n - 2] =  (  h[n - 3] * (2.0 * h[n - 2] + h[n - 3]) * y[n - 1]
                     - pow((h[n - 2] + h[n - 3]), 2.0) * y[n - 2]
                     + pow (h[n - 2], 2.0) * y[n - 3])
                  / (h[n - 3] * (h[n - 2] + h[n - 3]) * (y[n - 1] - y[n - 2]));
    }

    /* Correct values for al and be */
    for (int i = 0; i < nm1; i++) {
      if (   al[i] + be[i] > 2.0
          && 2.0 * al[i] + be[i] > 3.0
          && al[i] + 2.0 * be[i] > 3.0) {

        cs_real_t phi = al[i] - pow((2.0 * al[i] + be[i] - 3.0), 2.0) / (al[i] + be[i] - 2.0) / 3.0;

        if (phi < 0.0) {
          cs_real_t ti = 3.0 / sqrt (pow (al[i], 2.0) + pow (be[i], 2.0));
          al[i] = ti * al[i];
          be[i] = ti * be[i];
        }

      }
    }

    /* Calculate spline coefficients */

    for (int i = 0; i < nm1; i++) {
      d[i] = (al[i] + be[i] - 2.0) * delta[i] / (pow (h[i], 2.0));
      c[i] = (3.0 - 2.0 * al[i] - be[i]) * delta[i] / (h[i]);
      b[i] = al[i] * delta[i];
      if (b[i] * delta[i] < 0.0) {
        b[i] = delta[i];
        c[i] = 0.0;
        d[i] = 0.0;
      }
    }

    b[n - 1] = be[n - 2] * delta[n - 2];
    c[n - 1] = 0.0;
    d[n - 1] = 0.0;

  }
  else if (n == 2) {

    b[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c[0] = 0.0;
    d[0] = 0.0;
    b[1] = b[0];
    c[1] = 0.0;
    d[1] = 0.0;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief simple_interpg
 *
 * \param[in]     nxy
 * \param[in]     xx
 * \param[in]     yy
 * \param[in]     ni
 * \param[in]     xi
 * \param[out]    yi
 */
/*----------------------------------------------------------------------------*/

static inline void
_simple_interpg (int        nxy,
                 cs_real_t  xx[],
                 cs_real_t  yy[],
                 int        ni,
                 cs_real_t  xi[],
                 cs_real_t  yi[])
{
  int ibgn = 0;

  for (int iq = 0; iq < ni; iq++) {

    for (int i = ibgn; i < nxy; i++) {

      if (xi[iq] < xx[0]) {
        yi[iq] = yy[0] * xi[iq] / CS_MAX(1e-09, xx[0]);
        break;
      }
      else if (xi[iq] > xx[nxy])
        yi[iq] = yy[nxy];


      /* interpolate */
      if (CS_ABS(xi[iq] - xx[i]) / (xx[i] + 1e-15) < 0.001) {
        yi[iq] = yy[i];
        ibgn = i;
        break;
      }
      else if (xx[i] > xi[iq]) {
        yi[iq] =   yy[i - 1]
                + (yy[i] - yy[i - 1]) * (xi[iq] - xx[i - 1])
                                      / CS_MAX(1e-09, (xx[i] - xx[i - 1]));
        ibgn = i;
        break;
      }
    }
  }

  if (CS_ABS(xi[ni] - xx[nxy]) / (xx[nxy] + 1e-15) < 0.001)
    yi[ni] = yy[nxy];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief interpolation4d
 *
 * \param[in]     trad          Reference temperature
 * \param[in]     t             Local temperature
 * \param[in]     xco2          Reference CO2 volume fraction
 * \param[in]     xh2o          Reference H2O volume fraction
 * \param[in]     interp_method Interpolation method
 * \param[out]    gdb
 * \param[out]    kdb
 */
/*----------------------------------------------------------------------------*/

inline static void
_interpolation4d(cs_real_t trad,
                 cs_real_t t,
                 cs_real_t xco2,
                 cs_real_t xh2o,
                 int       interp_method,
                 cs_real_t gdb[],
                 cs_real_t kdb[])
{
  int     itx[4][4];
  int     nix, nit;

  cs_real_t *karray, *kint1, *kint2, *kint3;
  BFT_MALLOC(karray, ng*4*4*4*4, cs_real_t);
  BFT_MALLOC(kint1,  ng*4*4*4, cs_real_t);
  BFT_MALLOC(kint2,  ng*4*4, cs_real_t);
  BFT_MALLOC(kint3,  ng*4, cs_real_t);

  cs_real_t *b, *c, *d, *kg_t2, *kg_x2;
  BFT_MALLOC(b, 4, cs_real_t);
  BFT_MALLOC(c, 4, cs_real_t);
  BFT_MALLOC(d, 4, cs_real_t);
  BFT_MALLOC(kg_t2, 4, cs_real_t);
  BFT_MALLOC(kg_x2, 4, cs_real_t);

  /* Determine positions in x and T
   * in the NB database for interpolation. */
  for (int i = 0; i < ng; i++)
    gdb[i] = gi[i];

  /* number of interpolation points along t & x: 2 linear or 4 spline */
  _gridposnbsg1(trad,
                t,
                xco2,
                xh2o,
                interp_method,
                itx);

  /* spline over x */
  if (interp_method == im2dspline)
    nix  = 4;
  else
    nix  = 2;

  /* spline over t */
  if (interp_method != 0)
    nit  = 4;
  else
    nit  = 2;

  /* Attribute interpolation point indexes along T & x */

  int itrada[4] = {0, 0, 0, 0};
  int ita[4]    = {0, 0, 0, 0};
  int ico2a[4] = {0, 0, 0, 0};
  int ih2oa[4] = {0, 0, 0, 0};

  for (int i = 0; i < nit; i++) {
    itrada[i] = itx[0][i];
    ita[i] = itx[1][i];
  }

  for (int i = 0; i < nix; i++) {
    ico2a[i] = itx[2][i];
    ih2oa[i] = itx[3][i];
  }

  for (int ih2o = 0; ih2o < nix; ih2o++) {
    for (int ico2 = 0; ico2 < nix; ico2++) {
      for (int it = 0; it < nit; it++) {
        for (int itrad = 0; itrad < nit; itrad++) {
          for (int ig = 0; ig < ng; ig++) {
            karray[  ih2o
                   + ico2*4
                   + it*4*4
                   + itrad*4*4*4
                   + ig*4*4*4*4] =
              kmfs[  ih2oa[ih2o]
                   + ico2a[ico2]*nconc
                   + ita[it]*nconc*nconc
                   + itrada[itrad]*nconc*nconc*nt
                   + ig*nconc*nconc*nt*nt];
          }
        }
      }
    }
  }

  /* Interpolation on XH2O; spline over x */

  if (interp_method == im2dspline) {
    for (int ico2 = 0; ico2 < nix; ico2++) {
      for (int it = 0; it < nit; it++) {
        for (int itrad = 0; itrad < nit; itrad++) {
          for (int ig = 0; ig < ng; ig++) {
            kg_x2[0] = x_kg[ih2oa[0]];
            kg_x2[1] = x_kg[ih2oa[1]];
            kg_x2[2] = x_kg[ih2oa[2]];
            kg_x2[3] = x_kg[ih2oa[3]];
            int ih2o = 0;
            int nargs = 4;
            _splmi(nargs,
                   kg_x2,
                   &karray[ih2o + ico2*4 + it*4*4 + itrad*4*4*4 + ig*4*4*4*4],
                   b,
                   c,
                   d);

            kint1[ico2 + it*4 + itrad*4*4 + ig*4*4*4]
              = _seval(nargs,
                       xh2o,
                       kg_x2,
                       &karray[ih2o + ico2*4 + it*4*4 + itrad*4*4*4 + ig*4*4*4*4],
                       b,
                       c,
                       d);
          }
        }
      }
    }
  }
  else {
    cs_real_t wx = (xh2o - x_kg[ih2oa[0]]) / (x_kg[ih2oa[1]] - x_kg[ih2oa[0]]);
    for (int ico2 = 0; ico2 < nix; ico2++) {
      for (int it = 0; it < nit; it++) {
        for (int itrad = 0; itrad < nit; itrad++) {
          for (int ig = 0; ig < ng; ig++) {
            kint1[ico2 + it*4 + itrad*4*4 + ig*4*4*4]
              =          wx  * karray[1 + ico2*4 + it*4*4 + itrad*4*4*4 + ig*4*4*4*4]
                + (1.0 - wx) * karray[0 + ico2*4 + it*4*4 + itrad*4*4*4 + ig*4*4*4*4];
          }
        }
      }
    }
  }

  /* Interpolation on XCO2: spline over x */

  if (interp_method == im2dspline) {
    for (int it = 0; it < nit; it++) {
      for (int itrad = 0; itrad < nit; itrad++) {
        for (int ig = 0; ig < ng; ig++) {
          kg_x2[0] = x_kg[ico2a[0]];
          kg_x2[1] = x_kg[ico2a[1]];
          kg_x2[2] = x_kg[ico2a[2]];
          kg_x2[3] = x_kg[ico2a[3]];

          int ico2 = 0;
          int nargs = 4;
          _splmi(nargs,
                 kg_x2,
                 &kint1[ico2 + it*4 + itrad*4*4 + ig*4*4*4],
                 b,
                 c,
                 d);
          kint2[it + itrad*4 + ig*4*4]
            = _seval(nargs,
                     xco2,
                     kg_x2,
                     &kint1[ico2 + it*4 + itrad*4*4 + ig*4*4*4],
                     b,
                     c,
                     d);
        }
      }
    }
  }
  else {
    cs_real_t wx = (xco2 - x_kg[ico2a[0]]) / (x_kg[ico2a[1]] - x_kg[ico2a[0]]);
    for (int it = 0; it < nit; it++) {
      for (int itrad = 0; itrad < nit; itrad++) {
        for (int ig = 0; ig < ng; ig++) {
          kint2[it + itrad*4 + ig*4*4]
            =          wx  * kint1[1 + it*4 + itrad*4*4 + ig*4*4*4]
              + (1.0 - wx) * kint1[0 + it*4 + itrad*4*4 + ig*4*4*4];
        }
      }
    }
  }

  /* Interpolation on T: spline on t */

  if (interp_method != 0) {
    for (int itrad = 0; itrad < nit; itrad++) {
      for (int ig = 0; ig < ng; ig++) {
        kg_t2[0] = tt[ita[0]];
        kg_t2[1] = tt[ita[1]];
        kg_t2[2] = tt[ita[2]];
        kg_t2[3] = tt[ita[3]];
        int nargs = 4;
        _splmi (nargs,
                kg_t2,
                &kint2[0 + itrad*4 + ig*4*4],
                b,
                c,
                d);
        kint3[itrad + ig*4] = _seval(nargs,
                                     t,
                                     kg_t2,
                                     &kint2[0 + itrad*4 + ig*4*4],
                                     b,
                                     c,
                                     d);
      }
    }
  }
  else {
    cs_real_t wt = (t - tt[ita[0]]) / (tt[ita[1]] - tt[ita[0]]);
    for (int itrad = 0; itrad < nit; itrad++) {
      for (int ig = 0; ig < ng; ig++) {
        kint3[itrad + ig*4] =         wt  * kint2[1 + itrad*4 + ig*4*4]
                             + (1.0 - wt) * kint2[0 + itrad*4 + ig*4*4];
      }
    }
  }

  /* Interpolation on Trad: spline */

  if (interp_method != 0) {
    for (int ig = 0; ig < ng; ig++) {
      kg_t2[0] = tt[itrada[0]];
      kg_t2[1] = tt[itrada[1]];
      kg_t2[2] = tt[itrada[2]];
      kg_t2[3] = tt[itrada[3]];
      int nargs = 4;
      _splmi(nargs,
             kg_t2,
             &kint3[0 + ig*4],
             b,
             c,
             d);
      kdb[ig] = _seval(nargs,
                       trad,
                       kg_t2,
                       &kint3[0 + ig*4],
                       b,
                       c,
                       d);
    }
  }
  else {
    cs_real_t wt = (trad - tt[itrada[0]]) / (tt[itrada[1]] - tt[itrada[0]]);
    for (int ig = 0; ig < ng; ig++) {
      kdb[ig] = wt * kint3[1 + ig*4] + (1.0 - wt) * kint3[0 + ig*4];
    }
  }

  /* Free memory */
  BFT_FREE(karray);
  BFT_FREE(kint1);
  BFT_FREE(kint2);
  BFT_FREE(kint3);
  BFT_FREE(b);
  BFT_FREE(c);
  BFT_FREE(d);
  BFT_FREE(kg_t2);
  BFT_FREE(kg_x2);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Determine the radiation coefficients of the FSCK model
 *         as well as the corresponding weights.
 *
 * \param[in]     pco2        CO2 volume fraction
 * \param[in]     ph2o        H2O volume fraction
 * \param[in]     teloc       gas temperature
 * \param[out]    kloc        radiation coefficient of the i different gases
 * \param[out]    aloc        weights of the i different gases in cells
 * \param[out]    alocb       weights of the i different gases at boundaries
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_fsck(const cs_real_t  *restrict pco2,
                     const cs_real_t  *restrict ph2o,
                     const cs_real_t  *restrict teloc,
                     cs_real_t        *restrict kloc,
                     cs_real_t        *restrict aloc,
                     cs_real_t        *restrict alocb)
{
  /* Initialization */

  cs_real_t   *gfskref, *kfskref, *kfsk, *gfsk, *gg1;
  cs_real_t   *kg1, *as, *ag, *aw, *kloctmp;

  BFT_MALLOC(gfskref, ng, cs_real_t);
  BFT_MALLOC(kfskref, ng, cs_real_t);
  BFT_MALLOC(gfsk, ng, cs_real_t);
  BFT_MALLOC(kfsk, ng, cs_real_t);
  BFT_MALLOC(gg1, ng, cs_real_t);
  BFT_MALLOC(kg1, ng, cs_real_t);
  BFT_MALLOC(as, ng, cs_real_t);

  BFT_MALLOC(ag, cs_glob_rad_transfer_params->nwsgg, cs_real_t);
  BFT_MALLOC(aw, cs_glob_rad_transfer_params->nwsgg, cs_real_t);
  BFT_MALLOC(kloctmp, cs_glob_rad_transfer_params->nwsgg, cs_real_t);

  cs_field_t *f_bound_t = cs_field_by_name_try("boundary_temperature");
  cs_real_t *tpfsck = f_bound_t->val;

  /* Read the data base files */

  ipass++;
  if (ipass == 1) { /* Read parameters files */

    FILE *radfile = NULL;
    const char *pathdatadir = cs_base_get_pkgdatadir();
    char filepath[256];

    BFT_MALLOC(gi,    ng, cs_real_t);
    BFT_MALLOC(tt,    nt, cs_real_t);
    BFT_MALLOC(kpco2, nt, cs_real_t);
    BFT_MALLOC(kph2o, nt, cs_real_t);
    BFT_MALLOC(wv,    nband, cs_real_t);
    BFT_MALLOC(dwv,   nband, cs_real_t);

    BFT_MALLOC(kmfs, nconc * nconc * nt *nt *ng, cs_real_t);

    /* Read k-distributions */
    {
      snprintf(filepath, 256, "%s/data/thch/dp_radiat_MFS", pathdatadir);
      radfile = fopen(filepath, "r");
      char line[256];
      for (int cco2 = 0; cco2 < nconc; cco2++) {
        for (int ch2o = 0; ch2o < nconc; ch2o++) {
          for (int it = 0; it < nt; it++) {
            for (int itrad = 0; itrad < nt; itrad++) {
              fgets(line, 256, radfile);
              fgets(line, 256, radfile);
              for (int ig = 0; ig < ng; ig++) {
                cs_real_t temp[2] = {0.};
                int nvalues;
                _line_to_array(radfile, temp, &nvalues);
                assert(nvalues == 2);
                gi[ig] = temp[0];
                kmfs[  ch2o
                     + cco2*nconc
                     + it*nconc*nconc
                     + itrad*nconc*nconc*nt
                     + ig*nconc*nconc*nt*nt] = temp[1];
              }
            }
          }
        }
      }
      fclose(radfile);
    }

    /* Read the Planck coefficients */
    {
      snprintf(filepath, 256, "%s/data/thch/dp_radiat_Planck_CO2", pathdatadir);
      radfile = fopen(filepath, "r");
      for (int it = 0; it < nt; it++) {
        cs_real_t temp[2] = {0.};
        int nvalues;
        _line_to_array(radfile, temp, &nvalues);
        assert(nvalues == 2);
        tt[it]    = temp[0];
        kpco2[it] = temp[1];
      }
      fclose(radfile);
    }

    {
      snprintf(filepath, 256, "%s/data/thch/dp_radiat_Planck_H2O", pathdatadir);
      radfile = fopen(filepath, "r");
      for (int it = 0; it < nt; it++) {
        cs_real_t temp[2] = {0.};
        int nvalues;
        _line_to_array(radfile, temp, &nvalues);
        assert(nvalues == 2);
        tt[it]    = temp[0];
        kph2o[it] = temp[1];
      }
      fclose(radfile);
    }

    /* Read the wavelength interval */
    {
      snprintf(filepath, 256, "%s/data/thch/dp_radiat_wave", pathdatadir);
      radfile = fopen(filepath, "r");
      for (int iwvnb = 0; iwvnb < nband; iwvnb++) {
        cs_real_t temp[2] = {0.};
        int nvalues;
        _line_to_array(radfile, temp, &nvalues);
        assert(nvalues == 2);
        wv[iwvnb]  = temp[0];
        dwv[iwvnb] = temp[1];
      }
      fclose(radfile);
    }

    /* Gaussian quadrature */

    /* Allocation */
    BFT_MALLOC(gq, cs_glob_rad_transfer_params->nwsgg, cs_real_t);
    int m = (cs_glob_rad_transfer_params->nwsgg + 1) / 2;
    cs_real_t *p1, *p2, *p3, *pp;
    cs_real_t *z, *z1, *arth;
    bool *unfinished;
    BFT_MALLOC(p1, m, cs_real_t);
    BFT_MALLOC(p2, m, cs_real_t);
    BFT_MALLOC(p3, m, cs_real_t);
    BFT_MALLOC(pp, m, cs_real_t);
    BFT_MALLOC(z,  m, cs_real_t);
    BFT_MALLOC(z1, m, cs_real_t);
    BFT_MALLOC(arth, m, cs_real_t);
    BFT_MALLOC(unfinished, m, bool);

    /* Boundaries */
    cs_real_t x1 = 0.0;
    cs_real_t x2 = 1.0;

    /* Given the lower and upper limits of integration x1 and x2, this routine
     * returns arrays x[0..n-1]
     * and w[0..n-1] of length n, containing the abscissas and weights of the
     * Gauss-Legendre n-point
     * quadrature formula. The parameter EPS is the relative precision.
     * Note that internal computations are done in double precision. */

    int n = cs_glob_rad_transfer_params->nwsgg;

    /* The roots are symmetric in the interval, so we only have to find half of them. */

    cs_real_t xm = 0.5 * (x2 + x1);
    cs_real_t xl = 0.5 * (x2 - x1);

    for (int i = 0; i < m; i++)
      arth[i] = i;

    /* initial aproximation of the roots */

    for (int i = 0; i < m; i++) {
      z[i] = cos(cs_math_pi * (arth[i] - 0.25) / (n + 0.5));
      unfinished[i] = true;
    }

    /* Newton's method carried out simultaneously on the roots */
    int its;
    for (its = 0; its < maxit; its++) {

      /* initialisation */

      for (int j = 0; j < m; j++) {
        if (unfinished[j]) {
          p1[j]    = 1.0;
          p2[j]    = 0.0;
        }
      }

      for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
          p3[i] = p2[i];
          p2[i] = p1[i];
          p1[i] = ((2.0 * j - 1.0) * z[i] * p2[i] - (j - 1.0) * p3[i]) / j;
        }
      }

      /* p1 now contains the desired Legendre polynomials.
       * We next compute pp, its derivative, by a standard relation involving
       * also p2, the polynomial of one lower order. */

      for (int j = 0; j < m; j++) {
        pp[j] = n * (z[j] * p1[j] - p2[j]) / (z[j] * z[j] - 1.0);
        z1[j] = z[j];
        z[j]  = z1[j] - p1[j] / pp[j];
        unfinished[j] = (CS_ABS(z[j] - z1[j]) > eps);
      }

      /* check if we should continue */
      int cont = 0;
      for (int j = 0; j < m; j++) {
        cont = 1;
      }
      if (cont == 0)
        break;

    }

    if (its == maxit)
      cs_log_printf(CS_LOG_DEFAULT, "Maximum number of iterations during GAULEG");

    /* Scale the root to the desired interval, and put in its symmetric counterpart. */

    for (int j = 0; j < m; j++) {
      gq[j]         = xm - xl * z[j];
      gq[n - j + 1] = xm + xl * z[j];
    }

    /* Compute the weight and its symmetric counterpart. */

    for (int j = 0; j < m; j++) {
      cs_glob_rad_transfer_params->wq[j] = 2.0 * xl / ((1.0 - pow (z[j], 2.0)) * pow (pp[j], 2.0));
      cs_glob_rad_transfer_params->wq[n - j + 1] = cs_glob_rad_transfer_params->wq[j];
    }

    BFT_FREE(p1);
    BFT_FREE(p2);
    BFT_FREE(p3);
    BFT_FREE(pp);
    BFT_FREE(z);
    BFT_FREE(z1);
    BFT_FREE(arth);
    BFT_FREE(unfinished);

  }

  /* Compute the reference state */

  cs_real_t pco2ref = 0.0;
  cs_real_t ph2oref = 0.0;
  cs_real_t tref   = 0.0;
  cs_real_t sum1   = 0.0;
  cs_real_t sum2   = 0.0;

  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {

    /* Calculation of pco2ref and ph2oref */
    pco2ref += pco2[iel] * cs_glob_mesh_quantities->cell_vol[iel];
    ph2oref += ph2o[iel] * cs_glob_mesh_quantities->cell_vol[iel];

    /* Calculation of tref */
    /* Interplolation of Planck coefficient for mix */
    cs_real_t kp;
    if (teloc[iel] <= tt[0])
      kp = pco2[iel] * kpco2[0] + ph2o[iel] * kph2o[0];
    else if (teloc[iel] >= tt[nt - 1])
      kp = pco2[iel] * kpco2[nt - 1] + ph2o[iel] * kpco2[nt - 1];
    else {
      kp = 0.;
      for (int it = 0; it < nt - 1; it++) {
        if ((teloc[iel] >= tt[it]) && (teloc[iel] < tt[it + 1])) {
          cs_real_t kp1 = pco2[iel] * kpco2[it] + ph2o[iel] * kph2o[it];
          cs_real_t kp2 = pco2[iel] * kpco2[it + 1] + ph2o[iel] * kph2o[it + 1];
          kp = (kp2 - kp1) / (tt[it + 1] - tt[it]) * (teloc[iel] - tt[it]) + kp1;
          break;
        }
      }
    }
    cs_real_t kpt4dv =  kp
                      * pow(teloc[iel], 4.0)
                      * cs_glob_mesh_quantities->cell_vol[iel];
    sum1 += kpt4dv * teloc[iel];
    sum2 += kpt4dv;
  }

  if (cs_glob_rank_id >= 0) {
    cs_parall_sum(1, CS_DOUBLE, &pco2ref);
    cs_parall_sum(1, CS_DOUBLE, &ph2oref);
    cs_parall_sum(1, CS_DOUBLE, &sum1);
    cs_parall_sum(1, CS_DOUBLE, &sum2);
  }

  pco2ref /= cs_glob_mesh_quantities->tot_vol;
  ph2oref /= cs_glob_mesh_quantities->tot_vol;
  tref = sum1 / sum2;

  /* Interpolation */

  /* Determination of the k-distribution at the reference state */
  int interp_method = imlinear;
  _interpolation4d(tref,
                   tref,
                   pco2ref,
                   ph2oref,
                   interp_method,
                   gfskref,
                   kfskref);

  /* [m^-1] */
  for (int i = 0; i < ng; i++)
    kfskref[i] *= 100.0;

  for (cs_lnum_t iel = 0; iel < cs_glob_mesh->n_cells; iel++) {

    /* Determination of the local absorbtion coefficient */
    for (int i = 0; i < ng; i++) {
      kfsk[i] = 0.;
      gfsk[i] = 0.;
    }
    _interpolation4d (tref,
                      teloc[iel],
                      pco2[iel],
                      ph2o[iel],
                      interp_method,
                      gfsk,
                      kfsk);
    /* [m^-1] */
    for (int i = 0; i < ng; i++)
      kfsk[i] /= 100.0;
    _simple_interpg(ng,
                    gfsk,
                    kfsk,
                    cs_glob_rad_transfer_params->nwsgg,
                    gq,
                    kloctmp);
    for (int i = 0; i < cs_glob_rad_transfer_params->nwsgg; i++)
      kloc[i * cs_glob_mesh->n_cells + iel] = kloctmp[i];


    /* Determination of the local weights */
    for (int i = 0; i < ng; i++) {
      kfsk[i] = 0.;
      gfsk[i] = 0.;
    }
    _interpolation4d(teloc[iel],
                     tref,
                     pco2ref,
                     ph2oref,
                     interp_method,
                     gg1,
                     kg1);
    /* [m^-1] */
    for (int i = 0; i < ng; i++)
      kg1[i] *= 100.0;
    _simple_interpg(ng,
                    kg1,
                    gg1,
                    ng,
                    kfskref,
                    gfsk);
    as[0]  = (gfsk[1] - gfsk[0]) / (gfskref[1] - gfskref[0] + 1e-15);
    as[ng-1] = (gfsk[ng-1] - gfsk[ng - 2]) / (gfskref[ng-1] - gfskref[ng - 2] + 1e-15);
    for (int k = 1; k < ng - 1; k++)
      as[k] = (gfsk[k + 1] - gfsk[k - 1]) / (gfskref[k + 1] - gfskref[k - 1] + 1e-15);
    _simple_interpg(ng,
                    gfskref,
                    as,
                    cs_glob_rad_transfer_params->nwsgg,
                    gq,
                    ag);
    for (int i = 0; i < cs_glob_rad_transfer_params->nwsgg; i++)
      aloc[i * cs_glob_mesh->n_cells + iel] = ag[i];

  }

  for (cs_lnum_t ifac = 0; ifac < cs_glob_mesh->n_b_faces; ifac++) {
    _interpolation4d(tpfsck[ifac],
                     tref,
                     pco2ref,
                     ph2oref,
                     interp_method,
                     gg1,
                     kg1);
    for (int i = 0; i < ng; i++)
      kg1[i] *= 100.0;
    _simple_interpg(ng,
                    kg1,
                    gg1,
                    ng,
                    kfskref,
                    gfsk);
    as[0]  = (gfsk[1] - gfsk[0]) / (gfskref[1] - gfskref[0] + 1e-15);
    as[ng] = (gfsk[ng] - gfsk[ng - 1]) / (gfskref[ng] - gfskref[ng - 1] + 1e-15);
    for (int k = 1; k < ng - 1; k++)
      as[k] = (gfsk[k + 1] - gfsk[k - 1]) / (gfskref[k + 1] - gfskref[k - 1] + 1e-15);
    _simple_interpg(ng,
                    gfskref,
                    as,
                    cs_glob_rad_transfer_params->nwsgg,
                    gq,
                    aw);
    for (int i = 0; i < cs_glob_rad_transfer_params->nwsgg; i++)
      alocb[i * cs_glob_mesh->n_b_faces + ifac] = aw[i];
  }

  /* free memory */
  if (cs_glob_time_step->nt_cur == cs_glob_time_step->nt_max) {
    BFT_FREE(gi);
    BFT_FREE(tt);
    BFT_FREE(kpco2);
    BFT_FREE(kph2o);
    BFT_FREE(wv);
    BFT_FREE(dwv);
    BFT_FREE(kmfs);
    BFT_FREE(gq);
  }

  BFT_FREE(kloctmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
