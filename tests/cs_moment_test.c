/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/
/* Return Kahan sum of array values                                           */
/*----------------------------------------------------------------------------*/

static double
_sum_kahan(size_t         n,
           const double  *x)
{
  size_t  i;
  double  s = 0, c = 0;

  if (n < 1)
    return s;

  for (i = 0; i < n; i++) {
    double z = x[i] - c;
    double t = s + z;
    c = (t - s) - z;
    s = t;
  }

  return s;
}

/*----------------------------------------------------------------------------*/

# define NR 50

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  double m_ref, m_y_ref, v_ref, c_ref;
  double *xr = NULL, *yr = NULL, *wr = NULL;

  const size_t nr = 50;

  /* Initialization and environment */

  xr = malloc(nr*sizeof(double));
  yr = malloc(nr*sizeof(double));
  wr = malloc(nr*sizeof(double));

  printf("\n");

  srand(2);
  for (size_t i = 0; i < nr; i++) {
    double a = (double)(rand()) / RAND_MAX;
    double b = (double)(rand()) / RAND_MAX;
    wr[i] = (2. + a) * 0.0001;
    xr[i] = 1. + b;
    yr[i] = 1. + a + b;
  }

  /* Reference moments */

  {
    double *wx = malloc(nr*sizeof(double));
    double *wy = malloc(nr*sizeof(double));
    for (size_t i = 0; i < nr; i++) {
      wx[i] = wr[i]*xr[i];
      wy[i] = wr[i]*yr[i];
    }

    double swx = _sum_kahan(nr, wx);
    double swy = _sum_kahan(nr, wy);
    double sw = _sum_kahan(nr, wr);

    m_ref = swx / sw;
    m_y_ref = swy / sw;

    for (size_t i = 0; i < nr; i++) {
      wx[i] = wr[i]*(xr[i]-m_ref)*(xr[i]-m_ref);
      wy[i] = wr[i]*(yr[i]-m_y_ref)*(xr[i]-m_ref);
    }

    swx = _sum_kahan(nr, wx);
    swy = _sum_kahan(nr, wy);

    v_ref = swx / sw;
    c_ref = swy / sw;

    printf("Reference mean:      %12.5g\n"
           "Reference variance:  %12.5g\n"
           "Reference covariance:  %12.5g\n\n",
           m_ref, v_ref, c_ref);

    free(wy);
    free(wx);
  }

  /* Precision tests */
  /*-----------------*/

  for (int t_id = 0; t_id < 7; t_id++) {

    size_t n_test = nr;

    for (int tm = 0; tm < t_id; tm++)
      n_test *= 10;

    /* Variant 1a */

    {
      double ws = 0;
      double m = 0;
      double m_y = 0;
      double m2 = 0;
      double c2 = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double delta_y = yr[j] - m_y;

        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double r_y = delta_y * (wr[j] / (fmax(ws_n, 1e-100)));
        m2 += ws * delta * r;
        c2 += ws * delta * r_y;
        m += r;
        m_y += r_y;
        ws += wr[j];
      }

      double v = m2 / ws;
      double c = c2 / ws;

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      double ce = fabs(c_ref - c);
      printf("Algorithm 1a, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n"
             "  covariance error:  %12.5g\n\n", (int)n_test, me, ve, ce);
    }

    /* Variant 1b */

    {
      double ws = 0;
      double m = 0;
      double m_y = 0;
      double m2 = 0;
      double c2 = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double delta_y = yr[j] - m_y;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double r_y = delta_y * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        m2 += wr[j] * delta * (xr[j]-m_n);
        c2 += wr[j] * delta_y * (xr[j]-m_n);
        m += r;
        m_y += r_y;
        ws += wr[j];
      }

      double v = m2 / ws;
      double c = c2 / ws;

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      double ce = fabs(c_ref - c);
      printf("Algorithm 1b, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n"
             "  covariance error:  %12.5g\n\n", (int)n_test, me, ve, ce);
    }

    /* Variant 2a */

    {
      double ws = 0;
      double m = 0;
      double m_y = 0;
      double v = 0;
      double c = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double delta_y = yr[j] - m_y;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double r_y = delta_y * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        v = v*(ws/ws_n) + (wr[j] * delta * (xr[j]-m_n)/ws_n);
        c = c*(ws/ws_n) + (wr[j] * delta_y * (xr[j]-m_n)/ws_n);
        m += r;
        m_y += r_y;
        ws += wr[j];
      }

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      double ce = fabs(c_ref - c);
      printf("Algorithm 2a, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n"
             "  covariance error:  %12.5g\n\n", (int)n_test, me, ve, ce);
    }

    /* Variant 2b */

    {
      double ws = 0;
      double m = 0;
      double m_y = 0;
      double v = 0;
      double c = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double delta_y = yr[j] - m_y;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double r_y = delta_y * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        v = (v*ws + (wr[j] * delta * (xr[j]-m_n))) / ws_n;
        c = (c*ws + (wr[j] * delta_y * (xr[j]-m_n))) / ws_n;
        m += r;
        m_y += r_y;
        ws += wr[j];
      }

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      double ce = fabs(c_ref - c);
      printf("Algorithm 2b, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n"
             "  covariance error:  %12.5g\n\n", (int)n_test, me, ve, ce);
    }

    /* Variant 3 */

    {
      double ws = 0;
      double m = 0;
      double m_y = 0;
      double c = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double delta_y = yr[j] - m_y;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double r_y = delta_y * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        double m_y_n = m_y + r_y;
        c = (c*ws + (wr[j] * 0.5* (  delta_y * (xr[j]-m_n)
                                   + delta * (yr[j]-m_y_n)))) / ws_n;
        m += r;
        m_y += r_y;
        ws += wr[j];
      }

      double me = fabs(m_ref - m);
      double mye = fabs(m_y_ref - m_y);
      double ce = fabs(c_ref - c);
      printf("Algorithm 3, n_test = %d\n"
             "  mean error x:      %12.5g\n"
             "  mean error y:      %12.5g\n"
             "  covariance error:  %12.5g\n\n", (int)n_test, me, mye, ce);
    }
  }

  free(xr);
  free(wr);

  exit (EXIT_SUCCESS);
}
