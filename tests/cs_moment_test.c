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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

/*----------------------------------------------------------------------------*/

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
  int sub_id;
  double m_ref, v_ref;
  double *xr = NULL, *wr = NULL;

  const size_t nr = 50;

  /* Initialization and environment */

  xr = malloc(nr*sizeof(double));
  wr = malloc(nr*sizeof(double));

  printf("\n");

  srand(2);
  for (size_t i = 0; i < nr; i++) {
    double a = (double)(rand()) / RAND_MAX;
    double b = (double)(rand()) / RAND_MAX;
    wr[i] = (2. + a) * 0.0001;
    xr[i] = 1. + b;
  }

  /* Reference moments */

  {
    double *wx = malloc(nr*sizeof(double));
    for (size_t i = 0; i < nr; i++)
      wx[i] = wr[i]*xr[i];

    double swx = _sum_kahan(nr, wx);
    double sw = _sum_kahan(nr, wr);

    m_ref = swx / sw;

    for (size_t i = 0; i < nr; i++)
      wx[i] = wr[i]*(xr[i]-m_ref)*(xr[i]-m_ref);

    swx = _sum_kahan(nr, wx);

    v_ref = swx / sw;

    printf("Reference mean:      %12.5g\n"
           "Reference variance:  %12.5g\n\n", m_ref, v_ref);

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
      double m2 = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        m2 += ws * delta * r;
        m += r;
        ws += wr[j];
      }

      double v = m2 / ws;

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      printf("Algorithm 1a, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n\n", n_test, me, ve);
    }

    /* Variant 1b */

    {
      double ws = 0;
      double m = 0;
      double m2 = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        m2 += wr[j] * delta * (xr[j]-m_n);
        m += r;
        ws += wr[j];
      }

      double v = m2 / ws;

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      printf("Algorithm 1b, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n\n", n_test, me, ve);
    }

    /* Variant 2a */

    {
      double ws = 0;
      double m = 0;
      double v = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        v = v*(ws/ws_n) + (wr[j] * delta * (xr[j]-m_n)/ws_n);
        m += r;
        ws += wr[j];
      }

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      printf("Algorithm 2a, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n\n", n_test, me, ve);
    }

    /* Variant 2b */

    {
      double ws = 0;
      double m = 0;
      double v = 0;

      for (size_t i = 0; i < n_test; i++) {
        size_t j = i % nr;
        double ws_n = wr[j] + ws;
        double delta = xr[j] - m;
        double r = delta * (wr[j] / (fmax(ws_n, 1e-100)));
        double m_n = m + r;
        v = (v*ws + (wr[j] * delta * (xr[j]-m_n))) / ws_n;
        m += r;
        ws += wr[j];
      }

      double me = fabs(m_ref - m);
      double ve = fabs(v_ref - v);
      printf("Algorithm 2b, n_test = %d\n"
             "  mean error:        %12.5g\n"
             "  variance error:    %12.5g\n\n", n_test, me, ve);
    }
  }

  free(xr);
  free(wr);

  exit (EXIT_SUCCESS);
}
