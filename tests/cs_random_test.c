/*============================================================================
 * Unit test for random number generator.
 *============================================================================*/

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

#include "cs_defs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bft_error.h>
#include <bft_printf.h>

#include "cs_timer.h"

#include "cs_random.h"

/*---------------------------------------------------------------------------*/

#define NPTS 5000

/*---------------------------------------------------------------------------*/

static void
_uniform_test(cs_lnum_t   n,
              cs_real_t  *a)
{
  double diff;
  double b[607];
  double svblk[1634];
  double t1, t2, t3;
  int i, k, ii;
  int ia[20];

  /* number of iterations of test: nits */

  int nits = 128;

  for (k = 0; k < 20; ++k) {
    ia[k] = 0;
  }

  t1 = 100.;
  for (k = 0; k < nits; ++k) {
    t2 = cs_timer_wtime();
    cs_random_uniform(n, a);
    t3 = cs_timer_wtime();
    t2 = t3 - t2;
    t1 = CS_MIN(t2, t1);
    for (i = 0; i < n; ++i) {
      ii = (int) (a[i] * (double)20.);
      ++ia[ii];
    }

    /*  last time, save klotz0 for save/resore test */

    if (k == nits - 2) {
      cs_random_save(svblk);
    }
  }

  /*  test save/restore sequence */

  cs_random_restore(svblk);
  cs_random_uniform(607, b);
  diff = 0.;
  for (i = 0; i < (CS_MIN(n, 607)); ++i) {
    diff += (b[i] - a[i])*(b[i] - a[i]);
  }

  if (fabs(diff) > 1e-18)
    printf("ERROR in start/restart: diff = %e\n", diff);
  else
    printf("  cs_random save/restore test OK \n");

  t1 = t1 / ((float)n);
  printf("\n    time/uniform = %e \n",t1);
  printf("\n    Histogram of uniform distribution:\n");
  printf("    --------- -- ------- ------------ \n");
  for (k = 0; k < 20; ++k) {
    if(k<9) printf("    bin[%d]  = %d\n",k+1,ia[k]);
    else    printf("    bin[%d] = %d\n",k+1,ia[k]);
  }
}

static void
_normal_test(cs_lnum_t   n,
             cs_real_t  *x)
{
  double diff;
  int i, k;
  double y[128], boxsv[1634];
  cs_real_t t1, t2, t3;
  double  x1, x2, x3, x4, x5, x6;
  int kk;
  double xx2, xx4;
  int bin[21];

  /* number of iterations of test */

  int nits = 128;

  /* initialize moments */

  x1 = 0.;
  x2 = 0.;
  x3 = 0.;
  x4 = 0.;
  x5 = 0.;
  x6 = 0.;

  cs_random_normal(n, x);
  for (i = 0; i < 21; ++i) {
    bin[i] = 0;
  }

  t1 = 100.;
  for (k = 0; k < nits; ++k) {

    /*  save seeds and pointers for save/restore test */

    if (k == nits-1) {
      cs_random_save(boxsv);
    }

    t2 = cs_timer_wtime();
    cs_random_normal(n, x);
    t3 = cs_timer_wtime();
    /* t2 = t3 - t2; */
    t1 = CS_MIN(t1,t3-t2);

    for (i = 0; i < n; ++i) {
      kk = (int) ((x[i] + 5.25) * 2.);
      ++bin[kk];
    }
    for (i = 0; i < n; ++i) {
      x1 += x[i];
      xx2 = x[i] * x[i];
      x2 += xx2;
      x3 += xx2 * x[i];
      xx4 = xx2 * xx2;
      x4 += xx4;
      x5 += xx4 * x[i];
      x6 += xx4 * xx2;
    }

    /*  restore previous seeds and pointers for save/restore test */

    if (k == nits-1) {
      cs_random_restore(boxsv);
      cs_random_normal(128, y);
    }
  }

  /*  save/restore check: */

  diff = 0.;
  for (i = 0; i < (CS_MIN(n,128)); ++i) {
    diff += (y[i] - x[i])*(y[i] - x[i]);
  }
  if (fabs(diff) > 1e-18)
    printf("ERROR in normalsv/normalrs: diff = %e\n",diff);
  else
    printf("    normal law save/restore test OK\n");

  x1 /= (double) (n * nits);
  x2 /= (double) (n * nits);
  x3 /= (double) (n * nits);
  x4 /= (double) (n * nits);
  x5 /= (double) (n * nits);
  x6 /= (double) (n * nits);

  t1 = t1 / ((float) n);
  printf("\n    Time/normal = %e seconds \n",t1);
  printf("    Moments: \n");
  printf("      Compare to (0.0)               (1.0) \n");
  printf("              %e       %e \n",x1,x2);
  printf("      Compare to (0.0)               (3.0) \n");
  printf("              %e       %e \n",x3,x4);
  printf("      Compare to (0.0)              (15.0) \n");
  printf("              %e       %e \n",x5,x6);
  printf("\n    Histogram of gaussian distribution:\n");
  printf("    --------- -- --------- ------------ \n");
  for (k = 0; k < 21; ++k) {
    if (k<9) printf("    bin[%d]  = %d \n",k+1,bin[k]);
    else     printf("    bin[%d] = %d \n",k+1,bin[k]);
  }
}

static void
_poisson_test(cs_lnum_t   n,
              int        *p)
{
  int nits, i, k;
  double p1,p2,p3,p4,x1,x2,x3,x4,fp;
  float t1,t2,t3;
  int kk;
  double mu;
  int bin[20];

  mu = 2.;
  nits = 128;

  for (k = 0; k < 20; ++k) {
    bin[k] = 0;
  }

  /* moment comparison values */

  p1 = mu;
  p2 = mu + mu * mu;
  p3 = mu + mu * (double)3. * mu + mu * mu * mu;
  p4 = mu + 7.*mu*mu + 6.*mu*mu*mu + mu*mu*mu*mu;

  x1 = 0.;
  x2 = 0.;
  x3 = 0.;
  x4 = 0.;

  t1 = 10.;
  for (k = 0; k < nits; ++k) {

    t2 = cs_timer_wtime();
    cs_random_poisson(n, mu, p);
    t3 = cs_timer_wtime();
    t2 = t3 - t2;
    t1 = CS_MIN(t1,t2);

    for (i = 0; i < n; ++i) {
      kk = p[i];
      ++bin[kk];
    }

    for (i = 0; i < n; ++i) {
      fp = (double) p[i];
      x1 += fp;
      x2 += fp * fp;
      x3 += fp * fp * fp;
      x4 += fp * fp * fp * fp;
    }

  }

  x1 /= (double) (n * nits);
  x2 /= (double) (n * nits);
  x3 /= (double) (n * nits);
  x4 /= (double) (n * nits);

  t1 = t1 / (float) n;
  printf("\n    Time/poisson = %e seconds \n",t1);
  printf("    Moments: \n");
  printf("       Compare: (%e)           (%e)\n",p1,p2);
  printf("                 %e             %e \n",x1,x2);
  printf("       Compare: (%e)           (%e)\n",p3,p4);
  printf("                 %e             %e \n",x3,x4);
  printf("\n    Histogram of Poisson distribution: mu = %e\n",mu);
  printf("    --------- -- ------- ------------ \n");
  for (k = 0; k < 20; ++k) {
    if(k<9) printf("    bin[%d]  = %d \n",k+1,bin[k]);
    else     printf("    bin[%d] = %d \n",k+1,bin[k]);
  }
}

/*---------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  int p[NPTS];
  double a[NPTS];
  double wt0, wt1;

  /* Initialize the seeds for the uniform random number generator */

  int seed = 0;

  wt0 = cs_timer_wtime();

  cs_random_seed(seed);

  wt1 = cs_timer_wtime();

  printf("Random number generator initialized in %f seconds\n",
         wt1 - wt0);

  wt0 = cs_timer_wtime();

  _uniform_test(NPTS, a);

  wt1 = cs_timer_wtime();

  printf("Uniform distribution for %d values in %f seconds\n",
         NPTS, wt1 - wt0);

  wt0 = cs_timer_wtime();

  _normal_test(NPTS, a);

  wt1 = cs_timer_wtime();

  printf("Normal (Gaussian) distribution for %d values in %f seconds\n",
         NPTS, wt1 - wt0);

  wt0 = cs_timer_wtime();

  _poisson_test(NPTS, p);

  wt1 = cs_timer_wtime();

  printf("Fischer distribution for %d values in %f seconds\n",
         NPTS, wt1 - wt0);

  exit(EXIT_SUCCESS);
}
