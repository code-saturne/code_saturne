/*============================================================================
 * Low-level operator benchmarking
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(__STDC_VERSION__)      /* size_t */
#if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#include <stdlib.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_ESSL_H)
#include <essl.h>

#elif defined(HAVE_MKL)
#include <mkl_cblas.h>
#include <mkl_spblas.h>

#elif defined(HAVE_ACML)
#include <acml.h>

#elif defined(HAVE_CBLAS)
#include <cblas.h>

#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_matrix.h"
#include "cs_perio.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_benchmark.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Count number of operations.
 *
 * parameters:
 *   n_runs       <-- Local number of runs
 *   n_ops        <-- Local number of operations
 *   n_ops_single <-- Single-processor equivalent number of operations
 *                    (without ghosts); ignored if 0
 *   wt           <-- wall-clock time
 *----------------------------------------------------------------------------*/

static void
_print_stats(long    n_runs,
             long    n_ops,
             long    n_ops_single,
             double  wt)
{
  double fm = 1.0 * n_runs / (1.e9 * (CS_MAX(wt, 1)));

  if (cs_glob_n_ranks == 1)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  N ops:       %12ld\n"
                    "  Wall clock:  %12.5e\n"
                    "  GFLOPS:      %12.5e\n"),
                  n_ops, wt/n_runs, n_ops*fm);

#if defined(HAVE_MPI)

  else {

    long n_ops_min, n_ops_max, n_ops_tot;
    double loc_count[2], glob_sum[2], glob_min[2], glob_max[2], fmg;

    loc_count[0] = wt;
    loc_count[1] = n_ops*fm;

    MPI_Allreduce(&n_ops, &n_ops_min, 1, MPI_LONG, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_max, 1, MPI_LONG, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_tot, 1, MPI_LONG, MPI_SUM,
                  cs_glob_mpi_comm);

    MPI_Allreduce(loc_count, glob_min, 2, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_max, 2, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_sum, 2, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);

    /* global flops multiplier */
    fmg = n_runs / (1.e9 * CS_MAX(glob_max[0], 1));

    glob_sum[0] /= n_runs;
    glob_min[0] /= n_runs;
    glob_max[0] /= n_runs;

    if (n_ops_single == 0)
      cs_log_printf
        (CS_LOG_PERFORMANCE,
         _("               Mean         Min          Max          Total\n"
           "  N ops:       %12ld %12ld %12ld %12ld\n"
           "  Wall clock:  %12.5e %12.5e %12.5e\n"
           "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e\n"),
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

    else
      cs_log_printf
        (CS_LOG_PERFORMANCE,
         _("               Mean         Min          Max          Total"
           "        Single\n"
           "  N ops:       %12ld %12ld %12ld %12ld %12ld\n"
           "  Wall clock:  %12.5e %12.5e %12.5e\n"
           "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e %12.5e\n"),
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         n_ops_single,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1],
         n_ops_tot*fmg, n_ops_single*fmg);
  }

#endif

  cs_log_printf_flush(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------
 * Count number of operations.
 *
 * parameters:
 *   n_runs       <-- Local number of runs
 *   n_ops        <-- Local number of operations
 *   wt           <-- wall-clock time
 *----------------------------------------------------------------------------*/

static void
_print_mem_stats(long    n_runs,
                 long    n_ops,
                 double  wt)
{
  double fm = 1.0 * n_runs / (1.e9 * (CS_MAX(wt, 1)));

  if (cs_glob_n_ranks == 1)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  N ops:       %12ld\n"
                    "  Wall clock:  %12.5e\n"
                    "  GB/s:        %12.5e\n"),
                  n_ops, wt/n_runs, n_ops*fm);

#if defined(HAVE_MPI)

  else {

    long n_ops_min, n_ops_max, n_ops_tot;
    double loc_count[2], glob_sum[2], glob_min[2], glob_max[2], fmg;

    loc_count[0] = wt;
    loc_count[1] = n_ops*fm;

    MPI_Allreduce(&n_ops, &n_ops_min, 1, MPI_LONG, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_max, 1, MPI_LONG, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_tot, 1, MPI_LONG, MPI_SUM,
                  cs_glob_mpi_comm);

    MPI_Allreduce(loc_count, glob_min, 2, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_max, 2, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_sum, 2, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);

    /* global flops multiplier */
    fmg = n_runs / (8.e9 * CS_MAX(glob_max[0], 1));

    glob_sum[0] /= n_runs;
    glob_min[0] /= n_runs;
    glob_max[0] /= n_runs;

    cs_log_printf
      (CS_LOG_PERFORMANCE,
       _("               Mean         Min          Max          Total\n"
         "  N ops:       %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GB/s:        %12.5e %12.5e %12.5e %12.5e\n"),
       n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
       glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
       glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

  }

#endif

  cs_log_printf_flush(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------
 * Simple dot product.
 *
 * parameters:
 *   t_measure       <-- minimum time for each measure (< 0 for single pass)
 *   global          <-- 0 for local use, 1 for MPI sum
 *   n_cells         <-- number of cells (array size)
 *   x               <-- Vector
 *   y               <-- Vector
 *----------------------------------------------------------------------------*/

static void
_dot_product_1(double            t_measure,
               int               global,
               cs_lnum_t         n_cells,
               const cs_real_t  *x,
               const cs_real_t  *y)
{
  double wt0, wt1;
  int run_id, n_runs;
  long n_ops;
  cs_lnum_t ii;

  double test_sum = 0.0;
  int _global = global;

  char type_name[] = "X.Y";

  if (x == y)
    type_name[2] = 'X';

  if (cs_glob_n_ranks == 1)
    _global = 0;

  n_ops = CS_MAX(n_cells*2 - 1, 0);

  /* First simple local x.x version */

#if defined(HAVE_BLAS)

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      double s1 = cs_dot(n_cells, x, y);
#if defined(HAVE_MPI)
      if (_global) {
        double s1_glob = 0.0;
        MPI_Allreduce(&s1, &s1_glob, 1, MPI_DOUBLE, MPI_SUM,
                      cs_glob_mpi_comm);
        s1 = s1_glob;
      }
#endif
      test_sum += test_sum_mult*s1;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  if (_global == 0)
    cs_log_printf(CS_LOG_PERFORMANCE,
      _("\n"
         "Simple local dot product %s with BLAS\n"
         "-------------------------------------\n"),
       type_name);
  else
    cs_log_printf(CS_LOG_PERFORMANCE,
      _("\n"
         "Simple global dot product %s with BLAS\n"
         "---------------------------------------\n"),
       type_name);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      double s1 = 0.0;
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        s1 += x[ii] * y[ii];
#if defined(HAVE_MPI)
      if (_global) {
        double s1_glob = 0.0;
        MPI_Allreduce(&s1, &s1_glob, 1, MPI_DOUBLE, MPI_SUM,
                      cs_glob_mpi_comm);
        s1 = s1_glob;
      }
#endif
      test_sum += test_sum_mult*s1;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  if (_global == 0)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Simple local dot product %s\n"
                    "---------------------------\n"),
                  type_name);
  else
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Simple global dot product %s\n"
                    "----------------------------\n"),
                  type_name);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

 }

/*----------------------------------------------------------------------------
 * Double local dot product.
 *
 * parameters:
 *   t_measure       <-- minimum time for each measure (< 0 for single pass)
 *   n_cells         <-- number of cells (array size)
 *   x               <-- Vector
 *   y               <-- Vector
 *----------------------------------------------------------------------------*/

static void
_dot_product_2(double            t_measure,
               cs_lnum_t         n_cells,
               const cs_real_t  *x,
               const cs_real_t  *y)
{
  double wt0, wt1;
  int run_id, n_runs;
  long n_ops;
  cs_lnum_t ii;

  double test_sum = 0.0;

  n_ops = CS_MAX(n_cells*2 - 1, 0) * 2;

  /* First simple local x.x version */

#if defined(HAVE_BLAS)

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      double s1 = cs_dot(n_cells, x, x);
      double s2 = cs_dot(n_cells, x, y);
      test_sum += test_sum_mult*(s1+s2);
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Double local dot product X.X, Y.Y with BLAS\n"
                  "-------------------------------------------\n"));
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      double s1 = 0.0;
      double s2 = 0.0;
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++) {
        s1 += x[ii] * x[ii];
        s2 += x[ii] * y[ii];
      }
      test_sum += test_sum_mult*(s1+s2);
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Double local dot product X.X, Y.Y\n"
                  "---------------------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

 }

/*----------------------------------------------------------------------------
 * y -> ax + y test
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *   x             <-- Vector
 *   y             <-> Vector
 *----------------------------------------------------------------------------*/

static void
_axpy_(double             t_measure,
       cs_lnum_t          n_cells,
       const cs_real_t   *restrict x,
       cs_real_t         *restrict y)
{
  double wt0, wt1;
  int run_id, n_runs;
  long n_ops;
  cs_lnum_t ii;

  double test_sum = 0.0;

  n_ops = n_cells * 2;

  /* First simple local x.x version */

# pragma omp parallel for
  for (ii = 0; ii < n_cells; ii++)
    y[ii] = 0.0;

#if defined(HAVE_BLAS)

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      cs_axpy(n_cells, test_sum_mult, x, y);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- aX + Y with BLAS\n"
                  "---------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        y[ii] += test_sum_mult * x[ii];
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- aX + Y\n"
                  "-----------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  /* Variant with alpha = -1 */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

# pragma omp parallel for
  for (ii = 0; ii < n_cells; ii++)
    y[ii] = 0.0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        y[ii] -= x[ii];
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    if (run_id % 4096 && run_id > 0) {
      for (ii = 0; ii < n_cells; ii++)
        y[ii] = 0.0;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- Y - X\n"
                  "----------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);
}

/*----------------------------------------------------------------------------
 * Simple divisions on a vector.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *----------------------------------------------------------------------------*/

static void
_division_test(double     t_measure,
               cs_lnum_t  n_cells)
{
  double wt0, wt1;
  int run_id, n_runs;
  long n_ops;
  cs_lnum_t ii;

  cs_real_t  *x = NULL, *y = NULL, *z = NULL;

  BFT_MALLOC(x, n_cells, cs_real_t);
  BFT_MALLOC(y, n_cells, cs_real_t);
  BFT_MALLOC(z, n_cells, cs_real_t);

  n_ops = n_cells;

  /* Division */
  /*----------*/

# pragma omp parallel for
  for (ii = 0; ii < n_cells; ii++) {
    x[ii] = 2.0 + ii%3;
    y[ii] = 2.0 + (ii+1)%3;
  }

  /* Division of two vectors */

  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        z[ii] = x[ii] / y[ii];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Division z = x/y\n"
                  "----------------\n"));

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(z);

  /* Copy inverse of a vector */

  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        y[ii] = 1.0 / x[ii];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Division y = 1/x\n"
                  "----------------\n"));

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(y);

  /* Inverse of a vector */

  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        x[ii] = 1.0 / x[ii];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Division x <- 1/x\n"
                  "-----------------\n"));

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(x);
 }

/*----------------------------------------------------------------------------
 * Simple square root on a vector.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *----------------------------------------------------------------------------*/

static void
_sqrt_test(double     t_measure,
           cs_lnum_t  n_cells)
{
  double wt0, wt1;
  int run_id, n_runs;
  long n_ops;
  cs_lnum_t ii;

  cs_real_t  *x = NULL, *y = NULL;

  BFT_MALLOC(x, n_cells, cs_real_t);
  BFT_MALLOC(y, n_cells, cs_real_t);

  n_ops = n_cells;

# pragma omp parallel for
  for (ii = 0; ii < n_cells; ii++)
    x[ii] = 2.0 + ii%3;

  /* Copy square root of a vector */

  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        y[ii] = sqrt(x[ii]);
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "y = sqrt(x)\n"
                  "-----------\n"));

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(y);

  /* In place square root of a vector */

  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        x[ii] = sqrt(x[ii]);
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "x = sqrt(x)\n"
                  "-----------\n"));

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(x);
}

/*----------------------------------------------------------------------------
 * y -> da.x test
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *   da            <-- Vector
 *   x             <-- Vector
 *   y             <-> Vector
 *----------------------------------------------------------------------------*/

static void
_ad_x_test(double             t_measure,
           cs_lnum_t          n_cells,
           const cs_real_t   *restrict da,
           const cs_real_t   *restrict x,
           cs_real_t         *restrict y)
{
  double wt0, wt1;
  int    run_id, n_runs;
  long   n_ops;
  int    ii;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *da)
# endif

  n_ops = n_cells;

# pragma omp parallel for
  for (ii = 0; ii < n_cells; ii++)
    y[ii] = 0.0;

#if defined(HAVE_CBLAS) || defined(HAVE_MKL)

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      cblas_dgbmv(CblasRowMajor, CblasNoTrans,
                  n_cells,  /* n rows */
                  n_cells,  /* n columns */
                  0,        /* kl */
                  0,        /* ku */
                  1.0,      /* alpha */
                  da,       /* matrix */
                  1,        /* lda */
                  x, 1, 0.0, y, 1);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- DX with BLAS dgbmv\n"
                  "-----------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_CBLAS) || defined(HAVE_MKL) */

#if defined(HAVE_MKL)

  /* Diagonal variant */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    char transa = 'n';
    int n_rows = n_cells;
    int ndiag = 1;
    int idiag[1] = {0};
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      mkl_ddiagemv(&transa,
                   &n_rows,
                   da,
                   &n_rows,
                   idiag,
                   &ndiag,
                   x,
                   y);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- DX with MKL ddiagemv\n"
                  "-------------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_MKL) */

#if defined(HAVE_ESSL)

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      dndot(n_cells, 1, y, 1, 1, da, 1, 1, x, 1, 1);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- DX with ESSL dndot\n"
                  "-----------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

#endif /* defined(HAVE_ESSL) */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
#     pragma omp parallel for
      for (ii = 0; ii < n_cells; ii++)
        y[ii] = da[ii] * x[ii];
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Y <- DX\n"
                  "-------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

}

/*----------------------------------------------------------------------------
 * Compute matrix-vector product for dense blocks: y[i] = a[i].x[i]
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization,
 * and an inline function
 *
 * parameters:
 *   b_id   <-- block id
 *   a      <-- Pointer to block matrixes array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_ax(cs_lnum_t         b_id,
              const cs_real_t  *restrict a,
              const cs_real_t  *restrict x,
              cs_real_t        *restrict y)
{
  cs_lnum_t   jj, kk;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, * a)
# endif

  for (jj = 0; jj < 3; jj++) {
    y[b_id*3 + jj] = 0;
    for (kk = 0; kk < 3; kk++)
      y[b_id*3 + jj] +=   a[b_id*9 + jj*3 + kk] * x[b_id*3 + kk];
  }
}

static inline void
_3_3_diag_vec_p_l_a(const cs_real_t  *restrict da,
                    const cs_real_t  *restrict x,
                    cs_real_t        *restrict y,
                    cs_lnum_t         n_elts)
{
  cs_lnum_t  ii;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, *da)
# endif

# pragma omp parallel for
  for (ii = 0; ii < n_elts; ii++)
    _dense_3_3_ax(ii, da, x, y);
}

/*----------------------------------------------------------------------------
 * Block version of y[i] = da[i].x[i].
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization,
 * and no inline function
 *
 * parameters:
 *   da     <-- Pointer to coefficients array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *   n_elts <-- Array size
 *----------------------------------------------------------------------------*/

static inline void
_3_3_diag_vec_p_l_b(const cs_real_t  *restrict da,
                    const cs_real_t  *restrict x,
                    cs_real_t        *restrict y,
                    cs_lnum_t         n_elts)
{
  cs_lnum_t  ii, jj, kk;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, *da)
# endif

# pragma omp parallel for private(jj, kk)
  for (ii = 0; ii < n_elts; ii++) {
    for (jj = 0; jj < 3; jj++) {
      y[ii*3 + jj] = 0;
      for (kk = 0; kk < 3; kk++)
        y[ii*3 + jj] += da[ii*9 + jj*3 + kk] * x[ii*3 + kk];
    }
  }
}

/*----------------------------------------------------------------------------
 * Block version of y[i] = da[i].x[i], with unrolled loop.
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization,
 * and no inline function
 *
 * parameters:
 *   da     <-- Pointer to coefficients array (usually matrix diagonal)
 *   x      <-- Multipliying vector values
 *   y      --> Resulting vector
 *   n_elts <-- Array size
 *----------------------------------------------------------------------------*/

static inline void
_3_3_diag_vec_p_l_c(const cs_real_t  *restrict da,
                    const cs_real_t  *restrict x,
                    cs_real_t        *restrict y,
                    cs_lnum_t         n_elts)
{
  cs_lnum_t  ii;

# if defined(__xlc__) /* Tell IBM compiler not to alias */
# pragma disjoint(*x, *y, *da)
# endif

# pragma omp parallel for
  for (ii = 0; ii < n_elts; ii++) {
    y[ii*3]     =   da[ii*9]         * x[ii*3]
                  + da[ii*9 + 1]     * x[ii*3 + 1]
                  + da[ii*9 + 2]     * x[ii*3 + 2];
    y[ii*3 + 1] =   da[ii*9 + 3]     * x[ii*3]
                  + da[ii*9 + 3 + 1] * x[ii*3 + 1]
                  + da[ii*9 + 3 + 2] * x[ii*3 + 2];
    y[ii*3 + 2] =   da[ii*9 + 6]     * x[ii*3]
                  + da[ii*9 + 6 + 1] * x[ii*3 + 1]
                  + da[ii*9 + 6 + 2] * x[ii*3 + 2];
  }
}

/*----------------------------------------------------------------------------
 * y -> da.x test
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *----------------------------------------------------------------------------*/

static void
_block_ad_x_test(double             t_measure,
                 cs_lnum_t          n_cells)
{
  double wt0, wt1;
  int    run_id, n_runs;
  long   n_ops;
  int    ii;
  cs_real_t  *da, *x, *y;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *da)
# endif

  BFT_MALLOC(da, n_cells*15, cs_real_t);
  BFT_MALLOC(x, n_cells*3, cs_real_t);
  BFT_MALLOC(y, n_cells*3, cs_real_t);

  n_ops = n_cells * (3+2) * 3;

  /* First simple local x.x version */

  for (ii = 0; ii < n_cells*15; ii++)
    da[ii] = 1.0;
  for (ii = 0; ii < n_cells*3; ii++) {
    x[ii] = 0.0;
    y[ii] = 0.0;
  }

  /* Variant a */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _3_3_diag_vec_p_l_a(da, x, y, n_cells);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Block Y <- DX (compiler inline)\n"
                  "-------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  /* Variant b */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _3_3_diag_vec_p_l_b(da, x, y, n_cells);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Block Y <- DX (manual inline)\n"
                  "-------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  /* Variant c */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;

  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _3_3_diag_vec_p_l_c(da, x, y, n_cells);
      test_sum += test_sum_mult*y[run_id%n_cells];
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Block Y <- DX (3x3 manually unrolled)\n"
                  "-------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, 0, wt1 - wt0);

  BFT_FREE(da);
  BFT_FREE(x);
  BFT_FREE(y);
}

/*----------------------------------------------------------------------------
 * Memory copy test using different variants.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *   n_cells       <-- number of cells (array size)
 *   x             <-- Vector
 *   y             <-> Vector
 *----------------------------------------------------------------------------*/

static void
_copy_test(double             t_measure,
           cs_lnum_t          n_cells,
           const cs_real_t   *restrict x,
           cs_real_t         *restrict y)
{
  double     wt0, wt1;
  int        sub_id, run_id, n_runs;
  cs_lnum_t  _n_cells, n_div, ii;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y)
# endif

  /* Copy with memcpy */

  for (sub_id = 0, n_div = 1;
       sub_id < 4;
       sub_id++, n_div *= 8) {

    test_sum = 0.0;
    wt0 = cs_timer_wtime(), wt1 = wt0;

    _n_cells = n_cells / n_div;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        memcpy(y, x, _n_cells*sizeof(cs_real_t));
        y[0] += 0.1;
        test_sum += test_sum_mult*y[run_id%n_cells];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Copy with memcpy (n_cells/%d)\n"
                    "----------------\n"), n_div);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  (calls: %d;  test sum: %12.5f)\n"),
                  n_runs, test_sum);

    _print_mem_stats(n_runs, _n_cells, wt1 - wt0);

  }

  /* Copy with loop */

  for (sub_id = 0, n_div = 1;
       sub_id < 4;
       sub_id++, n_div *= 8) {

    test_sum = 0.0;
    wt0 = cs_timer_wtime(), wt1 = wt0;

    _n_cells = n_cells / n_div;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < _n_cells; ii++)
          y[ii] = x[ii];
        y[0] += 0.1;
        test_sum += test_sum_mult*y[run_id%n_cells];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Copy with loop (n_cells/%d)\n"
                    "--------------\n"), n_div);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  (calls: %d;  test sum: %12.5f)\n"),
                  n_runs, test_sum);

    _print_mem_stats(n_runs, _n_cells, wt1 - wt0);

  }
}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product related performance.
 *
 * parameters:
 *   t_measure   <-- minimum time for each measure (< 0 for single pass)
 *   m_variant   <-- matrix type
 *   sym_coeffs  <-- symmetric coefficients
 *   n_cells     <-- number of local cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   cell_num    <-- global cell numbers (1 to n)
 *   face_cell   <-- face -> cells connectivity (1 to n)
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   da          <-- diagonal values
 *   xa          <-- extradiagonal values
 *   x           <-> vector
 *   y           --> vector
 *----------------------------------------------------------------------------*/

static void
_matrix_vector_test(double                 t_measure,
                    cs_matrix_variant_t   *m_variant,
                    bool                   sym_coeffs,
                    cs_int_t               n_cells,
                    cs_int_t               n_cells_ext,
                    cs_int_t               n_faces,
                    const cs_gnum_t       *cell_num,
                    const cs_int_t        *face_cell,
                    const cs_halo_t       *halo,
                    const cs_numbering_t  *numbering,
                    const cs_real_t       *restrict da,
                    const cs_real_t       *restrict xa,
                    cs_real_t             *restrict x,
                    cs_real_t             *restrict y)
{
  cs_int_t ii;
  double wt0, wt1;
  int    run_id, n_runs;
  long   n_ops, n_ops_glob;

  double test_sum = 0.0;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;
  cs_matrix_type_t m_type = cs_matrix_variant_type(m_variant);

  /* n_cells + n_faces*2 nonzeroes,
     n_row_elts multiplications + n_row_elts-1 additions per row */

  n_ops = n_cells + n_faces*4;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (cs_glob_mesh->n_g_cells + cs_glob_mesh->n_g_i_faces*4);

  ms = cs_matrix_structure_create(m_type,
                                  true,
                                  n_cells,
                                  n_cells_ext,
                                  n_faces,
                                  cell_num,
                                  face_cell,
                                  halo,
                                  numbering);

  m = cs_matrix_create_tuned(ms, m_variant);

  cs_matrix_set_coefficients(m,
                             sym_coeffs,
                             NULL,
                             da,
                             xa);

  /* Matrix.vector product */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      cs_matrix_vector_multiply(CS_PERIO_ROTA_COPY, m, x, y);
      test_sum += y[n_cells-1]*test_sum_mult;
      run_id++;
#if 0
      for (ii = 0; ii < n_cells; ii++)
        cs_log_printf(CS_LOG_PERFORMANCE,
                      "y[%d] = %12.4f\n", ii, y[ii]);
#endif
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  if (sym_coeffs == true)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Matrix.vector product (symm coeffs)\n"
                    "---------------------\n"));
  else
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Matrix.vector product\n"
                    "---------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  /* Local timing in parallel mode */

  if (cs_glob_n_ranks > 1) {

    test_sum = 0.0;
    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        cs_matrix_vector_multiply_nosync(m,
                                         x,
                                         y);
        test_sum += y[n_cells-1]*test_sum_mult;
        if (run_id > 0 && run_id % 64) {
          for (ii = n_cells; ii < n_cells_ext; ii++)
            y[ii] = 0;
        }
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Local matrix.vector product\n"
                    "---------------------------\n"));

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  (calls: %d;  test sum: %12.5f)\n"),
                  n_runs, test_sum);

    _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  }

  /* (Matrix - diagonal).vector product */

  /* n_faces*2 nonzeroes,
     n_row_elts multiplications + n_row_elts-1 additions per row */

  n_ops = n_faces*4 - n_cells;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (cs_glob_mesh->n_g_i_faces*4 - cs_glob_mesh->n_g_cells);

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      cs_matrix_exdiag_vector_multiply(CS_PERIO_ROTA_COPY,
                                       m,
                                       x,
                                       y);
      test_sum += y[n_cells-1]*test_sum_mult;
      if (run_id > 0 && run_id % 64) {
        for (ii = n_cells; ii < n_cells_ext; ii++)
          y[ii] = 0;
      }
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "(Matrix-diagonal).vector product (%s)\n"
                  "--------------------------------\n"),
                _(cs_matrix_type_name[m_type]));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  cs_matrix_destroy(&m);
  cs_matrix_structure_destroy(&ms);

}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product extradiagonal terms related performance
 * (symmetric matrix case).
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity (1 to n)
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_exdiag_native(cs_int_t             n_faces,
                       const cs_int_t      *face_cell,
                       const cs_real_t     *restrict xa,
                       cs_real_t           *restrict x,
                       cs_real_t           *restrict y)
{
  cs_int_t  ii, jj, face_id;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *y, *xa)
#endif

  const cs_int_t *restrict face_cel_p = face_cell;

  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = *face_cel_p++ - 1;
    jj = *face_cel_p++ - 1;
    y[ii] += xa[face_id] * x[jj];
    y[jj] += xa[face_id] * x[ii];
  }
}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product extradiagonal terms related performance
 * (symmetric matrix case, variant 1).
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity (1 to n)
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_exdiag_native_v1(cs_int_t             n_faces,
                          const cs_int_t      *face_cell,
                          const cs_real_t     *restrict xa,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  cs_int_t  ii, ii_prev, kk, face_id, kk_max;
  cs_real_t y_it, y_it_prev;

  const int l1_cache_size = 508;

  /*
   * 1/ Split y[ii] and y[jj] computation into 2 loops to remove compiler
   *    data dependency assertion between y[ii] and y[jj].
   * 2/ keep index (*face_cel_p) in L1 cache from y[ii] loop to y[jj] loop
   *    and xa in L2 cache.
   * 3/ break high frequency occurence of data dependency from one iteration
   *    to another in y[ii] loop (nonzero matrix value on the same line ii).
   */

  const cs_int_t *restrict face_cel_p = face_cell;

  for (face_id = 0;
       face_id < n_faces;
       face_id += l1_cache_size) {

    kk_max = CS_MIN((n_faces - face_id), l1_cache_size);

    /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

    ii = face_cel_p[0] - 1;
    ii_prev = ii;
    y_it_prev = y[ii_prev] + xa[face_id] * x[face_cel_p[1] - 1];

    for (kk = 1; kk < kk_max; ++kk) {
      ii = face_cel_p[2*kk] - 1;
      /* y[ii] += xa[face_id+kk] * x[jj]; */
      if(ii == ii_prev) {
        y_it = y_it_prev;
      }
      else {
        y_it = y[ii];
        y[ii_prev] = y_it_prev;
      }
      ii_prev = ii;
      y_it_prev = y_it + xa[face_id+kk] * x[face_cel_p[2*kk+1] - 1];
    }
    y[ii] = y_it_prev;

    /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

    for (kk = 0; kk < kk_max; ++kk) {
      y[face_cel_p[2*kk+1] - 1]
        += xa[face_id+kk] * x[face_cel_p[2*kk] - 1];
    }
    face_cel_p += 2 * l1_cache_size;
  }
}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product extradiagonal terms related performance
 * with contribution to face-based array instead of cell-based array
 * (symmetric matrix case).
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity (1 to n)
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   ya              <-> vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_exdiag_part_p1(cs_int_t             n_faces,
                        const cs_int_t      *face_cell,
                        const cs_real_t     *restrict xa,
                        cs_real_t           *restrict x,
                        cs_real_t           *restrict ya)
{
  cs_int_t  ii, jj, face_id;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*x, *xa, *ya)
#endif

  const cs_int_t *restrict face_cel_p = face_cell;

  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = *face_cel_p++ - 1;
    jj = *face_cel_p++ - 1;
    ya[face_id] += xa[face_id] * x[ii];
    ya[face_id] += xa[face_id] * x[jj];
  }
}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product local extradiagonal part related performance.
 *
 * parameters:
 *   t_measure   <-- minimum time for each measure (< 0 for single pass)
 *   n_cells     <-- number of cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   face_cell   <-- face -> cells connectivity (1 to n)
 *   xa          <-- extradiagonal values
 *   x           <-> vector
 *   y           --> vector
 *----------------------------------------------------------------------------*/

static void
_sub_matrix_vector_test(double               t_measure,
                        cs_int_t             n_cells,
                        cs_int_t             n_cells_ext,
                        cs_int_t             n_faces,
                        const cs_int_t      *face_cell,
                        const cs_real_t     *restrict xa,
                        cs_real_t           *restrict x,
                        cs_real_t           *restrict y)
{
  cs_int_t  jj;
  double wt0, wt1;
  int    run_id, n_runs;
  long   n_ops, n_ops_glob;
  double *ya = NULL;

  double test_sum = 0.0;

  /* n_faces*2 nonzeroes,
     n_row_elts multiplications + n_row_elts-1 additions per row */

  n_ops = n_faces*4 - n_cells;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (cs_glob_mesh->n_g_i_faces*4 - cs_glob_mesh->n_g_cells);

  for (jj = 0; jj < n_cells_ext; jj++)
    y[jj] = 0.0;

  /* Matrix.vector product, variant 0 */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _mat_vec_exdiag_native(n_faces, face_cell, xa, x, y);
      test_sum += y[n_cells-1]*test_sum_mult;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Matrix.vector product, extradiagonal part, variant 0\n"
                  "---------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  for (jj = 0; jj < n_cells_ext; jj++)
    y[jj] = 0.0;

  test_sum = 0.0;

  /* Matrix.vector product, variant 1 */

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _mat_vec_exdiag_native_v1(n_faces, face_cell, xa, x, y);
      test_sum += y[n_cells-1]*test_sum_mult;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Matrix.vector product, extradiagonal part, variant 1\n"
                  "---------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  /* Matrix.vector product, contribute to faces only */

  /* n_faces*2 nonzeroes, n_row_elts multiplications */

  n_ops = n_faces*2;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (cs_glob_mesh->n_g_i_faces*2);

  BFT_MALLOC(ya, n_faces, cs_real_t);
  for (jj = 0; jj < n_faces; jj++)
    ya[jj] = 0.0;

  test_sum = 0.0;
  wt0 = cs_timer_wtime(), wt1 = wt0;
  if (t_measure > 0)
    n_runs = 8;
  else
    n_runs = 1;
  run_id = 0;
  while (run_id < n_runs) {
    double test_sum_mult = 1.0/n_runs;
    while (run_id < n_runs) {
      _mat_vec_exdiag_part_p1(n_faces, face_cell, xa, x, ya);
      test_sum += y[n_cells-1]*test_sum_mult;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  BFT_FREE(ya);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Matrix.vector product, face values only\n"
                  "---------------------\n"));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  (calls: %d;  test sum: %12.5f)\n"),
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Run simple benchmarks.
 *
 * parameters:
 *   mpi_trace_mode <-- indicates if timing mode (0) or MPI trace-friendly
 *                      mode (1) is to be used
 *----------------------------------------------------------------------------*/

void
cs_benchmark(int  mpi_trace_mode)
{
  /* Local variable definitions */
  /*----------------------------*/

  size_t ii;

  double t_measure = (mpi_trace_mode) ? -1.0 : 3.0;

  cs_real_t *x1 = NULL, *x2 = NULL;
  cs_real_t *y = NULL;
  cs_real_t *da = NULL, *xa = NULL;

  cs_matrix_variant_t *mv = NULL;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_v = cs_glob_mesh_quantities;

  size_t n_cells = mesh->n_cells;
  size_t n_cells_ext = mesh->n_cells_with_ghosts;
  size_t n_faces = mesh->n_i_faces;

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  BFT_MALLOC(x1, n_cells_ext, cs_real_t);
  BFT_MALLOC(x2, n_cells_ext, cs_real_t);

  for (ii = 0; ii < n_cells_ext; ii++) {
    x1[ii] = mesh_v->cell_cen[ii*3];
    x2[ii] = mesh_v->cell_cen[ii*3 + 1];
  }

  if (CS_MEM_ALIGN > 0)
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n_cells_ext, cs_real_t);
  else
    BFT_MALLOC(y, n_cells_ext, cs_real_t);

  BFT_MALLOC(da, n_cells_ext, cs_real_t);
  BFT_MALLOC(xa, n_faces*2, cs_real_t);

  for (ii = 0; ii < n_cells_ext; ii++)
    da[ii] = 1.0;

  for (ii = 0; ii < n_faces; ii++) {
    xa[ii*2] = 0.5;
    xa[ii*2 + 1] = -0.5;
  }

  /* Run tests */
  /*-----------*/

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Benchmark mode activated\n"
                  "========================\n"));

  /* Dot product test */
  /*------------------*/

  _dot_product_1(t_measure, 0, n_cells, x1, x2);
  _dot_product_2(t_measure, n_cells, x1, x2);
  _axpy_(t_measure, n_cells, x1, y);

  _division_test(t_measure, n_cells);
  _sqrt_test(t_measure, n_cells);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    _dot_product_1(t_measure, 1, n_cells, x1, x1);

#endif /* HAVE_MPI */

  _ad_x_test(t_measure, n_cells, da, x1, y);

  _block_ad_x_test(t_measure, n_cells);

  _copy_test(t_measure, n_cells, x1, y);

  /* Call matrix tuning */
  /*--------------------*/

  /* Test local matrix.vector product operations. */

  cs_matrix_variant_test(n_cells,
                         n_cells_ext,
                         n_faces,
                         mesh->global_cell_num,
                         mesh->i_face_cells,
                         mesh->halo,
                         mesh->i_face_numbering);

  /* Enter tuning phase */

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "General tuning for matrices\n"
                  "=====================================\n"));

  mv = cs_matrix_variant_tuned(t_measure,
                               0.000001, /* sym_weight */
                               0.000001, /* block_weight */
                               10,       /* min expected SpMV products */
                               n_cells,
                               n_cells_ext,
                               n_faces,
                               mesh->global_cell_num,
                               mesh->i_face_cells,
                               mesh->halo,
                               mesh->i_face_numbering);

  _matrix_vector_test(t_measure,
                      mv, false,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y);

  cs_matrix_variant_destroy(&mv);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Tuning for symmetric matrices\n"
                  "=============================\n"));

  mv = cs_matrix_variant_tuned(t_measure,
                               1.0, /* sym_weight */
                               0.0, /* block_weight */
                               10,  /* min expected SpMV products */
                               n_cells,
                               n_cells_ext,
                               n_faces,
                               mesh->global_cell_num,
                               mesh->i_face_cells,
                               mesh->halo,
                               mesh->i_face_numbering);

  _matrix_vector_test(t_measure,
                      mv, true,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y);

  cs_matrix_variant_destroy(&mv);

  _sub_matrix_vector_test(t_measure,
                          n_cells,
                          n_cells_ext,
                          n_faces,
                          mesh->i_face_cells,
                          xa,
                          x1,
                          y);

  cs_log_separator(CS_LOG_PERFORMANCE);

  /* Free working arrays */
  /*---------------------*/

  BFT_FREE(x1);
  BFT_FREE(x2);

  BFT_FREE(y);

  BFT_FREE(da);
  BFT_FREE(xa);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
