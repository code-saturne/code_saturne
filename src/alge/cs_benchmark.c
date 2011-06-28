/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Low-level operator benchmarking
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_matrix.h"
#include "cs_perio.h"

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
 * Start timer.
 *
 * parameters:
 *   wt       --> wall-clock time (start in, stop - start out)
 *   cpu      --> CPU time (start in, stop - start out)
 *----------------------------------------------------------------------------*/

static void
_timer_start(double  *wt,
             double  *cpu)
{
  *wt = bft_timer_wtime();
  *cpu = bft_timer_cpu_time();
}

/*----------------------------------------------------------------------------
 * Stop timer.
 *
 * parameters:
 *   n_runs     <-- Number of timing runs
 *   wt         <-> wall-clock time (start in, stop - start out)
 *   cpu        <-> CPU time (start in, stop - start out)
 *----------------------------------------------------------------------------*/

static void
_timer_stop(int      n_runs,
            double  *wt,
            double  *cpu)
{
  double wt_s, cpu_s;

  wt_s = bft_timer_wtime();
  cpu_s = bft_timer_cpu_time();

  *wt = (wt_s - *wt) / (double)n_runs;
  *cpu = (cpu_s - *cpu) / (double)n_runs;
 }

/*----------------------------------------------------------------------------
 * Print overhead.
 *
 * parameters:
 *   wt       <-- wall-clock time
 *   cpu      <-- CPU time
 *----------------------------------------------------------------------------*/

static void
_print_overhead(double  wt,
                double  cpu)
{
  if (cs_glob_n_ranks == 1)
    bft_printf(_("  Wall clock : %12.5e\n"
                 "  CPU :        %12.5e\n"),
               wt, cpu);

  else {

    double loc_count[2], glob_min[2], glob_max[2];

    loc_count[0] = wt;
    loc_count[1] = cpu;

#if defined(HAVE_MPI)
    MPI_Allreduce(loc_count, glob_min, 2, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_max, 2, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);
#else
    { /* We should never enter here unless we have an alternative to MPI */
      int i;
      for (i = 0; i < 3; i++) {
        glob_min[i] = loc_count[i];
        glob_max[i] = loc_count[i];
      }
    }
#endif

    bft_printf(_("               Min          Max\n"
                 "  Wall clock : %12.5e %12.5e\n"
                 "  CPU :        %12.5e %12.5e\n"),
               glob_min[0], glob_max[0], glob_min[1], glob_max[1]);
  }
}

/*----------------------------------------------------------------------------
 * Count number of operations.
 *
 * parameters:
 *   n_ops        <-- Local number of operations
 *   n_ops_single <-- Single-processor equivalent number of operations
 *                    (without ghosts); ignored if 0
 *   wt           <-- wall-clock time
 *   cpu          <-- CPU time
 *----------------------------------------------------------------------------*/

static void
_print_stats(long    n_ops,
             long    n_ops_single,
             double  wt,
             double  cpu)
{
  double fm = 1.0 / (1.e9 * wt);

  if (cs_glob_n_ranks == 1)
    bft_printf(_("  N ops :      %12ld\n"
                 "  Wall clock : %12.5e\n"
                 "  CPU :        %12.5e\n"
                 "  GFLOPS :     %12.5e\n"),
               n_ops, wt, cpu, n_ops*fm);

  else {

    long n_ops_min, n_ops_max, n_ops_tot;
    double loc_count[3], glob_min[3], glob_max[3], fmg;

    loc_count[0] = wt;
    loc_count[1] = cpu;
    loc_count[2] = n_ops*fm;

#if defined(HAVE_MPI)

    MPI_Allreduce(&n_ops, &n_ops_min, 1, MPI_LONG, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_max, 1, MPI_LONG, MPI_MAX,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&n_ops, &n_ops_tot, 1, MPI_LONG, MPI_SUM,
                  cs_glob_mpi_comm);

    MPI_Allreduce(loc_count, glob_min, 3, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);
    MPI_Allreduce(loc_count, glob_max, 3, MPI_DOUBLE, MPI_MAX,
                  cs_glob_mpi_comm);

#else
    { /* We should never enter here unless we have an alternative to MPI */
      int i;
      n_ops_min = n_ops; n_ops_max = n_ops; n_ops_tot = n_ops;
      for (i = 0; i < 3; i++) {
        glob_min[i] = loc_count[i];
        glob_max[i] = loc_count[i];
      }
    }
#endif

    fmg = 1.0 / (1.e9 * glob_max[0]); /* global flops multiplier */

    if (n_ops_single == 0)
      bft_printf
        (_("               Min          Max          Total\n"
           "  N ops :      %12ld %12ld %12ld\n"
           "  Wall clock : %12.5e %12.5e\n"
           "  CPU :        %12.5e %12.5e\n"
           "  GFLOPS :     %12.5e %12.5e %12.5e\n"),
         n_ops_min, n_ops_max, n_ops_tot,
         glob_min[0], glob_max[0],
         glob_min[1], glob_max[1],
         glob_min[2], glob_max[2], n_ops_tot*fmg);

    else
      bft_printf
        (_("               Min          Max          Total        Single\n"
           "  N ops :      %12ld %12ld %12ld %12ld\n"
           "  Wall clock : %12.5e %12.5e\n"
           "  CPU :        %12.5e %12.5e\n"
           "  GFLOPS :     %12.5e %12.5e %12.5e %12.5e\n"),
         n_ops_min, n_ops_max, n_ops_tot, n_ops_single,
         glob_min[0], glob_max[0],
         glob_min[1], glob_max[1],
         glob_min[2], glob_max[2], n_ops_tot*fmg, n_ops_single*fmg);
  }
}

/*----------------------------------------------------------------------------
 * Simple dot product.
 *
 * parameters:
 *   global          <-- 0 for local use, 1 for MPI sum
 *   n_runs          <-- number of operation runs
 *   n_cells         <-- number of cells (array size)
 *   x               <-- Vector
 *   y               <-- Vector
 *----------------------------------------------------------------------------*/

static void
_dot_product_1(int                  global,
               int                  n_runs,
               size_t               n_cells,
               const cs_real_t     *x,
               const cs_real_t     *y)
{
  double wt, cpu;
  int    run_id;
  long   n_ops;
  size_t ii;

  double test_sum_mult = 1.0/n_runs;
  double test_sum = 0.0;
  int _global = global;

  char type_name[] = "X.Y";

  if (n_runs < 1)
    return;

  if (x == y)
    type_name[2] = 'X';

  if (cs_glob_n_ranks == 1)
    _global = 0;

  n_ops = n_cells;

  /* First simple local x.x version */

  _timer_start(&wt, &cpu);

#if defined(HAVE_BLAS)

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {
    double s1 = cblas_ddot(n_cells, x, 1, y, 1);
#if defined(HAVE_MPI)
    if (_global) {
      double s1_glob = 0.0;
      MPI_Allreduce(&s1, &s1_glob, 1, MPI_DOUBLE, MPI_SUM,
                    cs_glob_mpi_comm);
      s1 = s1_glob;
    }
#endif
    test_sum += test_sum_mult*s1;
  }

  _timer_stop(n_runs, &wt, &cpu);

  if (_global == 0)
    bft_printf(_("\n"
                 "Simple local dot product %s with BLAS\n"
                 "-------------------------------------\n"),
               type_name);
  else
    bft_printf(_("\n"
                 "Simple global dot product %s with BLAS\n"
                 "---------------------------------------\n"),
               type_name);

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {
    double s1 = 0.0;
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
  }

  _timer_stop(n_runs, &wt, &cpu);

  if (_global == 0)
    bft_printf(_("\n"
                 "Simple local dot product %s\n"
                 "---------------------------\n"),
               type_name);
  else
    bft_printf(_("\n"
                 "Simple global dot product %s\n"
                 "----------------------------\n"),
               type_name);

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

 }

/*----------------------------------------------------------------------------
 * Double local dot product.
 *
 * parameters:
 *   n_runs          <-- number of operation runs
 *   n_cells         <-- number of cells (array size)
 *   x               <-- Vector
 *   y               <-- Vector
 *----------------------------------------------------------------------------*/

static void
_dot_product_2(int                  n_runs,
               size_t               n_cells,
               const cs_real_t     *x,
               const cs_real_t     *y)
{
  double wt, cpu;
  int    run_id;
  long   n_ops;
  size_t ii;

  double test_sum_mult = 1.0/n_runs;
  double test_sum = 0.0;

  if (n_runs < 1)
    return;

  n_ops = n_cells * 2;

  /* First simple local x.x version */

  _timer_start(&wt, &cpu);

#if defined(HAVE_BLAS)

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {
    double s1 = cblas_ddot(n_cells, x, 1, x, 1);
    double s2 = cblas_ddot(n_cells, x, 1, y, 1);
    test_sum += test_sum_mult*(s1+s2);
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Double local dot product X.X, Y.Y with BLAS\n"
               "-------------------------------------------\n"));
  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {
    double s1 = 0.0;
    double s2 = 0.0;
    for (ii = 0; ii < n_cells; ii++) {
      s1 += x[ii] * x[ii];
      s2 += x[ii] * y[ii];
    }
    test_sum += test_sum_mult*(s1+s2);
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Double local dot product X.X, Y.Y\n"
               "---------------------------------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

 }

/*----------------------------------------------------------------------------
 * y -> ax + y test
 *
 * parameters:
 *   n_runs        <-- number of operation runs
 *   n_cells       <-- number of cells (array size)
 *   x             <-- Vector
 *   y             <-> Vector
 *----------------------------------------------------------------------------*/

static void
_axpy_(int                n_runs,
       size_t             n_cells,
       const cs_real_t   *restrict x,
       cs_real_t         *restrict y)
{
  double wt, cpu;
  int    run_id;
  long   n_ops;
  size_t ii;

  double test_sum_mult = 1.0/n_runs;
  double test_sum = 0.0;

  if (n_runs < 1)
    return;

  n_ops = n_cells * 2;

  /* First simple local x.x version */

  for (ii = 0; ii < n_cells; ii++)
    y[ii] = 0.0;

  _timer_start(&wt, &cpu);

#if defined(HAVE_BLAS)

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {

    cblas_daxpy(n_cells, test_sum_mult, x, 1, y, 1);

    test_sum += test_sum_mult*y[run_id%n_cells];

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Y <- aX + Y with BLAS\n"
               "---------------------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

#endif /* defined(HAVE_BLAS) */

  test_sum = 0.0;

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      y[ii] += test_sum_mult * x[ii];

    test_sum += test_sum_mult*y[run_id%n_cells];

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Y <- aX + Y\n"
               "-----------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, 0, wt, cpu);

}

/*----------------------------------------------------------------------------
 * Simple divisions on a vector.
 *
 * parameters:
 *   n_runs          <-- number of operation runs
 *   n_cells         <-- number of cells (array size)
 *----------------------------------------------------------------------------*/

static void
_division_test(int     n_runs,
               size_t  n_cells)
{
  double wt, cpu;
  int    run_id;
  long   n_ops;
  size_t ii;

  cs_real_t  *x = NULL, *y = NULL, *z = NULL;

  if (n_runs < 1)
    return;

  BFT_MALLOC(x, n_cells, cs_real_t);
  BFT_MALLOC(y, n_cells, cs_real_t);
  BFT_MALLOC(z, n_cells, cs_real_t);

  n_ops = n_cells;

  /* Division */
  /*----------*/

  for (ii = 0; ii < n_cells; ii++) {
    x[ii] = 2.0 + ii%3;
    y[ii] = 2.0 + (ii+1)%3;
  }

  /* Division of two vectors */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      z[ii] = x[ii] / y[ii];

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Division z = x/y\n"
               "----------------\n"));

  _print_stats(n_ops, 0, wt, cpu);

  BFT_FREE(z);

  /* Copy inverse of a vector */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      y[ii] = 1.0 / x[ii];

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Division y = 1/x\n"
               "----------------\n"));

  _print_stats(n_ops, 0, wt, cpu);

  BFT_FREE(y);

  /* Inverse of a vector */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      x[ii] = 1.0 / x[ii];

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Division x <- 1/x\n"
               "-----------------\n"));

  _print_stats(n_ops, 0, wt, cpu);

  BFT_FREE(x);

 }

/*----------------------------------------------------------------------------
 * Simple square root on a vector.
 *
 * parameters:
 *   n_runs          <-- number of operation runs
 *   n_cells         <-- number of cells (array size)
 *----------------------------------------------------------------------------*/

static void
_sqrt_test(int     n_runs,
           size_t  n_cells)
{
  double wt, cpu;
  int    run_id;
  long   n_ops;
  size_t ii;

  cs_real_t  *x = NULL, *y = NULL;

  if (n_runs < 1)
    return;

  BFT_MALLOC(x, n_cells, cs_real_t);
  BFT_MALLOC(y, n_cells, cs_real_t);

  n_ops = n_cells;

  for (ii = 0; ii < n_cells; ii++)
    x[ii] = 2.0 + ii%3;

  /* Copy square root of a vector */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      y[ii] = sqrt(x[ii]);

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "y = sqrt(x)\n"
               "-----------\n"));

  _print_stats(n_ops, 0, wt, cpu);

  BFT_FREE(y);

  /* In place square root of a vector */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {

    for (ii = 0; ii < n_cells; ii++)
      x[ii] = sqrt(x[ii]);

  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "x = sqrt(x)\n"
               "-----------\n"));

  _print_stats(n_ops, 0, wt, cpu);

  BFT_FREE(x);

}

/*----------------------------------------------------------------------------
 * Measure matrix creation + destruction related performance.
 *
 * parameters:
 *   n_runs      <-- number of operation runs
 *   type_name   <-- type name
 *   type        <-- matrix type
 *   n_cells     <-- number of local cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   cell_num    <-- global cell numbers (1 to n)
 *   face_cell   <-- face -> cells connectivity (1 to n)
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *----------------------------------------------------------------------------*/

static void
_matrix_creation_test(int                    n_runs,
                      const char            *type_name,
                      cs_matrix_type_t       type,
                      cs_int_t               n_cells,
                      cs_int_t               n_cells_ext,
                      cs_int_t               n_faces,
                      const fvm_gnum_t      *cell_num,
                      const cs_int_t        *face_cell,
                      const cs_halo_t       *halo,
                      const cs_numbering_t  *numbering)
{
  double wt, cpu;
  int    run_id;

  cs_matrix_structure_t  *ms = NULL;

  if (n_runs < 1)
    return;

  /* First count creation/destruction overhead */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    ms = cs_matrix_structure_create(type,
                                    true,
                                    n_cells,
                                    n_cells_ext,
                                    n_faces,
                                    cell_num,
                                    face_cell,
                                    halo,
                                    numbering);
    cs_matrix_structure_destroy(&ms);
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Matrix structure construction / destruction (%s)\n"
               "------------------------------------------------\n"),
             _(type_name));

  bft_printf(_("  (calls: %d)\n"), n_runs);

  _print_overhead(wt, cpu);

}

/*----------------------------------------------------------------------------
 * Measure matrix assignment performance.
 *
 * parameters:
 *   n_runs      <-- number of operation runs
 *   type_name   <-- type name
 *   type        <-- matrix type
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
 *----------------------------------------------------------------------------*/

static void
_matrix_assignment_test(int                    n_runs,
                        const char            *type_name,
                        cs_matrix_type_t       type,
                        cs_bool_t              sym_coeffs,
                        cs_int_t               n_cells,
                        cs_int_t               n_cells_ext,
                        cs_int_t               n_faces,
                        const fvm_gnum_t      *cell_num,
                        const cs_int_t        *face_cell,
                        const cs_halo_t       *halo,
                        const cs_numbering_t  *numbering,
                        const cs_real_t       *restrict da,
                        const cs_real_t       *restrict xa)
{
  double wt, cpu;
  int    run_id;

  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;

  if (n_runs < 1)
    return;

  /* Count assignment overhead */

  ms = cs_matrix_structure_create(type,
                                  true,
                                  n_cells,
                                  n_cells_ext,
                                  n_faces,
                                  cell_num,
                                  face_cell,
                                  halo,
                                  numbering);

  m = cs_matrix_create(ms);

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    cs_matrix_set_coefficients(m,
                               sym_coeffs,
                               NULL,
                               da,
                               xa);
  }

  _timer_stop(n_runs, &wt, &cpu);

  cs_matrix_destroy(&m);
  cs_matrix_structure_destroy(&ms);

  bft_printf(_("\n"
               "Matrix value assignment (%s)\n"
               "-----------------------\n"), _(type_name));

  bft_printf(_("  (calls: %d)\n"), n_runs);

  _print_overhead(wt, cpu);

}

/*----------------------------------------------------------------------------
 * Measure matrix.vector product related performance.
 *
 * parameters:
 *   n_runs      <-- number of operation runs
 *   type_name   <-- type name
 *   type        <-- matrix type
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
_matrix_vector_test(int                    n_runs,
                    const char            *type_name,
                    cs_matrix_type_t       type,
                    cs_bool_t              sym_coeffs,
                    cs_int_t               n_cells,
                    cs_int_t               n_cells_ext,
                    cs_int_t               n_faces,
                    const fvm_gnum_t      *cell_num,
                    const cs_int_t        *face_cell,
                    const cs_halo_t       *halo,
                    const cs_numbering_t  *numbering,
                    const cs_real_t       *restrict da,
                    const cs_real_t       *restrict xa,
                    cs_real_t             *restrict x,
                    cs_real_t             *restrict y)
{
  cs_int_t ii;
  double wt, cpu;
  int    run_id;
  long   n_ops, n_ops_glob;

  double test_sum = 0.0;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;

  if (n_runs < 1)
    return;

  n_ops = n_cells + n_faces*2;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (  cs_glob_mesh->n_g_cells
                  + cs_glob_mesh->n_g_i_faces*2);

  ms = cs_matrix_structure_create(type,
                                  true,
                                  n_cells,
                                  n_cells_ext,
                                  n_faces,
                                  cell_num,
                                  face_cell,
                                  halo,
                                  numbering);

  m = cs_matrix_create(ms);

  cs_matrix_set_coefficients(m,
                             sym_coeffs,
                             NULL,
                             da,
                             xa);

  /* Matrix.vector product */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    cs_matrix_vector_multiply(CS_PERIO_ROTA_COPY,
                              m,
                              x,
                              y);
    test_sum += y[n_cells-1];
#if 0
    for (ii = 0; ii < n_cells; ii++)
      bft_printf("y[%d] = %12.4f\n", ii, y[ii]);
#endif
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Matrix.vector product (%s)\n"
               "---------------------\n"), _(type_name));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

  /* Local timing in parallel mode */

  if (cs_glob_n_ranks > 1) {

    test_sum = 0.0;

    _timer_start(&wt, &cpu);

    for (run_id = 0; run_id < n_runs; run_id++) {
      cs_matrix_vector_multiply_nosync(m,
                                       x,
                                       y);
      test_sum += y[n_cells-1];
    }

    _timer_stop(n_runs, &wt, &cpu);

    bft_printf(_("\n"
                 "Local matrix.vector product (%s)\n"
                 "---------------------------\n"), _(type_name));

    bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
               n_runs, test_sum);

    _print_stats(n_ops, n_ops_glob, wt, cpu);

  }

  /* Combined matrix.vector product: alpha.A.x + beta.y */

  test_sum = 0.0;
  for (ii = 0; ii < n_cells_ext; y[ii++] = 0.0);

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    cs_matrix_alpha_a_x_p_beta_y(CS_PERIO_ROTA_COPY,
                                 0.5,
                                 0.0,
                                 m,
                                 x,
                                 y);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Matrix.vector product alpha.A.x + beta.y (%s)\n"
               "----------------------------------------\n"), _(type_name));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

  /* (Matrix - diagonal).vector product (with diagonal structure) */

  cs_matrix_set_coefficients(m,
                             sym_coeffs,
                             NULL,
                             NULL,
                             xa);

  test_sum = 0.0;

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    cs_matrix_vector_multiply(CS_PERIO_ROTA_COPY,
                              m,
                              x,
                              y);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "(Matrix-diagonal).vector product (%s)\n"
               "--------------------------------\n"), _(type_name));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

  /* (Matrix - diagonal).vector product */

  cs_matrix_destroy(&m);
  cs_matrix_structure_destroy(&ms);

  n_ops = n_faces*2;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (  cs_glob_mesh->n_g_cells
                  + cs_glob_mesh->n_g_i_faces*2);

  ms = cs_matrix_structure_create(type,
                                  false,
                                  n_cells,
                                  n_cells_ext,
                                  n_faces,
                                  cell_num,
                                  face_cell,
                                  halo,
                                  numbering);

  m = cs_matrix_create(ms);

  cs_matrix_set_coefficients(m,
                             sym_coeffs,
                             NULL,
                             NULL,
                             xa);

  test_sum = 0.0;

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    cs_matrix_vector_multiply(CS_PERIO_ROTA_COPY,
                              m,
                              x,
                              y);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "(Matrix without diagonal).vector product (%s)\n"
               "----------------------------------------\n"), _(type_name));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

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
 *   n_runs      <-- number of operation runs
 *   n_cells     <-- number of cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   face_cell   <-- face -> cells connectivity (1 to n)
 *   xa          <-- extradiagonal values
 *   x           <-> vector
 *   y           --> vector
 *----------------------------------------------------------------------------*/

static void
_sub_matrix_vector_test(int                  n_runs,
                        cs_int_t             n_cells,
                        cs_int_t             n_cells_ext,
                        cs_int_t             n_faces,
                        const cs_int_t      *face_cell,
                        const cs_real_t     *restrict xa,
                        cs_real_t           *restrict x,
                        cs_real_t           *restrict y)
{
  cs_int_t  jj;
  double wt, cpu;
  int    run_id;
  long   n_ops, n_ops_glob;
  double *ya = NULL;

  double test_sum = 0.0;

  if (n_runs < 1)
    return;

  n_ops = n_faces*2;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (  cs_glob_mesh->n_g_cells
                  + cs_glob_mesh->n_g_i_faces*2);

  for (jj = 0; jj < n_cells_ext; jj++)
    y[jj] = 0.0;

  /* Matrix.vector product, variant 0 */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    _mat_vec_exdiag_native(n_faces, face_cell, xa, x, y);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Matrix.vector product, extradiagonal part, variant 0\n"
               "---------------------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

  for (jj = 0; jj < n_cells_ext; jj++)
    y[jj] = 0.0;

  test_sum = 0.0;

  /* Matrix.vector product, variant 1 */

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    _mat_vec_exdiag_native_v1(n_faces, face_cell, xa, x, y);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  bft_printf(_("\n"
               "Matrix.vector product, extradiagonal part, variant 1\n"
               "---------------------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

  /* Matrix.vector product, contribute to faces only */

  BFT_MALLOC(ya, n_faces, cs_real_t);
  for (jj = 0; jj < n_faces; jj++)
    ya[jj] = 0.0;

  test_sum = 0.0;

  _timer_start(&wt, &cpu);

  for (run_id = 0; run_id < n_runs; run_id++) {
    _mat_vec_exdiag_part_p1(n_faces, face_cell, xa, x, ya);
    test_sum += y[n_cells-1];
  }

  _timer_stop(n_runs, &wt, &cpu);

  BFT_FREE(ya);

  bft_printf(_("\n"
               "Matrix.vector product, face values only\n"
               "---------------------\n"));

  bft_printf(_("  (calls: %d;  test sum: %12.5f)\n"),
             n_runs, test_sum);

  _print_stats(n_ops, n_ops_glob, wt, cpu);

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

  int n_runs;
  size_t ii;

  cs_real_t *x1 = NULL, *x2 = NULL;
  cs_real_t *y1 = NULL, *y2 = NULL;
  cs_real_t *da = NULL, *xa = NULL;

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

  if (CS_MEM_ALIGN > 0) {

    BFT_MEMALIGN(y1, CS_MEM_ALIGN, n_cells_ext, cs_real_t);
    BFT_MEMALIGN(y2, CS_MEM_ALIGN, n_cells_ext, cs_real_t);

  }
  else {

    BFT_MALLOC(y1, n_cells_ext, cs_real_t);
    BFT_MALLOC(y2, n_cells_ext, cs_real_t);

  }

  BFT_MALLOC(da, n_cells_ext, cs_real_t);
  BFT_MALLOC(xa, n_faces*2, cs_real_t);

  for (ii = 0; ii < n_cells_ext; ii++) {
    da[ii] = 1.0 + ii/n_cells_ext;
    xa[ii] = mesh_v->cell_cen[ii*3 + 1];
  }

  for (ii = 0; ii < n_faces; ii++) {
    xa[ii*2] = 0.5*(1.0 + ii/n_faces);
    xa[ii*2 + 1] = -0.5*(1.0 + ii/n_faces);
  }

  /* Run tests */
  /*-----------*/

  bft_printf(_("\n"
               "Benchmark mode activated\n"
               "========================\n"));

  /* Dot product test */
  /*------------------*/

  n_runs = (mpi_trace_mode) ? 1 : 10000;

  _dot_product_1(0, n_runs, n_cells, x1, x1);
  _dot_product_1(0, n_runs, n_cells, x1, x2);
  _dot_product_2(n_runs, n_cells, x1, x2);
  _axpy_(n_runs, n_cells, x1, y1);

  _division_test(n_runs/5, n_cells);
  _sqrt_test(n_runs/10, n_cells);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    n_runs = (mpi_trace_mode) ? 1 : 10000;

    _dot_product_1(1, n_runs, n_cells, x1, x1);
    _dot_product_1(1, n_runs, n_cells, x1, x2);

  }

#endif /* HAVE_MPI */

  /* Matrix test */
  /*-------------*/

  n_runs = (mpi_trace_mode) ? 0 : 500;

  if (!mpi_trace_mode) {

    /* Creation tests */

    n_runs = 2000;

    _matrix_creation_test(n_runs,
                          _("native"),
                          CS_MATRIX_NATIVE,
                          n_cells, n_cells_ext, n_faces,
                          mesh->global_cell_num, mesh->i_face_cells,
                          mesh->halo,
                          mesh->i_face_numbering);

    n_runs = 300;

    _matrix_creation_test(n_runs,
                          _("CSR"),
                          CS_MATRIX_CSR,
                          n_cells, n_cells_ext, n_faces,
                          mesh->global_cell_num, mesh->i_face_cells,
                          mesh->halo,
                          mesh->i_face_numbering);

    _matrix_creation_test(n_runs,
                          _("symmetric CSR"),
                          CS_MATRIX_CSR_SYM,
                          n_cells, n_cells_ext, n_faces,
                          mesh->global_cell_num, mesh->i_face_cells,
                          mesh->halo,
                          mesh->i_face_numbering);

    /* Assignment tests */

    n_runs = 2000;

    _matrix_assignment_test(n_runs,
                            _("native"),
                            CS_MATRIX_NATIVE, false,
                            n_cells, n_cells_ext, n_faces,
                            mesh->global_cell_num, mesh->i_face_cells,
                            mesh->halo, mesh->i_face_numbering, da, xa);

    _matrix_assignment_test(n_runs,
                            _("native, sym coeffs"),
                            CS_MATRIX_NATIVE, true,
                            n_cells, n_cells_ext, n_faces,
                            mesh->global_cell_num, mesh->i_face_cells,
                            mesh->halo, mesh->i_face_numbering, da, xa);

    n_runs = 300;

    _matrix_assignment_test(n_runs,
                            _("CSR"),
                            CS_MATRIX_CSR, false,
                            n_cells, n_cells_ext, n_faces,
                            mesh->global_cell_num, mesh->i_face_cells,
                            mesh->halo, mesh->i_face_numbering, da, xa);

    _matrix_assignment_test(n_runs,
                            _("CSR, sym coeffs"),
                            CS_MATRIX_CSR, true,
                            n_cells, n_cells_ext, n_faces,
                            mesh->global_cell_num, mesh->i_face_cells,
                            mesh->halo, mesh->i_face_numbering, da, xa);

    _matrix_assignment_test(n_runs,
                            _("symmetric CSR"),
                            CS_MATRIX_CSR_SYM, true,
                            n_cells, n_cells_ext, n_faces,
                            mesh->global_cell_num, mesh->i_face_cells,
                            mesh->halo, mesh->i_face_numbering, da, xa);

  }

  /* Matrix.vector tests */

  n_runs = (mpi_trace_mode) ? 1 : 1000;

  _matrix_vector_test(n_runs,
                      _("native"),
                      CS_MATRIX_NATIVE, false,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y1);

  _matrix_vector_test(n_runs,
                      _("native, sym coeffs"),
                      CS_MATRIX_NATIVE, true,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y1);

  _matrix_vector_test(n_runs,
                      _("CSR"),
                      CS_MATRIX_CSR, false,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y1);

  _matrix_vector_test(n_runs,
                      _("symmetric CSR"),
                      CS_MATRIX_CSR_SYM, true,
                      n_cells, n_cells_ext, n_faces,
                      mesh->global_cell_num, mesh->i_face_cells, mesh->halo,
                      mesh->i_face_numbering, da, xa, x1, y1);

  _sub_matrix_vector_test(n_runs,
                          n_cells,
                          n_cells_ext,
                          n_faces,
                          mesh->i_face_cells,
                          xa,
                          x1,
                          y1);

  /* Free working arrays */
  /*---------------------*/

  BFT_FREE(x1);
  BFT_FREE(x2);

  BFT_FREE(y1);
  BFT_FREE(y2);

  BFT_FREE(da);
  BFT_FREE(xa);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
