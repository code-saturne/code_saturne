/*============================================================================
 * Low-level operator benchmarking.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined(HAVE_MKL)
#include <mkl_cblas.h>
#include <mkl_spblas.h>

#elif defined(HAVE_CBLAS)
#include <cblas.h>

#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_matrix.h"
#include "cs_matrix_assembler.h"
#include "cs_matrix_default.h"
#include "cs_matrix_tuning.h"
#include "cs_timer.h"

#if defined(HAVE_HYPRE)
#include "cs_matrix_hypre.h"
#include "cs_sles_hypre.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_matrix_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_benchmark.h"
#include "cs_benchmark_matrix.h"

#if defined(HAVE_CUDA)
#include "cs_benchmark_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#if defined(HAVE_CBLAS) || defined (HAVE_MKL)
#define HAVE_BLAS 1
#endif

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static const char *_matrix_operation_name[CS_MATRIX_N_FILL_TYPES][2]
  = {{"y <- A.x",
      "y <- (A-D).x"},
     {"Symmetric y <- A.x",
      "Symmetric y <- (A-D).x"},
     {"Block diagonal y <- A.x",
      "Block diagonal y <- (A-D).x"},
     {"Block 6 diagonal y <- A.x",
      "Block 6 diagonal y <- (A-D).x"},
     {"Block diagonal symmetric y <- A.x",
      "Block diagonal symmetric y <- (A-D).x"},
     {"Block y <- A.x",
      "Block y <- (A-D).x"}};

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
                  "  N ops:       %12ld\n"
                  "  Wall clock:  %12.5e\n"
                  "  GFLOPS:      %12.5e\n",
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
         "               Mean         Min          Max          Total\n"
         "  N ops:       %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e\n",
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

    else
      cs_log_printf
        (CS_LOG_PERFORMANCE,
         "               Mean         Min          Max          Total"
         "        Single\n"
         "  N ops:       %12ld %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e %12.5e\n",
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
 * Measure matrix.vector product extradiagonal terms related performance
 * (symmetric matrix case).
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_exdiag_native(cs_lnum_t            n_faces,
                       const cs_lnum_2_t   *face_cell,
                       const cs_real_t     *restrict xa,
                       cs_real_t           *restrict x,
                       cs_real_t           *restrict y)
{
  cs_lnum_t  ii, jj, face_id;

  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)face_cell;

  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = *face_cel_p++;
    jj = *face_cel_p++;
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
_mat_vec_exdiag_native_v1(cs_lnum_t            n_faces,
                          const cs_lnum_2_t   *face_cell,
                          const cs_real_t     *restrict xa,
                          cs_real_t           *restrict x,
                          cs_real_t           *restrict y)
{
  cs_lnum_t  ii, ii_prev, kk, face_id, kk_max;
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

  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)face_cell;

  for (face_id = 0;
       face_id < n_faces;
       face_id += l1_cache_size) {

    kk_max = CS_MIN((n_faces - face_id), l1_cache_size);

    /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

    ii = face_cel_p[0];
    ii_prev = ii;
    y_it_prev = y[ii_prev] + xa[face_id] * x[face_cel_p[1]];

    for (kk = 1; kk < kk_max; ++kk) {
      ii = face_cel_p[2*kk];
      /* y[ii] += xa[face_id+kk] * x[jj]; */
      if(ii == ii_prev) {
        y_it = y_it_prev;
      }
      else {
        y_it = y[ii];
        y[ii_prev] = y_it_prev;
      }
      ii_prev = ii;
      y_it_prev = y_it + xa[face_id+kk] * x[face_cel_p[2*kk+1]];
    }
    y[ii] = y_it_prev;

    /* sub-loop to compute y[ii] += xa[face_id] * x[jj] */

    for (kk = 0; kk < kk_max; ++kk) {
      y[face_cel_p[2*kk+1]]
        += xa[face_id+kk] * x[face_cel_p[2*kk]];
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
 *   face_cell       <-- face -> cells connectivity
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   ya              <-> vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_exdiag_part_p1(cs_lnum_t            n_faces,
                        const cs_lnum_2_t   *face_cell,
                        const cs_real_t     *restrict xa,
                        cs_real_t           *restrict x,
                        cs_real_t           *restrict ya)
{
  cs_lnum_t  ii, jj, face_id;

  const cs_lnum_t *restrict face_cel_p
    = (const cs_lnum_t *restrict)face_cell;

  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = *face_cel_p++;
    jj = *face_cel_p++;
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
 *   face_cell   <-- face -> cells connectivity
 *   xa          <-- extradiagonal values
 *   x           <-> vector
 *   y           --> vector
 *----------------------------------------------------------------------------*/

static void
_sub_matrix_vector_test(double               t_measure,
                        cs_lnum_t            n_cells,
                        cs_lnum_t            n_cells_ext,
                        cs_lnum_t            n_faces,
                        const cs_lnum_2_t   *face_cell,
                        const cs_real_t     *restrict xa,
                        cs_real_t           *restrict x,
                        cs_real_t           *restrict y)
{
  cs_lnum_t  jj;
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
                "\n"
                "Matrix.vector product, extradiagonal part, variant 0\n"
                "---------------------\n");

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  (calls: %d;  test sum: %12.5f)\n",
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
                "\n"
                "Matrix.vector product, extradiagonal part, variant 1\n"
                "---------------------\n");

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  (calls: %d;  test sum: %12.5f)\n",
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

  /* Matrix.vector product, CUDA variant */

#if (HAVE_CUDA)

  const cs_lnum_2_t *__restrict__ d_face_cell
    = (const cs_lnum_2_t *)cs_get_device_ptr((void *)face_cell);
  const cs_real_t *__restrict__ d_xa
    = (const cs_real_t *)cs_get_device_ptr((void *)xa);
  const cs_real_t *__restrict__ d_x
    = (const cs_real_t *)cs_get_device_ptr((void *)x);
  cs_real_t *__restrict__ d_y
    = (cs_real_t *)cs_get_device_ptr((void *)y);

  //cs_sync_h2d(face_cell);
  cs_sync_h2d(xa);
  cs_sync_h2d(x);
  cs_sync_h2d(y);

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
      cs_mat_vec_exdiag_native_sym_cuda(n_faces, d_face_cell, d_xa, d_x, d_y);
      test_sum += y[n_cells-1]*test_sum_mult;
      run_id++;
    }
    wt1 = cs_timer_wtime();
    if (wt1 - wt0 < t_measure)
      n_runs *= 2;
  }

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Matrix.vector product, extradiagonal part, CUDA variant\n"
                "---------------------\n");

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  (calls: %d;  test sum: %12.5f)\n",
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);

#endif /* (HAVE_CUDA) */

  /* Matrix.vector product, contribute to faces only */

  /* n_faces*2 nonzeroes, n_row_elts multiplications */

  n_ops = n_faces*2;

  if (cs_glob_n_ranks == 1)
    n_ops_glob = n_ops;
  else
    n_ops_glob = (cs_glob_mesh->n_g_i_faces*2);

  CS_MALLOC_HD(ya, n_faces, cs_real_t, cs_alloc_mode);
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

  CS_FREE_HD(ya);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Matrix.vector product, face values only\n"
                "---------------------\n");

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  (calls: %d;  test sum: %12.5f)\n",
                n_runs, test_sum);

  _print_stats(n_runs, n_ops, n_ops_glob, wt1 - wt0);
}

/*----------------------------------------------------------------------------
 * Copy array to reference for matrix computation check.
 *
 * parameters:
 *   n_elts      <-- number values to compare
 *   y           <-- array to copare or copy
 *   yr          <-- reference array
 *
 * returns:
 *   maximum difference between values
 *----------------------------------------------------------------------------*/

static double
_matrix_check_compare(cs_lnum_t        n_elts,
                      const cs_real_t  y[],
                      cs_real_t        yr[])
{
  double dmax = 0.0;

  for (cs_lnum_t ii = 0; ii < n_elts; ii++) {
    double d = CS_ABS(y[ii] - yr[ii]);
    if (d > dmax)
      dmax = d;
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double dmaxg;
    MPI_Allreduce(&dmax, &dmaxg, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
    dmax = dmaxg;
  }

#endif

  return dmax;
}

/*----------------------------------------------------------------------------
 * Check local matrix.vector product operations using matrix assembler
 *
 * parameters:
 *   n_rows      <-- local number of rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of (undirected) graph edges
 *   cell_num    <-- optional global cell numbers (1 to n), or NULL
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   halo        <-- cell halo structure
 *----------------------------------------------------------------------------*/

static void
_matrix_check_asmb(cs_lnum_t              n_rows,
                   cs_lnum_t              n_cols_ext,
                   cs_lnum_t              n_edges,
                   const cs_lnum_2_t     *edges,
                   const cs_halo_t       *halo)
{
  cs_real_t  *da = NULL, *xa = NULL, *x = NULL, *y = NULL;
  cs_real_t  *yr0 = NULL;
  cs_lnum_t a_block_size = 3;
  cs_lnum_t a_block_stride = a_block_size*a_block_size;

  cs_matrix_fill_type_t f_type[]
    = {CS_MATRIX_SCALAR,           /* Simple scalar matrix */
       CS_MATRIX_BLOCK_D};

  const char *t_name[] = {"general assembly",
                          "local rows assembly",
                          "assembly from shared"};

  /* Allocate and initialize  working arrays */

  if (CS_MEM_ALIGN > 0) {
    BFT_MEMALIGN(x, CS_MEM_ALIGN, n_cols_ext*a_block_size, cs_real_t);
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n_cols_ext*a_block_size, cs_real_t);
    BFT_MEMALIGN(yr0, CS_MEM_ALIGN, n_cols_ext*a_block_size, cs_real_t);
  }
  else {
    BFT_MALLOC(x, n_cols_ext*a_block_size, cs_real_t);
    BFT_MALLOC(y, n_cols_ext*a_block_size, cs_real_t);
    BFT_MALLOC(yr0, n_cols_ext*a_block_size, cs_real_t);
  }

  BFT_MALLOC(da, n_cols_ext*a_block_stride, cs_real_t);
  BFT_MALLOC(xa, n_edges*2*a_block_stride, cs_real_t);

  cs_gnum_t *cell_gnum = NULL;
  BFT_MALLOC(cell_gnum, n_cols_ext, cs_gnum_t);
  if (cs_glob_mesh->global_cell_num != NULL) {
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      cell_gnum[ii] = cs_glob_mesh->global_cell_num[ii];
  }
  else {
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      cell_gnum[ii] = ii+1;
  }
  if (halo != NULL)
    cs_halo_sync_untyped(halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_gnum_t),
                         cell_gnum);

  /* Global cell ids, based on range/scan */

  cs_gnum_t l_range[2] = {0, n_rows};
  cs_gnum_t n_g_rows = n_rows;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t g_shift;
    cs_gnum_t l_shift = n_rows;
    MPI_Scan(&l_shift, &g_shift, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&l_shift, &n_g_rows, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    l_range[0] = g_shift - l_shift;
    l_range[1] = g_shift;
  }
#endif

  cs_gnum_t *r_g_id;
  BFT_MALLOC(r_g_id, n_cols_ext, cs_gnum_t);
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    r_g_id[ii] = ii + l_range[0];
  if (halo != NULL)
    cs_halo_sync_untyped(halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_gnum_t),
                         r_g_id);

  /* Loop on fill options */

  for (int f_id = 0; f_id < 2; f_id++) {

    const cs_lnum_t d_block_size
      = (f_type[f_id] >= CS_MATRIX_BLOCK_D) ? a_block_size : 1;
    const cs_lnum_t stride = d_block_size;
    const cs_lnum_t sd = stride*stride; /* for current fill types */
    const cs_lnum_t se = 1;             /* for current fill types */

    /* Initialize arrays; we need to be careful here, so that
       the array values are consisten across MPI ranks.
       This requires using a specific initialiation for each fill type. */

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cols_ext; ii++) {
      cs_gnum_t jj = (cell_gnum[ii] - 1)*sd;
      for (cs_lnum_t kk = 0; kk < sd; kk++) {
        da[ii*sd+kk] = 1.0 + cos(jj*sd+kk);
      }
    }

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
      cs_lnum_t i0 = edges[ii][0];
      cs_lnum_t i1 = edges[ii][1];
      cs_gnum_t j0 = (cell_gnum[i0] - 1)*se;
      cs_gnum_t j1 = (cell_gnum[i1] - 1)*se;
      for (cs_lnum_t kk = 0; kk < se; kk++) {
        xa[(ii*se+kk)*2]
          = 0.5*(0.45 + cos(j0*se+kk) + cos(j1*se+kk));
        xa[(ii*se+kk)*2 + 1]
          = -0.5*(0.45 + cos(j0*se+kk) + cos(j1*se+kk));
      }
    }

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cols_ext; ii++) {
      cs_gnum_t jj = (cell_gnum[ii] - 1)*stride;
      for (cs_lnum_t kk = 0; kk < stride; kk++)
        x[ii*stride+kk] = sin(jj*stride+kk);
    }

    /* Reference */

    cs_matrix_vector_native_multiply(false,  /* symmetric */
                                     d_block_size,
                                     1,      /* extra diag block size */
                                     -1,     /* field id or -1 */
                                     da,
                                     xa,
                                     x,
                                     yr0);

    /* Test for matrix assembler (for MSR case) */

    const cs_lnum_t block_size = 800;
    cs_gnum_t g_row_id[800];
    cs_gnum_t g_col_id[800];
    cs_real_t val[1600];

    cs_matrix_assembler_t *ma = NULL;

    /* 3 variant construction methods, 2 coefficient methods */

    const char *ma_name[] = {"distributed contribution assember",
                             "local rows assembler",
                             "shared index assembler"};

    for (int c_id = 0; c_id < 3; c_id++) {

      /* Matrix created from shared index may not always handle
         periodic elements in the same manner */

      if (halo != NULL) {
        if (halo->n_transforms > 0 && c_id == 2)
          continue;
      }

      if (c_id < 2) {

        ma = cs_matrix_assembler_create(l_range, true);

        /* Add connectivities */

        cs_matrix_assembler_add_g_ids(ma, n_rows, r_g_id, r_g_id);

        if (c_id == 0) { /* Rank contributions through global edges */
          cs_lnum_t jj = 0;
          for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
            cs_lnum_t i0 = edges[ii][0];
            cs_lnum_t i1 = edges[ii][1];
            g_row_id[jj] = r_g_id[i0];
            g_col_id[jj] = r_g_id[i1];
            jj++;
            g_row_id[jj] = r_g_id[i1];
            g_col_id[jj] = r_g_id[i0];
            jj++;
            if (jj >= block_size - 1) {
              cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
              jj = 0;
            }
          }
          cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
        }
        else { /* Rank contributions are local */
          cs_lnum_t jj = 0;
          for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
            cs_lnum_t i0 = edges[ii][0];
            cs_lnum_t i1 = edges[ii][1];
            if (i0 < n_rows) {
              g_row_id[jj] = r_g_id[i0];
              g_col_id[jj] = r_g_id[i1];
              jj++;
            }
            if (i1 < n_rows) {
              g_row_id[jj] = r_g_id[i1];
              g_col_id[jj] = r_g_id[i0];
              jj++;
            }
            if (jj >= block_size - 1) {
              cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
              jj = 0;
            }
          }
          cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
        }

        cs_matrix_assembler_compute(ma);
      }
      else {
        const cs_mesh_adjacencies_t  *madj = cs_glob_mesh_adjacencies;
        ma = cs_matrix_assembler_create_from_shared(n_rows,
                                                    true,
                                                    madj->cell_cells_idx,
                                                    madj->cell_cells,
                                                    halo);
      }

      if (f_id == 0) /* no need to log multiple identical assemblers */
        cs_matrix_assembler_log_rank_counts(ma, CS_LOG_DEFAULT, ma_name[c_id]);

      cs_matrix_structure_t  *ms
        = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      for (int m_type_idx = 0; m_type_idx < 3; m_type_idx++) {

        cs_matrix_t  *m = NULL;

        switch (m_type_idx) {
        case 0:
          m = cs_matrix_create(ms);
          break;
        case 1:
#if defined(HAVE_HYPRE)
          {
            int device_id = cs_get_device_id();
            int use_device = (device_id < 0) ? 0 : 1;
            m = cs_matrix_create(ms);
            cs_matrix_set_type_hypre(m, use_device);
          }
#else
          continue;
#endif
          break;
        case 2:
#if defined(HAVE_PETSC)
          {
            m = cs_matrix_create(ms);
            cs_matrix_set_type_petsc(m, 0);
          }
#else
          continue;
#endif
          break;
        default:
          continue;
        }

        cs_matrix_assembler_values_t *mav = NULL;

        if (f_type[f_id] == CS_MATRIX_SCALAR)
          mav = cs_matrix_assembler_values_init(m, 1, 1);
        else if (f_type[f_id] == CS_MATRIX_BLOCK_D)
          mav = cs_matrix_assembler_values_init(m, d_block_size, 1);

        cs_matrix_assembler_values_add_g(mav, n_rows,
                                         r_g_id, r_g_id, da);

        if (c_id == 0) { /* Rank contributions through global edges */
          cs_lnum_t jj = 0;
          for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
            cs_lnum_t i0 = edges[ii][0];
            cs_lnum_t i1 = edges[ii][1];
            g_row_id[jj] = r_g_id[i0];
            g_col_id[jj] = r_g_id[i1];
            if (i0 < n_rows && i1 < n_rows)
              val[jj] = xa[ii*2];
            else
              val[jj] = xa[ii*2]*0.5; /* count half contribution twice */
            jj++;
            g_row_id[jj] = r_g_id[i1];
            g_col_id[jj] = r_g_id[i0];
            if (i0 < n_rows && i1 < n_rows)
              val[jj] = xa[ii*2+1];
            else
              val[jj] = xa[ii*2+1]*0.5;
            jj++;
            if (jj >= block_size - 1) {
              cs_matrix_assembler_values_add_g(mav, jj,
                                               g_row_id, g_col_id, val);
              jj = 0;
            }
          }
          for (cs_lnum_t ii = 0; ii < jj; ii++) {
          }

          cs_matrix_assembler_values_add_g(mav, jj, g_row_id, g_col_id, val);
        }
        else { /* Rank contributions are local */
          cs_lnum_t jj = 0;
          for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
            cs_lnum_t i0 = edges[ii][0];
            cs_lnum_t i1 = edges[ii][1];
            if (i0 < n_rows) {
              g_row_id[jj] = r_g_id[i0];
              g_col_id[jj] = r_g_id[i1];
              val[jj] = xa[ii*2];
              jj++;
            }
            if (i1 < n_rows) {
              g_row_id[jj] = r_g_id[i1];
              g_col_id[jj] = r_g_id[i0];
              val[jj] = xa[ii*2+1];
              jj++;
            }
            if (jj >= block_size - 1) {
              cs_matrix_assembler_values_add_g(mav, jj,
                                               g_row_id, g_col_id, val);
              jj = 0;
            }
          }
          cs_matrix_assembler_values_add_g(mav, jj, g_row_id, g_col_id, val);
        }

        cs_matrix_assembler_values_finalize(&mav);

        cs_matrix_vector_multiply(m, x, y);

        cs_matrix_release_coefficients(m);

        double dmax = _matrix_check_compare(n_rows*stride, y, yr0);
        bft_printf("\n%s (%s)\n",
                   _matrix_operation_name[f_type[f_id]][0],
                   cs_matrix_get_type_name(m));
        bft_printf("  %-32s : %12.5e\n",
                   t_name[c_id],
                   dmax);
        bft_printf_flush();

        cs_matrix_destroy(&m);

      }  /* End of loop on matrix types */

      cs_matrix_structure_destroy(&ms);

      cs_matrix_assembler_destroy(&ma);

    } /* end of loop on construction method */

  } /* end of loop on fill types */

  BFT_FREE(r_g_id);
  BFT_FREE(cell_gnum);

  BFT_FREE(yr0);

  BFT_FREE(y);
  BFT_FREE(x);

  BFT_FREE(xa);
  BFT_FREE(da);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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

  double t_measure = (mpi_trace_mode) ? -1.0 : 0.8;

  cs_real_t *x = NULL, *y = NULL;
  cs_real_t *da = NULL, *xa = NULL;

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mesh_v = cs_glob_mesh_quantities;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);

  size_t n_cells = mesh->n_cells;
  size_t n_cells_ext = mesh->n_cells_with_ghosts;
  size_t n_faces = mesh->n_i_faces;

  int                    n_fill_types_nsym = 4;
  int                    n_fill_types_sym = 2;
  cs_matrix_fill_type_t  fill_types_nsym[] = {CS_MATRIX_SCALAR,
                                              CS_MATRIX_BLOCK_D,
                                              CS_MATRIX_BLOCK_D_66,
                                              CS_MATRIX_BLOCK};
  cs_matrix_fill_type_t  fill_types_sym[] = {CS_MATRIX_SCALAR_SYM,
                                             CS_MATRIX_BLOCK_D_SYM};

  cs_matrix_initialize();

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Benchmark mode activated\n"
                "========================\n");

#if defined(HAVE_HYPRE)
  /* Force HYPRE initialization */
  void *hypre_sles
    = cs_sles_hypre_create(CS_SLES_HYPRE_NONE, CS_SLES_HYPRE_NONE, NULL, NULL);
#endif

  /* Run some feature tests */
  /*------------------------*/

  _matrix_check_asmb(n_cells,
                     n_cells_ext,
                     n_faces,
                     i_face_cells,
                     mesh->halo);

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  CS_MALLOC_HD(x, n_cells_ext, cs_real_t, cs_alloc_mode);

  for (ii = 0; ii < n_cells_ext; ii++)
    x[ii] = mesh_v->cell_cen[ii*3];

  CS_MALLOC_HD(y, n_cells_ext, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(da, n_cells_ext, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xa, n_faces*2, cs_real_t, cs_alloc_mode);

  for (ii = 0; ii < n_cells_ext; ii++)
    da[ii] = 1.0;

  for (ii = 0; ii < n_faces; ii++) {
    xa[ii*2] = 0.5;
    xa[ii*2 + 1] = -0.5;
  }

  /* Call matrix tuning */
  /*--------------------*/

  /* Enter tuning phase */

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "General timing for matrices\n"
                "===========================\n");

  cs_benchmark_matrix(t_measure,
                      0,
                      n_fill_types_nsym,
                      NULL,
                      fill_types_nsym,
                      n_cells,
                      n_cells_ext,
                      n_faces,
                      i_face_cells,
                      mesh->halo,
                      mesh->i_face_numbering);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Timing for symmetric matrices\n"
                "=============================\n");

  cs_benchmark_matrix(t_measure,
                      0,
                      n_fill_types_sym,
                      NULL,
                      fill_types_sym,
                      n_cells,
                      n_cells_ext,
                      n_faces,
                      i_face_cells,
                      mesh->halo,
                      mesh->i_face_numbering);

  _sub_matrix_vector_test(t_measure,
                          n_cells,
                          n_cells_ext,
                          n_faces,
                          i_face_cells,
                          xa,
                          x,
                          y);

  cs_matrix_finalize();

  cs_mesh_adjacencies_finalize();

  cs_log_separator(CS_LOG_PERFORMANCE);

#if defined(HAVE_HYPRE)
  cs_sles_hypre_destroy(&hypre_sles);
#endif

  /* Free working arrays */
  /*---------------------*/

  CS_FREE_HD(x);
  CS_FREE_HD(y);

  CS_FREE_HD(da);
  CS_FREE_HD(xa);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
