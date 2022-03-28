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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <sys/time.h>
#include <unistd.h>

/* For the Intel MKL library, function prototypes are defined in mkl_cblas.h,
   with standard legacy C BLAS names */

#if defined(HAVE_MKL)
#include <mkl_cblas.h>

#elif defined(HAVE_ATLAS) || defined(HAVE_CBLAS)
#include <cblas.h>

#endif

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

#include "cs_system_info.h"

#include "cs_blas.h"
#include "cs_defs.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/* Block size for superblock algorithm */
#define CS_SBLOCK_BLOCK_SIZE 60

#if defined(HAVE_ATLAS)
const char ext_blas_type[] = "ATLAS";
#elif defined(HAVE_MKL)
const char ext_blas_type[] = "MKL";
#elif defined(HAVE_BLAS)
const char ext_blas_type[] = "external";
#else
const char ext_blas_type[] = "";
#endif

/* Reduction algorithm name */
const char *reduce_name[] = {"superblock", "Kahan compensated sum"};

/* Size of arrays for tests */
static const int _n_sizes = 9;
static const cs_lnum_t _n_elts[]
  = {1000000, 200000, 50000, 20000, 10000, 3000, 1000, 500, 100};

static double _pi = 3.14159265358979323846;

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * False print of a message to standard output for discarded logs
 *----------------------------------------------------------------------------*/

static int
_bft_printf_null(const char  *format,
                 va_list      arg_ptr)
{
  CS_UNUSED(format);
  CS_UNUSED(arg_ptr);

  return 0;
}

/*----------------------------------------------------------------------------
 * Analysis of environment variables to determine
 * if we require MPI, and initialization if necessary.
 *----------------------------------------------------------------------------*/

static void
_mpi_init(void)
{
  int flag = 0;
  bool use_mpi = false;

#if   defined(__bg__) || defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  use_mpi = true;

#elif defined(MPICH2)
  if (getenv("PMI_RANK") != NULL)
    use_mpi = true;

#elif defined(MPICH_NAME)

  if (getenv("GMPI_ID") != NULL) /* In case we are using MPICH-GM */
    use_mpi = true;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL)         /* OpenMPI 1.2 */
    use_mpi = true;
  else if (getenv("OMPI_COMM_WORLD_RANK") != NULL)    /* OpenMPI 1.3 and above */
    use_mpi = true;

#endif /* Tests for known MPI variants */

  /* If we have determined from known MPI environment variables
     of command line arguments that we are running under MPI,
     initialize MPI */

  if (use_mpi == true) {

    MPI_Initialized(&flag);

    if (!flag) {
#if defined(MPI_VERSION) && (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
      int mpi_threads;
      MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &mpi_threads);
#else
      MPI_Init(NULL, NULL);
#endif
    }

    cs_glob_mpi_comm = MPI_COMM_WORLD;
    MPI_Comm_size(cs_glob_mpi_comm, &cs_glob_n_ranks);
    MPI_Comm_rank(cs_glob_mpi_comm, &cs_glob_rank_id);

    if (cs_glob_rank_id > 0)
      bft_printf_proxy_set(_bft_printf_null);

  }

}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/* Return the dot product of 2 vectors: x.y                                   */
/*----------------------------------------------------------------------------*/

static double
_ddot_canonical(cs_lnum_t      n,
                const double  *x,
                const double  *y)
{
  cs_lnum_t  i;
  double     s = 0;

  if (n < 1)
    return s;

# pragma omp parallel for reduction(+: s)
  for (i = 0; i < n; i++)
    s += (x[i] * y[i]);

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute blocks sizes for superblock algorithm.
 *
 * \param[in]   n                 size of array
 * \param[in]   block_size        block size
 * \param[out]  n_sblocks         number of superblocks
 * \param[out]  blocks_in_sblocks number of blocks per superblock
 */
/*----------------------------------------------------------------------------*/

static inline void
_sbloc_sizes(cs_lnum_t   n,
             cs_lnum_t   block_size,
             cs_lnum_t  *n_sblocks,
             cs_lnum_t  *blocks_in_sblocks)
{
  cs_lnum_t n_blocks = (n + block_size - 1) / block_size;
  *n_sblocks = (n_blocks > 1) ? sqrt(n_blocks) : 1;

  cs_lnum_t n_b = block_size * *n_sblocks;
  *blocks_in_sblocks = (n + n_b - 1) / n_b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the dot product of 2 vectors: x.y
 *        using a superblock algorithm with different x and y precision.
 *
 * \param[in]  n  size of arrays x and y
 * \param[in]  x  array of floating-point values
 * \param[in]  y  array of floating-point values
 *
 * \return  dot product
 */
/*----------------------------------------------------------------------------*/

static double
_cs_dot_superblock_m(cs_lnum_t         n,
                     const double     *x,
                     const float      *y)
{
  double dot = 0.0;

# pragma omp parallel reduction(+:dot) if (n > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    cs_parall_thread_range(n, sizeof(cs_real_t), &s_id, &e_id);

    const cs_lnum_t _n = e_id - s_id;
    const double *_x = x + s_id;
    const float  *_y = y + s_id;

    const cs_lnum_t block_size = CS_SBLOCK_BLOCK_SIZE;
    cs_lnum_t n_sblocks, blocks_in_sblocks;

    _sbloc_sizes(_n, block_size, &n_sblocks, &blocks_in_sblocks);

    for (cs_lnum_t sid = 0; sid < n_sblocks; sid++) {

      double sdot = 0.0;

      for (cs_lnum_t bid = 0; bid < blocks_in_sblocks; bid++) {
        cs_lnum_t start_id = block_size * (blocks_in_sblocks*sid + bid);
        cs_lnum_t end_id = block_size * (blocks_in_sblocks*sid + bid + 1);
        if (end_id > _n)
          end_id = _n;
        double cdot = 0.0;
        for (cs_lnum_t i = start_id; i < end_id; i++)
        cdot += _x[i]*_y[i];
        sdot += cdot;
      }

      dot += sdot;

    }

  }

  return dot;
}

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
    bft_printf("  N ops:       %12ld\n"
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
      bft_printf
        ("               Mean         Min          Max          Total\n"
         "  N ops:       %12ld %12ld %12ld %12ld\n"
         "  Wall clock:  %12.5e %12.5e %12.5e\n"
         "  GFLOPS:      %12.5e %12.5e %12.5e %12.5e\n",
         n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
         glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
         glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

    else
      bft_printf
        ("               Mean         Min          Max          Total"
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

  bft_printf_flush();
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
    bft_printf("  N ops:       %12ld\n"
               "  Wall clock:  %12.5e\n"
               "  GB/s:        %12.5e\n",
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

    bft_printf
      ("               Mean         Min          Max          Total\n"
       "  N ops:       %12ld %12ld %12ld %12ld\n"
       "  Wall clock:  %12.5e %12.5e %12.5e\n"
       "  GB/s:        %12.5e %12.5e %12.5e %12.5e\n",
       n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
       glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
       glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

  }

#endif

  bft_printf_flush();
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
_print_time_stats(long    n_runs,
                 long    n_ops,
                 double  wt)
{
  double fm = 1.0 * n_runs / (1.e9 * (CS_MAX(wt, 1)));

  if (cs_glob_n_ranks == 1)
    bft_printf("  N ops:       %12ld\n"
               "  Wall clock:  %12.5e\n"
               "  Ops/s:       %12.5e\n",
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

    bft_printf
      ("               Mean         Min          Max          Total\n"
       "  N ops:       %12ld %12ld %12ld %12ld\n"
       "  Wall clock:  %12.5e %12.5e %12.5e\n"
       "  Ops/s:       %12.5e %12.5e %12.5e %12.5e\n",
       n_ops_tot/cs_glob_n_ranks, n_ops_min, n_ops_max, n_ops_tot,
       glob_sum[0]/cs_glob_n_ranks, glob_min[0], glob_max[0],
       glob_sum[1]/cs_glob_n_ranks, glob_min[1], glob_max[1], n_ops_tot*fmg);

  }

#endif

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Simple dot product.
 *
 * parameters:
 *   t_measure <-- minimum time for each measure (< 0 for single pass)
 *   global    <-- 0 for local use, 1 for MPI sum
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_dot_product_1(double   t_measure,
               int      global)
{
  double wt0, wt1;
  int type_id, sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  double test_sum = 0.0;
  int _global = global;

  const char *type_name[] = {"X.Y", "X.X"};

  cs_real_t *x = NULL;
  cs_real_t *y = NULL;

  if (cs_glob_n_ranks == 1)
    _global = 0;

  for (type_id = 0; type_id < 2; type_id++) {

    /* First simple local x.x version */

#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) ||defined(HAVE_MKL)

    for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

      n = _n_elts[sub_id];
      n_ops = n*2 - 1;

      /* Realloc and initialize arrays for each test, as
         first touch may affect memory locality on some systems */

      BFT_MALLOC(x, n, double);
      if (type_id == 0)
        BFT_MALLOC(y, n, double);
      else
        y = x;

#     pragma omp parallel for
      for (ii = 0; ii < n; ii++) {
        x[ii] = (ii%10 + 1);
        y[ii] = (ii%10 + 1)*_pi;
      }

      test_sum = test_sum - floor(test_sum);

      wt0 = cs_timer_wtime(), wt1 = wt0;
      if (t_measure > 0)
        n_runs = 8;
      else
        n_runs = 1;
      run_id = 0;
      while (run_id < n_runs) {
        double test_sum_mult = 1.0/n_runs;
        while (run_id < n_runs) {
#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) ||defined(HAVE_MKL)
          double s1 = cblas_ddot(n, x, 1, y, 1);
#endif
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
        bft_printf("\n"
                   "Local dot product %s with %s BLAS (%d elts.)\n"
                   "-----------------\n",
                   type_name[type_id], ext_blas_type, (int)n);
      else
        bft_printf("\n"
                   "Global dot product %s with %s BLAS (%d elts/rank)\n"
                   "------------------\n",
                   type_name[type_id], ext_blas_type, (int)n);

      bft_printf("  (calls: %d)\n", n_runs);

      _print_stats(n_runs, n_ops, 0, wt1 - wt0);

      if (type_id == 0)
        BFT_FREE(y);
      BFT_FREE(x);

  }

#endif /* external BLAS */

    /* Local dot product */

    for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

      n = _n_elts[sub_id];
      n_ops = n*2 - 1;

      /* Realloc and initialize arrays for each test, as
         first touch may affect memory locality on some systems */

      BFT_MALLOC(x, n, double);
      if (type_id == 0)
        BFT_MALLOC(y, n, double);
      else
        y = x;

#     pragma omp parallel for
      for (ii = 0; ii < n; ii++) {
        x[ii] = (ii%10 + 1)*_pi;
        y[ii] = (ii%10 + 1);
      }

      test_sum = test_sum - floor(test_sum);

      for (cs_blas_reduce_t j = 0; j < 2; j++) {

        cs_blas_set_reduce_algorithm(j);

        wt0 = cs_timer_wtime(), wt1 = wt0;
        if (t_measure > 0)
          n_runs = 8;
        else
          n_runs = 1;
        run_id = 0;
        while (run_id < n_runs) {
          double test_sum_mult = 1.0/n_runs;
          while (run_id < n_runs) {
            double s1 = cs_dot(n, x, y);
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
          bft_printf("\n"
                     "Local dot product %s with %s (%d elts.)\n"
                     "-----------------\n",
                     type_name[type_id], reduce_name[j], (int)n);
        else
          bft_printf("\n"
                     "Global dot product %s with %s (%d elts/rank)\n"
                     "------------------\n",
                     type_name[type_id], reduce_name[j], (int)n);

        bft_printf("  (calls: %d)\n", n_runs);

        _print_stats(n_runs, n_ops, 0, wt1 - wt0);

      }

      if (type_id == 0)
        BFT_FREE(y);
      BFT_FREE(x);
    }

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Double local dot product.
 *
 * parameters:
 *   t_measure <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_dot_product_2(double  t_measure)
{
  double wt0, wt1;
  int sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  double test_sum = 0.0;

  cs_real_t *restrict x = NULL;
  cs_real_t *restrict y = NULL;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y)
# endif

  /* First simple local x.x version */

#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) ||defined(HAVE_MKL)

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = 2 * (n*2 - 1);

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) || defined(HAVE_MKL)
        double s1 = cblas_ddot(n, x, 1, x, 1);
        double s2 = cblas_ddot(n, x, 1, y, 1);
#endif
        test_sum += test_sum_mult*(s1+s2);
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Double local dot product X.X, X.Y with %s BLAS (%d elts.)\n"
               "---------------------------------\n",
               ext_blas_type, (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);
    BFT_FREE(y);
  }

#endif /* external BLAS */

  /* Local dot product */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = 2 * (n*2 - 1);

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);

    for (cs_blas_reduce_t j = 0; j < 2; j++) {

      cs_blas_set_reduce_algorithm(j);

      wt0 = cs_timer_wtime(), wt1 = wt0;

      if (t_measure > 0)
        n_runs = 8;
      else
        n_runs = 1;
      run_id = 0;
      while (run_id < n_runs) {
        double s1, s2;
        double test_sum_mult = 1.0/n_runs;
        while (run_id < n_runs) {
          cs_dot_xx_xy(n, x, y, &s1, &s2);
          test_sum += test_sum_mult*s1;
          run_id++;
        }
        wt1 = cs_timer_wtime();
        if (wt1 - wt0 < t_measure)
          n_runs *= 2;
      }

      bft_printf("\n"
                 "Double local dot product X.X, X.Y with %s (%d elts.)\n"
                 "---------------------------------\n", reduce_name[j], (int)n);

      bft_printf("  (calls: %d)\n", n_runs);

      _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    }

    BFT_FREE(x);
    BFT_FREE(y);
  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Simple dot product with mixed precision (double/float)
 *
 * parameters:
 *   t_measure <-- minimum time for each measure (< 0 for single pass)
 *   global    <-- 0 for local use, 1 for MPI sum
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_dot_product_m(double   t_measure,
               int      global)
{
  double wt0, wt1;
  int sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  double test_sum = 0.0;
  int _global = global;

  double *x = NULL;
  float  *y = NULL;

  if (cs_glob_n_ranks == 1)
    _global = 0;

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n*2 - 1;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, float);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        double s1 = _cs_dot_superblock_m(n, x, y);
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
      bft_printf("\n"
                 "Local dot product X.Y (sb, %d elts, mixed precision)\n"
                 "-----------------\n",
                 (int)n);
    else
      bft_printf("\n"
                 "Global dot product X.Y (sb, %d elts/rank, mixed precision)\n"
                 "------------------\n",
                 (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * y -> ax + y test
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_axpy_test(double  t_measure)
{
  double wt0, wt1;
  int sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  cs_real_t *restrict x = NULL;
  cs_real_t *restrict y = NULL;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y)
# endif

  /* First simple local x.x version */

#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) ||defined(HAVE_MKL)

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n*2;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);
    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
    n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#if defined(HAVE_ATLAS) || defined(HAVE_CBLAS) || defined(HAVE_MKL)
        cblas_daxpy(n, test_sum_mult, x, 1, y, 1);
#endif
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- aX + Y with %s BLAS (%d elts.)\n"
               "-----------------------------\n",
               ext_blas_type, (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);
    BFT_FREE(y);
  }

#endif /* external BLAS */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n*2;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);
    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          y[ii] += test_sum_mult * x[ii];
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- aX + Y (%d elts.)\n"
               "-----------\n",
               (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);
    BFT_FREE(y);
  }

  /* Variant with alpha = -1 */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n*2;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 + 1)*_pi;
      y[ii] = (ii%10 + 1);
    }

    test_sum = test_sum - floor(test_sum);
    wt0 = cs_timer_wtime(), wt1 = wt0;

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      y[ii] = 0.0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          y[ii] -= x[ii];
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      if (run_id % 4096 && run_id > 0) {
        for (ii = 0; ii < n; ii++)
          y[ii] = 0.0;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- -X + Y (%d elts.)\n"
               "-----------\n",
               (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);
    BFT_FREE(y);
  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Simple divisions on a vector.
 *
 * parameters:
 *   t_measure <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_division_test(double  t_measure)
{
  double wt0, wt1;
  int sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  double test_sum = 0.0;

  cs_real_t  *x = NULL, *y = NULL, *z = NULL;

  /* Division of 2 vectors */
  /*-----------------------*/

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);
    BFT_MALLOC(z, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = 2.0 + ii%3;
      y[ii] = 2.0 + (ii+1)%3;
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          z[ii] = x[ii] / y[ii];
        test_sum += z[n-1];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Division z = x/y (%d elts.)\n"
               "----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(z);
    BFT_FREE(y);
    BFT_FREE(x);

  }

  /* Copy inverse of a vector */
  /*--------------------------*/

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = 2.0 + (ii+1)%3;

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          y[ii] = 1.0 / x[ii];
        test_sum += y[n-1];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Division y = 1/x (%d elts.)\n"
               "----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);

  }

  /* Inverse of a vector */
  /*---------------------*/

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = 2.0 + (ii+1)%3;

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          x[ii] = 1.0 / x[ii];
        test_sum += x[n-1];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Division x = 1/x (%d elts.)\n"
               "----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Simple square root on a vector.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_sqrt_test(double  t_measure)
{
  double wt0, wt1;
  int sub_id, run_id, n_runs;
  long n_ops;
  cs_lnum_t n, ii;

  double test_sum = 0.0;

  cs_real_t  *restrict x = NULL, *restrict y = NULL;

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = 2.0 + ii%3;

    test_sum = test_sum - floor(test_sum);

    /* Copy square root of a vector */

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          y[ii] = sqrt(x[ii]);
        test_sum += y[n-1];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "y = sqrt(x) (%d elts.)\n"
               "-----------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
  }

  /* In place square root of a vector */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = 2.0 + ii%3;

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;
    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          x[ii] = sqrt(x[ii]);
        test_sum += x[n-1];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "x = sqrt(x) (%d elts.)\n"
               "-----------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(x);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * y -> da.x test
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_ad_x_test(double  t_measure)
{
  double wt0, wt1;
  int    sub_id, run_id, n_runs;
  long   n_ops;
  int    n, ii;

  double test_sum = 0.0;

  cs_real_t  *restrict da = NULL;
  cs_real_t  *restrict x = NULL;
  cs_real_t  *restrict y = NULL;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *da)
# endif

#if defined(HAVE_CBLAS) || defined(HAVE_MKL)

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n, cs_real_t);
    BFT_MALLOC(x, n, cs_real_t);
    BFT_MALLOC(y, n, cs_real_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      da[ii] = 1.0 + ii%3;
      x[ii] = 2.0 + ii%3;
    }

    test_sum = test_sum - floor(test_sum);

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
                    n,        /* n rows */
                    n,        /* n columns */
                    0,        /* kl */
                    0,        /* ku */
                    1.0,      /* alpha */
                    da,       /* matrix */
                    1,        /* lda */
                    x, 1, 0.0, y, 1);
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- DX with BLAS dgbmv (%d elts.)\n"
               "-----------------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

#endif /* defined(HAVE_CBLAS) || defined(HAVE_MKL) */

#if defined(HAVE_MKL)

  /* Diagonal variant */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n, cs_real_t);
    BFT_MALLOC(x, n, cs_real_t);
    BFT_MALLOC(y, n, cs_real_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      da[ii] = 1.0 + ii%3;
      x[ii] = 2.0 + ii%3;
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      char transa = 'n';
      int n_rows = n;
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
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- DX with MKL ddiagemv (%d elts.)\n"
               "-------------------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

#endif /* defined(HAVE_MKL) */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n, cs_real_t);
    BFT_MALLOC(x, n, cs_real_t);
    BFT_MALLOC(y, n, cs_real_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      da[ii] = 1.0 + ii%3;
      x[ii] = 2.0 + ii%3;
    }

    test_sum = test_sum - floor(test_sum);

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
        for (ii = 0; ii < n; ii++)
          y[ii] = da[ii] * x[ii];
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Y <- DX (%d elts.)\n"
               "-------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

  return test_sum;
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
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_block_ad_x_test(double  t_measure)
{
  double wt0, wt1;
  int    sub_id, run_id, n_runs;
  long   n_ops;
  int    n, ii;
  cs_real_t  *restrict da = NULL, *restrict x = NULL, *restrict y = NULL;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y, *da)
# endif

  /* Variant a */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n*15, cs_real_t);
    BFT_MALLOC(x, n*3, cs_real_t);
    BFT_MALLOC(y, n*3, cs_real_t);

    n_ops = n * (3+2) * 3;

    /* First simple local x.x version */

#   pragma omp parallel for
    for (ii = 0; ii < n*15; ii++)
      da[ii] = 1.0;
#   pragma omp parallel for
    for (ii = 0; ii < n*3; ii++) {
      x[ii] = 0.0;
      y[ii] = 0.0;
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        _3_3_diag_vec_p_l_a(da, x, y, n);
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Block Y <- DX (compiler inline, %d elts.)\n"
               "-------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

  /* Variant b */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n*15, cs_real_t);
    BFT_MALLOC(x, n*3, cs_real_t);
    BFT_MALLOC(y, n*3, cs_real_t);

    n_ops = n * (3+2) * 3;

    /* First simple local x.x version */

#   pragma omp parallel for
    for (ii = 0; ii < n*15; ii++)
      da[ii] = 1.0;
#   pragma omp parallel for
    for (ii = 0; ii < n*3; ii++) {
      x[ii] = 0.0;
      y[ii] = 0.0;
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        _3_3_diag_vec_p_l_b(da, x, y, n);
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Block Y <- DX (manual inline, %d elts.)\n"
               "-------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

  /* Variant c */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];
    n_ops = n;

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(da, n*15, cs_real_t);
    BFT_MALLOC(x, n*3, cs_real_t);
    BFT_MALLOC(y, n*3, cs_real_t);

    n_ops = n * (3+2) * 3;

    /* First simple local x.x version */

#   pragma omp parallel for
    for (ii = 0; ii < n*15; ii++)
      da[ii] = 1.0;
#   pragma omp parallel for
    for (ii = 0; ii < n*3; ii++) {
      x[ii] = 0.0;
      y[ii] = 0.0;
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        _3_3_diag_vec_p_l_c(da, x, y, n);
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Block Y <- DX (3x3 manually unrolled, %d elts.)\n"
               "-------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_stats(n_runs, n_ops, 0, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);
    BFT_FREE(da);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Memory copy test using different variants.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_copy_test(double  t_measure)
{
  double     wt0, wt1;
  int        sub_id, run_id, n_runs;
  cs_lnum_t  n, ii;

  cs_real_t  *restrict x = NULL, *restrict y = NULL;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *y)
# endif

  /* Copy with memcpy */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = (ii%10 + 1);

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        memcpy(y, x, n*sizeof(cs_real_t));
        y[0] += 0.1;
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Copy with memcpy (%d elts.)\n"
               "----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_mem_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);

  }

  /* Copy with loop */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++)
      x[ii] = (ii%10 + 1);

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          y[ii] = x[ii];
        y[0] += 0.1;
        test_sum += test_sum_mult*y[run_id%n];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Copy with loop (%d elts.)\n"
               "--------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_mem_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(y);
    BFT_FREE(x);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------
 * Solve linear system for 3x3 dense blocks: x = a.b
 *
 * This variant uses elimination without pivoting, and copies
 * the matrix to allow further use; the right-hand side may be  modified,
 * and is considered a work array.
 *
 * This assumes the matrix is invertible
 *
 * parameters:
 *   a      <-- block matrix
 *   b      <-> right-hand side
 *   x      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_sv_pc(const cs_real_t   a[restrict 3][3],
                 cs_real_t         b[restrict 3],
                 cs_real_t         x[restrict 3])
{
  int i, j;

  double f[2];
  double _a[3][3];

  /* Copy array */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      _a[i][j] = a[i][j];
    }
  }

  /* Forward elimination (unrolled) */

  f[0] = _a[1][0] / _a[0][0];
  f[1] = _a[2][0] / _a[0][0];

  _a[1][0] = 0.0;
  _a[1][1] -= _a[0][1]*f[0];
  _a[1][2] -= _a[0][2]*f[0];

  _a[2][0] = 0.0;
  _a[2][1] -= _a[0][1]*f[1];
  _a[2][2] -= _a[0][2]*f[1];

  b[1] -= b[0]*f[0];
  b[2] -= b[0]*f[1];

  /* Eliminate values */

  f[0] = _a[2][1] / _a[1][1];

  _a[2][1] = 0.0;
  _a[2][2] -= _a[1][2]*f[0];
  b[2] -= b[1]*f[0];

  /* Solve system */

  x[2] =  b[2]                                   / _a[2][2];
  x[1] = (b[1] - _a[1][2]*x[2])                  / _a[1][1];
  x[0] = (b[0] - _a[0][2]*x[2]  - _a[0][1]*x[1]) / _a[0][0];
}

/*----------------------------------------------------------------------------
* LU factorization of dense 3*3 matrices.
*
* parameters:
*   a  <-> linear equation matrix in, lu decomposition out
*----------------------------------------------------------------------------*/

static inline void
_fact_lu33(cs_real_t  a[restrict 3][3])
{
  cs_real_t lu[3][3];

  lu[1][0] = a[1][0] / a[0][0];
  lu[1][1] = a[1][1] - a[1][0]*a[0][1];
  lu[1][2] = a[1][2] - a[1][0]*a[0][2];

  lu[2][0] = a[2][0]/a[0][0];
  lu[2][1] = (a[2][1] - lu[2][0]*a[0][1])/lu[1][1];
  lu[2][2] = a[2][2] - lu[2][0]*a[0][2] - lu[2][1]*lu[1][2];

  a[1][0] = lu[1][0];
  a[1][1] = lu[1][1];
  a[1][2] = lu[1][2];

  a[2][0] = lu[2][0];
  a[2][1] = lu[2][1];
  a[2][2] = lu[2][2];
}

/*----------------------------------------------------------------------------
 * Compute forward and backward to solve an LU 3*3 system.
 *
 * parameters:
 *   mat   <-- 3*3 matrix
 *   b     <-- RHS
 *   x     --> solution
 *----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_lu33(const cs_real_t  a[restrict 3][3],
                const cs_real_t  b[restrict 3],
                cs_real_t        x[restrict 3])
{
  cs_real_t aux = b[1] - b[0]*a[1][0];

  x[2] = (b[2] - b[0]*a[2][0] - aux*a[2][1])/a[2][2];
  x[1] = (aux - a[1][2]*x[2])/a[1][1];
  x[0] = (b[0] - a[0][1]*x[1] - a[0][2]*x[2])/a[0][0];
}

/*----------------------------------------------------------------------------
 * Solve linear system for 3x3 dense blocks: x = a.b
 *
 * This variant uses LU decomposition.
 *
 * This assumes the matrix is invertible
 *
 * parameters:
 *   a      <-- block matrix
 *   b      <-> right-hand side
 *   x      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_sv_lu(const cs_real_t   a[restrict 3][3],
                 cs_real_t         b[restrict 3],
                 cs_real_t         x[restrict 3])
{
  cs_real_t lu[3][3];

  lu[1][0] = a[1][0] / a[0][0];
  lu[1][1] = a[1][1] - a[1][0]*a[0][1];
  lu[1][2] = a[1][2] - a[1][0]*a[0][2];

  lu[2][0] = a[2][0]/a[0][0];
  lu[2][1] = (a[2][1] - lu[2][0]*a[0][1])/lu[1][1];
  lu[2][2] = a[2][2] - lu[2][0]*a[0][2] - lu[2][1]*lu[1][2];

  cs_real_t aux = b[1] - b[0]*lu[1][0];

  x[2] = (b[2] - b[0]*lu[2][0] - aux*lu[2][1])/lu[2][2];
  x[1] = (aux - lu[1][2]*x[2])/lu[1][1];
  x[0] = (b[0] - a[0][1]*x[1] - a[0][2]*x[2])/a[0][0];
}

/*----------------------------------------------------------------------------
 * Compute matrix-vector product for dense blocks: y[i] = a[i].x[i]
 *
 * This variant uses a fixed 3x3 block, for better compiler optimization,
 * and an inline function
 *
 * parameters:
 *   a      <-- Pointer to block matrixes array (usually matrix diagonal)
 *   b      <-- Multipliying vector values
 *   x      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_axl(const cs_real_t  a[restrict 3][3],
               const cs_real_t  b[restrict 3],
               cs_real_t        x[restrict 3])
{
  x[0] = a[0][0]*b[0] + a[0][1]*b[1] + a[0][2]*b[2];
  x[1] = a[1][0]*b[0] + a[1][1]*b[1] + a[1][2]*b[2];
  x[2] = a[2][0]*b[0] + a[2][1]*b[1] + a[2][2]*b[2];
}

/*----------------------------------------------------------------------------
 * Solve linear system for 3x3 dense blocks: x = a.b
 *
 * This variant uses Cramer's rule.
 *
 * This assumes the matrix is invertible
 *
 * parameters:
 *   a      <-- block matrix
 *   b      <-> right-hand side
 *   x      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_sv_c(const cs_real_t   a[restrict 3][3],
                const cs_real_t   b[restrict 3],
                cs_real_t         x[restrict 3])
{
  cs_real_t t[12], di;

  t[0]  = a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]);
  t[1]  = a[1][0] * (a[2][1]*a[0][2] - a[0][1]*a[2][2]);
  t[2]  = a[2][0] * (a[0][1]*a[1][2] - a[1][1]*a[0][2]);

  t[3]  = b[0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]);
  t[4]  = b[1] * (a[2][1]*a[0][2] - a[0][1]*a[2][2]);
  t[5]  = b[2] * (a[0][1]*a[1][2] - a[1][1]*a[0][2]);

  t[6]  = a[0][0] * (b[1]*a[2][2] - b[2]*a[1][2]);
  t[7]  = a[1][0] * (b[2]*a[0][2] - b[0]*a[2][2]);
  t[8]  = a[2][0] * (b[0]*a[2][1] - b[1]*a[0][2]);

  t[9]  = a[0][0] * (a[1][1]*b[2] - a[2][1]*b[2]);
  t[10] = a[1][0] * (a[2][1]*b[2] - a[0][1]*b[2]);
  t[11] = a[2][0] * (a[0][1]*b[2] - a[1][1]*b[2]);

  di = 1.0 / (t[0] + t[1]  + t[2]);

  x[0] = (t[3] + t[4]  + t[5])*di;
  x[1] = (t[6] + t[7]  + t[8])*di;
  x[2] = (t[9] + t[10] + t[11])*di;
}

/*----------------------------------------------------------------------------
 * Solve linear system for 3x3 dense blocks: x = a.b
 *
 * This variant uses elimination without pivoting, and may modify both
 * the matrix and the right-hand side, which are considered work arrays.
 *
 * This assumes the matrix is invertible
 *
 * parameters:
 *   a      <-- block matrix
 *   b      <-> right-hand side
 *   x      --> Resulting vector
 *----------------------------------------------------------------------------*/

static inline void
_dense_3_3_sv_p(cs_real_t   a[restrict 3][3],
                cs_real_t   b[restrict 3],
                cs_real_t   x[restrict 3])
{
  double factor;

  /* Forward elimination (unrolled) */

  factor = a[1][0] / a[0][0];

  a[1][0] = 0.0;
  a[1][1] -= a[0][1]*factor;
  a[1][2] -= a[0][2]*factor;
  b[1] -= b[0]*factor;

  factor = a[2][0] / a[0][0];

  a[2][0] = 0.0;
  a[2][1] -= a[0][1]*factor;
  a[2][2] -= a[0][2]*factor;
  b[2] -= b[0]*factor;

  /* Eliminate values */

  factor = a[2][1] / a[1][1];

  a[2][1] = 0.0;
  a[2][2] -= a[1][2]*factor;
  b[2] -= b[1]*factor;

  /* Solve system */

  x[2] =  b[2]                                 / a[2][2];
  x[1] = (b[1] - a[1][2]*x[2])                 / a[1][1];
  x[0] = (b[0] - a[0][2]*x[2]  - a[0][1]*x[1]) / a[0][0];
}

/*----------------------------------------------------------------------------
 * Test solving of 3x3 linear systems using different variants.
 *
 * parameters:
 *   t_measure     <-- minimum time for each measure (< 0 for single pass)
 *
 * returns:
 *   test sum (ensures compiler may not optimize loops out)
 *----------------------------------------------------------------------------*/

static double
_solve_33_test(double  t_measure)
{
  double     wt0, wt1, wt2, wt3;
  int        sub_id, run_id, n_runs;
  cs_lnum_t  n, ii;

  cs_real_3_t  *restrict x = NULL, *restrict b = NULL;
  cs_real_33_t  *restrict a = NULL;

  double test_sum = 0.0;

  /* Tell IBM compiler not to alias */
# if defined(__xlc__)
# pragma disjoint(*x, *a, *b)
# endif

  /* Gaussian elimination */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, cs_real_3_t);
    BFT_MALLOC(b, n, cs_real_3_t);
    BFT_MALLOC(a, n, cs_real_33_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++)
          a[ii][jj][kk] = 1.;
      }
      for (int jj = 0; jj < 3; jj++) {
        b[ii][jj] = (ii%10 + 1);
        a[ii][jj][jj] += 2.;
      }
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _dense_3_3_sv_pc(a[ii], b[ii], x[ii]);
        test_sum += test_sum_mult*x[run_id%n][0];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Solve 3x3 systems with Gaussian elimination (%d elts.)\n"
               "-----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_time_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(b);
    BFT_FREE(x);
    BFT_FREE(a);

  }

  /* LU */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, cs_real_3_t);
    BFT_MALLOC(b, n, cs_real_3_t);
    BFT_MALLOC(a, n, cs_real_33_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++)
          a[ii][jj][kk] = 1.;
      }
      for (int jj = 0; jj < 3; jj++) {
        b[ii][jj] = (ii%10 + 1);
        a[ii][jj][jj] += 2.;
      }
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _dense_3_3_sv_lu(a[ii], b[ii], x[ii]);
        test_sum += test_sum_mult*x[run_id%n][0];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Solve 3x3 systems with LU (%d elts.)\n"
               "-----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_time_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(b);
    BFT_FREE(x);
    BFT_FREE(a);

  }

  /* Cramer */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, cs_real_3_t);
    BFT_MALLOC(b, n, cs_real_3_t);
    BFT_MALLOC(a, n, cs_real_33_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++)
          a[ii][jj][kk] = 1.;
      }
      for (int jj = 0; jj < 3; jj++) {
        b[ii][jj] = (ii%10 + 1);
        a[ii][jj][jj] += 2.;
      }
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _dense_3_3_sv_c(a[ii], b[ii], x[ii]);
        test_sum += test_sum_mult*x[run_id%n][0];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Solve 3x3 systems with Cramer's rule (%d elts.)\n"
               "-----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_time_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(b);
    BFT_FREE(x);
    BFT_FREE(a);

  }

  /* 2-stage LU */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, cs_real_3_t);
    BFT_MALLOC(b, n, cs_real_3_t);
    BFT_MALLOC(a, n, cs_real_33_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++)
          a[ii][jj][kk] = 1.;
      }
      for (int jj = 0; jj < 3; jj++) {
        b[ii][jj] = (ii%10 + 1);
        a[ii][jj][jj] += 2.;
      }
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        wt2 = cs_timer_wtime();
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++)
              a[ii][jj][kk] = 1.;
          }
          for (int jj = 0; jj < 3; jj++) {
            a[ii][jj][jj] += 2.;
          }
        }
        wt3 = cs_timer_wtime();
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _fact_lu33(a[ii]);
        test_sum += test_sum_mult*a[run_id%n][0][0];
        run_id++;
      }
      wt1 = cs_timer_wtime() - (wt3-wt2);
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Prepare 3x3 systems LU decomposition (%d elts.)\n"
               "-------------------\n", (int)n);
    bft_printf("  (calls: %d)\n", n_runs);
    _print_time_stats(n_runs, n, wt1 - wt0);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _fw_and_bw_lu33(a[ii], b[ii], x[ii]);
        test_sum += test_sum_mult*x[run_id%n][0];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Apply LU decomposition to 3x3 systems (%d elts.)\n"
               "-----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_time_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(b);
    BFT_FREE(x);
    BFT_FREE(a);

  }

  /* 2-stage Cramer */

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    n = _n_elts[sub_id];

    /* Realloc and initialize arrays for each test, as
       first touch may affect memory locality on some systems */

    BFT_MALLOC(x, n, cs_real_3_t);
    BFT_MALLOC(b, n, cs_real_3_t);
    BFT_MALLOC(a, n, cs_real_33_t);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < 3; kk++)
          a[ii][jj][kk] = 1.;
      }
      for (int jj = 0; jj < 3; jj++) {
        b[ii][jj] = (ii%10 + 1);
        a[ii][jj][jj] += 2.;
      }
    }

    test_sum = test_sum - floor(test_sum);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
        wt2 = cs_timer_wtime();
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            for (int kk = 0; kk < 3; kk++)
              a[ii][jj][kk] = 1.;
          }
          for (int jj = 0; jj < 3; jj++) {
            a[ii][jj][jj] += 2.;
          }
        }
        wt3 = cs_timer_wtime();
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          cs_math_33_inv_cramer_in_place(a[ii]);
        test_sum += test_sum_mult*a[run_id%n][0][0];
        run_id++;
      }
      wt1 = cs_timer_wtime() - (wt3-wt2);
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Prepare 3x3 systems Cramer (%d elts.)\n"
               "-------------------\n", (int)n);
    bft_printf("  (calls: %d)\n", n_runs);
    _print_time_stats(n_runs, n, wt1 - wt0);

    wt0 = cs_timer_wtime(), wt1 = wt0;

    if (t_measure > 0)
      n_runs = 8;
    else
      n_runs = 1;
    run_id = 0;
    while (run_id < n_runs) {
      double test_sum_mult = 1.0/n_runs;
      while (run_id < n_runs) {
#       pragma omp parallel for
        for (ii = 0; ii < n; ii++)
          _dense_3_3_axl(a[ii], b[ii], x[ii]);
        test_sum += test_sum_mult*x[run_id%n][0];
        run_id++;
      }
      wt1 = cs_timer_wtime();
      if (wt1 - wt0 < t_measure)
        n_runs *= 2;
    }

    bft_printf("\n"
               "Apply inverse matrix to 3x3 systems (%d elts.)\n"
               "-----------------\n", (int)n);

    bft_printf("  (calls: %d)\n", n_runs);

    _print_time_stats(n_runs, n, wt1 - wt0);

    BFT_FREE(b);
    BFT_FREE(x);
    BFT_FREE(a);

  }

  return test_sum;
}

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

  int sub_id;
  double t_measure = 1.0;
  double test_sum = 0.0;

  /* Internationalization */

#ifdef HAVE_SETLOCALE
  if (!setlocale (LC_ALL,"")) {
#if defined (DEBUG)
     bft_printf("locale not supported by C library"
                " or bad LANG environment variable");
#endif
  }
#endif /* HAVE_SETLOCALE */

  /* Initialization and environment */

#if defined(HAVE_MPI)
  _mpi_init();
#endif

  _pi = 4 * atan(1);

  (void)cs_timer_wtime();

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif

  bft_printf("\n");

  /* Precision tests */
  /*-----------------*/

  for (sub_id = 0; sub_id < _n_sizes; sub_id++) {

    cs_lnum_t ii;
    double s;

    cs_lnum_t n = _n_elts[sub_id];
    double ref_s = 165 * (n/10) * _pi;
    double *x = NULL, *y = NULL;

    /* Initialize arrays */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (ii = 0; ii < n; ii++) {
      x[ii] = (ii%10 - 3)*_pi;
      y[ii] = (ii%10 + 1);
    }

    s = _ddot_canonical(n, x, y);
    bft_printf("  %-22s dot product error for n = %7d: %12.5e\n",
               "canonical", (int)n, ref_s - s);

    for (cs_blas_reduce_t j = 0; j < 2; j++) {
      cs_blas_set_reduce_algorithm(j);
      s = cs_dot(n, x, y);
      bft_printf("  %-22s dot product error for n = %7d: %12.5e\n",
                 reduce_name[j], (int)n, ref_s - s);
    }

#if defined(HAVE_MKL)
    s = cblas_ddot(n, x, 1, y, 1);
    bft_printf("  %-22s dot product error for n = %7d: %12.5e\n",
               "MKL", (int)n, ref_s - s);
#elif defined(HAVE_ATLAS)
    s = cblas_ddot(n, x, 1, y, 1);
    bft_printf("  %-22s dot product error for n = %7d: %12.5e\n",
               "ATLAS", (int)n, ref_s - s);
#elif defined(HAVE_CBLAS)
    s = cblas_ddot(n, x, 1, y, 1);
    bft_printf("  %-22s dot product error for n = %7d: %12.5e\n",
               "C BLAS", (int)n, ref_s - s);
#endif

    bft_printf("\n");

    BFT_FREE(y);
    BFT_FREE(x);
  }

  /* Performance tests */
  /*-------------------*/

  if (cs_glob_n_ranks > 1)
    test_sum += _dot_product_1(t_measure, 1);

  test_sum += _dot_product_1(t_measure, 0);

  test_sum += _dot_product_2(t_measure);

  test_sum += _dot_product_m(t_measure, 0);

  test_sum += _axpy_test(t_measure);

  test_sum += _division_test(t_measure);

  test_sum += _sqrt_test(t_measure);

  test_sum += _ad_x_test(t_measure);

  test_sum += _block_ad_x_test(t_measure);

  test_sum += _copy_test(t_measure);

  test_sum += _solve_33_test(t_measure);

  if (isnan(test_sum))
    bft_printf("test_sum is NaN\n");

  /* Finalize */

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Finalize();
  }
#endif /* HAVE_MPI */

  exit (EXIT_SUCCESS);
}
