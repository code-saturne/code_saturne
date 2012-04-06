/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_mem_usage.h"

#include "cs_system_info.h"

#include "cs_blas.h"
#include "cs_defs.h"
#include "cs_timer.h"

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * False print of a message to standard output for discarded logs
 *----------------------------------------------------------------------------*/

static int
_bft_printf_null(const char  *format,
                 va_list      arg_ptr)
{
  return 0;
}

/*----------------------------------------------------------------------------
 * Analysis of environment variables to determine
 * if we require MPI, and initialization if necessary.
 *----------------------------------------------------------------------------*/

static void
_mpi_init(void)
{

  char *s;

  int arg_id = 0, flag = 0;
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
_ddot_simple(cs_lnum_t      n,
             const double  *x,
             const double  *y)
{
  cs_lnum_t  i;
  double     s = 0;

  if (n < 1)
    return s;

# pragma omp parallel for reduction(+: s)
  for (i = 0; i < n; i++) {
    s += (x[i] * y[i]);
  }

  return s;
}

/*----------------------------------------------------------------------------
 * Return the dot product of 2 vectors: x.y
 *
 * The algorithm used is l3superblock60, based on the article:
 * "Reducing Floating Point Error in Dot Product Using the Superblock Family
 * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
 * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156–1174
 * 2008 Society for Industrial and Applied Mathematics
 *
 * parameters:
 *   n <-- size of arrays x and y
 *   x <-- array of floating-point values
 *   y<-- array of floating-point values
 *
 * returns:
 *   dot product
 *----------------------------------------------------------------------------*/

static double
_ddot_l3superblock60(cs_lnum_t      n,
                     const double  *x,
                     const double  *y)
{
  const cs_lnum_t block_size = 60;

  cs_lnum_t sid, bid, i;
  cs_lnum_t start_id, end_id;
  double sdot, cdot;

  cs_lnum_t n_blocks = n / 60;
  cs_lnum_t n_sblocks = sqrt(n_blocks);
  cs_lnum_t blocks_in_sblocks = (n_sblocks > 0) ? n_blocks / n_sblocks : 0;

  double dot = 0.0;

# pragma omp parallel for reduction(+:dot) private(bid, start_id, end_id, i, \
                                                   cdot, sdot) if (n > THR_MIN)
  for (sid = 0; sid < n_sblocks; sid++) {

    sdot = 0.0;

    for (bid = 0; bid < blocks_in_sblocks; bid++) {
      start_id = block_size * bid*sid;
      end_id = block_size * (bid*sid + 1);
      cdot = 0.0;
      for (i = start_id; i < end_id; i++)
        cdot += x[i]*y[i];
      sdot += cdot;
    }

    dot += sdot;

  }

  cdot = 0.0;
  start_id = block_size * n_sblocks*blocks_in_sblocks;
  end_id = n;
  for (i = start_id; i < end_id; i++)
    cdot += x[i] * y[i];
  dot += cdot;

  return dot;
}

/*----------------------------------------------------------------------------*/

int
main (int argc, char *argv[])
{
  cs_lnum_t n;
  double walltime, cputime;

  /* Internationalization */

#ifdef HAVE_SETLOCALE
  if (!setlocale (LC_ALL,"")) {
#if defined (DEBUG)
     printf("locale not supported by C library"
            " or bad LANG environment variable");
#endif
  }
#endif /* HAVE_SETLOCALE */

  /* Initialization and environment */

  (void)cs_timer_wtime();

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif

  printf("\n");

  walltime = cs_timer_wtime();
  cputime  = cs_timer_cpu_time();

  /* Loop on array sizes */

  for (n = 10; n < 10000000; n *= 10) {

    cs_lnum_t i;
    double s;
    double pi = 3.14159265358979323846;
    double ref_s = 385 * (n/10) * pi;
    double *x = NULL, *y = NULL;

    /* Initialize arrays */

    BFT_MALLOC(x, n, double);
    BFT_MALLOC(y, n, double);

#   pragma omp parallel for
    for (i = 0; i < n; i++) {
      x[i] = (i%10 + 1)*pi;
      y[i] = (i%10 + 1);
    }

    s = _ddot_simple(n, x, y);
    printf("  Simple     dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);

    s = _ddot_l3superblock60(n, x, y);
    printf("  Superblock dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);

#if defined(HAVE_ESSL)
    s = ddot(n, (double *)x, 1, (double *)y, 1);
    printf("  ESSL       dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);
#elif defined(HAVE_ACML)
    s = ddot(n, (double *)x, 1, (double *)y, 1);
    printf("  ACML       dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);
#elif defined(HAVE_MKL)
    s = cblas_ddot(n, x, 1, y, 1);
    printf("  MKL        dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);
#elif defined(HAVE_ATLAS)
    s = cblas_ddot(n, x, 1, y, 1);
    printf("  ATLAS      dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);
#elif defined(HAVE_CBLAS)
    s = cblas_ddot(n, x, 1, y, 1);
    printf("  BLAS       dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);
#endif

    s = cs_dot(n, x, y);
    printf("  Default    dot product error for n = %7d: %12.5e\n",
           (int)n, ref_s - s);

    printf("\n");

    BFT_FREE(x);
    BFT_FREE(y);
  }

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
