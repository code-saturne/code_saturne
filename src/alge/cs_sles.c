/*============================================================================
 * Sparse Linear Equation Solvers
 *============================================================================*/

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_log.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_matrix.h"
#include "cs_matrix_util.h"
#include "cs_post.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define EPZERO  1.E-12
#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures
 * through at least 2012) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_sles_info_t {

  char                *name;               /* System name */
  cs_sles_type_t       type;               /* Solver type */

  unsigned             n_calls;            /* Number of times system solved */

  unsigned             n_iterations_last;  /* Number of iterations for last
                                              system resolution */
  unsigned             n_iterations_min;   /* Minimum number ot iterations
                                              in system resolution history */
  unsigned             n_iterations_max;   /* Maximum number ot iterations
                                              in system resolution history */
  unsigned long long   n_iterations_tot;   /* Total accumulated number of
                                              iterations */

  cs_timer_counter_t   t_tot;              /* Total time used */

} cs_sles_info_t;

/* Convergence testing and tracking */
/*----------------------------------*/

typedef struct _cs_sles_convergence_t {

  int                  verbosity;          /* Verbosity level */

  unsigned             n_iterations;       /* Current number of iterations */
  unsigned             n_iterations_max;   /* Maximum number of iterations */

  double               precision;          /* Precision limit */
  double               r_norm;             /* Residue normalization */
  double               residue;            /* Current residue */
  double               initial_residue;    /* Initial residue */

} cs_sles_convergence_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

static int cs_glob_sles_n_systems = 0;      /* Current number of systems */
static int cs_glob_sles_n_max_systems = 0;  /* Max. number of sytems for
                                               cs_glob_sles_systems. */

static cs_sles_info_t **cs_glob_sles_systems = NULL; /* System info array */

/*
  Matrix structures re-used for various resolutions.

  These structures are kept throughout the whole run, to avoid paying the
  CPU overhead for their construction at each system resolution
  (at the cost of extra memory use, depending on the chosen structure).

  Two simultaneous matrixes may be needed for some solvers: one for the
  linear system, one for the preconditionner.
  We always have at least one structure of "native" matrix type, as
  this type incurs negligible memory and assignment cpu overhead;
  we use it for Jacobi (where the number of iterations done is often
  small, and an assignment cost equivalent to a few matrix.vector
  products may not be amortized).
*/


cs_matrix_structure_t *cs_glob_sles_native_matrix_struct = NULL;
cs_matrix_t *cs_glob_sles_native_matrix = NULL;

/* Sparse linear equation solver type names */

/* Short names for solver types */

const char *cs_sles_type_name[] = {N_("Conjugate gradient"),
                                   N_("Conjugate gradient, single reduction"),
                                   N_("Jacobi"),
                                   N_("Bi-CGstab"),
                                   N_("GMRES")};

/* Communicator used for reduction operations */

#if defined(HAVE_MPI)
MPI_Comm _cs_sles_mpi_reduce_comm = MPI_COMM_NULL;
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to new linear system info structure.
 *
 * parameters:
 *   name <-- system name
 *   type <-- resolution method
 *
 * returns:
 *   pointer to newly created linear system info structure
 *----------------------------------------------------------------------------*/

static cs_sles_info_t *
_sles_info_create(const char      *name,
                  cs_sles_type_t   type)
{
  cs_sles_info_t *new_info = NULL;

  BFT_MALLOC(new_info, 1, cs_sles_info_t);
  BFT_MALLOC(new_info->name, strlen(name) + 1, char);

  strcpy(new_info->name, name);
  new_info->type = type;

  new_info->n_calls = 0;
  new_info->n_iterations_min = 0;
  new_info->n_iterations_max = 0;
  new_info->n_iterations_last = 0;
  new_info->n_iterations_tot = 0;

  CS_TIMER_COUNTER_INIT(new_info->t_tot);

  return new_info;
}

/*----------------------------------------------------------------------------
 * Destroy linear system info structure.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure pointer
 *----------------------------------------------------------------------------*/

static void
_sles_info_destroy(cs_sles_info_t  **this_info)
{
  if (*this_info != NULL) {
    BFT_FREE((*this_info)->name);
    BFT_FREE(*this_info);
  }
}

/*----------------------------------------------------------------------------
 * Output information regarding linear system resolution.
 *
 * parameters:
 *   this_info <-> pointer to linear system info structure
 *----------------------------------------------------------------------------*/

static void
_sles_info_dump(cs_sles_info_t *this_info)
{
  int n_calls = this_info->n_calls;
  int n_it_min = this_info->n_iterations_min;
  int n_it_max = this_info->n_iterations_max;
  int n_it_mean = (int)(  this_info->n_iterations_tot
                        / ((unsigned long long)n_calls));

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Summary of resolutions for %s (%s):\n"
                  "\n"
                  "  Number of calls:               %12d\n"
                  "  Minimum number of iterations:  %12d\n"
                  "  Maximum number of iterations:  %12d\n"
                  "  Mean number of iterations:     %12d\n"
                  "  Total elapsed time:            %12.3f\n"),
                this_info->name, cs_sles_type_name[this_info->type],
                n_calls, n_it_min, n_it_max, n_it_mean,
                this_info->t_tot.wall_nsec*1e-9);
}

/*----------------------------------------------------------------------------
 * Return pointer to linear system info.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name <-- system name
 *   type <-- resolution method
 *----------------------------------------------------------------------------*/

static cs_sles_info_t *
_find_or_add_system(const char      *name,
                    cs_sles_type_t   type)
{
  int ii, start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find system */

  start_id = 0;
  end_id = cs_glob_sles_n_systems - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp((cs_glob_sles_systems[mid_id])->name, name);
    if (cmp_ret == 0)
      cmp_ret = (cs_glob_sles_systems[mid_id])->type - type;
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If found, return */

  if (cmp_ret == 0)
    return cs_glob_sles_systems[mid_id];

  /* Reallocate global array if necessary */

  if (cs_glob_sles_n_systems >= cs_glob_sles_n_max_systems) {

    if (cs_glob_sles_n_max_systems == 0)
      cs_glob_sles_n_max_systems = 10;
    else
      cs_glob_sles_n_max_systems *= 2;
    BFT_REALLOC(cs_glob_sles_systems,
                cs_glob_sles_n_max_systems,
                cs_sles_info_t*);

  }

  /* Insert in sorted list */

  for (ii = cs_glob_sles_n_systems; ii > mid_id; ii--)
    cs_glob_sles_systems[ii] = cs_glob_sles_systems[ii - 1];

  cs_glob_sles_systems[mid_id] = _sles_info_create(name,
                                                   type);
  cs_glob_sles_n_systems += 1;

  return cs_glob_sles_systems[mid_id];
}

/*----------------------------------------------------------------------------
 * Initialize or reset convergence info structure.
 *
 * parameters:
 *   solver_name <-- solver name
 *   var_name    <-- variable name
 *   convergence <-> Convergence info structure
 *   verbosity   <-- Verbosity level
 *   n_iter_max  <-- Maximum number of iterations
 *   precision   <-- Precision limit
 *   r_norm      <-- Residue normalization
 *   residue     <-- Initial residue
 *----------------------------------------------------------------------------*/

static void
_convergence_init(cs_sles_convergence_t  *convergence,
                  const char             *solver_name,
                  const char             *var_name,
                  int                     verbosity,
                  unsigned                n_iter_max,
                  double                  precision,
                  double                  r_norm,
                  double                  residue)
{
  convergence->verbosity = verbosity;

  convergence->n_iterations = 0;
  convergence->n_iterations_max = n_iter_max;

  convergence->precision = precision;
  convergence->r_norm = r_norm;
  convergence->residue = residue;
  convergence->initial_residue = residue;

  if (verbosity > 1) {
    bft_printf("%s [%s]:\n", solver_name, var_name);
    if (verbosity > 2)
      bft_printf(_("  n_iter     res_abs     res_nor\n"));
  }

}

/*----------------------------------------------------------------------------
 * Convergence test.
 *
 * parameters:
 *   solver_name <-- solver name
 *   var_name    <-- variable name
 *   n_iter      <-- Number of iterations done
 *   residue     <-- Non normalized residue
 *   convergence <-> Convergence information structure
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_convergence_test(const char             *solver_name,
                  const char             *var_name,
                  unsigned                n_iter,
                  double                  residue,
                  cs_sles_convergence_t  *convergence)
{
  const int verbosity = convergence->verbosity;

  const char final_fmt[]
    = N_("  n_iter : %5d, res_abs : %11.4e, res_nor : %11.4e\n");

  /* Update conversion info structure */

  convergence->n_iterations = n_iter;
  convergence->residue = residue;

  /* Print convergence values if high verbosity */

  if (verbosity > 2)
    bft_printf("   %5u %11.4e %11.4e\n",
               n_iter, residue, residue/convergence->r_norm);

  /* If not converged */

  if (residue > convergence->precision * convergence->r_norm) {

    if (n_iter < convergence->n_iterations_max) {
      int diverges = 0;
      if (residue > convergence->initial_residue * 10000.0 && residue > 100.)
        diverges = 1;
#if (__STDC_VERSION__ >= 199901L)
      else if (isnan(residue) || isinf(residue))
        diverges = 1;
#endif
      if (diverges) {
        bft_printf(_("\n\n"
                     "%s [%s]: divergence after %u iterations:\n"
                     "  initial residual: %11.4e; current residual: %11.4e\n"),
                   solver_name, var_name,
                   convergence->n_iterations,
                   convergence->initial_residue, convergence->residue);
        return -2;
      }
      else
        return 0;
    }
    else {
      if (verbosity > 0) {
        if (verbosity == 1) /* Already output if verbosity > 1 */
          bft_printf("%s [%s]:\n", solver_name, var_name);
        if (verbosity <= 2) /* Already output if verbosity > 2 */
          bft_printf(_(final_fmt),
                     n_iter, residue, residue/convergence->r_norm);
        bft_printf(_(" @@ Warning: non convergence\n"));
      }
      return -1;
    }

  }

  /* If converged */

  else {
    if (verbosity == 2) /* Already output if verbosity > 2 */
      bft_printf(final_fmt, n_iter, residue, residue/convergence->r_norm);
    return 1;
  }

}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- first vector in s = x.y
 *   y      <-- second vector in s = x.y
 *
 * returns:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(cs_int_t          n_elts,
             const cs_real_t  *x,
             const cs_real_t  *y)
{
  double s = cs_dot(n_elts, x, y);

#if defined(HAVE_MPI)

  if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                  _cs_sles_mpi_reduce_comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- vector in s1 = x.x and s2 = x.y
 *   y      <-- vector in s2 = x.y
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy(cs_int_t          n_elts,
                    const cs_real_t  *x,
                    const cs_real_t  *y,
                    double           *s1,
                    double           *s2)
{
  double s[2];

  cs_dot_xx_xy(n_elts, x, y, s, s+1);

#if defined(HAVE_MPI)

  if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM,
                  _cs_sles_mpi_reduce_comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- vector in s1 = x.y
 *   y      <-- vector in s1 = x.y and s2 = y.z
 *   z      <-- vector in s2 = y.z
 *   s1     --> result of s1 = x.y
 *   s2     --> result of s2 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xy_yz(cs_int_t          n_elts,
                    const cs_real_t  *x,
                    const cs_real_t  *y,
                    const cs_real_t  *z,
                    double           *s1,
                    double           *s2)
{
  double s[2];

  cs_dot_xy_yz(n_elts, x, y, z, s, s+1);

#if defined(HAVE_MPI)

  if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM,
                  _cs_sles_mpi_reduce_comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute 3 dot products, summing result over all ranks.
 *
 * parameters:
 *   n_elts <-- Local number of elements
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *   s3     --> result of s3 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy_yz(cs_lnum_t         n_elts,
                       const cs_real_t  *x,
                       const cs_real_t  *y,
                       const cs_real_t  *z,
                       double           *s1,
                       double           *s2,
                       double           *s3)
{
  double s[3];

  cs_dot_xx_xy_yz(n_elts, x, y, z, s, s+1, s+2);

#if defined(HAVE_MPI)

  if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
    double _sum[3];

    MPI_Allreduce(s, _sum, 3, MPI_DOUBLE, MPI_SUM, _cs_sles_mpi_reduce_comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
    s[2] = _sum[2];
  }

#endif /* defined(HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
  *s3 = s[2];
}

/*----------------------------------------------------------------------------
 * Residue preconditioning Gk = C.Rk
 *
 * To increase the polynomial order with no major overhead, we may
 * "implicit" the resolution using a red-black mesh coloring.
 *
 * poly_degree:
 *   0: Gk = (1/ad).Rk
 *   1: Gk = (1/ad - (1/ad).ax.(1/ad)).Rk
 *   2: Gk = (1/ad - (1/ad).ax.(1/ad) + (1/ad).ax.(1/ad).ax.(1/ad)).Rk
 *
 * parameters:
 *   n_cells       <--  Local number of cells
 *   poly_degree   <--  preconditioning polynomial degree (0: diagonal)
 *   rotation_mode <--  Halo update option for rotational periodicity
 *   ad_inv        <--  Inverse of matrix diagonal
 *   a             <--  Linear equation matrix
 *   rk            <--  Residue vector
 *   gk            -->  Result vector
 *   wk            ---  Working array
 *----------------------------------------------------------------------------*/

static void
_polynomial_preconditionning(cs_int_t             n_cells,
                             int                  poly_degree,
                             cs_halo_rotation_t   rotation_mode,
                             const cs_real_t     *ad_inv,
                             const cs_matrix_t   *a,
                             const cs_real_t     *rk,
                             cs_real_t           *restrict gk,
                             cs_real_t           *restrict wk)
{
  int deg_id;
  cs_int_t ii;

#if defined(__xlc__)
#pragma disjoint(*ad_inv, *rk, *gk, *wk)
#endif

  /* Polynomial of degree 0 (diagonal)
   *-----------------------------------*/

# pragma omp parallel for if(n_cells > THR_MIN)
  for (ii = 0; ii < n_cells; ii++)
    gk[ii] = rk[ii] * ad_inv[ii];

  /* Polynomial of degree n
   *-----------------------
   *
   *                  n=1                    n=2
   * gk = ((1/ad) - (1/ad).ax.(1/ad) + (1/ad).ax.(1/ad).ax.(1/ad) + ... ).rk
   */

  for (deg_id = 1; deg_id <= poly_degree; deg_id++) {

    /* Compute Wk = (A-diag).Gk */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, gk, wk);

#   pragma omp parallel for if(n_cells > THR_MIN)
    for (ii = 0; ii < n_cells; ii++)
      gk[ii] = (rk[ii] - wk[ii]) * ad_inv[ii];

  }

}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned conjugate gradient.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      <-- Variable name
 *   a             <-- Matrix
 *   poly_degree   <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   convergence   <-- Convergence information structure
 *   rhs           <-- Right hand side
 *   vx            --> System solution
 *   aux_size      <-- Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_conjugate_gradient(const char             *var_name,
                    const cs_matrix_t      *a,
                    int                     diag_block_size,
                    int                     poly_degree,
                    cs_halo_rotation_t      rotation_mode,
                    cs_sles_convergence_t  *convergence,
                    const cs_real_t        *rhs,
                    cs_real_t              *restrict vx,
                    size_t                  aux_size,
                    void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_rows, ii;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, beta, residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk;
  cs_real_t *restrict zk, *restrict wk, *restrict ad_inv;

  unsigned n_iter = 1;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *rk, *dk, *gk, *zk, *wk, *ad_inv)
#endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_PCG]);

  n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
  n_rows = cs_matrix_get_n_rows(a) * diag_block_size;

  /* Allocate work arrays */

  {
    size_t  n_wa = 5;
    size_t  wa_size = CS_SIMD_SIZE(n_cols);

    if (poly_degree > 0)
      n_wa += 1;

    if (aux_vectors == NULL || aux_size < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;

    rk = _aux_vectors + wa_size;
    dk = _aux_vectors + wa_size*2;
    gk = _aux_vectors + wa_size*3;
    zk = _aux_vectors + wa_size*4;

    if (poly_degree > 0)
      wk = _aux_vectors + wa_size*5;
    else
      wk = NULL;
  }

  cs_matrix_copy_diagonal(a, ad_inv);

  /* Initialize arrays */

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] = rk[ii] - rhs[ii];

  /* Polynomial preconditionning of order poly_degre */

  _polynomial_preconditionning(n_rows,
                               poly_degree,
                               rotation_mode,
                               ad_inv,
                               a,
                               rk,
                               gk,
                               wk);

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  rk_gkm1 = _dot_product(n_rows, rk, gk);

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

  /* Descent parameter */

  _dot_products_xy_yz(n_rows, rk, dk, zk, &ro_0, &ro_1);

  alpha =  - ro_0 / ro_1;

# pragma omp parallel if(n_rows > THR_MIN)
  {
#   pragma omp for nowait
    for (ii = 0; ii < n_rows; ii++)
      vx[ii] += (alpha * dk[ii]);

#   pragma omp for nowait
    for (ii = 0; ii < n_rows; ii++)
      rk[ii] += (alpha * zk[ii]);
  }

  /* Convergence test */

  residue = sqrt(_dot_product(n_rows, rk, rk));

  cvg = _convergence_test(sles_name, var_name,
                          n_iter, residue, convergence);

  /* Current Iteration */
  /*-------------------*/

  while (cvg == 0) {

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 a,
                                 rk,
                                 gk,
                                 wk);

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy(n_rows, rk, gk, &residue, &rk_gk);

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(sles_name, var_name,
                              n_iter, residue, convergence);

    if (cvg != 0)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = gk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    _dot_products_xy_yz(n_rows, rk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        rk[ii] += (alpha * zk[ii]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using preconditioned conjugate gradient
 * with single reduction.
 *
 * For more information, see Lapack Working note 56, at
 * http://www.netlib.org/lapack/lawnspdf/lawn56.pdf)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name        <-- Variable name
 *   a               <-- Matrix
 *   diag_block_size <-- Block size of element ii, ii
 *   poly_degree     <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode   <-- Halo update option for rotational periodicity
 *   convergence     <-- Convergence information structure
 *   rhs             <-- Right hand side
 *   vx              --> System solution
 *   aux_size        <-- Number of elements in aux_vectors
 *   aux_vectors     --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_conjugate_gradient_sr(const char             *var_name,
                       const cs_matrix_t      *a,
                       int                     diag_block_size,
                       int                     poly_degree,
                       cs_halo_rotation_t      rotation_mode,
                       cs_sles_convergence_t  *convergence,
                       const cs_real_t        *rhs,
                       cs_real_t              *restrict vx,
                       size_t                  aux_size,
                       void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_rows, ii;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, gk_sk, beta, residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk, *restrict sk;
  cs_real_t *restrict zk, *restrict wk, *restrict ad_inv;

  unsigned n_iter = 1;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *rk, *dk, *gk, *zk, *wk, *ad_inv)
#endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_PCG_SR]);

  n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
  n_rows = cs_matrix_get_n_rows(a) * diag_block_size;

  /* Allocate work arrays */

  {
    size_t  n_wa = 6;
    size_t  wa_size = CS_SIMD_SIZE(n_cols);

    if (poly_degree > 0)
      n_wa += 1;

    if (aux_vectors == NULL || aux_size < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;

    rk = _aux_vectors + wa_size;
    dk = _aux_vectors + wa_size*2;
    gk = _aux_vectors + wa_size*3;
    zk = _aux_vectors + wa_size*4;
    sk = _aux_vectors + wa_size*5;

    if (poly_degree > 0)
      wk = _aux_vectors + wa_size*6;
    else
      wk = NULL;
  }

  cs_matrix_copy_diagonal(a, ad_inv);

  /* Initialize arrays */

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] = rk[ii] - rhs[ii];

  /* Polynomial preconditionning of order poly_degre */
  /* gk = c_1 * rk  (zk = c_1 * rk) */

  _polynomial_preconditionning(n_rows,
                               poly_degree,
                               rotation_mode,
                               ad_inv,
                               a,
                               rk,
                               gk,
                               wk);

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk); /* zk = A.dk */

  /* Descent parameter */

  _dot_products_xy_yz(n_rows, rk, dk, zk, &ro_0, &ro_1);

  alpha =  - ro_0 / ro_1;

  rk_gkm1 = ro_0;

# pragma omp parallel if(n_rows > THR_MIN)
  {
#   pragma omp for nowait
    for (ii = 0; ii < n_rows; ii++)
      vx[ii] += (alpha * dk[ii]);

#   pragma omp for nowait
    for (ii = 0; ii < n_rows; ii++)
      rk[ii] += (alpha * zk[ii]);
  }

  /* Convergence test */

  residue = sqrt(_dot_product(n_rows, rk, rk));

  cvg = _convergence_test(sles_name, var_name,
                          n_iter, residue, convergence);

  /* Current Iteration */
  /*-------------------*/

  while (cvg == 0) {

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 a,
                                 rk,
                                 gk,
                                 wk);

    cs_matrix_vector_multiply(rotation_mode, a, gk, sk);  /* sk = A.zk */

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy_yz(n_rows, rk, gk, sk, &residue, &rk_gk, &gk_sk);

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(sles_name, var_name,
                              n_iter, residue, convergence);

    if (cvg != 0)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

    ro_1 = gk_sk - beta*beta*ro_1;
    ro_0 = rk_gk;

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        dk[ii] = gk[ii] + (beta * dk[ii]);
        vx[ii] = vx[ii] + (alpha * dk[ii]);
      }
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        zk[ii] = sk[ii] + (beta * zk[ii]);
        rk[ii] = rk[ii] + (alpha * zk[ii]);
      }
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      <-- Variable name
 *   a             <-- Linear equation matrix
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   convergence   <-- Convergence information structure
 *   rhs           <-- Right hand side
 *   vx            --> System solution
 *   aux_size      <-- Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_jacobi(const char             *var_name,
        const cs_matrix_t      *a,
        cs_halo_rotation_t      rotation_mode,
        cs_sles_convergence_t  *convergence,
        const cs_real_t        *rhs,
        cs_real_t              *restrict vx,
        size_t                  aux_size,
        void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_rows, ii;
  double  res2, residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict ad_inv, *restrict rk;
  const cs_real_t  *restrict ad;

  unsigned n_iter = 0;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *ad, *ad_inv)
#endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_JACOBI]);

  n_cols = cs_matrix_get_n_columns(a);
  n_rows = cs_matrix_get_n_rows(a);

  ad = cs_matrix_get_diagonal(a);

  /* Allocate work arrays */

  {
    size_t  wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size < (wa_size * 2))
      BFT_MALLOC(_aux_vectors, wa_size * 2, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;
    rk     = _aux_vectors + wa_size;
  }

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad[ii];

  cvg = 0;

  /* Current iteration */
  /*-------------------*/

  while (cvg == 0) {

    register double r;

    n_iter += 1;

#if defined(HAVE_OPENMP)

#   pragma omp parallel for if(n_rows > THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      rk[ii] = vx[ii];

#else

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));   /* rk <- vx */

#endif

    /* Compute Vx <- Vx - (A-diag).Rk and residue. */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, rk, vx);

    res2 = 0.0;

#   pragma omp parallel for private(r) reduction(+:res2) if(n_rows > THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
      r = ad[ii] * (vx[ii]-rk[ii]);
      res2 += (r*r);
    }

#if defined(HAVE_MPI)

    if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                    _cs_sles_mpi_reduce_comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residue = sqrt(res2);

    /* Convergence test */

    cvg = _convergence_test(sles_name, var_name,
                            n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
* Block Jacobi utilities.
* Used to factorize dense 3*3 matrices.
*
* parameters:
*   ad            <-- Diagonal part of linear equation matrix
*   ad_inv        --> Inverse of the diagonal part of linear equation matrix
*   n_blocks      <-- Number of blocks
*
* returns:
*----------------------------------------------------------------------------*/


static void
_fact_lu33(const cs_real_t   *ad,
           cs_real_t         *ad_inv,
           cs_lnum_t          n_blocks)
{
  cs_lnum_t i;
  cs_real_t * _ad_inv = NULL;
  const cs_real_t * _ad = NULL;

# pragma omp parallel for private(_ad, _ad_inv) if(n_blocks > THR_MIN)
  for (i = 0; i < n_blocks; i++) {
    _ad_inv = &ad_inv[9*i];
    _ad = &ad[9*i];

    _ad_inv[0] = _ad[0];
    _ad_inv[1] = _ad[1];
    _ad_inv[2] = _ad[2];

    _ad_inv[3] = _ad[3]/_ad[0];
    _ad_inv[4] = _ad[4] - _ad_inv[3]*_ad[1];
    _ad_inv[5] = _ad[5] - _ad_inv[3]*_ad[2];

    _ad_inv[6] = _ad[6]/_ad[0];
    _ad_inv[7] = (_ad[7] - _ad_inv[6]*_ad[1])/_ad_inv[4];
    _ad_inv[8] = _ad[8] - _ad_inv[6]*_ad[2] - _ad_inv[7]*_ad_inv[5];
  }
}

/*----------------------------------------------------------------------------
* Block Jacobi utilities.
* Compute forward and backward to solve an LU 3*3 system.
*
* parameters:
*   mat           <-- 3*3*dim matrix
*   x             --> Solution
*   b             --> 1st part of RHS (c - b)
*   c             --> 2nd part of RHS (c - b)
*   aux           --> Work array
*
* returns:
*----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_lu33(const cs_real_t     mat[],
                cs_real_t           x[],
                const cs_real_t     b[],
                const cs_real_t     c[],
                cs_real_t           aux[])
{
  aux[0] = (c[0] - b[0]);
  aux[1] = (c[1] - b[1]) - aux[0]*mat[3];
  aux[2] = (c[2] - b[2]) - aux[0]*mat[6] - aux[1]*mat[7];

  x[2] = aux[2]/mat[8];
  x[1] = (aux[1] - mat[5]*x[2])/mat[4];
  x[0] = (aux[0] - mat[1]*x[1] - mat[2]*x[2])/mat[0];
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      <-- Variable name
 *   a             <-- Linear equation matrix
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   convergence   <-- Convergence information structure
 *   rhs           <-- Right hand side
 *   vx            --> System solution
 *   aux_size      <-- Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_block_3_jacobi(const char             *var_name,
                const cs_matrix_t      *a,
                cs_halo_rotation_t      rotation_mode,
                cs_sles_convergence_t  *convergence,
                const cs_real_t        *rhs,
                cs_real_t              *restrict vx,
                size_t                  aux_size,
                void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_blocks, n_rows, ii, kk, jj;
  double  res2, residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  temp[3];
  cs_real_t  *restrict ad_inv, *restrict rk, *restrict vxx;
  const cs_real_t  *restrict ad;

  unsigned n_iter = 0;

  /* Tell IBM compiler not to alias */
  #if defined(__xlc__)
  #pragma disjoint(*rhs, *vx, *ad, *ad_inv)
  #endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_JACOBI]);

  n_cols = cs_matrix_get_n_columns(a) * 3;
  n_blocks = cs_matrix_get_n_rows(a);
  n_rows = n_blocks * 3;

  ad = cs_matrix_get_diagonal(a);

  /* Allocate work arrays */

  {
    size_t  wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL
        || aux_size < (wa_size * 2) + 3*wa_size)
      BFT_MALLOC(_aux_vectors,
                 (wa_size * 2) + 3*wa_size, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;
    rk     = _aux_vectors + wa_size*3;
    vxx    = _aux_vectors + wa_size*3 + wa_size;

  }

  _fact_lu33(ad, ad_inv, n_blocks);

  cvg = 0;

  /* Current iteration */
  /*-------------------*/

  while (cvg == 0) {

    register double r;

    n_iter += 1;
    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute Vx <- Vx - (A-diag).Rk and residue. */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, rk, vxx);

    res2 = 0.0;

#   pragma omp parallel for private(jj, kk, r, temp) \
                            reduction(+:res2) if(n_blocks > THR_MIN)
    for (ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu33(ad_inv + 9*ii,
                      vx + 3*ii,
                      vxx + 3*ii,
                      rhs + 3*ii,
                      temp);
      for (jj = 0; jj < 3; jj++) {
        r = 0.0;
        for (kk = 0; kk < 3; kk++)
          r +=    ad[ii*9 + jj*3 + kk]
               * (vx[ii*3 + kk] - rk[ii*3 + kk]);
        res2 += (r*r);
      }
    }

#if defined(HAVE_MPI)

    if (cs_glob_n_ranks > 1) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                    cs_glob_mpi_comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residue = sqrt(res2);

    /* Convergence test */

    cvg = _convergence_test(sles_name, var_name,
                            n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using preconditioned Bi-CGSTAB.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name         <-- Variable name
 *   a                <-- Matrix
 *   diag_block_size  <-- Block size of element ii, ii
 *   poly_degree      <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode    <-- Halo update option for rotational periodicity
 *   convergence      <-- Convergence information structure
 *   rhs              <-- Right hand side
 *   vx               --> System solution
 *   aux_size         <-- Number of elements in aux_vectors
 *   aux_vectors      --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_bi_cgstab(const char             *var_name,
           const cs_matrix_t      *a,
           int                     diag_block_size,
           int                     poly_degree,
           cs_halo_rotation_t      rotation_mode,
           cs_sles_convergence_t  *convergence,
           const cs_real_t        *rhs,
           cs_real_t              *restrict vx,
           size_t                  aux_size,
           void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_rows, ii;
  double  _epzero = 1.e-30; /* smaller than epzero */
  double  ro_0, ro_1, alpha, beta, betam1, gamma, omega, ukres0;
  double  residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict res0, *restrict rk, *restrict pk, *restrict zk;
  cs_real_t  *restrict uk, *restrict vk, *restrict ad_inv;

  unsigned n_iter = 0;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *res0, *rk, *pk, *zk, *uk, *vk, *ad_inv)
#endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_BICGSTAB]);

  n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
  n_rows = cs_matrix_get_n_rows(a) * diag_block_size;

  /* Allocate work arrays */

  {
    size_t  n_wa = 7;
    size_t  wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;

    res0 = _aux_vectors + wa_size;
    rk = _aux_vectors + wa_size*2;
    pk = _aux_vectors + wa_size*3;
    zk = _aux_vectors + wa_size*4;
    uk = _aux_vectors + wa_size*5;
    vk = _aux_vectors + wa_size*6;

  }

  cs_matrix_copy_diagonal(a, ad_inv);
# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    pk[ii] = 0.0;
    uk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(rotation_mode, a, vx, res0);

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    res0[ii] = -res0[ii] + rhs[ii];
    rk[ii] = res0[ii];
  }

  alpha = 1.0;
  betam1 = 1.0;
  gamma = 1.0;

  cvg = 0;

  /* Current Iteration */
  /*-------------------*/

  while (cvg == 0) {

    /* Compute beta and omega;
       group dot products for new iteration's beta
       and previous iteration's residue to reduce total latency */

    if (n_iter == 0) {

      beta = _dot_product(n_rows, res0, rk);
      residue = sqrt(beta);

    }
    else {

      _dot_products_xx_xy(n_rows, rk, res0, &residue, &beta);
      residue = sqrt(residue);

      /* Convergence test */
      cvg = _convergence_test(sles_name, var_name,
                              n_iter, residue, convergence);
      if (cvg != 0)
        break;

    }

    n_iter += 1;

    if (CS_ABS(beta) < _epzero) {

      if (convergence->verbosity == 2)
        bft_printf(_("  n_iter : %5u, res_abs : %11.4e, res_nor : %11.4e\n"),
                   n_iter, residue, residue/convergence->r_norm);
      else if (convergence->verbosity > 2)
        bft_printf("   %5u %11.4e %11.4e\n",
                   n_iter, residue, residue/convergence->r_norm);

      cvg = 0;
      break;
    }

    if (CS_ABS(alpha) < _epzero) {
      bft_printf
        (_("\n\n"
           "%s [%s]:\n"
           " @@ Warning: non convergence and abort\n"
           "\n"
           "    Alpha coefficient is lower than %12.4e\n"
           "\n"
           "    The matrix cannot be considered as invertible anymore."),
         sles_name, var_name, alpha);
      cvg = -2;
      break;
    }

    omega = beta*gamma / (alpha*betam1);
    betam1 = beta;

    /* Compute pk */

#   pragma omp parallel for if(n_rows > THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      pk[ii] = rk[ii] + omega*(pk[ii] - alpha*uk[ii]);

    /* Compute zk = c.pk */

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 a,
                                 pk,
                                 zk,
                                 uk);

    /* Compute uk = A.zk */

    cs_matrix_vector_multiply(rotation_mode, a, zk, uk);

    /* Compute uk.res0 and gamma */

    ukres0 = _dot_product(n_rows, uk, res0);

    gamma = beta / ukres0;

    /* First update of vx and rk */

#   pragma omp parallel if(n_rows > THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (gamma * zk[ii]);

#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        rk[ii] -= (gamma * uk[ii]);
    }

    /* Compute zk = C.rk (zk is overwritten, vk is a working array */

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 a,
                                 rk,
                                 zk,
                                 vk);

    /* Compute vk = A.zk and alpha */

    cs_matrix_vector_multiply(rotation_mode, a, zk, vk);

    _dot_products_xx_xy(n_rows, vk, rk, &ro_1, &ro_0);

    if (ro_1 < _epzero) {
      bft_printf
        (_("\n\n"
           "%s [%s]:\n"
           " @@ Warning: non convergence and abort\n"
           "\n"
           "    The square of the norm of the descent vector\n"
           "    is lower than %12.4e\n"
           "\n"
           "    The resolution does not progress anymore."),
         sles_name, var_name, _epzero);
      cvg = -2;
      break;
    }

    alpha = ro_0 / ro_1;

    /* Final update of vx and rk */

#   pragma omp parallel if(n_rows > THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * zk[ii]);

#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        rk[ii] -= (alpha * vk[ii]);
    }

    /* Convergence test at beginning of next iteration so
       as to group dot products for better parallel performance */
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Transform using Givens rotations a system Ax=b where A is an upper
 * triangular matrix of size n*(n-1) with a lower diagonal to an equivalent
 * system A'x=b' where A' is an upper triangular matrix.
 *
 * parameters:
 *   a              <-> dense matrix (size: size a_size*a_size):
 *                      a(i,j) = a[i + j*a_size]
 *                      input: A; output: A'
 *   a_size         <-- matrix dim
 *   b              <-> pre-allocated vector of a_size elements
 *                      input: b; output: b'
 *   givens_coeff   <-> pre-allocated vector of a_size elements
 *                      input: previous Givens coefficients
 *                      output: updated Givens coefficients
 *   update_rank    <-- rank of first non null coefficient on lower diagonal
 *   end_update     <-- rank of last non null coefficient on lower diagonal
 *----------------------------------------------------------------------------*/

static void
_givens_rot_update(cs_real_t    *restrict a,
                   int                    a_size,
                   cs_real_t    *restrict b,
                   cs_real_t    *restrict givens_coeff,
                   int                    update_rank,
                   int                    end_update)
{
  int i, j;
  cs_real_t _aux;
  cs_real_t norm;

  for (i = 0; i < update_rank; ++i) {
    for (j = update_rank; j < end_update; ++j) {

      _aux =   givens_coeff[i]*a[j*a_size + i]
             + givens_coeff[i + a_size] * a[j*a_size + i+1];

      a[j*a_size + i+1] =   givens_coeff[i] * a[i+1 + j*a_size]
                          - givens_coeff[i + a_size] * a[j*a_size + i];

      a[j*a_size + i] = _aux;
    }
  }

  for (i = update_rank; i < end_update; ++i) {

    norm = pow(a[i*a_size + i], 2) + pow(a[i*a_size + i+1], 2);
    norm = sqrt(norm);

    givens_coeff[a_size + i] = a[i*a_size + i+1]/norm;
    givens_coeff[i] = a[i*a_size + i]/norm;

    b[i+1] = -b[i]*givens_coeff[a_size + i];
    b[i] = b[i]*givens_coeff[i];

    for (j = i; j < end_update; j++) {

      _aux =   givens_coeff[i] * a[j*a_size + i]
             + givens_coeff[a_size+i] * a[j*a_size + i+1];
      if (j == i)
        a[j*a_size+i+1] = 0;
      else
        a[i+1 + j*a_size] =   givens_coeff[i]*a[i+1 + j*a_size]
                            - givens_coeff[a_size + i]*a[i + j*a_size];

      a[j*a_size + i] = _aux;
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute solution of Ax = b where A is an upper triangular matrix.
 *
 * As the system solved by GMRES will grow with iteration number, we
 * preallocate a allocated size, and solve only the useful part:
 *
 *   |       |       |   |x1|   |b1|
 *   |   A   |  pad  |   |x2|   |b2|
 *   |_______|       |   |x3|   |b3|
 *   |               | * |p | = |p |
 *   |     pad       |   |a |   |a |
 *   |               |   |d |   |d |
 *
 * parameters:
 *   a          <-- dense upper triangular matrix A (size: glob_size*glob_size)
 *                 a(i,j) = a[i + j*glob_size]
 *   a_size     <-- system size
 *   alloc_size <-- a_size + halo size
 *   b          <-- pre-allocated vector of a_size elements
 *   x          --> system solution, pre-allocated vector of a_size elements
 *
 * returns:
 *   0 in case of success, 1 in case of zero-pivot.
 *----------------------------------------------------------------------------*/

static int
_solve_diag_sup_halo(cs_real_t  *restrict a,
                     int                  a_size,
                     int                  alloc_size,
                     cs_real_t  *restrict b,
                     cs_real_t  *restrict x)
{
  int i, j;

  for (i = a_size - 1; i > -1; i--) {

    x[i] = b[i];

    for (j = i + 1; j < a_size; j++)
      x[i] = x[i] - a[j*alloc_size + i]*x[j];

    x[i] /= a[i*alloc_size + i];
  }

  return 0; /* We should add a check for zero-pivot */
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned GMRES.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      <-- Variable name
 *   a             <-- Matrix
 *   poly_degree   <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   convergence   <-- Convergence information structure
 *   rhs           <-- Right hand side
 *   vx            --> System solution
 *   aux_size      <-- Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_gmres(const char             *var_name,
       const cs_matrix_t      *a,
       int                     poly_degree,
       cs_halo_rotation_t      rotation_mode,
       cs_sles_convergence_t  *convergence,
       const cs_real_t        *rhs,
       cs_real_t              *restrict vx,
       size_t                  aux_size,
       void                   *aux_vectors)
{
  int cvg;
  const char *sles_name;
  int check_freq, l_iter, l_old_iter, scaltest;
  int krylov_size, _krylov_size;
  cs_lnum_t  n_cols, n_rows, ii, kk, jj;
  double    beta, dot_prod, residue, _residue, epsi;
  cs_real_t  *_aux_vectors;
  cs_real_t *restrict _krylov_vectors, *restrict _h_matrix;
  cs_real_t *restrict _givens_coeff, *restrict _beta;
  cs_real_t *restrict dk, *restrict gk, *restrict ad_inv;
  cs_real_t *restrict bk, *restrict fk, *restrict krk;

  int krylov_size_max = 75;
  unsigned n_iter = 0;

#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *_krylov_vectors, *_h_matrix, \
  *_givens_coeff, *_beta, *dk, *gk, *bk, *fk, *ad_inv, *krk)
#endif

  sles_name = _(cs_sles_type_name[CS_SLES_GMRES]);

  /* Preliminary calculations */
  /*--------------------------*/

  n_cols = cs_matrix_get_n_columns(a);
  n_rows = cs_matrix_get_n_rows(a);

  /* Allocate work arrays */

  krylov_size =  (krylov_size_max < (int)sqrt(n_rows)*1.5) ?
                  krylov_size_max : (int)sqrt(n_rows)*1.5 + 1;

#if defined(HAVE_MPI)
  if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
    MPI_Allreduce(&krylov_size,
                  &_krylov_size,
                  1,
                  MPI_INT,
                  MPI_MIN,
                  _cs_sles_mpi_reduce_comm);
    krylov_size = _krylov_size;
  }
#endif

  check_freq = (int)(krylov_size/10) + 1;
  epsi = 1.e-15;
  scaltest = 0;

  {
    size_t _aux_size;
    size_t  n_wa = 5;
    size_t  wa_size = n_cols < krylov_size? krylov_size:n_cols;

    wa_size = CS_SIMD_SIZE(wa_size);
    _aux_size =   wa_size*n_wa
                + (krylov_size-1)*(n_rows + krylov_size) + 3*krylov_size;

    if (aux_vectors == NULL || aux_size < _aux_size)
      BFT_MALLOC(_aux_vectors, _aux_size, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    dk = _aux_vectors;
    gk = _aux_vectors + wa_size;
    ad_inv = _aux_vectors + 2*wa_size;
    bk = _aux_vectors + 3*wa_size;
    fk = _aux_vectors + 4*wa_size;
    _krylov_vectors = _aux_vectors + 5*wa_size;
    _h_matrix = _aux_vectors + 5*wa_size + (krylov_size - 1)*n_rows;
    _givens_coeff =   _aux_vectors + 5*wa_size
                    + (krylov_size - 1)*(n_rows + krylov_size);
    _beta =   _aux_vectors + 5*wa_size
            + (krylov_size - 1)*(n_rows + krylov_size) + 2*krylov_size;
  }

  for (ii = 0; ii < krylov_size*(krylov_size - 1); ii++)
    _h_matrix[ii] = 0.;

  cs_matrix_copy_diagonal(a, ad_inv);

# pragma omp parallel for if(n_rows > THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1./ad_inv[ii];

  cvg = 0;

  while (cvg == 0) {

    /* compute  rk <- a*vx (vx = x0) */

    cs_matrix_vector_multiply(rotation_mode, a, vx, dk);

    /* compute  rk <- rhs - rk (r0 = b-A*x0) */

#   pragma omp parallel for if(n_rows > THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = rhs[ii] - dk[ii];

    /* beta = ||r0|| */
    beta = sqrt(_dot_product(n_rows, dk, dk));
    dot_prod = beta;

    _beta[0] = beta;
    for (ii = 1; ii < krylov_size; ii++)
      _beta[ii] = 0.;

    /* Lap */

    l_iter = 0;
    l_old_iter = 0;

    for (ii = 0; ii < krylov_size - 1; ii++) {

      /* krk = k-ieth col of _krylov_vector = vi */
      krk = _krylov_vectors + ii*n_rows;

#     pragma omp parallel for if(n_rows > THR_MIN)
      for (jj = 0; jj < n_rows; jj++)
        krk[jj] = dk[jj]/dot_prod;

      _polynomial_preconditionning(n_rows,
                                   poly_degree,
                                   rotation_mode,
                                   ad_inv,
                                   a,
                                   krk,
                                   gk,
                                   dk);

      /* compute w = dk <- A*vj */

      cs_matrix_vector_multiply(rotation_mode, a, gk, dk);

      for (jj = 0; jj < ii + 1; jj++) {

        /* compute h(k,i) = <w,vi> = <dk,vi> */
        _h_matrix[ii*krylov_size + jj]
          = _dot_product(n_rows, dk, (_krylov_vectors + jj*n_rows));

        /* compute w = dk <- w - h(i,k)*vi */

        cs_axpy(n_rows,
                -_h_matrix[ii*krylov_size+jj],
                (_krylov_vectors + jj*n_rows),
                dk);

      }

      /* compute h(i+1,i) = sqrt<w,w> */
      dot_prod = sqrt(_dot_product(n_rows, dk, dk));
      _h_matrix[ii*krylov_size + ii + 1] = dot_prod;

      if (dot_prod < epsi) scaltest = 1;

      if (   (l_iter + 1)%check_freq == 0
          || l_iter == krylov_size - 2
          || scaltest == 1) {

          /* H matrix to diagonal sup matrix */

        _givens_rot_update(_h_matrix,
                           krylov_size,
                           _beta,
                           _givens_coeff,
                           l_old_iter,
                           l_iter + 1);

        l_old_iter = l_iter + 1;

        /* solve diag sup system */
        _solve_diag_sup_halo(_h_matrix, l_iter + 1, krylov_size, _beta, gk);

#       pragma omp parallel for private(kk) if(n_rows > THR_MIN)
        for (jj = 0; jj < n_rows; jj++) {
          fk[jj] = 0.0;
          for (kk = 0; kk <= l_iter; kk++)
            fk[jj] += _krylov_vectors[kk*n_rows + jj] * gk[kk];
        }

        _polynomial_preconditionning(n_rows,
                                     poly_degree,
                                     rotation_mode,
                                     ad_inv,
                                     a,
                                     fk,
                                     gk,
                                     bk);


#       pragma omp parallel for if(n_rows > THR_MIN)
        for (jj = 0; jj < n_rows; jj++)
          fk[jj] = vx[jj] + gk[jj];

        cs_matrix_vector_multiply(rotation_mode, a, fk, bk);

        /* compute residue = | Ax - b |_1 */

        residue = 0.;
#       pragma omp parallel for reduction(+:residue) if(n_rows > THR_MIN)
        for (jj = 0; jj < n_rows; jj++)
          residue += pow(rhs[jj] - bk[jj], 2);

#if defined(HAVE_MPI)

        if (_cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
          MPI_Allreduce(&residue, &_residue, 1, MPI_DOUBLE, MPI_SUM,
                        _cs_sles_mpi_reduce_comm);
          residue = _residue;
        }

#endif

        residue = sqrt(residue);

        cvg = _convergence_test(sles_name,
                                var_name,
                                n_iter,
                                residue,
                                convergence);

      }

      n_iter++;
      l_iter++;

      if (cvg == 1 || l_iter == krylov_size - 1 || scaltest == 1) {
#       pragma omp parallel for if(n_rows > THR_MIN)
        for (jj = 0; jj < n_rows; jj++)
          vx[jj] = fk[jj];
        break;
      }
    }
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Output post-processing data for failed system convergence.
 *
 * parameters:
 *   n_vals        <-- Size of val and val_type array
 *   val           <-> Values to post-process (set to 0 on output if not
 *                     normal floating-point values)
 *   val_type      --> 0: normal values, 1: infinite, 2: Nan
 *
 * returns:
 *   number of non-normal values
 *----------------------------------------------------------------------------*/

static size_t
_value_type(size_t      n_vals,
            cs_real_t   val[],
            cs_real_t   val_type[])
{
  size_t ii;
  size_t retval = 0;

#if (__STDC_VERSION__ >= 199901L)

  for (ii = 0; ii < n_vals; ii++) {

    int v_type = fpclassify(val[ii]);

    if (v_type == FP_INFINITE) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else if (v_type == FP_NAN) {
      val[ii] = 0.;
      val_type[ii] = 2;
      retval += 1;
    }

    else if (val[ii] > 1.e38 || val[ii] < -1.e38) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else
      val_type[ii] = 0;
  }

#else

  for (ii = 0; ii < n_vals; ii++) {

    if (val[ii] != val[ii]) { /* Test for NaN with IEEE 754 arithmetic */
      val[ii] = 0.;
      val_type[ii] = 2;
      retval += 1;
    }

    else if (val[ii] > 1.e38 || val[ii] < -1.e38) {
      val[ii] = 0.;
      val_type[ii] = 1;
      retval += 1;
    }

    else
      val_type[ii] = 0;
  }

#endif

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute per-cell residual for Ax = b.
 *
 * parameters:
 *   diag_block_size  <-- Block sizes for diagonal
 *   rotation_mode    <-- Halo update option for rotational periodicity
 *   a                <-- Linear equation matrix
 *   rhs              <-- Right hand side
 *   vx               <-> Current system solution
 *   res              --> Residual
 *----------------------------------------------------------------------------*/

static void
_cell_residual(const int           *diag_block_size,
               cs_halo_rotation_t   rotation_mode,
               const cs_matrix_t   *a,
               const cs_real_t      rhs[],
               cs_real_t            vx[],
               cs_real_t            res[])
{
  cs_int_t ii;

  const cs_int_t n_vals = cs_glob_mesh->n_cells * diag_block_size[1];

  cs_matrix_vector_multiply(rotation_mode, a, vx, res);

# pragma omp parallel for if(n_vals > THR_MIN)
  for (ii = 0; ii < n_vals; ii++)
    res[ii] = fabs(res[ii] - rhs[ii]);
}

/*----------------------------------------------------------------------------
 * Sparse linear system resolution for non-interleaved non-symmetric matrixes.
 *
 * This function does not have the full functionnality of cs_sles_solve.c,
 * as it should only be used for "legacy" non-interleaved non-symmetric
 * matrixes: it may not be used inside a multigrid solver, and always updates
 * statistics.
 *
 * Note that in most cases (if the right-hand side is not already zero
 * within convergence criteria), coefficients are assigned to matrixes
 * then released by this function, so coefficients need not be assigned
 * prior to this call, and will have been released upon returning.
 *
 * parameters:
 *   var_name      <-- Variable name
 *   solver_type   <-- Type of solver (PCG, Jacobi, ...)
 *   a             <-- Matrix
 *   poly_degree   <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   verbosity     <-- Verbosity level
 *   n_max_iter    <-- Maximum number of iterations
 *   precision     <-- Precision limit
 *   r_norm        <-- Residue normalization
 *   n_iter        --> Number of iterations
 *   residue       <-> Residue
 *   rhs           <-- Right hand side
 *   vx            --> System solution
 *   aux_size      <-- Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

static int
_solve_ni(const char          *var_name,
          cs_sles_type_t       solver_type,
          const cs_matrix_t   *a,
          int                  poly_degree,
          cs_halo_rotation_t   rotation_mode,
          int                  verbosity,
          int                  n_max_iter,
          double               precision,
          double               r_norm,
          int                 *n_iter,
          double              *residue,
          const cs_real_t     *rhs,
          cs_real_t           *vx,
          size_t               aux_size,
          void                *aux_vectors)
{
  cs_int_t  n_rows;
  cs_timer_t t0, t1;

  unsigned _n_iter = 0;
  int cvg = 0;
  cs_sles_convergence_t  convergence;

  cs_sles_info_t *sles_info = NULL;

  n_rows = cs_matrix_get_n_rows(a);

  t0 = cs_timer_time();
  sles_info = _find_or_add_system(var_name, solver_type);

  /* Initialize number of iterations and residue,
     check for immediate return,
     and solve sparse linear system */

  *n_iter = 0;

  if (cs_sles_needs_solving(var_name,
                            _(cs_sles_type_name[solver_type]),
                            n_rows,
                            verbosity,
                            r_norm,
                            residue,
                            rhs) != 0) {

    /* Solve sparse linear system */

    _convergence_init(&convergence,
                      _(cs_sles_type_name[solver_type]),
                      var_name,
                      verbosity,
                      n_max_iter,
                      precision,
                      r_norm,
                      *residue);

    switch (solver_type) {
    case CS_SLES_PCG:
    case CS_SLES_PCG_SR:
      cvg = _conjugate_gradient(var_name,
                                a,
                                1,
                                poly_degree,
                                rotation_mode,
                                &convergence,
                                rhs,
                                vx,
                                aux_size,
                                aux_vectors);
      break;
    case CS_SLES_JACOBI:
      cvg = _jacobi(var_name,
                    a,
                    rotation_mode,
                    &convergence,
                    rhs,
                    vx,
                    aux_size,
                    aux_vectors);
      break;
    case CS_SLES_BICGSTAB:
      cvg = _bi_cgstab(var_name,
                       a,
                       1,
                       poly_degree,
                       rotation_mode,
                       &convergence,
                       rhs,
                       vx,
                       aux_size,
                       aux_vectors);
      break;
    case CS_SLES_GMRES:
      cvg = _gmres(var_name,
                   a,
                   poly_degree,
                   rotation_mode,
                   &convergence,
                   rhs,
                   vx,
                   aux_size,
                   aux_vectors);
      break;
    default:
      break;
    }

    /* Update return values */

    _n_iter = convergence.n_iterations;

    *n_iter = convergence.n_iterations;
    *residue = convergence.residue;
  }

  t1 = cs_timer_time();

  if (sles_info->n_calls == 0)
    sles_info->n_iterations_min = _n_iter;

  sles_info->n_calls += 1;

  if (sles_info->n_iterations_min > _n_iter)
    sles_info->n_iterations_min = _n_iter;
  if (sles_info->n_iterations_max < _n_iter)
    sles_info->n_iterations_max = _n_iter;

  sles_info->n_iterations_last = _n_iter;
  sles_info->n_iterations_tot += _n_iter;

  cs_timer_counter_add_diff(&(sles_info->t_tot), &t0, &t1);

  return cvg;
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(reslin, RESLIN)
(
 const char       *cname,     /* <-- variable name */
 const cs_int_t   *lname,     /* <-- variable name length */
 const cs_int_t   *ncelet,    /* <-- Number of cells, halo included */
 const cs_int_t   *ncel,      /* <-- Number of local cells */
 const cs_int_t   *nfac,      /* <-- Number of faces */
 const cs_int_t   *isym,      /* <-- Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ilved,     /* <-- Interleaved indicator  */
                              /*     1: interleaved; 2: not interleaved */
 const cs_int_t   *ibsize,    /* <-- Block size of element ii, ii */
 const cs_int_t   *iesize,    /* <-- Block size of element ij */
 const cs_int_t   *ireslp,    /* <-- Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab,
                                     3: gmres, 10: pcg_sr */
 const cs_int_t   *ipol,      /* <-- Preconditioning polynomial degree
                                     (0: diagonal) */
 const cs_int_t   *nitmap,    /* <-- Number of max iterations */
 const cs_int_t   *iinvpe,    /* <-- Indicator to cancel increments
                                     in rotational periodicty (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* <-- Verbosity level */
 cs_int_t         *niterf,    /* --> Number of iterations done */
 const cs_real_t  *epsilp,    /* <-- Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* <-- Residue normalization */
 cs_real_t        *residu,    /* --> Final non normalized residue */
 const cs_int_t   *ifacel,    /* <-- Face -> cell connectivity  */
 const cs_real_t  *dam,       /* <-- Matrix diagonal */
 const cs_real_t  *xam,       /* <-- Matrix extra-diagonal terms */
 const cs_real_t  *rhs,       /* <-- System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
)
{
  char *var_name;
  cs_sles_type_t type;

  int cvg = 0;
  int n_iter = *niterf;
  int diag_block_size[4] = {1, 1, 1, 1};
  int extra_diag_block_size[4] = {1, 1, 1, 1};
  bool symmetric = (*isym == 1) ? true : false;
  bool interleaved = (*ilved == 1) ? true : false;
  cs_halo_rotation_t rotation_mode = CS_HALO_ROTATION_COPY;

  cs_matrix_t *a = cs_glob_matrix_default;

  assert(*ncelet >= *ncel);
  assert(*nfac > 0);
  assert(ifacel != NULL);

  if (*iinvpe == 2)
    rotation_mode = CS_HALO_ROTATION_ZERO;
  else if (*iinvpe == 3)
    rotation_mode = CS_HALO_ROTATION_IGNORE;

  if (*ibsize > 1) {
    diag_block_size[0] = *ibsize;
    diag_block_size[1] = *ibsize;
    diag_block_size[2] = *ibsize;
    diag_block_size[3] = (*ibsize)*(*ibsize);
  }

  if (*iesize > 1) {
    extra_diag_block_size[0] = *iesize;
    extra_diag_block_size[1] = *iesize;
    extra_diag_block_size[2] = *iesize;
    extra_diag_block_size[3] = (*iesize)*(*iesize);
  }

  var_name = cs_base_string_f_to_c_create(cname, *lname);

  switch ((int)(*ireslp)) {
  case 0:
    type = CS_SLES_PCG;
    break;
  case 1:
    type = CS_SLES_JACOBI;
    break;
  case 2:
    type = CS_SLES_BICGSTAB;
    break;
  case 3:
    type = CS_SLES_GMRES;
    break;
  case 10:
    type = CS_SLES_PCG_SR;
    break;
  default:
    type = CS_SLES_N_TYPES;
    assert(0);
  }

  if (interleaved || symmetric) {
    cs_matrix_set_coefficients(a,
                               symmetric,
                               diag_block_size,
                               extra_diag_block_size,
                               dam,
                               xam);
    cvg = cs_sles_solve(var_name,
                        type,
                        true,
                        a,
                        *ipol,
                        rotation_mode,
                        *iwarnp,
                        *nitmap,
                        *epsilp,
                        *rnorm,
                        &n_iter,
                        residu,
                        rhs,
                        vx,
                        0,
                        NULL);
  }
  else {
    cs_matrix_set_coefficients_ni(a, false, dam, xam);
    cvg = _solve_ni(var_name,
                    type,
                    a,
                    *ipol,
                    rotation_mode,
                    *iwarnp,
                    *nitmap,
                    *epsilp,
                    *rnorm,
                    &n_iter,
                    residu,
                    rhs,
                    vx,
                    0,
                    NULL);
  }

  *niterf = n_iter;

  /* If divergence is detected, try diagnostics and abort */

  if (cvg == -2) {

    int mesh_id = cs_post_init_error_writer_cells();

    cs_sles_post_error_output_def(var_name,
                                  mesh_id,
                                  rotation_mode,
                                  a,
                                  rhs,
                                  vx);

    cs_post_finalize();

    bft_error(__FILE__, __LINE__, 0,
              _("%s: error (divergence) solving for %s"),
              _(cs_sles_type_name[type]), var_name);
  }

  cs_matrix_release_coefficients(a);

  cs_base_string_f_to_c_free(&var_name);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_initialize(void)
{
  cs_mesh_t  *mesh = cs_glob_mesh;

  assert(mesh != NULL);

   cs_glob_sles_native_matrix_struct
    = cs_matrix_structure_create(CS_MATRIX_NATIVE,
                                 true,
                                 mesh->n_cells,
                                 mesh->n_cells_with_ghosts,
                                 mesh->n_i_faces,
                                 mesh->global_cell_num,
                                 mesh->i_face_cells,
                                 mesh->halo,
                                 mesh->i_face_numbering);

  cs_glob_sles_native_matrix
    = cs_matrix_create(cs_glob_sles_native_matrix_struct);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    _cs_sles_mpi_reduce_comm = cs_glob_mpi_comm;
  else
    _cs_sles_mpi_reduce_comm = MPI_COMM_NULL;
#endif
}

/*----------------------------------------------------------------------------
 * Finalize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_finalize(void)
{
  int ii;

  /* Free system info */

  for (ii = 0; ii < cs_glob_sles_n_systems; ii++) {
    _sles_info_dump(cs_glob_sles_systems[ii]);
    _sles_info_destroy(&(cs_glob_sles_systems[ii]));
  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  BFT_FREE(cs_glob_sles_systems);

  cs_glob_sles_n_systems = 0;
  cs_glob_sles_n_max_systems = 0;

  /* Free matrix structures */

  cs_matrix_destroy(&cs_glob_sles_native_matrix);

  cs_matrix_structure_destroy(&cs_glob_sles_native_matrix_struct);
 }

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Set MPI communicator for dot products.
 *----------------------------------------------------------------------------*/

void
cs_sles_set_mpi_reduce_comm(MPI_Comm comm)
{
  _cs_sles_mpi_reduce_comm = comm;

  if (comm != cs_glob_mpi_comm)
    cs_halo_set_use_barrier(0);
  else {
    cs_halo_set_use_barrier(1);
    if (cs_glob_n_ranks < 2)
      _cs_sles_mpi_reduce_comm = MPI_COMM_NULL;
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Test if a general sparse linear system needs solving or if the right-hand
 * side is already zero within convergence criteria.
 *
 * The computed residue is also updated;
 *
 * parameters:
 *   var_name      <-- Variable name
 *   solver_name   <-- Name of solver
 *   n_rows        <-- Number of (non ghost) rows in rhs
 *   verbosity     <-- Verbosity level
 *   r_norm        <-- Residue normalization
 *   residue       <-> Residue
 *   rhs           <-- Right hand side
 *
 * returns:
 *   1 if solving is required, 0 if the rhs is already zero within tolerance
 *   criteria (precision of residue normalization)
 *----------------------------------------------------------------------------*/

int
cs_sles_needs_solving(const char        *var_name,
                      const char        *solver_name,
                      cs_int_t           n_rows,
                      int                verbosity,
                      double             r_norm,
                      double            *residue,
                      const cs_real_t   *rhs)
{
  int retval = 1;

  /* Initialize residue, check for immediate return */

  *residue = sqrt(_dot_product(n_rows, rhs, rhs));

#if defined(HAVE_MPI)
  if (_cs_sles_mpi_reduce_comm != cs_glob_mpi_comm )
    MPI_Bcast(residue, 1, MPI_DOUBLE, 0, cs_glob_mpi_comm);
#endif

  if (r_norm <= EPZERO || *residue <= EPZERO) {
    if (verbosity > 1)
      bft_printf(_("%s [%s]:\n"
                   "  immediate exit; r_norm = %11.4e, residual = %11.4e\n"),
                 solver_name, var_name, r_norm, *residue);
    retval = 0;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * General sparse linear system resolution.
 *
 * parameters:
 *   var_name          <-- Variable name
 *   solver_type       <-- Type of solver (PCG, Jacobi, ...)
 *   update_stats      <-- Automatic solver statistics indicator
 *   a                 <-- Matrix
 *   poly_degree       <-- Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode     <-- Halo update option for rotational periodicity
 *   verbosity         <-- Verbosity level
 *   n_max_iter        <-- Maximum number of iterations
 *   precision         <-- Precision limit
 *   r_norm            <-- Residue normalization
 *   n_iter            --> Number of iterations
 *   residue           <-> Residue
 *   rhs               <-- Right hand side
 *   vx                --> System solution
 *   aux_size          <-- Number of elements in aux_vectors
 *   aux_vectors       --- Optional working area (allocation otherwise)
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached, -2 if divergence is detected.
 *----------------------------------------------------------------------------*/

int
cs_sles_solve(const char          *var_name,
              cs_sles_type_t       solver_type,
              bool                 update_stats,
              const cs_matrix_t   *a,
              int                  poly_degree,
              cs_halo_rotation_t   rotation_mode,
              int                  verbosity,
              int                  n_max_iter,
              double               precision,
              double               r_norm,
              int                 *n_iter,
              double              *residue,
              const cs_real_t     *rhs,
              cs_real_t           *vx,
              size_t               aux_size,
              void                *aux_vectors)
{
  cs_int_t  n_rows;
  cs_timer_t t0, t1;

  unsigned _n_iter = 0;
  int cvg = 0;
  int _diag_block_size = 1;
  cs_sles_convergence_t  convergence;

  cs_sles_info_t *sles_info = NULL;

  const int *diag_block_size = cs_matrix_get_diag_block_size(a);

  n_rows = cs_matrix_get_n_rows(a);

  if (update_stats == true) {
    t0 = cs_timer_time();
    sles_info = _find_or_add_system(var_name, solver_type);
  }

  if (diag_block_size != NULL) {
    assert(diag_block_size[0] == diag_block_size[1]);
    _diag_block_size = diag_block_size[0];
  }

  /* Initialize number of iterations and residue,
     check for immediate return,
     and solve sparse linear system */

  *n_iter = 0;

  if (cs_sles_needs_solving(var_name,
                            _(cs_sles_type_name[solver_type]),
                            _diag_block_size*n_rows,
                            verbosity,
                            r_norm,
                            residue,
                            rhs) != 0) {

    /* Solve sparse linear system */

    _convergence_init(&convergence,
                      _(cs_sles_type_name[solver_type]),
                      var_name,
                      verbosity,
                      n_max_iter,
                      precision,
                      r_norm,
                      *residue);

    /* Only call solver for "active" ranks */

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks < 2 || _cs_sles_mpi_reduce_comm != MPI_COMM_NULL) {
#endif

      switch (solver_type) {
      case CS_SLES_PCG:
        cvg = _conjugate_gradient(var_name,
                                  a,
                                  _diag_block_size,
                                  poly_degree,
                                  rotation_mode,
                                  &convergence,
                                  rhs,
                                  vx,
                                  aux_size,
                                  aux_vectors);
        break;
      case CS_SLES_PCG_SR:
        cvg = _conjugate_gradient_sr(var_name,
                                     a,
                                     _diag_block_size,
                                     poly_degree,
                                     rotation_mode,
                                       &convergence,
                                     rhs,
                                     vx,
                                     aux_size,
                                     aux_vectors);
        break;
      case CS_SLES_JACOBI:
        if (_diag_block_size == 1)
          cvg = _jacobi(var_name,
                        a,
                        rotation_mode,
                        &convergence,
                        rhs,
                        vx,
                        aux_size,
                        aux_vectors);
        else if (_diag_block_size == 3)
          cvg = _block_3_jacobi(var_name,
                                a,
                                rotation_mode,
                                &convergence,
                                rhs,
                                vx,
                                aux_size,
                                aux_vectors);
        break;
      case CS_SLES_BICGSTAB:
        cvg = _bi_cgstab(var_name,
                         a,
                         _diag_block_size,
                         poly_degree,
                         rotation_mode,
                         &convergence,
                         rhs,
                         vx,
                         aux_size,
                         aux_vectors);
        break;
      case CS_SLES_GMRES:
        if (_diag_block_size == 1)
          cvg = _gmres(var_name,
                       a,
                       poly_degree,
                       rotation_mode,
                       &convergence,
                       rhs,
                       vx,
                       aux_size,
                       aux_vectors);
        else
          bft_error
            (__FILE__, __LINE__, 0,
             _("GMRES not supported with block_size > 1 (velocity coupling)."));
        break;
      default:
        bft_error
          (__FILE__, __LINE__, 0,
           _("Resolution of linear equation on \"%s\"\n"
             "with solver type %d, which is not defined)."),
           var_name, (int)solver_type);
        break;
      }

#if defined(HAVE_MPI)
    }
#endif

    /* Broadcast convergence info from "active" ranks to others*/

#if defined(HAVE_MPI)
    if (_cs_sles_mpi_reduce_comm != cs_glob_mpi_comm) {
      /* cvg is signed, so shift (with some margin) before copy to unsigned. */
      unsigned buf[2] = {cvg+10, convergence.n_iterations};
      MPI_Bcast(buf, 2, MPI_UNSIGNED, 0, cs_glob_mpi_comm);
      MPI_Bcast(&convergence.residue, 1, MPI_DOUBLE, 0,
                cs_glob_mpi_comm);
      cvg = buf[0] - 10;
      convergence.n_iterations = buf[1];
    }
#endif

    /* Update return values */

    _n_iter = convergence.n_iterations;

    *n_iter = convergence.n_iterations;
    *residue = convergence.residue;
  }

  if (update_stats == true) {

    t1 = cs_timer_time();

    if (sles_info->n_calls == 0)
      sles_info->n_iterations_min = _n_iter;

    sles_info->n_calls += 1;

    if (sles_info->n_iterations_min > _n_iter)
      sles_info->n_iterations_min = _n_iter;
    if (sles_info->n_iterations_max < _n_iter)
      sles_info->n_iterations_max = _n_iter;

    sles_info->n_iterations_last = _n_iter;
    sles_info->n_iterations_tot += _n_iter;

    cs_timer_counter_add_diff(&(sles_info->t_tot), &t0, &t1);
  }

  return cvg;
}

/*----------------------------------------------------------------------------
 * Output default post-processing data for failed system convergence.
 *
 * parameters:
 *   var_name         <-- Variable name
 *   mesh_id          <-- id of error output mesh, or 0 if none
 *   rotation_mode    <-- Halo update option for rotational periodicity
 *   a                <-- Linear equation matrix
 *   ax               <-- Non-diagonal part of linear equation matrix
 *   rhs              <-- Right hand side
 *   vx               <-> Current system solution
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_def(const char          *var_name,
                              int                  mesh_id,
                              cs_halo_rotation_t   rotation_mode,
                              const cs_matrix_t   *a,
                              const cs_real_t     *rhs,
                              cs_real_t           *vx)
{
  if (mesh_id != 0) {

    const cs_mesh_t *mesh = cs_glob_mesh;

    char base_name[32], val_name[32];

    int val_id;
    const cs_int_t n_cells = mesh->n_cells;
    const int *diag_block_size = cs_matrix_get_diag_block_size(a);

    cs_real_t *val;

    BFT_MALLOC(val, mesh->n_cells_with_ghosts*diag_block_size[1], cs_real_t);

    for (val_id = 0; val_id < 5; val_id++) {

      switch(val_id) {

      case 0:
        strcpy(base_name, "Diag");
        cs_matrix_copy_diagonal(a, val);
        break;

      case 1:
        strcpy(base_name, "RHS");
        memcpy(val, rhs, n_cells*diag_block_size[1]*sizeof(cs_real_t));
        break;

      case 2:
        strcpy(base_name, "X");
        memcpy(val, vx, n_cells*diag_block_size[1]*sizeof(cs_real_t));
        break;

      case 3:
        strcpy(base_name, "Residual");
        _cell_residual(diag_block_size,
                       rotation_mode,
                       a,
                       rhs,
                       vx,
                       val);
        break;

      case 4:
        strcpy(base_name, "Diag_Dom");
        cs_matrix_diag_dominance(a, val);
        break;

      }

      if (strlen(var_name) + strlen(base_name) < 31) {
        strcpy(val_name, base_name);
        strcat(val_name, "_");
        strcat(val_name, var_name);
      }
      else {
        strncpy(val_name, base_name, 31);
        val_name[31] = '\0';
      }

      assert(diag_block_size[0] == diag_block_size[1]); /* for now */

      cs_sles_post_error_output_var(val_name,
                                    mesh_id,
                                    diag_block_size[1],
                                    val);
    }

    BFT_FREE(val);
  }
}

/*----------------------------------------------------------------------------
 * Output post-processing variable for failed system convergence.
 *
 * parameters:
 *   var_name        <-- Variable name
 *   mesh_id         <-- id of error output mesh, or 0 if none
 *   diag_block_size <-- Block size for diagonal
 *   var             <-- Variable values
 *----------------------------------------------------------------------------*/

void
cs_sles_post_error_output_var(const char  *var_name,
                              int          mesh_id,
                              int          diag_block_size,
                              cs_real_t   *var)
{
  if (mesh_id != 0) {

    int _diag_block_size[4] = {1, 1, 1, 1};

    const cs_mesh_t *mesh = cs_glob_mesh;

    size_t n_non_norm;
    const cs_int_t n_cells = mesh->n_cells;

    cs_real_t *val_type;

    assert(_diag_block_size[0] == _diag_block_size[1]); /* no padding */

    if (diag_block_size > 1) {
      int i;
      for (i = 0; i < 3; i++) /* will become false if padding is used */
        _diag_block_size[i] = diag_block_size;
      _diag_block_size[3] = diag_block_size*diag_block_size;
    }

    BFT_MALLOC(val_type, _diag_block_size[1]*n_cells, cs_real_t);

    n_non_norm = _value_type(_diag_block_size[1]*n_cells, var, val_type);

    cs_post_write_var(mesh_id,
                      var_name,
                      _diag_block_size[0],
                      true, /* interlace */
                      true, /* use parents */
                      CS_POST_TYPE_cs_real_t,
                      -1,
                      0.0,
                      var,
                      NULL,
                      NULL);

    if (n_non_norm > 0) {

      char type_name[32];
      size_t l = strlen(var_name);

      if (l > 31)
        l = 31;

      l -= strlen("_fp_type");

      strncpy(type_name, var_name, l);
      type_name[l] = '\0';

      strcat(type_name, "_fp_type");

      cs_post_write_var(mesh_id,
                        type_name,
                        _diag_block_size[0],
                        true, /* interlace */
                        true, /* use parents */
                        CS_POST_TYPE_cs_real_t,
                        -1,
                        0.0,
                        val_type,
                        NULL,
                        NULL);

    }

    BFT_FREE(val_type);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
