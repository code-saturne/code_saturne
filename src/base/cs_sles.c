/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
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
 * Sparse Linear Equation Solvers
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(_CS_HAVE_MPI)
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
#include "cs_prototypes.h"
#include "cs_mesh.h"
#include "cs_matrix.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define EPZERO  1.E-12
#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

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

  double               wt_tot;             /* Total wall-clock time used */
  double               cpu_tot;            /* Total (local) CPU used */

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

cs_matrix_t *cs_glob_sles_base_matrix = NULL;
cs_matrix_t *cs_glob_sles_native_matrix = NULL;

/* Sparse linear equation solver type names */

/* Short names for matrix types */

const char *cs_sles_type_name[] = {N_("Gradient conjugué"),
                                   N_("Jacobi"),
                                   N_("Bi-CGstab")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to new linear system info structure.
 *
 * parameters:
 *   name --> system name
 *   type --> resolution method
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

  new_info->wt_tot = 0.0;
  new_info->cpu_tot = 0.0;

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

  bft_printf(_("\n"
               "Bilan des résolutions pour \"%s\" (%s) :\n\n"
               "  Nombre d'appels :                 %d\n"
               "  Nombre d'itérations minimal :     %d\n"
               "  Nombre d'itérations maximal :     %d\n"
               "  Nombre d'itérations moyen :       %d\n"
               "  Temps écoulé cumulé :             %12.3f\n"),
             this_info->name, cs_sles_type_name[this_info->type],
             n_calls, n_it_min, n_it_max, n_it_mean,
             this_info->wt_tot);

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_nbr > 1) {

    double cpu_min, cpu_max, cpu_tot;
    double cpu_loc = this_info->cpu_tot;

    MPI_Allreduce(&cpu_loc, &cpu_min, 1, MPI_DOUBLE, MPI_MIN,
                  cs_glob_base_mpi_comm);
    MPI_Allreduce(&cpu_loc, &cpu_max, 1, MPI_DOUBLE, MPI_MAX,
                  cs_glob_base_mpi_comm);
    MPI_Allreduce(&cpu_loc, &cpu_tot, 1, MPI_DOUBLE, MPI_SUM,
                  cs_glob_base_mpi_comm);

    bft_printf(_("  Temps CPU cumulé local min :      %12.3f\n"
                 "  Temps CPU cumulé local max :      %12.3f\n"
                 "  Temps CPU cumulé total :          %12.3f\n"),
               cpu_min, cpu_max, cpu_tot);

  }

#endif

  if (cs_glob_base_nbr == 1)
    bft_printf(_("  Temps CPU cumulé :                %12.3f\n"),
               this_info->cpu_tot);
}

/*----------------------------------------------------------------------------
 * Return pointer to linear system info.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems.
 *
 * parameters:
 *   name --> system name
 *   type --> resolution method
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
 *   solver_name --> solver name
 *   var_name    --> variable name
 *   convergence <-> Convergence info structure
 *   verbosity   --> Verbosity level
 *   n_iter_max  --> Maximum number of iterations
 *   precision   --> Precision limit
 *   r_norm      --> Residue normalization
 *   residue     --> Initial residue
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
    bft_printf(_("%s [%s]:\n"), solver_name, var_name);
    if (verbosity > 2)
      bft_printf(_("  n_iter     res_abs     res_nor\n"));
  }

}

/*----------------------------------------------------------------------------
 * Convergence test.
 *
 * parameters:
 *   solver_name --> solver name
 *   var_name    --> variable name
 *   n_iter      --> Number of iterations done
 *   residue     --> Non normalized residue
 *   convergence <-> Convergence information structure
 *
 * returns:
 *   1 if converged, 0 if not converged, -1 if not converged and maximum
 *   iteration number reached.
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
    bft_printf("   %5d %11.4e %11.4e\n",
               n_iter, residue, residue/convergence->r_norm);

  /* If not converged */

  if (residue > convergence->precision * convergence->r_norm) {

    if (n_iter < convergence->n_iterations_max) {
      int diverges = 0;
      if (residue > convergence->initial_residue * 10000.0 && residue > 100.)
        diverges = 1;
#if (_CS_STDC_VERSION >= 199901L)
      else if (isnan(residue) || isinf(residue))
        diverges = 1;
#endif
      if (diverges)
        bft_error(__FILE__, __LINE__, 0,
                  _("%s [%s]: divergence après %u itérations :\n"
                    "  résidu initial : %11.4e ; résidu courant : %11.4e"),
                  solver_name, var_name,
                  convergence->n_iterations,
                  convergence->initial_residue, residue);
      return 0;
    }
    else {
      if (verbosity > 0) {
        if (verbosity == 1) /* Already output if verbosity > 1 */
          bft_printf(_("%s [%s]:\n"), solver_name, var_name);
        if (verbosity <= 2) /* Already output if verbosity > 2 */
          bft_printf(_(final_fmt),
                     n_iter, residue, residue/convergence->r_norm);
        bft_printf(_(" @@ Attention : non convergence\n"));
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
 *   n_elts --> Local number of elements
 *   x      --> first vector in s = x.y
 *   y      --> second vector in s = x.y
 *
 * returns:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(cs_int_t          n_elts,
             const cs_real_t  *x,
             const cs_real_t  *y)
{
  double s = cblas_ddot(n_elts, x, 1, y, 1);

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_nbr > 1) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, cs_glob_base_mpi_comm);
    s = _sum;
  }

#endif /* defined(_CS_HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products, summing result over all ranks.
 *
 * parameters:
 *   n_elts --> Local number of elements
 *   x1     --> first vector in s1 = x1.y1
 *   y1     --> second vector in s1 = x1.y1
 *   x2     --> first vector in s2 = x2.y2
 *   y2     --> second vector in s2 = x2.y2
 *   s1     <-- result of s1 = x1.y1
 *   s2     <-- result of s2 = x2.y2
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_2(cs_int_t          n_elts,
                const cs_real_t  *x1,
                const cs_real_t  *y1,
                const cs_real_t  *x2,
                const cs_real_t  *y2,
                double           *s1,
                double           *s2)
{
  cs_int_t ii;
  double s[2];

  /* If a term appears in both dot products, we do not use the BLAS,
     as grouping both dot products in a same loop allows for better
     cache use and often better performance than separate dot products
     with even optimized BLAS */

  if (x1 == x2 || x1 == y2 || y1 == x2 || y1 == y2) {

    s[0] = 0.0; s[1] = 0.0;
    for (ii = 0; ii < n_elts; ii++) {
      s[0] += x1[ii] * y1[ii];
      s[1] += x2[ii] * y2[ii];
    }

  }
  else {

    s[0] = cblas_ddot(n_elts, x1, 1, y1, 1);
    s[1] = cblas_ddot(n_elts, x2, 1, y2, 1);

  }

#if defined(_CS_HAVE_MPI)

  if (cs_glob_base_nbr > 1) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM, cs_glob_base_mpi_comm);
    s[0] = _sum[0];
    s[1] = _sum[1];
  }

#endif /* defined(_CS_HAVE_MPI) */

  *s1 = s[0];
  *s2 = s[1];
}

/*----------------------------------------------------------------------------
 * Compute y <- x + alpha.y
 *
 * parameters:
 *   n     --> Number of elements in vectors x, y, z
 *   alpha --> Scalar alpha
 *   x     --> Vector x (size: n)
 *   y     <-> Vector y (size: n)
 *----------------------------------------------------------------------------*/

inline static void
_y_aypx(cs_int_t          n,
        cs_real_t         alpha,
        const cs_real_t  *x,
        cs_real_t        *restrict y)
{
#if defined(_CS_HAVE_ESSL)

  dzaxpy(n, alpha, y, 1, x, 1, y, 1);

#else

 {
   cs_int_t ii;

#if defined(__xlc__)
#pragma disjoint(alpha, *x, *y)
#endif

   for (ii = 0; ii < n; ii++)
     y[ii] = x[ii] + (alpha * y[ii]);
 }

#endif

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
 *   n_cells       -->  Local number of cells
 *   poly_degree   -->  preconditioning polynomial degree (0: diagonal)
 *   rotation_mode -->  Halo update option for rotational periodicity
 *   ad_inv        -->  Inverse of matrix diagonal
 *   ax            -->  Non-diagonal part of linear equation matrix
 *   rk            -->  Residue vector
 *   gk            <--  Result vector
 *   wk            ---  Working array
 *----------------------------------------------------------------------------*/

static void
_polynomial_preconditionning(cs_int_t            n_cells,
                             int                 poly_degree,
                             cs_perio_rota_t     rotation_mode,
                             const cs_real_t    *ad_inv,
                             const cs_matrix_t  *ax,
                             const cs_real_t    *rk,
                             cs_real_t          *restrict gk,
                             cs_real_t          *restrict wk)
{
  int deg_id;
  cs_int_t ii;

#if defined(__xlc__)
#pragma disjoint(*ad_inv, *ax, *rk, *gk, *wk)
#endif

  /* Polynomial of degree 0 (diagonal)
   *-----------------------------------*/

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

    cs_matrix_vector_multiply(rotation_mode, ax, gk, wk);

    for (ii = 0; ii < n_cells; ii++)
      gk[ii] = (rk[ii] - wk[ii]) * ad_inv[ii];

  }

}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using preconditioned conjugate gradient.
 *
 * Single-processor-optimized version.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      --> Variable name
 *   a             --> Matrix
 *   ax            --> Non-diagonal part of linear equation matrix
 *                     (only necessary if poly_degree > 0)
 *   n_cells_ext   --> Local number of cells + ghost cells sharing a face
 *   n_cells       --> Local number of cells
 *   n_faces       --> Local number of internal faces
 *   symmetric     --> 1: symmetric; 2: non-symmetric
 *   poly_degree   --> Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode --> Halo update option for rotational periodicity
 *   convergence   --> Convergence information structure
 *   rhs           --> Right hand side
 *   aux_size      --> Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *----------------------------------------------------------------------------*/

static void
_conjugate_gradient_sp(const char             *var_name,
                       const cs_matrix_t      *a,
                       const cs_matrix_t      *ax,
                       int                     poly_degree,
                       cs_perio_rota_t         rotation_mode,
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

  n_cols = cs_matrix_get_n_columns(a);
  n_rows = cs_matrix_get_n_rows(a);

  /* Allocate work arrays */

  {
    size_t  n_wa = 5;
    size_t  wa_size = n_cols;

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

  cs_matrix_get_diagonal(a, ad_inv);
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

  /* Initialize work arrays (not necessary, just for debugging) */

  for (ii = 0; ii < n_rows; ii++) {
    rk[ii] = 0.0;
    dk[ii] = 0.0;
    gk[ii] = 0.0;
    zk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);

  for (ii = 0; ii < n_rows; ii++) {
    rk[ii] = rk[ii] - rhs[ii];
    dk[ii] = rk[ii];
  }

  /* Polynomial preconditionning of order poly_degre */

  _polynomial_preconditionning(n_rows,
                               poly_degree,
                               rotation_mode,
                               ad_inv,
                               ax,
                               rk,
                               gk,
                               wk);

  /* Descent direction */
  /*-------------------*/

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

  rk_gkm1 = _dot_product(n_rows, rk, gk);

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

  /* Descent parameter */

  _dot_products_2(n_rows, rk, dk, dk, zk, &ro_0, &ro_1);

  alpha =  - ro_0 / ro_1;

  cblas_daxpy(n_rows, alpha, dk, 1, vx, 1);
  cblas_daxpy(n_rows, alpha, zk, 1, rk, 1);

  /* Convergence test */

  residue = sqrt(_dot_product(n_rows, rk, rk));

  cvg = _convergence_test(sles_name, var_name,
                          n_iter, residue, convergence);

  /* Current Iteration */
  /*-------------------*/

  while (cvg == 0) {

    n_iter += 1;

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 ax,
                                 rk,
                                 gk,
                                 wk);

    /* Descent parameter */

    rk_gk = _dot_product(n_rows, rk, gk);

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

    _y_aypx(n_rows, beta, gk, dk);  /* dk <- gk + (beta.dk) */

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    _dot_products_2(n_rows, rk, dk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

    cblas_daxpy(n_rows, alpha, dk, 1, vx, 1);
    cblas_daxpy(n_rows, alpha, zk, 1, rk, 1);

    /* Convergence test */

    residue = sqrt(_dot_product(n_rows, rk, rk));

    cvg = _convergence_test(sles_name, var_name,
                            n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using preconditioned conjugate gradient.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      --> Variable name
 *   a             --> Matrix
 *   ax            --> Non-diagonal part of linear equation matrix
 *                     (only necessary if poly_degree > 0)
 *   n_cells_ext   --> Local number of cells + ghost cells sharing a face
 *   n_cells       --> Local number of cells
 *   n_faces       --> Local number of internal faces
 *   symmetric     --> 1: symmetric; 2: non-symmetric
 *   poly_degree   --> Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode --> Halo update option for rotational periodicity
 *   convergence   --> Convergence information structure
 *   rhs           --> Right hand side
 *   aux_size      --> Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *----------------------------------------------------------------------------*/

static void
_conjugate_gradient_mp(const char             *var_name,
                       const cs_matrix_t      *a,
                       const cs_matrix_t      *ax,
                       int                     poly_degree,
                       cs_perio_rota_t         rotation_mode,
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

  n_cols = cs_matrix_get_n_columns(a);
  n_rows = cs_matrix_get_n_rows(a);

  /* Allocate work arrays */

  {
    size_t  n_wa = 5;
    size_t  wa_size = n_cols;

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

  cs_matrix_get_diagonal(a, ad_inv);
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

  /* Initialize work arrays (not necessary, just for debugging) */

  for (ii = 0; ii < n_rows; ii++) {
    rk[ii] = 0.0;
    dk[ii] = 0.0;
    gk[ii] = 0.0;
    zk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);

  for (ii = 0; ii < n_rows; ii++) {
    rk[ii] = rk[ii] - rhs[ii];
    dk[ii] = rk[ii];
  }

  /* Polynomial preconditionning of order poly_degre */

  _polynomial_preconditionning(n_rows,
                               poly_degree,
                               rotation_mode,
                               ad_inv,
                               ax,
                               rk,
                               gk,
                               wk);

  /* Descent direction */
  /*-------------------*/

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

  rk_gkm1 = _dot_product(n_rows, rk, gk);

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

  /* Descent parameter */

  _dot_products_2(n_rows, rk, dk, dk, zk, &ro_0, &ro_1);

  alpha =  - ro_0 / ro_1;

  cblas_daxpy(n_rows, alpha, dk, 1, vx, 1);
  cblas_daxpy(n_rows, alpha, zk, 1, rk, 1);

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
                                 ax,
                                 rk,
                                 gk,
                                 wk);

    /* compute residue and prepare descent parameter */

    _dot_products_2(n_rows, rk, rk, rk, gk, &residue, &rk_gk);

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

    _y_aypx(n_rows, beta, gk, dk);  /* dk <- gk + (beta.dk) */

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    _dot_products_2(n_rows, rk, dk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

    cblas_daxpy(n_rows, alpha, dk, 1, vx, 1);
    cblas_daxpy(n_rows, alpha, zk, 1, rk, 1);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using Jacobi.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   var_name      --> Variable name
 *   ad            --> Diagonal part of linear equation matrix
 *                     (only necessary if poly_degree > 0)
 *   ax            --> Non-diagonal part of linear equation matrix
 *                     (only necessary if poly_degree > 0)
 *   n_cells_ext   --> Local number of cells + ghost cells sharing a face
 *   n_cells       --> Local number of cells
 *   n_faces       --> Local number of internal faces
 *   symmetric     --> 1: symmetric; 2: non-symmetric
 *   rotation_mode --> Halo update option for rotational periodicity
 *   convergence   --> Convergence information structure
 *   rhs           --> Right hand side
 *   aux_size      --> Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *----------------------------------------------------------------------------*/

static void
_jacobi(const char             *var_name,
        const cs_real_t        *restrict ad,
        const cs_matrix_t      *ax,
        cs_perio_rota_t         rotation_mode,
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

  unsigned n_iter = 0;

  /* Tell IBM compiler not to alias */
#if defined(__xlc__)
#pragma disjoint(*rhs, *vx, *ad, *ad_inv)
#endif

  /* Preliminary calculations */
  /*--------------------------*/

  sles_name = _(cs_sles_type_name[CS_SLES_JACOBI]);

  n_cols = cs_matrix_get_n_columns(ax);
  n_rows = cs_matrix_get_n_rows(ax);

  /* Allocate work arrays */

  {
    size_t  wa_size = n_cols;

    if (aux_vectors == NULL || aux_size < (wa_size * 2))
      BFT_MALLOC(_aux_vectors, wa_size * 2, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    ad_inv = _aux_vectors;
    rk     = _aux_vectors + wa_size;
  }

  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad[ii];

  cvg = 0;

  /* Current iteration */
  /*-------------------*/

  while (cvg == 0) {

    n_iter += 1;

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    memcpy(vx, rhs, n_rows * sizeof(cs_real_t));  /* vx <- rhs */
    for (ii = n_rows; ii < n_cols; ii++)
      vx[ii] = 0.0;

    /* Compute Vx <- Vx - (A-diag).Rk */

    cs_matrix_alpha_a_x_p_beta_y(rotation_mode, -1.0, 1.0, ax, rk, vx);

    for (ii = 0; ii < n_rows; ii++)
      vx[ii] = vx[ii] * ad_inv[ii];

    /* Compute residue */

    res2 = 0.0;

    for (ii = 0; ii < n_rows; ii++) {
      register double r = ad[ii] * (vx[ii]-rk[ii]);
      res2 += (r*r);
    }

#if defined(_CS_HAVE_MPI)

    if (cs_glob_base_nbr > 1) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                    cs_glob_base_mpi_comm);
      res2 = _sum;
    }

#endif /* defined(_CS_HAVE_MPI) */

    residue = sqrt(res2);

    /* Convergence test */

    cvg = _convergence_test(sles_name, var_name,
                            n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);
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
 *   var_name      --> Variable name
 *   a             --> Matrix
 *   ax            --> Non-diagonal part of linear equation matrix
 *                     (only necessary if poly_degree > 0)
 *   n_cells_ext   --> Local number of cells + ghost cells sharing a face
 *   n_cells       --> Local number of cells
 *   n_faces       --> Local number of internal faces
 *   symmetric     --> 1: symmetric; 2: non-symmetric
 *   poly_degree   --> Preconditioning polynomial degree (0: diagonal)
 *   rotation_mode --> Halo update option for rotational periodicity
 *   convergence   --> Convergence information structure
 *   rhs           --> Right hand side
 *   aux_size      --> Number of elements in aux_vectors
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *----------------------------------------------------------------------------*/

static void
_bi_cgstab(const char             *var_name,
           const cs_matrix_t      *a,
           const cs_matrix_t      *ax,
           int                     poly_degree,
           cs_perio_rota_t         rotation_mode,
           cs_sles_convergence_t  *convergence,
           const cs_real_t        *rhs,
           cs_real_t              *restrict vx,
           size_t                  aux_size,
           void                   *aux_vectors)
{
  const char *sles_name;
  int cvg;
  cs_int_t  n_cols, n_rows, ii;
  double  _epzero = EPZERO*EPZERO; /* smaller than epzero */
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

  n_cols = cs_matrix_get_n_columns(a);
  n_rows = cs_matrix_get_n_rows(a);

  /* Allocate work arrays */

  {
    size_t  n_wa = 7;
    size_t  wa_size = n_cols;

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

  cs_matrix_get_diagonal(a, ad_inv);
  for (ii = 0; ii < n_rows; ii++)
    ad_inv[ii] = 1.0 / ad_inv[ii];

  /* Initialize work arrays (not necessary, just for debugging) */

  for (ii = 0; ii < n_rows; ii++) {
    res0[ii] = 0.0;
    rk[ii] = 0.0;
    pk[ii] = 0.0;
    zk[ii] = 0.0;
    uk[ii] = 0.0;
    vk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(rotation_mode, a, vx, res0);

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

    if (n_iter == 0)
      beta = _dot_product(n_rows, res0, rk);

    else {

      _dot_products_2(n_rows, res0, rk, rk, rk, &beta, &residue);
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
        bft_printf(_("  n_iter : %5d, res_abs : %11.4e, res_nor : %11.4e\n"),
                   n_iter, residue, residue/convergence->r_norm);
      else if (convergence->verbosity > 2)
        bft_printf("   %5d %11.4e %11.4e\n",
                   n_iter, residue, residue/convergence->r_norm);
      return;
    }

    if (CS_ABS(alpha) < _epzero) {
      bft_error(__FILE__, __LINE__, 0,
                _("%s [%s]:\n"
                  " @@ Attention : non convergence et arrêt\n\n"
                  "    Le coefficient alpha est inférieur à %12.4e\n\n"
                  "    La matrice ne peut plus être considérée come "
                  " inversible."),
                sles_name, var_name, alpha);
    }

    omega = beta*gamma / (alpha*betam1);
    betam1 = beta;

    /* Compute pk */

    for (ii = 0; ii < n_rows; ii++) {
      pk[ii] = rk[ii] + omega*(pk[ii] - alpha*uk[ii]);
    }

    /* Compute zk = c.pk */

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 ax,
                                 pk,
                                 zk,
                                 uk);

    /* Compute uk = A.zk */

    cs_matrix_vector_multiply(rotation_mode, a, zk, uk);

    /* Compute uk.res0 and gamma */

    ukres0 = _dot_product(n_rows, uk, res0);

    gamma = beta / ukres0;

    /* First update of vx and rk */

    cblas_daxpy(n_rows,  gamma, zk, 1, vx, 1);
    cblas_daxpy(n_rows, -gamma, uk, 1, rk, 1);

    /* Compute zk = C.rk (zk is overwritten, vk is a working array */

    _polynomial_preconditionning(n_rows,
                                 poly_degree,
                                 rotation_mode,
                                 ad_inv,
                                 ax,
                                 rk,
                                 zk,
                                 vk);

    /* Compute vk = A.zk and alpha */

    cs_matrix_vector_multiply(rotation_mode, a, zk, vk);

    _dot_products_2(n_rows, vk, rk, vk, vk, &ro_0, &ro_1);

    if (ro_1 < _epzero) {
      bft_error(__FILE__, __LINE__, 0,
                _("%s [%s]:\n"
                  " @@ Attention : non convergence et arrêt\n\n"
                  "    Le carre de la norme du vecteur de descente\n"
                  "    est inférieur à %12.4e\n\n"
                  "    La résolution ne progresse plus."),
                sles_name, var_name, _epzero);
    }

    alpha = ro_0 / ro_1;

    /* Final update of vx and rk */

    cblas_daxpy(n_rows,  alpha, zk, 1, vx, 1);
    cblas_daxpy(n_rows, -alpha, vk, 1, rk, 1);

    /* Convergence test at beginning of next iteration so
       as to group dot products for better parallel performance */
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);
}

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * General sparse linear system resolution
 *----------------------------------------------------------------------------*/

void CS_PROCF(reslin, RESLIN)
(
 const char       *cname,     /* --> variable name */
 const cs_int_t   *lname,     /* --> variable name length */
 const cs_int_t   *ncelet,    /* --> Number of cells, halo included */
 const cs_int_t   *ncel,      /* --> Number of local cells */
 const cs_int_t   *nfac,      /* --> Number of faces */
 const cs_int_t   *isym,      /* --> Symmetry indicator:
                                     1: symmetric; 2: not symmetric */
 const cs_int_t   *ireslp,    /* --> Resolution type:
                                     0: pcg; 1: Jacobi; 2: cg-stab */
 const cs_int_t   *ipol,      /* --> Preconditioning polynomial degree
                                     (0: diagonal) */
 const cs_int_t   *nitmap,    /* --> Number of max iterations */
 const cs_int_t   *iinvpe,    /* --> Indicator to cancel increments
                                     in rotational periodicty (2) or
                                     to exchange them as scalars (1) */
 const cs_int_t   *iwarnp,    /* --> Verbosity level */
 cs_int_t         *niterf,    /* <-- Number of iterations done */
 const cs_real_t  *epsilp,    /* --> Precision for iterative resolution */
 const cs_real_t  *rnorm,     /* --> Residue normalization */
 cs_real_t        *residu,    /* <-- Final non normalized residue */
 const cs_int_t   *ifacel,    /* --> Face -> cell connectivity  */
 const cs_real_t  *dam,       /* --> Matrix diagonal */
 const cs_real_t  *xam,       /* --> Matrix extra-diagonal terms */
 const cs_real_t  *rhs,       /* --> System right-hand side */
 cs_real_t        *vx         /* <-> System solution */
)
{
  char *var_name;
  cs_sles_type_t type;
  double  wt_start, wt_stop, cpu_start, cpu_stop;

  cs_sles_info_t *sles_info = NULL;
  cs_bool_t symmetric = (*isym == 1) ? true : false;
  cs_perio_rota_t rotation_mode = CS_PERIO_ROTA_COPY;

  assert(*ncelet >= *ncel);
  assert(*nfac > 0);
  assert(ifacel != NULL);

  wt_start =bft_timer_wtime();
  cpu_start =bft_timer_cpu_time();

  if (*iinvpe == 2)
    rotation_mode = CS_PERIO_ROTA_RESET;
  else if (*iinvpe == 3)
    rotation_mode = CS_PERIO_ROTA_IGNORE;

  var_name = cs_base_chaine_f_vers_c_cree(cname, *lname);

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
  default:
    type = CS_SLES_N_TYPES;
    assert(0);
  }

  sles_info = _find_or_add_system(var_name, type);

  /* Initialize number of iterations and residue, check for immediate return */

  *niterf = 0;
  *residu = sqrt(_dot_product(*ncel, rhs, rhs));

  if (*rnorm <= EPZERO || *residu <= EPZERO) {
    if (*iwarnp > 1)
      bft_printf(_("%s [%s]:\n"
                   "  sortie immédiate ; r_norm = %11.4e, residu = %11.4e\n"),
                 _(cs_sles_type_name[type]), var_name, *rnorm, *residu);
  }

  /* Solve sparse linear system */

  else {

    cs_sles_convergence_t  convergence;
    cs_matrix_t  *a = NULL;
    cs_matrix_t  *ax = NULL;

    if (type == CS_SLES_JACOBI) {

      ax = cs_glob_sles_base_matrix;
      cs_matrix_set_coefficients(ax, symmetric, NULL, xam);

    }
    else { /* if (type != CS_SLES_JACOBI) */

      a = cs_glob_sles_base_matrix;
      cs_matrix_set_coefficients(a, symmetric, dam, xam);

      if (*ipol > 0) {
        ax = cs_glob_sles_native_matrix;
        cs_matrix_set_coefficients(ax, symmetric, NULL, xam);
      }
    }

    _convergence_init(&convergence,
                      _(cs_sles_type_name[type]),
                      var_name,
                      *iwarnp,
                      *nitmap,
                      *epsilp,
                      *rnorm,
                      *residu);

    switch (type) {
    case CS_SLES_PCG:
      if (cs_glob_base_nbr == 1)
        _conjugate_gradient_sp(var_name,
                               a,
                               ax,
                               *ipol,
                               rotation_mode,
                               &convergence,
                               rhs,
                               vx,
                               0,
                               NULL);
      else
        _conjugate_gradient_mp(var_name,
                               a,
                               ax,
                               *ipol,
                               rotation_mode,
                               &convergence,
                               rhs,
                               vx,
                               0,
                               NULL);
      break;
    case CS_SLES_JACOBI:
      _jacobi(var_name,
              dam,
              ax,
              rotation_mode,
              &convergence,
              rhs,
              vx,
              0,
              NULL);
      break;
    case CS_SLES_BICGSTAB:
      _bi_cgstab(var_name,
                 a,
                 ax,
                 *ipol,
                 rotation_mode,
                 &convergence,
                 rhs,
                 vx,
                 0,
                 NULL);
      break;
    default:
      break;
    }

    if (a != NULL)
      cs_matrix_release_coefficients(a);

    if (ax != NULL)
      cs_matrix_release_coefficients(ax);

    *niterf = convergence.n_iterations;
    *residu = convergence.residue;

  }

  cs_base_chaine_f_vers_c_detruit(var_name);

  wt_stop =bft_timer_wtime();
  cpu_stop =bft_timer_cpu_time();

  if (sles_info->n_calls == 1)
    sles_info->n_iterations_min = *niterf;

  sles_info->n_calls += 1;

  if (sles_info->n_iterations_min > (unsigned)(*niterf))
    sles_info->n_iterations_min = *niterf;
  if (sles_info->n_iterations_max < (unsigned)(*niterf))
    sles_info->n_iterations_max = *niterf;

  sles_info->n_iterations_last = *niterf;
  sles_info->n_iterations_tot += *niterf;

  sles_info->wt_tot += (wt_stop - wt_start);
  sles_info->cpu_tot += (cpu_stop - cpu_start);
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
  cs_bool_t  periodic = CS_FALSE;

  assert(mesh != NULL);

  if (mesh->n_init_perio > 0)
    periodic = CS_TRUE;

  cs_glob_sles_base_matrix = cs_matrix_create(CS_MATRIX_NATIVE,
                                              CS_FALSE,
                                              CS_TRUE,
                                              periodic,
                                              mesh->n_cells,
                                              mesh->n_cells_with_ghosts,
                                              mesh->n_i_faces,
                                              mesh->global_cell_num,
                                              mesh->i_face_cells);

  cs_glob_sles_native_matrix = cs_matrix_create(CS_MATRIX_NATIVE,
                                                CS_FALSE,
                                                CS_TRUE,
                                                periodic,
                                                mesh->n_cells,
                                                mesh->n_cells_with_ghosts,
                                                mesh->n_i_faces,
                                                mesh->global_cell_num,
                                                mesh->i_face_cells);
}

/*----------------------------------------------------------------------------
 * Finalize sparse linear equation solver API.
 *----------------------------------------------------------------------------*/

void
cs_sles_finalize(void)
{
  int ii;

  /* Free sytem info */

  for (ii = 0; ii < cs_glob_sles_n_systems; ii++) {
    _sles_info_dump(cs_glob_sles_systems[ii]);
    _sles_info_destroy(&(cs_glob_sles_systems[ii]));
  }

  BFT_FREE(cs_glob_sles_systems);

  cs_glob_sles_n_systems = 0;
  cs_glob_sles_n_max_systems = 0;

  /* Free matrix structures */

  cs_matrix_destroy(&cs_glob_sles_native_matrix);
  cs_matrix_destroy(&cs_glob_sles_base_matrix);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
