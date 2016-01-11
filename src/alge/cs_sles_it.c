/*============================================================================
 * Sparse Linear Equation Solvers
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_blas.h"
#include "cs_log.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_sles_it.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_it.c
        Iterative linear solvers

  \enum cs_sles_it_type_t

  \brief Iterative solver types

  \var CS_SLES_PCG
       Preconditioned conjugate gradient
  \var CS_SLES_JACOBI
       Jacobi
  \var CS_SLES_BICGSTAB
       Preconditioned BiCGstab (biconjugate gradient stabilized)
  \var CS_SLES_BICGSTAB2
       Preconditioned BiCGstab2 (biconjugate gradient stabilized)
  \var CS_SLES_GMRES
       Preconditioned GMRES (generalized minimum residual)
  \var CS_SLES_B_GAUSS_SEIDEL
       Block Gauss-Seidel

 \page sles_it Iterative linear solvers.

 For Krylov space solvers (all here except Jacobi), reconditioning is based
 on a Neumann polynomial of degree \a poly_degree, with a negative value
 meaning no preconditionig, and 0 diagonal preconditioning.

 For positive values of \a poly_degree, the preconditioning is explained here:
 \a D being the diagonal part of matrix \a A and \a X its extra-diagonal
 part, it can be written \f$A=D(Id+D^{-1}X)\f$. Therefore
 \f$A^{-1}=(Id+D^{-1}X)^{-1}D^{-1}\f$. A series development of
 \f$Id+D^{-1}X\f$ can then be used which yields, symbolically,
 \f[
 Id+\sum\limits_{I=1}^{poly\_degree}\left(-D^{-1}X\right)^{I}
 \f]

 The efficiency of the polynomial preconditioning will vary depending
 on the system type. In most cases, diagonal or degree 1 provide
 best results. Each polynomial precoditioning degree above 0 adds one
 matrix-vector product per inital matrix-vector product of the algorithm.
 Switching from diagonal to polynomial degree 1 often divides the number of
 required iterations by approximately 2, but each iteration then costs
 close to 2 times that of diagonal preconditoning (other vector operations
 are not doubled), so the net gain is often about 10%. Higher degree
 polynomials usually lead to diminishing returns.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define RINFIN  1.E+30

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures
 * through at least 2012) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Solver setup data */
/*-------------------*/

typedef struct _cs_sles_it_setup_t {

  bool                 single_reduce;    /* single reduction mode for PCG */

  double               initial_residue;  /* Last initial residue value */

  cs_lnum_t            n_rows;           /* Number of associated rows */

  const cs_real_t     *ad_inv;           /* pointer to diagonal inverse */
  cs_real_t           *_ad_inv;          /* private pointer to
                                            diagonal inverse */

} cs_sles_it_setup_t;

/* Solver additional data */
/*------------------------*/

typedef struct _cs_sles_it_add_t {

  cs_lnum_t           *order;            /* ordering */

} cs_sles_it_add_t;

/* Basic per linear system options and logging */
/*---------------------------------------------*/

struct _cs_sles_it_t {

  /* Base settings */

  cs_sles_it_type_t    type;               /* Solver type */

  int                  poly_degree;        /* preconditioning polynomial
                                              degree (0: diagonal,
                                              -1: non-preconditioned) */

  bool                 update_stats;       /* do stats need to be updated ? */

  int                  n_max_iter;         /* maximum number of iterations */

  /* Performance data */

  unsigned             n_setups;           /* Number of times system setup */
  unsigned             n_solves;           /* Number of times system solved */

  unsigned             n_iterations_last;  /* Number of iterations for last
                                              system resolution */
  unsigned             n_iterations_min;   /* Minimum number ot iterations
                                              in system resolution history */
  unsigned             n_iterations_max;   /* Maximum number ot iterations
                                              in system resolution history */
  unsigned long long   n_iterations_tot;   /* Total accumulated number of
                                              iterations */

  cs_timer_counter_t   t_setup;            /* Total setup */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /* Communicator used for reduction operations
     (if left at NULL, main communicator will be used) */

# if defined(HAVE_MPI)
  MPI_Comm comm;
# endif

  /* Solver setup */

  const struct _cs_sles_it_t  *shared;     /* pointer to context sharing some
                                              setup data, or NULL */

  cs_sles_it_add_t            *add_data;   /* additional data */

  cs_sles_it_setup_t          *setup_data; /* setup data */

};

/* Convergence testing and tracking */
/*----------------------------------*/

typedef struct _cs_sles_it_convergence_t {

  const char          *name;               /* Pointer to name string */

  int                  verbosity;          /* Verbosity level */

  unsigned             n_iterations;       /* Current number of iterations */
  unsigned             n_iterations_max;   /* Maximum number of iterations */

  double               precision;          /* Precision limit */
  double               r_norm;             /* Residue normalization */
  double               residue;            /* Current residue */

} cs_sles_it_convergence_t;

/*============================================================================
 *  Global variables
 *=================================================================x===========*/

/* Mean system rows threshold under which we use single-reduce version of PCG */

static cs_lnum_t _pcg_sr_threshold = 512;

/* Sparse linear equation solver type names */

const char *cs_sles_it_type_name[] = {N_("Conjugate gradient"),
                                      N_("Jacobi"),
                                      N_("BiCGstab"),
                                      N_("BiCGstab2"),
                                      N_("GMRES"),
                                      N_("Block Gauss-Seidel")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize or reset convergence info structure.
 *
 * At this stage, the initial residue is set to HUGE_VAL, as it is unknown.
 *
 * parameters:
 *   solver_name <-- solver name
 *   var_name    <-- variable name
 *   convergence <-> Convergence info structure
 *   verbosity   <-- Verbosity level
 *   n_iter_max  <-- Maximum number of iterations
 *   precision   <-- Precision limit
 *   r_norm      <-- Residue normalization
 *   residue     <-> Initial residue
 *----------------------------------------------------------------------------*/

static void
_convergence_init(cs_sles_it_convergence_t  *convergence,
                  const char                *solver_name,
                  const char                *var_name,
                  int                        verbosity,
                  unsigned                   n_iter_max,
                  double                     precision,
                  double                     r_norm,
                  double                    *residue)
{
  *residue = HUGE_VAL;  /* Unknown at this stage */

  convergence->name = solver_name;
  convergence->verbosity = verbosity;

  convergence->n_iterations = 0;
  convergence->n_iterations_max = n_iter_max;

  convergence->precision = precision;
  convergence->r_norm = r_norm;
  convergence->residue = *residue;

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
 *   c           <-- pointer to solver context info
 *   n_iter      <-- Number of iterations done
 *   residue     <-- Non normalized residue
 *   convergence <-> Convergence information structure
 *
 * returns:
 *   convergence status.
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_convergence_test(cs_sles_it_t              *c,
                  unsigned                   n_iter,
                  double                     residue,
                  cs_sles_it_convergence_t  *convergence)
{
  const int verbosity = convergence->verbosity;
  const cs_sles_it_setup_t  *s = c->setup_data;

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
      if (residue > s->initial_residue * 10000.0 && residue > 100.)
        diverges = 1;
#if (__STDC_VERSION__ >= 199901L)
      else if (isnan(residue) || isinf(residue))
        diverges = 1;
#endif
      if (diverges) {
        bft_printf(_("\n\n"
                     "%s [%s]: divergence after %u iterations:\n"
                     "  initial residual: %11.4e; current residual: %11.4e\n"),
                   cs_sles_it_type_name[c->type], convergence->name,
                   convergence->n_iterations,
                   s->initial_residue, convergence->residue);
        return CS_SLES_DIVERGED;
      }
      else
        return 0;
    }
    else {
      if (verbosity > 0) {
        if (verbosity == 1) /* Already output if verbosity > 1 */
          bft_printf("%s [%s]:\n", cs_sles_it_type_name[c->type],
                     convergence->name);
        if (verbosity <= 2) /* Already output if verbosity > 2 */
          bft_printf(_(final_fmt),
                     n_iter, residue, residue/convergence->r_norm);
        bft_printf(_(" @@ Warning: non convergence\n"));
      }
      return CS_SLES_MAX_ITERATION;
    }

  }

  /* If converged */

  else {
    if (verbosity == 2) /* Already output if verbosity > 2 */
      bft_printf(final_fmt, n_iter, residue, residue/convergence->r_norm);
    return CS_SLES_CONVERGED;
  }

}

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- first vector in s = x.y
 *   y      <-- second vector in s = x.y
 *
 * returns:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

inline static double
_dot_product(const cs_sles_it_t  *c,
             const cs_real_t     *x,
             const cs_real_t     *y)
{
  double s = cs_dot(c->setup_data->n_rows, x, y);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM, c->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute dot product x.x, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s = x.x
 *
 * returns:
 *   result of s = x.x
 *----------------------------------------------------------------------------*/

inline static double
_dot_product_xx(const cs_sles_it_t  *c,
                const cs_real_t     *x)
{
  double s;

  s = cs_dot_xx(c->setup_data->n_rows, x);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum;
    MPI_Allreduce(&s, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                  c->comm);
    s = _sum;
  }

#endif /* defined(HAVE_MPI) */

  return s;
}

/*----------------------------------------------------------------------------
 * Compute 2 dot products x.x and x.y, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s1 = x.x and s2 = x.y
 *   y      <-- vector in s2 = x.y
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy(const cs_sles_it_t  *c,
                    const cs_real_t     *x,
                    const cs_real_t     *y,
                    double              *s1,
                    double              *s2)
{
  double s[2];

  cs_dot_xx_xy(c->setup_data->n_rows, x, y, s, s+1);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM,
                  c->comm);
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
 *   c      <-- pointer to solver context info
 *   x      <-- vector in s1 = x.y
 *   y      <-- vector in s1 = x.y and s2 = y.z
 *   z      <-- vector in s2 = y.z
 *   s1     --> result of s1 = x.y
 *   s2     --> result of s2 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xy_yz(const cs_sles_it_t  *c,
                    const cs_real_t     *x,
                    const cs_real_t     *y,
                    const cs_real_t     *z,
                    double              *s1,
                    double              *s2)
{
  double s[2];

  cs_dot_xy_yz(c->setup_data->n_rows, x, y, z, s, s+1);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[2];
    MPI_Allreduce(s, _sum, 2, MPI_DOUBLE, MPI_SUM,
                  c->comm);
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
 *   c      <-- pointer to solver context info
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   s1     --> result of s1 = x.x
 *   s2     --> result of s2 = x.y
 *   s3     --> result of s3 = y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_xy_yz(const cs_sles_it_t  *c,
                       const cs_real_t     *x,
                       const cs_real_t     *y,
                       const cs_real_t     *z,
                       double              *s1,
                       double              *s2,
                       double              *s3)
{
  double s[3];

  cs_dot_xx_xy_yz(c->setup_data->n_rows, x, y, z, s, s+1, s+2);

#if defined(HAVE_MPI)

  if (c->comm != MPI_COMM_NULL) {
    double _sum[3];

    MPI_Allreduce(s, _sum, 3, MPI_DOUBLE, MPI_SUM, c->comm);
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
 * Compute 5 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   xx     --> result of x.x
 *   yy     --> result of y.y
 *   xy     --> result of x.y
 *   xz     --> result of x.z
 *   yz     --> result of y.z
 *----------------------------------------------------------------------------*/

inline static void
_dot_products_xx_yy_xy_xz_yz(const cs_sles_it_t  *c,
                             const cs_real_t     *x,
                             const cs_real_t     *y,
                             const cs_real_t     *z,
                             double              *xx,
                             double              *yy,
                             double              *xy,
                             double              *xz,
                             double              *yz)
{
  double s[5];

  cs_dot_xx_yy_xy_xz_yz(c->setup_data->n_rows, x, y, z, s, s+1, s+2, s+3, s+4);

  cs_parall_sum(5, CS_DOUBLE, s);

  *xx = s[0];
  *yy = s[1];
  *xy = s[2];
  *xz = s[3];
  *yz = s[4];
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
 *   c             <-- pointer to solver context info
 *   rotation_mode <-- halo update option for rotational periodicity
 *   a             <-- linear equation matrix
 *   rk            <-- residue vector
 *   gk            --> result vector
 *   wk            --- Working array
 *----------------------------------------------------------------------------*/

static void
_polynomial_preconditionning(const cs_sles_it_t  *c,
                             cs_halo_rotation_t   rotation_mode,
                             const cs_matrix_t   *a,
                             const cs_real_t     *rk,
                             cs_real_t           *restrict gk,
                             cs_real_t           *restrict wk)
{
  int deg_id;
  cs_lnum_t ii;

  const cs_sles_it_setup_t *s = c->setup_data;
  const cs_lnum_t n_rows = s->n_rows;
  const cs_real_t *restrict ad_inv = s->ad_inv;

  if (c->poly_degree < 0) {
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      gk[ii] = rk[ii];
    return;
  }

  /* Polynomial of degree 0 (diagonal)
   *-----------------------------------*/

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    gk[ii] = rk[ii] * ad_inv[ii];

  /* Polynomial of degree n
   *-----------------------
   *
   *                  n=1                    n=2
   * gk = ((1/ad) - (1/ad).ax.(1/ad) + (1/ad).ax.(1/ad).ax.(1/ad) + ... ).rk
   */

  for (deg_id = 1; deg_id <= c->poly_degree; deg_id++) {

    /* Compute Wk = (A-diag).Gk */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, gk, wk);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      gk[ii] = (rk[ii] - wk[ii]) * ad_inv[ii];

  }
}

/*----------------------------------------------------------------------------
* Compute inverses of dense 3*3 matrices.
*
* parameters:
*   n_blocks <-- number of blocks
*   ad       <-- diagonal part of linear equation matrix
*   ad_inv   --> inverse of the diagonal part of linear equation matrix
*----------------------------------------------------------------------------*/

static void
_fact_lu33(cs_lnum_t         n_blocks,
           const cs_real_t  *ad,
           cs_real_t        *ad_inv)
{
# pragma omp parallel for if(n_blocks > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_blocks; i++) {

    cs_real_t *restrict _ad_inv = &ad_inv[9*i];
    const cs_real_t *restrict  _ad = &ad[9*i];

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
 * Setup context for iterative linear solver.
 *
 * This function is common to most solvers
 *
 * parameters:
 *   c                <-> pointer to solver context info
 *   a                <-- matrix
 *   diag_block_size  <-- diagonal block size
 *   block_33_inverse <-- if diagonal block is of size 3, compute inverse of
 *                        block if true, inverse of block diagonal otherwise
 *----------------------------------------------------------------------------*/

static void
_setup_sles_it(cs_sles_it_t       *c,
               const cs_matrix_t  *a,
               int                 diag_block_size,
               bool                block_33_inverse)
{
  cs_timer_t t0;
  if (c->update_stats == true)
    t0 = cs_timer_time();

  cs_sles_it_setup_t *sd = c->setup_data;

  if (sd == NULL) {
    BFT_MALLOC(c->setup_data, 1, cs_sles_it_setup_t);
    sd = c->setup_data;
    sd->ad_inv = NULL;
    sd->_ad_inv = NULL;
  }

  sd->n_rows = cs_matrix_get_n_rows(a) * diag_block_size;

  sd->single_reduce = false;
  sd->initial_residue = -1;

  /* Setup diagonal inverse */

  if (c->poly_degree >= 0) {

    const cs_sles_it_t  *s = c->shared;
    if (s != NULL) {
      if (s->setup_data == NULL)
        s = NULL;
      else if (s->setup_data->ad_inv == NULL)
        s = NULL;
    }

    if (s != NULL) {
      sd->ad_inv = s->setup_data->ad_inv;
      BFT_FREE(sd->_ad_inv);
    }
    else {

      if (diag_block_size != 3 || block_33_inverse == false) {

        const cs_lnum_t n_rows = sd->n_rows;

        BFT_REALLOC(sd->_ad_inv, sd->n_rows, cs_real_t);
        sd->ad_inv = sd->_ad_inv;

        cs_matrix_copy_diagonal(a, sd->_ad_inv);

#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_rows; i++)
          sd->_ad_inv[i] = 1.0 / sd->_ad_inv[i];

      }

      else {

        BFT_REALLOC(sd->_ad_inv, sd->n_rows*diag_block_size, cs_real_t);
        sd->ad_inv = sd->_ad_inv;

        const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);
        const cs_lnum_t  n_blocks = sd->n_rows / diag_block_size;

        _fact_lu33(n_blocks, ad, sd->_ad_inv);

      }

    }

  }

  /* Check for single-reduction */

#if defined(HAVE_MPI)

  if (c->type == CS_SLES_PCG) {
    cs_gnum_t n_m_rows = sd->n_rows;
    cs_parall_sum(1, CS_GNUM_TYPE, &n_m_rows);
    if (c->comm != MPI_COMM_NULL && c->comm != cs_glob_mpi_comm) {
      int size;
      MPI_Comm_size(c->comm, &size);
      n_m_rows /= (cs_gnum_t)size;
    }
    else
      n_m_rows /= (cs_gnum_t)cs_glob_n_ranks;
    if (cs_glob_n_ranks > 1 && c->comm != cs_glob_mpi_comm)
      MPI_Bcast(&n_m_rows, 1, CS_MPI_GNUM, 0, cs_glob_mpi_comm);

    if (n_m_rows < (cs_gnum_t)_pcg_sr_threshold)
      sd->single_reduce = true;
  }

#endif

  /* Now finish */

  if (c->update_stats == true) {
    cs_timer_t t1 = cs_timer_time();
    c->n_setups += 1;
    cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
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
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    int                        diag_block_size,
                    cs_halo_rotation_t         rotation_mode,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *restrict vx,
                    size_t                     aux_size,
                    void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii, jj;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, beta, residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 4;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    gk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Polynomial preconditionning of order poly_degre */

  _polynomial_preconditionning(c,
                               rotation_mode,
                               a,
                               rk,
                               gk,
                               dk); /* dk used as work array here */

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  _dot_products_xx_xy(c, rk, gk, &residue, &rk_gkm1);
  residue = sqrt(residue);

  /* If no solving required, finish here */

  c->setup_data->initial_residue = residue;
  cvg = _convergence_test(c, n_iter, residue, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

    /* Convergence test */

    residue = sqrt(_dot_product_xx(c, rk));

    cvg = _convergence_test(c, n_iter, residue, convergence);

    /* Current Iteration */
    /*-------------------*/

  }

  while (cvg == CS_SLES_ITERATING) {

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 rk,
                                 gk,
                                 zk); /* zk used as work array here */

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy(c, rk, gk, &residue, &rk_gk);

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residue, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = gk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned conjugate gradient
 * with single reduction.
 *
 * For more information, see Lapack Working note 56, at
 * http://www.netlib.org/lapack/lawnspdf/lawn56.pdf)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_sr(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       int                        diag_block_size,
                       cs_halo_rotation_t         rotation_mode,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *restrict vx,
                       size_t                     aux_size,
                       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii, jj;
  double  ro_0, ro_1, alpha, rk_gkm1, rk_gk, gk_sk, beta, residue;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict gk, *restrict sk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 5;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    gk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
    sk = _aux_vectors + wa_size*4;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Polynomial preconditionning of order poly_degre */
  /* gk = c_1 * rk  (zk = c_1 * rk) */

  _polynomial_preconditionning(c,
                               rotation_mode,
                               a,
                               rk,
                               gk,
                               dk); /* dk used as work array here */

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = gk[ii];

#else

  memcpy(dk, gk, n_rows * sizeof(cs_real_t));

#endif

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk); /* zk = A.dk */

  /* Descent parameter */

  _dot_products_xx_xy_yz(c, rk, dk, zk, &residue, &ro_0, &ro_1);
  residue = sqrt(residue);

  c->setup_data->initial_residue = residue;

  /* If no solving required, finish here */

  cvg = _convergence_test(c, n_iter, residue, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    alpha =  - ro_0 / ro_1;

    rk_gkm1 = ro_0;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

    /* Convergence test */

    residue = sqrt(_dot_product_xx(c, rk));

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 rk,
                                 gk,
                                 sk); /* sk used as work array here */

    cs_matrix_vector_multiply(rotation_mode, a, gk, sk);  /* sk = A.zk */

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy_yz(c, rk, gk, sk, &residue, &rk_gk, &gk_sk);

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residue, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_gk / rk_gkm1;
    rk_gkm1 = rk_gk;

    ro_1 = gk_sk - beta*beta*ro_1;
    ro_0 = rk_gk;

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        dk[ii] = gk[ii] + (beta * dk[ii]);
        vx[ii] += alpha * dk[ii];
      }
#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++) {
        zk[jj] = sk[jj] + (beta * zk[jj]);
        rk[jj] += alpha * zk[jj];
      }
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using non-preconditioned conjugate gradient.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_npc(cs_sles_it_t              *c,
                        const cs_matrix_t         *a,
                        int                        diag_block_size,
                        cs_halo_rotation_t         rotation_mode,
                        cs_sles_it_convergence_t  *convergence,
                        const cs_real_t           *rhs,
                        cs_real_t                 *restrict vx,
                        size_t                     aux_size,
                        void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii, jj;
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, beta, residue;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict zk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 3;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    zk = _aux_vectors + wa_size*2;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] -= rhs[ii];

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = rk[ii];

#else

  memcpy(dk, rk, n_rows * sizeof(cs_real_t));

#endif

  rk_rkm1 = _dot_product_xx(c, rk);
  residue = sqrt(rk_rkm1);

  /* If no solving required, finish here */

  c->setup_data->initial_residue = residue;
  cvg = _convergence_test(c, n_iter, residue, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    /* Descent parameter */

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

    /* Convergence test */

    residue = sqrt(_dot_product(c, rk, rk));

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy(c, rk, rk, &residue, &rk_rk);

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residue, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_rk / rk_rkm1;
    rk_rkm1 = rk_rk;

#   pragma omp parallel for firstprivate(alpha) if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = rk[ii] + (beta * dk[ii]);

    cs_matrix_vector_multiply(rotation_mode, a, dk, zk);

    _dot_products_xy_yz(c, rk, dk, zk, &ro_0, &ro_1);

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using non-preconditioned conjugate gradient
 * with single reduction.
 *
 * For more information, see Lapack Working note 56, at
 * http://www.netlib.org/lapack/lawnspdf/lawn56.pdf)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_conjugate_gradient_npc_sr(cs_sles_it_t              *c,
                           const cs_matrix_t         *a,
                           int                        diag_block_size,
                           cs_halo_rotation_t         rotation_mode,
                           cs_sles_it_convergence_t  *convergence,
                           const cs_real_t           *rhs,
                           cs_real_t                 *restrict vx,
                           size_t                     aux_size,
                           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii, jj;
  double  ro_0, ro_1, alpha, rk_rkm1, rk_rk, rk_sk, beta, residue;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict dk, *restrict sk;
  cs_real_t *restrict zk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 4;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
    dk = _aux_vectors + wa_size;
    zk = _aux_vectors + wa_size*2;
    sk = _aux_vectors + wa_size*3;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  /* Residue and descent direction */

  cs_matrix_vector_multiply(rotation_mode, a, vx, rk);  /* rk = A.x0 */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    rk[ii] = rk[ii] - rhs[ii];

  /* Descent direction */
  /*-------------------*/

#if defined(HAVE_OPENMP)

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++)
    dk[ii] = rk[ii];

#else

  memcpy(dk, rk, n_rows * sizeof(cs_real_t));

#endif

  cs_matrix_vector_multiply(rotation_mode, a, dk, zk); /* zk = A.dk */

  /* Descent parameter */

  _dot_products_xx_xy_yz(c, rk, dk, zk, &residue, &ro_0, &ro_1);
  residue = sqrt(residue);

  /* If no solving required, finish here */

  c->setup_data->initial_residue = residue;

  cvg = _convergence_test(c, n_iter, residue, convergence);

  if (cvg == CS_SLES_ITERATING) {

    n_iter = 1;

    alpha =  - ro_0 / ro_1;

    rk_rkm1 = ro_0;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * dk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] += (alpha * zk[jj]);
    }

    /* Convergence test */

    residue = sqrt(_dot_product_xx(c, rk));

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    cs_matrix_vector_multiply(rotation_mode, a, rk, sk);  /* sk = A.zk */

    /* compute residue and prepare descent parameter */

    _dot_products_xx_xy(c, rk, sk, &residue, &rk_sk);

    rk_rk = residue;

    residue = sqrt(residue);

    /* Convergence test for end of previous iteration */

    if (n_iter > 1)
      cvg = _convergence_test(c, n_iter, residue, convergence);

    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    /* Complete descent parameter computation and matrix.vector product */

    beta = rk_rk / rk_rkm1;
    rk_rkm1 = rk_rk;

    ro_1 = rk_sk - beta*beta*ro_1;
    ro_0 = rk_rk;

    alpha =  - ro_0 / ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++) {
        dk[ii] = rk[ii] + (beta * dk[ii]);
        vx[ii] += alpha * dk[ii];
      }
#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++) {
        zk[jj] = sk[jj] + (beta * zk[jj]);
        rk[jj] += alpha * zk[jj];
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
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_jacobi(cs_sles_it_t              *c,
        const cs_matrix_t         *a,
        int                        diag_block_size,
        cs_halo_rotation_t         rotation_mode,
        cs_sles_it_convergence_t  *convergence,
        const cs_real_t           *rhs,
        cs_real_t                 *restrict vx,
        size_t                     aux_size,
        void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii;
  double  res2, residue;
  cs_real_t *_aux_vectors;
  cs_real_t *restrict rk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 1;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk = _aux_vectors;
  }

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    register double r;

    n_iter += 1;

#if defined(HAVE_OPENMP)

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      rk[ii] = vx[ii];

#else

    memcpy(rk, vx, n_rows * sizeof(cs_real_t));   /* rk <- vx */

#endif

    /* Compute Vx <- Vx - (A-diag).Rk and residue. */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, rk, vx);

    res2 = 0.0;

#   pragma omp parallel for private(r) reduction(+:res2) if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      vx[ii] = (rhs[ii]-vx[ii])*ad_inv[ii];
      r = ad[ii] * (vx[ii]-rk[ii]);
      res2 += (r*r);
    }

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                    c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residue = sqrt(res2); /* Actually, residue of previous iteration */

    /* Convergence test */

    if (n_iter == 1)
      c->setup_data->initial_residue = residue;

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Block Jacobi utilities.
 * Compute forward and backward to solve an LU 3*3 system.
 *
 * parameters:
 *   mat   <-- 3*3*dim matrix
 *   x     --> solution
 *   b     --> 1st part of RHS (c - b)
 *   c     --> 2nd part of RHS (c - b)
 *----------------------------------------------------------------------------*/

inline static void
_fw_and_bw_lu33(const cs_real_t  mat[],
                cs_real_t        x[],
                const cs_real_t  b[],
                const cs_real_t  c[])
{
  cs_real_t  aux[3];

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
 *   c             <-- pointer to solver context info
 *   a             <-- linear equation matrix
 *   rotation_mode <-- halo update option for rotational periodicity
 *   convergence   <-- convergence information structure
 *   rhs           <-- right hand side
 *   vx            --> system solution
 *   aux_size      <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors   --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_block_3_jacobi(cs_sles_it_t              *c,
                const cs_matrix_t         *a,
                cs_halo_rotation_t         rotation_mode,
                cs_sles_it_convergence_t  *convergence,
                const cs_real_t           *rhs,
                cs_real_t                 *restrict vx,
                size_t                     aux_size,
                void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  res2, residue;
  cs_real_t *_aux_vectors;
  cs_real_t  *restrict rk, *restrict vxx;

  const int diag_block_size = 3;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, true);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_lnum_t n_rows = c->setup_data->n_rows;
  const cs_lnum_t n_blocks = c->setup_data->n_rows / 3;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 2;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    rk  = _aux_vectors;
    vxx = _aux_vectors + wa_size;
  }

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    n_iter += 1;
    memcpy(rk, vx, n_rows * sizeof(cs_real_t));  /* rk <- vx */

    /* Compute Vx <- Vx - (A-diag).Rk and residue. */

    cs_matrix_exdiag_vector_multiply(rotation_mode, a, rk, vxx);

    res2 = 0.0;

#   pragma omp parallel for reduction(+:res2) if(n_blocks > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_blocks; ii++) {
      _fw_and_bw_lu33(ad_inv + 9*ii,
                      vx + 3*ii,
                      vxx + 3*ii,
                      rhs + 3*ii);
      for (cs_lnum_t jj = 0; jj < 3; jj++) {
        register double r = 0.0;
        for (cs_lnum_t kk = 0; kk < 3; kk++)
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

    residue = sqrt(res2); /* Actually, residue of previous iteration */

    if (n_iter == 1)
      c->setup_data->initial_residue = residue;

    /* Convergence test */

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Test for (and eventually report) breakdown.
 *
 * parameters:
 *   c           <-- pointer to solver context info
 *   convergence <-- convergence information structure
 *   coeff       <-- coefficient name
 *   coeff       <-- coefficient to test
 *   epsilon     <-- value to test against
 *   n_iter      <-- current number of iterations
 *   cvg         <-> convergence status
 *
 * returns:
 *   true in case of breakdown, false otherwise
 *----------------------------------------------------------------------------*/

static inline bool
_breakdown(cs_sles_it_t                 *c,
           cs_sles_it_convergence_t     *convergence,
           const char                   *coeff_name,
           double                        coeff,
           double                        epsilon,
           double                        residue,
           int                           n_iter,
           cs_sles_convergence_state_t  *cvg)
{
  bool retval = false;

  if (CS_ABS(coeff) < epsilon) {

    bft_printf
      (_("\n\n"
         "%s [%s]:\n"
         " @@ Warning: non convergence\n"
         "\n"
         "    norm of coefficient \"%s\" is lower than %12.4e\n"
         "\n"
         "    The resolution does not progress anymore."),
       cs_sles_it_type_name[c->type], convergence->name, coeff_name, epsilon);
    bft_printf(_("  n_iter : %5u, res_abs : %11.4e, res_nor : %11.4e\n"),
               n_iter, residue, residue/convergence->r_norm);

    *cvg = CS_SLES_BREAKDOWN;
    retval = true;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using preconditioned Bi-CGSTAB.
 *
 * Parallel-optimized version, groups dot products, at the cost of
 * computation of the preconditionning for n+1 iterations instead of n.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_bi_cgstab(cs_sles_it_t              *c,
           const cs_matrix_t         *a,
           int                        diag_block_size,
           cs_halo_rotation_t         rotation_mode,
           cs_sles_it_convergence_t  *convergence,
           const cs_real_t           *rhs,
           cs_real_t                 *restrict vx,
           size_t                     aux_size,
           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t  ii, jj;
  double  _epzero = 1.e-30; /* smaller than epzero */
  double  ro_0, ro_1, alpha, beta, betam1, gamma, omega, ukres0;
  double  residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict res0, *restrict rk, *restrict pk, *restrict zk;
  cs_real_t  *restrict uk, *restrict vk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    const size_t n_wa = 6;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    res0 = _aux_vectors;
    rk = _aux_vectors + wa_size;
    pk = _aux_vectors + wa_size*2;
    zk = _aux_vectors + wa_size*3;
    uk = _aux_vectors + wa_size*4;
    vk = _aux_vectors + wa_size*5;
  }

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    pk[ii] = 0.0;
    uk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(rotation_mode, a, vx, res0);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    res0[ii] = -res0[ii] + rhs[ii];
    rk[ii] = res0[ii];
  }

  alpha = 1.0;
  betam1 = 1.0;
  gamma = 1.0;

  cvg = CS_SLES_ITERATING;

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Compute beta and omega;
       group dot products for new iteration's beta
       and previous iteration's residue to reduce total latency */

    if (n_iter == 0) {
      beta = _dot_product_xx(c, rk); /* rk == res0 here */
      residue = sqrt(beta);
      c->setup_data->initial_residue = residue;
    }
    else {
      _dot_products_xx_xy(c, rk, res0, &residue, &beta);
      residue = sqrt(residue);
    }

    /* Convergence test */
    cvg = _convergence_test(c, n_iter, residue, convergence);
    if (cvg != CS_SLES_ITERATING)
      break;

    n_iter += 1;

    if (_breakdown(c, convergence, "beta", beta, _epzero,
                   residue, n_iter, &cvg))
      break;

    if (_breakdown(c, convergence, "alpha", alpha, _epzero,
                   residue, n_iter, &cvg))
      break;

    omega = beta*gamma / (alpha*betam1);
    betam1 = beta;

    /* Compute pk */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      pk[ii] = rk[ii] + omega*(pk[ii] - alpha*uk[ii]);

    /* Compute zk = c.pk */

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 pk,
                                 zk,
                                 uk);

    /* Compute uk = A.zk */

    cs_matrix_vector_multiply(rotation_mode, a, zk, uk);

    /* Compute uk.res0 and gamma */

    ukres0 = _dot_product(c, uk, res0);

    gamma = beta / ukres0;

    /* First update of vx and rk */

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (gamma * zk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] -= (gamma * uk[jj]);
    }

    /* Compute zk = C.rk (zk is overwritten, vk is a working array */

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 rk,
                                 zk,
                                 vk);

    /* Compute vk = A.zk and alpha */

    cs_matrix_vector_multiply(rotation_mode, a, zk, vk);

    _dot_products_xx_xy(c, vk, rk, &ro_1, &ro_0);

    if (_breakdown(c, convergence, "rho1", ro_1, _epzero,
                   residue, n_iter, &cvg))
      break;

    alpha = ro_0 / ro_1;

    /* Final update of vx and rk */

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vx[ii] += (alpha * zk[ii]);

#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        rk[jj] -= (alpha * vk[jj]);
    }

    /* Convergence test at beginning of next iteration so
       as to group dot products for better parallel performance */
  }

  if (_aux_vectors != aux_vectors)
    BFT_FREE(_aux_vectors);

  return cvg;
}

/*----------------------------------------------------------------------------
 * Solution of (ad+ax).vx = Rhs using (not yet preconditioned) Bi-CGSTAB2.
 *
 * This Krylov method is based on an implemantation of J.P. Caltagirone
 * in his file bibpac6.f90 for Aquillon. He refers to it as Bi-CGSTAB2,
 * but a proper reference for the method has yet to be found.
 * (Gutknecht's BiCGstab2?)
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- block size of element ii, ii
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_bicgstab2(cs_sles_it_t              *c,
           const cs_matrix_t         *a,
           int                        diag_block_size,
           cs_halo_rotation_t         rotation_mode,
           cs_sles_it_convergence_t  *convergence,
           const cs_real_t           *rhs,
           cs_real_t                 *restrict vx,
           size_t                     aux_size,
           void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  cs_lnum_t ii, jj;
  double  _epzero = 1.e-30;/* smaller than epzero */
  double  ro_0, ro_1, alpha, beta, gamma;
  double  omega_1, omega_2, mu, nu, tau;
  double  residue;
  cs_real_t  *_aux_vectors;
  cs_real_t  *restrict res0, *restrict qk, *restrict rk, *restrict sk;
  cs_real_t  *restrict tk, *restrict uk, *restrict vk, *restrict wk;
  cs_real_t  *restrict zk;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;
    size_t  n_wa = 9;
    const size_t wa_size = CS_SIMD_SIZE(n_cols);

    if (aux_vectors == NULL || aux_size < (wa_size * n_wa))
      BFT_MALLOC(_aux_vectors, wa_size * n_wa, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    res0 = _aux_vectors;
    zk = _aux_vectors + wa_size;
    qk = _aux_vectors + wa_size*2;
    rk = _aux_vectors + wa_size*3;
    sk = _aux_vectors + wa_size*4;
    tk = _aux_vectors + wa_size*5;
    uk = _aux_vectors + wa_size*6;
    vk = _aux_vectors + wa_size*7;
    wk = _aux_vectors + wa_size*8;
  }

  /* Initialize work arrays */

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    uk[ii] = 0.0;
  }

  /* Initialize iterative calculation */
  /*----------------------------------*/

  cs_matrix_vector_multiply(rotation_mode, a, vx, res0);

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (ii = 0; ii < n_rows; ii++) {
    res0[ii] = -res0[ii] + rhs[ii];
    rk[ii] = res0[ii];
    qk[ii] = rk[ii];
  }

  ro_0    = 1.0;
  alpha   = 0.0;
  omega_2 = 1.0;

  cvg = CS_SLES_ITERATING;

  /* Current Iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    /* Compute beta and omega;
       group dot products for new iteration's beta
       and previous iteration's residue to reduce total latency */

    double mprec = 1.0e-60;

    if (n_iter == 0) {
      residue = sqrt(_dot_product_xx(c, rk)); /* rk == res0 here */
      c->setup_data->initial_residue = residue;
    }
    else
      residue = sqrt(_dot_product_xx(c, rk));

    /* Convergence test */
    cvg = _convergence_test(c, n_iter, residue, convergence);
    if (cvg != CS_SLES_ITERATING)
        break;

    n_iter += 1;

    ro_0 = -omega_2*ro_0;
    ro_1 = _dot_product(c, qk, rk);

    if (_breakdown(c, convergence, "rho0", ro_0, 1.e-60,
                   residue, n_iter, &cvg))
      break;

    if (_breakdown(c, convergence, "rho1", ro_1, _epzero,
                   residue, n_iter, &cvg))
      break;

    beta = alpha*ro_1/ro_0;
    ro_0 = ro_1;

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      uk[ii] = rk[ii] - beta* uk[ii];

    /* Compute vk =  A*uk */

    cs_matrix_vector_multiply(rotation_mode, a, uk, vk);

    /* Compute zk = c.vk */

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 vk,
                                 zk,
                                 sk);

    /* Compute gamma and alpha */

    gamma = _dot_product(c, qk, vk);

    if (_breakdown(c, convergence, "gamma", gamma, 1.e-60,
                   residue, n_iter, &cvg))
      break;

    alpha = ro_0/gamma;

    /* Update rk */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      rk[ii] -= alpha*vk[ii];
      vx[ii] += alpha*uk[ii];
    }

    /* p = A*r */

    cs_matrix_vector_multiply(rotation_mode, a, rk, sk);

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 sk,
                                 zk,
                                 wk);

    ro_1 = _dot_product(c, qk, sk);
    beta = alpha*ro_1/ro_0;
    ro_0 = ro_1;

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        vk[ii] = sk[ii] - beta*vk[ii];
#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        uk[jj] = rk[jj] - beta*uk[jj];
    }

    /* wk = A*vk */

    cs_matrix_vector_multiply(rotation_mode, a, vk, wk);

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 wk,
                                 zk,
                                 tk);

    gamma = _dot_product(c, qk, wk);
    alpha = (ro_0+mprec)/(gamma+mprec);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++) {
      rk[ii] -= alpha*vk[ii];
      sk[ii] -= alpha*wk[ii];
    }

    /* tk = A*sk */

    cs_matrix_vector_multiply(rotation_mode, a, sk, tk);

    _polynomial_preconditionning(c,
                                 rotation_mode,
                                 a,
                                 tk,
                                 zk,
                                 vk);

    _dot_products_xx_yy_xy_xz_yz(c, sk, tk, rk,
                                 &mu, &tau, &nu, &omega_1, &omega_2);

    tau = tau - (nu*nu/mu);
    omega_2 = (omega_2 - ((nu+mprec)*(omega_1+mprec)/(mu+mprec)))/(tau+mprec);

    omega_1 = (omega_1 - nu*omega_2)/mu;

    /* sol <- sol + omega_1*r + omega_2*s + alpha*u */
#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      vx[ii] += omega_1*rk[ii] + omega_2*sk[ii] + alpha*uk[ii];

#   pragma omp parallel if(n_rows > CS_THR_MIN)
    {
      /* r <- r - omega_1*s - omega_2*t */
#     pragma omp for nowait
      for (ii = 0; ii < n_rows; ii++)
        rk[ii] += - omega_1*sk[ii] - omega_2*tk[ii];
      /* u <- u - omega_1*v - omega_2*w */
#     pragma omp for nowait
      for (jj = 0; jj < n_rows; jj++)
        uk[jj] += - omega_1*vk[jj] - omega_2*wk[jj];
    }

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
 *   rotation_mode <-- Halo update option for rotational periodicity
 *   convergence   <-- Convergence information structure
 *   rhs           <-- Right hand side
 *   vx            <-> System solution
 *   aux_size      <-- Number of elements in aux_vectors (in bytes)
 *   aux_vectors   --- Optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_gmres(cs_sles_it_t              *c,
       const cs_matrix_t         *a,
       cs_halo_rotation_t         rotation_mode,
       cs_sles_it_convergence_t  *convergence,
       const cs_real_t           *rhs,
       cs_real_t                 *restrict vx,
       size_t                     aux_size,
       void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  int check_freq, l_iter, l_old_iter, scaltest;
  int krylov_size, _krylov_size;
  cs_lnum_t  ii, kk, jj;
  double    beta, dot_prod, residue, _residue, epsi;
  cs_real_t  *_aux_vectors;
  cs_real_t *restrict _krylov_vectors, *restrict _h_matrix;
  cs_real_t *restrict _givens_coeff, *restrict _beta;
  cs_real_t *restrict dk, *restrict gk;
  cs_real_t *restrict bk, *restrict fk, *restrict krk;

  int krylov_size_max = 75;
  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  const int diag_block_size = 1;

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = c->setup_data->n_rows;

  /* Allocate work arrays */

  krylov_size =  (krylov_size_max < (int)sqrt(n_rows)*1.5) ?
                  krylov_size_max : (int)sqrt(n_rows)*1.5 + 1;

#if defined(HAVE_MPI)
  if (c->comm != MPI_COMM_NULL) {
    MPI_Allreduce(&krylov_size,
                  &_krylov_size,
                  1,
                  MPI_INT,
                  MPI_MIN,
                  c->comm);
    krylov_size = _krylov_size;
  }
#endif

  check_freq = (int)(krylov_size/10) + 1;
  epsi = 1.e-15;
  scaltest = 0;

  {
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * diag_block_size;

    size_t _aux_r_size;
    size_t  n_wa = 4;
    size_t  wa_size = n_cols < krylov_size? krylov_size:n_cols;

    wa_size = CS_SIMD_SIZE(wa_size);
    _aux_r_size =   wa_size*n_wa
                  + (krylov_size-1)*(n_rows + krylov_size) + 3*krylov_size;

    if (aux_vectors == NULL || aux_size/sizeof(cs_real_t) < _aux_r_size)
      BFT_MALLOC(_aux_vectors, _aux_r_size, cs_real_t);
    else
      _aux_vectors = aux_vectors;

    dk = _aux_vectors;
    gk = _aux_vectors + wa_size;
    bk = _aux_vectors + 2*wa_size;
    fk = _aux_vectors + 3*wa_size;
    _krylov_vectors = _aux_vectors + n_wa*wa_size;
    _h_matrix = _aux_vectors + n_wa*wa_size + (krylov_size - 1)*n_rows;
    _givens_coeff =   _aux_vectors + n_wa*wa_size
                    + (krylov_size - 1)*(n_rows + krylov_size);
    _beta =   _aux_vectors + n_wa*wa_size
            + (krylov_size - 1)*(n_rows + krylov_size) + 2*krylov_size;
  }

  for (ii = 0; ii < krylov_size*(krylov_size - 1); ii++)
    _h_matrix[ii] = 0.;

  cvg = CS_SLES_ITERATING;

  while (cvg == CS_SLES_ITERATING) {

    /* compute  rk <- a*vx (vx = x0) */

    cs_matrix_vector_multiply(rotation_mode, a, vx, dk);

    /* compute  rk <- rhs - rk (r0 = b-A*x0) */

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (ii = 0; ii < n_rows; ii++)
      dk[ii] = rhs[ii] - dk[ii];

    if (n_iter == 0) {
      residue = sqrt(_dot_product_xx(c, dk));
      c->setup_data->initial_residue = residue;
      cvg = _convergence_test(c, n_iter, residue, convergence);
      if (cvg != CS_SLES_ITERATING)
        break;
    }

    /* beta = ||r0|| */
    beta = sqrt(_dot_product(c, dk, dk));
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

#     pragma omp parallel for if(n_rows > CS_THR_MIN)
      for (jj = 0; jj < n_rows; jj++)
        krk[jj] = dk[jj]/dot_prod;

      _polynomial_preconditionning(c,
                                   rotation_mode,
                                   a,
                                   krk,
                                   gk,
                                   dk);

      /* compute w = dk <- A*vj */

      cs_matrix_vector_multiply(rotation_mode, a, gk, dk);

      for (jj = 0; jj < ii + 1; jj++) {

        /* compute h(k,i) = <w,vi> = <dk,vi> */
        _h_matrix[ii*krylov_size + jj]
          = _dot_product(c, dk, (_krylov_vectors + jj*n_rows));

        /* compute w = dk <- w - h(i,k)*vi */

        cs_axpy(n_rows,
                -_h_matrix[ii*krylov_size+jj],
                (_krylov_vectors + jj*n_rows),
                dk);

      }

      /* compute h(i+1,i) = sqrt<w,w> */
      dot_prod = sqrt(_dot_product(c, dk, dk));
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

#       pragma omp parallel for private(kk) if(n_rows > CS_THR_MIN)
        for (jj = 0; jj < n_rows; jj++) {
          fk[jj] = 0.0;
          for (kk = 0; kk <= l_iter; kk++)
            fk[jj] += _krylov_vectors[kk*n_rows + jj] * gk[kk];
        }

        _polynomial_preconditionning(c,
                                     rotation_mode,
                                     a,
                                     fk,
                                     gk,
                                     bk);


#       pragma omp parallel for if(n_rows > CS_THR_MIN)
        for (jj = 0; jj < n_rows; jj++)
          fk[jj] = vx[jj] + gk[jj];

        cs_matrix_vector_multiply(rotation_mode, a, fk, bk);

        /* compute residue = | Ax - b |_1 */

        residue = 0.;
#       pragma omp parallel for reduction(+:residue) if(n_rows > CS_THR_MIN)
        for (jj = 0; jj < n_rows; jj++)
          residue += pow(rhs[jj] - bk[jj], 2);

#if defined(HAVE_MPI)

        if (c->comm != MPI_COMM_NULL) {
          MPI_Allreduce(&residue, &_residue, 1, MPI_DOUBLE, MPI_SUM,
                        c->comm);
          residue = _residue;
        }

#endif

        residue = sqrt(residue);

        cvg = _convergence_test(c, n_iter, residue, convergence);

      }

      n_iter++;
      l_iter++;

      if (cvg == 1 || l_iter == krylov_size - 1 || scaltest == 1) {
#       pragma omp parallel for if(n_rows > CS_THR_MIN)
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
 * Solution of A.vx = Rhs using Block Gauss-Seidel.
 *
 * On entry, vx is considered initialized.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

static cs_sles_convergence_state_t
_b_gauss_seidel(cs_sles_it_t              *c,
                const cs_matrix_t         *a,
                int                        diag_block_size,
                cs_halo_rotation_t         rotation_mode,
                cs_sles_it_convergence_t  *convergence,
                const cs_real_t           *rhs,
                cs_real_t                 *restrict vx,
                size_t                     aux_size,
                void                      *aux_vectors)
{
  cs_sles_convergence_state_t cvg;
  double  res2, residue;

  unsigned n_iter = 0;

  /* Call setup if not already done, allocate or map work arrays */
  /*-------------------------------------------------------------*/

  if (c->setup_data == NULL)
    _setup_sles_it(c, a, diag_block_size, false);

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(a);

  const cs_halo_t *halo = cs_matrix_get_halo(a);

  const cs_real_t  *restrict ad_inv = c->setup_data->ad_inv;

  const cs_real_t  *restrict ad = cs_matrix_get_diagonal(a);

  const cs_lnum_t  *order = NULL;

  const cs_lnum_t  *a_row_index, *a_col_id;
  const cs_real_t  *a_d_val, *a_x_val;

  const int *db_size = cs_matrix_get_diag_block_size(a);
  cs_matrix_get_msr_arrays(a, &a_row_index, &a_col_id, &a_d_val, &a_x_val);

  if (c->add_data != NULL)
    order = c->add_data->order;

  cvg = CS_SLES_ITERATING;

  /* Current iteration */
  /*-------------------*/

  while (cvg == CS_SLES_ITERATING) {

    register double r;

    n_iter += 1;

    /* Synchronize ghost cells first */

    if (halo != NULL) {

      if (db_size[3] == 1)
        cs_halo_sync_component(halo,
                               CS_HALO_STANDARD,
                               rotation_mode,
                               vx);

      else { /* if (matrix->db_size[3] > 1) */

        cs_halo_sync_var_strided(halo,
                                 CS_HALO_STANDARD,
                                 vx,
                                 db_size[1]);

        /* Synchronize periodic values */

        if (halo->n_transforms > 0 && db_size[0] == 3)
          cs_halo_perio_sync_var_vect(halo,
                                      CS_HALO_STANDARD,
                                      vx,
                                      db_size[1]);

      }

    }

    /* Compute Vx <- Vx - (A-diag).Rk and residue. */

    res2 = 0.0;

    if (order == NULL) {

      if (diag_block_size == 1) {

#       pragma omp parallel for private(r) reduction(+:res2)      \
                            if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

          const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
          const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
          const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

          cs_real_t vxm1 = vx[ii];

          vx[ii] = rhs[ii];

          for (cs_lnum_t jj = 0; jj < n_cols; jj++)
            vx[ii] -= (m_row[jj]*vx[col_id[jj]]);

          vx[ii] *= ad_inv[ii];

          r = ad[ii] * (vx[ii]-vxm1);
          res2 += (r*r);
        }

      }
      else {

#       pragma omp parallel for private(r) reduction(+:res2) \
                          if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

          const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
          const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
          const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

          cs_real_t vxm1;

          for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {

            vxm1 = vx[ii*db_size[1] + kk];

            vx[ii*db_size[1] + kk] = rhs[ii*db_size[1] + kk];

            for (cs_lnum_t jj = 0; jj < kk; jj++)
              vx[ii*db_size[1] + kk]
                -=   ad[ii*db_size[3] + kk*db_size[2] + jj]
                   * vx[ii*db_size[1] + jj];
            for (cs_lnum_t jj = kk+1; jj < db_size[0]; jj++)
              vx[ii*db_size[1] + kk]
                -=   ad[ii*db_size[3] + kk*db_size[2] + jj]
                   * vx[ii*db_size[1] + jj];

            for (cs_lnum_t jj = 0; jj < n_cols; jj++)
              vx[ii*db_size[1] + kk]
                -= (m_row[jj]*vx[col_id[jj]*db_size[1] + kk]);

            vx[ii*db_size[1] + kk] *= ad_inv[ii*db_size[1] + kk];

            r = ad[ii*db_size[1] + kk] * (vx[ii*db_size[1] + kk]-vxm1);
            res2 += (r*r);

          }

        }

      }

    }

    /* Ordered variant */

    else {

      if (diag_block_size == 1) {

#       pragma omp parallel for private(r) reduction(+:res2)      \
                            if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

          cs_lnum_t ii = order[ll];

          const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
          const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
          const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

          cs_real_t vxm1 = vx[ii];

          vx[ii] = rhs[ii];

          for (cs_lnum_t jj = 0; jj < n_cols; jj++)
            vx[ii] -= (m_row[jj]*vx[col_id[jj]]);

          vx[ii] *= ad_inv[ii];

          r = ad[ii] * (vx[ii]-vxm1);
          res2 += (r*r);
        }

      }
      else {

#       pragma omp parallel for private(r) reduction(+:res2) \
                            if(n_rows > CS_THR_MIN)
        for (cs_lnum_t ll = 0; ll < n_rows; ll++) {

          cs_lnum_t ii = order[ll];

          const cs_lnum_t *restrict col_id = a_col_id + a_row_index[ii];
          const cs_real_t *restrict m_row = a_x_val + a_row_index[ii];
          const cs_lnum_t n_cols = a_row_index[ii+1] - a_row_index[ii];

          cs_real_t vxm1;

          for (cs_lnum_t kk = 0; kk < db_size[0]; kk++) {

            vxm1 = vx[ii*db_size[1] + kk];

            vx[ii*db_size[1] + kk] = rhs[ii*db_size[1] + kk];

            for (cs_lnum_t jj = 0; jj < kk; jj++)
              vx[ii*db_size[1] + kk]
                -=   ad[ii*db_size[3] + kk*db_size[2] + jj]
                   * vx[ii*db_size[1] + jj];
            for (cs_lnum_t jj = kk+1; jj < db_size[0]; jj++)
              vx[ii*db_size[1] + kk]
                -=   ad[ii*db_size[3] + kk*db_size[2] + jj]
                   * vx[ii*db_size[1] + jj];

            for (cs_lnum_t jj = 0; jj < n_cols; jj++)
              vx[ii*db_size[1] + kk]
                -= (m_row[jj]*vx[col_id[jj]*db_size[1] + kk]);

            vx[ii*db_size[1] + kk] *= ad_inv[ii*db_size[1] + kk];

            r = ad[ii*db_size[1] + kk] * (vx[ii*db_size[1] + kk]-vxm1);
            res2 += (r*r);

          }

        }

      }

    }

#if defined(HAVE_MPI)

    if (c->comm != MPI_COMM_NULL) {
      double _sum;
      MPI_Allreduce(&res2, &_sum, 1, MPI_DOUBLE, MPI_SUM,
                    c->comm);
      res2 = _sum;
    }

#endif /* defined(HAVE_MPI) */

    residue = sqrt(res2); /* Actually, residue of previous iteration */

    /* Convergence test */

    if (n_iter == 1)
      c->setup_data->initial_residue = residue;

    cvg = _convergence_test(c, n_iter, residue, convergence);

  }

  return cvg;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate an iterative sparse linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_it_create.
 *
 * Note that this function returns a pointer directly to the iterative solver
 * management structure. This may be used to set further options,
 * for example using \ref cs_sles_it_set_pcg_single_reduction. If needed,
 * \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]  f_id          associated field id, or < 0
 * \param[in]  name          associated name if f_id < 0, or NULL
 * \param[in]  solver_type   type of solver (PCG, Jacobi, ...)
 * \param[in]  poly_degree   preconditioning polynomial degree
 *                           (0: diagonal; -1: non-preconditioned)
 * \param[in]  n_max_iter    maximum number of iterations
 *
 * \return  pointer to newly created iterative solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_define(int                 f_id,
                  const char         *name,
                  cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter)
{
  cs_sles_it_t *
    c = cs_sles_it_create(solver_type,
                          poly_degree,
                          n_max_iter,
                          true); /* update stats */

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_it_t",
                                 cs_sles_it_setup,
                                 cs_sles_it_solve,
                                 cs_sles_it_free,
                                 cs_sles_it_log,
                                 cs_sles_it_copy,
                                 cs_sles_it_destroy);

  cs_sles_set_error_handler(sc,
                            cs_sles_it_error_post_and_abort);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context.
 *
 * \param[in]  solver_type   type of solver (PCG, Jacobi, ...)
 * \param[in]  poly_degree   preconditioning polynomial degree
 *                           (0: diagonal; -1: non-preconditioned;
 *                           see \ref sles_it for details)
 * \param[in]  n_max_iter    maximum number of iterations
 * \param[in]  update_stats  automatic solver statistics indicator
 *
 * \return  pointer to newly created solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_t *
cs_sles_it_create(cs_sles_it_type_t   solver_type,
                  int                 poly_degree,
                  int                 n_max_iter,
                  bool                update_stats)
{
  cs_sles_it_t *c;

  BFT_MALLOC(c, 1, cs_sles_it_t);

  c->type = solver_type;

  c->poly_degree = poly_degree;
  if (c->type == CS_SLES_JACOBI)
    c->poly_degree = 0;

  c->update_stats = update_stats;

  c->n_max_iter = n_max_iter;

  c->n_setups = 0;
  c->n_solves = 0;

  c->n_iterations_min = 0;
  c->n_iterations_max = 0;
  c->n_iterations_last = 0;
  c->n_iterations_tot = 0;

  CS_TIMER_COUNTER_INIT(c->t_setup);
  CS_TIMER_COUNTER_INIT(c->t_solve);

#if defined(HAVE_MPI)
  c->comm = cs_glob_mpi_comm;
#endif

  c->setup_data = NULL;
  c->add_data = NULL;
  c->shared = NULL;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy iterative sparse linear system solver info and context.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_it_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_destroy(void **context)
{
  cs_sles_it_t *c = (cs_sles_it_t *)(*context);
  if (c != NULL) {
    cs_sles_it_free(c);
    if (c->add_data != NULL) {
      BFT_FREE(c->add_data->order);
      BFT_FREE(c->add_data);
    }
    BFT_FREE(c);
    *context = c;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create iterative sparse linear system solver info and context
 *        based on existing info and context.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_it_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_it_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_it_copy(const void  *context)
{
  cs_sles_it_t *d = NULL;

  if (context != NULL) {
    const cs_sles_it_t *c = context;
    d = cs_sles_it_create(c->type,
                          c->poly_degree,
                          c->n_max_iter,
                          c->update_stats);
#if defined(HAVE_MPI)
    d->comm = c->comm;
#endif
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to iterative solver info and context
 *                       (actual type: cs_sles_it_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_log(const void  *context,
               cs_log_t     log_type)
{
  const cs_sles_it_t  *c = context;

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       %s\n"
                    "  Preconditioning:                   "),
                  _(cs_sles_it_type_name[c->type]));
    if (c->poly_degree < 0)
      cs_log_printf(log_type, _("none\n"));
    else if (c->poly_degree == 0)
      cs_log_printf(log_type, _("diagonal\n"));
    else
      cs_log_printf(log_type, _("polynomial, degree %d\n"),
                    c->poly_degree);
    cs_log_printf(log_type,
                  _("  Maximum number of iterations:      %d\n"),
                  c->n_max_iter);

  }

  else if (log_type == CS_LOG_PERFORMANCE) {

    int n_calls = c->n_solves;
    int n_it_min = c->n_iterations_min;
    int n_it_max = c->n_iterations_max;
    int n_it_mean = 0;

    if (n_calls > 0)
      n_it_mean = (int)(  c->n_iterations_tot
                         / ((unsigned long long)n_calls));

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   %s\n"
                    "  Number of setups:              %12d\n"
                    "  Number of calls:               %12d\n"
                    "  Minimum number of iterations:  %12d\n"
                    "  Maximum number of iterations:  %12d\n"
                    "  Mean number of iterations:     %12d\n"
                    "  Total setup time:              %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  _(cs_sles_it_type_name[c->type]),
                  c->n_setups, n_calls, n_it_min, n_it_max, n_it_mean,
                  c->t_setup.wall_nsec*1e-9,
                  c->t_solve.wall_nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup iterative sparse linear equation solver.
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 *                             (actual type: cs_sles_it_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_setup(void               *context,
                 const char         *name,
                 const cs_matrix_t  *a,
                 int                 verbosity)
{
  cs_sles_it_t  *c = context;

  const int diag_block_size = (cs_matrix_get_diag_block_size(a))[0];

  if (c->type == CS_SLES_JACOBI)
    _setup_sles_it(c, a, diag_block_size, true);

  else
    _setup_sles_it(c, a, diag_block_size, false);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call iterative sparse linear equation solver.
 *
 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_it_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       rotation_mode  halo update option for rotational periodicity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residue normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residue        residue
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if NULL)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_solve(void                *context,
                 const char          *name,
                 const cs_matrix_t   *a,
                 int                  verbosity,
                 cs_halo_rotation_t   rotation_mode,
                 double               precision,
                 double               r_norm,
                 int                 *n_iter,
                 double              *residue,
                 const cs_real_t     *rhs,
                 cs_real_t           *vx,
                 size_t               aux_size,
                 void                *aux_vectors)
{
  cs_sles_it_t  *c = context;

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_timer_t t0, t1;

  unsigned _n_iter = 0;
  cs_sles_it_convergence_t  convergence;

  if (c->update_stats == true)
    t0 = cs_timer_time();

  const int *diag_block_size = cs_matrix_get_diag_block_size(a);
  const int _diag_block_size = diag_block_size[0];

  assert(diag_block_size[0] == diag_block_size[1]);

  /* Initialize number of iterations and residue,
     check for immediate return,
     and solve sparse linear system */

  *n_iter = 0;

  /* Setup if not already done */

  if (c->setup_data == NULL) {

    if (c->update_stats) { /* Stop solve timer to switch to setup timer */
      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);
    }

    cs_sles_it_setup(c, name, a, verbosity);

    if (c->update_stats) /* Restart solve timer */
      t0 = cs_timer_time();

  }

  /* Solve sparse linear system */

  _convergence_init(&convergence,
                    _(cs_sles_it_type_name[c->type]),
                    name,
                    verbosity,
                    c->n_max_iter,
                    precision,
                    r_norm,
                    residue);

  c->setup_data->initial_residue = -1;

  /* Only call solver for "active" ranks */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks < 2 || c->comm != MPI_COMM_NULL) {
#endif

    switch (c->type) {
    case CS_SLES_PCG:
      if (! c->setup_data->single_reduce) {
        if (c->poly_degree > -1)
          cvg = _conjugate_gradient(c,
                                    a,
                                    _diag_block_size,
                                    rotation_mode,
                                    &convergence,
                                    rhs,
                                    vx,
                                    aux_size,
                                    aux_vectors);
        else
          cvg = _conjugate_gradient_npc(c,
                                        a,
                                        _diag_block_size,
                                        rotation_mode,
                                        &convergence,
                                        rhs,
                                        vx,
                                        aux_size,
                                        aux_vectors);
        break;
      }
      else {
        if (c->poly_degree > -1)
          cvg = _conjugate_gradient_sr(c,
                                       a,
                                       _diag_block_size,
                                       rotation_mode,
                                       &convergence,
                                       rhs,
                                       vx,
                                       aux_size,
                                       aux_vectors);
        else
          cvg = _conjugate_gradient_npc_sr(c,
                                           a,
                                           _diag_block_size,
                                           rotation_mode,
                                           &convergence,
                                           rhs,
                                           vx,
                                           aux_size,
                                           aux_vectors);
      }
      break;
    case CS_SLES_JACOBI:
      if (_diag_block_size == 1)
        cvg = _jacobi(c,
                      a,
                      _diag_block_size,
                      rotation_mode,
                      &convergence,
                      rhs,
                      vx,
                      aux_size,
                      aux_vectors);
      else if (_diag_block_size == 3)
        cvg = _block_3_jacobi(c,
                              a,
                              rotation_mode,
                              &convergence,
                              rhs,
                              vx,
                              aux_size,
                              aux_vectors);
      break;
    case CS_SLES_BICGSTAB:
      cvg = _bi_cgstab(c,
                       a,
                       _diag_block_size,
                       rotation_mode,
                       &convergence,
                       rhs,
                       vx,
                       aux_size,
                       aux_vectors);
      break;
    case CS_SLES_BICGSTAB2:
      cvg = _bicgstab2(c,
                       a,
                       _diag_block_size,
                       rotation_mode,
                       &convergence,
                       rhs,
                       vx,
                       aux_size,
                       aux_vectors);
      break;
    case CS_SLES_GMRES:
      if (_diag_block_size == 1)
        cvg = _gmres(c,
                     a,
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
    case CS_SLES_B_GAUSS_SEIDEL:
      cvg = _b_gauss_seidel(c,
                            a,
                            _diag_block_size,
                            rotation_mode,
                            &convergence,
                            rhs,
                            vx,
                            aux_size,
                            aux_vectors);
      break;
    default:
      bft_error
        (__FILE__, __LINE__, 0,
         _("Resolution of linear equation on \"%s\"\n"
           "with solver type %d, which is not defined)."),
         name, (int)c->type);
      break;
    }

#if defined(HAVE_MPI)
  }
#endif

  /* Broadcast convergence info from "active" ranks to others*/

#if defined(HAVE_MPI)
  if (c->comm != cs_glob_mpi_comm) {
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

  if (c->update_stats == true) {

    t1 = cs_timer_time();

    if (c->n_solves == 0)
      c->n_iterations_min = _n_iter;

    c->n_solves += 1;

    if (c->n_iterations_min > _n_iter)
      c->n_iterations_min = _n_iter;
    if (c->n_iterations_max < _n_iter)
      c->n_iterations_max = _n_iter;

    c->n_iterations_last = _n_iter;
    c->n_iterations_tot += _n_iter;

    cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  }

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free iterative sparse linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 *                           (actual type: cs_sles_it_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_free(void  *context)
{
  cs_sles_it_t  *c = context;

  cs_timer_t t0;
  if (c->update_stats == true)
    t0 = cs_timer_time();

  if (c->setup_data != NULL) {
    BFT_FREE(c->setup_data->_ad_inv);
    BFT_FREE(c->setup_data);
  }

  if (c->update_stats == true) {
    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(c->t_setup), &t0, &t1);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return iterative solver type.
 *
 * \param[in]  context  pointer to iterative solver info and context
 *
 * \return  selected solver type
 */
/*----------------------------------------------------------------------------*/

cs_sles_it_type_t
cs_sles_it_get_type(const cs_sles_it_t  *context)
{
  return context->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the initial residue for the previous solve operation
 *        with a solver.
 *
 * This is useful for convergence tests when this solver is used as
 * a preconditioning smoother.
 *
 * This operation is only valid between calls to \ref cs_sles_it_setup
 * (or \ref cs_sles_it_solve) and \ref cs_sles_it_free.
 * It returns -1 otherwise.
 *
 * \param[in]  context  pointer to iterative solver info and context
 *
 * \return initial residue from last call to \ref cs_sles_solve with this
 *         solver
 */
/*----------------------------------------------------------------------------*/

double
cs_sles_it_get_last_initial_residue(const cs_sles_it_t  *context)
{
  double retval = 1;
  if (context->setup_data != NULL)
    retval = context->setup_data->initial_residue;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a similar info and context object with which some setup
 *        data may be shared.
 *
 * This is especially useful for sharing preconditioning data between
 * similar solver contexts (for example ascending and descending multigrid
 * smoothers based on the same matrix).
 *
 * For preconditioning data to be effectively shared, \ref cs_sles_it_setup
 * (or \ref cs_sles_it_solve) must be called on \p shareable before being
 * called on \p context (without \ref cs_sles_it_free being called in between,
 * of course).
 *
 * It is the caller's responsibility to ensure the context is not used
 * for a \ref cs_sles_it_setup or \ref cs_sles_it_solve operation after the
 * shareable object has been destroyed (normally by \ref cs_sles_it_destroy).
 *
 * \param[in, out]  context    pointer to iterative solver info and context
 * \param[in]       shareable  pointer to iterative solver info and context
 *                             whose context may be shared
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_shareable(cs_sles_it_t        *context,
                         const cs_sles_it_t  *shareable)
{
  cs_sles_it_t  *c = context;

  c->shared = shareable;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set MPI communicator for dot products.
 *
 * \param[in, out]  context   pointer to iterative solver info and context
 * \param[in]       comm      MPI communicator
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_mpi_reduce_comm(cs_sles_it_t  *context,
                               MPI_Comm       comm)
{
  cs_sles_it_t  *c = context;

  static int flag = -1;

  if (flag < 0)
    flag = cs_halo_get_use_barrier();

  c->comm = comm;

  if (comm != cs_glob_mpi_comm)
    cs_halo_set_use_barrier(0);
  else {
    cs_halo_set_use_barrier(flag);
    if (cs_glob_n_ranks < 2)
      c->comm = MPI_COMM_NULL;
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign ordering to iterative solver.
 *
 * The solver context takes ownership of the order array (i.e. it will
 * handle its later deallocation).
 *
 * This is useful only for Block Gauss-Seidel.
 *
 * \param[in, out]  context  pointer to iterative solver info and context
 * \param[in, out]  order    pointer to ordering array
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_assign_order(cs_sles_it_t   *context,
                        cs_lnum_t     **order)
{
  if (context->type != CS_SLES_B_GAUSS_SEIDEL)
    BFT_FREE(*order);

  else {

    if (context->add_data == NULL) {
      BFT_MALLOC(context->add_data, 1, cs_sles_it_add_t);
      context->add_data->order = NULL;
    }

    BFT_FREE(context->add_data->order);

    context->add_data->order = *order;

    *order = NULL;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query mean number of rows under which Conjugate Gradient algorithm
 *        uses the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs, so 0 is returned.
 *
 * \returns  mean number of rows per active rank under which the
 *           single-reduction variant will be used
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_sles_it_get_pcg_single_reduction(void)
{
#if defined(HAVE_MPI)
  return _pcg_sr_threshold;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set mean number of rows under which Conjugate Gradient algorithm
 *        should use the single-reduction variant.
 *
 * The single-reduction variant requires only one parallel sum per
 * iteration (instead of 2), at the cost of additional vector operations,
 * so it tends to be more expensive when the number of matrix rows per
 * MPI rank is high, then becomes cheaper when the MPI latency cost becomes
 * more significant.
 *
 * This option is ignored for non-parallel runs.
 *
 * \param[in]  threshold  mean number of rows per active rank under which the
 *             single-reduction variant will be used
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_set_pcg_single_reduction(cs_lnum_t  threshold)
{
#if defined(HAVE_MPI)
  _pcg_sr_threshold = threshold;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the current global settings relative to parallelism.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_log_parallel_options(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "Iterative linear solvers parallel parameters:\n"
                    "  PCG single-reduction threshold:     %d\n"),
                 _pcg_sr_threshold);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler for iterative sparse linear equation solver.
 *
 * In case of divergence or breakdown, this error handler outputs
 * postprocessing data to assist debugging, then aborts the run.
 * It does nothing in case the maximum iteration count is reached.

 * \param[in, out]  context        pointer to iterative solver info and context
 *                                 (actual type: cs_sles_it_t  *)
 * \param[in]       state          convergence state
 * \param[in]       name           pointer to name of linear system
 * \param[in]       a              matrix
 * \param[in]       rotation_mode  halo update option for rotational periodicity
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_it_error_post_and_abort(void                         *context,
                                cs_sles_convergence_state_t   state,
                                const char                   *name,
                                const cs_matrix_t            *a,
                                cs_halo_rotation_t            rotation_mode,
                                const cs_real_t              *rhs,
                                cs_real_t                    *vx)
{
  if (state >= CS_SLES_BREAKDOWN)
    return;

  cs_sles_it_t  *c = context;

  int mesh_id = cs_post_init_error_writer_cells();

  cs_sles_post_error_output_def(name,
                                mesh_id,
                                rotation_mode,
                                a,
                                rhs,
                                vx);

  cs_post_finalize();

  const char *error_type[] = {N_("divergence"), N_("breakdown")};
  int err_id = (state == CS_SLES_BREAKDOWN) ? 1 : 0;

  bft_error(__FILE__, __LINE__, 0,
            _("%s: error (%s) solving for %s"),
            _(cs_sles_it_type_name[c->type]),
            _(error_type[err_id]),
            name);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
