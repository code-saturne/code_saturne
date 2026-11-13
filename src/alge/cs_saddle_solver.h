#ifndef CS_SADDLE_SOLVER_H
#define CS_SADDLE_SOLVER_H

/*============================================================================
 * Solvers for saddle-point systems arising from CDO discretizations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cdo/cs_cdo_blas.h"
#include "cdo/cs_cdo_system.h"
#include "cdo/cs_iter_algo.h"
#include "cdo/cs_property.h"
#include "alge/cs_param_saddle.h"
#include "alge/cs_sles.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Generic function prototype to perform a matrix vector operation
 *        This operation takes place between an unassembled matrix and a vector
 *
 * \param[in]      n2_dofs  number of (scatter) DoFs for the (2,2)-block
 * \param[in]      vec      array of values
 * \param[in]      mat_adj  adjacency related to the matrix operator
 * \param[in]      mat_op   values associated to the matrix operator
 * \param[in, out] matvec   resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_saddle_solver_matvec_t)(cs_lnum_t              n2_dofs,
                            const cs_real_t       *vec,
                            const cs_adjacency_t  *mat_adj,
                            const cs_real_t       *mat_op,
                            cs_real_t             *matvec);

/*
 * Main structure to solve a saddle-point problem described as follows
 *
 *      | M11 | M12 | | x1 |   |rhs1 |
 *  M = |-----------| |----| = |-----|
 *      | M21 |  0  | | x2 |   |rhs2 |
 *
 *  One assumes that M12 = M21^T
 */

typedef struct {

  /* Set of parameters to know how to solve the saddle-point system (shared) */

  const cs_param_saddle_t  *param;

  /* Structure to handle iterative algorithms */

  cs_iter_algo_t           *algo;

  /* Description of the saddle-point system to solve (shared) */

  cs_cdo_system_helper_t   *system_helper;

  bool                      do_setup;

  /* Main SLES structure associated to the saddle-point problem. According to
   * the type of saddle-point solver, additional SLES can be allocated, for
   * instance to compute an approximation of the Schur complement. */

  cs_sles_t                *main_sles;

  /* Scatter viewpoint */

  cs_lnum_t                 n1_elts;
  int                       n1_dofs_by_elt;
  cs_lnum_t                 n2_elts;
  int                       n2_dofs_by_elt;

  cs_lnum_t                 n1_scatter_dofs; /* = n1_elts * n1_dofs_by_elt */
  cs_lnum_t                 n2_scatter_dofs; /* = n2_elts * n2_dofs_by_elt */

  /* Additional members according to the type of saddle-point solver. This
   * structure is cast on-the-fly */

  void                     *context;

  /* Monitoring */

  unsigned                  n_calls;
  unsigned                  n_iter_min;
  unsigned                  n_iter_max;
  unsigned                  n_iter_tot;

} cs_saddle_solver_t;


/* ==================================================== */
/* List of context structures dedicated to some solvers */
/* ==================================================== */

/* Context structure for the ALU algorithm */
/* --------------------------------------- */

typedef struct {

  /* Auxiliary buffers */

  cs_real_t  *inv_m22;  /* reciprocal of the mass matrix for the (2,2)
                         * block; buffer of size n2_dofs */
  cs_real_t  *res2;     /* buffer of size n2_dofs */
  cs_real_t  *m21x1;    /* buffer of size n2_dofs */

  cs_real_t  *b1_tilda; /* Modified RHS (size n1_dofs) */
  cs_real_t  *rhs;      /* buffer of size n1_dofs */

  /* SLES associated to the transformation of the right-hand side */

  cs_sles_t  *init_sles;

  /* Function pointers for computing operations needed in the algorithm */
  /* ------------------------------------------------------------------ */

  cs_cdo_blas_square_norm_t  *square_norm_b11;
  cs_saddle_solver_matvec_t  *m12_vector_multiply;
  cs_saddle_solver_matvec_t  *m21_vector_multiply;

  /* Shared pointers */

  const cs_real_t            *m21_val;

  /* Indexed list used to scan the unassembled m21 operator */

  const cs_adjacency_t       *m21_adj;

} cs_saddle_solver_context_alu_t;


/* Context structure for block preconditioner used in combination with a Krylov
 * solver such as MINRES or GCR.
 * ----------------------------------------------------------------------------
 *
 * A hybrid storage is used to represent the system described below
 *
 *      | M11 | M12 | | x1 |   |rhs1 |
 *  M = |-----------| |----| = |-----|
 *      | M21 |  0  | | x2 |   |rhs2 |
 *
 *  One assumes that M12 = M21^T
 */

typedef struct {

  /* Shortcut on pointers available through the system helper */

  cs_matrix_t      *m11;
  cs_range_set_t   *b11_range_set;

  /* Max. size of the blocks (scatter or gather view). This size takes into
     account the size needed for synchronization. */

  cs_lnum_t         b11_max_size;
  cs_lnum_t         b22_max_size;

  /* SLES structure associated to the Schur complement. It depends on the type
     of Schur complement approximation used */

  cs_matrix_t      *schur_matrix;
  cs_sles_t        *schur_sles;
  double            schur_scaling;

  /* Native arrays for the Schur matrix (optional) */

  cs_real_t        *schur_diag;
  cs_real_t        *schur_xtra;

  /* Diagonal approximations of block matrices (optional) */

  cs_real_t        *m11_inv_diag;
  cs_real_t        *m22_mass_diag;

  cs_sles_t        *xtra_sles;

  /* Function pointers */

  cs_saddle_solver_matvec_t  *m12_vector_multiply;
  cs_saddle_solver_matvec_t  *m21_vector_multiply;

  /* Shared pointers */

  const cs_property_t   *pty_22;  /* Property related to the (2,2) block */
  const cs_real_t       *m21_val;
  const cs_adjacency_t  *m21_adj; /* Indexed list used to scan the unassembled
                                     m21 operator */

} cs_saddle_solver_context_block_pcd_t;


/* Context structure for the GKB algorithm */
/* --------------------------------------- */
/* This structure follows the algorithm given in the article entitled "An
 * iterative generalized Golub-Kahan algorithm for problems in structural
 * mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 */

typedef struct {

  /* Orthogonalization coefficients */

  cs_real_t   alpha;
  cs_real_t   beta;
  cs_real_t   zeta;

  /* Energy norm estimation */

  int         zeta_size;
  cs_real_t  *zeta_array;
  cs_real_t   zeta_square_sum;

  /* Auxiliary buffers */

  cs_real_t  *q;         /* buffer of size n2_dofs */
  cs_real_t  *d;         /* buffer of size n2_dofs */
  cs_real_t  *m21v;      /* buffer of size n2_dofs */
  cs_real_t  *inv_m22;   /* reciprocal of the mass matrix for the (2,2)
                          * block; buffer of size n2_dofs */
  const cs_real_t  *m22; /* mass matrix for the (2,2) block (shared pointer) */

  cs_real_t  *w;         /* buffer of size n1_dofs */
  cs_real_t  *m12q;      /* buffer of size n1_dofs */
  cs_real_t  *v;         /* buffer of size n1_dofs */
  cs_real_t  *x1_tilda;  /* buffer of size n1_dofs */

  cs_real_t  *rhs_tilda;  /* buffer of size MAX(n1_dofs, n2_dofs) */

  /* Optional SLES to perform the transformation of the right hand-side */

  cs_sles_t  *init_sles;

  /* Function pointers */

  cs_cdo_blas_square_norm_t  *square_norm_b11;
  cs_saddle_solver_matvec_t  *m12_vector_multiply;
  cs_saddle_solver_matvec_t  *m21_vector_multiply;

  /* Shared pointers */

  const cs_real_t            *m21_val;

  /* Indexed list used to scan the unassembled m21 operator */

  const cs_adjacency_t       *m21_adj;

} cs_saddle_solver_context_gkb_t;


/* Context structure for the Notay's algebraic transformation */
/* ---------------------------------------------------------- */

typedef struct {

  /* Function pointer */

  cs_saddle_solver_matvec_t  *m12_vector_multiply;

  /* Shared pointers */

  const cs_real_t            *m21_val;

  /* Indexed list used to scan the unassembled m21 operator */

  const cs_adjacency_t       *m21_adj;

} cs_saddle_solver_context_notay_t;


/* Context structure for the Uzawa-CG algorithm */
/* -------------------------------------------- */

typedef struct {

  /* Auxiliary scaling coefficient */

  cs_real_t   alpha;

  /* Auxiliary buffers */

  cs_real_t  *res2;     /* buffer of size n2_dofs */
  cs_real_t  *m21x1;    /* buffer of size n2_dofs */
  cs_real_t  *gk;       /* buffer of size n2_dofs */

  cs_real_t  *b1_tilda; /* Modified RHS (size n1_dofs) */
  cs_real_t  *dzk;      /* buffer of size n1_dofs */
  cs_real_t  *rhs;      /* buffer of size n1_dofs */

  /* Function pointers */

  cs_cdo_blas_square_norm_t  *square_norm_b11;
  cs_saddle_solver_matvec_t  *m12_vector_multiply;
  cs_saddle_solver_matvec_t  *m21_vector_multiply;

  /* Shortcut on pointers available through the system helper */

  cs_matrix_t           *m11;
  cs_range_set_t        *b11_range_set;

  /* Max. size of the blocks (scatter or gather view). This size takes into
     account the size needed for synchronization. */

  cs_lnum_t              b11_max_size;
  cs_lnum_t              b22_max_size;

  /* SLES structure associated to the Schur complement. It depends on the type
     of Schur complement approximation used */

  cs_real_t             *inv_m22;  /* reciprocal of the mass matrix for the
                                    * (2,2) block; buffer of size n2_dofs */
  cs_matrix_t           *schur_matrix;
  cs_sles_t             *schur_sles;

  /* Native arrays for the Schur matrix (optional) */

  cs_real_t             *schur_diag;
  cs_real_t             *schur_xtra;

  /* Diagonal approximations of block matrices (optional) */

  cs_real_t             *m11_inv_diag;

  cs_sles_t             *xtra_sles;
  cs_sles_t             *init_sles;

  /* Shared pointers */

  const cs_property_t   *pty_22;  /* Property related to the (2,2) block */
  const cs_real_t       *m21_val;
  const cs_adjacency_t  *m21_adj; /* Indexed list used to scan the unassembled
                                     m21 operator */

} cs_saddle_solver_context_uzawa_cg_t;

/* Context structure for the SIMPLE-lile algorithm */
/* ----------------------------------------------- */

typedef struct {

  /* Auxiliary buffers */

  cs_real_t  *res2;     /* buffer of size n2_dofs */
  cs_real_t  *m21x1;    /* buffer of size n2_dofs */

  cs_real_t  *b1_tilda; /* Modified RHS (size n1_dofs) */
  cs_real_t  *rhs;      /* buffer of size n1_dofs */

  /* Function pointers */

  cs_cdo_blas_square_norm_t  *square_norm_b11;
  cs_saddle_solver_matvec_t  *m12_vector_multiply;
  cs_saddle_solver_matvec_t  *m21_vector_multiply;

  /* Shortcut on pointers available through the system helper */

  cs_matrix_t           *m11;
  cs_range_set_t        *b11_range_set;

  /* Max. size of the blocks (scatter or gather view). This size takes into
     account the size needed for synchronization. */

  cs_lnum_t              b11_max_size;
  cs_lnum_t              b22_max_size;

  /* SLES structure associated to the Schur complement. It depends on the type
     of Schur complement approximation used */

  cs_real_t             *inv_m22;  /* reciprocal of the mass matrix for the
                                    * (2,2) block; buffer of size n2_dofs */
  cs_matrix_t           *schur_matrix;
  cs_sles_t             *schur_sles;

  /* Native arrays for the Schur matrix (optional) */

  cs_real_t             *schur_diag;
  cs_real_t             *schur_xtra;

  /* Diagonal approximations of block matrices (optional) */

  cs_real_t             *m11_inv_diag;

  cs_sles_t             *xtra_sles;
  cs_sles_t             *init_sles;

  /* Shared pointers */

  const cs_property_t   *pty_22;  /* Property related to the (2,2) block */
  const cs_real_t       *m21_val;
  const cs_adjacency_t  *m21_adj; /* Indexed list used to scan the unassembled
                                     m21 operator */

} cs_saddle_solver_context_simple_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Retrieve the number of saddle-point systems which have been added
 */
/*----------------------------------------------------------------------------*/

int
cs_saddle_solver_get_n_systems(void);

/*----------------------------------------------------------------------------*/
/*
 * Get a pointer to saddle-point solver from its id.
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_by_id(int id);

/*----------------------------------------------------------------------------*/
/*
 * Add a new solver for solving a saddle-point problem.
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_add(cs_lnum_t                n1_elts,
                     int                      n1_dofs_by_elt,
                     cs_lnum_t                n2_elts,
                     int                      n2_dofs_by_elt,
                     const cs_param_saddle_t *saddlep,
                     cs_cdo_system_helper_t  *sh,
                     cs_sles_t               *main_sles);

/*----------------------------------------------------------------------------*/
/*
 * Free all remaining structures related to saddle-point solvers.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_finalize(void);

/*----------------------------------------------------------------------------*/
/*
 * Free a cs_saddle_solver_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_free(cs_saddle_solver_t **p_solver);

/*----------------------------------------------------------------------------*/
/*
 * Free/reset only a part of a cs_saddle_solver_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_clean(cs_saddle_solver_t *solver);

/*----------------------------------------------------------------------------*/
/*
 * Update the current monitoring state with n_iter.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_update_monitoring(cs_saddle_solver_t *solver,
                                   unsigned            n_iter);

/*----------------------------------------------------------------------------*/
/*
 * Log the monitoring performance of all saddle-point systems
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_log_monitoring(void);

/*----------------------------------------------------------------------------*/
/*
 * Retrieve the lumped matrix the inverse of the diagonal of the (1,1)-block
 * matrix. The storage of a matrix is in a gather view and the resulting array
 * is in scatter view.
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_saddle_solver_m11_inv_lumped(cs_saddle_solver_t   *solver,
                                const cs_real_t       scaling_coef,
                                const cs_matrix_t    *m11,
                                const cs_range_set_t *b11_rset,
                                cs_sles_t            *xtra_sles,
                                int                  *n_iter);

/*----------------------------------------------------------------------------*/
/*
 * Compute the resulting vector of the operation m12*x2
 * The stride is equal to 3 for the operator m21 (unassembled) x2 corresponds
 * to a "scatter" view
 *
 * This is an update operation. Be careful that the resulting array has to be
 * initialized.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_vector(cs_lnum_t             n2_dofs,
                                     const cs_real_t      *x2,
                                     const cs_adjacency_t *m21_adj,
                                     const cs_real_t      *m21_val,
                                     cs_real_t            *m12x2);

/*----------------------------------------------------------------------------*/
/*
 * Compute the resulting vector of the operation m12*x2
 * The stride is equal to 1 for the operator m21 (unassembled) x2 corresponds
 * to a "scatter" view.
 *
 * This is an update operation. Be careful that the resulting array has to be
 * initialized.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_scalar(cs_lnum_t             n2_elts,
                                     const cs_real_t      *x2,
                                     const cs_adjacency_t *m21_adj,
                                     const cs_real_t      *m21_val,
                                     cs_real_t            *m12x2);

/*----------------------------------------------------------------------------*/
/*
 * Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 3 for the operator m21 operator
 *        x1 corresponds to a "scatter" view
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_vector(cs_lnum_t             n2_dofs,
                                     const cs_real_t      *x1,
                                     const cs_adjacency_t *m21_adj,
                                     const cs_real_t      *m21_val,
                                     cs_real_t            *m21x1);

/*----------------------------------------------------------------------------*/
/*
 * Compute the resulting vector of the operation m21*x1
 * The stride is equal to 1 for the operator m21 operator x1 corresponds to a
 * "scatter" view
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_scalar(cs_lnum_t             n2_dofs,
                                     const cs_real_t      *x1,
                                     const cs_adjacency_t *m21_adj,
                                     const cs_real_t      *m21_val,
                                     cs_real_t            *m21x1);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for an algorithm related to the
 * ALU algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_create(cs_saddle_solver_t *solver);

/*----------------------------------------------------------------------------*/
/*
 * Free main memory consuming part of the context structure associated to an
 * ALU algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_clean(cs_saddle_solver_context_alu_t *ctx);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure associated to an ALU algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_free(cs_saddle_solver_context_alu_t **p_ctx);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for a block preconditioner used
 * in combination with a Krylov solver such as MINRES or GCR.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_create(cs_lnum_t           b22_max_size,
                                          cs_saddle_solver_t *solver);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure for a block preconditioner used in combination
 * with a Krylov solver such as MINRES or GCR.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_clean
(
 cs_saddle_solver_context_block_pcd_t *ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure for a block preconditioner used in combination
 * with a Krylov solver such as MINRES or GCR.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_free
(
 cs_saddle_solver_context_block_pcd_t **p_ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for the GKB algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_create(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*
 * Free the main memory consuming part of the context structure associated to
 * the GKB algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_clean(cs_saddle_solver_context_gkb_t *ctx);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure associated to the GKB algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_free
(
 cs_saddle_solver_context_gkb_t **p_ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for the algorithm relying on the
 * Notay's algebraic transformation.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_notay_create(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for an algorithm related to the
 * Uzawa-CG algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_create(cs_lnum_t           b22_max_size,
                                         cs_saddle_solver_t *solver);

/*----------------------------------------------------------------------------*/
/*
 * Free main memory consuming part of the context structure associated to a
 * Uzawa-CG algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_clean
(
 cs_saddle_solver_context_uzawa_cg_t *ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure associated to a Uzawa-CG algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_free
(
 cs_saddle_solver_context_uzawa_cg_t **p_ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Create and initialize the context structure for an algorithm related to the
 * SIMPLE algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_create(cs_lnum_t           b22_max_size,
                                       cs_saddle_solver_t *solver);

/*----------------------------------------------------------------------------*/
/*
 * Free main memory consuming part of the context structure associated to a
 * SIMPLE algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_clean(cs_saddle_solver_context_simple_t *ctx);

/*----------------------------------------------------------------------------*/
/*
 * Free the context structure associated to a SIMPLE algorithm.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_free
(
 cs_saddle_solver_context_simple_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*
 * Apply the Augmented Lagrangian-Uzawa algorithm to a saddle point problem
 * (the system is stored in a hybrid way). The stride is equal to 1 for the
 * matrix (db_size[3] = 1) and the vector.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_alu_incr(cs_saddle_solver_t *solver,
                          cs_real_t          *x1,
                          cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the GCR algorithm to a saddle point problem (the system is stored in a
 * hybrid way).  The stride is equal to 1 for the matrix in the (1,1) block
 * (db_size[3] = 1).
 *
 * This algorithm is taken from 2010 Notay's paper:
 * "An aggregation-based algebraic multigrid method" ETNA (vol. 37)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gcr(cs_saddle_solver_t *solver,
                     cs_real_t          *x1,
                     cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the GKB algorithm to a saddle point problem (the system is stored in a
 * hybrid way).  The stride is equal to 1 for the matrix in the (1,1) block
 * (db_size[3] = 1).
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gkb_inhouse(cs_saddle_solver_t *solver,
                             cs_real_t          *x1,
                             cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply an (external) solver to solve a saddle point problem (the system is
 * stored in a monolithic way).
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_sles_full_system(cs_saddle_solver_t *solver,
                                  cs_real_t          *x1,
                                  cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the MINRES algorithm to a saddle point problem (the system is stored
 * in a hybrid way).  The stride is equal to 1 for the matrix (db_size[3] = 1)
 * and the vector.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_minres(cs_saddle_solver_t *solver,
                        cs_real_t          *x1,
                        cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the Notay's transformation algorithm to solve a saddle point problem
 * (the system is stored in a monolithic way).
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_notay(cs_saddle_solver_t *solver,
                       cs_real_t          *x1,
                       cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the Uzawa-CG algorithm to solve a saddle point problem (the system is
 * stored in a hybrid way). This algorithm is based on Koko's paper "Uzawa
 * conjugate gradient method for the Stokes problem: Matlab implementation with
 * P1-iso-P2/ P1 finite element".
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_uzawa_cg(cs_saddle_solver_t *solver,
                          cs_real_t          *x1,
                          cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the SIMPLE-like algorithm to solve a saddle point problem.
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_simple(cs_saddle_solver_t *solver,
                        cs_real_t          *x1,
                        cs_real_t          *x2);

/*----------------------------------------------------------------------------*/
/*
 * Apply the AFS (Algebraic Fractional Step) algorithm to solve a saddle point
 * problem cf. Bonelle & Ma (HAL preprint available online).
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_afs(cs_saddle_solver_t *solver,
                     cs_real_t          *x1,
                     cs_real_t          *x2);

/*----------------------------------------------------------------------------*/

#endif /* CS_SADDLE_SOLVER_H */
