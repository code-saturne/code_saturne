#ifndef __CS_SADDLE_SOLVER_H__
#define __CS_SADDLE_SOLVER_H__

/*============================================================================
 * Solvers for saddle-point systems arising from CDO discretizations
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_blas.h"
#include "cs_cdo_system.h"
#include "cs_iter_algo.h"
#include "cs_param_saddle.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * \param[in]      mat_val  values associated to the matrix operator
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

  const cs_param_saddle_t       *param;

  /* Structure to handle iterative algorithms */

  cs_iter_algo_t                *algo;

  /* Description of the saddle-point system to solve (shared) */

  cs_cdo_system_helper_t        *system_helper;

  bool                           do_setup;

  /* Main SLES structure associated to the saddle-point problem. According to
   * the type of saddle-point solver, additional SLES can be allocated, for
   * instance to compute an approximation of the Schur complement. */

  cs_sles_t                     *main_sles;

  /* Scatter viewpoint */

  cs_lnum_t                      n1_elts;
  int                            n1_dofs_by_elt;
  cs_lnum_t                      n2_elts;
  int                            n2_dofs_by_elt;

  cs_lnum_t                      n1_scatter_dofs;
  cs_lnum_t                      n2_scatter_dofs;

  /* Additional members according to the type of saddle-point solver. This
   * structure is cast on-the-fly */

  void                          *context;

  /* Monitoring */

  unsigned                       n_calls;
  unsigned                       n_iter_min;
  unsigned                       n_iter_max;
  unsigned                       n_iter_tot;

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
/*!
 * \brief Retrieve the number of saddle-point systems which have been added
 *
 * \return the current number of systems
 */
/*----------------------------------------------------------------------------*/

int
cs_saddle_solver_get_n_systems(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a pointer to saddle-point solver from its id
 *
 * \param[in] id  id of the saddle-point system
 *
 * \return a pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new solver for solving a saddle-point problem.
 *
 * \param[in] n1_elts         number of elements associated to the (1,1)-block
 * \param[in] n1_dofs_by_elt  number of DoFs by elements in the (1,1)-block
 * \param[in] n2_elts         number of elements associated to the (2,2)-block
 * \param[in] n2_dofs_by_elt  number of DoFs by elements in the (2,2)-block
 * \param[in] saddlep         set of parameters for the saddle-point solver
 * \param[in] sh              pointer to a system helper structure
 * \param[in] main_sles       pointer to the main SLES structure related to
 *                            this saddle-point problem
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_add(cs_lnum_t                 n1_elts,
                     int                       n1_dofs_by_elt,
                     cs_lnum_t                 n2_elts,
                     int                       n2_dofs_by_elt,
                     const cs_param_saddle_t  *saddlep,
                     cs_cdo_system_helper_t   *sh,
                     cs_sles_t                *main_sles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all remaining structures related to saddle-point solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_saddle_solver_t structure
 *
 * \param[in, out] p_solver  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_free(cs_saddle_solver_t  **p_solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free/reset only a part of a cs_saddle_solver_t structure
 *
 * \param[in, out] solver  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_clean(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the current monitoring state with n_iter
 *
 * \param[in, out] solver  pointer to a saddle solver structure
 * \param[in]      n_iter  number of iterations needed for a new call
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_update_monitoring(cs_saddle_solver_t  *solver,
                                   unsigned             n_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the monitoring performance of all saddle-point systems
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_log_monitoring(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the lumped matrix the inverse of the diagonal of the
 *        (1,1)-block matrix. The storage of a matrix is in a gather view and
 *        the resulting array is in scatter view.
 *
 * \param[in]      solver     solver for saddle-point problems
 * \param[in]      m11        matrix related to the (1,1) block
 * \param[in]      b11_rset   range set structure for the (1,1) block
 * \param[in, out] xtra_sles  pointer to an extra SLES structure
 * \param[out]     n_iter     number of iterations for this operation
 *
 * \return a pointer to the computed array (scatter view)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_saddle_solver_m11_inv_lumped(cs_saddle_solver_t     *solver,
                                const cs_matrix_t      *m11,
                                const cs_range_set_t   *b11_rset,
                                cs_sles_t              *xtra_sles,
                                int                    *n_iter);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 3 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view
 *
 *        This is an update operation. Be careful that the resulting array has
 *        to be initialized.
 *
 * \param[in]      n2_dofs  number of DoFs for x2
 * \param[in]      x2       array for the second set
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  array associated to the M21 operator (unassembled)
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_vector(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x2,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m12x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 1 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view.
 *
 *        This is an update operation. Be careful that the resulting array has
 *        to be initialized.
 *
 * \param[in]      n2_elts  number of elements for x2 (not DoFs)
 * \param[in]      x2       array for the second set
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  array associated to the M21 operator (unassembled)
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m12_multiply_scalar(cs_lnum_t              n2_elts,
                                     const cs_real_t       *x2,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m12x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 3 for the operator m21 operator
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      n2_dofs  number of (scatter) DoFs for (2,2)-block
 * \param[in]      x1       array for the first part
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  values associated to the M21 operator (unassembled)
 * \param[in, out] m21x1    resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_vector(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x1,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m21x1);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 1 for the operator m21 operator
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      n2_dofs  number of (scatter) DoFs for (2,2)-block
 * \param[in]      x1       array for the first part
 * \param[in]      m21_adj  adjacency related to the M21 operator
 * \param[in]      m21_val  values associated to the M21 operator (unassembled)
 * \param[in, out] m21x1    resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_m21_multiply_scalar(cs_lnum_t              n2_dofs,
                                     const cs_real_t       *x1,
                                     const cs_adjacency_t  *m21_adj,
                                     const cs_real_t       *m21_val,
                                     cs_real_t             *m21x1);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the ALU algorithm
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_create(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to an ALU algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_clean
(
 cs_saddle_solver_context_alu_t  *ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to an ALU algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_alu_free
(
 cs_saddle_solver_context_alu_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for a block preconditioner
 *        used in combination with a Krylov solver such as MINRES or GCR.
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_create(cs_lnum_t            b22_max_size,
                                          cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure for a block preconditioner used in
 *        combination with a Krylov solver such as MINRES or GCR.
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_clean
(
 cs_saddle_solver_context_block_pcd_t  *ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure for a block preconditioner used in
 *        combination with a Krylov solver such as MINRES or GCR.
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_block_pcd_free
(
 cs_saddle_solver_context_block_pcd_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for the GKB algorithm
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_create(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the main memory consuming part of the context structure
 *        associated to the GKB algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_clean(cs_saddle_solver_context_gkb_t  *ctx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to the GKB algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_gkb_free
(
 cs_saddle_solver_context_gkb_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for the algorithm relying
 *        on the Notay's algebraic transformation
 *
 * \param[in, out] solver  pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_notay_create(cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the Uzawa-CG algorithm
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_create(cs_lnum_t            b22_max_size,
                                         cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to a Uzawa-CG algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_clean
(
 cs_saddle_solver_context_uzawa_cg_t  *ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to a Uzawa-CG algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_uzawa_cg_free
(
 cs_saddle_solver_context_uzawa_cg_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize the context structure for an algorithm related
 *        to the SIMPLE algorithm
 *
 * \param[in]      b22_max_size  max. size for the second part of unknows
 * \param[in, out] solver        pointer to a saddle-point solver structure
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_create(cs_lnum_t            b22_max_size,
                                       cs_saddle_solver_t  *solver);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free main memory consuming part of the context structure associated
 *        to a SIMPLE algorithm
 *
 * \param[in, out] ctx  pointer to the context structure to clean
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_clean
(
 cs_saddle_solver_context_simple_t  *ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the context structure associated to a Uzawa-CG algorithm
 *
 * \param[in, out] p_ctx  double pointer to the context structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_context_simple_free
(
 cs_saddle_solver_context_simple_t  **p_ctx
);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Augmented Lagrangian-Uzawa algorithm to a saddle point
 *        problem (the system is stored in a hybrid way). The stride is equal
 *        to 1 for the matrix (db_size[3] = 1) and the vector.
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_alu_incr(cs_saddle_solver_t  *solver,
                          cs_real_t           *x1,
                          cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the GCR algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix in the (1,1) block
 *        (db_size[3] = 1).
 *
 *        This algorithm is taken from 2010 Notay's paper:
 *        "An aggregation-based algebraic multigrid method" ETNA (vol. 37)
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gcr(cs_saddle_solver_t  *solver,
                     cs_real_t           *x1,
                     cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the GKB algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix in the (1,1) block
 *        (db_size[3] = 1).
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_gkb_inhouse(cs_saddle_solver_t  *solver,
                             cs_real_t           *x1,
                             cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply an (external) solver to solve a saddle point problem (the
 *        system is stored in a monolithic way).
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_sles_full_system(cs_saddle_solver_t  *solver,
                                  cs_real_t           *x1,
                                  cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the MINRES algorithm to a saddle point problem (the system is
 *        stored in a hybrid way).
 *        The stride is equal to 1 for the matrix (db_size[3] = 1) and the
 *        vector.
 *
 * \param[in, out] solver  pointer to a saddle_point solver structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_minres(cs_saddle_solver_t  *solver,
                        cs_real_t           *x1,
                        cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Notay's transformation algorithm to solve a saddle point
 *        problem (the system is stored in a monolithic way).
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_notay(cs_saddle_solver_t  *solver,
                       cs_real_t           *x1,
                       cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the Uzawa-CG algorithm to solve a saddle point problem (the
 *        system is stored in a hybrid way). This algorithm is based on Koko's
 *        paper "Uzawa conjugate gradient method for the Stokes problem: Matlab
 *        implementation with P1-iso-P2/ P1 finite element"
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_uzawa_cg(cs_saddle_solver_t  *solver,
                          cs_real_t           *x1,
                          cs_real_t           *x2);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the SIMPLE-like algorithm to solve a saddle point problem
 *
 * \param[in, out] solver  pointer to a cs_saddle_solver_t structure
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_simple(cs_saddle_solver_t  *solver,
                        cs_real_t           *x1,
                        cs_real_t           *x2);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SADDLE_SOLVER_H__ */
