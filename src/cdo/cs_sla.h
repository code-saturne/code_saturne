#ifndef __CS_SLA_H__
#define __CS_SLA_H__

/*============================================================================
 *  Sparse linear algebra routines
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Matrix flag */
#define  CS_SLA_MATRIX_SYM        (1 <<  0)  /*    1: symmetric */
#define  CS_SLA_MATRIX_SORTED     (1 <<  1)  /*    2: sorted */
#define  CS_SLA_MATRIX_LU         (1 <<  2)  /*    4: lower/upper separation */
#define  CS_SLA_MATRIX_SHARED     (1 <<  3)  /*    8: share pattern */

/*============================================================================
 * Type definitions for matrices
 *============================================================================*/

typedef enum {

  CS_SLA_MAT_NONE,
  CS_SLA_MAT_DEC,
  CS_SLA_MAT_CSR,
  CS_SLA_MAT_MSR,
  CS_SLA_MAT_N_TYPES

} cs_sla_matrix_type_t;

typedef struct {

  /* Stencil */
  int     stencil_min;
  int     stencil_max;
  double  stencil_mean;

  /* Fill-in */
  size_t  nnz;
  double  fillin;

} cs_sla_matrix_info_t;

typedef struct _sla_mat_prop_t cs_sla_mat_prop_t; /* Private definition */

typedef struct {

  cs_sla_matrix_type_t   type;
  cs_sla_mat_prop_t     *properties;
  int                    flag;       /* Symmetric, sorted, LU, shared... */

  int     stride;   /* Number of entries in "val" for each couple A(i,j) */
  int     n_rows;
  int     n_cols;

  cs_lnum_t   *idx;
  cs_lnum_t   *col;

  short int   *sgn;  /* Only for DEC type: -1 or 1 according to orientation */
  double      *val;  /* Only for CSR and MSR type */

  cs_lnum_t   *didx;  /* Diagonal index: used to flag diagonal index and to
                         separate lower and upper part */
  double      *diag;  /* Used in MSR type but can also be used for a
                         given preconditioning technique */

} cs_sla_matrix_t;

/*============================================================================
 * Type definitions for iterative solvers
 *============================================================================*/

typedef enum {

  CS_SLA_CODE_CONTINUE,   /* Need one more iteration */
  CS_SLA_CODE_CVG,        /* Convergence criterion is fullfilled */
  CS_SLA_CODE_CVG_EPS,    /* Convergence is reached to the machine precision
                             we can not go further */
  CS_SLA_CODE_MAX_ITER,   /* Max iteration is reached. Stop solver iterations */
  CS_SLA_CODE_DIVZERO,    /* Division by zero is encountered. Stop algorithm */
  CS_SLA_CODE_STAG,       /* No evolution in the residual norm since several
                             iterations. Stop algorithm */
  CS_SLA_CODE_DVG,        /* Residual norm is now greater than the initial
                             residual by a big factor */
  CS_SLA_CODE_STOP        /* Algorithm encountered a problem. Stop algorithm */

} cs_sla_code_t;

typedef struct {

  cs_sla_code_t    code;       // Convergence code
  int              iter;       // Current iteration
  double           residual;   // Current residual norm computed

} cs_sla_sumup_t;


/*============================================================================
 * Public function prototypes for SLA matrices
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write in binary format a matrix in CSR format, its righ hand side
 *         and the solution
 *
 * \param[in]     name      name of the output file
 * \param[in]     A         system to solve
 * \param[in]     rhs       right hand side
 * \param[in]     x         solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_bwrite(const char            *name,
              const cs_sla_matrix_t *m,
              const double          *rhs,
              const double          *sol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read from a binary file a matrix in CSR format, its righ hand side
 *         and the solution. Matrix must have a stride equal to 1.
 *
 * \param[in]     name      name of the output file
 * \param[inout]  p_mat     system to solve
 * \param[inout]  p_rhs     right hand side
 * \param[inout]  p_sol       solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_bread(const char        *name,
             cs_sla_matrix_t  **p_mat,
             double            *p_rhs[],
             double            *p_sol[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Change matrix representation from MSR to CSR
 *
 * \param[inout] a     matrix to transform
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_msr2csr(cs_sla_matrix_t  *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Change matrix representation from CSR to MSR
 *
 * \param[inout] a     matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_csr2msr(cs_sla_matrix_t  *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate its own pattern if shared
 *
 * \param[inout] a     matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_share2own(cs_sla_matrix_t  *a);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a new matrix resulting from the extraction of some listed
 *          rows and columns. The output is a new matrix packed (or zipped)
 *          w.r.t the initial matrix.
 *
 * \param[in]  n_final_rows   number of rows to extract
 * \param[in]  n_final_cols   number of columns to extract
 * \param[in]  init           init matrix to work with
 * \param[in]  row_z2i_ids    zipped -> initial ids for rows
 * \param[in]  col_i2z_ids    initial-> zipped ids for columns (-1 ==> remove)
 * \param[in]  keep_sym       true or false
 *
 * \return a pointer to the new (packed) matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_pack(cs_lnum_t                n_final_rows,
                   cs_lnum_t                n_final_cols,
                   const cs_sla_matrix_t   *init,
                   const cs_lnum_t         *row_z2i_ids,
                   const cs_lnum_t         *col_i2z_ids,
                   bool                     keep_sym);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build diagonal index
 *
 * \param[inout]  m   matrix to work with
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_diag_idx(cs_sla_matrix_t    *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the diagonal entries of a given matrix.
 *
 * \param[in]    m        matrix to work with
 * \param[inout] p_diag   pointer to diag array to define (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_get_diag(const cs_sla_matrix_t  *m,
                       double                 *p_diag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Sort each row by increasing colomn number
 *
 * \param[in]  mat     matrix to sort
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_sort(cs_sla_matrix_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a Lower/Upper pattern by sorting entries in each row by
 *          increasing number and by defining diagonal index
 *          Matrix should stored in CSR format and be square
 *
 * \param[inout]  m   matrix to work with
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_lu_pattern(cs_sla_matrix_t    *m);

/*----------------------------------------------------------------------------
 * Create a new matrix from a block description (not done for all matrix
 * combinations, see below).
 *
 *     |A | B|
 * M = |--|--|
 *     |C | D|
 *----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_block2mat(cs_sla_matrix_t         *A,
                 const cs_sla_matrix_t   *B,
                 const cs_sla_matrix_t   *C,
                 const cs_sla_matrix_t   *D,
                 bool                     sym);

/* Operations on matrices */
void cs_sla_matvec(const cs_sla_matrix_t *m,
                   const double           v[],
                   double                *inout[],
                   bool                   reset);

void cs_sla_amxby(double                  alpha,
                  const cs_sla_matrix_t  *m,
                  const double            x[],
                  double                  beta,
                  const double            y[],
                  double                 *inout[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Matrix block 2x2 multiply by a vector
 *
 *              | A  B | |X| = |F|= |AX + BY|
 *              | C  D | |Y|   |G|  |CX + DY|
 *
 * \param[in]    A     pointer to a cs_sla_matrix_t block (1,1)
 * \param[in]    B     pointer to a cs_sla_matrix_t block (1,2)
 * \param[in]    C     pointer to a cs_sla_matrix_t block (2,1)
 * \param[in]    D     pointer to a cs_sla_matrix_t block (2,2)
 * \param[in]    X     upper vector
 * \param[in]    Y     lower vector
 * \param[inout] F     upper result
 * \param[inout] G     lower result
 * \param[in]    reset reset before computation (true/false)
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matvec_block2(const cs_sla_matrix_t  *A,
                     const cs_sla_matrix_t  *B,
                     const cs_sla_matrix_t  *C,
                     const cs_sla_matrix_t  *D,
                     const double            X[],
                     const double            Y[],
                     double                 *F[],
                     double                 *G[],
                     bool                    reset);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two sparse matrices a and b. c = alpha*a + beta*b
 *
 * \param[in]  alpha   first coef.
 * \param[in]  a       first matrix to add
 * \param[in]  beta    second coef.
 * \param[in]  b       second matrix to add
 *
 * \return  a pointer to the resulting matrix
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_add(double                  alpha,
                  const cs_sla_matrix_t  *a,
                  double                  beta,
                  const cs_sla_matrix_t  *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main subroutine to multiply two sparse matrices a and b. c= a*b
 *
 * \param[in]  a     first matrix to multiply
 * \param[in]  b     second matrix to multiply
 *
 * \return  a pointer to the resulting matrix
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_multiply(const cs_sla_matrix_t   *a,
                       const cs_sla_matrix_t   *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Specific matrix multiplication. Compute Bt * beta * B + alpha * A
 *          alpha and beta are scalar
 *
 * \param[in] alpha     real coefficient
 * \param[in] a         square sym. matrix
 * \param[in] beta      real coefficient
 * \param[in] b         matrix (CSR or DEC)
 * \param[in] bt        adjoint matrix of b
 *
 * \return the matrix corresponding to this operation
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_combine(double                  alpha,
                      const cs_sla_matrix_t  *a,
                      double                  beta,
                      const cs_sla_matrix_t  *bt,
                      const cs_sla_matrix_t  *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the product C = At * Diag * A
 *
 * \param[in]     At   pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[in]     D    array standing for a diagonal operator
 * \param[in]     A    pointer to a cs_sla_matrix_t struct. (DEC type)
 * \param[inout]  w    work buffer
 *
 * \return a pointer to a new cs_sla_matrix_t
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t * cs_sla_multiply_AtDA(const cs_sla_matrix_t   *At,
                                       const double             D[],
                                       const cs_sla_matrix_t   *A,
                                       int                     *w);

int cs_sla_get_mat_spectral_radius(const cs_sla_matrix_t  *matrix,
                                   int                    iter_max,
                                   double                 epsilon,
                                   double                *lambda);

double cs_sla_get_matrix_norm(const cs_sla_matrix_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_sla_matrix_t structure
 *
 * \param[in]  n_rows     number of rows
 * \param[in]  n_cols     number of columns
 * \param[in]  stride     number of values related to each entry
 * \param[in]  type       kind of matrix
 * \param[in]  sym        true or false
 *
 * \return  a pointer to new allocated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_create(int                   n_rows,
                     int                   n_cols,
                     int                   stride,
                     cs_sla_matrix_type_t  type,
                     bool                  sym);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a cs_sla_matrix_t structure from an existing one
 *
 * \param[in]  ref      pointer to a reference matrix with the same pattern
 * \param[in]  type     type of the matrix to create
 * \param[in]  stride   number of values associated to each entry
 *
 * \return  a pointer to new allocated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_create_from_pattern(const cs_sla_matrix_t  *ref,
                                  cs_sla_matrix_type_t    type,
                                  int                     stride);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Create a new matrix structure from the copy of an existing one.
 *
 * \param[in] a       matrix to copy
 * \param[in] shared  true: only pointer are copied other allocated
 *
 * \return the matrix corresponding to this operation
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_copy(const cs_sla_matrix_t  *a,
                   bool                    shared);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Transpose a cs_sla_matrix_t structure
 *
 * \param[in]  mat     matrix to transpose
 *
 * \return  a pointer to new allocated matrix which is transposed to mat
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_transpose(const cs_sla_matrix_t  *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_sla_matrix_t structure
 *
 * \param[in]  m     matrix to free
 *
 *\return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_sla_matrix_free(cs_sla_matrix_t  *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Remove entries in a cs_sla_matrix_t structure below a given
 *          threshold. |a(i,j)| < eps * max|a(i,j)|
 *
 * \param[inout]  mat     matrix to clean
 * \param[in]     eps     value of the threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_clean(cs_sla_matrix_t   *m,
                    double             eps);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the number of non-zeros (nnz) elements in a matrix
 *
 * \param[in]  mat     matrix
 *
 * \return the number of nnz
 */
/*----------------------------------------------------------------------------*/

size_t
cs_sla_matrix_get_nnz(const cs_sla_matrix_t   *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute general information related to a cs_sla_matrix_t structure
 *
 * \param[in]  m           matrix to analyse
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_info_t
cs_sla_matrix_analyse(const cs_sla_matrix_t   *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Synthesis of a cs_sla_matrix_t structure
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_resume(const char              *name,
                     FILE                    *f,
                     const cs_sla_matrix_t   *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_sla_matrix_t structure
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_matrix_dump(const char              *name,
                   FILE                    *f,
                   const cs_sla_matrix_t   *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a cs_sla_matrix_t structure and its related right-hand side
 *
 * \param[in]  name        either name of the file if f is NULL or description
 * \param[in]  f           pointer to a FILE struct.
 * \param[in]  m           matrix to dump
 * \param[in]  rhs         right-hand side to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_system_dump(const char              *name,
                   FILE                    *f,
                   const cs_sla_matrix_t   *m,
                   const double            *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assemble a MSR matrix from local contributions
 *          --> We assume that the local matrices are symmetric
 *          --> We assume that the assembled matrix has its columns sorted
 *
 * \param[in]     loc     pointer to a local matrix
 * \param[inout]  ass     pointer to a cs_sla_matrix_t struct. collecting data
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_assemble_msr(const cs_toolbox_locmat_t  *loc,
                    cs_sla_matrix_t            *ass);

/*============================================================================
 * Public function prototypes for iterative solvers
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Evaluate the sizes of temporary buffers used in iterative solvers
 *
 * \param[in]  refsize        reference size to evaluate buffer dimension
 * \param[in]  info           information about the iterative solver to use
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_update_buffer_sizes(size_t                   refsize,
                                 const cs_param_itsol_t   info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate temporary buffers used in iterative solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_alloc_buffers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free temporary buffers useful for iterative solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_sla_itsol_free_buffers(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a sparse linear system A*x = rhs
 *         Matrix A can re-order during the solving process.
 *
 * \param[in]     itsol_info   set of parameters for the iterative solver
 * \param[in]     A            matrix to invert
 * \param[in]     rhs          right hand side
 * \param[inout]  x            initial guess (in) / solution (out)
 *
 * \return  a summary of the convergence monitoring
 */
/*----------------------------------------------------------------------------*/

cs_sla_sumup_t
cs_sla_solve(const cs_param_itsol_t    itsol_info,
             const cs_sla_matrix_t    *A,
             const double              rhs[],
             double                   *x);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name of a monitoring code
 *
 * \param[in] code     type of code
 *
 * \return the associated code name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_sla_get_code_name(cs_sla_code_t  code);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLA_H__ */
