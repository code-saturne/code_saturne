/*============================================================================
 * Set of operations to handle Small Dense Matrices (SDM)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo.h"
#include "cs_blas.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_sort.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sdm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static definitions
 *============================================================================*/

static const char  _msg_small_p[] =
  " %s: Very small or null pivot.\n Stop inversion.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic way to allocate a cs_sdm_t structure
 *
 * \param[in]  flag         metadata related to a cs_sdm_t structure
 * \param[in]  n_max_rows   max number of rows
 * \param[in]  n_max_cols   max number of columns
 * \param[in]  n_max_rows   max number of rows
 * \param[in]  n_max_cols   max number of columns
 *
 *
 * \return  a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sdm_t *
_create_sdm(cs_flag_t   flag,
            int         n_max_rows,
            int         n_max_cols)
{
  cs_sdm_t  *mat = NULL;

  BFT_MALLOC(mat, 1, cs_sdm_t);

  mat->flag = flag;
  mat->n_max_rows = n_max_rows;
  mat->n_max_cols = n_max_cols;
  mat->n_rows = n_max_rows;
  mat->n_cols = n_max_cols;

  BFT_MALLOC(mat->val, mat->n_max_rows*mat->n_max_cols, cs_real_t);
  memset(mat->val, 0, sizeof(cs_real_t)*mat->n_max_rows*mat->n_max_cols);

  if (flag & CS_SDM_BY_BLOCK) {

    cs_sdm_block_t  *bd = NULL;

    BFT_MALLOC(bd, 1, cs_sdm_block_t);
    bd->blocks = NULL;
    bd->n_max_blocks_by_row = bd->n_max_blocks_by_col = 0;
    bd->n_row_blocks = bd->n_col_blocks = 0;
    mat->block_desc = bd;

  }
  else
    mat->block_desc = NULL;

  return mat;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure
 *          Most generic function to create a cs_sdm_t structure without block
 *          description
 *
 * \param[in]  flag         metadata related to a cs_sdm_t structure
 * \param[in]  n_max_rows   max number of rows
 * \param[in]  n_max_cols   max number of columns
 *
 * \return  a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_create(cs_flag_t   flag,
              int         n_max_rows,
              int         n_max_cols)
{
  return _create_sdm(flag, n_max_rows, n_max_cols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure
 *          Case of a square matrix
 *
 * \param[in]  n_max_rows   max number of rows
 *
 * \return  a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_square_create(int   n_max_rows)
{
  return _create_sdm(0, n_max_rows, n_max_rows);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate a cs_sdm_t structure and initialized it with the
 *          copy of the matrix m in input
 *
 * \param[in]     m     pointer to a cs_sdm_t structure to copy
 *
 * \return  a pointer to a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_create_copy(const cs_sdm_t   *m)
{
  cs_sdm_t  *c = _create_sdm(m->flag, m->n_max_rows, m->n_max_cols);

  c->n_rows = m->n_rows;
  c->n_cols = m->n_cols;
  memcpy(c->val, m->val, sizeof(cs_real_t)*m->n_rows*m->n_cols);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix by adding the given matrix with its transpose.
 *          Keep the transposed matrix for a future use.
 *
 * \param[in] mat   local matrix to transpose
 *
 * \return a pointer to the new allocated transposed matrix
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_create_transpose(cs_sdm_t  *mat)
{
  /* Sanity check */
  assert(mat != NULL);

  cs_sdm_t  *tr = _create_sdm(mat->flag, mat->n_max_cols, mat->n_max_rows);

  tr->n_rows = mat->n_cols;
  tr->n_cols = mat->n_rows;

  for (short int i = 0; i < mat->n_rows; i++) {

    const cs_real_t  *mval_i = mat->val + i*mat->n_cols;
    for (short int j = 0; j < mat->n_cols; j++)
      tr->val[j*tr->n_cols+i] = mval_i[j];

  }

  return tr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure
 *
 * \param[in]  n_max_blocks_by_row    max number of blocks in a row
 * \param[in]  n_max_blocks_by_col    max number of blocks in a column
 * \param[in]  max_row_block_sizes    max number of rows by block in a column
 * \param[in]  max_col_block_sizes    max number of columns by block in a row
 *
 * \return  a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_block_create(short int          n_max_blocks_by_row,
                    short int          n_max_blocks_by_col,
                    const short int    max_row_block_sizes[],
                    const short int    max_col_block_sizes[])
{
  cs_sdm_t  *m = NULL;

  if (n_max_blocks_by_row < 1 || n_max_blocks_by_col < 1)
    return m;

  int  row_size = 0, col_size = 0;
  for (int i = 0; i < n_max_blocks_by_row; i++)
    row_size += max_row_block_sizes[i];
  for (int j = 0; j < n_max_blocks_by_col; j++)
    col_size += max_col_block_sizes[j];

  m = _create_sdm(CS_SDM_BY_BLOCK, row_size, col_size);

  /* Define the block description */
  m->block_desc->n_max_blocks_by_row = n_max_blocks_by_row;
  m->block_desc->n_max_blocks_by_col = n_max_blocks_by_col;
  m->block_desc->n_row_blocks = n_max_blocks_by_row;
  m->block_desc->n_col_blocks = n_max_blocks_by_col;
  BFT_MALLOC(m->block_desc->blocks,
             n_max_blocks_by_row*n_max_blocks_by_col, cs_sdm_t);

  cs_real_t  *p_val = m->val;
  int  shift = 0;
  for (int i = 0; i < n_max_blocks_by_row; i++) {
    const short int  n_rows_i = max_row_block_sizes[i];
    for (int j = 0; j < n_max_blocks_by_col; j++) {
      const short int  n_cols_j = max_col_block_sizes[j];

      /* Set the block (i,j) */
      cs_sdm_t  *b_ij = m->block_desc->blocks + shift;
      int  _size = n_rows_i*n_cols_j;

      cs_sdm_map_array(n_rows_i, n_cols_j, b_ij, p_val);
      shift++;
      p_val += _size;

    }
  }

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_sdm_t structure
 *
 * \param[in]  mat    pointer to a cs_sdm_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_free(cs_sdm_t  *mat)
{
  if (mat == NULL)
    return mat;

  if ((mat->flag & CS_SDM_SHARED_VAL) == 0)
    BFT_FREE(mat->val);

  if (mat->flag & CS_SDM_BY_BLOCK) {
    BFT_FREE(mat->block_desc->blocks);
    BFT_FREE(mat->block_desc);
  }

  BFT_FREE(mat);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the pattern of cs_sdm_t structure defined by block
 *          The matrix should have been allocated before calling this function
 *
 * \param[in, out] m
 * \param[in]      n_blocks_by_row    number of blocks in a row
 * \param[in]      n_blocks_by_col    number of blocks in a column
 * \param[in]      row_block_sizes    number of rows by block in a column
 * \param[in]      col_block_sizes    number of columns by block in a row
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_init(cs_sdm_t          *m,
                  short int          n_blocks_by_row,
                  short int          n_blocks_by_col,
                  const short int    row_block_sizes[],
                  const short int    col_block_sizes[])
{
  /* Sanity checks */
  assert(m != NULL && row_block_sizes != NULL && col_block_sizes != NULL);
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);

  cs_sdm_block_t  *bd = m->block_desc;

  assert(n_blocks_by_row <= bd->n_max_blocks_by_row);
  assert(n_blocks_by_col <= bd->n_max_blocks_by_col);
  bd->n_row_blocks = n_blocks_by_row;
  bd->n_col_blocks = n_blocks_by_col;
  m->n_rows = 0;
  for (int i = 0; i < n_blocks_by_row; i++)
    m->n_rows += row_block_sizes[i];
  m->n_cols = 0;
  for (int i = 0; i < n_blocks_by_col; i++)
    m->n_cols += col_block_sizes[i];

  memset(m->val, 0, m->n_rows*m->n_cols*sizeof(cs_real_t));

  cs_real_t  *p_val = m->val;
  int  shift = 0;
  for (int i = 0; i < bd->n_row_blocks; i++) {
    const short int  n_rows_i = row_block_sizes[i];
    for (int j = 0; j < bd->n_col_blocks; j++) {
      const short int  n_cols_j = col_block_sizes[j];

      /* Set the block (i,j) */
      cs_sdm_t  *b_ij = bd->blocks + shift;

      cs_sdm_map_array(n_rows_i, n_cols_j, b_ij, p_val);
      p_val += n_rows_i*n_cols_j;
      shift++;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a local dense matrix-product c = a*b
 *          c has been previously allocated
 *
 * \param[in]      a     local dense matrix to use
 * \param[in]      b     local dense matrix to use
 * \param[in, out] c     result of the local matrix-product
 *                       is updated
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_multiply(const cs_sdm_t   *a,
                const cs_sdm_t   *b,
                cs_sdm_t         *c)
{
  /* Sanity checks */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->n_cols == b->n_rows &&
         a->n_rows == c->n_rows &&
         c->n_cols == b->n_cols);

  const cs_real_t *bv = b->val;

  for (short int i = 0; i < a->n_rows; i++) {

    const cs_real_t *av_i = a->val + i*a->n_cols;
    cs_real_t  *cv_i = c->val + i*b->n_cols;

    for (short int j = 0; j < b->n_cols; j++){
      cs_real_t p = 0.0;
      for (short int k = 0; k < a->n_cols; k++)
        p += av_i[k] * bv[k*b->n_cols + j];
      cv_i[j] += p;

    } /* Loop on b columns */
  } /* Loop on a rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Generic version: all compatible sizes
 *
 * \param[in]      a    local matrix to use
 * \param[in]      b    local matrix to use
 * \param[in, out] c    result of the local matrix-product (already allocated)
 *                      is updated.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_multiply_rowrow(const cs_sdm_t   *a,
                       const cs_sdm_t   *b,
                       cs_sdm_t         *c)
{
  /* Sanity check */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->n_cols == b->n_cols &&
         a->n_rows == c->n_rows &&
         c->n_cols == b->n_rows);

  for (short int i = 0; i < a->n_rows; i++) {

    const cs_real_t  *av_i = a->val + i*a->n_cols;

    cs_real_t  *cv_i = c->val + i*b->n_rows;

    for (short int j = 0; j < b->n_rows; j++) {

      const cs_real_t  *bv_j = b->val + j * b->n_cols;

      cs_real_t  dp = 0;
      for (short int k = 0; k < a->n_cols; k++)
        dp += av_i[k] * bv_j[k];
      cv_i[j] += dp;

    } /* Loop on b rows */
  } /* Loop on a rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Generic version: all compatible sizes
 *           Result is known to be symmetric.
 *
 * \param[in]      a    local matrix to use
 * \param[in]      b    local matrix to use
 * \param[in, out] c    result of the local matrix-product (already allocated)
 *                      is updated.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_multiply_rowrow_sym(const cs_sdm_t   *a,
                           const cs_sdm_t   *b,
                           cs_sdm_t         *c)
{
  /* Sanity check */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->n_cols == b->n_cols &&
         a->n_rows == c->n_rows &&
         c->n_cols == b->n_rows);

  for (short int i = 0; i < a->n_rows; i++) {

    const cs_real_t  *av_i = a->val + i*a->n_cols;

    cs_real_t  *cv_i = c->val + i*b->n_rows;

    for (short int j = i; j < b->n_rows; j++) {

      const cs_real_t  *bv_j = b->val + j * b->n_cols;

      cs_real_t  dp = 0;
      for (short int k = 0; k < a->n_cols; k++)
        dp += av_i[k] * bv_j[k];
      cv_i[j] += dp;

      if (j > i)
        (c->val + j*b->n_rows)[i] += dp;

    } /* Loop on b rows */
  } /* Loop on a rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Case of matrices defined by block.
 *
 * \param[in]      a    local matrix to use
 * \param[in]      b    local matrix to use
 * \param[in, out] c    result of the local matrix-product (already allocated)
 *                      is updated
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_multiply_rowrow(const cs_sdm_t   *a,
                             const cs_sdm_t   *b,
                             cs_sdm_t         *c)
{
  /* Sanity check */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->flag & CS_SDM_BY_BLOCK);
  assert(b->flag & CS_SDM_BY_BLOCK);
  assert(c->flag & CS_SDM_BY_BLOCK);
  assert(a->n_cols == b->n_cols &&
         a->n_rows == c->n_rows &&
         c->n_cols == b->n_rows);

  const cs_sdm_block_t  *a_desc = a->block_desc;
  const cs_sdm_block_t  *b_desc = b->block_desc;
  const cs_sdm_block_t  *c_desc = c->block_desc;

  assert(a_desc->n_col_blocks == b_desc->n_col_blocks &&
         a_desc->n_row_blocks == c_desc->n_row_blocks &&
         c_desc->n_col_blocks == b_desc->n_row_blocks);

  for (short int i = 0; i < a_desc->n_row_blocks; i++) {
    for (short int j = 0; j < b_desc->n_row_blocks; j++) {

      cs_sdm_t  *cIJ = cs_sdm_get_block(c, i, j);

      for (short int k = 0; k < a_desc->n_col_blocks; k++) {

        cs_sdm_t  *aIK = cs_sdm_get_block(a, i, k);
        cs_sdm_t  *bJK = cs_sdm_get_block(b, j, k);

        cs_sdm_multiply_rowrow(aIK, bJK, cIJ);

      } /* Loop on common blocks between a and b */
    } /* Loop on b row blocks */
  } /* Loop on a row blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Case of matrices defined by block.
 *           Results is known to be symmetric.
 *
 * \param[in]      a    local matrix to use
 * \param[in]      b    local matrix to use
 * \param[in, out] c    result of the local matrix-product (already allocated)
 *                      is updated
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_multiply_rowrow_sym(const cs_sdm_t   *a,
                                 const cs_sdm_t   *b,
                                 cs_sdm_t         *c)
{
  /* Sanity check */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->flag & CS_SDM_BY_BLOCK);
  assert(b->flag & CS_SDM_BY_BLOCK);
  assert(c->flag & CS_SDM_BY_BLOCK);
  assert(a->n_cols == b->n_cols &&
         a->n_rows == c->n_rows &&
         c->n_cols == b->n_rows);

  const cs_sdm_block_t  *a_desc = a->block_desc;
  const cs_sdm_block_t  *b_desc = b->block_desc;
  const cs_sdm_block_t  *c_desc = c->block_desc;

  assert(a_desc->n_col_blocks == b_desc->n_col_blocks &&
         a_desc->n_row_blocks == c_desc->n_row_blocks &&
         c_desc->n_col_blocks == b_desc->n_row_blocks);

  for (short int i = 0; i < a_desc->n_row_blocks; i++) {

    for (short int j = i; j < b_desc->n_row_blocks; j++) {

      cs_sdm_t  *cIJ = cs_sdm_get_block(c, i, j);

      for (short int k = 0; k < a_desc->n_col_blocks; k++) {

        cs_sdm_t  *aIK = cs_sdm_get_block(a, i, k);
        cs_sdm_t  *bJK = cs_sdm_get_block(b, j, k);

        cs_sdm_multiply_rowrow_sym(aIK, bJK, cIJ);

      } /* Loop on common blocks between a and b */

      if (j > i)
        cs_sdm_transpose_and_update(cIJ, cs_sdm_get_block(c, j, i));

    } /* Loop on b row blocks */
  } /* Loop on a row blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a small square matrix
 *          mv has been previously allocated
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] mv result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_matvec(const cs_sdm_t    *mat,
                     const cs_real_t   *vec,
                     cs_real_t         *mv)
{
  /* Sanity checks */
  assert(mat != NULL && vec != NULL && mv != NULL);
  assert(mat->n_rows == mat->n_cols);

  const int  n = mat->n_rows;

  /* Initialize mv */
  const cs_real_t  v = vec[0];
  for (short int i = 0; i < n; i++)
    mv[i] = v*mat->val[i*n];

  /* Increment mv */
  for (short int i = 0; i < n; i++) {
    cs_real_t *m_i = mat->val + i*n;
    for (short int j = 1; j < n; j++)
      mv[i] +=  m_i[j] * vec[j];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a rectangular matrix
 *          mv has been previously allocated
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use (size = mat->n_cols)
 * \param[in, out] mv     result of the operation (size = mat->n_rows)
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_matvec(const cs_sdm_t    *mat,
              const cs_real_t   *vec,
              cs_real_t         *mv)
{
  /* Sanity checks */
  assert(mat != NULL && vec != NULL && mv != NULL);

  if (mat->n_rows == mat->n_cols) {
    cs_sdm_square_matvec(mat, vec, mv);
    return;
  }

  const short int  nr = mat->n_rows;
  const short int  nc = mat->n_cols;

  /* Initialize mv with the first column */
  const cs_real_t  v = vec[0];
  for (short int i = 0; i < nr; i++)
    mv[i] = v*mat->val[i*nc];

  /* Increment mv */
  for (short int i = 0; i < nr; i++) {
    cs_real_t *m_i = mat->val + i*nc;
    for (short int j = 1; j < nc; j++)
      mv[i] +=  m_i[j] * vec[j];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two matrices defined by block: loc += add
 *
 * \param[in, out] mat   local matrix storing the result
 * \param[in]      add   values to add to mat
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_add(cs_sdm_t        *mat,
                 const cs_sdm_t  *add)
{
  if (mat == NULL || add == NULL)
    return;

  const cs_sdm_block_t  *add_desc = add->block_desc;
  const cs_sdm_block_t  *mat_desc = mat->block_desc;

  assert(add_desc != NULL && mat_desc != NULL);
  assert(add_desc->n_row_blocks == mat_desc->n_row_blocks);
  assert(add_desc->n_col_blocks == mat_desc->n_col_blocks);

  for (short int bi = 0; bi < mat_desc->n_row_blocks; bi++) {
    for (short int bj = 0; bj < mat_desc->n_col_blocks; bj++) {

      cs_sdm_t  *mat_ij = cs_sdm_get_block(mat, bi, bj);
      const cs_sdm_t  *add_ij = cs_sdm_get_block(add, bi, bj);

      cs_sdm_add(mat_ij, add_ij);

    } /* Loop on column blocks */
  } /* Loop on row blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two small dense matrices: loc += add
 *
 * \param[in, out] mat   local matrix storing the result
 * \param[in]      add   values to add to mat
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_add(cs_sdm_t        *mat,
           const cs_sdm_t  *add)
{
  /* Sanity checks */
  assert(mat != NULL && add != NULL);
  assert(mat->n_rows == add->n_rows);
  assert(mat->n_cols == add->n_cols);

  for (int i = 0; i < mat->n_rows*mat->n_cols; i++)
    mat->val[i] += add->val[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two small dense matrices: loc += alpha*add
 *
 * \param[in, out] mat     local matrix storing the result
 * \param[in]      alpha   multiplicative coefficient
 * \param[in]      add     values to add to mat
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_add_mult(cs_sdm_t        *mat,
                cs_real_t        alpha,
                const cs_sdm_t  *add)
{
  /* Sanity checks */
  assert(mat != NULL && add != NULL);
  assert(mat->n_rows == add->n_rows);
  assert(mat->n_cols == add->n_cols);

  if (fabs(alpha) < FLT_MIN)
    return;

  for (int i = 0; i < mat->n_rows*mat->n_cols; i++)
    mat->val[i] += alpha * add->val[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix by adding the given matrix with its transpose.
 *          Keep the transposed matrix for a future use.
 *
 * \param[in, out] mat   local matrix to transpose and add
 * \param[in, out] tr    transposed of the local matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_add_transpose(cs_sdm_t  *mat,
                            cs_sdm_t  *tr)
{
  /* Sanity check */
  assert(mat != NULL && tr != NULL && tr->n_max_rows == mat->n_max_rows);
  assert(mat->n_rows == mat->n_cols);

  if (mat->n_rows < 1 || mat->n_cols < 1)
    return;

  tr->n_rows = mat->n_cols;
  tr->n_cols = mat->n_rows;

  for (short int i = 0; i < mat->n_rows; i++) {

    const int  ii = i*mat->n_cols + i;
    tr->val[ii] = mat->val[ii];
    mat->val[ii] *= 2;

    for (short int j = i+1; j < mat->n_cols; j++) {

      const int  ij = i*mat->n_cols + j;
      const int  ji = j*mat->n_cols + i;

      tr->val[ji] = mat->val[ij];
      tr->val[ij] = mat->val[ji];
      mat->val[ij] += tr->val[ij];
      mat->val[ji] += tr->val[ji];

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given matrix into its anti-symmetric part
 *
 * \param[in, out] mat   small dense matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_asymm(cs_sdm_t   *mat)
{
  /* Sanity check */
  assert(mat != NULL);
  assert(mat->n_rows == mat->n_cols);

  if (mat->n_rows < 1)
    return;

  for (short int i = 0; i < mat->n_rows; i++) {

    cs_real_t  *mi = mat->val + i*mat->n_cols;

    mi[i] = 0;

    for (short int j = i+1; j < mat->n_cols; j++) {

      int  ji = j*mat->n_rows + i;

      mi[j] = 0.5*(mi[j] - mat->val[ji]);
      mat->val[ji] = mi[j];

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 3x3 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \Note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[6])
{
  /* Sanity check */
  assert(m != NULL && facto != NULL);
  assert(m->n_cols == m->n_rows && m->n_cols == 3);

  // j=0: first row
  const cs_real_t  d00 = m->val[0];
  if (fabs(d00) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);

  facto[0] = 1. / d00;
  const cs_real_t  l10 = facto[1] = m->val[1] * facto[0];
  const cs_real_t  l20 = facto[3] = m->val[2] * facto[0];

  // j=1: second row
  const cs_real_t  d11 = m->val[4] - l10*l10 * d00;
  if (fabs(d11) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[2] = 1. / d11;
  const cs_real_t l21 = facto[4] = (m->val[5] - l20*d00*l10) * facto[2];

  // j=2: third row
  const cs_real_t  d22 = m->val[8] - l20*l20*d00 - l21*l21*d11;
  if (fabs(d22) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[5] = 1. / d22;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \Note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_44_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[10])
{
  /* Sanity check */
  assert(m != NULL && facto != NULL);
  assert(m->n_cols == m->n_rows && m->n_cols == 4);

  // j=0: first row
  const cs_real_t  d00 = m->val[0];
  if (fabs(d00) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);

  facto[0] = 1. / d00;
  const cs_real_t  l10 = facto[1] = m->val[1] * facto[0];
  const cs_real_t  l20 = facto[3] = m->val[2] * facto[0];
  const cs_real_t  l30 = facto[6] = m->val[3] * facto[0];

  // j=1: second row
  const cs_real_t  d11 = m->val[5] - l10*l10 * d00;
  if (fabs(d11) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[2] = 1. / d11;
  const cs_real_t  l21 = facto[4] = (m->val[6] - l20*d00*l10) * facto[2];
  const cs_real_t  l31 = facto[7] = (m->val[7] - l30*d00*l10) * facto[2];

  // j=2: third row
  const cs_real_t  d22 = m->val[10] - l20*l20*d00 - l21*l21*d11;
  if (fabs(d22) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[5] = 1. / d22;
  const cs_real_t  l32 = facto[8] =
    (m->val[11] - l30*d00*l20 - l31*d11*l21) * facto[5];

  // j=3: row 4
  const cs_real_t  d33 = m->val[15] - l30*l30*d00 - l31*l31*d11 - l32*l32*d22;
  if (fabs(d33) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[9] = 1. / d33;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 6x6 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \Note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_66_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[21])
{
  /* Sanity check */
  assert(m != NULL && facto != NULL);
  assert(m->n_cols == m->n_rows && m->n_cols == 6);

  // j=0: first row
  const cs_real_t  d00 = m->val[0];
  if (fabs(d00) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);

  facto[0] = 1. / d00;
  const cs_real_t  l10 = facto[ 1] = m->val[1] * facto[0];
  const cs_real_t  l20 = facto[ 3] = m->val[2] * facto[0];
  const cs_real_t  l30 = facto[ 6] = m->val[3] * facto[0];
  const cs_real_t  l40 = facto[10] = m->val[4] * facto[0];
  const cs_real_t  l50 = facto[15] = m->val[5] * facto[0];

  // j=1: second row
  const cs_real_t  d11 = m->val[7] - l10*l10 * d00;
  if (fabs(d11) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[2] = 1. / d11;
  const cs_real_t  d0l10 = d00*l10;
  const cs_real_t  l21 = facto[ 4] = (m->val[ 8] - l20*d0l10) * facto[2];
  const cs_real_t  l31 = facto[ 7] = (m->val[ 9] - l30*d0l10) * facto[2];
  const cs_real_t  l41 = facto[11] = (m->val[10] - l40*d0l10) * facto[2];
  const cs_real_t  l51 = facto[16] = (m->val[11] - l50*d0l10) * facto[2];

  // j=2: third row
  const cs_real_t  d22 = m->val[14] - l20*l20*d00 - l21*l21*d11;
  if (fabs(d22) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[5] = 1. / d22;
  const cs_real_t  d1l21 = d11*l21, d0l20 = d00*l20;
  const cs_real_t  l32 = facto[ 8] =
    (m->val[15] - l30*d0l20 - l31*d1l21) * facto[5];
  const cs_real_t  l42 = facto[12] =
    (m->val[16] - l30*d0l20 - l31*d1l21) * facto[5];
  const cs_real_t  l52 = facto[17] =
    (m->val[17] - l30*d0l20 - l31*d1l21) * facto[5];

  // j=3: row 4
  const cs_real_t  d33 = m->val[21] - l30*l30*d00 - l31*l31*d11 - l32*l32*d22;
  if (fabs(d33) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[9] = 1. / d33;
  const cs_real_t  d1l31 = d11*l31, d0l30 = d00*l30, d2l32 = d22*l32;
  const cs_real_t  l43 = facto[13] =
    (m->val[22] - l40*d0l30 - l41*d1l31 - l42*d2l32) * facto[9];
  const cs_real_t  l53 = facto[18] =
    (m->val[23] - l50*d0l30 - l51*d1l31 - l52*d2l32) * facto[9];

  // j=4: row 5
  const cs_real_t  d44 =
    m->val[28] - l40*l40*d00 - l41*l41*d11 - l42*l42*d22 - l43*l43*d33;
  if (fabs(d44) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[14] = 1. / d44;
  const cs_real_t  l54 = facto[19] = facto[14] *
    (m->val[29] - l50*d00*l40 - l51*d11*l41 - l52*d22*l42 - l53*d33*l43);

  // j=5: row 6
  const cs_real_t  d55 = m->val[35]
    - l50*l50*d00 - l51*l51*d11 - l52*l52*d22 - l53*l53*d33 - l54*l54*d44;
  if (fabs(d55) < cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
  facto[20] = 1. / d55;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 * \param[in, out] dkk      store temporary the diagonal (size = n_rows)
 *
 * \Note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_ldlt_compute(const cs_sdm_t     *m,
                    cs_real_t          *facto,
                    cs_real_t          *dkk)
{
  /* Sanity checks */
  assert(m != NULL && facto != NULL);
  assert(m->n_cols == m->n_rows);

  const short int n = m->n_cols;

  if (n == 1) {
    facto[0] = 1. / m->val[0];
    return;
  }

  int  rowj_idx = 0;

  /* Factorization (column-major algorithm) */
  for (short int j = 0; j < n; j++) {

    rowj_idx += j;
    const int  djj_idx = rowj_idx + j;

    switch (j) {

    case 0:  /* Optimization for the first colum */
      {
        dkk[0] = m->val[0]; // d00

        if (fabs(dkk[0]) < cs_math_zero_threshold)
          bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);
        const cs_real_t  inv_d00 = facto[0] = 1. / dkk[0];

        // l_i0 = a_i0 / d_00
        short int rowi_idx = rowj_idx;
        const cs_real_t  *a_0 = m->val;  /* a_ij = a_ji */
        for (short int i = j+1; i < n; i++) { /* Loop on rows */

          rowi_idx += i;
          cs_real_t  *l_i = facto + rowi_idx;
          l_i[0] = a_0[i] * inv_d00;

        }

      }
      break;

    case 1:  /* Optimization for the second colum */
      {
        // d_11 = a_11 - l_10^2 * d_00
        cs_real_t  *l_1 = facto + rowj_idx;

        const cs_real_t  d11 = dkk[1] = m->val[n+1] - l_1[0]*l_1[0]*dkk[0];
        if (fabs(d11) < cs_math_zero_threshold)
          bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);

        const cs_real_t  inv_d11 = facto[djj_idx] = 1. / d11;

        // l_i1 = (a_i1 - l_i0 * d_00 * l_10 ) / d_11
        short int rowi_idx = rowj_idx;
        const cs_real_t  *a_1 = m->val + n;  /* a_i1 = a_1i */
        for (short int i = 2; i < n; i++) { /* Loop on rows */

          rowi_idx += i;
          cs_real_t  *l_i = facto + rowi_idx;
          l_i[1] = (a_1[i] - l_i[0] *  dkk[0] * l_1[0]) * inv_d11;

        }

      }
      break;

    default:
      {
        // d_jj = a_jj - \sum_{k=0}^{j-1} l_jk^2 * d_kk
        cs_real_t  *l_j = facto + rowj_idx;

        cs_real_t  sum = 0.;
        for (short int k = 0; k < j; k++)
          sum += l_j[k]*l_j[k] * dkk[k];
        const cs_real_t  djj = dkk[j] = m->val[j*n+j] - sum;

        if (fabs(djj) < cs_math_zero_threshold)
          bft_error(__FILE__, __LINE__, 0, _msg_small_p, __func__);

        const cs_real_t  inv_djj = facto[djj_idx] = 1. / djj;

        // l_ij = (a_ij - \sum_{k=1}^{j-1} l_ik * d_kk * l_jk ) / d_jj
        short int rowi_idx = rowj_idx;
        const cs_real_t  *a_j = m->val + j*n;  /* a_ij = a_ji */
        for (short int i = j+1; i < n; i++) { /* Loop on rows */

          rowi_idx += i;
          cs_real_t  *l_i = facto + rowi_idx;
          sum = 0.;
          for (short int k = 0; k < j; k++)
            sum += l_i[k] *  dkk[k] * l_j[k];
          l_i[j] = (a_j[i] - sum) * inv_djj;

        }

      }
      break;
    } /* End of switch */

  } /* Loop on column j */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a 3x3 matrix with a modified Cholesky decomposition (L.D.L^T)
 *         The solution should be already allocated
 * Ref. http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      facto   vector of the coefficients of the decomposition
 * \param[in]      rhs     right-hand side
 * \param[in,out]  sol     solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_ldlt_solve(const cs_real_t    facto[6],
                     const cs_real_t    rhs[3],
                     cs_real_t          sol[3])
{
  /* Sanity check */
  assert(facto != NULL && rhs != NULL && sol != NULL);

  sol[0] = rhs[0];
  sol[1] = rhs[1] - sol[0]*facto[1];
  sol[2] =(rhs[2] - sol[0]*facto[3] - sol[1]*facto[4]) * facto[5];

  sol[1] = sol[1] * facto[2] - facto[4]*sol[2];
  sol[0] = sol[0] * facto[0] - facto[1]*sol[1] - facto[3]*sol[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a 4x4 matrix with a modified Cholesky decomposition (L.D.L^T)
 *         The solution should be already allocated
 * Ref. http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      facto   vector of the coefficients of the decomposition
 * \param[in]      rhs     right-hand side
 * \param[in,out]  x       solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_44_ldlt_solve(const cs_real_t    facto[10],
                     const cs_real_t    rhs[4],
                     cs_real_t          x[4])
{
  /* Sanity check */
  assert(facto != NULL && rhs != NULL && x != NULL);

  x[0] = rhs[0];
  x[1] = rhs[1] - x[0]*facto[1];
  x[2] = rhs[2] - x[0]*facto[3] - x[1]*facto[4];
  x[3] = rhs[3] - x[0]*facto[6] - x[1]*facto[7] - x[2]*facto[8];

  x[3] = x[3]*facto[9];
  x[2] = x[2]*facto[5] - facto[8]*x[3];
  x[1] = x[1]*facto[2] - facto[7]*x[3] - facto[4]*x[2];
  x[0] = x[0]*facto[0] - facto[6]*x[3] - facto[3]*x[2] - facto[1]*x[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a 6x6 matrix with a modified Cholesky decomposition (L.D.L^T)
 *         The solution should be already allocated
 * Ref. http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      f    vector of the coefficients of the decomposition
 * \param[in]      b    right-hand side
 * \param[in,out]  x    solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_66_ldlt_solve(const cs_real_t    f[21],
                     const cs_real_t    b[3],
                     cs_real_t          x[3])
{
  /* Sanity check */
  assert(f != NULL && b != NULL && x != NULL);

  x[0] = b[0];
  x[1] = b[1] - x[0]*f[1];
  x[2] = b[2] - x[0]*f[3]  - x[1]*f[4];
  x[3] = b[3] - x[0]*f[6]  - x[1]*f[7]  - x[2]*f[8];
  x[4] = b[4] - x[0]*f[10] - x[1]*f[11] - x[2]*f[12] - x[3]*f[13];
  x[5] = b[5] - x[0]*f[15] - x[1]*f[16] - x[2]*f[17] - x[3]*f[18] - x[4]*f[19];

  x[5] = x[5]*f[20];
  x[4] = x[4]*f[14] -f[19]*x[5];
  x[3] = x[3]*f[9]  -f[18]*x[5] -f[13]*x[4];
  x[2] = x[2]*f[5]  -f[17]*x[5] -f[12]*x[4] -f[8]*x[3];
  x[1] = x[1]*f[2]  -f[16]*x[5] -f[11]*x[4] -f[7]*x[3] -f[4]*x[2];
  x[0] = x[0]*f[0]  -f[15]*x[5] -f[10]*x[4] -f[6]*x[3] -f[3]*x[2] -f[1]*x[1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a SPD matrix with a L.D.L^T (Modified Cholesky decomposition)
 *         The solution should be already allocated
 * Reference:
 * http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]       n_rows   dimension of the system to solve
 * \param[in]       facto    vector of the coefficients of the decomposition
 * \param[in]       rhs      right-hand side
 * \param[in, out]  sol      solution
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_ldlt_solve(int                n_rows,
                  const cs_real_t   *facto,
                  const cs_real_t   *rhs,
                  cs_real_t         *sol)
{
  /* Sanity check */
  assert(facto != NULL && rhs != NULL && sol != NULL);

  if (n_rows == 1) {
    sol[0] = rhs[0] * facto[0];
    return;
  }

  /* 1 - Solving Lz = b with forward substitution :
   *     z_i = b_i - \sum_{k=0}^{i-1} l_ik * z_k
   */

  sol[0] = rhs[0]; /* case i = 0 */

  short int rowi_idx = 0;
  for (short int i = 1; i < n_rows; i++){

    rowi_idx += i;

    const cs_real_t  *l_i = facto + rowi_idx;
    cs_real_t  sum = 0.;
    for (short int k = 0; k < i; k++)
      sum += sol[k] * l_i[k];
    sol[i] = rhs[i] - sum;

  } /* forward substitution */

  /* 2 - Solving Dy = z and facto^Tx=y with backwards substitution
   *     x_i = z_i/d_ii - \sum_{k=i+1}^{n} l_ki * x_k
   */

  const short int  last_row_id = n_rows - 1;
  const int  shift = n_rows*(last_row_id)/2;   // idx with n_rows - 1
  int  diagi_idx = shift + last_row_id;        // last entry of the facto.
  sol[last_row_id] *= facto[diagi_idx];         // 1 / d_nn

  for (short int i = last_row_id - 1; i >= 0; i--) {

    diagi_idx -= (i+2);
    sol[i] *= facto[diagi_idx];

    short int  rowk_idx = shift;
    cs_real_t  sum = 0.0;
    for (short int k = last_row_id; k > i; k--) {
      /*sol[i] -= facto[k*(k+1)/2+i] * sol[k];*/
      const cs_real_t  *l_k = facto + rowk_idx;
      sum += l_k[i] * sol[k];
      rowk_idx -= k;
    }
    sol[i] -= sum;

  } /* backward substitution */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a small dense matrix
 *
 * \param[in]  mat         pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_simple_dump(const cs_sdm_t     *mat)
{
  if (mat == NULL)
    return;

  if (mat->n_rows < 1 || mat->n_cols < 1) {
    cs_log_printf(CS_LOG_DEFAULT, " No value.\n");
    return;
  }

  for (short int i = 0; i < mat->n_rows; i++) {
    for (short int j = 0; j < mat->n_cols; j++)
      cs_log_printf(CS_LOG_DEFAULT, " % .4e", mat->val[i*mat->n_cols+j]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a small dense matrix
 *
 * \param[in]  parent_id   id of the related parent entity
 * \param[in]  row_ids     list of ids related to associated entities (or NULL)
 * \param[in]  col_ids     list of ids related to associated entities (or NULL)
 * \param[in]  mat         pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_dump(cs_lnum_t           parent_id,
            const cs_lnum_t    *row_ids,
            const cs_lnum_t    *col_ids,
            const cs_sdm_t     *mat)
{
  if (mat == NULL)
    return;

  cs_log_printf(CS_LOG_DEFAULT, "<< MATRIX parent id: %d >>\n", parent_id);

  if (mat->n_rows < 1 || mat->n_cols < 1) {
    cs_log_printf(CS_LOG_DEFAULT, " No value.\n");
    return;
  }

  if (row_ids == NULL || col_ids == NULL)
    cs_sdm_simple_dump(mat);

  else {

    cs_log_printf(CS_LOG_DEFAULT, " %8s %11d", " ", col_ids[0]);
    for (short int i = 1; i < mat->n_cols; i++)
      cs_log_printf(CS_LOG_DEFAULT, " %11d", col_ids[i]);
    cs_log_printf(CS_LOG_DEFAULT, "\n");

    for (short int i = 0; i < mat->n_rows; i++) {
      cs_log_printf(CS_LOG_DEFAULT, " %8d ", row_ids[i]);
      for (short int j = 0; j < mat->n_cols; j++)
        cs_log_printf(CS_LOG_DEFAULT, " % .4e", mat->val[i*mat->n_cols+j]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a small dense matrix defined by blocks
 *
 * \param[in]  parent_id   id of the related parent entity
 * \param[in]  mat         pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_dump(cs_lnum_t           parent_id,
                  const cs_sdm_t     *mat)
{
  if (mat == NULL)
    return;

  if ((mat->flag & CS_SDM_BY_BLOCK) == 0) {
    cs_sdm_simple_dump(mat);
    return;
  }
  assert(mat->block_desc != NULL);

  cs_log_printf(CS_LOG_DEFAULT, "\n << BLOCK MATRIX parent id: %d >>\n",
                parent_id);

  int  n_b_rows = mat->block_desc->n_row_blocks;
  int  n_b_cols = mat->block_desc->n_col_blocks;

  if (n_b_rows < 1 || n_b_cols < 1) {
    cs_log_printf(CS_LOG_DEFAULT, " No block\n");
    return;
  }
  cs_log_printf(CS_LOG_DEFAULT, " n_row_blocks: %d; n_col_blocks: %d\n",
                n_b_rows, n_b_cols);

  cs_sdm_t  *blocks = mat->block_desc->blocks;

  for (short int bi = 0; bi < n_b_rows; bi++) {
    for (short int bj = 0; bj < n_b_cols; bj++) {

      cs_sdm_t  *bij = blocks + bi*n_b_cols + bj;
      cs_log_printf(CS_LOG_DEFAULT, "<< BLOCK (%2d, %2d) >>\n", bi, bj);
      cs_sdm_simple_dump(bij);

    } /* Block j */
  } /* Block i */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
