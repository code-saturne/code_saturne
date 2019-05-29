#ifndef __CS_SDM_H__
#define __CS_SDM_H__

/*============================================================================
 * Set of operations to handle Small Dense Matrices (SDM)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_SDM_BY_BLOCK    (1 << 0) /* Matrix is defined by block */
#define CS_SDM_SYMMETRIC   (1 << 1) /* Matrix is symmetric by construction */
#define CS_SDM_SHARED_VAL  (1 << 2) /* Matrix is not owner of its values */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_sdm_t cs_sdm_t;

typedef struct {

  int    n_max_blocks_by_row;
  int    n_row_blocks;
  int    n_max_blocks_by_col;
  int    n_col_blocks;

  /* Allocated to n_max_blocks_by_row*n_max_blocks_by_col
     Cast in cs_sdm_t where values are shared with the master cs_sdm_t struct.
  */
  cs_sdm_t    *blocks;

} cs_sdm_block_t;

/* Structure enabling the repeated usage of Small Dense Matrices (SDM) during
   the building of the linear system by a cellwise process */
struct _cs_sdm_t {

  cs_flag_t   flag;        /* Metadata */

  /* Row-related members */
  int         n_max_rows;  // max number of entities by primal cells
  int         n_rows;      // current number of entities

  /* Column-related members. Not useful if the matrix is square */
  int         n_max_cols;  // max number of entries in a column
  int         n_cols;      // current number of columns

  cs_real_t  *val;         // small dense matrix (size: n_max_rows*n_max_cols)

  /* Structure describing the matrix in terms of blocks */
  cs_sdm_block_t   *block_desc;

};

/*============================================================================
 * Prototypes for pointer of functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic prototype for computing a local dense matrix-product
 *          c = a*b where c has been previously allocated
 *
 * \param[in]      a     local dense matrix to use
 * \param[in]      b     local dense matrix to use
 * \param[in, out] c     result of the local matrix-product
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_sdm_product_t) (const cs_sdm_t   *a,
                    const cs_sdm_t   *b,
                    cs_sdm_t         *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Generic prototype for computing a dense matrix-vector product
 *          mv has been previously allocated
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] mv     result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_sdm_matvec_t) (const cs_sdm_t    *mat,
                   const cs_real_t   *vec,
                   cs_real_t         *mv);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Basic dot product for a small vector
 *         For very small array sizes (3, 4, 6) prefer functions in cs_math
 *         For large array sizes ( from 10^3 to ..) prefer functions in cs_blas
 *
 * \param[in]  n      size of arrays x and y (small)
 * \param[in]  x      first array
 * \param[in]  y      second array
 *
 * \return  the dot product
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_sdm_dot(int               n,
           const cs_real_t   x[],
           const cs_real_t   y[])
{
  cs_real_t  dp = 0;

  if (x == NULL || y == NULL)
    return dp;

  for (int i = 0; i < n; i++)
    dp += x[i]*y[i];

  return dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Multiply a small vector by a scalar coefficient: \a y = \a s \a x
 *         For very small array sizes (3, 4, 6) prefer functions in cs_math
 *         For large array sizes ( from 10^3 to ..) prefer functions in cs_blas
 *
 * \param[in]      n      size of arrays x and y (small)
 * \param[in]      s      scalar coefficient
 * \param[in]      x      in array
 * \param[in,out]  y      out array
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_scalvect(int               n,
                const cs_real_t   s,
                const cs_real_t   x[],
                cs_real_t         y[])
{
  if (x == NULL || y == NULL)
    return;

  for (int i = 0; i < n; i++)
    y[i] = s * x[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Multiply a small vector by a scalar coefficient: \a y += \a s \a x
 *         For very small array sizes (3, 4, 6) prefer functions in cs_math
 *         For large array sizes ( from 10^3 to ..) prefer functions in cs_blas
 *
 * \param[in]      n      size of arrays x and y (small)
 * \param[in]      s      scalar coefficient
 * \param[in]      x      in array
 * \param[in,out]  y      out array
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_add_scalvect(int               n,
                    const cs_real_t   s,
                    const cs_real_t   x[],
                    cs_real_t         y[])
{
  if (x == NULL || y == NULL)
    return;

  for (int i = 0; i < n; i++)
    y[i] += s * x[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure
 *          Most generic function to create a cs_sdm_t structure
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
              int         n_max_cols);

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
cs_sdm_square_create(int   n_max_rows);

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
cs_sdm_create_copy(const cs_sdm_t   *m);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix which is its transpose.
 *
 * \param[in] mat   local matrix to transpose
 *
 * \return a pointer to the new allocated transposed matrix
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_create_transpose(cs_sdm_t  *mat);

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
cs_sdm_block_create(int          n_max_blocks_by_row,
                    int          n_max_blocks_by_col,
                    const int    max_row_block_sizes[],
                    const int    max_col_block_sizes[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure by block when the
 *          block size is constant and equal to 3
 *
 * \param[in]  n_max_blocks_by_row    max number of blocks in a row
 * \param[in]  n_max_blocks_by_col    max number of blocks in a column
 *
 * \return  a new allocated cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_block33_create(int      n_max_blocks_by_row,
                      int      n_max_blocks_by_col);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Map an array into a predefined cs_sdm_t structure. This array is
 *          shared and the lifecycle of this array is not managed by the
 *          cs_sdm_t structure
 *
 * \param[in]      n_max_rows   max. number of rows
 * \param[in]      n_max_cols   max. number of columns
 * \param[in, out] m            pointer to a cs_sdm_t structure to set
 * \param[in, out] array        pointer to an array of values of size equal to
 *                              n_max_rows x n_max_cols
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_map_array(int          n_max_rows,
                 int          n_max_cols,
                 cs_sdm_t    *m,
                 cs_real_t   *array)
{
  assert(array != NULL && m != NULL);  /* Sanity check */

  m->flag = CS_SDM_SHARED_VAL;
  m->n_rows = m->n_max_rows = n_max_rows;
  m->n_cols = m->n_max_cols = n_max_cols;
  m->val = array;
  m->block_desc = NULL;
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
cs_sdm_free(cs_sdm_t  *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a cs_sdm_t structure
 *          Case of a square matrix
 *
 * \param[in]      n_rows   current number of rows
 * \param[in]      n_cols   current number of columns
 * \param[in, out] mat      pointer to the cs_sdm_t structure to initialize
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_init(int         n_rows,
            int         n_cols,
            cs_sdm_t   *mat)
{
  assert(mat != NULL);

  mat->n_rows = n_rows;
  mat->n_cols = n_cols;
  memset(mat->val, 0, n_rows*n_cols*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a cs_sdm_t structure
 *          Case of a square matrix
 *
 * \param[in]      n_rows   current number of rows
 * \param[in, out] mat      pointer to the cs_sdm_t structure to initialize
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_square_init(int         n_rows,
                   cs_sdm_t   *mat)
{
  assert(mat != NULL);

  mat->n_rows = mat->n_cols = n_rows; /* square matrix */
  memset(mat->val, 0, n_rows*n_rows*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the pattern of cs_sdm_t structure defined by block
 *          The matrix should have been allocated before calling this function
 *
 * \param[in, out] m
 * \param[in]      n_row_blocks      number of blocks in a row
 * \param[in]      n_col_blocks      number of blocks in a column
 * \param[in]      row_block_sizes   number of rows by block in a column
 * \param[in]      col_block_sizes   number of columns by block in a row
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_init(cs_sdm_t      *m,
                  int            n_row_blocks,
                  int            n_col_blocks,
                  const int      row_block_sizes[],
                  const int      col_block_sizes[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the pattern of cs_sdm_t structure defined by 3x3 block
 *          The matrix should have been allocated before calling this function
 *
 * \param[in, out] m
 * \param[in]      n_row_blocks      number of blocks in a row
 * \param[in]      n_col_blocks      number of blocks in a column
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block33_init(cs_sdm_t     *m,
                    int           n_row_blocks,
                    int           n_col_blocks);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Copy a cs_sdm_t structure into another cs_sdm_t structure
 *          which has been already allocated
 *
 * \param[in, out]  recv    pointer to a cs_sdm_t struct.
 * \param[in]       send    pointer to a cs_sdm_t struct.
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_copy(cs_sdm_t        *recv,
            const cs_sdm_t  *send)
{
  /* Sanity check */
  assert(recv->n_max_rows >= send->n_rows);
  assert(recv->n_max_cols >= send->n_cols);

  recv->flag = send->flag;
  recv->n_rows = send->n_rows;
  recv->n_cols = send->n_cols;

  /* Copy values */
  memcpy(recv->val, send->val, sizeof(cs_real_t)*send->n_rows*send->n_cols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_sdm_t structure w.r.t. to a given
 *          matrix
 *
 * \param[in]  mref       pointer to a matrix to copy
 *
 * \return  a new allocated cs_sdm_t structure which is a copy of mref
 */
/*----------------------------------------------------------------------------*/

cs_sdm_t *
cs_sdm_block_create_copy(const cs_sdm_t   *mref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get a specific block in a cs_sdm_t structure defined by block
 *
 * \param[in]       m              pointer to a cs_sdm_t struct.
 * \param[in]       row_block_id   id of the block row, zero-based.
 * \param[in]       col_block_id   id of the block column, zero-based.
 *
 * \return a pointer to a cs_sdm_t structure corresponfing to a block
 */
/*----------------------------------------------------------------------------*/

static inline cs_sdm_t *
cs_sdm_get_block(const cs_sdm_t    *const m,
                 int                      row_block_id,
                 int                      col_block_id)
{
  /* Sanity checks */
  assert(m != NULL);
  assert(m->flag & CS_SDM_BY_BLOCK && m->block_desc != NULL);
  assert(col_block_id < m->block_desc->n_col_blocks);
  assert(row_block_id < m->block_desc->n_row_blocks);

  return  m->block_desc->blocks
    + row_block_id*m->block_desc->n_col_blocks + col_block_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get a copy of a column in a preallocated vector
 *
 * \param[in]       m         pointer to a cs_sdm_t struct.
 * \param[in]       col_id    id of the column, zero-based.
 * \param[in, out]  col_vals  extracted values
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_get_col(const cs_sdm_t    *m,
               int                col_id,
               cs_real_t         *col_vals)
{
  /* Sanity checks */
  assert(m != NULL && col_vals != NULL);
  assert(col_id < m->n_cols);

  const cs_real_t  *_col = m->val + col_id;
  for(int i = 0; i < m->n_rows; i++, _col += m->n_cols)
    col_vals[i] = *_col;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  copy a block of a matrix into a sub-matrix starting from (r_id, c_id)
 *  with a size of nr rows and nc cols
 *
 * \param[in]      m      pointer to cs_sdm_t structure
 * \param[in]      r_id   row index
 * \param[in]      c_id   column index
 * \param[in]      nr     number of rows to extract
 * \param[in]      nc     number of column to extract
 * \param[in,out]  b      submatrix
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_copy_block(const cs_sdm_t       *m,
                  const short int       r_id,
                  const short int       c_id,
                  const short int       nr,
                  const short int       nc,
                  cs_sdm_t             *b)
{
  /* Sanity checks */
  assert(m != NULL && b != NULL);
  assert(r_id >= 0 && c_id >= 0);
  assert((r_id + nr) <= m->n_rows);
  assert((c_id + nc) <= m->n_cols);
  assert(nr == b->n_rows && nc == b->n_cols);

  const cs_real_t *_start = m->val + c_id + r_id*m->n_cols;
  for (short int i = 0; i < nr; i++, _start += m->n_cols)
    memcpy(b->val + i*nc, _start, sizeof(cs_real_t)*nc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  transpose and copy a matrix into another one already shaped
 *         sub-matrix starting from (r_id, c_id)
 *
 * \param[in]      m      pointer to cs_sdm_t structure
 * \param[in, out] mt     matrix to update with the transposed of m
 *
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_transpose_and_update(const cs_sdm_t   *m,
                            cs_sdm_t         *mt)
{
  assert(m != NULL && mt != NULL);
  assert(m->n_rows == mt->n_cols && m->n_cols == mt->n_rows);

  for (short int i = 0; i < m->n_rows; i++) {
    const cs_real_t  *m_i = m->val + i*m->n_cols;
    for (short int j = 0; j < m->n_cols; j++)
      mt->val[j*mt->n_cols + i] += m_i[j];
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
                cs_sdm_t         *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Specific version: a (1x3) and b (3x1)
 *
 * \param[in]      a    local matrix to use
 * \param[in]      b    local matrix to use
 * \param[in, out] c    result of the local matrix-product (already allocated)
 *                      is updated
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_multiply_r1c3_rowrow(const cs_sdm_t   *a,
                            const cs_sdm_t   *b,
                            cs_sdm_t         *c)
{
  /* Sanity check */
  assert(a != NULL && b != NULL && c != NULL);
  assert(a->n_cols == 3 && b->n_cols == 3 &&
         a->n_rows == 1 && c->n_rows == 1 &&
         c->n_cols == 1 && b->n_rows == 1);

  c->val[0] += a->val[0]*b->val[0] + a->val[1]*b->val[1] + a->val[2]*b->val[2];
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
 *                      is updated
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_multiply_rowrow(const cs_sdm_t   *a,
                       const cs_sdm_t   *b,
                       cs_sdm_t         *c);

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
                           cs_sdm_t         *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Compute a row-row matrix product of a and b. It is basically equal
 *           to the classical a*b^T. It is a fast (matrices are row-major) way
 *           of computing a*b if b is symmetric or if b^T is given.
 *           Case of a matrix defined by block
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
                             cs_sdm_t         *c);

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
                                 cs_sdm_t         *c);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a small square matrix
 *          mv has been previously allocated
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] mv     result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_matvec(const cs_sdm_t    *mat,
                     const cs_real_t   *vec,
                     cs_real_t         *mv);

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
              cs_real_t         *mv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a rectangular matrix
 *          mv has been previously allocated and initialized
 *          Thus mv is updated inside this function
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use (size = mat->n_cols)
 * \param[in, out] mv     result of the operation (size = mat->n_rows)
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_update_matvec(const cs_sdm_t    *mat,
                     const cs_real_t   *vec,
                     cs_real_t         *mv);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a rectangular matrix
 *          which is transposed.
 *          mv has been previously allocated. mv is updated inside this
 *          function. Don't forget to initialize mv if needed.
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use (size = mat->n_cols)
 * \param[in, out] mv     result of the operation (size = mat->n_rows)
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_matvec_transposed(const cs_sdm_t    *mat,
                         const cs_real_t   *vec,
                         cs_real_t         *mv);

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
                 const cs_sdm_t  *add);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two matrices defined by block: loc += mult_coef * add
 *
 * \param[in, out] mat         local matrix storing the result
 * \param[in]      mult_coef   multiplicative coefficient
 * \param[in]      add         values to add to mat
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_add_mult(cs_sdm_t        *mat,
                      cs_real_t        mult_coef,
                      const cs_sdm_t  *add);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a dense matrix-vector product for a rectangular matrix
 *          defined by block
 *          mv has been previously allocated
 *
 * \param[in]      mat    local matrix to use
 * \param[in]      vec    local vector to use (size = mat->n_cols)
 * \param[in, out] mv     result of the operation (size = mat->n_rows)
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_matvec(const cs_sdm_t    *mat,
                    const cs_real_t   *vec,
                    cs_real_t         *mv);

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
           const cs_sdm_t  *add);

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
                const cs_sdm_t  *add);

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
                            cs_sdm_t  *tr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given matrix to two times its symmetric part
 *          mat --> mat + mat_tr = 2*symm(mat)
 *
 * \param[in, out] mat   small dense matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_2symm(cs_sdm_t   *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given matrix into its anti-symmetric part
 *
 * \param[in, out] mat   small dense matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_square_asymm(cs_sdm_t   *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given block matrix into its anti-symmetric part
 *
 * \param[in, out] mat   small dense matrix defined by block to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_square_asymm(cs_sdm_t   *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Decompose a matrix into the matrix product Q.R
 *         Case of a 3x3 symmetric matrix
 *
 * \param[in]      m      matrix values
 * \param[in, out] Qt     transposed of matrix Q
 * \param[in, out] R      vector of the coefficient of the decomposition
 *
 * \note: R is an upper triangular matrix. Stored in a compact way.
 *
 *    j= 0, 1, 2
 *  i=0| 0| 1| 2|
 *  i=1   | 4| 5|
 *  i=2        6|
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_sym_qr_compute(const cs_real_t   m[9],
                         cs_real_t         Qt[9],
                         cs_real_t         R[6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LU factorization of a small dense 3x3 matrix.
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    compact storage of coefficients for the LU
 *                          factorization
 *
 * \note: facto stores L the lower triangular matrix (without its diagonal
 *        entries assumed to be equal to 1) and U the upper triangular matrix.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_lu_compute(const cs_sdm_t   *m,
                     cs_real_t         facto[9]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LU factorization of a small dense matrix. Small means that the
 *         number m->n_rows is less than 100 for instance.
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    compact storage of coefficients for the LU
 *                          factorization (should be allocated to the right
 *                          size, i.e. m->n_rows*m->n_rows)
 *
 * \note: facto stores L the lower triangular matrix (without its diagonal
 *        entries assumed to be equal to 1) and U the upper triangular matrix.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_lu_compute(const cs_sdm_t   *m,
                  cs_real_t         facto[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a system A.sol = rhs using a LU factorization of A (a small
 *         3x3 dense matrix).
 *
 * \param[in]       facto    compact storage of coefficients for the LU
 *                           factorization (should be allocated to the right
 *                           size, i.e. n_rows*n_rows)
 * \param[in]       rhs      right-hand side
 * \param[in, out]  sol      solution
 *
 * \note: facto stores L the lower triangular matrix (without its diagonal
 *        entries assumed to be equal to 1) and U the upper triangular matrix.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_lu_solve(const cs_real_t    facto[9],
                   const cs_real_t    rhs[3],
                   cs_real_t          sol[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a system A.sol = rhs using a LU factorization of A (a small
 *         dense matrix).
 *
 * \param[in]       n_rows   dimension of the system to solve
 * \param[in]       facto    compact storage of coefficients for the LU
 *                           factorization (should be allocated to the right
 *                           size, i.e. n_rows*n_rows)
 * \param[in]       rhs      right-hand side
 * \param[in, out]  sol      solution
 *
 * \note: facto stores L the lower triangular matrix (without its diagonal
 *        entries assumed to be equal to 1) and U the upper triangular matrix.
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_lu_solve(cs_lnum_t          n_rows,
                const cs_real_t    facto[],
                const cs_real_t   *rhs,
                cs_real_t         *sol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 3x3 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_33_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 4x4 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_44_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[10]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  LDL^T: Modified Cholesky decomposition of a 6x6 SPD matrix.
 *         For more reference, see for instance
 *   http://mathforcollege.com/nm/mws/gen/04sle/mws_gen_sle_txt_cholesky.pdf
 *
 * \param[in]      m        pointer to a cs_sdm_t structure
 * \param[in, out] facto    vector of the coefficient of the decomposition
 *
 * \note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_66_ldlt_compute(const cs_sdm_t   *m,
                       cs_real_t         facto[21]);

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
 * \note: facto is a lower triangular matrix. The first value of the
 *        j-th (zero-based) row is easily accessed: its index is j*(j+1)/2
 *        (cf sum of the first j natural numbers). Instead of 1 on the diagonal
 *        we store the inverse of D mat in the L.D.L^T decomposition
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_ldlt_compute(const cs_sdm_t     *m,
                    cs_real_t          *facto,
                    cs_real_t          *dkk);

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
                     cs_real_t          sol[3]);

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
                     cs_real_t          x[4]);

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
                     const cs_real_t    b[6],
                     cs_real_t          x[6]);

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
                  cs_real_t         *sol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test if a matrix is symmetric. Return 0. if the extradiagonal
 *          differences are lower thann the machine precision.
 *
 * \param[in]  mat         pointer to the cs_sdm_t structure to test
 *
 * \return  0 if the matrix is symmetric at the machine tolerance otherwise
 *          the absolute max. value between two transposed terms
 */
/*----------------------------------------------------------------------------*/

double
cs_sdm_test_symmetry(const cs_sdm_t     *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a small dense matrix
 *
 * \param[in]  mat         pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_simple_dump(const cs_sdm_t     *mat);

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
            const cs_sdm_t     *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a cs_sdm_t structure not defined by block
 *          Print into the file f if given otherwise open a new file named
 *          fname if given otherwise print into the standard output
 *          The usage of threshold allows one to compare more easier matrices
 *          without taking into account numerical roundoff.
 *
 * \param[in]  fp         pointer to a file structure or NULL
 * \param[in]  fname      filename or NULL
 * \param[in]  thd        threshold (below this value --> set 0)
 * \param[in]  m          pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_fprintf(FILE             *fp,
               const char       *fname,
               cs_real_t         thd,
               const cs_sdm_t   *m);

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
                  const cs_sdm_t     *mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a cs_sdm_t structure which is defined by block
 *          Print into the file f if given otherwise open a new file named
 *          fname if given otherwise print into the standard output
 *          The usage of threshold allows one to compare more easier matrices
 *          without taking into account numerical roundoff.
 *
 * \param[in]  fp         pointer to a file structure or NULL
 * \param[in]  fname      filename or NULL
 * \param[in]  thd        threshold (below this value --> set 0)
 * \param[in]  m          pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_sdm_block_fprintf(FILE             *fp,
                     const char       *fname,
                     cs_real_t         thd,
                     const cs_sdm_t   *m);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SDM_H__ */
