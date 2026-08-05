#ifndef CS_SDM_H
#define CS_SDM_H

/*============================================================================
 * Set of operations to handle Small Dense Matrices (SDM)
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

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <cassert>
#include <cstring>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

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

struct cs_sdm_block_t {

  int    n_max_blocks_by_row;
  int    n_row_blocks;
  int    n_max_blocks_by_col;
  int    n_col_blocks;

  /* Allocated to n_max_blocks_by_row*n_max_blocks_by_col
     Cast in cs_sdm_t where values are shared with the master cs_sdm_t struct.
  */
  cs_sdm_t    *blocks;

} ;

/* Structure enabling the repeated usage of Small Dense Matrices (SDM) during
   the building of the linear system by a cellwise process */

struct _cs_sdm_t {

  cs_flag_t   flag;        /* Metadata */

  /* Row-related members */
  int         n_max_rows;  // max number of entries in a row
  int         n_rows;      // current number of entities

  /* Column-related members. Not useful if the matrix is square */
  int         n_max_cols;  // max number of entries in a column
  int         n_cols;      // current number of columns

  cs_real_t *val; // small dense matrix (size: n_max_rows*n_max_cols)
                  // storage is row-major

  /* Structure describing the matrix in terms of blocks */
  cs_sdm_block_t   *block_desc;

  /*===========================================================================
   * Constructors/Destructors
   *==========================================================================*/

  _cs_sdm_t(const cs_flag_t flag_,
            const int       n_max_rows_,
            const int       n_max_cols_);

  _cs_sdm_t(const int n_max_rows_);

  _cs_sdm_t(const _cs_sdm_t &other);

  ~_cs_sdm_t();

  /*===========================================================================
   * Operators
   *==========================================================================*/

  /* Deleted for the moment */
  _cs_sdm_t &
  operator=(const _cs_sdm_t &) = delete;

  _cs_sdm_t &
  operator=(cs_real_t value);

  cs_real_t &
  operator()(const int i,
             const int j);

  const cs_real_t &
  operator()(const int i,
             const int j) const;

  _cs_sdm_t &
  operator*=(cs_real_t scaling);

  _cs_sdm_t &
  operator+=(const _cs_sdm_t &add);

  /*===========================================================================
   * Inline methods
   *==========================================================================*/

  /*-------------------------------------------------------------------------*/
  /*!
   * \brief Return the size of the matrix
   */
  /*-------------------------------------------------------------------------*/
  inline int
  size() const
  {
    return n_rows * n_cols;
  }

  /*-------------------------------------------------------------------------*/
  /*!
   * \brief Return a pointer on the ith row
   *
   * \param[in] row_i  index of the row
   *
   * \return pointer on row_i
   */
  /*-------------------------------------------------------------------------*/

  inline cs_real_t *
  row(const int row_i) const
  {
    assert(0 <= row_i && row_i < n_rows);
    return val + row_i * n_cols;
  }

  /*-------------------------------------------------------------------------*/
  /*!
   * \brief Return a pointer on the matrix values
   */
  /*--------------------------------------------------------------------------*/
  inline cs_real_t *
  data() const
  {
    return val;
  }

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Get a specific block in a cs_sdm_t structure defined by block
   *
   * \param[in] row_block_id  id of the block row, zero-based.
   * \param[in] col_block_id  id of the block column, zero-based.
   *
   * \return a pointer to a cs_sdm_t structure corresponfing to a block
   */
  /*--------------------------------------------------------------------------*/
  inline cs_sdm_t *
  get_block(const int row_block_id, const int col_block_id) const
  {
    /* Sanity checks */
    assert(flag & CS_SDM_BY_BLOCK && block_desc != nullptr);
    assert(col_block_id < block_desc->n_col_blocks);
    assert(row_block_id < block_desc->n_row_blocks);

    return block_desc->blocks + row_block_id * block_desc->n_col_blocks +
      col_block_id;
  }

  /*===========================================================================
   * Public methods
   *==========================================================================*/

  void
  init(const int nrows,
       const int ncols);

  void
  init(const int nrows);

  void
  set_zero();

  void
  fill(cs_real_t value);

  _cs_sdm_t
  transpose() const;

  void
  symmetrize_ur();

  void
  matvec(const cs_real_t vec[],
         cs_real_t       mv[]) const;

  void
  dump() const;

  void
  dump(cs_lnum_t parent_id) const;

  void
  dump(cs_lnum_t        parent_id,
       const cs_lnum_t *row_ids,
       const cs_lnum_t *col_ids) const;

  void
  dump(FILE       *fp,
       const char *fname,
       cs_real_t   thd) const;

  void
  map_array(int        n_max_rows_,
            int        n_max_cols_,
            cs_real_t *array);

}; // structure _sdm_t

namespace cs {
  using sdm_t = cs_sdm_t;
}

/*============================================================================
 * Public inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_sdm_t structure into another cs_sdm_t structure
 *        which has been already allocated
 *
 * \param[in, out] recv  pointer to a cs_sdm_t struct.
 * \param[in]      send  pointer to a cs_sdm_t struct.
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_copy(cs_sdm_t       *recv,
            const cs_sdm_t *send)
{
  // Sanity check
  assert(recv->n_max_rows >= send->n_rows);
  assert(recv->n_max_cols >= send->n_cols);

  recv->flag = send->flag;
  recv->n_rows = send->n_rows;
  recv->n_cols = send->n_cols;

  // Copy values
  std::memcpy(recv->val,
              send->val,
              sizeof(cs_real_t) * send->n_rows * send->n_cols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a block of a matrix into a sub-matrix starting from (r_id, c_id)
 *        with a size of nr rows and nc cols
 *
 * \param[in]     m     pointer to cs_sdm_t structure
 * \param[in]     r_id  row index
 * \param[in]     c_id  column index
 * \param[in]     nr    number of rows to extract
 * \param[in]     nc    number of column to extract
 * \param[in,out] b     submatrix
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_copy_block(const cs_sdm_t  *m,
                  const int  r_id,
                  const int  c_id,
                  const int  nr,
                  const int  nc,
                  cs_sdm_t        *b)
{
  // Sanity checks
  assert(m != NULL && b != NULL);
  assert(r_id >= 0 && c_id >= 0);
  assert((r_id + nr) <= m->n_rows);
  assert((c_id + nc) <= m->n_cols);
  assert(nr == b->n_rows && nc == b->n_cols);

  const cs_real_t *_start = m->val + c_id + r_id*m->n_cols;
  for (int i = 0; i < nr; i++, _start += m->n_cols)
    std::memcpy(b->val + i * nc, _start, sizeof(cs_real_t) * nc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief transpose and copy a matrix into another one already shaped
 *        sub-matrix starting from (r_id, c_id)
 *
 * \param[in]      m   pointer to cs_sdm_t structure
 * \param[in, out] mt  matrix to update with the transposed of m
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_transpose_and_update(const cs_sdm_t *m,
                            cs_sdm_t       *mt)
{
  assert(m != NULL && mt != NULL);
  assert(m->n_rows == mt->n_cols && m->n_cols == mt->n_rows);

  for (int i = 0; i < m->n_rows; i++) {
    const cs_real_t  *m_i = m->val + i*m->n_cols;
    for (int j = 0; j < m->n_cols; j++)
      mt->val[j*mt->n_cols + i] += m_i[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Basic dot product for a small vector.
 *        For very small array sizes (3, 4, 6) prefer functions in cs_math
 *        For large array sizes (from 10^3 and higher) prefer functions defined
 *        in cs_blas.cpp file
 *
 * \param[in] n  size of arrays x and y (small)
 * \param[in] x  first array
 * \param[in] y  second array
 *
 * \return the dot product
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
cs_sdm_dot(int             n,
           const cs_real_t x[],
           const cs_real_t y[])
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
 * \brief Multiply a small vector by a scalar coefficient: \a y = \a s \a x
 *        For very small array sizes (3, 4, 6) prefer functions in cs_math
 *        For large array sizes ( from 10^3 to ..) prefer functions in cs_blas
 *
 * \param[in]     n  size of arrays x and y (small)
 * \param[in]     s  scalar coefficient
 * \param[in]     x  input array
 * \param[in,out] y  output array
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_scalvect(int             n,
                const cs_real_t s,
                const cs_real_t x[],
                cs_real_t       y[])
{
  if (x == NULL || y == NULL)
    return;

  for (int i = 0; i < n; i++)
    y[i] = s * x[i];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Multiply a small vector by a scalar coefficient: \a y += \a s \a x
 *        For very small array sizes (3, 4, 6) prefer functions in cs_math
 *        For large array sizes ( from 10^3 to ..) prefer functions in cs_blas
 *
 * \param[in]     n  size of arrays x and y (small)
 * \param[in]     s  scalar coefficient
 * \param[in]     x  input array
 * \param[in,out] y  output array
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_sdm_add_scalvect(int             n,
                    const cs_real_t s,
                    const cs_real_t x[],
                    cs_real_t       y[])
{
  if (x == NULL || y == NULL)
    return;

  for (int i = 0; i < n; i++)
    y[i] += s * x[i];
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

cs_sdm_t *
cs_sdm_create(cs_flag_t flag,
              int       n_max_rows,
              int       n_max_cols);

cs_sdm_t *
cs_sdm_square_create(int n_max_rows);

cs_sdm_t *
cs_sdm_create_copy(const cs_sdm_t *m);

cs_sdm_t *
cs_sdm_create_transpose(cs_sdm_t *mat);

cs_sdm_t *
cs_sdm_block_create(int       n_max_blocks_by_row,
                    int       n_max_blocks_by_col,
                    const int max_row_block_sizes[],
                    const int max_col_block_sizes[]);

cs_sdm_t *
cs_sdm_block33_create(int n_max_blocks_by_row,
                      int n_max_blocks_by_col);

cs_sdm_t *
cs_sdm_free(cs_sdm_t *mat);

void
cs_sdm_block_init(cs_sdm_t  *m,
                  int        n_row_blocks,
                  int        n_col_blocks,
                  const int  row_block_sizes[],
                  const int  col_block_sizes[]);

void
cs_sdm_block33_init(cs_sdm_t *m,
                    int       n_row_blocks,
                    int       n_col_blocks);

void
cs_sdm_block_33_to_xyz(const cs_sdm_t *mb33,
                       cs_sdm_t       *mbxyz);

cs_sdm_t *
cs_sdm_block_create_copy(const cs_sdm_t *mref);

void
cs_sdm_multiply(const cs_sdm_t *a,
                const cs_sdm_t *b,
                cs_sdm_t       *c);

void
cs_sdm_multiply_rowrow(const cs_sdm_t *a,
                       const cs_sdm_t *b,
                       cs_sdm_t       *c);

void
cs_sdm_multiply_rowrow_sym(const cs_sdm_t *a,
                           const cs_sdm_t *b,
                           cs_sdm_t       *c);

void
cs_sdm_block_multiply_rowrow(const cs_sdm_t *a,
                             const cs_sdm_t *b,
                             cs_sdm_t       *c);

void
cs_sdm_block_multiply_rowrow_sym(const cs_sdm_t *a,
                                 const cs_sdm_t *b,
                                 cs_sdm_t       *c);
void
cs_sdm_update_matvec(const cs_sdm_t  *mat,
                     const cs_real_t *vec,
                     cs_real_t       *mv);

void
cs_sdm_matvec_transposed(const cs_sdm_t  *mat,
                         const cs_real_t *vec,
                         cs_real_t       *mv);

void
cs_sdm_block_add(cs_sdm_t       *mat,
                 const cs_sdm_t *add);

void
cs_sdm_block_add_mult(cs_sdm_t       *mat,
                      cs_real_t       mult_coef,
                      const cs_sdm_t *add);

void
cs_sdm_block_matvec(const cs_sdm_t  *mat,
                    const cs_real_t *vec,
                    cs_real_t       *mv);

void
cs_sdm_add_mult(cs_sdm_t       *mat,
                cs_real_t       alpha,
                const cs_sdm_t *add);

void cs_sdm_add_block(cs_sdm_t       *mat,
                      const int       r_id,
                      const int       c_id,
                      const int       nr,
                      const int       nc,
                      const cs_sdm_t *add);

void cs_sdm_add_block_topleft(cs_sdm_t       *mat,
                              const int       nr,
                              const int       nc,
                              const cs_sdm_t *add);

void
cs_sdm_square_add_transpose(cs_sdm_t *mat,
                            cs_sdm_t *tr);

void
cs_sdm_square_2symm(cs_sdm_t *mat);

void
cs_sdm_square_asymm(cs_sdm_t *mat);

void
cs_sdm_block_square_asymm(cs_sdm_t *mat);

void
cs_sdm_33_sym_qr_compute(const cs_real_t m[9],
                         cs_real_t       q_tr[9],
                         cs_real_t       r[6]);

void
cs_sdm_33_lu_compute(const cs_sdm_t *m,
                     cs_real_t       facto[9]);

void
cs_sdm_lu_compute(const cs_sdm_t *m,
                  cs_real_t       facto[]);

void
cs_sdm_33_lu_solve(const cs_real_t facto[9],
                   const cs_real_t rhs[3],
                   cs_real_t       sol[3]);

void
cs_sdm_lu_solve(cs_lnum_t        n_rows,
                const cs_real_t  facto[],
                const cs_real_t *rhs,
                cs_real_t       *sol);

void
cs_sdm_33_ldlt_compute(const cs_sdm_t *m,
                       cs_real_t       facto[6]);

void
cs_sdm_44_ldlt_compute(const cs_sdm_t *m,
                       cs_real_t       facto[10]);

void
cs_sdm_66_ldlt_compute(const cs_sdm_t *m,
                       cs_real_t       facto[21]);

void
cs_sdm_ldlt_compute(const cs_sdm_t     *m,
                    cs_real_t          *facto,
                    cs_real_t          *dkk);

void
cs_sdm_33_ldlt_solve(const cs_real_t facto[6],
                     const cs_real_t rhs[3],
                     cs_real_t       sol[3]);

void
cs_sdm_44_ldlt_solve(const cs_real_t facto[10],
                     const cs_real_t rhs[4],
                     cs_real_t       x[4]);

void
cs_sdm_66_ldlt_solve(const cs_real_t f[21],
                     const cs_real_t b[6],
                     cs_real_t       x[6]);

void
cs_sdm_ldlt_solve(int              n_rows,
                  const cs_real_t *facto,
                  const cs_real_t *rhs,
                  cs_real_t       *sol);

double
cs_sdm_test_symmetry(const cs_sdm_t *mat);

void
cs_sdm_block_dump(cs_lnum_t       parent_id,
                  const cs_sdm_t *mat);
void
cs_sdm_block_fprintf(FILE           *fp,
                     const char     *fname,
                     cs_real_t       thd,
                     const cs_sdm_t *m);

/*----------------------------------------------------------------------------*/

#endif /* CS_SDM_H */
