/*============================================================================
 * Operations on cs_cdo_system_t structures: matrix/vector multiplications
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

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_array.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_saddle_system.h"

/*----------------------------------------------------------------------------*/
BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_saddle_system.c
 *
 * \brief Set of functions to manipulate saddle-point systems defined by blocks
 *        relying on the cs_cdo_system_t structures: matrix/vector
 *        multiplications
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_SADDLE_SYSTEM_DBG        0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Static inline private function prototypes
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 3 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view
 *
 *        This is an update operation. Be careful that the resulting array has
 *        been initialized.
 *
 * \param[in]      x2       array for the second set
 * \param[in]      ub       unassembled block related to the (2,1) block
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m12_multiply_vector(const cs_cdo_system_ublock_t  *ub,
                     const cs_real_t               *x2,
                     cs_real_t                     *m12x2)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_lnum_t  n2_dofs =  m21_adj->n_elts; /* adjacency should be related
                                                  to the (2,1) block */
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {

    const cs_real_t  _x2 = x2[i2];
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  *m21_vals = m21_val + 3*j;
      cs_real_t  *_m12x2 = m12x2 + 3*m21_adj->ids[j];

#     pragma omp atomic
      _m12x2[0] += m21_vals[0] * _x2;
#     pragma omp atomic
      _m12x2[1] += m21_vals[1] * _x2;
#     pragma omp atomic
      _m12x2[2] += m21_vals[2] * _x2;

    } /* Loop on x1 elements associated to a given x2 DoF */

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The stride is equal to 1 for the operator m21 (unassembled)
 *        x2 corresponds to a "scatter" view
 *
 *        This is an update operation. Be careful that the resulting array has
 *        been initialized.
 *
 * \param[in]      x2       array for the second set
 * \param[in]      ub       unassembled block related to the (2,1) block
 * \param[in, out] m12x2    resulting array (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m12_multiply_scalar(const cs_cdo_system_ublock_t  *ub,
                     const cs_real_t               *x2,
                     cs_real_t                     *m12x2)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_lnum_t  n2_dofs =  m21_adj->n_elts; /* adjacency should be related
                                                  to the (2,1) block */
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (n2_dofs > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < n2_dofs; i2++) {

    const cs_real_t  _x2 = x2[i2];
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
#     pragma omp atomic
      m12x2[m21_adj->ids[j]] += _x2 * m21_val[j];

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 3 for the operator m21 operator (unassembled)
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first part
 * \param[in, out] m21x1  resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m21_multiply_vector(const cs_cdo_system_ublock_t  *ub,
                     const cs_real_t               *x1,
                     cs_real_t                     *m21x1)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
      _m21x1 += cs_math_3_dot_product(m21_val + 3*j, x1 + 3*m21_adj->ids[j]);

    m21x1[i2] = _m21x1;

  } /* Loop on x2 elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The stride is equal to 1 for the operator m21 operator (unassembled)
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first part
 * \param[in, out] m21x1  resulting vector (have to be allocated)
 */
/*----------------------------------------------------------------------------*/

static void
_m21_multiply_scalar(const cs_cdo_system_ublock_t  *ub,
                     const cs_real_t               *x1,
                     cs_real_t                     *m21x1)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    cs_real_t  _m21x1 = 0.;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++)
      _m21x1 += m21_val[j] * x1[m21_adj->ids[j]];

    m21x1[i2] = _m21x1;

  } /* Loop on x2 elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a part of the matrix-vector operation for a saddle-point
 *        system in the case of a scalar-valued M11 block
 *        x1 and x2 are in scatter view
 *
 *        a) Resulting vector of the operation m12*x2
 *        b) r2 = m21x1
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first set
 * \param[in]      x2     array for the second set
 * \param[in, out] m12x2  resulting array for the first set
 * \param[in, out] m21x1  resulting array for the second set
 */
/*----------------------------------------------------------------------------*/

static void
_m12x2_m21x1_scalar(const cs_cdo_system_ublock_t  *ub,
                    const cs_real_t               *x1,
                    const cs_real_t               *x2,
                    cs_real_t                     *m12x2,
                    cs_real_t                     *m21x1)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    const cs_real_t _x2 = x2[i2];
    cs_real_t  _m21x1 = 0;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  m21_value = m21_val[j];
      const cs_lnum_t  elt1_id = m21_adj->ids[j];

#     pragma omp atomic
      m12x2[elt1_id] += _x2 * m21_value;

      _m21x1 += m21_value * x1[elt1_id];

    }

    m21x1[i2] = _m21x1;

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a part of the matrix-vector operation for a saddle-point
 *        system in the case of a vector-valued M11 block
 *        x1 and x2 are in scatter view
 *
 *        a) Resulting vector of the operation m12*x2
 *        b) r2 = m21x1
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first set
 * \param[in]      x2     array for the second set
 * \param[in, out] m12x2  resulting array for the first set
 * \param[in, out] m21x1  resulting array for the second set
 */
/*----------------------------------------------------------------------------*/

static void
_m12x2_m21x1_vector(const cs_cdo_system_ublock_t  *ub,
                    const cs_real_t               *x1,
                    const cs_real_t               *x2,
                    cs_real_t                     *m12x2,
                    cs_real_t                     *m21x1)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    const cs_real_t _x2 = x2[i2];
    cs_real_t  _m21x1 = 0;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  *m21_values = m21_val + 3*j;
      const cs_lnum_t  shift = 3*m21_adj->ids[j];

      cs_real_t *_m12x2 = m12x2 + shift;

#     pragma omp atomic
      _m12x2[0] += _x2 * m21_values[0];
#     pragma omp atomic
      _m12x2[1] += _x2 * m21_values[1];
#     pragma omp atomic
      _m12x2[2] += _x2 * m21_values[2];

      _m21x1 += cs_math_3_dot_product(m21_values, x1 + shift);

    }

    m21x1[i2] = _m21x1;

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a part of the residual of a saddle-point system in the case
 *        of a scalar-valued M11 block
 *        x1 and x2 are in scatter view
 *
 *        a) Resulting vector of the operation m12*x2
 *        b) res2 = rhs2 - m21x1
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first set
 * \param[in]      x2     array for the second set
 * \param[in]      rhs2   rhs for the second set
 * \param[in, out] m12x2  resulting array for the first set
 * \param[in, out] res2   computed residual for the second set
 */
/*----------------------------------------------------------------------------*/

static void
_partial_residual_scalar(const cs_cdo_system_ublock_t  *ub,
                         const cs_real_t               *x1,
                         const cs_real_t               *x2,
                         const cs_real_t               *rhs2,
                         cs_real_t                     *m12x2,
                         cs_real_t                     *res2)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    const cs_real_t _x2 = x2[i2];
    cs_real_t  m21x1 = 0;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  m21_value = m21_val[j];
      const cs_lnum_t  elt1_id = m21_adj->ids[j];

#     pragma omp atomic
      m12x2[elt1_id] += _x2 * m21_value;

      m21x1 += m21_value * x1[elt1_id];

    }

    res2[i2] = rhs2[i2] - m21x1;

  } /* Loop on x2 DoFs */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a part of the residual of a saddle-point system in the case
 *        of a scalar-valued M11 block
 *        x1 and x2 are in scatter view
 *
 *        a) Resulting vector of the operation m12*x2
 *        b) res2 = rhs2 - m21x1
 *
 * \param[in]      ub     unassembled block related to the (2,1) block
 * \param[in]      x1     array for the first set
 * \param[in]      x2     array for the second set
 * \param[in]      rhs2   rhs for the second set
 * \param[in, out] m12x2  resulting array for the first set
 * \param[in, out] res2   computed residual for the second set
 */
/*----------------------------------------------------------------------------*/

static void
_partial_residual_vector(const cs_cdo_system_ublock_t  *ub,
                         const cs_real_t               *x1,
                         const cs_real_t               *x2,
                         const cs_real_t               *rhs2,
                         cs_real_t                     *m12x2,
                         cs_real_t                     *res2)
{
  const cs_adjacency_t  *m21_adj = ub->adjacency;
  const cs_real_t  *m21_val = ub->values;

# pragma omp parallel for if (m21_adj->n_elts > CS_THR_MIN)
  for (cs_lnum_t i2 = 0; i2 < m21_adj->n_elts; i2++) {

    const cs_real_t _x2 = x2[i2];
    cs_real_t  m21x1 = 0;
    for (cs_lnum_t j = m21_adj->idx[i2]; j < m21_adj->idx[i2+1]; j++) {

      const cs_real_t  *m21_values = m21_val + 3*j;
      const cs_lnum_t  shift = 3*m21_adj->ids[j];

      cs_real_t *_m12x2 = m12x2 + shift;

#     pragma omp atomic
      _m12x2[0] += _x2 * m21_values[0];
#     pragma omp atomic
      _m12x2[1] += _x2 * m21_values[1];
#     pragma omp atomic
      _m12x2[2] += _x2 * m21_values[2];

      m21x1 += cs_math_3_dot_product(m21_values, x1 + shift);

    }

    res2[i2] = rhs2[i2] - m21x1;

  } /* Loop on x2 DoFs */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the inverse of the diagonal of the (1,1)-block matrix
 *        The storage of a matrix is in a gather view and the resulting array is
 *        in a scatter view.
 *
 * \param[in] b11_max_size  max size related to the (1,1) block
 * \param[in] sh            pointer to a system helper structure
 *
 * \return a pointer to the computed array (scatter view)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_saddle_system_b11_inv_diag(cs_lnum_t                b11_max_size,
                              cs_cdo_system_helper_t  *sh)
{
  if (sh == nullptr)
    return nullptr;

  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);
  const cs_lnum_t  n11_rows = cs_matrix_get_n_rows(m11); /* gather view */

  cs_real_t *inv_diag_m11 = nullptr;
  assert(n11_rows <= b11_max_size);
  BFT_MALLOC(inv_diag_m11, b11_max_size, cs_real_t);

  switch (cs_cdo_system_get_matrix_class(sh, 0)) {

  case CS_CDO_SYSTEM_MATRIX_CS:
    {
      const cs_real_t  *diag_m11 = cs_matrix_get_diagonal(m11);

#     pragma omp parallel for if (n11_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n11_rows; i++)
        inv_diag_m11[i] = 1./diag_m11[i];
    }
    break;

  case CS_CDO_SYSTEM_MATRIX_HYPRE:
  case CS_CDO_SYSTEM_MATRIX_PETSC:
    {
      /* Get diagonal is not available in this case. Avoid ending with an
         error */

      cs_matrix_copy_diagonal(m11, inv_diag_m11);

#     pragma omp parallel for if (n11_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n11_rows; i++)
        inv_diag_m11[i] = 1./inv_diag_m11[i];
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid case.", __func__);
    break;

  }

  /* Switch to a scatter view */

  cs_range_set_scatter(cs_cdo_system_get_range_set(sh, 0),
                       CS_REAL_TYPE, 1, /* treated as scalar-valued up to now */
                       inv_diag_m11,    /* gathered view */
                       inv_diag_m11);   /* scatter view */

  return inv_diag_m11;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Perform a matrix-vector multiplication for the (1,1) block when the
 *        input vector is in a scatter state.

 *        Thus, one performs a scatter --> gather (before the multiplication)
 *        and a gather --> scatter operation after the multiplication. One
 *        assumes that m11x1 is allocated to the right size. No check is done.
 *
 * \param[in]      sh        pointer to a system helper structure
 * \param[in, out] vec       vector
 * \param[in, out] matvec    resulting vector for the matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b11_matvec(const cs_cdo_system_helper_t  *sh,
                            cs_real_t                     *vec,
                            cs_real_t                     *matvec)

{
  if (sh == nullptr || vec == nullptr || matvec == nullptr)
    return;

  /* Remark:
   * n_rows = n_gather_elts <= n_scatter_elts = n_dofs (mesh view) <= n_cols
   */

#if defined(DEBUG) /* Only case handled up to now */
  cs_cdo_system_block_t  *blk11 = sh->blocks[0];
  assert(blk11->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
#endif

  const cs_range_set_t  *b11_rset = cs_cdo_system_get_range_set(sh, 0);
  const cs_matrix_t  *m11 = cs_cdo_system_get_matrix(sh, 0);

  /* scatter view to gather view for the input vector */

  cs_range_set_gather(b11_rset,
                      CS_REAL_TYPE, 1, /* type and stride */
                      vec,             /* in:  size=n1_scatter_elts */
                      vec);            /* out: size=n1_gather_elts */

  cs_matrix_vector_multiply(m11, vec, matvec);

  /* gather view to scatter view (i.e. algebraic to mesh view) */

  cs_range_set_scatter(b11_rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       vec, vec);

  cs_range_set_scatter(b11_rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       matvec, matvec);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m12*x2
 *        The block (1,2) is stored in an unassembled way and corresponds to
 *        the (2,1) block. Thus, one considers a transposition of this block.
 *        x2 corresponds to a "scatter" view
 *
 * \param[in]      sh         pointer to a system helper structure
 * \param[in]      x2         array to be multiplied
 * \param[in, out] res        result array storing the matrix.vector operation
 * \param[in]      reset_res  false --> this is an update
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b12_matvec(const cs_cdo_system_helper_t  *sh,
                            const cs_real_t               *x2,
                            cs_real_t                     *res,
                            bool                           reset_res)
{
  if (sh == nullptr || x2 == nullptr)
    return;
  assert(res != nullptr);

  /* Get information from the system helper */

  assert(sh->type == CS_CDO_SYSTEM_SADDLE_POINT);
  assert(sh->n_blocks == 2);

  /* n2_dofs = n2_elts since stride is equal to 1 in all our cases of
     saddle-point systems up to now (these are cell unknowns) */

  cs_cdo_system_block_t  *blk12 = sh->blocks[1];

  if (reset_res)
    cs_array_real_fill_zero(sh->col_block_sizes[0], res);

  /* Handled only unassembled blocks when dealing with saddle-point systems */

  assert(blk12->type == CS_CDO_SYSTEM_BLOCK_UNASS);
  cs_cdo_system_ublock_t *ub12
    = static_cast<cs_cdo_system_ublock_t *>(blk12->block_pointer);

  switch (blk12->info.stride) {
  case 1:
    _m12_multiply_scalar(ub12, x2, res);
    break;

  case 3:
    _m12_multiply_vector(ub12, x2, res);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid stride (%d)\n",
              __func__, blk12->info.stride);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the resulting vector of the operation m21*x1
 *        The block (2,1) is stored in an unassembled way.
 *        x1 corresponds to a "scatter" view
 *
 * \param[in]      sh   pointer to a system helper structure
 * \param[in]      x1   array to be multiplied
 * \param[in, out] res  result array storing the matrix.vector operation
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_b21_matvec(const cs_cdo_system_helper_t  *sh,
                            const cs_real_t               *x1,
                            cs_real_t                     *res)
{
  if (sh == nullptr || x1 == nullptr)
    return;
  assert(res != nullptr);

  /* Get information from the system helper */

  assert(sh->type == CS_CDO_SYSTEM_SADDLE_POINT);
  assert(sh->n_blocks == 2);

  cs_cdo_system_block_t  *blk12 = sh->blocks[1];

  /* Handled only unassembled blocks when dealing with saddle-point systems */

  assert(blk12->type == CS_CDO_SYSTEM_BLOCK_UNASS);
  cs_cdo_system_ublock_t *ub12
    = static_cast<cs_cdo_system_ublock_t *>(blk12->block_pointer);

  switch (blk12->info.stride) {
  case 1:
    _m21_multiply_scalar(ub12, x1, res);
    break;

  case 3:
    _m21_multiply_vector(ub12, x1, res);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid stride (%d)\n",
              __func__, blk12->info.stride);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the matrix-vector operation for a saddle-point system
 *        r1 = M11.x1 + M12.x2 (result for the first set)
 *        r2 = M21.x1          (result for the second set)
 *
 * \param[in]      sh  pointer to a system helper structure
 * \param[in, out] x1  array for the first part
 * \param[in, out] x2  array for the second part
 * \param[in, out] r1  result array for the first set of DoFs
 * \param[in, out] r2  result array for the second set of DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_matvec(const cs_cdo_system_helper_t  *sh,
                        cs_real_t                     *x1,
                        cs_real_t                     *x2,
                        cs_real_t                     *r1,
                        cs_real_t                     *r2)
{
  if (sh == nullptr)
    return;
  if (r1 == nullptr || r2 == nullptr)
    return;
  assert(x1 != nullptr && x2 != nullptr);

  /* Get information from the system helper */

  assert(sh->type == CS_CDO_SYSTEM_SADDLE_POINT);
  assert(sh->n_blocks == 2);

  const cs_lnum_t  n1_dofs = sh->col_block_sizes[0]; /* scatter view */

  /* 2. r1 <-- M11.x1 */

  cs_saddle_system_b11_matvec(sh, x1, r1);

  /* 1. m12x2 = M12.x2
   *    r2 = M21.x1
   *
   * One assumes an unassembled block for the (2,1) and its transposed (1,2)
   * block
   */

  cs_cdo_system_block_t  *blk12 = sh->blocks[1];
  assert(blk12->type == CS_CDO_SYSTEM_BLOCK_UNASS);
  cs_cdo_system_ublock_t *ub12
    = static_cast<cs_cdo_system_ublock_t *>(blk12->block_pointer);

  cs_real_t *m12x2 = nullptr;
  BFT_MALLOC(m12x2, n1_dofs, cs_real_t);
  cs_array_real_fill_zero(n1_dofs, m12x2);

  switch (blk12->info.stride) {
  case 1:
    _m12x2_m21x1_scalar(ub12, x1, x2, m12x2, r2);
    break;

  case 3:
    _m12x2_m21x1_vector(ub12, x1, x2, m12x2, r2);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid stride (%d)\n",
              __func__, blk12->info.stride);
  }

  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         n1_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         m12x2);

  /* 3. r1 <-- r1 + m12x2 */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    r1[i1] += m12x2[i1];

  BFT_FREE(m12x2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the residual of the saddle-point system
 *        res1 = rhs1 - M11.x1 - M12.x2 (residual for the first set)
 *        res2 = rhs2 - M21.x1          (residual for the second set)
 *
 * \param[in]      sh    pointer to a system helper structure
 * \param[in, out] x1    array for the first part
 * \param[in, out] x2    array for the second part
 * \param[in, out] res1  resulting residual vector for the first set of DoFs
 * \param[in, out] res2  resulting residual vector for the second set of DoFs
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_system_residual(const cs_cdo_system_helper_t  *sh,
                          cs_real_t                     *x1,
                          cs_real_t                     *x2,
                          cs_real_t                     *res1,
                          cs_real_t                     *res2)
{
  if (sh == nullptr)
    return;
  if (res1 == nullptr || res2 == nullptr)
    return;
  assert(x1 != nullptr && x2 != nullptr);

  /* Get information from the system helper */

  assert(sh->type == CS_CDO_SYSTEM_SADDLE_POINT);
  assert(sh->n_blocks == 2);

  const cs_lnum_t  n1_dofs = sh->col_block_sizes[0]; /* scatter view */

  /* Two parts for each set of DoFs:
   * a) res1 = rhs1 - M11.x1 - M12.x2
   * b) res2 = rhs2 - M21.x1
   */

  const cs_real_t  *rhs1 = sh->rhs_array[0];
  const cs_real_t  *rhs2 = sh->rhs_array[1];

  /* 1. m12x2 = M12.x2 (the sign will be changed in the last operation)
   *    res2 = rhs2 - M21.x1
   *
   * One assumes an unassembled block for the (2,1) and its transposed (1,2)
   * block
   */

  cs_cdo_system_block_t  *blk12 = sh->blocks[1];
  assert(blk12->type == CS_CDO_SYSTEM_BLOCK_UNASS);
  cs_cdo_system_ublock_t *ub12
    = static_cast<cs_cdo_system_ublock_t *>(blk12->block_pointer);

  cs_real_t *m12x2 = nullptr;
  BFT_MALLOC(m12x2, n1_dofs, cs_real_t);
  cs_array_real_fill_zero(n1_dofs, m12x2);

  switch (blk12->info.stride) {
  case 1:
    _partial_residual_scalar(ub12, x1, x2, rhs2, m12x2, res2);
    break;

  case 3:
    _partial_residual_vector(ub12, x1, x2, rhs2, m12x2, res2);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid stride (%d)\n",
              __func__, blk12->info.stride);
  }

  const cs_range_set_t  *rset = cs_cdo_system_get_range_set(sh, 0);

  if (rset->ifs != nullptr)
    cs_interface_set_sum(rset->ifs,
                         n1_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         m12x2);

  /* 2. res1 <-- M11.x1 */

  cs_saddle_system_b11_matvec(sh, x1, res1);

  /* 3. res1 <-- rhs1 - res1 - m12x2 */

# pragma omp parallel for if (n1_dofs > CS_THR_MIN)
  for (cs_lnum_t i1 = 0; i1 < n1_dofs; i1++)
    res1[i1] = rhs1[i1] - res1[i1] - m12x2[i1];

  BFT_FREE(m12x2);
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
