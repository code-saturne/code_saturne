/*============================================================================
 * Assembly of local cellwise system into a cs_matrix_t structure through
 * the cs_matrix_assembler_t and its related structures
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_assembler_priv.h"
#include "cs_matrix_assembler.h"
#include "cs_param_cdo.h"
#include "cs_parall.h"
#include "cs_sort.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_assembly.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdo_assembly.c

  \brief Assembly of local cellwise system into a cs_matrix_t structure through
  the cs_matrix_assembler_t and its related structures.

  This function are specific to CDO schemes. Thus one can assume a more specific
  behavior in order to get a more optimzed version of the standard assembly
  process.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_CDO_ASSEMBLY_DBG    0  /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

static cs_cdo_assembly_t  **cs_cdo_assembly = NULL;

/*=============================================================================
 * Local function pointer definitions
 *============================================================================*/

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

typedef struct {

  /* Row numberings */

  cs_gnum_t   g_id;             /* Global row numbering */
  cs_lnum_t   l_id;             /* Local range set numbering */
  int         i;                /* Cellwise system numbering */

  /* Column-related members */

  int                 n_cols;    /* Number of columns (cellwise system) */
  cs_gnum_t          *col_g_id;  /* Global numbering of columns */
  int                *col_idx;   /* Array to build and to give as parameter
                                    for add_vals() function */

  const cs_real_t    *val;       /* Row values */
  cs_real_t          *expval;    /* Expanded row values (when unrolling non
                                    scalar-valued block) */

} cs_cdo_assembly_row_t;

struct _cs_cdo_assembly_t {

  int      n_cw_dofs;   /* Number of DoFs in a cell */
  int      ddim;        /* Number of real values related to each diagonal
                            entry */
  int      edim;        /* Number of real values related to each
                           extra-diagonal entry */

  /* When working with matrix build by scalar-valued blocks, one may need to
     shift the row and/or the column local ids */

  cs_lnum_t   l_col_shift;
  cs_lnum_t   l_row_shift;

  cs_cdo_assembly_row_t    *row;

};

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given local id in a given array of
 *        ordered local ids, when the id might not be present
 *
 * We assume the id is present in the array.
 *
 * \param[in]  start_id         begin search with this id
 * \param[in]  l_id_array size  array_size
 * \param[in]  l_id             local id to search for
 * \param[in]  l_id_array       ordered unique local ids array
 *
 * \return  index of l_id in l_id_array, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline int
_l_binary_search(int              start_id,
                 int              l_id_array_size,
                 cs_lnum_t        l_id,
                 const cs_lnum_t  l_id_array[])
{
  assert(start_id > -1);
  int  end_id = l_id_array_size - 1;

  while (start_id <= end_id) {

    const int  mid_id = (start_id + end_id)/2;
    const cs_lnum_t test_val = l_id_array[mid_id];

    if (test_val < l_id)
      start_id = mid_id + 1;
    else if (test_val > l_id)
      end_id = mid_id - 1;
    else
      return mid_id;

  }

  return -1; /* Not found */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given global id in a given array of
 *        ordered global ids
 *
 * We assume the id is present in the array.
 *
 * \param[in]  g_id_array size  array_size
 * \param[in]  g_id             global id to search for
 * \param[in]  g_id_array       ordered unique global ids array
 *
 * \return  index of g_id in g_id_array or -1 if not found.
 */
/*----------------------------------------------------------------------------*/

static inline int
_g_binary_search(int              g_id_array_size,
                 cs_gnum_t        g_id,
                 const cs_gnum_t  g_id_array[])
{
  int  start_id = 0;
  int  end_id = g_id_array_size - 1;

  while (start_id <= end_id) {

    const int  mid_id = (end_id + start_id) / 2;
    const cs_gnum_t mid_val = g_id_array[mid_id];

    if (mid_val < g_id)
      start_id = mid_id + 1;
    else if (mid_val > g_id)
      end_id = mid_id - 1;
    else
      return mid_id;

  }

  return -1; /* Not found */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with no openMP and scalar-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_single(const cs_cdo_assembly_row_t    *row,
                        void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;

  /* Update the diagonal value */

  mc->_d_val[row->l_id] += row->val[row->i];

  /* Update the extra-diagonal values */

  cs_real_t  *xvals = mc->_e_val + ms->e.row_index[row->l_id];
  for (int j = 0; j < row->i; j++) /* Lower part */
    xvals[row->col_idx[j]] += row->val[j];
  for (int j = row->i+1; j < row->n_cols; j++) /* Upper part */
    xvals[row->col_idx[j]] += row->val[j];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with openMP atomic section and scalar-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_atomic(const cs_cdo_assembly_row_t    *row,
                        void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;

  /* Update the diagonal value */

# pragma omp atomic
  mc->_d_val[row->l_id] += row->val[row->i];

  /* Update the extra-diagonal values */

  cs_real_t  *xvals = mc->_e_val + ms->e.row_index[row->l_id];
  for (int j = 0; j < row->n_cols; j++) {
    if (j != row->i) {
#     pragma omp atomic
      xvals[row->col_idx[j]] += row->val[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with openMP critical section and scalar-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_critical(const cs_cdo_assembly_row_t    *row,
                          void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;

  /* Update the diagonal value */

# pragma omp critical
  {
    mc->_d_val[row->l_id] += row->val[row->i];

    /* Update the extra-diagonal values */

    cs_real_t  *xvals = mc->_e_val + ms->e.row_index[row->l_id];
    for (int j = 0; j < row->n_cols; j++)
      if (j != row->i)
        xvals[row->col_idx[j]] += row->val[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the column index for each entry of a row
 *        Case where the row belong to the local rank and all its colums too.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  ma     pointer to matrix assembler values structure
 * \param[in, out]  row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_set_col_idx_scal_l(const cs_matrix_assembler_t     *ma,
                    cs_cdo_assembly_row_t           *row)
{
  const cs_lnum_t  l_r_id = row->l_id; /* g_r_id - ma->l_range[0]; */
  const cs_lnum_t  l_start = ma->r_idx[l_r_id], l_end = ma->r_idx[l_r_id+1];
  const int  n_l_cols = l_end - l_start;
  const cs_lnum_t  *col_ids = ma->c_id + l_start;

  /* Loop on columns to fill col_idx for extra-diag entries
   * Diagonal is treated separately */

  for (int j = 0; j < row->i; j++) { /* Lower part */
    row->col_idx[j] = _l_binary_search(0,
                                       n_l_cols,
                         /* l_c_id */  row->col_g_id[j]-ma->l_range[0],
                                       col_ids);
    assert(row->col_idx[j] > -1);
  }

  for (int j = row->i + 1; j < row->n_cols; j++) { /* Upper part */
    row->col_idx[j] = _l_binary_search(0,
                                       n_l_cols,
                         /* l_c_id */  row->col_g_id[j]-ma->l_range[0],
                                       col_ids);
    assert(row->col_idx[j] > -1);
  }
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the column index for each entry of a row
 *        Case where the row belong to the local rank but a part of the columns
 *        belong to a distant rank. Hence the naming *_ld
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  ma     pointer to matrix assembler values structure
 * \param[in, out]  row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_set_col_idx_scal_ld(const cs_matrix_assembler_t    *ma,
                     cs_cdo_assembly_row_t          *row)
{
  assert(ma->d_r_idx != NULL); /* local-id-based function, need to adapt */

  const cs_lnum_t l_r_id = row->l_id; /* g_r_id - ma->l_range[0]; */
  const cs_lnum_t l_start = ma->r_idx[l_r_id], l_end = ma->r_idx[l_r_id+1];
  const cs_lnum_t d_start = ma->d_r_idx[l_r_id], d_end = ma->d_r_idx[l_r_id+1];
  const int n_d_cols = d_end - d_start;
  const int n_l_cols = l_end - l_start - n_d_cols;

  /* Loop on columns to fill col_idx for extra-diag entries */

  for (int j = 0; j < row->i; j++) {

    const cs_gnum_t g_c_id = row->col_g_id[j];

    if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) { /* Local part */

      row->col_idx[j] = _l_binary_search(0,
                                         n_l_cols,
                                         g_c_id - ma->l_range[0], /* l_c_id */
                                         ma->c_id + l_start);
      assert(row->col_idx[j] > -1);

    }
    else { /* Distant part */

      /* column ids start and end of local row, so add n_l_cols */

      row->col_idx[j] = n_l_cols + _g_binary_search(n_d_cols,
                                                    g_c_id,
                                                    ma->d_g_c_id + d_start);

    }

  } /* Loop on columns of the row */

  for (int j = row->i + 1; j < row->n_cols; j++) {

    const cs_gnum_t g_c_id = row->col_g_id[j];

    if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) { /* Local part */

      row->col_idx[j] = _l_binary_search(0,
                                         n_l_cols,
                                         g_c_id-ma->l_range[0], /* l_c_id */
                                         ma->c_id + l_start);
      assert(row->col_idx[j] > -1);

    }
    else { /* Distant part */

      /* column ids start and end of local row, so add n_l_cols */

      row->col_idx[j] = n_l_cols + _g_binary_search(n_d_cols,
                                                    g_c_id,
                                                    ma->d_g_c_id + d_start);

    }

  } /* Loop on columns of the row */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row does not belong to the local rank.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  mav    pointer to matrix assembler values structure
 * \param[in]       ma     pointer to matrix assembler values structure
 * \param[in]       row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_dt(cs_matrix_assembler_values_t       *mav,
                      const cs_matrix_assembler_t        *ma,
                      const cs_cdo_assembly_row_t        *row)
{
  /* Case where coefficient is handled by other rank. No need to call
     add_values() function in this case */

  assert(row->g_id < ma->l_range[0] || row->g_id >= ma->l_range[1]);

  const cs_lnum_t  e_r_id = _g_binary_search(ma->coeff_send_n_rows,
                                             row->g_id,
                                             ma->coeff_send_row_g_id);
  const cs_lnum_t  r_start = ma->coeff_send_index[e_r_id];
  const int  n_e_rows = ma->coeff_send_index[e_r_id+1] - r_start;
  const cs_gnum_t  *coeff_send_g_id = ma->coeff_send_col_g_id + r_start;

  /* Diagonal term */

  const cs_lnum_t  e_diag_id = r_start + _g_binary_search(n_e_rows,
                                                          row->g_id,
                                                          coeff_send_g_id);

  /* Now add values to send coefficients */

# pragma omp atomic
  mav->coeff_send[e_diag_id] += row->val[row->i];

  /* Loop on extra-diagonal entries */

  for (int j = 0; j < row->i; j++) { /* Lower-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

#   pragma omp atomic
    mav->coeff_send[e_id] += row->val[j];

  }

  for (int j = row->i + 1; j < row->n_cols; j++) { /* Upper-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

#   pragma omp atomic
    mav->coeff_send[e_id] += row->val[j];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row does not belong to the local rank. No openMP.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  mav    pointer to matrix assembler values structure
 * \param[in]       ma     pointer to matrix assembler values structure
 * \param[in]       row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_ds(cs_matrix_assembler_values_t       *mav,
                      const cs_matrix_assembler_t        *ma,
                      const cs_cdo_assembly_row_t        *row)
{
  /* Case where coefficient is handled by other rank. No need to call
     add_values() function in this case */

  assert(row->g_id < ma->l_range[0] || row->g_id >= ma->l_range[1]);

  const cs_lnum_t  e_r_id = _g_binary_search(ma->coeff_send_n_rows,
                                             row->g_id,
                                             ma->coeff_send_row_g_id);
  const cs_lnum_t  r_start = ma->coeff_send_index[e_r_id];
  const int  n_e_rows = ma->coeff_send_index[e_r_id+1] - r_start;
  const cs_gnum_t  *coeff_send_g_id = ma->coeff_send_col_g_id + r_start;

  /* Diagonal term */

  const cs_lnum_t  e_diag_id = r_start + _g_binary_search(n_e_rows,
                                                          row->g_id,
                                                          coeff_send_g_id);

  /* Now add values to send coefficients */

  mav->coeff_send[e_diag_id] += row->val[row->i];

  /* Loop on extra-diagonal entries */

  for (int j = 0; j < row->i; j++) { /* Lower-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

    mav->coeff_send[e_id] += row->val[j];

  }

  for (int j = row->i + 1; j < row->n_cols; j++) { /* Upper-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

    mav->coeff_send[e_id] += row->val[j];

  }
}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with no openMP and vector-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_vect_values_single(const cs_cdo_assembly_row_t    *row,
                        void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;

  cs_lnum_t stride;

  /* db_size != ebsize has not been handled by
   * cs_matrix_assembler_t in mpi case yet */

  assert(db_size == 3);
  assert(eb_size == 3);

  /* Update the diagonal value */

  stride = db_size * db_size;

  for (int j = 0; j < stride; j++)
    mc->_d_val[row->l_id*stride + j] += row->val[9*row->i + j];

  /* Update the extra-diagonal values */

  stride = eb_size * eb_size;

  cs_real_t  *xvals = mc->_e_val + stride*ms->e.row_index[row->l_id];
  for (int j = 0; j < row->i; j++) /* Lower part */
    for (int k = 0; k < stride; k++)
      xvals[row->col_idx[j]*stride + k] += row->val[9*j + k];

  for (int j = row->i+1; j < row->n_cols; j++) /* Upper part */
    for (int k = 0; k < stride; k++)
      xvals[row->col_idx[j]*stride + k] += row->val[9*j + k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with openMP atomic section and vector-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_vect_values_atomic(const cs_cdo_assembly_row_t    *row,
                        void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;

  cs_lnum_t stride;

  /* db_size != ebsize has not been handled by
   * cs_matrix_assembler_t in mpi case yet */

  assert(db_size == 3);
  assert(eb_size == 3);

  /* Update the diagonal value */

  stride = db_size * db_size;

  for (int j = 0; j < stride; j++) {
# pragma omp atomic
    mc->_d_val[row->l_id*stride + j] += row->val[9*row->i + j];
  }

  /* Update the extra-diagonal values */

  stride = eb_size * eb_size;

  cs_real_t  *xvals = mc->_e_val + stride*ms->e.row_index[row->l_id];
  for (int j = 0; j < row->n_cols; j++) {
    if (j != row->i) {
      for (int k = 0; k < stride; k++) {
#     pragma omp atomic
        xvals[row->col_idx[j]*stride + k] += row->val[9*j + k];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function pointer for adding values to a MSR matrix.
 *
 *  Specific case:
 *        CDO schemes with openMP critical section and vector-valued quantities
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \remark  Note that we pass column indexes (not ids) here; as the
 *          caller is already assumed to have identified the index
 *          matching a given column id.
 *
 * \param[in]      row         pointer to a cs_cdo_assembly_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_vect_values_critical(const cs_cdo_assembly_row_t    *row,
                          void                           *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_dist_t  *mc = matrix->coeffs;

  const cs_matrix_struct_dist_t  *ms = matrix->structure;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;

  cs_lnum_t stride;

  /* db_size != ebsize has not been handled by
   * cs_matrix_assembler_t in mpi case yet */

  assert(db_size == 3);
  assert(eb_size == 3);

  /* Update the diagonal value */

# pragma omp critical
  {
    stride = db_size*db_size;

    for (int j = 0; j < stride; j++)
      mc->_d_val[row->l_id*stride + j] += row->val[9*row->i + j];

    /* Update the extra-diagonal values */

    stride = eb_size*eb_size;

    cs_real_t  *xvals = mc->_e_val + stride*ms->e.row_index[row->l_id];
    for (int j = 0; j < row->n_cols; j++)
      if (j != row->i)
        for (int k = 0; k < stride; k++)
          xvals[row->col_idx[j]*stride + k] += row->val[9*j + k];
    }
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row does not belong to the local rank.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  mav    pointer to matrix assembler values structure
 * \param[in]       ma     pointer to matrix assembler values structure
 * \param[in]       row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_vect_dt(cs_matrix_assembler_values_t       *mav,
                      const cs_matrix_assembler_t        *ma,
                      const cs_cdo_assembly_row_t        *row)
{
  /* Case where coefficient is handled by other rank. No need to call
     add_values() function in this case */

  assert(row->g_id < ma->l_range[0] || row->g_id >= ma->l_range[1]);

  const cs_lnum_t  e_r_id = _g_binary_search(ma->coeff_send_n_rows,
                                             row->g_id,
                                             ma->coeff_send_row_g_id);
  const cs_lnum_t  r_start = ma->coeff_send_index[e_r_id];
  const int  n_e_rows = ma->coeff_send_index[e_r_id+1] - r_start;
  const cs_gnum_t  *coeff_send_g_id = ma->coeff_send_col_g_id + r_start;

  cs_matrix_t  *matrix = (cs_matrix_t *)mav->matrix;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;

  cs_lnum_t stride;

  /* db_size != ebsize has not been handled by
   * cs_matrix_assembler_t in mpi case yet */

  assert(db_size == 3);
  assert(eb_size == 3);

  /* Diagonal term */

  const cs_lnum_t  e_diag_id = r_start + _g_binary_search(n_e_rows,
                                                          row->g_id,
                                                          coeff_send_g_id);

  /* Now add values to send coefficients */

  stride = db_size*db_size;

  for (int k = 0; k < stride; k++) {
# pragma omp atomic
    mav->coeff_send[e_diag_id*stride + k] += row->val[9*row->i + k];
  }

  /* Loop on extra-diagonal entries */

  stride = eb_size*eb_size;

  for (int j = 0; j < row->i; j++) { /* Lower-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

    for (int k = 0; k < stride; k++) {
#   pragma omp atomic
      mav->coeff_send[e_id*stride + k] += row->val[9*row->i + k];
    }
  }

  for (int j = row->i + 1; j < row->n_cols; j++) { /* Upper-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */

    for (int k = 0; k < stride; k++) {
#   pragma omp atomic
      mav->coeff_send[e_id*stride + k] += row->val[9*row->i + k];
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row does not belong to the local rank. No openMP.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  mav    pointer to matrix assembler values structure
 * \param[in]       ma     pointer to matrix assembler values structure
 * \param[in]       row    pointer to a cs_cdo_assembly_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_vect_ds(cs_matrix_assembler_values_t       *mav,
                      const cs_matrix_assembler_t        *ma,
                      const cs_cdo_assembly_row_t        *row)
{
  /* Case where coefficient is handled by other rank. No need to call
     add_values() function in this case */

  assert(row->g_id < ma->l_range[0] || row->g_id >= ma->l_range[1]);

  const cs_lnum_t  e_r_id = _g_binary_search(ma->coeff_send_n_rows,
                                             row->g_id,
                                             ma->coeff_send_row_g_id);
  const cs_lnum_t  r_start = ma->coeff_send_index[e_r_id];
  const int  n_e_rows = ma->coeff_send_index[e_r_id+1] - r_start;
  const cs_gnum_t  *coeff_send_g_id = ma->coeff_send_col_g_id + r_start;

  cs_matrix_t  *matrix = (cs_matrix_t *)mav->matrix;
  const cs_lnum_t db_size = matrix->db_size;
  const cs_lnum_t eb_size = matrix->eb_size;

  cs_lnum_t stride;

  /* db_size != ebsize has not been handled by
   * cs_matrix_assembler_t in mpi case yet */

  assert(db_size == 3);
  assert(eb_size == 3);

  /* Diagonal term */

  const cs_lnum_t  e_diag_id = r_start + _g_binary_search(n_e_rows,
                                                          row->g_id,
                                                          coeff_send_g_id);

  /* Now add values to send coefficients */

  stride = db_size*db_size;

  for (int k = 0; k < stride; k++)
    mav->coeff_send[e_diag_id*stride + k] += row->val[9*row->i + k];

  /* Loop on extra-diagonal entries */

  stride = eb_size*eb_size;

  for (int j = 0; j < row->i; j++) { /* Lower-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */
    for (int k = 0; k < stride; k++)
      mav->coeff_send[e_id*stride + k] += row->val[9*j + k];

  }

  for (int j = row->i + 1; j < row->n_cols; j++) { /* Upper-part */

    const cs_lnum_t e_id = r_start + _g_binary_search(n_e_rows,
                                                      row->col_g_id[j],
                                                      coeff_send_g_id);

    /* Now add values to send coefficients */
    for (int k = 0; k < stride; k++)
      mav->coeff_send[e_id*stride + k] += row->val[9*j + k];
  }

}
#endif /* defined(HAVE_MPI) */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_cdo_assembly_t structure
 *
 * \param[in]      ddim         dim of the diagonal entries
 * \param[in]      edim         dim of the extra-diagonal entries
 * \param[in]      n_cw_dofs    number of DoF to be handled
 * \param[in, out] p_asb        double pointer to an assembly structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_assembly_struct(int                   ddim,
                      int                   edim,
                      int                   n_cw_dofs,
                      cs_cdo_assembly_t   **p_asb)
{
  cs_cdo_assembly_t  *asb = *p_asb;
  bool  create = false, reallocate = false;

  if (asb == NULL) {

    BFT_MALLOC(asb, 1, cs_cdo_assembly_t);
    create = true;

    /* Diagonal and extra-diagonal max. number of entries */

    asb->n_cw_dofs = n_cw_dofs;
    asb->ddim = ddim;
    asb->edim = edim;

  }
  else { /* Already allocated */

    if (asb->n_cw_dofs < n_cw_dofs)
      asb->n_cw_dofs = n_cw_dofs, reallocate = true;

    /* Diagonal and extra-diagonal max. number of entries */

    if (asb->ddim < ddim)
      asb->ddim = ddim, reallocate = true;
    if (asb->edim < edim)
      asb->edim = edim, reallocate = true;

  }

  /* When working with matrix build by scalar-valued blocks, one may need to
     shift the row and/or the column local ids */

  asb->l_row_shift = 0;
  asb->l_col_shift = 0;

  if (create) {

    /* Allocate the row structure used in the assembly process */

    BFT_MALLOC(asb->row, 1, cs_cdo_assembly_row_t);

    cs_cdo_assembly_row_t  *row = asb->row;

    if (asb->ddim < 2) {

      BFT_MALLOC(row->col_g_id, asb->n_cw_dofs, cs_gnum_t);
      BFT_MALLOC(row->col_idx, asb->n_cw_dofs, int);
      row->expval = NULL;

    }
    else {

      /* Temporary (until the global matrix is not defined by block) */

      int _size = asb->n_cw_dofs * asb->ddim;

      BFT_MALLOC(row->col_g_id, _size, cs_gnum_t);
      BFT_MALLOC(row->col_idx, _size, int);
      BFT_MALLOC(row->expval, asb->ddim*_size, cs_real_t);

    }

  }
  else { /* Update only */

    if (reallocate) {

      cs_cdo_assembly_row_t  *row = asb->row;
      assert(row != NULL);

      /* Re-allocate the row structure given the new sizes */

      if (asb->ddim < 2) {

        BFT_REALLOC(row->col_g_id, asb->n_cw_dofs, cs_gnum_t);
        BFT_REALLOC(row->col_idx, asb->n_cw_dofs, int);
        assert(row->expval == NULL);

      }
      else {

        /* Temporary (until the global matrix is not defined by block) */

        int _size = asb->n_cw_dofs * asb->ddim;

        BFT_REALLOC(row->col_g_id, _size, cs_gnum_t);
        BFT_REALLOC(row->col_idx, _size, int);
        BFT_REALLOC(row->expval, asb->ddim * _size, cs_real_t);

      }

    } /* reallocate = true */

  } /* update = true */

  /* Return the pointer */

  *p_asb = asb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_assembly_t structure
 *
 * \param[in, out]  p_asb    pointer to a structure pointer to be freed
 */
/*----------------------------------------------------------------------------*/

static void
_free_assembly_struct(cs_cdo_assembly_t  **p_asb)
{
  if (*p_asb == NULL)
    return;

  cs_cdo_assembly_t  *asb = *p_asb;

  if (asb->ddim > 1)
    BFT_FREE(asb->row->expval);

  BFT_FREE(asb->row->col_g_id);
  BFT_FREE(asb->row->col_idx);
  BFT_FREE(asb->row);

  BFT_FREE(asb);
  *p_asb = NULL;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_cdo_assembly_t structure related
 *         to a given thread
 *
 * \param[in]  t_id    id in the array of pointer
 *
 * \return a pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_assembly_t *
cs_cdo_assembly_get(int    t_id)
{
  if (t_id < 0 || t_id >= cs_glob_n_threads)
    return NULL;

  return cs_cdo_assembly[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate cs_cdo_assembly_t structure (shared among schemes). Each
 *         thread has its own copy of this structure to enable a multithreaded
 *         assembly process.
 *
 * \param[in]  ddim          max number of dof values on the diagonal part
 * \param[in]  edim          max number of dof values on the extra-diag. part
 * \param[in]  n_cw_dofs     max number of DoFs in a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_init(int     ddim,
                     int     edim,
                     int     n_cw_dofs)
{
  /* Common buffers for assemble usage */

  if (cs_cdo_assembly == NULL) {

    const int  n_threads = cs_glob_n_threads;
    BFT_MALLOC(cs_cdo_assembly, n_threads, cs_cdo_assembly_t *);
    for (int i = 0; i < n_threads; i++)
      cs_cdo_assembly[i] = NULL;

  }

#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
#pragma omp parallel
  {
    /* Each thread allocate its part. This yields a better memory affinity */

    int  t_id = omp_get_thread_num();
    _init_assembly_struct(ddim, edim, n_cw_dofs, &(cs_cdo_assembly[t_id]));
  }
#else
  _init_assembly_struct(ddim, edim, n_cw_dofs, &(cs_cdo_assembly[0]));
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free matrix-related structures used during the simulation.
 *         Display overall statistic about the assembly stage for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_finalize(void)
{
  /* Free shared buffers for the assembly process */

  for (int t_id = 0; t_id < cs_glob_n_threads; t_id++)
    _free_assembly_struct(&(cs_cdo_assembly[t_id]));

  BFT_FREE(cs_cdo_assembly);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current shift values to consider during the assembly stage
 *
 * \param[in, out] asb          pointer to a cs_cdo_assembly_t to update
 * \param[in]      l_row_shift  shift to apply to local row ids
 * \param[in]      l_col_shift  shift to apply to local col ids
 *
 * \return a function pointer cs_cdo_assembly_func_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_set_shift(cs_cdo_assembly_t    *asb,
                          cs_lnum_t             l_row_shift,
                          cs_lnum_t             l_col_shift)
{
  if (asb == NULL)
    return;

  assert(l_row_shift > -1 && l_col_shift > -1);

  asb->l_row_shift = l_row_shift;
  asb->l_col_shift = l_col_shift;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Rely on the generic cs_matrix_assembler_values_add_g() function
 *         Case of scalar-valued matrices.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_scal_generic(const cs_sdm_t                   *m,
                                    const cs_lnum_t                  *dof_ids,
                                    const cs_range_set_t             *rset,
                                    cs_cdo_assembly_t                *asb,
                                    cs_matrix_assembler_values_t     *mav)
{
  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->n_rows == m->n_cols);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < m->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  cs_gnum_t  r_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_gnum_t  c_gids[CS_CDO_ASSEMBLE_BUF_SIZE];

  if (CS_CDO_ASSEMBLE_BUF_SIZE > m->n_rows*m->n_cols) {

    int  bufsize = 0;
    for (int r = 0; r < m->n_rows; r++) {

      const cs_gnum_t r_gid = row->col_g_id[r];

      for (int c = 0; c < m->n_cols; c++) {

        r_gids[bufsize] = r_gid;
        c_gids[bufsize] = row->col_g_id[c];
        bufsize++;

      } /* Loop on columns */
    }   /* Loop on rows */

#   pragma omp critical
    cs_matrix_assembler_values_add_g(mav, bufsize, r_gids, c_gids, m->val);

  }
  else {  /* Buffer size is too small for this matrix. One needs several calls
             to *_values_add_g() */

    int  bufsize = 0;
    for (int r = 0; r < m->n_rows; r++) {

      const cs_gnum_t r_gid = row->col_g_id[r];

      for (int c = 0; c < m->n_cols; c++) {

        r_gids[bufsize] = r_gid;
        c_gids[bufsize] = row->col_g_id[c];
        bufsize++;

        if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#         pragma omp critical
          cs_matrix_assembler_values_add_g(mav,
                                           bufsize, r_gids, c_gids, m->val);

          bufsize = 0;
        }

      } /* Loop on columns */
    } /* Loop on rows */

    if (bufsize > 0) {
#     pragma omp critical
      cs_matrix_assembler_values_add_g(mav, bufsize, r_gids, c_gids, m->val);
      bufsize = 0;
    }

  } /* Test on the size of the local buffers */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Rely on the generic cs_matrix_assembler_values_add_g() function
 *         Case of vector-valued matrices with an expanded 33 block
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_e33_generic(const cs_sdm_t                  *m,
                                   const cs_lnum_t                 *dof_ids,
                                   const cs_range_set_t            *rset,
                                   cs_cdo_assembly_t               *asb,
                                   cs_matrix_assembler_values_t    *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;

  assert(m->n_rows == m->n_cols);
  assert(m->n_rows == 3*bd->n_row_blocks);
  assert(m->n_cols == 3*bd->n_col_blocks);

  cs_gnum_t  r_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_cdo_assembly_row_t  *row = asb->row;
  cs_real_t  *row_vals[3] = {row->expval,
                             row->expval + m->n_rows,
                             row->expval + 2*m->n_rows };

  if (CS_CDO_ASSEMBLE_BUF_SIZE < m->n_rows)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Increase the size of CS_CDO_ASSEMBLE_BUF_SIZE\n",
              __func__);

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < m->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  /* Fill rows 3 by 3 */

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {
    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;

      for (int k = 0; k < 3; k++) {
        row_vals[0][3*bj+k] = mvals[  k];
        row_vals[1][3*bj+k] = mvals[3+k];
        row_vals[2][3*bj+k] = mvals[6+k];
      }

    } /* Loop on column blocks */

    for (int k = 0; k < 3; k++) {

      row->g_id = row->col_g_id[3*bi+k]; /* global row id is constant for
                                            all columns */
      for (int i = 0; i < m->n_cols; i++)
        r_gids[i] = row->g_id;

#     pragma omp critical
      cs_matrix_assembler_values_add_g(mav,
                                       m->n_cols, r_gids, row->col_g_id,
                                       row_vals[k]);

    } /* Loop on the three rows in a block */

  } /* Loop on row blocks */
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case. Parallel and with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_mpit(const cs_sdm_t                   *m,
                            const cs_lnum_t                  *dof_ids,
                            const cs_range_set_t             *rset,
                            cs_cdo_assembly_t                *asb,
                            cs_matrix_assembler_values_t     *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_scal_dt(mav, ma, row);

    else {

      _set_col_idx_scal_ld(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_scal_values_critical(row, mav->matrix);
#else
      _add_scal_values_atomic(row, mav->matrix);
#endif

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_mpis(const cs_sdm_t                   *m,
                            const cs_lnum_t                  *dof_ids,
                            const cs_range_set_t             *rset,
                            cs_cdo_assembly_t                *asb,
                            cs_matrix_assembler_values_t     *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */
  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_scal_ds(mav, ma, row);

    else {

      _set_col_idx_scal_ld(ma, row);
      _add_scal_values_single(row, mav->matrix);

    }

  }
}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case. Sequential and with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_seqt(const cs_sdm_t                  *m,
                            const cs_lnum_t                 *dof_ids,
                            const cs_range_set_t            *rset,
                            cs_cdo_assembly_t               *asb,
                            cs_matrix_assembler_values_t    *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
    _add_scal_values_critical(row, mav->matrix);
#else
    _add_scal_values_atomic(row, mav->matrix);
#endif

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case.
 *         Sequential and without openMP.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_seqs(const cs_sdm_t                  *m,
                            const cs_lnum_t                 *dof_ids,
                            const cs_range_set_t            *rset,
                            cs_cdo_assembly_t               *asb,
                            cs_matrix_assembler_values_t    *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows; /* This is a square matrix */

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    _set_col_idx_scal_l(ma, row);
    _add_scal_values_single(row, mav->matrix);

  } /* Loop on rows */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise (no-block) matrix into the global matrix
 *         corresponding to a system of coupled equations.
 *         Scalar-valued case.
 *         Sequential and without openMP.
 *         Block matrices assembled from cellwise scalar-valued matrices
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_seqs(const cs_sdm_t                  *m,
                                const cs_lnum_t                 *dof_ids,
                                const cs_range_set_t            *rset,
                                cs_cdo_assembly_t               *asb,
                                cs_matrix_assembler_values_t    *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows; /* This is a square matrix */

  /* Switch to the global numbering */

  const cs_gnum_t  *_rset_g_id = rset->g_id + asb->l_col_shift;
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = _rset_g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  if (asb->l_col_shift == asb->l_row_shift) {

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                               /* cellwise numbering */
      row->g_id = row->col_g_id[i];             /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = m->val + i*row->n_cols;

      _set_col_idx_scal_l(ma, row);
      _add_scal_values_single(row, mav->matrix);

    } /* Loop on rows */

  }
  else {

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                                     /* cellwise numbering */
      row->g_id = rset->g_id[dof_ids[i] + asb->l_row_shift]; /* global num. */
      row->l_id = row->g_id - rset->l_range[0];      /* range set numbering */
      row->val = m->val + i*row->n_cols;

      _set_col_idx_scal_l(ma, row);
      _add_scal_values_single(row, mav->matrix);

    } /* Loop on rows */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise (no-block) matrix into the global matrix
 *         corresponding to a system of coupled equations.
 *         Scalar-valued case.
 *         Sequential and with openMP.
 *         Block matrices assembled from cellwise scalar-valued matrices
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_seqt(const cs_sdm_t                  *m,
                                const cs_lnum_t                 *dof_ids,
                                const cs_range_set_t            *rset,
                                cs_cdo_assembly_t               *asb,
                                cs_matrix_assembler_values_t    *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows; /* This is a square matrix */

  /* Switch to the global numbering */

  const cs_gnum_t  *_rset_g_id = rset->g_id + asb->l_col_shift;
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = _rset_g_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  if (asb->l_col_shift == asb->l_row_shift) {

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                               /* cellwise numbering */
      row->g_id = row->col_g_id[i];             /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = m->val + i*row->n_cols;

      _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_scal_values_critical(row, mav->matrix);
#else
      _add_scal_values_atomic(row, mav->matrix);
#endif

    } /* Loop on rows */

  }
  else {

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                                     /* cellwise numbering */
      row->g_id = rset->g_id[dof_ids[i] + asb->l_row_shift]; /* global num. */
      row->l_id = row->g_id - rset->l_range[0];      /* range set numbering */
      row->val = m->val + i*row->n_cols;

      _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_scal_values_critical(row, mav->matrix);
#else
      _add_scal_values_atomic(row, mav->matrix);
#endif

    } /* Loop on rows */

  }
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_matrix_sys_mpis(const cs_sdm_t                   *m,
                                const cs_lnum_t                  *dof_ids,
                                const cs_range_set_t             *rset,
                                cs_cdo_assembly_t                *asb,
                                cs_matrix_assembler_values_t     *mav)
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  const cs_gnum_t  *_rset_gc_id = rset->g_id + asb->l_col_shift;
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = _rset_gc_id[dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */

  if (asb->l_col_shift == asb->l_row_shift) {

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                               /* cellwise numbering */
      row->g_id = row->col_g_id[i];             /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = m->val + i*row->n_cols;

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Loop on rows */

  }
  else {

    const cs_gnum_t  *_rset_gr_id = rset->g_id + asb->l_row_shift;

    for (int i = 0; i < row->n_cols; i++) {

      row->i = i;                                 /* cellwise numbering */
      row->g_id = _rset_gr_id[dof_ids[i]];        /* global num. */
      row->l_id = row->g_id - rset->l_range[0];   /* range set numbering */
      row->val = m->val + i*row->n_cols;

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Loop on rows */

  } /* col_shift != row_shift */

}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_seqs(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + m->n_rows,
                          row->expval + 2*m->n_rows };

  assert(m->n_rows == m->n_rows);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;
      for (int k = 0; k < 3; k++) {
        _vxyz[0][3*bj+k] = mvals[  k];
        _vxyz[1][3*bj+k] = mvals[3+k];
        _vxyz[2][3*bj+k] = mvals[6+k];
      }

    } /* Loop on column-wise blocks */

    /* dof_ids is an interlaced array (get access to the next 3 values) */

    for (int k = 0; k < 3; k++) {
      row->i = 3*bi+k;                          /* cellwise numbering */
      row->g_id = row->col_g_id[3*bi+k];        /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vxyz[k];

      _set_col_idx_scal_l(ma, row);
      _add_scal_values_single(row, mav->matrix);

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_seqt(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + m->n_rows,
                          row->expval + 2*m->n_rows };

  assert(m->n_rows == m->n_rows);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;
      for (int k = 0; k < 3; k++) {
        _vxyz[0][3*bj+k] = mvals[  k];
        _vxyz[1][3*bj+k] = mvals[3+k];
        _vxyz[2][3*bj+k] = mvals[6+k];
      }

    } /* Loop on column-wise blocks */

    /* dof_ids is an interlaced array (get access to the next 3 values */

    for (int k = 0; k < 3; k++) {
      row->i = 3*bi+k;                          /* cellwise numbering */
      row->g_id = row->col_g_id[3*bi+k];        /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vxyz[k];

      _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_scal_values_critical(row, mav->matrix);
#else
      _add_scal_values_atomic(row, mav->matrix);
#endif

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_mpis(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + m->n_rows,
                          row->expval + 2*m->n_rows };

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;
      for (int k = 0; k < 3; k++) {
        _vxyz[0][3*bj+k] = mvals[  k];
        _vxyz[1][3*bj+k] = mvals[3+k];
        _vxyz[2][3*bj+k] = mvals[6+k];
      }

    } /* Loop on column-wise blocks */

    /* dof_ids is an interlaced array (get access to the next 3 values */

    for (int k = 0; k < 3; k++) {

      row->i = 3*bi+k;                          /* cellwise numbering */
      row->g_id = row->col_g_id[3*bi+k];        /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vxyz[k];

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock33_matrix_mpit(const cs_sdm_t               *m,
                                     const cs_lnum_t              *dof_ids,
                                     const cs_range_set_t         *rset,
                                     cs_cdo_assembly_t            *asb,
                                     cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + m->n_rows,
                          row->expval + 2*m->n_rows};

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;
      for (int k = 0; k < 3; k++) {
        _vxyz[0][3*bj+k] = mvals[  k];
        _vxyz[1][3*bj+k] = mvals[3+k];
        _vxyz[2][3*bj+k] = mvals[6+k];
      }

    } /* Loop on column-wise blocks */

    /* dof_ids is an interlaced array (get access to the next 3 values */

    for (int k = 0; k < 3; k++) {

      row->i = 3*bi+k;                          /* cellwise numbering */
      row->g_id = row->col_g_id[3*bi+k];        /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vxyz[k];

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
        _add_scal_values_critical(row, mav->matrix);
#else
        _add_scal_values_atomic(row, mav->matrix);
#endif

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_seqs(const cs_sdm_t               *m,
                                    const cs_lnum_t              *dof_ids,
                                    const cs_range_set_t         *rset,
                                    cs_cdo_assembly_t            *asb,
                                    cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz = row->expval;

  assert(m->n_rows == m->n_cols);

  row->n_cols = bd->n_row_blocks;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[dim*i]/dim];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;

      for (int k = 0; k < 9; k++) {
        _vxyz[9*bj+k] = mvals[k];
      }
    } /* Loop on column-wise blocks */

    row->i = bi;                              /* cellwise numbering */
    row->g_id = row->col_g_id[bi];            /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = _vxyz;

    /*All entries within one block share the same row and column */
    _set_col_idx_scal_l(ma, row);
    _add_vect_values_single(row, mav->matrix);

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_seqt(const cs_sdm_t               *m,
                                    const cs_lnum_t              *dof_ids,
                                    const cs_range_set_t         *rset,
                                    cs_cdo_assembly_t            *asb,
                                    cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz = row->expval;

  assert(m->n_rows == m->n_cols);

  row->n_cols = bd->n_row_blocks;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[dim*i]/dim];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;

      for (int k = 0; k < 9; k++) {
        _vxyz[9*bj+k] = mvals[k];
      }
    } /* Loop on column-wise blocks */

    row->i = bi;                              /* cellwise numbering */
    row->g_id = row->col_g_id[bi];            /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = _vxyz;

    /*All entries within one block share the same row and column */
    _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
    _add_vect_values_critical(row, mav->matrix);
#else
    _add_vect_values_atomic(row, mav->matrix);
#endif

  } /* Loop on row-wise blocks */
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_mpis(const cs_sdm_t               *m,
                                    const cs_lnum_t              *dof_ids,
                                    const cs_range_set_t         *rset,
                                    cs_cdo_assembly_t            *asb,
                                    cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz = row->expval;

  assert(m->n_rows == m->n_cols);

  row->n_cols = bd->n_row_blocks;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[dim*i]/dim];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;

      for (int k = 0; k < 9; k++) {
        _vxyz[9*bj+k] = mvals[k];
      }
    } /* Loop on column-wise blocks */

    row->i = bi;                              /* cellwise numbering */
    row->g_id = row->col_g_id[bi];            /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = _vxyz;

    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_vect_ds(mav, ma, row);

    else {

      /*All entries within one block share the same row and column */
      _set_col_idx_scal_ld(ma, row);
      _add_vect_values_single(row, mav->matrix);

    }

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_block33_matrix_mpit(const cs_sdm_t               *m,
                                    const cs_lnum_t              *dof_ids,
                                    const cs_range_set_t         *rset,
                                    cs_cdo_assembly_t            *asb,
                                    cs_matrix_assembler_values_t *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim >= 3);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vxyz = row->expval;

  assert(m->n_rows == m->n_cols);

  row->n_cols = bd->n_row_blocks;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[dim*i]/dim];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size 3 */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);
      const cs_real_t  *const mvals = mIJ->val;

      for (int k = 0; k < 9; k++) {
        _vxyz[9*bj+k] = mvals[k];
      }
    } /* Loop on column-wise blocks */

    row->i = bi;                              /* cellwise numbering */
    row->g_id = row->col_g_id[bi];            /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = _vxyz;

    /*All entries within one block share the same row and column */
    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_vect_dt(mav, ma, row);

    else {

      _set_col_idx_scal_ld(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_vect_values_critical(row, mav->matrix);
#else
      _add_vect_values_atomic(row, mav->matrix);
#endif

    }

  } /* Loop on row-wise blocks */
}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_seqs(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim < 19);
  assert(asb->edim == asb->ddim);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*m->n_rows;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size "dim" */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);

      for (int ki = 0; ki < dim; ki++) {
        const cs_real_t  *const mvals = mIJ->val + ki*dim;
        cs_real_t  *buf = _vpointer[ki] + dim*bj;
        for (int kj = 0; kj < dim; kj++)
          buf[kj] = mvals[kj];
      }

    } /* Loop on column blocks */

    /* dof_ids is an interlaced array (get access to the next "dim" values */

    for (int ki = 0; ki < dim; ki++) {
      row->i = dim*bi+ki;                       /* cellwise numbering */
      row->g_id = row->col_g_id[dim*bi+ki];     /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vpointer[ki];

      _set_col_idx_scal_l(ma, row);
      _add_scal_values_single(row, mav->matrix);
    }

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_seqt(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim < 19);
  assert(asb->edim == asb->ddim);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*m->n_rows;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size "dim" */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);

      for (int ki = 0; ki < dim; ki++) {
        const cs_real_t  *const mvals = mIJ->val + ki*dim;
        cs_real_t  *buf = _vpointer[ki] + dim*bj;
        for (int kj = 0; kj < dim; kj++)
          buf[kj] = mvals[kj];
      }

    } /* Loop on column blocks */

    /* dof_ids is an interlaced array (get access to the next "dim" values */

    for (int ki = 0; ki < dim; ki++) {
      row->i = dim*bi+ki;                       /* cellwise numbering */
      row->g_id = row->col_g_id[dim*bi+ki];     /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vpointer[ki];

      _set_col_idx_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
      _add_scal_values_critical(row, mav->matrix);
#else
      _add_scal_values_atomic(row, mav->matrix);
#endif
    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}

#if defined(HAVE_MPI)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_mpis(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim < 19);
  assert(asb->edim == asb->ddim);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*m->n_rows;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size "dim" */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);

      for (int ki = 0; ki < dim; ki++) {
        const cs_real_t  *const mvals = mIJ->val + ki*dim;
        cs_real_t  *buf = _vpointer[ki] + dim*bj;
        for (int kj = 0; kj < dim; kj++)
          buf[kj] = mvals[kj];
      }

    } /* Loop on column blocks */

    /* dof_ids is an interlaced array (get access to the next "dim" values */

    for (int ki = 0; ki < dim; ki++) {

      row->i = dim*bi+ki;                       /* cellwise numbering */
      row->g_id = row->col_g_id[dim*bi+ki];     /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vpointer[ki];

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise matrix into the global matrix
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      m        cellwise view of the algebraic system
 * \param[in]      dof_ids  local DoF numbering
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] asb      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_assembly_eblock_matrix_mpit(const cs_sdm_t                *m,
                                   const cs_lnum_t               *dof_ids,
                                   const cs_range_set_t          *rset,
                                   cs_cdo_assembly_t             *asb,
                                   cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_cdo_assembly_row_t  *row = asb->row;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(asb->ddim < 19);
  assert(asb->edim == asb->ddim);
  assert(row->expval != NULL);

  const int  dim = asb->ddim;

  /* Expand the values for a bundle of rows */

  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*m->n_rows;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */

  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[dof_ids[i]];

  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* Expand all the blocks for this row */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      /* mIJ matrices are small square matrices of size "dim" */

      const cs_sdm_t  *const mIJ = cs_sdm_get_block(m, bi, bj);

      for (int ki = 0; ki < dim; ki++) {
        const cs_real_t  *const mvals = mIJ->val + ki*dim;
        cs_real_t  *buf = _vpointer[ki] + dim*bj;
        for (int kj = 0; kj < dim; kj++)
          buf[kj] = mvals[kj];
      }

    } /* Loop on column blocks */

    /* dof_ids is an interlaced array (get access to the next "dim" values */

    for (int ki = 0; ki < dim; ki++) {

      row->i = dim*bi+ki;                       /* cellwise numbering */
      row->g_id = row->col_g_id[dim*bi+ki];     /* global numbering */
      row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
      row->val = _vpointer[ki];

      if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
        _assemble_row_scal_ds(mav, ma, row);

      else {

        _set_col_idx_scal_ld(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
        _add_scal_values_critical(row, mav->matrix);
#else
        _add_scal_values_atomic(row, mav->matrix);
#endif

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */
}
#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
