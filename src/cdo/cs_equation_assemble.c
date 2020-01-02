/*============================================================================
 * Assembly of local cellwise system into a cs_matrix_t structure through
 * the cs_matrix_assembler_t and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_assemble.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_equation_assemble.c

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

#define CS_EQUATION_ASSEMBLE_DBG          0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

/* Store the matrix structure and its assembler structures for each family
   of space discretizations */
static cs_matrix_assembler_t  **cs_equation_assemble_ma = NULL;
static cs_matrix_structure_t  **cs_equation_assemble_ms = NULL;
static cs_equation_assemble_t  **cs_equation_assemble = NULL;

static cs_timer_counter_t  cs_equation_ms_time;

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
  cs_real_t          *expval;    /* Expanded row values */

} cs_equation_assemble_row_t;

struct _cs_equation_assemble_t {

  int         ddim;         /* Number of real values related to each diagonal
                               entry */
  int         edim;         /* Number of real values related to each
                               extra-diagonal entry */

  cs_equation_assemble_row_t    *row;
};

/*============================================================================
 * Static inline private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of scalar-valued matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_equation_assembly_t *
_set_scalar_assembly_func(void)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_matrix_mpis;
    else                       /* With OpenMP */
      return cs_equation_assemble_matrix_mpit;

  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) { /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_matrix_seqs;
    else                       /* With OpenMP */
      return cs_equation_assemble_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of block 3x3 matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_equation_assembly_t *
_set_block33_assembly_func(void)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_eblock33_matrix_mpis;
    else                      /* With OpenMP */
      return cs_equation_assemble_eblock33_matrix_mpit;

  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) {  /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_eblock33_matrix_seqs;
    else                      /* With OpenMP */
      return cs_equation_assemble_eblock33_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Choose which function will be used to perform the matrix assembly
 *         Case of block NxN matrices.
 *
 * \return  a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_equation_assembly_t *
_set_block_assembly_func(void)
{
#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {  /* Parallel */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_eblock_matrix_mpis;
    else                      /* With OpenMP */
      return cs_equation_assemble_eblock_matrix_mpit;

  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks <= 1) {  /* Sequential */

    if (cs_glob_n_threads < 2) /* Without OpenMP */
      return cs_equation_assemble_eblock_matrix_seqs;
    else                      /* With OpenMP */
      return cs_equation_assemble_eblock_matrix_seqt;

  }

  return NULL; /* Case not handled */
}

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
 * \param[in]      row         pointer to a cs_equation_assemble_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_single(const cs_equation_assemble_row_t   *row,
                        void                               *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Update the diagonal value */
  mc->_d_val[row->l_id] += row->val[row->i];

  /* Update the extra-diagonal values */
  cs_real_t  *xvals = mc->_x_val + ms->row_index[row->l_id];
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
 * \param[in]      row         pointer to a cs_equation_assemble_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_atomic(const cs_equation_assemble_row_t   *row,
                        void                               *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Update the diagonal value */
# pragma omp atomic
  mc->_d_val[row->l_id] += row->val[row->i];

  /* Update the extra-diagonal values */
  cs_real_t  *xvals = mc->_x_val + ms->row_index[row->l_id];
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
 * \param[in]      row         pointer to a cs_equation_assemble_row_t type
 * \param[in, out] matrix_p    untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_add_scal_values_critical(const cs_equation_assemble_row_t   *row,
                          void                               *matrix_p)
{
  assert(row->l_id > -1);

  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Update the diagonal value */
# pragma omp critical
  {
    mc->_d_val[row->l_id] += row->val[row->i];

    /* Update the extra-diagonal values */
    cs_real_t  *xvals = mc->_x_val + ms->row_index[row->l_id];
    for (int j = 0; j < row->n_cols; j++)
      if (j != row->i)
        xvals[row->col_idx[j]] += row->val[j];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row belong to the local rank and all its colums too.
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  ma     pointer to matrix assembler values structure
 * \param[in, out]  row    pointer to a cs_equation_assemble_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_l(const cs_matrix_assembler_t     *ma,
                     cs_equation_assemble_row_t      *row)
{
  assert(ma->d_r_idx == NULL); /* local-id-based function, need to adapt */

  const cs_lnum_t  l_r_id = row->l_id; /* g_r_id - ma->l_range[0]; */
  const cs_lnum_t  l_start = ma->r_idx[l_r_id], l_end = ma->r_idx[l_r_id+1];
  const int  n_l_cols = l_end - l_start;
  const cs_lnum_t  *col_ids = ma->c_id + l_start;

  /* Loop on columns to fill col_idx for extra-diag entries
   *  Diagonal is treated separately */
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
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *        Case where the row belong to the local rank but a part of the columns
 *        belong to a distant rank. Hence the naming *_ld
 *
 * See \ref cs_matrix_assembler_values_add_g which performs the same operations
 * In the specific case of CDO system, one assumes predefined choices in
 * order to get a more optimized version of this function
 *
 * \param[in, out]  ma     pointer to matrix assembler values structure
 * \param[in, out]  row    pointer to a cs_equation_assemble_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_ld(const cs_matrix_assembler_t      *ma,
                      cs_equation_assemble_row_t       *row)
{
  assert(ma->d_r_idx != NULL); /* local-id-based function, need to adapt */

  const cs_lnum_t l_r_id = row->l_id; // g_r_id - ma->l_range[0];
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
 * \param[in]       row    pointer to a cs_equation_assemble_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_dt(cs_matrix_assembler_values_t         *mav,
                      const cs_matrix_assembler_t          *ma,
                      const cs_equation_assemble_row_t     *row)
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
 * \param[in]       row    pointer to a cs_equation_assemble_row_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble_row_scal_ds(cs_matrix_assembler_values_t         *mav,
                      const cs_matrix_assembler_t          *ma,
                      const cs_equation_assemble_row_t     *row)
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_equation_assemble_t structure
 *
 * \param[in]  max_ddim        max. dim of the diagonal entries
 * \param[in]  max_edim        max. dim of the extra-diagonal entries
 * \param[in]  n_max_cw_dofs   max. number of DoF to be handled
 *
 * \return a pointer to a new allocated cs_equation_assemble_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_equation_assemble_t *
_init_equation_assembler_struct(int      max_ddim,
                                int      max_edim,
                                int      n_max_cw_dofs)
{
  cs_equation_assemble_t  *eqa = NULL;

  BFT_MALLOC(eqa, 1, cs_equation_assemble_t);

  eqa->ddim = max_ddim;
  eqa->edim = max_edim;

  BFT_MALLOC(eqa->row, 1, cs_equation_assemble_row_t);
  if (max_ddim < 2) {
    BFT_MALLOC(eqa->row->col_g_id, n_max_cw_dofs, cs_gnum_t);
    BFT_MALLOC(eqa->row->col_idx, n_max_cw_dofs, int);
  }
  else {
    n_max_cw_dofs *= max_ddim; /* Temporary (until the global matrix is not
                                  defined by block) */

    BFT_MALLOC(eqa->row->col_g_id, n_max_cw_dofs, cs_gnum_t);
    BFT_MALLOC(eqa->row->col_idx, n_max_cw_dofs, int);
    BFT_MALLOC(eqa->row->expval, max_ddim*n_max_cw_dofs, cs_real_t);
  }

  return eqa;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_assemble_t structure
 *
 * \param[in, out]  p_eqa    pointer to a structure pointer to be freed
 */
/*----------------------------------------------------------------------------*/

static void
_free_equation_assembler_struct(cs_equation_assemble_t  **p_eqa)
{
  if (*p_eqa == NULL)
    return;

  cs_equation_assemble_t  *eqa = *p_eqa;

  if (eqa->ddim > 1)
    BFT_FREE(eqa->row->expval);
  BFT_FREE(eqa->row->col_g_id);
  BFT_FREE(eqa->row->col_idx);
  BFT_FREE(eqa->row);

  BFT_FREE(eqa);
  *p_eqa = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and define a cs_matrix_assembler_t structure
 *
 * \param[in]  n_elts     number of elements
 * \param[in]  n_dofbyx   number of DoFs by element
 * \param[in]  x2x        pointer to a cs_adjacency_t structure
 * \param[in]  rs         pointer to a range set or NULL if sequential
 *
 * \return a pointer to a new allocated cs_matrix_assembler_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_assembler_t *
_build_matrix_assembler(cs_lnum_t                n_elts,
                        int                      n_dofbyx,
                        const cs_adjacency_t    *x2x,
                        const cs_range_set_t    *rs)
{
  cs_gnum_t  *grows = NULL, *gcols = NULL;

  /* The second paramter is set to "true" meaning that the diagonal is stored
     separately --> MSR storage */
  cs_matrix_assembler_t  *ma = cs_matrix_assembler_create(rs->l_range, true);

  /* First loop to count max size of the buffer */
  cs_lnum_t  max_size = 0;
  for (cs_lnum_t id = 0; id < n_elts; id++)
    max_size = CS_MAX(max_size, x2x->idx[id+1] - x2x->idx[id]);

  /* We increment max_size to take into account the diagonal entry */
  int  buf_size = n_dofbyx * n_dofbyx * (max_size + 1);
  BFT_MALLOC(grows, buf_size, cs_gnum_t);
  BFT_MALLOC(gcols, buf_size, cs_gnum_t);

  if (n_dofbyx == 1)  { /* Simplified version */

    for (cs_lnum_t row_id = 0; row_id < n_elts; row_id++) {

      const cs_gnum_t  grow_id = rs->g_id[row_id];
      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];

      /* Diagonal term is excluded in this connectivity. Add it "manually" */
      grows[0] = grow_id, gcols[0] = grow_id;

      /* Extra diagonal couples */
      for (cs_lnum_t j = start, i = 1; j < end; j++, i++) {
        grows[i] = grow_id;
        gcols[i] = rs->g_id[x2x->ids[j]];
      }

      cs_matrix_assembler_add_g_ids(ma, end - start + 1, grows, gcols);

    } /* Loop on entities */

  }
  else {

    for (cs_lnum_t row_id = 0; row_id < n_elts; row_id++) {

      const cs_lnum_t  start = x2x->idx[row_id];
      const cs_lnum_t  end = x2x->idx[row_id+1];
      const int  n_entries = (end - start + 1) * n_dofbyx * n_dofbyx;
      const cs_gnum_t  *grow_ids = rs->g_id + row_id*n_dofbyx;

      int shift = 0;

      /* Diagonal term is excluded in this connectivity. Add it "manually" */
      for (int dof_i = 0; dof_i < n_dofbyx; dof_i++) {
        const cs_gnum_t  grow_id = grow_ids[dof_i];
        for (int dof_j = 0; dof_j < n_dofbyx; dof_j++) {
          grows[shift] = grow_id;
          gcols[shift] = grow_ids[dof_j];
          shift++;
        }
      }

      /* Extra diagonal couples */
      for (cs_lnum_t j = start; j < end; j++) {

        const cs_lnum_t  col_id = x2x->ids[j];
        const cs_gnum_t  *gcol_ids = rs->g_id + col_id*n_dofbyx;

        for (int dof_i = 0; dof_i < n_dofbyx; dof_i++) {
          const cs_gnum_t  grow_id = grow_ids[dof_i];
          for (int dof_j = 0; dof_j < n_dofbyx; dof_j++) {
            grows[shift] = grow_id;
            gcols[shift] = gcol_ids[dof_j];
            shift++;
          }
        }

      } /* Loop on number of DoFs by entity */

      assert(shift == n_entries);
      cs_matrix_assembler_add_g_ids(ma, n_entries, grows, gcols);

    } /* Loop on entities */

  }

  /* Now compute structure */
  cs_matrix_assembler_compute(ma);

  /* Free temporary buffers */
  BFT_FREE(grows);
  BFT_FREE(gcols);

  return ma;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the pointer to a requested \ref cs_matrix_structure_t
 *         structure
 *
 * \param[in]  flag_id       id in the array of matrix structures
 *
 * \return  a pointer to a cs_matrix_structure_t
 */
/*----------------------------------------------------------------------------*/

cs_matrix_structure_t *
cs_equation_get_matrix_structure(int  flag)
{
  if (cs_equation_assemble_ms == NULL || flag < 0)
    return NULL;

  if (flag < CS_CDO_CONNECT_N_CASES)
    return cs_equation_assemble_ms[flag];
  else
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_equation_assemble_t structure related
 *         to a given thread
 *
 * \param[in]  t_id    id in the array of pointer
 *
 * \return a pointer to a cs_equation_assemble_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_assemble_t *
cs_equation_assemble_get(int    t_id)
{
  if (t_id < 0 || t_id >= cs_glob_n_threads)
    return NULL;

  return cs_equation_assemble[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize matrix-related structures according to
 *         the type of discretization used for this simulation
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  eb_flag      metadata for Edge-based schemes
 * \param[in]  fb_flag      metadata for Face-based schemes
 * \param[in]  vb_flag      metadata for Vertex-based schemes
 * \param[in]  vcb_flag     metadata for Vertex+Cell-basde schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_init(const cs_cdo_connect_t       *connect,
                          cs_flag_t                     eb_flag,
                          cs_flag_t                     fb_flag,
                          cs_flag_t                     vb_flag,
                          cs_flag_t                     vcb_flag,
                          cs_flag_t                     hho_flag)
{
  assert(connect != NULL); /* Sanity check */

  cs_timer_t t0, t1;
  CS_TIMER_COUNTER_INIT(cs_equation_ms_time);

  /* Two types of matrix assembler are considered:
   *  - The one related to matrix based on vertices
   *  - The one related to matrix based on faces
   */
  BFT_MALLOC(cs_equation_assemble_ma,
             CS_CDO_CONNECT_N_CASES, cs_matrix_assembler_t *);
  for (int i = 0; i < CS_CDO_CONNECT_N_CASES; i++)
    cs_equation_assemble_ma[i] = NULL;

  BFT_MALLOC(cs_equation_assemble_ms,
             CS_CDO_CONNECT_N_CASES, cs_matrix_structure_t *);
  for (int i = 0; i < CS_CDO_CONNECT_N_CASES; i++)
    cs_equation_assemble_ms[i] = NULL;

  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_lnum_t  n_edges = connect->n_edges;

  /* Allocate shared buffer and initialize shared structures
     Number max of DoFs at the cell level: n_max_cw_dofs
     Greatest dimension of the diagonal block: max_ddim
     Greatest dimension of the extra-diagonal block: max_edim
   */
  int  n_max_cw_dofs = 0, max_ddim = 1, max_edim = 1;

  /* Allocate and initialize matrix assembler and matrix structures */
  if (vb_flag > 0 || vcb_flag > 0) {

    const cs_adjacency_t  *v2v = connect->v2v;

    assert(v2v != NULL);
    n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_vbyc);

    if (vb_flag & CS_FLAG_SCHEME_SCALAR || vcb_flag & CS_FLAG_SCHEME_SCALAR) {

      t0 = cs_timer_time();

      /* Build the matrix structure and the matrix assembler structure */
      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

      cs_matrix_assembler_t  *ma = _build_matrix_assembler(n_vertices, 1, v2v,
                                                           rs);
      cs_matrix_structure_t  *ms
        = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      cs_equation_assemble_ma[CS_CDO_CONNECT_VTX_SCAL] = ma;
      cs_equation_assemble_ms[CS_CDO_CONNECT_VTX_SCAL] = ms;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

    } /* scalar-valued DoFs */

    if (vb_flag & CS_FLAG_SCHEME_VECTOR || vcb_flag & CS_FLAG_SCHEME_VECTOR) {

      t0 = cs_timer_time();

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_VECT];

      cs_matrix_assembler_t  *ma = _build_matrix_assembler(n_vertices, 3, v2v,
                                                           rs);
      cs_matrix_structure_t  *ms
        = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      cs_equation_assemble_ma[CS_CDO_CONNECT_VTX_VECT] = ma;
      cs_equation_assemble_ms[CS_CDO_CONNECT_VTX_VECT] = ms;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

      max_ddim = CS_MAX(max_ddim, 3);
      max_edim = CS_MAX(max_edim, 3);

    } /* vector-valued DoFs */

  } /* Vertex-based schemes and related ones */

  /* Allocate and initialize matrix assembler and matrix structures */
  if (eb_flag > 0) {

    const cs_adjacency_t  *e2e = connect->e2e;

    assert(e2e != NULL);
    n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_ebyc);

    if (eb_flag & CS_FLAG_SCHEME_SCALAR) {

      t0 = cs_timer_time();

      /* Build the matrix structure and the matrix assembler structure */
      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_EDGE_SCAL];

      cs_matrix_assembler_t  *ma = _build_matrix_assembler(n_edges, 1, e2e, rs);
      cs_matrix_structure_t  *ms
        = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

      cs_equation_assemble_ma[CS_CDO_CONNECT_EDGE_SCAL] = ma;
      cs_equation_assemble_ms[CS_CDO_CONNECT_EDGE_SCAL] = ms;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

    } /* scalar-valued DoFs */

  } /* Edge-based schemes and related ones */

  if (fb_flag > 0 || hho_flag > 0) {

    cs_matrix_structure_t  *ms0 = NULL, *ms1 = NULL, *ms2 = NULL;
    cs_matrix_assembler_t  *ma0 = NULL, *ma1 = NULL, *ma2 = NULL;

    n_max_cw_dofs = CS_MAX(n_max_cw_dofs, connect->n_max_fbyc);

    const cs_adjacency_t  *f2f = connect->f2f;

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR)) {

      t0 = cs_timer_time();

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

      ma0 = _build_matrix_assembler(n_faces, 1, f2f, rs);
      ms0 = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma0);

      cs_equation_assemble_ma[CS_CDO_CONNECT_FACE_SP0] = ma0;
      cs_equation_assemble_ms[CS_CDO_CONNECT_FACE_SP0] = ms0;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

    } /* Scalar-valued CDO-Fb or HHO-P0 */

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY1 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR)) {

      t0 = cs_timer_time();

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP1];

      ma1 = _build_matrix_assembler(n_faces, CS_N_FACE_DOFS_1ST, f2f, rs);
      ms1 = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma1);

      assert((CS_CDO_CONNECT_FACE_SP1 == CS_CDO_CONNECT_FACE_VP0) &&
             (CS_CDO_CONNECT_FACE_SP1 == CS_CDO_CONNECT_FACE_VHP0));

      cs_equation_assemble_ma[CS_CDO_CONNECT_FACE_SP1] = ma1;
      cs_equation_assemble_ms[CS_CDO_CONNECT_FACE_SP1] = ms1;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

      max_ddim = CS_MAX(max_ddim, 3);
      max_edim = CS_MAX(max_edim, 3);

    } /* Vector CDO-Fb or HHO-P1 or vector HHO-P0 */

    if (cs_flag_test(hho_flag,
                     CS_FLAG_SCHEME_POLY2 | CS_FLAG_SCHEME_SCALAR)) {

      t0 = cs_timer_time();

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP2];

      ma2 = _build_matrix_assembler(n_faces, CS_N_FACE_DOFS_2ND, f2f, rs);
      ms2 = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma2);

      cs_equation_assemble_ma[CS_CDO_CONNECT_FACE_SP2] = ma2;
      cs_equation_assemble_ms[CS_CDO_CONNECT_FACE_SP2] = ms2;

      t1 = cs_timer_time();
      cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

      max_ddim = CS_MAX(max_ddim, CS_N_FACE_DOFS_2ND);
      max_edim = CS_MAX(max_edim, CS_N_FACE_DOFS_2ND);

    }

    /* For vector equations and HHO */
    if (cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY1) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY2)) {

      if  (hho_flag & CS_FLAG_SCHEME_POLY1) {

        t0 = cs_timer_time();

        const cs_range_set_t  *rs
          = connect->range_sets[CS_CDO_CONNECT_FACE_VHP1];

        ma1 = _build_matrix_assembler(n_faces, 3*CS_N_FACE_DOFS_1ST, f2f, rs);
        ms1 = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma1);

        cs_equation_assemble_ma[CS_CDO_CONNECT_FACE_VHP1] = ma1;
        cs_equation_assemble_ms[CS_CDO_CONNECT_FACE_VHP1] = ms1;

        t1 = cs_timer_time();
        cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

        max_ddim = CS_MAX(max_ddim, 3*CS_N_FACE_DOFS_1ST);
        max_edim = CS_MAX(max_edim, 3*CS_N_FACE_DOFS_1ST);

      }
      else if  (hho_flag & CS_FLAG_SCHEME_POLY2) {

        t0 = cs_timer_time();

        const cs_range_set_t  *rs
          = connect->range_sets[CS_CDO_CONNECT_FACE_VHP2];

        ma2 = _build_matrix_assembler(n_faces, 3*CS_N_FACE_DOFS_2ND, f2f, rs);
        ms2 = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma2);

        cs_equation_assemble_ma[CS_CDO_CONNECT_FACE_VHP2] = ma2;
        cs_equation_assemble_ms[CS_CDO_CONNECT_FACE_VHP2] = ms2;

        t1 = cs_timer_time();
        cs_timer_counter_add_diff(&cs_equation_ms_time, &t0, &t1);

        /* 18 = 3*6 =  3*CS_N_FACE_DOFS_2ND */
        max_ddim = CS_MAX(max_ddim, 3*CS_N_FACE_DOFS_2ND);
        max_edim = CS_MAX(max_edim, 3*CS_N_FACE_DOFS_2ND);

      }

    }

  } /* Face-based schemes (CDO or HHO) */

  /* Common buffers for assemble usage */
  const int  n_threads = cs_glob_n_threads;
  BFT_MALLOC(cs_equation_assemble, n_threads, cs_equation_assemble_t *);
  for (int i = 0; i < n_threads; i++)
    cs_equation_assemble[i] = NULL;

#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
#pragma omp parallel
  {
    /* Each thread allocate its part. This yields a better memory affinity */
    int  t_id = omp_get_thread_num();
    cs_equation_assemble[t_id] = _init_equation_assembler_struct(max_ddim,
                                                                 max_edim,
                                                                 n_max_cw_dofs);
  }
#else
  cs_equation_assemble[0] = _init_equation_assembler_struct(max_ddim,
                                                            max_edim,
                                                            n_max_cw_dofs);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free matrix-related structures used during the simulation.
 *         Display overall statistic about the assembly stage for CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_finalize(void)
{
  /* Display profiling/performance information */
  cs_log_printf(CS_LOG_PERFORMANCE, " <CDO/Assembly> structure: %5.3e\n",
                cs_equation_ms_time.wall_nsec*1e-9);

  /* Free common assemble buffers */
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
#pragma omp parallel
  {
    int  t_id = omp_get_thread_num();
    cs_equation_assemble_t  *eqa = cs_equation_assemble[t_id];
    _free_equation_assembler_struct(&eqa);
  }
#else
  cs_equation_assemble_t  *eqa = cs_equation_assemble[0];
  _free_equation_assembler_struct(&eqa);
#endif
  BFT_FREE(cs_equation_assemble);

  /* Free matrix structures */
  for (int i = 0; i < CS_CDO_CONNECT_N_CASES; i++)
    cs_matrix_structure_destroy(&(cs_equation_assemble_ms[i]));
  BFT_FREE(cs_equation_assemble_ms);

  /* Free matrix assemblers */
  for (int i = 0; i < CS_CDO_CONNECT_N_CASES; i++)
    cs_matrix_assembler_destroy(&(cs_equation_assemble_ma[i]));
  BFT_FREE(cs_equation_assemble_ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the function pointer used to assemble the algebraic system
 *
 * \param[in] scheme     space discretization scheme
 * \param[in] ma_id      id in the array of matrix assembler
 *
 * \return a function pointer cs_equation_assembly_t
 */
/*----------------------------------------------------------------------------*/

cs_equation_assembly_t *
cs_equation_assemble_set(cs_param_space_scheme_t    scheme,
                         int                        ma_id)
{
  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    if (ma_id == CS_CDO_CONNECT_VTX_SCAL)
      return _set_scalar_assembly_func();
    else if (ma_id == CS_CDO_CONNECT_VTX_VECT)
      return _set_block33_assembly_func();
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    if (ma_id == CS_CDO_CONNECT_VTX_SCAL)
      return _set_scalar_assembly_func();
    break;

  case CS_SPACE_SCHEME_HHO_P0:
  case CS_SPACE_SCHEME_CDOFB:
    if (ma_id == CS_CDO_CONNECT_FACE_SP0)
      return _set_scalar_assembly_func();
    else if (ma_id == CS_CDO_CONNECT_FACE_VP0)
      return _set_block33_assembly_func();
    break;

  case CS_SPACE_SCHEME_HHO_P1:
  case CS_SPACE_SCHEME_HHO_P2:
    if (ma_id == CS_CDO_CONNECT_FACE_SP1)
      return _set_block33_assembly_func();
    else
      return _set_block_assembly_func();
    break;

  case CS_SPACE_SCHEME_CDOEB:
    if (ma_id == CS_CDO_CONNECT_EDGE_SCAL)
      return _set_scalar_assembly_func();
    break;

  default:
    return NULL; /* Case not handle */
  }

  return NULL;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Parallel and with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_mpit(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav)
{
  const cs_sdm_t  *const m = csys->mat;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */
  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_scal_dt(mav, ma, row);

    else {

      _assemble_row_scal_ld(ma, row);

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
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Parallel without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_mpis(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav)
{
  const cs_sdm_t  *const m = csys->mat;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */
  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    if (row->l_id < 0 || row->l_id >= rset->n_elts[0])
      _assemble_row_scal_ds(mav, ma, row);

    else {

      _assemble_row_scal_ld(ma, row);
      _add_scal_values_single(row, mav->matrix);

    }

  }

}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Sequential and with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_seqt(const cs_cell_sys_t             *csys,
                                 const cs_range_set_t            *rset,
                                 cs_equation_assemble_t          *eqa,
                                 cs_matrix_assembler_values_t    *mav)
{
  const cs_sdm_t  *const m = csys->mat;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */
  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    _assemble_row_scal_l(ma, row);

#if CS_CDO_OMP_SYNC_SECTIONS > 0 /* OpenMP with critical section */
    _add_scal_values_critical(row, mav->matrix);
#else
    _add_scal_values_atomic(row, mav->matrix);
#endif

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system
 *         Scalar-valued case. Sequential and without openMP.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to a matrix assembler buffers
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_matrix_seqs(const cs_cell_sys_t              *csys,
                                 const cs_range_set_t             *rset,
                                 cs_equation_assemble_t           *eqa,
                                 cs_matrix_assembler_values_t     *mav)
{
  const cs_sdm_t  *const m = csys->mat;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

  /* Push each row of the cellwise matrix into the assembler */
  for (int i = 0; i < row->n_cols; i++) {

    row->i = i;                               /* cellwise numbering */
    row->g_id = row->col_g_id[i];             /* global numbering */
    row->l_id = row->g_id - rset->l_range[0]; /* range set numbering */
    row->val = m->val + i*row->n_cols;

    _assemble_row_scal_l(ma, row);
    _add_scal_values_single(row, mav->matrix);

  } /* Loop on rows */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_seqs(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + csys->n_dofs,
                          row->expval + 2*csys->n_dofs };

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

      _assemble_row_scal_l(ma, row);
      _add_scal_values_single(row, mav->matrix);

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Sequential run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_seqt(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + csys->n_dofs,
                          row->expval + 2*csys->n_dofs };

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

      _assemble_row_scal_l(ma, row);

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
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_mpis(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + csys->n_dofs,
                          row->expval + 2*csys->n_dofs };

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

        _assemble_row_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block 3x3 entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock33_matrix_mpit(const cs_cell_sys_t           *csys,
                                          const cs_range_set_t          *rset,
                                          cs_equation_assemble_t        *eqa,
                                          cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim >= 3);
  assert(row->expval != NULL);

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vxyz[3] = {row->expval,
                          row->expval + csys->n_dofs,
                          row->expval + 2*csys->n_dofs};

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

        _assemble_row_scal_ld(ma, row);

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
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Sequential run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_seqs(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim < 19);
  assert(eqa->edim == eqa->ddim);
  assert(row->expval != NULL);

  const int  dim = eqa->ddim;

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*csys->n_dofs;

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

      _assemble_row_scal_l(ma, row);
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
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_seqt(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim < 19);
  assert(eqa->edim == eqa->ddim);
  assert(row->expval != NULL);

  const int  dim = eqa->ddim;

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*csys->n_dofs;

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

      _assemble_row_scal_l(ma, row);

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
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run without openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_mpis(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim < 19);
  assert(eqa->edim == eqa->ddim);
  assert(row->expval != NULL);

  const int  dim = eqa->ddim;

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*csys->n_dofs;

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

        _assemble_row_scal_ld(ma, row);
        _add_scal_values_single(row, mav->matrix);

      }

    } /* Push each row of the block */

  } /* Loop on row-wise blocks */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system into the global algebraic system.
 *         Case of a block NxN entries. Expand each row.
 *         Parallel run with openMP threading.
 *
 * \param[in]      csys     cellwise view of the algebraic system
 * \param[in]      rset     pointer to a cs_range_set_t structure
 * \param[in, out] eqa      pointer to an equation assembly structure
 * \param[in, out] mav      pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_eblock_matrix_mpit(const cs_cell_sys_t           *csys,
                                        const cs_range_set_t          *rset,
                                        cs_equation_assemble_t        *eqa,
                                        cs_matrix_assembler_values_t  *mav)
{
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_equation_assemble_row_t  *row = eqa->row;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(eqa->ddim < 19);
  assert(eqa->edim == eqa->ddim);
  assert(row->expval != NULL);

  const int  dim = eqa->ddim;

  /* Expand the values for a bundle of rows */
  cs_real_t  *_vpointer[18];

  for (int k = 0; k < dim; k++)
    _vpointer[k] = row->expval + k*csys->n_dofs;

  assert(m->n_rows == csys->n_dofs);
  row->n_cols = m->n_rows;

  /* Switch to the global numbering */
  for (int i = 0; i < row->n_cols; i++)
    row->col_g_id[i] = rset->g_id[csys->dof_ids[i]];

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

        _assemble_row_scal_ld(ma, row);

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
