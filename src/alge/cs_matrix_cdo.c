/*============================================================================
 * Sparse Matrix Representation and Operations.
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

#if defined (HAVE_MKL)
#include <mkl_spblas.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_sort.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_cdo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */


/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Create a CSR matrix structure from an index and an array related
 *        to column id
 *
 * \param[in] have_diag  indicates if the diagonal structure contains nonzeroes
 * \param[in] owner      deallocate members in structure only if owner=true
 * \param[in] n_rows     local number of rows
 * \param[in] n_cols     local number of cols
 * \param[in] n_cols_ext local number of cols + ghost cols implied in sync.
 * \param[in] idx        index on rows
 * \param[inout] col_id  array of colum ids related to the row index
 *
 * \return  a pointer to a created CSR matrix structure
 *----------------------------------------------------------------------------*/

static cs_matrix_struct_csr_t *
_create_struct_csr(bool            have_diag,
                   bool            owner,
                   cs_lnum_t       n_rows,
                   cs_lnum_t       n_cols,
                   cs_lnum_t      *idx,
                   cs_lnum_t      *col_id)
{
  cs_lnum_t ii, _n_cols, n_cols_max;

  cs_matrix_struct_csr_t  *ms = NULL;

  /* Allocate and map */

  BFT_MALLOC(ms, 1, cs_matrix_struct_csr_t);

  ms->n_rows = n_rows;
  ms->n_cols = n_cols;

  ms->direct_assembly = false; // Always the case with CDO schemes
  ms->have_diag = have_diag;

  ms->row_index = idx;

  /* Count the max. number of nonzero elements per row */
  n_cols_max = 0;
  for (ii = 0; ii < ms->n_rows; ii++) {
    _n_cols = idx[ii+1] - idx[ii];
    if (_n_cols > n_cols_max)
      n_cols_max = _n_cols;
  }

  ms->n_cols_max = n_cols_max;
  if (have_diag == false) // diagonal term is stored elsewhere
    ms->n_cols_max += 1;

  /* Sort line elements by column id (for better access patterns) */
  if (owner == true)
    if (n_cols_max > 1)
      for (ii = 0; ii < ms->n_rows; ii++)
        cs_sort_shell(idx[ii], idx[ii+1], col_id);

  ms->col_id = col_id;

  return ms;
}

/*----------------------------------------------------------------------------
 * Set CSR matrix coefficients.
 * If da and xa are equal to NULL, then initialize val with zeros.
 *
 * parameters:
 *   matrix           <-> Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values (NULL if all zero)
 *   xa               <-- Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_csr(cs_matrix_t      *matrix,
                bool              symmetric,
                bool              copy,
                const cs_real_t  *restrict da,
                const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii;
  cs_matrix_coeff_csr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Sanity check */
  assert(ms->direct_assembly == false);

  if (mc->val == NULL)
    BFT_MALLOC(mc->val, ms->row_index[ms->n_rows], cs_real_t);

  /* Mark diagonal values as not queried (mc->_d_val not changed) */

  mc->d_val = NULL;

  if (da != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" da != NULL. This case is not handled within the CDO"
                " framework."));

  if (xa != NULL)
    memcpy(mc->val, xa, ms->row_index[ms->n_rows]*sizeof(cs_real_t));
  else
    for (ii = 0; ii < ms->row_index[ms->n_rows]; ii++)
      mc->val[ii] = 0.0;

}

/*----------------------------------------------------------------------------
 * Set MSR matrix coefficients.
 *
 * parameters:
 *   matrix           <-> Pointer to matrix structure
 *   symmetric        <-- Indicates if extradiagonal values are symmetric
 *   copy             <-- Indicates if coefficients should be copied
 *   da               <-- Diagonal values (NULL if all zero)
 *   xa               <-- Extradiagonal values (NULL if all zero)
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_msr(cs_matrix_t      *matrix,
                bool              symmetric,
                bool              copy,
                const cs_real_t  *restrict da,
                const cs_real_t  *restrict xa)
{
  cs_lnum_t  ii, jj;
  cs_matrix_coeff_msr_t  *mc = matrix->coeffs;

  const cs_matrix_struct_csr_t  *ms = matrix->structure;

  /* Sanity check */
  assert(ms->direct_assembly == false);

  /* Map or copy diagonal values */

  if (da != NULL) {

    if (copy) {
      if (mc->_d_val == NULL || mc->max_db_size < matrix->db_size[3]) {
        BFT_REALLOC(mc->_d_val, matrix->db_size[3]*ms->n_rows, cs_real_t);
        mc->max_db_size = matrix->db_size[3];
      }
      memcpy(mc->_d_val, da, matrix->db_size[3]*sizeof(cs_real_t) * ms->n_rows);
      mc->d_val = mc->_d_val;
    }
    else
      mc->d_val = da;

  }
  else
    mc->d_val = NULL;

  /* Extradiagonal values */ //TODO with matrix->eb_size[3] > 1

  if (mc->x_val == NULL)
    BFT_MALLOC(mc->x_val, ms->row_index[ms->n_rows], cs_real_t);

  if (xa != NULL)
    memcpy(mc->x_val, xa, ms->row_index[ms->n_rows]*sizeof(cs_real_t));
  else
    for (ii = 0; ii < ms->row_index[ms->n_rows]; ii++)
      mc->x_val[ii] = 0.0;

}

/*----------------------------------------------------------------------------
 * Select the sparse matrix-vector product function to be used by a
 * matrix or variant for a given fill type.
 *
 * Currently, possible variant functions are:
 *
 *   CS_MATRIX_NATIVE  (all fill types)
 *     default
 *     standard
 *     3_3_diag        (for CS_MATRIX_33_BLOCK_D or CS_MATRIX_33_BLOCK_D_SYM)
 *     omp             (for OpenMP with compatible numbering)
 *     vector          (For vector machine with compatible numbering)
 *
 *   CS_MATRIX_CSR     (for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     prefetch
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_CSR_SYM (for CS_MATRIX_SCALAR_SYM)
 *     default
 *     standard
 *     mkl             (with MKL)
 *
 *   CS_MATRIX_MSR     (all fill types except CS_MATRIX_33_BLOCK)
 *     standard
 *     prefetch
 *     mkl             (with MKL, for CS_MATRIX_SCALAR or CS_MATRIX_SCALAR_SYM)
 *
 * parameters:
 *   m_type          <-- Matrix type
 *   numbering       <-- mesh numbering type, or NULL
 *   fill type       <-- matrix fill type to merge from
 *   ed_flag         <-- 0: with diagonal only, 1 exclude only; 2; both
 *   func_name       <-- function type name, or NULL for default
 *   vector_multiply <-> multiplication function array
 *   loop_length     <-> loop length array
 *
 * returns:
 *   0 for success, 1 for incompatible function, 2 for compatible
 *   function not available in current build
 *----------------------------------------------------------------------------*/

static int
_set_spmv_func(cs_matrix_type_t             m_type,
               const cs_numbering_t        *numbering,
               cs_matrix_fill_type_t        fill_type,
               int                          ed_flag,
               const char                  *func_name,
               cs_matrix_vector_product_t  *vector_multiply[][2],
               int                          loop_length[][2])
{
  int retcode = 1;
  int standard = 0;

  cs_matrix_vector_product_t *spmv[2] = {NULL, NULL};
  int l_length[2] = {0, 0};

  if (func_name == NULL)
    standard = 2;
  else if (!strcmp(func_name, "default"))
    standard = 2;
  else if (!strcmp(func_name, "standard"))
    standard = 1;

  switch(m_type) {

  case CS_MATRIX_CSR:

    switch(fill_type) {
    case CS_MATRIX_SCALAR:
      if (standard > 0) {
        spmv[0] = cs_matrix_vec_p_l_csr;
        spmv[1] = cs_matrix_vec_p_l_csr;
      }
      else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
        spmv[0] = cs_matrix_vec_p_l_csr_mkl;
        spmv[1] = cs_matrix_vec_p_l_csr_mkl;
#else
        retcode = 2;
#endif
      }
      break;
    default:
      break;
    }

    break;

 case CS_MATRIX_MSR:

    if (standard > 0) {
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
        spmv[0] = cs_matrix_vec_p_l_msr;
        spmv[1] = cs_matrix_vec_p_l_msr;
        break;
      default:
        break;
      }
    }
    else if (!strcmp(func_name, "mkl")) {
#if defined(HAVE_MKL)
      switch(fill_type) {
      case CS_MATRIX_SCALAR:
        spmv[0] = cs_matrix_vec_p_l_msr_mkl;
        spmv[1] = cs_matrix_vec_p_l_msr_mkl;
        break;
      default:
        break;
      }
#else
      retcode = 2;
#endif
    }
    break;

  default:
    break;

  } /* end of switch */

  if (ed_flag != 1 && spmv[0] != NULL) {
    vector_multiply[fill_type][0] = spmv[0];
    loop_length[fill_type][0] = l_length[0];
    retcode = 0;
  }
  if (ed_flag != 0 && spmv[0] != NULL) {
    vector_multiply[fill_type][1] = spmv[1];
    loop_length[fill_type][1] = l_length[1];
    retcode = 0;
  }

  return retcode;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * \brief Create a CDO matrix structure.
 *
 * Only CSR and MSR formats are handled.
 * col_id is sorted row by row during the creation of this structure.
 *
 * \param[in] type       type of matrix considered
 * \param[in] owner      deallocate members in structure only if owner=true
 * \param[in] have_diag  indicates if the diagonal structure contains nonzeroes
 * \param[in] n_rows     local number of rows
 * \param[in] n_rows_ext local number of rows + ghost rows implied in sync.
 * \param[in] n_cols     local number of cols
 * \param[in] idx        index on rows
 * \param[inout] col_id  array of colum ids related to the row index
 * \param[in] halo       halo structure for synchronization in parallel, or NULL
 * \param[in] numbering  vectorization or thread-related numbering info, or NULL
 *
 * \return  a pointer to a created CDO matrix structure
 *----------------------------------------------------------------------------*/

cs_matrix_cdo_structure_t *
cs_matrix_cdo_structure_create(cs_matrix_type_t       type,
                               bool                   owner,
                               bool                   have_diag,
                               cs_lnum_t              n_rows,
                               cs_lnum_t              n_rows_ext,
                               cs_lnum_t              n_cols,
                               cs_lnum_t             *idx,
                               cs_lnum_t             *col_id,
                               const cs_halo_t       *halo,
                               const cs_numbering_t  *numbering)
{
  cs_matrix_cdo_structure_t *ms = NULL;

  BFT_MALLOC(ms, 1, cs_matrix_cdo_structure_t);

  ms->type = type;

  ms->n_rows = n_rows;
  ms->n_rows_ext = n_rows_ext;
  ms->n_cols = n_cols;

  ms->owner = owner;

  /* Define Structure */

  switch(ms->type) {
  case CS_MATRIX_CSR:
    ms->structure = _create_struct_csr(have_diag,
                                       owner,
                                       n_rows,
                                       n_cols,
                                       idx,
                                       col_id);

    break;
  case CS_MATRIX_MSR:
    ms->structure = _create_struct_csr(false,
                                       owner,
                                       n_rows,
                                       n_cols,
                                       idx,
                                       col_id);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrices in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[type]));
    break;
  }

  /* Set pointers to structures shared from mesh here */

  ms->halo = halo;
  ms->numbering = numbering;

  return ms;
}

/*----------------------------------------------------------------------------
 * \brief Destroy a CDO matrix structure.
 *
 * \param[inout]  ms  pointer to a CDO matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_matrix_cdo_structure_destroy(cs_matrix_cdo_structure_t  **ms)
{
  if (ms != NULL && *ms != NULL) {

    cs_matrix_cdo_structure_t *_ms = *ms;

    if (_ms->owner) {
      switch(_ms->type) {

      case CS_MATRIX_CSR:
      case CS_MATRIX_MSR:
        {
          cs_matrix_struct_csr_t *structure = _ms->structure;
          cs_matrix_destroy_struct_csr(&structure);
        }
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  _("Handling of matrices in %s format\n"
                    "is not operational yet."),
                  _(cs_matrix_type_name[_ms->type]));
        break;

      }

      _ms->structure = NULL;

      /* Now free main structure */

      BFT_FREE(*ms);
      *ms = NULL;

    } /* structure type */
  } /* owner */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix container using a given CDO matrix structure.
 *
 * Note that the matrix container maps to the assigned structure,
 * so it must be destroyed before that structure.
 *
 * \param[in]  ms  an associated CDO matrix structure
 *
 * \return  pointer to created matrix structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_cdo_create(const cs_matrix_cdo_structure_t  *ms)
{
  int i;
  cs_matrix_fill_type_t mft;

  cs_matrix_t *m = NULL;

  BFT_MALLOC(m, 1, cs_matrix_t);

  m->type = ms->type;

  /* Map shared structure */

  m->n_cells = ms->n_rows;
  m->n_cells_ext = ms->n_rows_ext;
  m->n_faces = ms->n_cols;

  for (i = 0; i < 4; i++) {
    m->db_size[i] = 1;
    m->eb_size[i] = 1;
  }
  m->fill_type = CS_MATRIX_N_FILL_TYPES;

  m->structure = ms->structure;

  m->face_cell = NULL;
  m->cell_num = NULL;
  m->halo = ms->halo;
  m->numbering = ms->numbering;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
    for (i = 0; i < 2; i++) {
      m->loop_length[mft][i] = 0;
      m->vector_multiply[mft][i] = NULL;
    }
  }

  /* Define coefficients */

  m->xa = NULL;

  switch(m->type) {
  case CS_MATRIX_CSR:
    m->coeffs = cs_matrix_create_coeff_csr();
    break;
  case CS_MATRIX_MSR:
    m->coeffs = cs_matrix_create_coeff_msr();
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Handling of matrixes in %s format\n"
                "is not operational yet."),
              _(cs_matrix_type_name[m->type]));
    break;
  }

  /* Set function pointers here */

  m->set_coefficients = NULL;

  for (mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++)
    _set_spmv_func(m->type,
                   m->numbering,
                   mft,
                   2,    /* ed_flag */
                   NULL, /* func_name */
                   m->vector_multiply,
                   m->loop_length);

  switch(m->type) {

  case CS_MATRIX_CSR:
    m->set_coefficients = _set_coeffs_csr;
    m->release_coefficients = cs_matrix_release_coeffs_csr;
    m->copy_diagonal = cs_matrix_copy_diagonal_csr;
    break;

  case CS_MATRIX_MSR:
    m->set_coefficients = _set_coeffs_msr;
    m->release_coefficients = cs_matrix_release_coeffs_msr;
    m->copy_diagonal = cs_matrix_copy_diagonal_separate;
    break;

  default:
    assert(0);
    break;

  }

  for (i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    if (m->vector_multiply[i][1] == NULL)
      m->vector_multiply[i][1] = m->vector_multiply[i][0];
  }

  /* Additional query buffers here */

  m->row_buffer_size = 0;
  m->row_buffer_id = NULL;
  m->row_buffer_val = NULL;

  return m;
}

