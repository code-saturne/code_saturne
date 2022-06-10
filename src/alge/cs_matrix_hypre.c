/*============================================================================
 * Sparse Matrix Representation and Operations using HYPRE library.
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * HYPRE headers
 *----------------------------------------------------------------------------*/

/* Avoid warnings due to previous values */
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_URL
#undef PACKAGE_VERSION

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_utilities.h>

#if !defined(HYPRE_RELEASE_NUMBER)
#define HYPRE_RELEASE_NUMBER 0
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_base_accel.h"
#include "cs_matrix_default.h"
#include "cs_matrix_hypre.h"
#include "cs_matrix_hypre_priv.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*! \file cs_matrix_hypre.c
 *
 * \brief Sparse Matrix Representation and Operations using HYPRE.
 */
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static const char _hypre_ij_type_name[] = "HYPRE_PARCSR";
static const char _hypre_ij_type_fullname[] = "HYPRE IJ (HYPRE_ParCSR)";

static int _device_is_setup = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ensure GPU device is setup if needed.
 *
 * parameters:
 *   use_device  <-- 1 if device is used.
 *----------------------------------------------------------------------------*/

static void
_ensure_device_setup(int  use_device)
{
  if (_device_is_setup != use_device) {

#if HYPRE_RELEASE_NUMBER >=  22100
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_GPU)

    if (use_device == 1) {

      HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE); /* setup AMG on GPUs */
      HYPRE_SetSpGemmUseCusparse(0);               /* use hypre's SpGEMM
                                                      instead of cuSPARSE */
      HYPRE_SetUseGpuRand(1);                      /* use GPU RNG */

#     if defined(HYPRE_USING_CUDA) && defined(HYPRE_USING_DEVICE_POOL)
#     if defined(HYPRE_USING_UMPIRE)
      HYPRE_SetUmpireUMPoolName("HYPRE_UM_POOL_CODE_SATURNE");
      HYPRE_SetUmpireDevicePoolName("HYPRE_DEVICE_POOL_CODE_SATURNE");
      #else
      /* HYPRE_SetGPUMemoryPoolSize(bin_growth,
         min_bin, max_bin, max_bytes); */
#     endif
#     endif

    }
    else {
      HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
    }

#endif /* defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_GPU) */
#endif /* HYPRE_RELEASE_NUMBER >=  22100 */

    _device_is_setup = use_device;

  }
}

/*----------------------------------------------------------------------------
 * Matrix.vector product y = A.x with HYPRE matrix
 *
 * Note that since this function creates vectors, it may have significant
 * overhead. Since it is present for completednes, this should not be an issue.
 * If used more often, vectors could be maintained in the coeffs structure
 * along with the matrix.
 *
 * parameters:
 *   matrix       <-- pointer to matrix structure
 *   exclude_diag <-- exclude diagonal if true
 *   sync         <-- synchronize ghost cells if true
 *   x            <-> multipliying vector values
 *   y            --> resulting vector
 *----------------------------------------------------------------------------*/

static void
_mat_vec_p_parcsr(const cs_matrix_t  *matrix,
                  bool                exclude_diag,
                  bool                sync,
                  cs_real_t          *restrict x,
                  cs_real_t          *restrict y)
{
  assert(exclude_diag == false);

  const cs_lnum_t  n_rows = matrix->n_rows * matrix->db_size;
  const cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Get pointers to structures through coefficients,
     and copy input values */

  HYPRE_ParCSRMatrix p_a;
  HYPRE_IJMatrixGetObject(coeffs->hm, (void **)&p_a);

  HYPRE_Real *_t = NULL;

  cs_alloc_mode_t  amode = CS_ALLOC_HOST;
  if (coeffs->memory_location != HYPRE_MEMORY_HOST)
    amode = CS_ALLOC_HOST_DEVICE_SHARED;

  if (sizeof(cs_real_t) == sizeof(HYPRE_Real) && amode == CS_ALLOC_HOST) {
    HYPRE_IJVectorSetValues(coeffs->hx, n_rows, NULL, x);
  }
  else {
    CS_MALLOC_HD(_t, n_rows, HYPRE_Real, amode);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      _t[ii] = x[ii];;
    }
    HYPRE_IJVectorSetValues(coeffs->hx, n_rows, NULL, _t);
  }
  if (sync)
    HYPRE_IJVectorAssemble(coeffs->hx);

  HYPRE_ParVector p_x, p_y;
  HYPRE_IJVectorGetObject(coeffs->hx, (void **)&p_x);
  HYPRE_IJVectorGetObject(coeffs->hy, (void **)&p_y);

  /* SpMv operation */

  HYPRE_ParCSRMatrixMatvec(1.0, p_a, p_x, 0, p_y);

  /* Copy data back */

  if (_t == NULL) {
    HYPRE_IJVectorGetValues(coeffs->hy, n_rows, NULL, y);
  }
  else {
    HYPRE_IJVectorGetValues(coeffs->hy, n_rows, NULL, _t);
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      y[ii] = _t[ii];
    }
    CS_FREE_HD(_t);
  }
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler(const cs_matrix_assembler_t   *ma,
                              HYPRE_Int                    **diag_sizes,
                              HYPRE_Int                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add = (cs_matrix_assembler_get_separate_diag(ma)) ? 1 : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  HYPRE_Int *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows, HYPRE_Int);
  BFT_MALLOC(_offdiag_sizes, n_rows, HYPRE_Int);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    HYPRE_Int n_r_diag = 0;
    HYPRE_Int n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows)
        n_r_diag++;
      else
        break;
    }

    _diag_sizes[i] = n_r_diag + n_diag_add;
    _offdiag_sizes[i] = n_cols - n_r_diag;
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 * with full diagonal blocks and A.I extradiangonal blocks fill.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   db_size       <-- diagonal block size
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler_db(const cs_matrix_assembler_t   *ma,
                                 cs_lnum_t                     db_size,
                                 HYPRE_Int                    **diag_sizes,
                                 HYPRE_Int                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add
    = (cs_matrix_assembler_get_separate_diag(ma)) ? db_size : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  HYPRE_Int *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows*db_size, HYPRE_Int);
  BFT_MALLOC(_offdiag_sizes, n_rows*db_size, HYPRE_Int);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    HYPRE_Int n_r_diag = 0;
    HYPRE_Int n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows) {
        n_r_diag++;
        if (col_ids[j] == i)
          n_r_diag += db_size -1;
      }
      else
        break;
    }

    for (cs_lnum_t j = 0; j < db_size; j++) {
      _diag_sizes[i*db_size + j] = n_r_diag + n_diag_add;
      _offdiag_sizes[i*db_size + j] = n_cols - n_r_diag;
    }
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of matrix entries using an assembler
 * with full blocks fill.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   ma            <-- pointer to matrix assembler structure
 *   b_size        <-- diagonal block size
 *   diag_sizes    --> diagonal values
 *   offdiag_sizes --> off-diagonal values
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_assembler_b(const cs_matrix_assembler_t   *ma,
                                cs_lnum_t                      b_size,
                                HYPRE_Int                    **diag_sizes,
                                HYPRE_Int                    **offdiag_sizes)
{
  const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);

  cs_lnum_t n_diag_add
    = (cs_matrix_assembler_get_separate_diag(ma)) ? 1 : 0;

  const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
  const cs_lnum_t *col_ids = cs_matrix_assembler_get_col_ids(ma);

  HYPRE_Int *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, n_rows*b_size, HYPRE_Int);
  BFT_MALLOC(_offdiag_sizes, n_rows*b_size, HYPRE_Int);

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t s_id = row_index[i];
    cs_lnum_t e_id = row_index[i+1];

    HYPRE_Int n_r_diag = 0;
    HYPRE_Int n_cols = e_id - s_id;

    for (cs_lnum_t j = s_id; j < e_id; j++) {
      if (col_ids[j] < n_rows) {
        n_r_diag++;
      }
      else
        break;
    }

    for (cs_lnum_t j = 0; j < b_size; j++) {
      _diag_sizes[i*b_size + j] = (n_r_diag + n_diag_add)*b_size;
      _offdiag_sizes[i*b_size + j] = (n_cols - n_r_diag)*b_size;
    }
  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------
 * Compute local and distant counts of a native matrix's entries.
 *
 * The caller is responsible for freeing the returned diag_sizes and
 * offdiag_sizes arrays.
 *
 * parameters:
 *   matrix        <-- pointer to matrix structure
 *   have_diag     <-- does the matrix include a diagonal ?
 *   n_edges       <-- local number of graph edges
 *   edges         <-- edges (symmetric row <-> column) connectivity
 *   g_e_id        <-- global element ids
 *   da            <-- diagonal values
 *   diag_sizes    --> diagonal sizes
 *   offdiag_sizes --> off-diagonal sizes
 *----------------------------------------------------------------------------*/

static void
_compute_diag_sizes_native(cs_matrix_t        *matrix,
                           bool                have_diag,
                           cs_lnum_t           n_edges,
                           const cs_lnum_t     edges[restrict][2],
                           const cs_gnum_t     g_e_id[],
                           const cs_real_t     da[restrict],
                           HYPRE_Int         **diag_sizes,
                           HYPRE_Int         **offdiag_sizes)
{
  cs_lnum_t  n_rows = matrix->n_rows;
  cs_lnum_t  b_size = matrix->db_size;
  cs_lnum_t  e_size = matrix->eb_size;
  cs_lnum_t  b_stride = b_size * b_size;

  cs_lnum_t _n_rows = n_rows*b_size;

  HYPRE_Int *_diag_sizes, *_offdiag_sizes;
  BFT_MALLOC(_diag_sizes, _n_rows, HYPRE_Int);
  BFT_MALLOC(_offdiag_sizes, _n_rows, HYPRE_Int);

  /* Case with b_size > e_size handled later */
  int n_diag = (have_diag && b_size == e_size) ? e_size : 0;

  cs_gnum_t g_id_lb = 0;
  cs_gnum_t g_id_ub = 0;

  if (_n_rows > 0) {
    g_id_lb = g_e_id[0];
    g_id_ub = g_e_id[n_rows-1] + 1;
  }

  /* Separate local and distant loops for better first touch logic */

# pragma omp parallel for  if(_n_rows > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < _n_rows; i++) {
    _diag_sizes[i] = n_diag;
    _offdiag_sizes[i] = 0;
  }

  const cs_lnum_t *group_index = NULL;
  if (matrix->numbering != NULL) {
    if (matrix->numbering->type == CS_NUMBERING_THREADS) {
      group_index = matrix->numbering->group_index;
    }
  }

  if (group_index != NULL) {

    const int n_threads = matrix->numbering->n_threads;
    const int n_groups = matrix->numbering->n_groups;

    for (int g_id = 0; g_id < n_groups; g_id++) {

#     pragma omp parallel for
      for (int t_id = 0; t_id < n_threads; t_id++) {

        for (cs_lnum_t edge_id = group_index[(t_id*n_groups + g_id)*2];
             edge_id < group_index[(t_id*n_groups + g_id)*2 + 1];
             edge_id++) {
          cs_lnum_t ii = edges[edge_id][0];
          cs_lnum_t jj = edges[edge_id][1];
          cs_gnum_t g_ii = g_e_id[ii];
          cs_gnum_t g_jj = g_e_id[jj];
          if (ii < n_rows) {
            if (g_jj >= g_id_lb && g_jj < g_id_ub)
              _diag_sizes[ii*b_size] += e_size;
            else
              _offdiag_sizes[ii*b_size] += e_size;
          }
          if (jj < n_rows) {
            if (g_ii >= g_id_lb && g_ii < g_id_ub)
              _diag_sizes[jj*b_size] += e_size;
            else
              _offdiag_sizes[jj*b_size] += e_size;
          }
        }
      }

    }

  }
  else {

    for (cs_lnum_t edge_id = 0; edge_id < n_edges; edge_id++) {
      cs_lnum_t ii = edges[edge_id][0];
      cs_lnum_t jj = edges[edge_id][1];
      cs_gnum_t g_ii = g_e_id[ii];
      cs_gnum_t g_jj = g_e_id[jj];
      if (ii < n_rows) {
        if (g_jj >= g_id_lb && g_jj < g_id_ub)
          _diag_sizes[ii*b_size] += e_size;
        else
          _offdiag_sizes[ii*b_size] += e_size;
      }
      if (jj < n_rows) {
        if (g_ii >= g_id_lb && g_ii < g_id_ub)
          _diag_sizes[jj*b_size] += e_size;
        else
          _offdiag_sizes[jj*b_size] += e_size;
      }
    }

  }

  /* Adjustment for "block" cases */

  if (b_size > 1) {

#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_rows; i++) {
      for (cs_lnum_t j = 1; j < b_size; j++) {
        _diag_sizes[i*b_size + j] = _diag_sizes[i*b_size];
        _offdiag_sizes[i*b_size + j] = _offdiag_sizes[i*b_size];
      }
    }

    /* Delayed diagonal terms for block diagonal */
    if (have_diag && e_size == 1) {
#     pragma omp parallel for  if(_n_rows > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_rows; i++) {
        for (cs_lnum_t j = 0; j < b_size; j++) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            cs_real_t a = da[i*b_stride + j*b_size + k];
            if (a < -0.0 || a > 0.0)
              _diag_sizes[i*b_size + j] += 1;
          }
        }
      }
    }

  }

  *diag_sizes = _diag_sizes;
  *offdiag_sizes = _offdiag_sizes;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize coefficients context structure.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       use_device  -1 for automatic, 0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

static void
_setup_coeffs(cs_matrix_t  *matrix,
              int           use_device)
{
  if (use_device < 0)
    use_device = (cs_get_device_id() < 0) ? 0 : 1;

  _ensure_device_setup(use_device);  /* Setup on first pass or device change */

  if (matrix->coeffs == NULL) {
    cs_matrix_coeffs_hypre_t  *coeffs;
    BFT_MALLOC(coeffs, 1, cs_matrix_coeffs_hypre_t);
    memset(coeffs, 0, sizeof(cs_matrix_coeffs_hypre_t));
    coeffs->matrix_state = 0;

    if (use_device == 1)
      coeffs->memory_location = HYPRE_MEMORY_DEVICE;
    else
      coeffs->memory_location = HYPRE_MEMORY_HOST;

    coeffs->max_chunk_size = 0; /* Defined later */

    matrix->coeffs = coeffs;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for initialization of HYPRE matrix coefficients using
 *        local row ids and column indexes.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 * \param[in]       db_size   optional diagonal block sizes
 * \param[in]       eb_size   optional extra-diagonal block sizes
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_init(void        *matrix_p,
                       cs_lnum_t    db_size,
                       cs_lnum_t    eb_size)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  if (matrix->coeffs == NULL)
    _setup_coeffs(matrix, -1);

  /* Associated matrix assembler */

  const cs_matrix_assembler_t  *ma = matrix->assembler;

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Create HYPRE matrix */

  HYPRE_IJMatrix hm = coeffs->hm;

  if (coeffs->matrix_state == 0) {

    cs_alloc_mode_t  amode = CS_ALLOC_HOST;
    if (coeffs->memory_location != HYPRE_MEMORY_HOST)
      amode = CS_ALLOC_HOST_DEVICE_SHARED;

    /* We seem to have memory issues when allocating by parts
       on GPU, so prepare buffers to transfer in a single step
       (we will also need to delay transfer of smaller pieces).
       On CPU, use smaller array to avoid excess duplicate memory */

    if (amode == CS_ALLOC_HOST) {
      coeffs->max_chunk_size = 32768;
    }
    else {
      const cs_lnum_t n_rows = cs_matrix_assembler_get_n_rows(ma);
      const cs_lnum_t *row_index = cs_matrix_assembler_get_row_index(ma);
      cs_lnum_t nnz =   row_index[n_rows] * eb_size*eb_size
                      + n_rows * db_size*db_size;
      coeffs->max_chunk_size = nnz;
    }
    CS_MALLOC_HD(coeffs->row_buf, coeffs->max_chunk_size, HYPRE_BigInt, amode);
    CS_MALLOC_HD(coeffs->col_buf, coeffs->max_chunk_size, HYPRE_BigInt, amode);
    CS_MALLOC_HD(coeffs->val_buf, coeffs->max_chunk_size, HYPRE_Real, amode);

    MPI_Comm comm = cs_glob_mpi_comm;
    if (comm == MPI_COMM_NULL)
      comm = MPI_COMM_WORLD;

    const cs_gnum_t *l_range = cs_matrix_assembler_get_l_range(ma);

    HYPRE_BigInt b_size = db_size;
    HYPRE_BigInt ilower = b_size*l_range[0];
    HYPRE_BigInt iupper = b_size*l_range[1] - 1;

    HYPRE_IJMatrixCreate(comm,
                         ilower,
                         iupper,
                         ilower,
                         iupper,
                         &hm);

    coeffs->l_range[0] = l_range[0];
    coeffs->l_range[1] = l_range[1];

    coeffs->hm = hm;

    HYPRE_IJMatrixSetObjectType(hm, HYPRE_PARCSR);

    HYPRE_Int *diag_sizes = NULL, *offdiag_sizes = NULL;

    if (db_size == 1)
      _compute_diag_sizes_assembler(ma,
                                    &diag_sizes,
                                    &offdiag_sizes);
    else if (eb_size == 1)
      _compute_diag_sizes_assembler_db(ma,
                                      db_size,
                                      &diag_sizes,
                                      &offdiag_sizes);
    else
      _compute_diag_sizes_assembler_b(ma,
                                      db_size,
                                      &diag_sizes,
                                      &offdiag_sizes);

    HYPRE_IJMatrixSetDiagOffdSizes(hm, diag_sizes, offdiag_sizes);
    HYPRE_IJMatrixSetMaxOffProcElmts(hm, 0);

    BFT_FREE(diag_sizes);
    BFT_FREE(offdiag_sizes);

    HYPRE_IJMatrixSetOMPFlag(hm, 0);

    HYPRE_IJMatrixInitialize_v2(hm, coeffs->memory_location);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add HYPRE matrix coefficients using global row ids
 *        and column indexes, using intermediate copy for indexes and values.
 *
 * \param[in, out]  coeffs    HYPRE Matrix coefficients handler
 * \param[in]       n         number of values to add
 * \param[in]       b_size    associated data block size
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 *
 * \return  HYPRE error code
 */
/*----------------------------------------------------------------------------*/

static HYPRE_Int
_assembler_values_add_block_cc(cs_matrix_coeffs_hypre_t  *coeffs,
                               HYPRE_Int                  n,
                               HYPRE_Int                  b_size,
                               HYPRE_Int                  stride,
                               const cs_gnum_t            row_g_id[],
                               const cs_gnum_t            col_g_id[],
                               const cs_real_t            vals[])
{
  HYPRE_IJMatrix hm = coeffs->hm;
  assert(hm != NULL);

  HYPRE_BigInt h_b_size = b_size;
  HYPRE_BigInt l_b = coeffs->l_range[0];
  HYPRE_BigInt u_b = coeffs->l_range[1];

  HYPRE_BigInt *rows = coeffs->row_buf;
  HYPRE_BigInt *cols = coeffs->col_buf;
  HYPRE_Real *values = coeffs->val_buf;

  HYPRE_Int block_step = coeffs->max_chunk_size / stride;

  for (HYPRE_Int s_id = 0; s_id < n; s_id += block_step) {

    HYPRE_Int n_group = block_step;
    if (s_id + n_group > n)
      n_group = n - s_id;

    HYPRE_Int l = 0;

    for (HYPRE_Int i = 0; i < n_group; i++) {
      HYPRE_Int r_g_id = row_g_id[s_id + i];
      if (r_g_id >= l_b && r_g_id < u_b) {
        for (cs_lnum_t j = 0; j < b_size; j++) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[l*stride + j*b_size + k]
              = row_g_id[s_id + i]*h_b_size + (HYPRE_BigInt)j;
            cols[l*stride + j*b_size + k]
              = col_g_id[s_id + i]*h_b_size + (HYPRE_BigInt)k;
            values[l*stride + j*b_size + k]
              = vals[(s_id + i)*stride + j*b_size + k];
          }
        }
        l++;
      }
    }

    HYPRE_Int hypre_ierr
      = HYPRE_IJMatrixAddToValues(hm, l*stride,
                                  NULL, rows, cols, values);

    if (hypre_ierr != 0)
      return hypre_ierr;

  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add extradiagonal HYPRE matrix coefficients using global row ids
 *        and column indexes, for fill type CS_MATRIX_BLOCK_D,
 *        CS_MATRIX_BLOCK_D_66, CS_MATRIX_BLOCK_D_SYM.
 *
 * \param[in, out]  coeffs    HYPRE Matrix coefficients handler
 * \param[in]       n         number of values to add
 * \param[in]       b_size    associated data block size
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 *
 * \return  HYPRE error code
 */
/*----------------------------------------------------------------------------*/

static HYPRE_Int
_assembler_values_add_block_d_e(cs_matrix_coeffs_hypre_t  *coeffs,
                                HYPRE_Int                  n,
                                HYPRE_Int                  b_size,
                                const cs_gnum_t            row_g_id[],
                                const cs_gnum_t            col_g_id[],
                                const cs_real_t            vals[])
{
  HYPRE_IJMatrix hm = coeffs->hm;
  assert(hm != NULL);

  HYPRE_BigInt h_b_size = b_size;

  HYPRE_BigInt l_b = coeffs->l_range[0];
  HYPRE_BigInt u_b = coeffs->l_range[1];

  HYPRE_BigInt *rows = coeffs->row_buf;
  HYPRE_BigInt *cols = coeffs->col_buf;
  HYPRE_Real *values = coeffs->val_buf;

  HYPRE_Int block_step = coeffs->max_chunk_size / b_size;

  for (HYPRE_Int s_id = 0; s_id < n; s_id += block_step) {

    HYPRE_Int n_group = block_step;
    if (s_id + n_group > n)
      n_group = n - s_id;

    HYPRE_Int l = 0;

    for (HYPRE_Int i = 0; i < n_group; i++) {
      HYPRE_Int r_g_id = row_g_id[s_id + i];
      if (r_g_id >= l_b && r_g_id < u_b) {
        for (cs_lnum_t j = 0; j < b_size; j++) {
          rows[l*b_size + j] = row_g_id[s_id + i]*h_b_size + (HYPRE_BigInt)j;
          cols[l*b_size + j] = col_g_id[s_id + i]*h_b_size + (HYPRE_BigInt)j;
          values[l*b_size + j] = vals[s_id + i];
        }
        l++;
      }
    }

    HYPRE_Int hypre_ierr
      = HYPRE_IJMatrixAddToValues(hm, l*b_size,
                                  NULL, rows, cols, values);

    if (hypre_ierr != 0)
      return hypre_ierr;

  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function for addition to HYPRE matrix coefficients using
 *        local row ids and column indexes.
 *
 * This function can be used in all cases, including when
 *  sizeof(HYPRE_BigInt) != sizeof(cs_gnum_t)
 *  sizeof(HYPRE_Real) != sizeof(cs_real_t)
 *
 * Values whose associated row index is negative should be ignored;
 * Values whose column index is -1 are assumed to be assigned to a
 * separately stored diagonal. Other indexes should be valid.
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
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 * \param[in]       n         number of values to add
 * \param[in]       stride    associated data block size
 * \param[in]       row_g_id  associated global row ids
 * \param[in]       col_g_id  associated global column ids
 * \param[in]       vals      pointer to values (size: n*stride)
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_add_g(void             *matrix_p,
                        cs_lnum_t         n,
                        cs_lnum_t         stride,
                        const cs_gnum_t   row_g_id[],
                        const cs_gnum_t   col_g_id[],
                        const cs_real_t   vals[])
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  HYPRE_Int hypre_ierr = 0;
  HYPRE_Int nrows = n;
  HYPRE_Int b_size = matrix->db_size;

  HYPRE_Int max_chunk_size = coeffs->max_chunk_size;

  /* Scalar matrix
     ------------- */

  if (b_size == 1) {

    HYPRE_IJMatrix hm = coeffs->hm;
    assert(hm != NULL);

    HYPRE_BigInt l_b = coeffs->l_range[0];
    HYPRE_BigInt u_b = coeffs->l_range[1];

    HYPRE_BigInt *rows = coeffs->row_buf;
    HYPRE_BigInt *cols = coeffs->col_buf;
    HYPRE_Real *values = coeffs->val_buf;

    for (HYPRE_Int s_id = 0; s_id < nrows; s_id += max_chunk_size) {

        HYPRE_Int n_group = max_chunk_size;
        if (s_id + n_group > n)
          n_group = nrows - s_id;

        HYPRE_Int l = 0;

        for (HYPRE_Int i = 0; i < n_group; i++) {
          HYPRE_Int r_g_id = row_g_id[s_id + i];
          if (r_g_id >= l_b && r_g_id < u_b) {
            rows[l] = row_g_id[s_id + i];
            cols[l] = col_g_id[s_id + i];
            values[l] = vals[s_id + i];
            l++;
          }
        }

        hypre_ierr
          = HYPRE_IJMatrixAddToValues(hm, l, NULL, rows, cols, values);

        if (hypre_ierr != 0)
          break;

    }
  }

  /* Block matrix
     ------------ */

  else {

    /* Full blocks (including diagonal terms for diagonal fill) */

    if (   matrix->fill_type >= CS_MATRIX_BLOCK
        || row_g_id[0] == col_g_id[0])
      hypre_ierr = _assembler_values_add_block_cc(coeffs,
                                                  nrows,
                                                  b_size,
                                                  stride,
                                                  row_g_id,
                                                  col_g_id,
                                                  vals);

    /* Diagonal bloc extra-diagonal terms only */

    else if (matrix->fill_type >= CS_MATRIX_BLOCK_D)
      hypre_ierr = _assembler_values_add_block_d_e(coeffs,
                                                   nrows,
                                                   b_size,
                                                   row_g_id,
                                                   col_g_id,
                                                   vals);

  }

  if (hypre_ierr != 0) {
    char err_desc_buffer[64];
    HYPRE_DescribeError(hypre_ierr, err_desc_buffer);
    err_desc_buffer[63] = '\0';
    bft_error(__FILE__, __LINE__, 0,
              _("%s: error in HYPRE matrix assembly:\n  %s"),
              __func__, err_desc_buffer);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function to start the final assembly of matrix coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_begin(void  *matrix_p)
{
  CS_UNUSED(matrix_p);

  /* Note: this function is called once all coefficients have
   *       been added, and before assembly is finalized.
   *       It could be used in a threading or tasking context to signify
   *       assembly finalization can start, returning immediately
   *       so the calling task can continue working during this finalization */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function to complete the final assembly of matrix
 *        coefficients.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure.
 *
 * \param[in, out]  matrix_p  untyped pointer to matrix description structure
 */
/*----------------------------------------------------------------------------*/

static void
_assembler_values_end(void  *matrix_p)
{
  cs_matrix_t  *matrix = (cs_matrix_t *)matrix_p;

  if (matrix == NULL)
    return;

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;
  HYPRE_IJMatrix hm = coeffs->hm;

  HYPRE_IJMatrixAssemble(hm);

  if (coeffs->matrix_state == 0) {

    CS_FREE_HD(coeffs->row_buf);
    CS_FREE_HD(coeffs->col_buf);

    MPI_Comm comm = cs_glob_mpi_comm;
    if (comm == MPI_COMM_NULL)
      comm = MPI_COMM_WORLD;

    /* Create associated vectors here also to avoid repeated creation
       (and possible overhead) where used */

    const HYPRE_Int  n_off_proc = matrix->n_cols_ext - matrix->n_rows;
    const HYPRE_BigInt b_size = matrix->db_size;

    HYPRE_BigInt ilower = b_size*(coeffs->l_range[0]);
    HYPRE_BigInt iupper = b_size*(coeffs->l_range[1]) - 1;

    HYPRE_IJVectorCreate(comm, ilower, iupper, &(coeffs->hx));
    HYPRE_IJVectorSetObjectType(coeffs->hx, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(coeffs->hx, n_off_proc);

    HYPRE_IJVectorCreate(comm, ilower, iupper, &(coeffs->hy));
    HYPRE_IJVectorSetObjectType(coeffs->hy, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(coeffs->hy, n_off_proc);

    HYPRE_IJVectorInitialize_v2(coeffs->hx, coeffs->memory_location);
    HYPRE_IJVectorInitialize_v2(coeffs->hy, coeffs->memory_location);
  }

  /* Set stat flag */

  coeffs->matrix_state = 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a CSR matrix assembler values structure.
 *
 * The associated matrix's structure must have been created using
 * \ref cs_matrix_structure_create_from_assembler.
 *
 * \param[in, out]  matrix                 pointer to matrix structure
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 *
 * \return  pointer to initialized matrix assembler values structure;
 */
/*----------------------------------------------------------------------------*/

static cs_matrix_assembler_values_t *
_assembler_values_create_hypre(cs_matrix_t      *matrix,
                               const cs_lnum_t   diag_block_size,
                               const cs_lnum_t   extra_diag_block_size)
{
  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_create(matrix->assembler,
                                        false,
                                        diag_block_size,
                                        extra_diag_block_size,
                                        (void *)matrix,
                                        _assembler_values_init,
                                        NULL,
                                        _assembler_values_add_g,
                                        _assembler_values_begin,
                                        _assembler_values_end);

  return mav;
}

/*----------------------------------------------------------------------------
 * Set HYPRE ParCSR matrix coefficients for block-diagonal cases.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_ij_db(cs_matrix_t        *matrix,
                  bool                symmetric,
                  cs_lnum_t           n_edges,
                  const cs_lnum_t     edges[restrict][2],
                  const cs_real_t     da[restrict],
                  const cs_real_t     xa[restrict])
{
  bool direct_assembly = false;

  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  assert(n_rows > 0);

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Sizes and buffers */

  HYPRE_IJMatrix hm = coeffs->hm;
  HYPRE_BigInt h_b_size = matrix->db_size;

  cs_lnum_t b_size = matrix->db_size;
  cs_lnum_t b_stride = b_size * b_size;

  assert(b_stride == b_size*b_size);

  HYPRE_Int max_chunk_size = coeffs->max_chunk_size;

  HYPRE_BigInt *rows = coeffs->row_buf;
  HYPRE_BigInt *cols = coeffs->col_buf;
  HYPRE_Real *aij = coeffs->val_buf;

  /* Diagonal part
     ------------- */

  if (have_diag) {

    HYPRE_Int ic = 0;

    for (cs_lnum_t ii = 0; ii < n_rows; ii++) {

      for (cs_lnum_t j = 0; j < b_size; j++) {
        for (cs_lnum_t k = 0; k < b_size; k++) {
          cs_real_t a = da[ii*b_stride + j*b_size + k];
          if (a < -0.0 || a > 0.0) {
            rows[ic] = g_id[ii]*h_b_size + (HYPRE_BigInt)j;
            cols[ic] = g_id[ii]*h_b_size + (HYPRE_BigInt)k;
            aij[ic] = a;
            ic++;
          }
        }
      }
      if (ic + b_stride > max_chunk_size) {
        if (direct_assembly)
          HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
        else
          HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
        ic = 0;
      }

    }

    if (ic > 0) {
      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
    }

  }

  /* Extradiagonal part */

  HYPRE_Int max_b_chunk_size = max_chunk_size / b_size / 2;

  if (symmetric) {

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0, ec = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_b_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[ic*b_size + k] = g_ii*h_b_size + (HYPRE_BigInt)k;
            cols[ic*b_size + k] = g_jj*h_b_size + (HYPRE_BigInt)k;
            aij[ic*b_size + k] = xa[e_id];
          }
          ic += 1;
        }
        if (jj < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[ic*b_size + k] = g_jj*h_b_size + (HYPRE_BigInt)k;
            cols[ic*b_size + k] = g_ii*h_b_size + (HYPRE_BigInt)k;
            aij[ic*b_size + k] = xa[e_id];
          }
          ic += 1;
        }
        ec++;
      }

      s_e_id += ec;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic*b_size, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic*b_size, NULL, rows, cols, aij);
    }

  }
  else { /* non-symmetric variant */

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0, ec = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_b_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[ic*b_size + k] = g_ii*h_b_size + (HYPRE_BigInt)k;
            cols[ic*b_size + k] = g_jj*h_b_size + (HYPRE_BigInt)k;
            aij[ic*b_size + k] = xa[e_id*2];
          }
          ic++;
        }
        if (jj < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[ic*b_size + k] = g_jj*h_b_size + (HYPRE_BigInt)k;
            cols[ic*b_size + k] = g_ii*h_b_size + (HYPRE_BigInt)k;
            aij[ic*b_size + k] = xa[e_id*2+1];
          }
          ic++;
        }
        ec++;
      }

      s_e_id += ec;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic*b_size, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic*b_size, NULL, rows, cols, aij);

    }

  }
}

/*----------------------------------------------------------------------------
 * Set HYPRE ParCSR matrix coefficients for full block cases.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_ij_b(cs_matrix_t        *matrix,
                 bool                symmetric,
                 cs_lnum_t           n_edges,
                 const cs_lnum_t     edges[restrict][2],
                 const cs_real_t     da[restrict],
                 const cs_real_t     xa[restrict])
{
  bool direct_assembly = false;

  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  assert(n_rows > 0);

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Sizes and buffers */

  HYPRE_IJMatrix hm = coeffs->hm;
  HYPRE_BigInt h_b_size = matrix->db_size;

  cs_lnum_t b_size = matrix->db_size;
  cs_lnum_t b_stride = b_size * b_size;

  assert(b_stride == b_size*b_size);

  HYPRE_Int max_chunk_size = coeffs->max_chunk_size;

  HYPRE_BigInt *rows = coeffs->row_buf;
  HYPRE_BigInt *cols = coeffs->col_buf;
  HYPRE_Real *aij = coeffs->val_buf;

  /* Diagonal part
     ------------- */

  if (have_diag) {

    cs_lnum_t s_id = 0;
    HYPRE_Int max_b_chunk_size = max_chunk_size / b_stride;

    while (s_id < n_rows) {

      HYPRE_Int ic = 0;

      for (cs_lnum_t ii = s_id;
           ii < n_rows && ic < max_b_chunk_size;
           ii++) {
        for (cs_lnum_t j = 0; j < b_size; j++) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            rows[ic*b_stride + j*b_size + k]
              = g_id[ii]*h_b_size + (HYPRE_BigInt)j;
            cols[ic*b_stride + j*b_size + k]
              = g_id[ii]*h_b_size + (HYPRE_BigInt)k;
            aij[ic*b_stride + j*b_size + k] = da[ii*b_stride + j*b_size + k];
          }
        }
        ic++;
      }

      s_id += ic;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic*b_stride, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic*b_stride, NULL, rows, cols, aij);
    }

  }  /* End of diagonal block addition */

  /* Extradiagonal part */

  HYPRE_Int max_b_chunk_size = max_chunk_size / b_stride / 2;

  if (symmetric) {

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0, ec = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_b_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            for (cs_lnum_t l = 0; l < b_size; l++) {
              rows[ic*b_stride + k*b_size + l] = g_ii*h_b_size + (HYPRE_BigInt)k;
              cols[ic*b_stride + k*b_size + l] = g_jj*h_b_size + (HYPRE_BigInt)l;
              aij[ic*b_stride + k*b_size + l] = xa[e_id*b_stride + k*b_size + l];
            }
          }
          ic += 1;
        }
        if (jj < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            for (cs_lnum_t l = 0; l < b_size; l++) {
              rows[ic*b_stride + k*b_size + l] = g_jj*h_b_size + (HYPRE_BigInt)k;
              cols[ic*b_stride + k*b_size + l] = g_ii*h_b_size + (HYPRE_BigInt)l;
              aij[ic*b_stride + k*b_size + l] = xa[e_id*b_stride + k*b_size + l];
            }
          }
          ic += 1;
        }
        ec++;
      }

      s_e_id += ec;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic*b_stride, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic*b_stride, NULL, rows, cols, aij);
    }

  }
  else { /* non-symmetric variant */

    cs_lnum_t s_e_id = 0;

    while (s_e_id < n_edges) {

      HYPRE_Int ic = 0, ec = 0;

      for (cs_lnum_t e_id = s_e_id;
           e_id < n_edges && ic < max_b_chunk_size;
           e_id++) {
        cs_lnum_t ii = edges[e_id][0];
        cs_lnum_t jj = edges[e_id][1];
        cs_gnum_t g_ii = g_id[ii];
        cs_gnum_t g_jj = g_id[jj];
        if (ii < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            for (cs_lnum_t l = 0; l < b_size; l++) {
              rows[ic*b_stride + k*b_size + l] = g_ii*h_b_size + (HYPRE_BigInt)k;
              cols[ic*b_stride + k*b_size + l] = g_jj*h_b_size + (HYPRE_BigInt)l;
              aij[ic*b_stride + k*b_size + l]
                = xa[e_id*2*b_stride + k*b_size + l];
            }
          }
          ic += 1;
        }
        if (jj < n_rows) {
          for (cs_lnum_t k = 0; k < b_size; k++) {
            for (cs_lnum_t l = 0; l < b_size; l++) {
              rows[ic*b_stride + k*b_size + l] = g_jj*h_b_size + (HYPRE_BigInt)k;
              cols[ic*b_stride + k*b_size + l] = g_ii*h_b_size + (HYPRE_BigInt)l;
              aij[ic*b_stride + k*b_size + l]
                = xa[(e_id*2+1)*b_stride + k*b_size + l];
            }
          }
          ic += 1;
        }
        ec++;
      }

      s_e_id += ec;

      if (direct_assembly)
        HYPRE_IJMatrixSetValues(hm, ic*b_stride, NULL, rows, cols, aij);
      else
        HYPRE_IJMatrixAddToValues(hm, ic*b_stride, NULL, rows, cols, aij);

    }

  }
}

/*----------------------------------------------------------------------------
 * Set HYPRE ParCSR matrix coefficients.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_set_coeffs_ij(cs_matrix_t        *matrix,
               bool                symmetric,
               bool                copy,
               cs_lnum_t           n_edges,
               const cs_lnum_t     edges[restrict][2],
               const cs_real_t     da[restrict],
               const cs_real_t     xa[restrict])
{
  CS_UNUSED(copy);

  bool direct_assembly = false;

  const cs_lnum_t  n_rows = matrix->n_rows;

  const cs_gnum_t *g_id = cs_matrix_get_block_row_g_id(matrix);
  const bool have_diag = (xa != NULL) ? true : false;

  MPI_Comm comm = cs_glob_mpi_comm;
  if (comm == MPI_COMM_NULL)
    comm = MPI_COMM_WORLD;

  assert(n_rows > 0);

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  /* Create HYPRE matrix */

  HYPRE_IJMatrix hm = coeffs->hm;

  cs_gnum_t  l_range[2] = {0, 0};

  if (n_rows > 0) {
    l_range[0] = g_id[0];
    l_range[1] = g_id[n_rows-1] + 1;
  }

  cs_lnum_t b_size = matrix->db_size;
  cs_lnum_t e_size = matrix->eb_size;
  cs_lnum_t b_stride = b_size * b_size;
  cs_lnum_t e_stride = e_size * e_size;

  assert(b_stride == b_size*b_size);

  if (coeffs->matrix_state == 0) {

    cs_alloc_mode_t  amode = CS_ALLOC_HOST;
    if (coeffs->memory_location != HYPRE_MEMORY_HOST)
      amode = CS_ALLOC_HOST_DEVICE_SHARED;

    /* We seem to have memory issues when allocating by parts
       on GPU, so prepare buffers to transfer in a single step.
       On CPU, use smaller array to avoid excess duplicate memory */

    if (amode == CS_ALLOC_HOST) {
      coeffs->max_chunk_size = 32768;
    }
    else {
      cs_lnum_t nnz = n_edges*2*e_stride + n_rows*b_stride;
      coeffs->max_chunk_size = nnz;
    }

    CS_MALLOC_HD(coeffs->row_buf, coeffs->max_chunk_size, HYPRE_BigInt, amode);
    CS_MALLOC_HD(coeffs->col_buf, coeffs->max_chunk_size, HYPRE_BigInt, amode);
    CS_MALLOC_HD(coeffs->val_buf, coeffs->max_chunk_size, HYPRE_Real, amode);

    HYPRE_BigInt ilower = b_size*l_range[0];
    HYPRE_BigInt iupper = b_size*l_range[1] - 1;

    coeffs->l_range[0] = l_range[0];
    coeffs->l_range[1] = l_range[1];

    HYPRE_IJMatrixCreate(comm,
                         ilower,
                         iupper,
                         ilower,
                         iupper,
                         &hm);

    coeffs->hm = hm;

    HYPRE_IJMatrixSetObjectType(hm, HYPRE_PARCSR);

    HYPRE_Int *diag_sizes = NULL, *offdiag_sizes = NULL;

    _compute_diag_sizes_native(matrix,
                               have_diag,
                               n_edges,
                               edges,
                               g_id,
                               da,
                               &diag_sizes,
                               &offdiag_sizes);

    HYPRE_IJMatrixSetDiagOffdSizes(hm, diag_sizes, offdiag_sizes);
    HYPRE_IJMatrixSetMaxOffProcElmts(hm, 0);

    BFT_FREE(diag_sizes);
    BFT_FREE(offdiag_sizes);

    HYPRE_IJMatrixSetOMPFlag(hm, 0);

    HYPRE_IJMatrixInitialize_v2(hm, coeffs->memory_location);
  }

  HYPRE_Int max_chunk_size = coeffs->max_chunk_size / 2;

  HYPRE_BigInt *rows = coeffs->row_buf;
  HYPRE_BigInt *cols = coeffs->col_buf;
  HYPRE_Real *aij = coeffs->val_buf;

  /* Scalar case */

  if (b_size == 1) {

    if (have_diag) {

      cs_lnum_t s_id = 0;

      while (s_id < n_rows) {
        HYPRE_Int ic = 0;

        for (cs_lnum_t ii = s_id;
             ii < n_rows && ic < max_chunk_size;
             ii++) {
          rows[ic] = g_id[ii];
          cols[ic] = g_id[ii];
          aij[ic] = da[ii];
          ic++;
        }

        s_id += ic;

        if (direct_assembly)
          HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
        else
          HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
      }

    }

    if (symmetric) {

      cs_lnum_t s_e_id = 0;

      while (s_e_id < n_edges) {

        HYPRE_Int ic = 0, ec = 0;

        for (cs_lnum_t e_id = s_e_id;
             e_id < n_edges && ic < max_chunk_size;
             e_id++) {
          cs_lnum_t ii = edges[e_id][0];
          cs_lnum_t jj = edges[e_id][1];
          cs_gnum_t g_ii = g_id[ii];
          cs_gnum_t g_jj = g_id[jj];
          if (ii < n_rows) {
            rows[ic] = g_ii;
            cols[ic] = g_jj;
            aij[ic] = xa[e_id];
            ic++;
          }
          if (jj < n_rows) {
            rows[ic] = g_jj;
            cols[ic] = g_ii;
            aij[ic] = xa[e_id];
            ic++;
          }
          ec++;
        }

        s_e_id += ec;

        if (direct_assembly)
          HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
        else
          HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);
      }

    }
    else { /* non-symmetric variant */

      cs_lnum_t s_e_id = 0;

      while (s_e_id < n_edges) {

        HYPRE_Int ic = 0, ec = 0;

        for (cs_lnum_t e_id = s_e_id;
             e_id < n_edges && ic < max_chunk_size;
             e_id++) {
          cs_lnum_t ii = edges[e_id][0];
          cs_lnum_t jj = edges[e_id][1];
          cs_gnum_t g_ii = g_id[ii];
          cs_gnum_t g_jj = g_id[jj];

          if (ii < n_rows) {
            rows[ic] = g_ii;
            cols[ic] = g_jj;
            aij[ic] = xa[e_id*2];
            ic++;
          }
          if (jj < n_rows) {
            rows[ic] = g_jj;
            cols[ic] = g_ii;
            aij[ic] = xa[e_id*2+1];
            ic++;
          }
          ec++;
        }

        s_e_id += ec;

        if (direct_assembly)
          HYPRE_IJMatrixSetValues(hm, ic, NULL, rows, cols, aij);
        else
          HYPRE_IJMatrixAddToValues(hm, ic, NULL, rows, cols, aij);

      }

    }

  }

  /* Block diagonal only */

  else if (b_size > 1 && e_size == 1)
    _set_coeffs_ij_db(matrix,
                      symmetric,
                      n_edges,
                      edges,
                      da,
                      xa);

  /* Full block */

  else /* if (b_size > 1 && e_size > 1) */
    _set_coeffs_ij_b(matrix,
                     symmetric,
                     n_edges,
                     edges,
                     da,
                     xa);

  /* Finalize asembly and update state */

  _assembler_values_end(matrix);
}

/*----------------------------------------------------------------------------
 * Release HYPRE ParCSR matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_release_coeffs_ij(cs_matrix_t  *matrix)
{
  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  if (matrix->coeffs != NULL) {

    if (coeffs->matrix_state > 0) {
      HYPRE_IJMatrixDestroy(coeffs->hm);

      HYPRE_IJVectorDestroy(coeffs->hx);
      HYPRE_IJVectorDestroy(coeffs->hy);

      CS_FREE_HD(coeffs->row_buf); /* precaution; usually done earlier */
      CS_FREE_HD(coeffs->col_buf); /* precaution; usually done earlier */
      CS_FREE_HD(coeffs->val_buf); /* only done here */

      coeffs->matrix_state = 0;
    }
  }
}

/*----------------------------------------------------------------------------
 * Release HYPRE ParCSR matrix coefficients.
 *
 * Depending on current options and initialization, values will be copied
 * or simply mapped.
 *
 * parameters:
 *   matrix    <-- pointer to matrix structure
 *   symmetric <-- indicates if extradiagonal values are symmetric
 *   copy      <-- indicates if coefficients should be copied
 *   n_edges   <-- local number of graph edges
 *   edges     <-- edges (symmetric row <-> column) connectivity
 *   da        <-- diagonal values
 *   xa        <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_destroy_coeffs_ij(cs_matrix_t  *matrix)
{
  if (matrix->coeffs != NULL) {
    _release_coeffs_ij(matrix);
    BFT_FREE(matrix->coeffs);
  }
}

/*----------------------------------------------------------------------------
 * Copy diagonal of HYPRE ParCSR matrix.
 *
 * parameters:
 *   matrix <-- pointer to matrix structure
 *   da     --> diagonal (pre-allocated, size: n_rows)
 *----------------------------------------------------------------------------*/

static void
_copy_diagonal_ij(const cs_matrix_t  *matrix,
                  cs_real_t          *restrict da)
{
  const cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  const HYPRE_BigInt b_size = matrix->db_size;

  HYPRE_BigInt ilower = b_size*(coeffs->l_range[0]);

  const HYPRE_BigInt n_rows = matrix->n_rows * b_size;

  HYPRE_BigInt *rows, *cols;
  HYPRE_Int *n_rcols;
  HYPRE_Real *aij;

  BFT_MALLOC(n_rcols, n_rows, HYPRE_Int);
  BFT_MALLOC(rows, n_rows, HYPRE_BigInt);
  BFT_MALLOC(cols, n_rows, HYPRE_BigInt);

  if (sizeof(HYPRE_Real) == sizeof(cs_real_t))
    aij = da;
  else
    BFT_MALLOC(aij, n_rows, HYPRE_Real);

# pragma omp parallel for  if(n_rows > CS_THR_MIN)
  for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
    n_rcols[ii] = 1;
    rows[ii] = ilower + ii;;
    cols[ii] = ilower + ii;;
  }

  HYPRE_IJMatrixGetValues(coeffs->hm,
                          matrix->n_rows,
                          n_rcols,
                          rows,
                          cols,
                          aij);

  if (aij != da) {
    for (HYPRE_BigInt ii = 0; ii < n_rows; ii++) {
      da[ii] = aij[ii];;
    }
    BFT_FREE(aij);
  }

  BFT_FREE(cols);
  BFT_FREE(rows);
  BFT_FREE(n_rcols);
}

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief return coefficients structure associated with HYPRE matrix.
 *
 * \param[in]  matrix  pointer to matrix structure
 *
 * \return  pointer to matrix coefficients handler structure for HYPRE matrix.
 */
/*----------------------------------------------------------------------------*/

cs_matrix_coeffs_hypre_t *
cs_matrix_hypre_get_coeffs(const cs_matrix_t  *matrix)
{
  return (cs_matrix_coeffs_hypre_t *)(matrix->coeffs);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Switch matrix type to HYPRE.
 *
 * This releases previous coefficients if present, so should be called
 * just after matrix creation, before assigning coefficients.
 *
 * \param[in, out]  matrix      pointer to matrix structure
 * \param[in]       use_device  0 for host, 1 for device (GPU)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_set_type_hypre(cs_matrix_t  *matrix,
                         int           use_device)
{
  _ensure_device_setup(use_device);  /* Setup on first pass */

  matrix->type = CS_MATRIX_N_BUILTIN_TYPES;

  matrix->type_name = _hypre_ij_type_name;
  matrix->type_fname = _hypre_ij_type_fullname;

  /* Release previous coefficients if present */

  if (matrix->coeffs != NULL)
    matrix->destroy_coefficients(matrix);

  _setup_coeffs(matrix, use_device);

  cs_matrix_coeffs_hypre_t  *coeffs = matrix->coeffs;

  matrix->coeffs = coeffs;

  /* Set function pointers here */

  matrix->set_coefficients = _set_coeffs_ij;
  matrix->release_coefficients = _release_coeffs_ij;
  matrix->destroy_coefficients = _destroy_coeffs_ij;
  matrix->assembler_values_create = _assembler_values_create_hypre;

  matrix->get_diagonal = NULL;

  /* Remark: block values are transformed into scalar values, so SpMv products
     should be possible, (and the function pointers updated). HYPRE also seems
     to have support for block matrixes (hypre_ParCSRBlockMatrix)  but the
     high-level documentation does not mention it. */

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    matrix->vector_multiply[i][0] = _mat_vec_p_parcsr;
    matrix->vector_multiply[i][1] = NULL;
  }

  matrix->copy_diagonal = _copy_diagonal_ij;

  /* Force MPI initialization if not already done.
   * The main code_saturne communicator is not modifed, as
   * this is purely for external libraries use. */

#if defined(HAVE_MPI)

  if (cs_glob_mpi_comm == MPI_COMM_NULL) {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag == 0)
      MPI_Init(NULL, NULL);
  }

#endif

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
