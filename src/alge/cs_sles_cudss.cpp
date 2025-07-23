/*============================================================================
 * Sparse Linear Equation Solvers using Cudss
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

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
 * cuDSS headers
 *----------------------------------------------------------------------------*/

#include <cudss.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_log.h"
#include "base/cs_fp_exception.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "base/cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_sles.h"
#include "alge/cs_sles_cudss.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_cudss.cpp

  \brief handling of cuDSS-based linear solvers

  \page sles_cudss cuDSS-based linear solvers.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic per linear system options and logging */
/*---------------------------------------------*/

typedef struct _cs_sles_cudss_setup_t {

  bool                 need_update;

  cudssData_t          data;
  cudssMatrix_t        matrix;     /* Linear system matrix */
  cudssMatrix_t        rhs;
  cudssMatrix_t        x;

  int32_t              *_row_index;
  int32_t              *_col_id;

} cs_sles_cudss_setup_t;

struct _cs_sles_cudss_t {

  /* Options */

  bool                 keep_data;        /* If true, keep data even
                                            when cs_sles_free is called,
                                            to amortize analysis */

  /* Performance data */

  int                  n_analysis;         /* Number of system analysis */
  int                  n_num_fact;         /* Number of system setups */
  int                  n_solves;           /* Number of system solves */

  int64_t              mem_estimates[16];

  cs_timer_counter_t   t_sym_fact;         /* Total symbolic factorization */
  cs_timer_counter_t   t_num_fact;         /* Total numerical factorization */
  cs_timer_counter_t   t_solve;            /* Total time used */

  /* Setup data */

  cudssConfig_t         config;             /* configuration */

  cs_sles_cudss_setup_t   *setup_data;

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static int  _n_cudss_systems = 0;

cudssHandle_t  _handle = nullptr;  /* cuDSS instance handle */

static const char *_status_str[] = {
  "CUDSS_STATUS_SUCCESS",
  "CUDSS_STATUS_NOT_INITIALIZED",
  "CUDSS_STATUS_ALLOC_FAILED",
  "CUDSS_STATUS_INVALID_VALUE",
  "CUDSS_STATUS_NOT_SUPPORTED",
  "CUDSS_STATUS_EXECUTION_FAILED",
  "CUDSS_STATUS_INTERNAL_ERROR"
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

//*---------------------------------------------------------------------------*/
/*!
 * \brief Setup cuDSS matrix
 *
 * \param[in, out]  c   pointer to cuDSS solver info and context
 * \param[in]       a   associated matrix
 */
/*----------------------------------------------------------------------------*/

static void
_setup_matrix(cs_sles_cudss_t    *c,
              const cs_matrix_t  *a)
{
  cs_sles_cudss_setup_t *sd = c->setup_data;

  if (sizeof(cs_lnum_t) != sizeof(int))
    bft_error
      (__FILE__, __LINE__, 0,
       _("cuDSS bindings are not currently handled for code_saturne builds\n"
         "using long cs_lnumt_ types (i.e. --enable-long-lnum)."));

  const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
  const int n_rows = cs_matrix_get_n_rows(a);
  const int n_cols = cs_matrix_get_n_columns(a);
  const int db_size = cs_matrix_get_diag_block_size(a);

  const cs_lnum_t *a_row_index, *a_col_id;
  const cs_real_t *a_val = nullptr;

  cs_matrix_get_csr_arrays(a, &a_row_index, &a_col_id, &a_val);

  const void  *row_index = a_row_index, *col_id = a_col_id;

  cs_alloc_mode_t amode = CS_ALLOC_HOST_DEVICE;
  cs_lnum_t nnz = a_row_index[n_rows];

  // const cs_gnum_t *grow_id = cs_matrix_get_block_row_g_id(a);

  if (sizeof(int32_t) != sizeof(cs_lnum_t)) {
    CS_REALLOC_HD(sd->_row_index, n_rows, int32_t, amode);
    int32_t  *_row_index = sd->_row_index;
    for (cs_lnum_t i = 0; i < n_rows; i++)
      _row_index[i] = a_row_index[i];
    row_index = _row_index;
  }

  if (sizeof(int32_t) != sizeof(cs_lnum_t)) {
    CS_REALLOC_HD(sd->_col_id, nnz, int32_t, amode);
    int32_t  *_col_id = sd->_col_id;
    for (cs_lnum_t i = 0; i < nnz; i++)
      _col_id[i] = a_col_id[i];
    col_id = _col_id;
  }

  /* Matrix */

  constexpr cudaDataType_t
    value_type = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;
  cudssMatrixType_t  mtype = CUDSS_MTYPE_GENERAL;
  if (cs_matrix_is_symmetric(a))
    mtype = CUDSS_MTYPE_SPD; // CUDSS_MTYPE_SYMMETRIC;

  cudssStatus_t retval
    = cudssMatrixCreateCsr(&(sd->matrix),
                           n_rows,
                           n_cols,
                           nnz,
                           const_cast<void *>(row_index),
                           nullptr,
                           const_cast<void *>(col_id),
                           const_cast<cs_real_t *>(a_val),
                           CUDA_R_32I, // indexType,
                           value_type,
                           mtype,
                           CUDSS_MVIEW_FULL,
                           CUDSS_BASE_ZERO);

  if (retval != CUDSS_STATUS_SUCCESS)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: error calling cudssMatrixCreateCsr: %s\n"),
              __func__, _status_str[retval]);

  /* Assume partitioning with continuous blocks */

#if defined(HAVE_MPI)

  if (cs_glob_mpi_comm != MPI_COMM_NULL) {
    cs_gnum_t local_shift = n_rows;
    cs_gnum_t global_shift = 0;
    MPI_Scan(&local_shift, &global_shift, 1, CS_MPI_GNUM, MPI_SUM,
             cs_glob_mpi_comm);

    int64_t first_row = global_shift - n_rows;
    int64_t last_row = global_shift - 1;

    retval
      = cudssMatrixSetDistributionRow1d(sd->matrix,
                                        first_row,
                                        last_row);

    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssMatrixSetDistributionRow1d: %s\n"),
                __func__, _status_str[retval]);
  }

#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and associate an cuDSS linear system solver
 *        for a given field or equation name.
 *
 * If this system did not previously exist, it is added to the list of
 * "known" systems. Otherwise, its definition is replaced by the one
 * defined here.
 *
 * This is a utility function: if finer control is needed, see
 * \ref cs_sles_define and \ref cs_sles_cudss_create.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * Note that this function returns a pointer directly to the solver
 * management structure. This may be used to set further options.
 * If needed, \ref cs_sles_find may be used to obtain a pointer to the matching
 * \ref cs_sles_t container.
 *
 * \param[in]      f_id          associated field id, or < 0
 * \param[in]      name          associated name if f_id < 0, or nullptr
 *
 * \return  pointer to newly created cuDSS solver info object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_cudss_t *
cs_sles_cudss_define(int           f_id,
                     const char  *name)
{
  cs_sles_cudss_t * c = cs_sles_cudss_create();

  cs_sles_t *sc = cs_sles_define(f_id,
                                 name,
                                 c,
                                 "cs_sles_cudss_t",
                                 cs_sles_cudss_setup,
                                 cs_sles_cudss_solve,
                                 cs_sles_cudss_free,
                                 cs_sles_cudss_log,
                                 cs_sles_cudss_copy,
                                 cs_sles_cudss_destroy);

  CS_NO_WARN_IF_UNUSED(sc);

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create cuDSS linear system solver info and context.
 *
 * In case of rotational periodicity for a block (non-scalar) matrix,
 * the matrix type will be forced to MATSHELL ("shell") regardless
 * of the option used.
 *
 * \return  pointer to associated linear system object.
 */
/*----------------------------------------------------------------------------*/

cs_sles_cudss_t *
cs_sles_cudss_create(void)

{
  cs_sles_cudss_t *c;

  // if (_n_cudss_systems == 0);
  if (_handle == nullptr) {
    cudssStatus_t retval = cudssCreate(&_handle);
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssCreate: %s\n"),
                __func__, _status_str[retval]);

    if (cs_glob_n_ranks > 1) {
      // With file name set to null, cuDSS will use the file
      // defined through the CUDSS_COMM_LIB environment variable */
      const char* comm_lib_file_name = nullptr;
      retval = cudssSetCommLayer(_handle, comm_lib_file_name);
      if (retval != CUDSS_STATUS_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: error calling cudssSetCommLayer: %s\n\n"
                    "The CUDSS_COMM_LIB environment variable should be set to\n"
                    "the absolute path of the cuDSS communicator layer.\n"),
                  __func__, _status_str[retval]);
    }

    if (cs_glob_n_threads > 1) {
      // With file name set to null, cuDSS will use the file
      // defined through the CUDSS_THREADING_LIB environment variable */
      const char* thread_lib_file_name = nullptr;
      retval = cudssSetThreadingLayer(_handle, thread_lib_file_name);
      if (retval != CUDSS_STATUS_SUCCESS)
        cs_log_warning
          (_("%s: error calling cudssSetThreadingLayer: %s\n\n"
             "The CUDSS_THREADING_LIB environment variable should be set to\n"
             "the absolute path of the cuDSS threading layer.\n"),
           __func__, _status_str[retval]);
    }
  }

  _n_cudss_systems += 1;

  CS_MALLOC(c, 1, cs_sles_cudss_t);

  c->keep_data = true;

  c->n_analysis = 0;
  c->n_num_fact = 0;
  c->n_solves = 0;

  for (int i = 0; i < 16; i++)
    c->mem_estimates[i] = 0;

  CS_TIMER_COUNTER_INIT(c->t_sym_fact);
  CS_TIMER_COUNTER_INIT(c->t_num_fact);
  CS_TIMER_COUNTER_INIT(c->t_solve);

  /* Setup data */

  cudssConfigCreate(&(c->config));

  /* (optional) use hybrid mode */
  {
    int mode = 1;
    cudssConfigSet(c->config,
                   CUDSS_CONFIG_HYBRID_MODE,
                   &mode,
                   sizeof(int));

    int n_threads = cs_glob_n_threads;
    cudssConfigSet(c->config,
                   CUDSS_CONFIG_HOST_NTHREADS,
                   &n_threads,
                   sizeof(int));
  }

  /* (optional) Modify solver settings, e.g., reordering algorithm */
  cudssAlgType_t reorder_alg = CUDSS_ALG_DEFAULT;
  cudssConfigSet(c->config,
                 CUDSS_CONFIG_REORDERING_ALG,
                 &reorder_alg,
                 sizeof(cudssAlgType_t));

  c->setup_data = nullptr;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create cuDSS linear system solver info and context
 *        based on existing info and context.
 *
 * Most configuration parameters will be copied from the existing context,
 * though not all, since cuDSS does not provide a comprehensive way to
 * do this.
 *
 * \param[in]  context  pointer to reference info and context
 *                     (actual type: cs_sles_cudss_t  *)
 *
 * \return  pointer to newly created solver info object.
 *          (actual type: cs_sles_cudss_t  *)
 */
/*----------------------------------------------------------------------------*/

void *
cs_sles_cudss_copy(const void  *context)
{
  cs_sles_cudss_t *d = nullptr;

  if (context != nullptr) {
    const cs_sles_cudss_t *c = (const cs_sles_cudss_t *)context;
    d = cs_sles_cudss_create();
    size_t sw;

    cudssConfigParam_t int_opts[] = {
      CUDSS_CONFIG_HYBRID_MODE,
      CUDSS_CONFIG_HOST_NTHREADS,
      CUDSS_CONFIG_USE_MATCHING,
      CUDSS_CONFIG_SOLVE_MODE,
      CUDSS_CONFIG_IR_N_STEPS,
      CUDSS_CONFIG_USE_CUDA_REGISTER_MEMORY,
      CUDSS_CONFIG_HYBRID_EXECUTE_MODE,
      CUDSS_CONFIG_ND_NLEVELS
    };
    for (int i = 0; i < 8; i++) {
      int val = 0;
      cudssConfigGet(c->config, int_opts[i], &val, sizeof(int), &sw);
      cudssConfigSet(d->config, int_opts[i], &val, sizeof(int));
    };

    cudssConfigParam_t int64_opts[] = {
      CUDSS_CONFIG_HYBRID_DEVICE_MEMORY_LIMIT
    };
    for (int i = 0; i < 1; i++) {
      int64_t val = 0;
      cudssConfigGet(c->config, int64_opts[i], &val, sizeof(int64_t), &sw);
      cudssConfigSet(d->config, int64_opts[i], &val, sizeof(int64_t));
    };

    cudssConfigParam_t alg_opts[] = {
      CUDSS_CONFIG_REORDERING_ALG,
      CUDSS_CONFIG_FACTORIZATION_ALG,
      CUDSS_CONFIG_SOLVE_ALG,
      CUDSS_CONFIG_PIVOT_EPSILON_ALG,
      CUDSS_CONFIG_MATCHING_ALG,
      CUDSS_CONFIG_IR_TOL,
      CUDSS_CONFIG_HYBRID_MODE
    };
    for (int i = 0; i < 7; i++) {
      cudssAlgType_t val;
      cudssConfigGet(c->config, alg_opts[i], &val, sizeof(cudssAlgType_t), &sw);
      cudssConfigSet(d->config, alg_opts[i], &val, sizeof(cudssAlgType_t));
    };


    cudssConfigParam_t double_opts[] = {
      CUDSS_CONFIG_IR_TOL
    };

    for (int i = 0; i < 1; i++) {
      double val;
      cudssConfigGet(c->config, double_opts[i], &val, sizeof(double), &sw);
      cudssConfigSet(d->config, double_opts[i], &val, sizeof(double));
    }
  }

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy cuDSS linear system solver info and context.
 *
 * \param[in, out]  context  pointer to cuDSS solver info and context
 *                           (actual type: cs_sles_cudss_t  **)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_destroy(void **context)
{
  cs_sles_cudss_t *c = (cs_sles_cudss_t *)(*context);
  if (c != nullptr) {

    cudssConfigDestroy(c->config);

    /* Free structure */

    c->keep_data = false;

    cs_sles_cudss_free(c);
    CS_FREE(c);
    *context = c;

    _n_cudss_systems -= 1;
    if (_n_cudss_systems == 0) {
      cudssDestroy(_handle);
      _handle = nullptr;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup cuDSS linear equation solver.
 *
 * \param[in, out]  context    pointer to cuDSS solver info and context
 *                             (actual type: cs_sles_cudss_t  *)
 * \param[in]       name       pointer to system name
 * \param[in]       a          associated matrix
 * \param[in]       verbosity  associated verbosity
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_setup(void               *context,
                    const char         *name,
                    const cs_matrix_t  *a,
                    int                 verbosity)
{
  CS_NO_WARN_IF_UNUSED(verbosity);

  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_cudss_t  *c = (cs_sles_cudss_t *)context;
  cs_sles_cudss_setup_t *sd = c->setup_data;

  cudssStatus_t retval = CUDSS_STATUS_SUCCESS;

  /* Case where setup data is not currently present */

  if (sd == nullptr) {

    CS_MALLOC(c->setup_data, 1, cs_sles_cudss_setup_t);
    sd = c->setup_data;

    sd->need_update = true;
    sd->data = nullptr;
    sd->matrix = nullptr;
    sd->rhs = nullptr;
    sd->x = nullptr;
    sd->_row_index = nullptr;
    sd->_col_id = nullptr;

    retval = cudssDataCreate(_handle, &(sd->data));
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssDataCreate: %s\n"),
                __func__, _status_str[retval]);

    const cs_matrix_type_t cs_mat_type = cs_matrix_get_type(a);
    const int db_size = cs_matrix_get_diag_block_size(a);
    const cs_halo_t *halo = cs_matrix_get_halo(a);

    /* Periodicity is not handled (at least not in serial mode), as the matrix
       is not square due to ghost cells */

    if (halo != nullptr) {
      bool have_perio = false;
      if (halo->n_transforms > 0)
        have_perio = true;
      assert(have_perio == false);
    }

    /* TODO: handle periodicity, by renumbering local periodic cells
       so as to use the main (and not ghost) cell id */

    if (   db_size > 1
        || (cs_mat_type != CS_MATRIX_CSR && cs_mat_type != CS_MATRIX_MSR)) {
      bft_error
        (__FILE__, __LINE__, 0,
         _("Matrix type %s with block size %d for system \"%s\" "
           "is not usable by cuDSS.\n"
           "Only block size 1 with CSR or MSR type "
           "is currently supported by cuDSS."),
         cs_matrix_get_type_name(a), db_size,
         name);
    }

    _setup_matrix(c, a);

  }

  /* Case where setup data is already present (simple update) */

  else {
    const cs_real_t *a_val = nullptr;
    cs_matrix_get_csr_arrays(a, nullptr, nullptr, &a_val);

    retval = cudssMatrixSetValues(sd->matrix,
                                  const_cast<cs_real_t *>(a_val));
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssMatrixSetValues: %s\n"),
                __func__, _status_str[retval]);
  }

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call cuDSS linear equation solver.
 *
 * \warning The precision, r_norm, and n_iter parameters are ignored here.
 *          the matching configuration options should be set earlier, using
 *          the \ref cs_sles_cudss_set_config function
 *
 * \param[in, out]  context        pointer to cuDSS solver info and context
 *                                 (actual type: cs_sles_cudss_t  *)
 * \param[in]       name           pointer to system name
 * \param[in]       a              matrix
 * \param[in]       verbosity      associated verbosity
 * \param[in]       precision      solver precision
 * \param[in]       r_norm         residual normalization
 * \param[out]      n_iter         number of "equivalent" iterations
 * \param[out]      residual       residual
 * \param[in]       rhs            right hand side
 * \param[in]       vx_ini         initial system solution
 *                                 (vx if nonzero, nullptr if zero)
 * \param[in, out]  vx             system solution
 * \param[in]       aux_size       number of elements in aux_vectors (in bytes)
 * \param           aux_vectors    optional working area
 *                                 (internal allocation if nullptr)
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_cudss_solve(void               *context,
                   const char          *name,
                   const cs_matrix_t   *a,
                   int                  verbosity,
                   double               precision,
                   double               r_norm,
                   int                 *n_iter,
                   double              *residual,
                   const cs_real_t     *rhs,
                   cs_real_t           *vx_ini,
                   cs_real_t           *vx,
                   size_t               aux_size,
                   void                *aux_vectors)
{
  CS_UNUSED(aux_size);
  CS_UNUSED(aux_vectors);

  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_cudss_t  *c = (cs_sles_cudss_t *)context;
  cs_sles_cudss_setup_t  *sd = c->setup_data;

  cudssSetStream(_handle, cs_cuda_get_stream(0));

  bool factorize = false;

  if (sd == nullptr) {
    cs_sles_cudss_setup(c, name, a, verbosity);
    sd = c->setup_data;
    factorize = true;
  }

  /* Vectors */

  const int n_rows = cs_matrix_get_n_rows(a);
  const int n_cols = cs_matrix_get_n_columns(a);

  constexpr cudaDataType_t
    value_type = (sizeof(cs_real_t) == 8) ? CUDA_R_64F : CUDA_R_32F;

  cudssStatus_t retval
    = cudssMatrixCreateDn(&(sd->rhs),
                          n_rows,
                          1,      // ncols
                          n_rows, // ld,
                          const_cast<cs_real_t *>(rhs),
                          value_type,
                          CUDSS_LAYOUT_COL_MAJOR);

  if (retval == CUDSS_STATUS_SUCCESS)
    retval
      = cudssMatrixCreateDn(&(sd->x),
                            n_rows,
                            1,      // ncols
                            n_rows, // ld,
                            const_cast<cs_real_t *>(vx),
                            value_type,
                            CUDSS_LAYOUT_COL_MAJOR);

  if (retval != CUDSS_STATUS_SUCCESS)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: error calling cudssMatrixCreateDn: %s\n"),
              __func__, _status_str[retval]);

  int       its = 1;
  double    _residual = -1;

  /* Resolution */

  cs_fp_exception_disable_trap();

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  if (factorize) {

    // Reordering & symbolic factorization
    retval = cudssExecute(_handle,
                          CUDSS_PHASE_ANALYSIS,
                          c->config,
                          sd->data,
                          sd->matrix,
                          sd->x,
                          sd->rhs);
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssExecute: %s for "
                  "analysis phase\n"),
                __func__, _status_str[retval]);

    c->n_analysis += 1;
    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(c->t_sym_fact), &t1, &t0);
    t1 = t0;

  }
  if (sd->need_update) {

    // Numerical factorization
    retval = cudssExecute(_handle,
                          CUDSS_PHASE_FACTORIZATION,
                          c->config,
                          sd->data,
                          sd->matrix,
                          sd->x,
                          sd->rhs);
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssExecute: %s for "
                  "numerical factorization phase\n"),
                __func__, _status_str[retval]);

    sd->need_update = false;

    c->n_num_fact += 1;
    t0 = cs_timer_time();
    cs_timer_counter_add_diff(&(c->t_num_fact), &t1, &t0);
    t1 = t0;

  }

  // Solve the system
  retval = cudssExecute(_handle,
                        CUDSS_PHASE_SOLVE,
                        c->config,
                        sd->data,
                        sd->matrix,
                        sd->x,
                        sd->rhs);
  if (retval != CUDSS_STATUS_SUCCESS)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: error calling cudssExecute: %s for "
                "solve phase\n"),
              __func__, _status_str[retval]);


  int64_t mem_estimates[16];
  size_t  size_written = 0;

  retval = cudssDataGet(_handle,
                        sd->data,
                        CUDSS_DATA_MEMORY_ESTIMATES,
                        mem_estimates,
                        sizeof(int64_t)*16,
                        &size_written);
  if (retval == CUDSS_STATUS_SUCCESS) {
    int n = cs::min(16, (int)(size_written/sizeof(int64_t)));
    for (int i = 0; i < n; i++)
      c->mem_estimates[i] = cs::max(c->mem_estimates[i], mem_estimates[i]);
  }

  cs_fp_exception_restore_trap();

  cvg = CS_SLES_CONVERGED;

  *residual = _residual;
  *n_iter = its;

  /* Update return values */

  c->n_solves += 1;
  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free cuDSS linear equation solver setup context.
 *
 * This function frees resolution-related data, such as
 * buffers and preconditioning but does not free the whole context,
 * as info used for logging (especially performance data) is maintained.
 *
 * \param[in, out]  context  pointer to cuDSS solver info and context
 *                           (actual type: cs_sles_cudss_t  *)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_free(void  *context)
{
  cs_timer_t t0;
  t0 = cs_timer_time();

  cs_sles_cudss_t  *c = (cs_sles_cudss_t *)context;
  cs_sles_cudss_setup_t *sd = c->setup_data;

  if (sd != nullptr) {

    if (c->keep_data == true) {
      sd->need_update = true;
      return;
    }

    cudssStatus_t retval
      = cudssDataDestroy(_handle, sd->data);
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssDataDestroy: %s\n"),
                __func__, _status_str[retval]);
    else
      sd->data = nullptr;

    retval  = cudssMatrixDestroy(sd->x); sd->x = nullptr;
    if (retval == CUDSS_STATUS_SUCCESS) {
      retval = cudssMatrixDestroy(sd->rhs); sd->rhs = nullptr;
    }
    if (retval == CUDSS_STATUS_SUCCESS) {
      retval = cudssMatrixDestroy(sd->matrix); sd->matrix = nullptr;
    }
    if (retval != CUDSS_STATUS_SUCCESS)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: error calling cudssMatrixDestroy: %s\n"),
                __func__, _status_str[retval]);
    else
      sd->matrix = nullptr;

    CS_FREE(sd->_row_index);
    CS_FREE(sd->_col_id);

  }
  if (c->setup_data != nullptr)
    CS_FREE(c->setup_data);

  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(c->t_solve), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log sparse linear equation solver info.
 *
 * \param[in]  context   pointer to cuDSS solver info and context
 *                       (actual type: cs_sles_cudss_t  *)
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_log(const void  *context,
                 cs_log_t     log_type)
{
  const cs_sles_cudss_t  *c = (const cs_sles_cudss_t *)context;

  const char m_type[] = "CSR";

  if (log_type == CS_LOG_SETUP) {

    cs_log_printf(log_type,
                  _("  Solver type:                       cuDSS\n"
                    "    Matrix format:                     %s\n"),
                  m_type);

  }
  else if (log_type == CS_LOG_PERFORMANCE) {

    cs_log_printf(log_type,
                  _("\n"
                    "  Solver type:                   cuDSS\n"
                    "  Number of analyses:            %12d\n"
                    "  Number of factorizations:      %12d\n"
                    "  Number of calls:               %12d\n"
                    "  Permanent device memory:       %12.1f Mb\n"
                    "  Peak device memory:            %12.1f Mb\n"
                    "  Permanent host memory:         %12.1f Mb\n"
                    "  Peak host memory:              %12.1f Mb\n"
                    "  Total analysis time:           %12.3f\n"
                    "  Total factorization time:      %12.3f\n"
                    "  Total solution time:           %12.3f\n"),
                  c->n_analysis, c->n_num_fact, c->n_solves,
                  (double)c->mem_estimates[0]/1e6,
                  (double)c->mem_estimates[1]/1e6,
                  (double)c->mem_estimates[2]/1e6,
                  (double)c->mem_estimates[3]/1e6,
                  c->t_sym_fact.nsec*1e-9,
                  c->t_num_fact.nsec*1e-9,
                  c->t_solve.nsec*1e-9);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on cuDSS library.
 *
 * \param[in]  log_type  log type
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_cudss_library_info(cs_log_t  log_type)
{
  int major = 0, minor = 0, patch = 0;

  cudssGetProperty(MAJOR_VERSION, &major);
  cudssGetProperty(MINOR_VERSION, &major);
  cudssGetProperty(PATCH_LEVEL, &patch);

  cs_log_printf(log_type,
                "    cuDSS v0.%d.%d-%d\n", major, minor, patch);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
