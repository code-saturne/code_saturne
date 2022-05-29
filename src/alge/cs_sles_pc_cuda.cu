/*============================================================================
 * Sparse Linear Equation Solver Preconditioners using CUDA
 * Sparse Linear Equation Solver Preconditioner driver
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

#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_base_cuda.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"
#include "cs_matrix_spmv_cuda.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_pc.h"
#include "cs_sles_pc_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_sles_pc.c

  CUDA implementation of some linear system solver preconditionners.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* SIMD unit size to ensure SIMD alignement (2 to 4 required on most
 * current architectures, so 16 should be enough on most architectures) */

#define CS_SIMD_SIZE(s) (((s-1)/16+1)*16)

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Apply Jacobi preconditionner.
 *
 * parameters:
 *   n      <-- number of rows
 *   ad_inv <-- diagonal inverse
 *   x_in   <-- input vector
 *   x_out  --> output vector
 *----------------------------------------------------------------------------*/

__global__ static void
_pc_jacobi(cs_lnum_t        n,
           const cs_real_t  ad_inv[restrict],
           const cs_real_t  x_in[restrict],
           cs_real_t        x_out[restrict])
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  while (ii < n) {
    x_out[ii] = ad_inv[ii] * x_in[ii];
    ii += grid_size;
  }
}

/*----------------------------------------------------------------------------
 * Apply Jacobi preconditionner, using same input and output vector.
 *
 * parameters:
 *   n      <-- number of rows
 *   ad_inv <-- diagonal inverse
 *   x      <-> input/output vector
 *----------------------------------------------------------------------------*/

__global__ static void
_pc_jacobi_self(cs_lnum_t        n,
                const cs_real_t  ad_inv[restrict],
                cs_real_t        x[restrict])
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  size_t grid_size = blockDim.x*gridDim.x;

  while (ii < n) {
    x[ii] *= ad_inv[ii];
    ii += grid_size;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function for application of a null-preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_none(void                *context,
                           const cs_real_t     *x_in,
                           cs_real_t           *x_out)
{
  if (x_in != NULL && x_out != x_in) {

    cs_sles_pc_poly_t  *c = (cs_sles_pc_poly_t *)context;
    const cs_lnum_t n_rows = c->n_rows;

    cudaMemcpy(x_out, x_in, n_rows * sizeof(cs_real_t),
               cudaMemcpyDeviceToDevice);
  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------
 * Function for application of a Jacobi preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_jacobi(void                *context,
                             const cs_real_t     *x_in,
                             cs_real_t           *x_out)
{
  cs_sles_pc_poly_t  *c = (cs_sles_pc_poly_t *)context;

  const cs_lnum_t n_rows = c->n_rows;
  const cs_real_t *restrict ad_inv = c->ad_inv;

  cudaStream_t stream = cs_matrix_spmv_cuda_get_stream();

  const unsigned int blocksize = 256;
  unsigned int gridsize = cs_cuda_grid_size(n_rows, blocksize);

  if (x_in != NULL) {
    _pc_jacobi<<<gridsize, blocksize, 0, stream>>>
      (n_rows, ad_inv, x_in, x_out);
  }
  else {
    _pc_jacobi_self<<<gridsize, blocksize, 0, stream>>>
      (n_rows, ad_inv, x_out);
  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------
 * Function for application of a polynomial preconditioner.
 *
 * In cases where it is desired that the preconditioner modify a vector
 * "in place", x_in should be set to NULL, and x_out contain the vector to
 * be modified (\f$x_{out} \leftarrow M^{-1}x_{out})\f$).
 *
 * parameters:
 *   context       <-> pointer to preconditioner context
 *   x_in          <-- input vector
 *   x_out         <-> input/output vector
 *
 * returns:
 *   preconditioner application status
 *----------------------------------------------------------------------------*/

cs_sles_pc_state_t
cs_sles_pc_cuda_apply_poly(void                *context,
                           const cs_real_t     *x_in,
                           cs_real_t           *x_out)
{
  cs_sles_pc_poly_t  *c = (cs_sles_pc_poly_t *)context;

  const cs_lnum_t n_rows = c->n_rows;
  const cs_lnum_t n_aux = (x_in == NULL) ?
    CS_SIMD_SIZE(c->n_cols) + c->n_cols : c->n_cols;

  if (c->n_aux < n_aux) {
    c->n_aux = n_aux;
    cs_alloc_mode_t amode = cs_check_device_ptr(c->ad_inv);
    CS_FREE_HD(c->aux);
    CS_MALLOC_HD(c->aux, c->n_aux, cs_real_t, amode);
  }

  cs_real_t *restrict w = c->aux;
  const cs_real_t *restrict r = x_in;
  const cs_real_t *restrict ad_inv = c->ad_inv;

  if (x_in == NULL) {

    cs_real_t *restrict _r = c->aux + CS_SIMD_SIZE(c->n_cols);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      _r[ii] = x_out[ii];

    r = _r;

  }

# pragma omp parallel for if(n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_rows; ii++)
    x_out[ii] = r[ii] * ad_inv[ii];

  for (int deg_id = 1; deg_id <= c->poly_degree; deg_id++) {

    /* Compute Wk = (A-diag).Gk */

    cs_matrix_vector_multiply_partial(c->a, CS_MATRIX_SPMV_E, x_out, w);

#   pragma omp parallel for if(n_rows > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      x_out[ii] = (r[ii] - w[ii]) * ad_inv[ii];

  }

  return CS_SLES_PC_CONVERGED;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
