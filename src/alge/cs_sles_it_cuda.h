#ifndef __CS_SLES_IT_CUDA_H__
#define __CS_SLES_IT_CUDA_H__

/*============================================================================
 * Sparse Linear Equation Solvers using CUDA
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "alge/cs_matrix.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_pc.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * User function prototypes
 *============================================================================*/

/*============================================================================
 * Semi-private function prototypes
 *============================================================================*/

#if defined(__CUDACC__)

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   stream <-- CUDA stream
 *   x      <-- first vector
 *   y      <-- vector
 *
 * return:
 *   result of s = x.x
 *----------------------------------------------------------------------------*/

double
cs_sles_it_dot_product
(
  const cs_sles_it_t  *c,
  cudaStream_t         stream,
  const cs_real_t     *x,
  const cs_real_t     *y
);

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   stream <-- CUDA stream
 *   x      <-- first vector
 *   y      <-- second vector
 *
 * return:
 *   result of s = x.y
 *----------------------------------------------------------------------------*/

double
cs_sles_it_dot_product_xx
(
  const cs_sles_it_t  *c,
  cudaStream_t         stream,
  const cs_real_t     *x,
  const cs_real_t     *y
);

/*----------------------------------------------------------------------------
 * Compute dot product, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   stream <-- CUDA stream
 *   x      <-- vector
 *
 * return:
 *   result of s = x.x
 *----------------------------------------------------------------------------*/

double
cs_sles_it_dot_product_xx
(
  const cs_sles_it_t  *c,
  cudaStream_t         stream,
  const cs_real_t     *x
);

/*----------------------------------------------------------------------------
 * Compute 2 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   stream <-- CUDA stream
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   xx     --> result of s1 = x.x
 *   xy     --> result of s2 = x.y
 *----------------------------------------------------------------------------*/

void
cs_sles_it_dot_products_xx_xy
(
  const cs_sles_it_t  *c,
  cudaStream_t         stream,
  const cs_real_t     *x,
  const cs_real_t     *y,
  double              *xx,
  double              *xy
);

/*----------------------------------------------------------------------------
 * Compute 3 dot products, summing result over all ranks.
 *
 * parameters:
 *   c      <-- pointer to solver context info
 *   stream <-- CUDA stream
 *   x      <-- first vector
 *   y      <-- second vector
 *   z      <-- third vector
 *   xx     --> result of s1 = x.x
 *   xy     --> result of s2 = x.y
 *   yz     --> result of s3 = y.z
 *----------------------------------------------------------------------------*/

void
cs_sles_it_dot_products_xx_xy_yz
(
  const cs_sles_it_t  *c,
  cudaStream_t         stream,
  const cs_real_t     *x,
  const cs_real_t     *y,
  const cs_real_t     *z,
  double              *xx,
  double              *xy,
  double              *yz
);

#endif // defined(__CUDACC__)

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_cuda_jacobi(cs_sles_it_t              *c,
                       const cs_matrix_t         *a,
                       cs_lnum_t                  diag_block_size,
                       cs_sles_it_convergence_t  *convergence,
                       const cs_real_t           *rhs,
                       cs_real_t                 *vx_ini,
                       cs_real_t                 *vx,
                       size_t                     aux_size,
                       void                      *aux_vectors);

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using block Jacobi.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- linear equation matrix
 *   diag_block_size <-- diagonal block size
 *   rotation_mode   <-- halo update option for rotational periodicity
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_cuda_block_jacobi(cs_sles_it_t              *c,
                             const cs_matrix_t         *a,
                             cs_lnum_t                  diag_block_size,
                             cs_sles_it_convergence_t  *convergence,
                             const cs_real_t           *rhs,
                             cs_real_t                 *vx_ini,
                             cs_real_t                 *vx,
                             size_t                     aux_size,
                             void                      *aux_vectors);

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using flexible preconditioned conjugate gradient.
 *
 * Compared to standard PCG, FCG supports variable preconditioners.
 *
 * This variant, described in \cite Notay:2015, allows computing the
 * required inner products with a single global communication.
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_cuda_fcg(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *vx_ini,
                    cs_real_t                 *vx,
                    size_t                     aux_size,
                    void                      *aux_vectors);

/*----------------------------------------------------------------------------
 * Solution of A.vx = Rhs using optimised preconditioned GCR (CUDA version).
 *
 * parameters:
 *   c               <-- pointer to solver context info
 *   a               <-- matrix
 *   diag_block_size <-- diagonal block size (unused here)
 *   convergence     <-- convergence information structure
 *   rhs             <-- right hand side
 *   vx_ini          <-- initial system solution
 *                       (vx if nonzero, nullptr if zero)
 *   vx              <-> system solution
 *   aux_size        <-- number of elements in aux_vectors (in bytes)
 *   aux_vectors     --- optional working area (allocation otherwise)
 *
 * returns:
 *   convergence state
 *----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_it_cuda_gcr(cs_sles_it_t              *c,
                    const cs_matrix_t         *a,
                    cs_lnum_t                  diag_block_size,
                    cs_sles_it_convergence_t  *convergence,
                    const cs_real_t           *rhs,
                    cs_real_t                 *vx_ini,
                    cs_real_t                 *vx,
                    size_t                     aux_size,
                    void                      *aux_vectors);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SLES_IT_CUDA_H__ */
