/*============================================================================
 * Low-level operator benchmarking with CUDA.
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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>
#include <list>

#include <stdio.h>

#if defined(HAVE_CUDA)
#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cusparse.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"

#include "cs_base_cuda.h"
// #include "cs_cuda_contrib.h"

#include "cs_log.h"
#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_numbering.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_benchmark_cuda.h"

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compatibility macro for __ldg (load from generic memory) intrinsic,
 * forcing load from read-only texture cache.
 *
 * This was not available in (very old) CUDA architectures.
 *----------------------------------------------------------------------------*/

#if __CUDA_ARCH__ < 350
#define __ldg(ptr) *(ptr);
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * SpMV extradiagonal terms using native to face-based array and scatter
 * approach, handling conflicts through atomic add.
 *
 * Symmetric matrix case.
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

__global__ static void
_mat_vec_exdiag_native_sym(cs_lnum_t        n_faces,
                           const cs_lnum_t  * __restrict__ face_cell,
                           const cs_real_t  *__restrict__ xa,
                           const cs_real_t  *__restrict__ x,
                           cs_real_t        *__restrict__ y)
{
  cs_lnum_t ii, jj;
  cs_lnum_t face_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (face_id < n_faces) {
    ii = face_cell[2 * face_id];
    jj = face_cell[2 * face_id + 1];
    atomicAdd(&y[jj], xa[face_id] * x[ii]);
    atomicAdd(&y[ii], xa[face_id] * x[jj]);
  }
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * SpMV extradiagonal terms using native to face-based array and scatter
 * approach, handling conflicts through atomic add.
 *
 * Symmetric matrix case.
 *
 * parameters:
 *   n_faces         <-- local number of internal faces
 *   face_cell       <-- face -> cells connectivity
 *   xa              <-- extradiagonal values
 *   x               <-- vector
 *   y               <-> vector
 *----------------------------------------------------------------------------*/

void
cs_mat_vec_exdiag_native_sym_cuda(cs_lnum_t           n_faces,
                                  const cs_lnum_2_t  *face_cell,
                                  const cs_real_t    *xa,
                                  const cs_real_t    *x,
                                  cs_real_t          *y)
{
  unsigned int blocksize = 512;
  unsigned int gridsize  = (unsigned int)ceil((double)n_faces / blocksize);

  _mat_vec_exdiag_native_sym<<<gridsize, blocksize>>>
    (n_faces, (const cs_lnum_t *)face_cell, xa, x, y);

  cudaDeviceSynchronize();
  CS_CUDA_CHECK(cudaGetLastError());
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
