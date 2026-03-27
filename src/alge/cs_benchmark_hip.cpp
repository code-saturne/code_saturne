/*============================================================================
 * Low-level operator benchmarking with HIP.
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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <algorithm>
#include <assert.h>
#include <list>

#include <stdio.h>

#if defined(HAVE_HIP)
#if defined(HAVE_HIPBLAS)
#include <hipblas.h>
#endif

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#if defined(HAVE_HIPSPARSE)
#include <hipsparse.h>
#endif
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_base_accel.h"
#include "base/cs_base_hip.h"

#include "base/cs_log.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_benchmark_hip.h"

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
    atomicAdd(&y[ii], xa[face_id] * x[jj]);
    atomicAdd(&y[jj], xa[face_id] * x[ii]);
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
cs_mat_vec_exdiag_native_sym_hip(cs_lnum_t           n_faces,
                                 const cs_lnum_2_t  *face_cell,
                                 const cs_real_t    *xa,
                                 const cs_real_t    *x,
                                 cs_real_t          *y)
{
  unsigned int blocksize = 512;
  unsigned int gridsize  = (unsigned int)ceil((double)n_faces / blocksize);
  cudaStream_t stream = cs_hip_get_stream(0);

  _mat_vec_exdiag_native_sym<<<gridsize, blocksize, 0, stream>>>
    (n_faces, (const cs_lnum_t *)face_cell, xa, x, y);

  CS_HIP_CHECK(hipStreamSynchronize(stream));
  CS_HIP_CHECK(hipGetLastError());
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
