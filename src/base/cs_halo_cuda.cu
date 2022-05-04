/*============================================================================
 * Functions dealing with ghost cells using CUDA.
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
#include <cuda.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_base_accel.h"
#include "cs_base_cuda.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_halo.h"
#include "cs_halo_cuda.h"

/*----------------------------------------------------------------------------*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Gather real values from a source buffer to another buffer.
 *
 * parameters:
 *   n    <-- number of elements to gather
 *   ids  <-- ids of values in src array
 *   src  <-- array of values to gather from (source)
 *   dest <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

__global__ static void
_gather_real(cs_lnum_t        n,
             const cs_lnum_t  ids[],
             const cs_real_t  src[],
             cs_real_t        dest[])
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if (ii < n)
    dest[ii] = src[ids[ii]];
}

/*----------------------------------------------------------------------------
 * Gather strided real values from a source buffer to another buffer.
 *
 * parameters:
 *   n      <-- number of elements to gather
 *   stride <-- stride (dimension) of elements to gather
 *   ids    <-- ids of values in src array
 *   src    <-- array of values to gather from (source)
 *   dest   <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

__global__ static void
_gather_real_strided(cs_lnum_t        n,
                     cs_lnum_t        stride,
                     const cs_lnum_t  ids[],
                     const cs_real_t  src[],
                     cs_real_t        dest[])
{
  cs_lnum_t i = blockIdx.x*blockDim.x + threadIdx.x;

  if (i < n) {
    cs_lnum_t j = ids[i];
    for (cs_lnum_t k = 0; k < stride; k++)
      dest[i*stride + k] = src[j*stride + k];
  }
}

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Pack cs_real_t halo data to send into dense buffer, using CUDA.
 *
 * A local state and/or buffer may be provided, or the default (global) state
 * and buffer will be used. If provided explicitely,
 * the buffer must be of sufficient size.
 *
 * \param[in]   halo         pointer to halo structure
 * \param[in]   sync_mode    synchronization mode (standard or extended)
 * \param[in]   stride       number of (interlaced) values by entity
 * \param[in]   var          pointer to value array (device)
 * \param[out]  send_buffer  pointer to send buffer, NULL for global buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_cuda_pack_send_buffer_real(const cs_halo_t  *halo,
                                   cs_halo_type_t    sync_mode,
                                   cs_lnum_t         stride,
                                   const cs_real_t   var[],
                                   cs_real_t         send_buffer[])
{
  const cs_lnum_t end_shift = (sync_mode == CS_HALO_STANDARD) ? 1 : 2;

  const cs_lnum_t *send_list = (cs_lnum_t *)cs_get_device_ptr(halo->send_list);

  cudaStream_t nstream[halo->n_c_domains];
  for (cs_lnum_t i = 0; i < halo->n_c_domains; i++)
    cudaStreamCreate(&nstream[i]);

  for (cs_lnum_t rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
    cs_lnum_t start = halo->send_index[2*rank_id];
    cs_lnum_t length =   halo->send_index[2*rank_id + end_shift]
                       - halo->send_index[2*rank_id];
    if (length == 0)
      continue;

    unsigned int blocksize = (length < 256)? length:256;
    unsigned int gridsize = (unsigned int)ceil((double)length/blocksize);

    if (stride == 1)
      _gather_real<<<gridsize, blocksize, 0, nstream[rank_id]>>>
        (length, send_list+start, var, send_buffer+start);
    else
      _gather_real_strided<<<gridsize, blocksize, 0, nstream[rank_id]>>>
        (length, stride, send_list+start, var, send_buffer+(start*stride));
  }

  for (cs_lnum_t i = 0; i < halo->n_c_domains; i++) {
    CS_CUDA_CHECK(cudaStreamSynchronize(nstream[i]));
    cudaStreamDestroy(nstream[i]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
