/*============================================================================
 * Functions dealing with ghost cells using CUDA.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Gather values from a source buffer to another buffer, by blocks
 *
 * parameters:
 *   n           <-- number of elements to gather
 *   send_blocks <-- ids of values in src array
 *   ids         <-- ids of values in src array
 *   src         <-- array of values to gather from (source)
 *   dest        <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_gather_block_cu(const cs_lnum_t   send_blocks[],
                 const cs_lnum_t   ids[],
                 const T           src[],
                 T                 dest[])
{
  cs_lnum_t i = send_blocks[blockIdx.x*2] + threadIdx.x;
  if (i < send_blocks[blockIdx.x*2 + 1])
    dest[i] = src[ids[i]];
}

/*----------------------------------------------------------------------------
 * Gather strided values from a source buffer to another buffer, by blocks
 *
 * parameters:
 *   n           <-- number of elements to gather
 *   stride      <-- number of elements to gather
 *   send_blocks <-- ids of values in src array
 *   ids         <-- ids of values in src array
 *   src         <-- array of values to gather from (source)
 *   dest        <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_gather_block_strided_cu(const cs_lnum_t   send_blocks[],
                         cs_lnum_t         stride,
                         const cs_lnum_t   ids[],
                         const T           src[],
                         T                 dest[])
{
  cs_lnum_t i = send_blocks[blockIdx.x*2] + threadIdx.x;
  if (i < send_blocks[blockIdx.x*2 + 1]) {
    cs_lnum_t j = ids[i];
    for (cs_lnum_t k = 0; k < stride; k++)
      dest[i*stride + k] = src[j*stride + k];
  }
}

/*----------------------------------------------------------------------------
 * Gather values from a source buffer to another buffer, by blocks
 *
 * parameters:
 *   n           <-- number of elements to gather
 *   ids         <-- ids of values in src array
 *   src         <-- array of values to gather from (source)
 *   dest        <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_gather_full_cu(cs_lnum_t         n,
                const cs_lnum_t   ids[],
                const T           src[],
                T                 dest[])
{
  cs_lnum_t ii = blockIdx.x*blockDim.x + threadIdx.x;
  if (ii < n)
    dest[ii] = src[ids[ii]];
}

/*----------------------------------------------------------------------------
 * Gather strided values from a source buffer to another buffer, by blocks
 *
 * parameters:
 *   n           <-- number of elements to gather
 *   stride      <-- number of elements to gather
 *   ids         <-- ids of values in src array
 *   src         <-- array of values to gather from (source)
 *   dest        <-- array of gathered values (destination)
 *----------------------------------------------------------------------------*/

template <typename T>
__global__ static void
_gather_full_strided_cu(cs_lnum_t         n,
                        cs_lnum_t         stride,
                        const cs_lnum_t   ids[],
                        const T           src[],
                        T                 dest[])
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
 * \param[in]   halo          pointer to halo structure
 * \param[in]   sync_mode     synchronization mode (standard or extended)
 * \param[in]   data_type     data type
 * \param[in]   stride        number of (interlaced) values by entity
 * \param[in]   var           pointer to value array (device)
 * \param[out]  send_buffer   pointer to send buffer, NULL for global buffer
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_cuda_pack_send_buffer(const cs_halo_t   *halo,
                              cs_halo_type_t     sync_mode,
                              cs_datatype_t      data_type,
                              cs_lnum_t          stride,
                              const void        *val,
                              void              *send_buffer)
{
  const size_t block_size = halo->std_send_block_size;
  const cs_lnum_t *send_list = halo->send_list;
  const cs_lnum_t *send_blocks = nullptr;

  size_t n_send = 0, n_blocks = 0;

  /* If we do not need to send the full halo, use blocks
     to possibly allow for threading */

  if (sync_mode == CS_HALO_STANDARD && halo->n_std_send_blocks > 0) {
    n_blocks = halo->n_std_send_blocks;
    send_blocks = halo->std_send_blocks;
  }
  else {
    n_send = halo->n_send_elts[1];
    n_blocks = (n_send % block_size) ? n_send/block_size + 1 : n_send/block_size;
  }

  cudaStream_t stream = cs_cuda_get_stream(0);

  if (data_type == CS_REAL_TYPE) {

    cs_real_t *buffer = (cs_real_t *)send_buffer;
    const cs_real_t *var = (const cs_real_t *)val;

    if (stride == 1) {
      if (send_blocks == nullptr)
        _gather_full_cu<<<n_blocks, block_size, 0, stream>>>
          (n_send, send_list, var, buffer);
      else
        _gather_block_cu<<<n_blocks, block_size, 0, stream>>>
          (send_blocks, send_list, var, buffer);
    }
    else {
      if (send_blocks == nullptr)
        _gather_full_strided_cu<<<n_blocks, block_size, 0, stream>>>
          (n_send, stride, send_list, var, buffer);
      else
        _gather_block_strided_cu<<<n_blocks, block_size, 0, stream>>>
          (send_blocks, stride, send_list, var, buffer);
    }

  }
  else {

    unsigned char *buffer = (unsigned char *)send_buffer;
    const unsigned char *var = (const unsigned char *)val;

    cs_lnum_t elt_size = cs_datatype_size[data_type] * stride;

    if (send_blocks == nullptr)
      _gather_full_strided_cu<<<n_blocks, block_size, 0, stream>>>
        (n_send, elt_size, send_list, var, buffer);
    else
      _gather_block_strided_cu<<<n_blocks, block_size, 0, stream>>>
        (send_blocks, elt_size, send_list, var, buffer);

  }

  CS_CUDA_CHECK(cudaStreamSynchronize(stream));
  CS_CUDA_CHECK(cudaGetLastError());
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
