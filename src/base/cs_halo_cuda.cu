/*============================================================================
 * Functions dealing with ghost cells using CUDA.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
 * CUDA kernel to gather sparse data to a dense array.
 *
 * parameters:
 *   length  <-- length of sub-list to gather
 *   src_ids <-- ids of elements to gather
 *   src     <-- source values
 *   var     <-> destination values
 *----------------------------------------------------------------------------*/

__global__ static void
_gather_from_var(cs_lnum_t         length,
                 const cs_lnum_t  *src_ids,
                 const cs_real_t  *src,
                 cs_real_t        *dest)
{
  cs_lnum_t i = blockIdx.x*blockDim.x + threadIdx.x;

  if (i < length)
    dest[i] = src[src_ids[i]];
}

/*----------------------------------------------------------------------------
 * CUDA kernel to gather strided sparse data to a dense array.
 *
 * parameters:
 *   stride    <-- data stride
 *   length    <-- length of sub-list to gather
 *   src_ids   <-- src_ids of elements to send
 *   src       <-- source values
 *   var       <-> destination values
 *----------------------------------------------------------------------------*/

__global__ static void
_gather_from_var_strided(cs_lnum_t   stride,
                         cs_lnum_t   length,
                         cs_lnum_t  *src_ids,
                         cs_real_t  *src,
                         cs_real_t  *dest)
{
  cs_lnum_t i = blockIdx.x*blockDim.x + threadIdx.x;

  if (i < length) {
    for (cs_lnum_t j = 0; j < stride; j++)
      dest[i*stride + j] = src[src_ids[i]*stride + j];
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
 * A local state handler may be provided, or the default state handler will
 * be used.
 *
 * A local state and/or buffer may be provided, or the default (global) state
 * and buffer will be used. If provided explicitely,
 * the buffer must be of sufficient size.
 *
 * \param[in]       halo        pointer to halo structure
 * \param[in]       sync_mode   synchronization mode (standard or extended)
 * \param[in]       data_type   data type
 * \param[in]       stride      number of (interlaced) values by entity
 * \param[in]       var         pointer to value array (device)
 * \param[out]      send_buf    pointer to send buffer, NULL for global buffer
 * \param[in, out]  hs          pointer to halo state, NULL for global state
 */
/*----------------------------------------------------------------------------*/

void
cs_halo_sync_pack_cuda_real(const cs_halo_t  *halo,
                            cs_halo_type_t    sync_mode,
                            cs_datatype_t     data_type,
                            int               stride,
                            cs_real_t        *var,
                            void             *send_buf,
                            cs_halo_state_t  *hs)
{
  if (halo == NULL)
    return;

  void *_send_buffer_p = cs_halo_sync_pack_init_state(halo,
                                                      sync_mode,
                                                      CS_REAL_TYPE,
                                                      stride,
                                                      var,
                                                      send_buf,
                                                      hs);

  cs_real_t *_send_buffer = (cs_real_t *)cs_get_device_ptr(_send_buffer_p);

  cs_lnum_t end_shift = 0;

  if (sync_mode == CS_HALO_STANDARD)
    end_shift = 1;

  else if (sync_mode == CS_HALO_EXTENDED)
    end_shift = 2;

  cs_lnum_t *send_list = (cs_lnum_t *)cs_get_device_ptr(halo->send_list);

  /* Assemble buffers for halo exchange. */

  cudaStream_t nstream[halo->n_c_domains];
  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++)
    cudaStreamCreate(&nstream[ii]);

  for (cs_lnum_t rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {
    cs_lnum_t start = halo->send_index[2*rank_id];
    cs_lnum_t length =   halo->send_index[2*rank_id + end_shift]
                       - halo->send_index[2*rank_id];
    if (length > 0) {
      unsigned int blocksize = (length < 256)? length:256;
      unsigned int gridsize = (unsigned int)ceil((double)length/blocksize);
      if (stride == 1)
        _gather_from_var<<<gridsize, blocksize, 0, nstream[rank_id]>>>
          (length, send_list+start, var, _send_buffer+start);
      else
        _gather_from_var_strided<<<gridsize, blocksize, 0, nstream[rank_id]>>>
          (stride, length, send_list+start, var, _send_buffer+stride*start);
    }
  }

  for (cs_lnum_t ii = 0; ii < halo->n_c_domains; ii++){
    CS_CUDA_CHECK(cudaStreamSynchronize(nstream[ii]));
    cudaStreamDestroy(nstream[ii]);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
