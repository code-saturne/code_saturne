#ifndef __CS_HALO_CUDA_H__
#define __CS_HALO_CUDA_H__

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

#include "cs_halo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
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
                                   cs_real_t         send_buffer[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HALO_CUDA_H__ */

