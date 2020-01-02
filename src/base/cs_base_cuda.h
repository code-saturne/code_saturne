#ifndef __CS_BASE_CUDA_H__
#define __CS_BASE_CUDA_H__

/*============================================================================
 * Definitions, global variables, and base functions for CUDA
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_log.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

#if defined(HAVE_CUDA)

extern int  cs_glob_cuda_device_id;

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_CUDA)

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on available CUDA devices.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_device_info(cs_log_t  log_id);

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_CUDA_H__ */
