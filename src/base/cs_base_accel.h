#ifndef __CS_BASE_ACCEL_H__
#define __CS_BASE_ACCEL_H__

/*============================================================================
 * Definitions, global variables, and base functions for accelerators.
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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_mem.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_FREE_HD(_ptr) \
cs_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

#if defined(HAVE_ACCEL)

extern int cs_mpi_device_support;

#else

#define cs_mpi_device_support 0;

#endif

/*! Default queue for SYCL */

#if defined(SYCL_LANGUAGE_VERSION) && !defined(CS_GLOB_SYCL_QUEUE_IS_DEFINED)
extern sycl::queue  cs_glob_sycl_queue;
#define CS_GLOB_SYCL_QUEUE_IS_DEFINED 1
#endif

/*=============================================================================
 * Public C function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return currently associated device id.
 *
 * \returns currently available device id, or -1 if none is available.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

int
cs_get_device_id(void);

#else

static inline int
cs_get_device_id(void)
{
  return -1;
}

#endif

#if defined(HAVE_OPENMP_TARGET)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set OpenMP Offload device based on MPI rank and number of devices.
 *
 * \param[in]  comm            associated MPI communicator
 * \param[in]  ranks_per_node  number of ranks per node (min and max)
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_omp_target_select_default_device(void);

#endif /* defined(HAVE_OPENMP_TARGET) */

#if defined(HAVE_SYCL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set SYCL device based on SYCL device selector.
 *
 * Should be based on based on MPI rank and number of devices in the future.
 *
 * \param[in]  comm            associated MPI communicator
 * \param[in]  ranks_per_node  number of ranks per node (min and max)
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_sycl_select_default_device(void);

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_ACCEL_H__ */
