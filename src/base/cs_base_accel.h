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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * If memory is not allocated on device yet at the call site, it will
 * be allocated automatically by this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void *
cs_get_device_ptr(void  *ptr);

#else

inline static void *
cs_get_device_ptr(void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const(T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr(reinterpret_cast<void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given constant pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * If memory is not allocated on device yet at the call site, it will
 * be allocated automatically by this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

const void *
cs_get_device_ptr_const(const void  *ptr);

#else

inline static const void *
cs_get_device_ptr_const(const void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const(const T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr_const(reinterpret_cast<const void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given constant pointer,
 *        prefetching if applicable.
 *
 * If separate pointers are used on the host and device, the host pointer
 * should be used with this function. In this case, it is assumed that
 * the host and device values have already been synchronized, unless
 * memory is not allocated on device yet at the call site, in which case
 * it will be allocated automatically by this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

const void *
cs_get_device_ptr_const_pf(const void  *ptr);

#else

inline static const void *
cs_get_device_ptr_const_pf(const void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const_pf(const T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr_const_pf(reinterpret_cast<const void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a pointer is associated with a device.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \returns allocation mode associated with pointer
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

cs_alloc_mode_t
cs_check_device_ptr(const void  *ptr);

#else

inline static cs_alloc_mode_t
cs_check_device_ptr(const void  *ptr)
{
  CS_UNUSED(ptr);
  return CS_ALLOC_HOST;
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate device memory with a given host memory pointer.
 *
 * If the host memory is already associated with the device, the existing
 * device pointer is returned. Otherwise, a new device allocation is
 * called and returned.
 *
 * \param [in]  host_ptr  host pointer
 * \param [in]  ni        number of elements
 * \param [in]  size      element size
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void *
cs_associate_device_ptr(void    *host_ptr,
                        size_t   ni,
                        size_t   size);

#else

#define cs_associate_device_ptr(_host_ptr, _ni, _size);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Detach device memory from a given host memory pointer.
 *
 * If the host memory is shared with the device (i.e. using CS_ALLOC_SHARED),
 * device memory stays shared.
 *
 * \param [in]  host_ptr  host pointer
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_disassociate_device_ptr(void  *host_ptr);

#else

#define cs_disassociate_device_ptr(_host_ptr);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set allocation mode for an already allocated pointer.
 *
 * If the allocation mode is different from the previous one,
 * the associated memory will be reallocated with the desired mode,
 * and the previous allocation freed.
 *
 * \param [in, out]  host_ptr   pointer to host pointer to modify
 * \param [in]       mode       desired allocation mode
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_set_alloc_mode(void             **host_ptr,
                  cs_alloc_mode_t     mode);

#else

#define cs_set_alloc_mode(_host_ptr, mode);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_mem_advise_set_read_mostly(void  *ptr);

#else

#define cs_mem_advise_set_read_mostly(ptr);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_mem_advise_unset_read_mostly(void  *ptr);

#else

#define cs_mem_advise_unset_read_mostly(ptr);

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

/*=============================================================================
 * Public C++ function prototypes and definitions.
 *============================================================================*/

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set allocation mode for an already allocated pointer using
 *        pass by reference semantics for the pointer (C++ only).
 *
 * In C++, this function is preferred to direct use of cs_set_alloc_mode,
 * as it allows a more concise syntax and does not require additional casting.
 *
 * \param [in, out]  host_ptr   reference to host pointer to modify
 * \param [in]       mode       desired allocation mode
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

template<typename T>
static inline void
cs_set_alloc_mode_r(T*                &host_ptr,
                    cs_alloc_mode_t    mode)
{
  void *p = host_ptr;
  cs_set_alloc_mode(&p, mode);
  host_ptr = (T *)p;
}

#else

#define cs_set_alloc_mode_r(_host_ptr, mode);

#endif

/*----------------------------------------------------------------------------*/

#endif /* defined(__cplusplus) */

#endif /* __CS_BASE_ACCEL_H__ */
