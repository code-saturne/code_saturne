#ifndef __CS_BASE_ACCEL_H__
#define __CS_BASE_ACCEL_H__

/*============================================================================
 * Definitions, global variables, and base functions for accelerators.
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

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls cs_malloc_hd(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer is returned.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 *   _mode <-- allocation mode.
 */

#define CS_MALLOC_HD(_ptr, _ni, _type, _mode) \
_ptr = (_type *) cs_malloc_hd(_mode, _ni, sizeof(_type), \
                                   #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls cs_realloc_hd(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * If the allocation parameters are unchanged, no actual reallocation
 * occurs.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 *   _mode <-- allocation mode.
 */

#define CS_REALLOC_HD(_ptr, _ni, _type, _mode) \
_ptr = (_type *) cs_realloc_hd(_ptr, _mode, _ni, sizeof(_type), \
                               #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls cs_mem_free_any(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be used with this
 * function.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#define CS_FREE_HD(_ptr) \
cs_free_hd(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

/*
 * Free allocated memory.
 *
 * This macro calls cs_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be used with this
 * function.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#define CS_FREE(_ptr) \
cs_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Variable value type.
 *----------------------------------------------------------------------------*/

/*!
 * Allocation modes for accelerated code.
 */

typedef enum {

  CS_ALLOC_HOST,                /*!< allocation on host only */
  CS_ALLOC_HOST_DEVICE,         /*!< allocation on host and device */
  CS_ALLOC_HOST_DEVICE_PINNED,  /*!< allocation on host and device,
                                  using page-locked memory on host
                                  if possible */
  CS_ALLOC_HOST_DEVICE_SHARED,  /*!< allocation on host and device,
                                  using mapped/shared memory */
  CS_ALLOC_DEVICE               /*!< allocation on device only */

} cs_alloc_mode_t;

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

#if defined(HAVE_ACCEL)

extern cs_alloc_mode_t  cs_alloc_mode;

#else

#define cs_alloc_mode CS_ALLOC_HOST

#endif

/*=============================================================================
 * Public function prototypes
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
 * \brief Allocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate allocation function based on
 * the requested mode, and allows introspection of the allocated memory.
 *
 * If separate pointers are used on the host and device,
 * the host pointer is returned.
 *
 * \param [in]  mode       allocation mode
 * \param [in]  ni         number of elements
 * \param [in]  size       element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void *
cs_malloc_hd(cs_alloc_mode_t   mode,
             size_t            ni,
             size_t            size,
             const char       *var_name,
             const char       *file_name,
             int               line_num);

#else

inline static void *
cs_malloc_hd(cs_alloc_mode_t   mode,
             size_t            ni,
             size_t            size,
             const char       *var_name,
             const char       *file_name,
             int               line_num)
{
  CS_UNUSED(mode);
  return bft_mem_malloc(ni, size, var_name, file_name, line_num);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reallocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate reallocation function based on
 * the requested mode, and allows introspection of the allocated memory.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * If the allocation parameters are unchanged, no actual reallocation
 * occurs on the host.
 *
 * If the device uses a separate allocation, it is freed, and a new
 * allocation is delayed (as per initial allocation) so as to invalidate copies
 * which will not be up to date anymore after the associated values
 * modification.
 *
 * \param [in]  ptr        pointer to previously allocated memory
 * \param [in]  mode       allocation mode
 * \param [in]  ni         number of elements
 * \param [in]  size       element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void *
cs_realloc_hd(void            *ptr,
              cs_alloc_mode_t  mode,
              size_t           ni,
              size_t           size,
              const char      *var_name,
              const char      *file_name,
              int              line_num);

#else

inline static void *
cs_realloc_hd(void             *ptr,
              cs_alloc_mode_t   mode,
              size_t            ni,
              size_t            size,
              const char       *var_name,
              const char       *file_name,
              int               line_num)
{
  CS_UNUSED(mode);
  return bft_mem_realloc(ptr, ni, size, var_name, file_name, line_num);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given host pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \param [in]  ptr        pointer to free
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_free_hd(void         *ptr,
           const char   *var_name,
           const char   *file_name,
           int           line_num);

#else

inline static void
cs_free_hd(void         *ptr,
           const char   *var_name,
           const char   *file_name,
           int           line_num)
{
  bft_mem_free(ptr, var_name, file_name, line_num);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given pointer.
 *
 * Compared to \cs_free_hd, this function also allows freeing memory
 * allocated through BFT_MEM_MALLOC / bft_mem_malloc.
 *
 * \param [in]  ptr        pointer to free
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_free(void        *ptr,
        const char  *var_name,
        const char  *file_name,
        int          line_num);

#else

inline static void
cs_free(void         *ptr,
        const char   *var_name,
        const char   *file_name,
        int           line_num)
{
  bft_mem_free(ptr, var_name, file_name, line_num);
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
 * \brief Synchronize data from host to device.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * Depending on the allocation type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * This function assumes the provided pointer was allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  host pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_sync_h2d(const void  *ptr);

#else

static inline void
cs_sync_h2d(const void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initiate synchronization of data from host to device for
 *        future access.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 * In this case, synchronization is started (asynchronously
 * if the allocation mode supports it).
 *
 * In other cases, synchronization will be delayed until actual use.
 * the host pointer should be used with this function.
 *
 * Depending on the allocation type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * This function assumes the provided pointer was allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  host pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_sync_h2d_future(const void  *ptr);

#else

static inline void
cs_sync_h2d_future(const void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize data from device to host.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be passed to this
 * function.
 *
 * Depending on the allocaton type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * This function assumes the provided pointer was allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_sync_d2h(void  *ptr);

#else

static inline void
cs_sync_d2h(void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prefetch data from host to device.
 *
 * This function should only be used on arrays using shared host and device
 * memory, shuch as those allocated using CS_ALLOC_HOST_DEVICE_SHARED.
 * It should be usable on a subset of such an array.
 *
 * \param [in, out]  ptr   pointer to data to prefetch
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_prefetch_h2d(void    *ptr,
                size_t   size);

#else

static inline void
cs_prefetch_h2d(void    *ptr,
                size_t   size)
{
  CS_UNUSED(ptr);
  CS_UNUSED(size);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prefetch data from device to host.
 *
 * This function should only be used on arrays using shared host and device
 * memory, shuch as those allocated using CS_ALLOC_HOST_DEVICE_SHARED.
 * It should be usable on a subset of such an array.
 *
 * \param [in, out]  ptr   pointer to data to prefetch
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_prefetch_d2h(void    *ptr,
                size_t   size);

#else

static inline void
cs_prefetch_d2h(void    *ptr,
                size_t   size)
{
  CS_UNUSED(ptr);
  CS_UNUSED(size);
}

#endif

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on device
 * \param [in, out]  src   pointer to source data on host
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_h2d(void        *dest,
            const void  *src,
            size_t       size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to host.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on host
 * \param [in, out]  src   pointer to source data on device
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_d2h(void        *dest,
            const void  *src,
            size_t       size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to device.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on host
 * \param [in, out]  src   pointer to source data on device
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_d2d(void        *dest,
            const void  *src,
            size_t       size);

#endif /* defined(HAVE_ACCEL) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of host-device allocations
 *
 * \returns current number of host-device allocations.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

int
cs_get_n_allocations_hd(void);

#else

static inline int
cs_get_n_allocations_hd(void)
{
  return 0;
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a given host pointer is allocated with associated with
 *        cs_alloc_hd or cs_realloc_hd.
 *
 * \returns allocated memory size, or zero if not allocated with this
 *          mechanism.
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

size_t
cs_get_allocation_hd_size(void  *host_ptr);

#else

static inline size_t
cs_get_allocation_hd_size(void  *host_ptr)
{
  CS_UNUSED(host_ptr);
  return 0;
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

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_ACCEL_H__ */
