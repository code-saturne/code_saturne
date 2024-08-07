#ifndef __CS_BASE_CUDA_H__
#define __CS_BASE_CUDA_H__

/*============================================================================
 * Definitions, global variables, and base functions for CUDA
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

#if defined(HAVE_CUDA)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base_accel.h"
#include "cs_log.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CUDA_CHECK(a) { \
    cudaError_t _l_ret_code = a; \
    if (cudaSuccess != _l_ret_code) { \
      bft_error(__FILE__, __LINE__, 0, "[CUDA error] %d: %s\n  running: %s", \
                _l_ret_code, ::cudaGetErrorString(_l_ret_code), #a); \
    } \
  }

#define CS_CUDA_CHECK_CALL(a, file_name, line_num) { \
    cudaError_t _l_ret_code = a; \
    if (cudaSuccess != _l_ret_code) { \
      bft_error(file_name, line_num, 0, "[CUDA error] %d: %s\n  running: %s", \
                _l_ret_code, ::cudaGetErrorString(_l_ret_code), #a); \
    } \
  }

/* For all current compute capabilities, the warp size is 32; If it ever
   changes, it can be obtained through cudaDeviceProp, so we could then
   replace this macro with a global variable */

#define CS_CUDA_WARP_SIZE 32

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

extern int  cs_glob_cuda_device_id;

/* Other device parameters */

extern int  cs_glob_cuda_max_threads_per_block;
extern int  cs_glob_cuda_max_block_size;
extern int  cs_glob_cuda_max_blocks;
extern int  cs_glob_cuda_n_mp;  /* Number of multiprocessors */

/* Allow graphs for kernel launches ? May interfere with profiling (nsys),
   so can be deactivated. */

extern bool cs_glob_cuda_allow_graph;

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_base_accel.c, and
 * not directly by the user.
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Allocate n bytes of CUDA device memory.
 *
 * This function simply wraps cudaMalloc, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
 * \param [in]  n          element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_cuda_mem_malloc_device(size_t        n,
                          const char   *var_name,
                          const char   *file_name,
                          int           line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Allocate n bytes of host memory using CUDA.
 *
 * This function simply wraps cudaMallocHost, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
 * \param [in]  n          element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_cuda_mem_malloc_host(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Allocate n bytes of CUDA managed memory.
 *
 * This function simply wraps cudaMallocManaged, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
 * \param [in]  n          element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_cuda_mem_malloc_managed(size_t        n,
                           const char   *var_name,
                           const char   *file_name,
                           int           line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Free CUDA memory associated with a given pointer.
 *
 * This function simply wraps cudaFree, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
 * \param [in]  p          pointer to device memory
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_mem_free(void         *p,
                 const char   *var_name,
                 const char   *file_name,
                 int           line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Free CUDA-allocated host memory associated with a given pointer.
 *
 * This function simply wraps cudaFreeHost, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
 * \param [in]  p          pointer to device memory
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_mem_free_host(void         *p,
                      const char   *var_name,
                      const char   *file_name,
                      int           line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over cudaMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_copy_h2d(void         *dst,
                 const void   *src,
                 size_t        size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy data from host to device, possibly returning on the host
 *        before the copy is finished.
 *
 * This is simply a wrapper over cudaMemcpyAsync.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_copy_h2d_async(void        *dst,
                       const void  *src,
                       size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy data from device to host.
 *
 * This is simply a wrapper over cudaMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_copy_d2h(void        *dst,
                 const void  *src,
                 size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over cudaMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_copy_d2h_async(void        *dst,
                       const void  *src,
                       size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Prefetch data from host to device.
 *
 * This is simply a wrapper over cudaMemPrefetchAsync.
 *
 * A safety check is added.
 *
 * \param [in]  dst   pointer to data
 * \param [in]  size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_prefetch_h2d(const void  *dst,
                     size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Prefetch data from device to host.
 *
 * This is simply a wrapper over cudaMemPrefetchAsync.
 *
 * A safety check is added.
 *
 * \param [in]  dst   pointer to data
 * \param [in]  size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_prefetch_d2h(const void  *dst,
                     size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Copy data from device to device.
 *
 * This is simply a wrapper over cudaMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_copy_d2d(void        *dst,
                 const void  *src,
                 size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Get host pointer for a managed or device pointer.
 *
 * This function can be called with a pointer inside an allocated block of
 * memory, so is not retricted to values returned by CS_ALLOC_HD.
 *
 * This makes it possible to check whether a pointer to an array inside
 * a larger array is shared or accessible from the device only
 * (for example when grouping allocations).
 *
 * \param [in]   ptr   pointer to device data
 *
 * \return  pointer to host data if shared or mapped at the CUDA level,
 *          NULL otherwise.
 */
/*----------------------------------------------------------------------------*/

void *
cs_cuda_get_host_ptr(const void  *ptr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a pointer is a device (or shared) pointer.
 *
 * \param [in]   ptr   pointer to device data
 *
 * \return  true if the pointer is usable from the device or null, false
 *          if available on host only or if query failed.
 */
/*----------------------------------------------------------------------------*/

bool
cs_cuda_is_device_ptr(const void  *ptr);

/*----------------------------------------------------------------------------*/
/*
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]    ptr   pointer to allocation
 * \param [size]  size  associated data size
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_mem_set_advise_read_mostly(const void  *ptr,
                                   size_t       size);

/*----------------------------------------------------------------------------*/
/*
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]    ptr    pointer to allocation
 * \param [size]  size   associated data size
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_mem_unset_advise_read_mostly(const void  *ptr,
                                   size_t       size);

/*=============================================================================
 * Inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute grid size for given array and block sizes.
 *
 * This assumes each thread on a given block handles a single array element.
 * For kernels in which each thread handles multiple elements, a grid size
 * divided by that multiple is sufficient.
 *
 * \param[in]  n           size of arrays
 * \param[in]  block_size  block size for kernels
 *
 * \return  grid size for kernels
 */
/*----------------------------------------------------------------------------*/

static inline unsigned int
cs_cuda_grid_size(cs_lnum_t     n,
                  unsigned int  block_size)
{
  return (n % block_size) ?  n/block_size + 1 : n/block_size;
}

END_C_DECLS

#if defined(__NVCC__)

/*----------------------------------------------------------------------------
 * Synchronize of copy a cs_real_t type array from the host to a device.
 *
 * parameters:
 *   val_h          <-- pointer to host data
 *   n_vals         <-- number of data values
 *   device_id      <-- associated device id
 *   stream         <-- associated stream (for async prefetch only)
 *   val_d          --> matching pointer on device
 *   buf_d          --> matching allocation pointer on device (should be freed
 *                      after use if non-NULL)
 *----------------------------------------------------------------------------*/

template <typename T>
void
cs_sync_or_copy_h2d(const T        *val_h,
                    cs_lnum_t       n_vals,
                    int             device_id,
                    cudaStream_t    stream,
                    const T       **val_d,
                    void          **buf_d)
{
  const T  *_val_d = NULL;
  void     *_buf_d = NULL;

  if (val_h != NULL) {

    cs_alloc_mode_t alloc_mode = cs_check_device_ptr(val_h);
    size_t size = n_vals * sizeof(T);

    if (alloc_mode == CS_ALLOC_HOST) {
      CS_CUDA_CHECK(cudaMalloc(&_buf_d, size));
      cs_cuda_copy_h2d(_buf_d, val_h, size);
      _val_d = (const T *)_buf_d;
    }
    else {
      _val_d = (const T *)cs_get_device_ptr_const((const void *)val_h);

      if (alloc_mode != CS_ALLOC_HOST_DEVICE_SHARED)
        cs_sync_h2d(val_h);
    }

  }

  *val_d = _val_d;
  *buf_d = _buf_d;
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Return stream handle from stream pool.
 *
 * If the requested stream id is higher than the current number of streams,
 * one or more new streams will be created, so that size of the stream pool
 * matches at least stream_id+1.
 *
 * By default, the first stream (with id 0) will be used for most operations,
 * while stream id 1 will be used for operations which can be done
 * concurrently, such as memory prefetching.
 *
 * Additional streams can be used for independent tasks, though opportunities
 * for this are limited in the current code (this would probably also require
 * associating different MPI communicators with each task).
 *
 * \param [in]  stream_id  id or requested stream
 *
 * \returns handle to requested stream
 */
/*----------------------------------------------------------------------------*/

cudaStream_t
cs_cuda_get_stream(int  stream_id);

#endif /* defined(__NVCC__) */

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * \brief  Log information on available CUDA devices.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_device_info(cs_log_t  log_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Log information on available CUDA version.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_version_info(cs_log_t  log_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Log information on CUDA compiler.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_compiler_info(cs_log_t  log_id);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set CUDA device based on MPI rank and number of devices.
 *
 * \param[in]  comm            associated MPI communicator
 * \param[in]  ranks_per_node  number of ranks per node (min and max)
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_base_cuda_select_default_device(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return currently selected CUDA devices.
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_base_cuda_get_device(void);

#endif  /* CS_HAVE_CUDA */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_CUDA_H__ */
