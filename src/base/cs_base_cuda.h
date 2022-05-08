#ifndef __CS_BASE_CUDA_H__
#define __CS_BASE_CUDA_H__

/*============================================================================
 * Definitions, global variables, and base functions for CUDA
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

#include "cs_log.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CUDA_CHECK(a) { \
    cudaError_t ret_code = a; \
    if (cudaSuccess != ret_code) { \
      bft_error(__FILE__, __LINE__, 0, "[CUDA error] %d: %s\n  running: %s", \
                ret_code, ::cudaGetErrorString(ret_code), #a); \
    } \
  }

#define CS_CUDA_CHECK_CALL(a, file_name, line_num) { \
    cudaError_t ret_code = a; \
    if (cudaSuccess != ret_code) { \
      bft_error(file_name, line_num, 0, "[CUDA error] %d: %s\n  running: %s", \
                ret_code, ::cudaGetErrorString(ret_code), #a); \
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

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_base_accel.c, and
 * not directly by the user.
 *============================================================================*/

#if defined(HAVE_CUDA)

/*----------------------------------------------------------------------------*/
/*!
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
/*!
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
/*!
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
/*!
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
/*!
 * \brief Free CUDA-allocated host memory associated with a given pointer.
 *
 * This function simply wraps cudaFreeHost, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in cuda_runtime.h requires
 * support of function attributes by compiler.
 *
 * A safety check is added.
 *
$ * \param [in]  p          pointer to device memory
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
/*!
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
/*!
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
/*!
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
/*!
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
/*!
 * \brief Prefetch data from host to device.
 *
 * This is simply a wrapper over cudaMemPrefetchAsync.
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
cs_cuda_prefetch_h2d(void    *dst,
                     size_t   size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prefetch data from device to host.
 *
 * This is simply a wrapper over cudaMemPrefetchAsync.
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
cs_cuda_prefetch_d2h(void    *dst,
                     size_t   size);

/*----------------------------------------------------------------------------*/
/*!
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

#endif

/*=============================================================================
 * Inline static function prototypes
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
  return (n % block_size) ?  n/block_size : n/block_size + 1;
}

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on available CUDA version.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_version_info(cs_log_t  log_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on CUDA compiler.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_compiler_info(cs_log_t  log_id);

/*----------------------------------------------------------------------------*/
/*!
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
/*!
 * \brief Return currently selected CUDA devices.
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_base_cuda_get_device(void);

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_CUDA_H__ */
