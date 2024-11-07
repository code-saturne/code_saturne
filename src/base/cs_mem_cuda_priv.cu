/*============================================================================
 * Memory allocation wrappers for CUDA
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "assert.h"
#include "bft_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mem_cuda_priv.h"

/*----------------------------------------------------------------------------*/

/*!
  \file cs_mem_cuda_priv.cu
        Memory allocation wrappers for CUDA.

   This file should only contain functions called in cs_mem.cpp.
   In most cases, CUDA functions could be called directly from that file,
   but in cases where the code is configured with "--disable_cuda_cpp",
   so that only functions explicitely in .cu files use CUDA, we still
   need a separate file (or to pass CUDA headers).
   If/when this mode is removed, the wrappers defined here can be
   replaced by direct calls in cs_mem.cpp.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

#undef CS_CUDA_CHECK
#undef CS_CUDA_CHECK_CALL

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

/*============================================================================
 * Static global variables
 *============================================================================*/

static cudaStream_t _cs_glob_stream_pf = 0;
static int          _cs_glob_cuda_device_id = -1;

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_mem.cpp, and
 * not directly by the user.
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define stream to use for prefetch
 *
 * \param [in]  stream  stream to use for prefetching
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_cuda_set_prefetch_stream(cudaStream_t  stream)
{
  _cs_glob_stream_pf = stream;
}

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
cs_mem_cuda_malloc_device(size_t        n,
                          const char   *var_name,
                          const char   *file_name,
                          int           line_num)
{
  void *ptr = nullptr;

  CS_CUDA_CHECK_CALL(cudaMalloc(&ptr, n), file_name, line_num);

  return ptr;
}

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
cs_mem_cuda_malloc_host(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num)
{
  void *ptr = nullptr;

  CS_CUDA_CHECK_CALL(cudaMallocHost(&ptr, n), file_name, line_num);

  return ptr;
}

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
cs_mem_cuda_malloc_managed(size_t        n,
                           const char   *var_name,
                           const char   *file_name,
                           int           line_num)
{
  void *ptr = nullptr;

  CS_CUDA_CHECK_CALL(cudaMallocManaged(&ptr, n), file_name, line_num);

  return ptr;
}

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
cs_mem_cuda_free(void         *p,
                 const char   *var_name,
                 const char   *file_name,
                 int           line_num)
{
  CS_CUDA_CHECK_CALL(cudaFree(p), file_name, line_num);
}

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
 * \param [in]  p          pointer to device memory
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_cuda_free_host(void         *p,
                      const char   *var_name,
                      const char   *file_name,
                      int           line_num)
{
  CS_CUDA_CHECK_CALL(cudaFreeHost(p), file_name, line_num);
}

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
cs_mem_cuda_copy_h2d(void        *dst,
                     const void  *src,
                     size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
}

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
cs_mem_cuda_copy_h2d_async(void        *dst,
                           const void  *src,
                           size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice));
}

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
cs_mem_cuda_copy_d2h(void        *dst,
                     const void  *src,
                     size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
}

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
cs_mem_cuda_copy_d2h_async(void        *dst,
                           const void  *src,
                           size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpyAsync(dst, src, size, cudaMemcpyDeviceToHost));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over cudaMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to data
 * \param [in]   size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_cuda_prefetch_h2d(const void  *dst,
                         size_t       size)
{
 if (_cs_glob_cuda_device_id < 0)
   cudaGetDevice(&_cs_glob_cuda_device_id);

  CS_CUDA_CHECK(cudaMemPrefetchAsync(dst, size, _cs_glob_cuda_device_id, \
                                     _cs_glob_stream_pf));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to host.
 *
 * This is simply a wrapper over cudaMemcpy.
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
cs_mem_cuda_prefetch_d2h(const void  *dst,
                         size_t       size)
{
  CS_CUDA_CHECK(cudaMemPrefetchAsync(dst, size, cudaCpuDeviceId, \
                                     _cs_glob_stream_pf));
}

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
cs_mem_cuda_copy_d2d(void        *dst,
                     const void  *src,
                     size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]    ptr   pointer to allocation
 * \param [size]  size   associated data size
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_cuda_set_advise_read_mostly(const void  *ptr,
                                   size_t       size)
{
 if (_cs_glob_cuda_device_id < 0)
   cudaGetDevice(&_cs_glob_cuda_device_id);

  CS_CUDA_CHECK(cudaMemAdvise(ptr,
                              size,
                              cudaMemAdviseSetReadMostly,
                              _cs_glob_cuda_device_id))
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]    ptr   pointer to allocation
 * \param [size]  size   associated data size
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_cuda_unset_advise_read_mostly(const void  *ptr,
                                     size_t       size)
{
 if (_cs_glob_cuda_device_id < 0)
   cudaGetDevice(&_cs_glob_cuda_device_id);

  CS_CUDA_CHECK(cudaMemAdvise(ptr,
                              size,
                              cudaMemAdviseUnsetReadMostly,
                              _cs_glob_cuda_device_id))
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
