/*============================================================================
 * Memory allocation wrappers for HIP
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "assert.h"
#include "bft/bft_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_mem_hip_priv.h"

/*----------------------------------------------------------------------------*/

/*!
  \file
  \brief Memory allocation wrappers for HIP.

   This file should only contain functions called in cs_mem.cpp.
   In most cases, HIP functions could be called directly from that file,
   but in cases where the code is configured with "--disable_hip_cpp",
   so that only functions explicitely in .cu files use HIP, we still
   need a separate file (or to pass HIP headers).
   If/when this mode is removed, the wrappers defined here can be
   replaced by direct calls in cs_mem.cpp.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

#undef CS_HIP_CHECK
#undef CS_HIP_CHECK_CALL

#define CS_HIP_CHECK(a) { \
  hipError_t _l_ret_code = a; \
  if (hipSuccess != _l_ret_code) { \
    bft_error(__FILE__, __LINE__, 0, "[HIP error] %d: %s\n  running: %s", \
              _l_ret_code, ::hipGetErrorString(_l_ret_code), #a); \
  } \
}

#define CS_HIP_CHECK_CALL(a, file_name, line_num) { \
  hipError_t _l_ret_code = a; \
  if (hipSuccess != _l_ret_code) { \
    bft_error(file_name, line_num, 0, "[HIP error] %d: %s\n  running: %s", \
              _l_ret_code, ::hipGetErrorString(_l_ret_code), #a); \
  } \
}

/*============================================================================
 * Static global variables
 *============================================================================*/

static hipStream_t _cs_glob_stream_pf = 0;
static int          _cs_glob_hip_device_id = -1;

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
cs_mem_hip_set_prefetch_stream(hipStream_t  stream)
{
  _cs_glob_stream_pf = stream;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of HIP device memory.
 *
 * This function simply wraps hipMalloc, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in hip/hip_runtime.h requires
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
cs_mem_hip_malloc_device(size_t        n,
                         const char   *var_name,
                         const char   *file_name,
                         int           line_num)
{
  void *ptr = nullptr;

  CS_HIP_CHECK_CALL(hipMalloc(&ptr, n), file_name, line_num);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of host memory using HIP.
 *
 * This function simply wraps hipHostMalloc, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in hip/hip_runtime.h requires
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
cs_mem_hip_malloc_host(size_t        n,
                       const char   *var_name,
                       const char   *file_name,
                       int           line_num)
{
  void *ptr = nullptr;

  CS_HIP_CHECK_CALL(hipHostMalloc(&ptr, n), file_name, line_num);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of HIP managed memory.
 *
 * This function simply wraps hipMallocManaged, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in hip/hip_runtime.h requires
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
cs_mem_hip_malloc_managed(size_t        n,
                          const char   *var_name,
                          const char   *file_name,
                          int           line_num)
{
  void *ptr = nullptr;

  CS_HIP_CHECK_CALL(hipMallocManaged(&ptr, n), file_name, line_num);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free HIP memory associated with a given pointer.
 *
 * This function simply wraps hipFree, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in hip/hip_runtime.h requires
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
cs_mem_hip_free(void         *p,
                const char   *var_name,
                const char   *file_name,
                int           line_num)
{
  CS_HIP_CHECK_CALL(hipFree(p), file_name, line_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free HIP-allocated host memory associated with a given pointer.
 *
 * This function simply wraps hipHostFree, which could probably be
 * directly called from C or C++, but whose use in such manner is not
 * well documented, and whose declaration in hip/hip_runtime.h requires
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
cs_mem_hip_free_host(void         *p,
                     const char   *var_name,
                     const char   *file_name,
                     int           line_num)
{
  CS_HIP_CHECK_CALL(hipHostFree(p), file_name, line_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over hipMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_hip_copy_h2d(void        *dst,
                    const void  *src,
                    size_t       size)
{
  CS_HIP_CHECK(hipMemcpy(dst, src, size, hipMemcpyHostToDevice));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device, possibly returning on the host
 *        before the copy is finished.
 *
 * This is simply a wrapper over hipMemcpyAsync.
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
cs_mem_hip_copy_h2d_async(void        *dst,
                          const void  *src,
                          size_t       size)
{
  CS_HIP_CHECK(hipMemcpyAsync(dst, src, size, hipMemcpyHostToDevice,
                _cs_glob_stream_pf));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to host.
 *
 * This is simply a wrapper over hipMemcpy.
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
cs_mem_hip_copy_d2h(void        *dst,
                    const void  *src,
                    size_t       size)
{
  CS_HIP_CHECK(hipMemcpy(dst, src, size, hipMemcpyDeviceToHost));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over hipMemcpy.
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
cs_mem_hip_copy_d2h_async(void        *dst,
                          const void  *src,
                          size_t       size)
{
  CS_HIP_CHECK(hipMemcpyAsync(dst, src, size, hipMemcpyDeviceToHost,
                                _cs_glob_stream_pf));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This is simply a wrapper over hipMemPrefetchAsync.
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
cs_mem_hip_prefetch_h2d(const void  *dst,
                        size_t       size)
{
 if (_cs_glob_hip_device_id < 0)
   CS_HIP_CHECK(hipGetDevice(&_cs_glob_hip_device_id));

  CS_HIP_CHECK(hipMemPrefetchAsync(dst, size, _cs_glob_hip_device_id,
                                     _cs_glob_stream_pf));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to host.
 *
 * This is simply a wrapper over hipMemPrefetchAsync.
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
cs_mem_hip_prefetch_d2h(const void  *dst,
                        size_t       size)
{
  CS_HIP_CHECK(hipMemPrefetchAsync(dst, size, hipCpuDeviceId,
                                     _cs_glob_stream_pf));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to device.
 *
 * This is simply a wrapper over hipMemcpy.
 *
 * A safety check is added.
 *
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_hip_copy_d2d(void        *dst,
                    const void  *src,
                    size_t       size)
{
  CS_HIP_CHECK(hipMemcpy(dst, src, size, hipMemcpyDeviceToDevice));
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
cs_mem_hip_set_advise_read_mostly(const void  *ptr,
                                  size_t       size)
{
 if (_cs_glob_hip_device_id < 0)
   CS_HIP_CHECK(hipGetDevice(&_cs_glob_hip_device_id));

  CS_HIP_CHECK(hipMemAdvise(ptr,
                              size,
                              hipMemAdviseSetReadMostly,
                              _cs_glob_hip_device_id));
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
cs_mem_hip_unset_advise_read_mostly(const void  *ptr,
                                    size_t       size)
{
 if (_cs_glob_hip_device_id < 0)
   CS_HIP_CHECK(hipGetDevice(&_cs_glob_hip_device_id));

  CS_HIP_CHECK(hipMemAdvise(ptr,
                              size,
                              hipMemAdviseUnsetReadMostly,
                              _cs_glob_hip_device_id));
}

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
cs_mem_hip_is_device_ptr(const void  *ptr)
{
  if (ptr == nullptr)
    return true;

  hipPointerAttribute_t attributes;

  int retcode = hipPointerGetAttributes(&attributes, ptr);

  if (retcode == hipSuccess) {
    if (ptr == attributes.devicePointer)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Get memory usage on HIP device (in kB)
 *
 * \return[in]  memory usage on current HIP device
 */
/*----------------------------------------------------------------------------*/

size_t
cs_mem_hip_get_device_memory_usage(void)
{
  size_t retval = 0;

  size_t mem_free = 0, mem_total = 0;
  int retcode = hipMemGetInfo(&mem_free, &mem_total);
  if (retcode == hipSuccess) {
    retval = mem_total - mem_free;
    retval /= 1024;
  }

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
