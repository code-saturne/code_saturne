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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "assert.h"
#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base_cuda.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Keep track of active device id; usually queried dynamically, but
   saving the value in this variable can be useful when debugging */

int  cs_glob_cuda_device_id = -1;

/* Other device parameters */

int  cs_glob_cuda_max_threads_per_block = -1;
int  cs_glob_cuda_max_block_size = -1;
int  cs_glob_cuda_max_blocks = -1;
int  cs_glob_cuda_n_mp = -1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_base_accel.c, and
 * not directly by the user.
 *============================================================================*/

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
                          int           line_num)
{
  void *ptr = NULL;

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
cs_cuda_mem_malloc_host(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num)
{
  void *ptr = NULL;

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
cs_cuda_mem_malloc_managed(size_t        n,
                           const char   *var_name,
                           const char   *file_name,
                           int           line_num)
{
  void *ptr = NULL;

  CS_CUDA_CHECK_CALL(cudaMallocManaged(&ptr, n), file_name, line_num);

#if 0
  CS_CUDA_CHECK_CALL(cudaMemPrefetchAsync (*pointer, size, cudaCpuDeviceId, 0),
                     file_name, line_num);
  CS_CUDA_CHECK_CALL(cudaDeviceSynchronize(), file_name, line_num);
#endif

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
cs_cuda_mem_free(void         *p,
                 const char   *var_name,
                 const char   *file_name,
                 int           line_num)
{
  CS_CUDA_CHECK_CALL(cudaFree(p), file_name, line_num);

#if 0
  CS_CUDA_CHECK_CALL((cudaDeviceSynchronize(), file_name, line_num);
#endif
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
cs_cuda_mem_free_host(void         *p,
                      const char   *var_name,
                      const char   *file_name,
                      int           line_num)
{
  CS_CUDA_CHECK_CALL(cudaFreeHost(p), file_name, line_num);

#if 0
  CS_CUDA_CHECK_CALL((cudaDeviceSynchronize(), file_name, line_num);
#endif
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
cs_cuda_copy_h2d(void        *dst,
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
cs_cuda_copy_h2d_async(void        *dst,
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
cs_cuda_copy_d2h(void        *dst,
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
cs_cuda_copy_d2h_async(void        *dst,
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
 * \param [out]  dst   pointer to destination data
 * \param [in]   src   pointer to source data
 * \param [in]   size  size of data to copy
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_prefetch_h2d(void    *dst,
                     size_t   size)
{
  CS_CUDA_CHECK(cudaMemPrefetchAsync(dst, size, cudaCpuDeviceId, 0));
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
cs_cuda_prefetch_d2h(void    *dst,
                     size_t   size)
{
  CS_CUDA_CHECK(cudaMemPrefetchAsync(dst, size, cs_glob_cuda_device_id, 0));
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
cs_cuda_copy_d2d(void        *dst,
                 const void  *src,
                 size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on available CUDA devices.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_device_info(cs_log_t  log_id)
{
  int n_devices = 0;

  cudaError_t retval = cudaGetDeviceCount(&n_devices);

  if (retval == cudaErrorNoDevice)
    cs_log_printf(log_id,
                  _("  CUDA device:         none available\n"));
  else if (retval)
    cs_log_printf(log_id,
                  _("  CUDA device:         %s\n"),
                  cudaGetErrorString(retval));

  char buffer[256] = "";

  for (int i = 0; i < n_devices; i++) {
    struct cudaDeviceProp prop;
    CS_CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
    unsigned long long mem = prop.totalGlobalMem / 1000000;
    char mode_name[32] = "";
    if (prop.computeMode == cudaComputeModeDefault)
      snprintf(mode_name, 31, "default");
    else if (prop.computeMode == cudaComputeModeExclusive)
      snprintf(mode_name, 31, "exclusive");
    else if (prop.computeMode == cudaComputeModeProhibited)
      snprintf(mode_name, 31, "prohibited");

    cs_log_printf
      (log_id,
       _("  CUDA device %d:       %s\n"),
       i, prop.name);

    if (strncmp(prop.name, buffer, 255) != 0)
      cs_log_printf
        (log_id,
         _("                       Compute capability: %d.%d\n"
           "                       Memory: %llu %s\n"
           "                       Multiprocessors: %d\n"
           "                       Integrated: %d\n"
           "                       Can map host memory: %d\n"
           "                       Compute mode: %s\n"),
         prop.major, prop.minor,
         mem, _("MB"),
         prop.multiProcessorCount,
         prop.integrated,
         prop.canMapHostMemory, mode_name);

    strncpy(buffer, prop.name, 255);
    buffer[255] = '\0';
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on available CUDA version.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_version_info(cs_log_t  log_id)
{
  int runtime_version = -1, driver_version = -1;

  if (cudaDriverGetVersion(&driver_version) == cudaSuccess)
    cs_log_printf(log_id,
                  "  %s%d\n", _("CUDA driver:         "), driver_version);
  if (cudaRuntimeGetVersion(&runtime_version) == cudaSuccess)
    cs_log_printf(log_id,
                  "  %s%d\n", _("CUDA runtime:        "), runtime_version);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on CUDA compiler.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_compiler_info(cs_log_t  log_id)
{
  cs_log_printf(log_id,
                "    %s%d.%d.%d\n", _("CUDA compiler:     "),
                __CUDACC_VER_MAJOR__,
                __CUDACC_VER_MINOR__,
                __CUDACC_VER_BUILD__);
}

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
cs_base_cuda_select_default_device(void)
{
  int device_id = 0, n_devices = 0;

  cudaError_t ret_code = cudaGetDeviceCount(&n_devices);

  if (ret_code == cudaErrorNoDevice)
    return -1;

  if (cudaSuccess != ret_code) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf("[CUDA error] %d: %s\n  running: %s\n  in: %s\n",
               ret_code, ::cudaGetErrorString(ret_code),
               "cudaGetDeviceCount", __func__);
    return -1;
  }

  if (cs_glob_rank_id > -1 && n_devices > 1) {

    device_id = cs_glob_node_rank_id*n_devices / cs_glob_node_n_ranks;

    assert(device_id > -1 && device_id < n_devices);

  }

  ret_code = cudaSetDevice(device_id);

  if (cudaSuccess != ret_code) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf("[CUDA error] %d: %s\n  running: %s\n  in: %s\n",
               ret_code, ::cudaGetErrorString(ret_code),
               "cudaSetDevice", __func__);
    return -1;
  }

  cs_glob_cuda_device_id = device_id;

  /* Also query some device properties */

  struct cudaDeviceProp prop;
  CS_CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id));
  cs_glob_cuda_max_threads_per_block = prop.maxThreadsPerBlock;
  cs_glob_cuda_max_block_size = prop.maxThreadsPerMultiProcessor;
  cs_glob_cuda_max_blocks
    =   prop.multiProcessorCount
      * (prop.maxThreadsPerMultiProcessor / prop.maxThreadsPerBlock);
  cs_glob_cuda_n_mp = prop.multiProcessorCount;

  return device_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return currently selected CUDA devices.
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

int
cs_base_cuda_get_device(void)
{
  int device_id = -1, n_devices = 0;

  cudaError_t ret_code = cudaGetDeviceCount(&n_devices);

  if (cudaSuccess == ret_code)
    ret_code = cudaGetDevice(&device_id);

  if (cudaSuccess != ret_code)
    device_id = -1;

  return device_id;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
