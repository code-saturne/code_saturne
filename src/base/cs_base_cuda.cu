/*============================================================================
 * Definitions, global variables, and base functions for CUDA
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

#if defined(HAVE_NCCL)
#include <nccl.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "assert.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "base/cs_mem_cuda_priv.h"

#if defined(HAVE_NCCL)
#include "base/cs_fp_exception.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_base_cuda.h"

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

/* Stream pool */

static int            _cs_glob_cuda_n_streams = -1;
static cudaStream_t  *_cs_glob_cuda_streams = nullptr;

#if defined(HAVE_NCCL)

ncclComm_t    cs_glob_nccl_comm = nullptr;

/* NCCL Datatypes associated with code_saturne datatypes */

ncclDataType_t  cs_datatype_to_nccl[] = {ncclChar,
                                         ncclChar,
                                         ncclFloat32,
                                         ncclFloat64,
                                         ncclFloat16,   /* No match (avoid) */
                                         ncclInt32,     /* CS_INT32 */
                                         ncclInt64,     /* CS_INT64 */
                                         ncclUint32,    /* CS_UINT32 */
                                         ncclUint64};   /* CS_UINT64 */

#endif

/* Reduce buffers associated with streams in pool */

static unsigned *_r_elt_size = nullptr;
static unsigned *_r_grid_size = nullptr;
static void  **_r_host = nullptr;
static void  **_r_reduce = nullptr;
static void  **_r_grid = nullptr;

static cudaStream_t _cs_glob_stream_pf = 0;

/* Allow graphs for kernel launches ? May interfere with profiling (nsys),
   so can be deactivated. */

bool cs_glob_cuda_allow_graph = false;

// Shared memory size per block (based on know GPUs, queried later).
int cs_glob_cuda_shared_mem_per_block = 48*1024;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy stream pool at exit.
 */
/*----------------------------------------------------------------------------*/

static void
finalize_streams_(void)
{
  cs_mem_cuda_set_prefetch_stream(0);
  cudaStreamDestroy(_cs_glob_stream_pf);

  for (int i = 0; i < _cs_glob_cuda_n_streams; i++) {
    cudaStreamDestroy(_cs_glob_cuda_streams[i]);
    CS_FREE(_r_host[i]);
    // _r_reduce is based on _r_host's device buffer
    CS_FREE(_r_grid[i]);
  }

  CS_FREE(_cs_glob_cuda_streams);

  CS_FREE(_r_elt_size);
  CS_FREE(_r_grid_size);
  CS_FREE(_r_host);
  CS_FREE(_r_reduce);
  CS_FREE(_r_grid);

  _cs_glob_cuda_n_streams = 0;
}

#if defined(HAVE_NCCL) && defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize NCCL.
 */
/*----------------------------------------------------------------------------*/

static void
_finalize_nccl(void)
{
  if (cs_glob_nccl_comm != nullptr) {
    ncclResult_t retcode = ncclCommFinalize(cs_glob_nccl_comm);
    if (retcode == ncclSuccess) {
      retcode = ncclCommDestroy(cs_glob_nccl_comm);
      cs_glob_nccl_comm = nullptr;
    }
    if (retcode != ncclSuccess) {
      const char * err_str = ncclGetErrorString(retcode);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "  Failed finalizing NCCL: %s;\n"),
         err_str);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize NCCL.
 *
 * \param [in]   n_devices     number of devices on local node
 * \param [in]   n_node_ranks  number of ranks per node
 */
/*----------------------------------------------------------------------------*/

static void
_initialize_nccl(int  n_devices,
                 int  n_node_ranks)
{
  if (cs_glob_n_ranks < 2)
    return;

  int may_use_nccl = (n_node_ranks <= n_devices) ? 1 : 0;

  const char *s = getenv("CS_DISABLE_NCCL");
  if (s != nullptr) {
    int i = atoi(s);
    if (i < 1)
      may_use_nccl = -1;
  }

  MPI_Allreduce(MPI_IN_PLACE, &may_use_nccl, 1, MPI_INT,
                MPI_MIN, cs_glob_mpi_comm);

  if (may_use_nccl > 0) {
    cs_fp_exception_disable_trap();

    /* Get NCCL identifier */
    ncclUniqueId uid;
    if (cs_glob_rank_id == 0)
      ncclGetUniqueId(&uid);
    MPI_Bcast(&uid, sizeof(uid), MPI_BYTE, 0, cs_glob_mpi_comm);

    /* Initialize NCCL communicator */

    ncclResult_t retcode = ncclCommInitRank(&cs_glob_nccl_comm,
                                            cs_glob_n_ranks,
                                            uid,
                                            cs_glob_rank_id);

    if (cs_glob_nccl_comm == nullptr || retcode != ncclSuccess) {
      const char * err_str = ncclGetErrorString(retcode);
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "  Failed initializing NCCL communicator: %s;\n"
           "    Fallback to MPI.\n"),
         err_str);
    }
    else
      cs_base_at_finalize(_finalize_nccl);

    cs_fp_exception_restore_trap();
  }

  else {
    if (may_use_nccl <  0)
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "  Not using NCCL (CS_DISABLE_NCCL=%s).\n"), s);
    else
      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n"
           "  Not using NCCL (requires 1 device per MPI rank, "
           "currently %d/%d).\n"),
         n_devices, n_node_ranks);
  }

}

#endif

/*============================================================================
 * Semi-private function prototypes
 *
 * The following functions are intended to be used by the common
 * host-device memory management functions from cs_base_accel.c, and
 * not directly by the user.
 *============================================================================*/

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
  CS_CUDA_CHECK(cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice,
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
  CS_CUDA_CHECK(cudaMemcpyAsync(dst, src, size, cudaMemcpyDeviceToHost,
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
cs_cuda_copy_d2d(void        *dst,
                 const void  *src,
                 size_t       size)
{
  CS_CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get host pointer for a managed or device pointer.
 *
 * This function can be called with a pointer inside an allocated block of
 * memory, so is not restricted to values returned by CS_MALLOC_HD.
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
cs_cuda_get_host_ptr(const void  *ptr)
{
  cudaPointerAttributes attributes;

  void *host_ptr = nullptr;
  int retcode = cudaPointerGetAttributes(&attributes, ptr);

  if (retcode == cudaSuccess) {
    if (ptr != attributes.devicePointer)
      bft_error(__FILE__, __LINE__, 0,
                _("%s: %p does not seem to be a managed or device pointer."),
                __func__, ptr);

    host_ptr = attributes.hostPointer;
  }

  return host_ptr;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

END_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

#ifdef __CUDACC__

/*----------------------------------------------------------------------------*/
/*!
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
cs_cuda_get_stream(int  stream_id)
{
  if (stream_id >= 0 && stream_id < _cs_glob_cuda_n_streams)
    return _cs_glob_cuda_streams[stream_id];
  else if (stream_id < 0)
    return nullptr;

  if (_cs_glob_cuda_n_streams < 0) {
    cs_base_at_finalize(finalize_streams_);
    _cs_glob_cuda_n_streams = 0;
  }

  CS_REALLOC(_cs_glob_cuda_streams, stream_id+1, cudaStream_t);
  CS_REALLOC(_r_elt_size, stream_id+1, unsigned);
  CS_REALLOC(_r_grid_size, stream_id+1, unsigned);
  CS_REALLOC(_r_reduce, stream_id+1, void *);
  CS_REALLOC(_r_host, stream_id+1, void *);
  CS_REALLOC(_r_grid, stream_id+1, void *);

  for (int i = _cs_glob_cuda_n_streams; i < stream_id+1; i++) {
    cudaStreamCreate(&_cs_glob_cuda_streams[i]);
    _r_elt_size[i] = 0;
    _r_grid_size[i] = 0;
    _r_host[i] = nullptr;
    _r_reduce[i] = nullptr;
    _r_grid[i] = nullptr;
  }

  _cs_glob_cuda_n_streams = stream_id+1;

  return _cs_glob_cuda_streams[stream_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream handle used for prefetching.
 *
 * By default, a single stream is created specifically for prefetching.
 *
 * \returns handle to prefetching stream
 */
/*----------------------------------------------------------------------------*/

cudaStream_t
cs_cuda_get_stream_prefetch(void)
{
  return _cs_glob_stream_pf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return stream id in stream pool matching a given CUDA stream.
 *
 * If the stream is not present in the stream pool, return -1.
 *
 * \param [in]  handle to given streams
 *
 * \returns if of stream in pool, or -1.
 */
/*----------------------------------------------------------------------------*/

int
cs_cuda_get_stream_id(cudaStream_t  stream)
{
  for (int i = 0; i < _cs_glob_cuda_n_streams; i++)
    if (stream == _cs_glob_cuda_streams[i])
      return i;

  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointers to reduction buffers needed for 2-stage reductions.
 *
 * These buffers are used internally by CUDA 2-stage operations, and are
 * allocated and resized updon demand.
 *
 * \param[in]   stream_id   stream id in pool
 * \param[in]   n_elts      size of arrays
 * \param[in]   n_elts      size of arrays
 * \param[in]   elt_size    size of element or structure simultaneously reduced
 * \param[in]   grid_size   associated grid size
 * \param[out]  r_grid      first stage reduce buffer
 * \param[out]  r_reduce    second stage (final result) reduce buffer
 * \param[out]  r_host      final result reduce buffer on host
 */
/*----------------------------------------------------------------------------*/

void
cs_cuda_get_2_stage_reduce_buffers(int            stream_id,
                                   cs_lnum_t      n_elts,
                                   size_t         elt_size,
                                   unsigned int   grid_size,
                                   void*         &r_grid,
                                   void*         &r_reduce,
                                   void*         &r_host)
{
  assert(stream_id > -1 && stream_id < _cs_glob_cuda_n_streams);

  unsigned int t_grid_size = grid_size * elt_size;

  if (_r_elt_size[stream_id] < elt_size) {
    _r_elt_size[stream_id] = elt_size;
    CS_FREE_HD(_r_host[stream_id]);
    unsigned char *b_ptr;
    CS_MALLOC_HD(b_ptr, elt_size, unsigned char, CS_ALLOC_HOST_DEVICE_PINNED);
    _r_host[stream_id] = b_ptr;
    _r_reduce[stream_id] = cs_get_device_ptr(b_ptr);
  }

  if (_r_grid_size[stream_id] < t_grid_size) {
    _r_grid_size[stream_id] = t_grid_size;
    CS_FREE(_r_grid[stream_id]);
    unsigned char *b_ptr;
    CS_MALLOC_HD(b_ptr, _r_grid_size[stream_id], unsigned char, CS_ALLOC_DEVICE);
    _r_grid[stream_id] = b_ptr;
  }

  r_grid = _r_grid[stream_id];
  r_reduce = _r_reduce[stream_id];
  r_host = _r_host[stream_id];
}

#endif /* defined(__CUDACC__) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log information on available CUDA devices.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

extern "C" void
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

    cs_glob_cuda_shared_mem_per_block = prop.sharedMemPerBlock;

    cs_log_printf
      (log_id,
       _("  CUDA device %d:       %s\n"),
       i, prop.name);

    if (strncmp(prop.name, buffer, 255) != 0) {
      cs_log_printf
        (log_id,
         _("                       Compute capability: %d.%d\n"
           "                       Memory: %llu %s\n"
           "                       Multiprocessors: %d\n"
           "                       Integrated: %d\n"
           "                       Unified addressing: %d\n"),
         prop.major, prop.minor,
         mem, _("MB"),
         prop.multiProcessorCount,
         prop.integrated,
         prop.unifiedAddressing);

#if (CUDART_VERSION >= 11000)
      cs_log_printf
        (log_id,
         _("                       Use host's page tables: %d\n"),
         prop.pageableMemoryAccessUsesHostPageTables);
#endif
    }

    strncpy(buffer, prop.name, 256);
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

extern "C" void
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

extern "C" void
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
 * \brief Log information on NCCL.
 *
 * \param[in]  log_id  id of log file in which to print information
 */
/*----------------------------------------------------------------------------*/

extern "C" void
cs_base_cuda_nccl_info(cs_log_t  log_id)
{
#if defined(HAVE_NCCL)
  int link_version;
  ncclGetVersion(&link_version);
  int compile_version = NCCL_VERSION(NCCL_MAJOR, NCCL_MINOR, NCCL_PATCH);
  if (compile_version != link_version)
    cs_log_printf(log_id,
                  _("  NCCL version %d (compiled with version %d)\n"),
                  link_version, compile_version);
  else
    cs_log_printf(log_id,
                  "  NCCL version %d.%d.%d\n",
                  NCCL_MAJOR, NCCL_MINOR, NCCL_PATCH);
#endif
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

extern "C" int
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

  cs_alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;
  cs_alloc_mode_read_mostly = CS_ALLOC_HOST_DEVICE_SHARED;
  cs_alloc_mode_device = CS_ALLOC_DEVICE;

  const char *s = getenv("CS_CUDA_ALLOC_DEVICE_UVM");
  if (s != nullptr) {
    int i = atoi(s);
    cs_alloc_mode_device = CS_ALLOC_HOST_DEVICE_SHARED;
  }

  /* Also query some device properties */

  struct cudaDeviceProp prop;
  CS_CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id));
  cs_glob_cuda_max_threads_per_block = prop.maxThreadsPerBlock;
  cs_glob_cuda_max_block_size = prop.maxThreadsPerMultiProcessor;
  cs_glob_cuda_max_blocks
    =   prop.multiProcessorCount
      * (prop.maxThreadsPerMultiProcessor / prop.maxThreadsPerBlock);
  cs_glob_cuda_n_mp = prop.multiProcessorCount;

  /* Create default stream for prefetching */
  if (_cs_glob_stream_pf == 0) {
    cudaStreamCreate(&_cs_glob_stream_pf);
    cs_mem_cuda_set_prefetch_stream(_cs_glob_stream_pf);
  }

  /* Initialize NCCL if available */
#if defined(HAVE_NCCL) && defined(HAVE_MPI)
  _initialize_nccl(n_devices, cs_glob_node_n_ranks);
#endif

  /* Finally, determine whether we may use graphs for some kernel launches. */

  s = getenv("CS_CUDA_ALLOW_GRAPH");
  if (s != nullptr) {
    int i = atoi(s);
    cs_glob_cuda_allow_graph = (i <= 0) ? false : true;
  }

  return device_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return currently selected CUDA devices.
 *
 * \return  selected device id, or -1 if no usable device is available
 */
/*----------------------------------------------------------------------------*/

extern "C" int
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
