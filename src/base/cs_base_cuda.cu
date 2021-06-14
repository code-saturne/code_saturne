/*============================================================================
 * Low-level functions and global variables definition for CUDA.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

#define CS_CUDA_CHECK(x)                                                       \
if (cudaError_t err = (x)) {                                                   \
  bft_error(__FILE__, __LINE__, 0, _("CUDA error: %s"), cudaGetErrorString(err)); \
}

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

int  cs_glob_cuda_device_id = -1;

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
           "                       Integrated: %d\n"
           "                       Can map host memory: %d\n"
           "                       Compute mode: %s\n"),
         prop.major, prop.minor,
         mem, _("MB"),
         prop.integrated,
         prop.canMapHostMemory, mode_name);

    strncpy(buffer, prop.name, 255);
    buffer[255] = '\0';
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set CUDA device based on MPI rank and number of devices.
 *
 * \param[in]  comm            associated MPI communicator
 * \param[in]  ranks_per_node  number of ranks per node (min and max)
 */
/*----------------------------------------------------------------------------*/

void
cs_base_cuda_set_default_device(void)
{
  int device_id = 0, n_devices = 0;

  cudaError_t ret_code = cudaGetDeviceCount(&n_devices);

  if (ret_code == cudaErrorNoDevice)
    return;

  if (cudaSuccess != ret_code)
    bft_error(__FILE__, __LINE__, 0, "[CUDA errror] %d: %s\n  running: %s",
              ret_code, ::cudaGetErrorString(ret_code), __func__);

#if defined(HAVE_MPI)

  if (cs_glob_rank_id > -1 && n_devices > 1) {

    MPI_Comm  comm = cs_glob_mpi_comm;
    int       max_ranks_per_node = -1;

    /* get local rank */

    MPI_Comm sh_comm;
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &sh_comm);
    int sh_rank;
    MPI_Comm_rank(sh_comm, &sh_rank);

    MPI_Allreduce(MPI_IN_PLACE, &sh_rank, 1, MPI_INT, MPI_MAX, sh_comm);
    MPI_Comm_free(&sh_comm);

    MPI_Allreduce(&sh_rank, &max_ranks_per_node, 1, MPI_INT, MPI_MAX, comm);
    max_ranks_per_node += 1;

    int n_ranks_per_device = max_ranks_per_node / n_devices;

    device_id = sh_rank / n_ranks_per_device;

  }

#endif

  CS_CUDA_CHECK_CALL(cudaSetDevice(device_id), __FILE__, __LINE__);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
