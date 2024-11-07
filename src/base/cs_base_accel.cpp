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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <map>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "cs_mem.h"

#if defined(HAVE_CUDA)
#include "cs_base_cuda.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base_accel.h"

/*----------------------------------------------------------------------------*/

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

/*! Is MPI device-aware ? */
/*------------------------*/

#if defined(OMPI_MAJOR_VERSION)
  #include <mpi-ext.h>
#endif

#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
int cs_mpi_device_support = 1;

#elif defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
/* We need better detection here, as OMPI_HAVE_MPI_EXT_CUDA = 1
   does not seem to guarantee device support is present or active
   (based on test on workstation). So do not activate yet.*/
int cs_mpi_device_support = 0;

#else
int cs_mpi_device_support = 0;
#endif

/*! Default queue for SYCL */

#if defined(SYCL_LANGUAGE_VERSION)
static int    _cs_glob_sycl_device_id = -1;
sycl::queue   cs_glob_sycl_queue;
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return currently associated device id.
 *
 * \returns currently available device id, or -1 if none is available.
 */
/*----------------------------------------------------------------------------*/

int
cs_get_device_id(void)
{
  int retval = -1;

#if defined(HAVE_CUDA)

  retval = cs_base_cuda_get_device();

#elif defined(SYCL_LANGUAGE_VERSION)

  retval = _cs_glob_sycl_device_id;

#elif defined (HAVE_OPENMP_TARGET)

  retval = omp_get_default_device();

#endif

  return retval;
}

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
cs_sycl_select_default_device(void)
{
  int device_id = -1;

#if defined(SYCL_LANGUAGE_VERSION)

  sycl::property_list q_prop{sycl::property::queue::in_order()};

  try {
    sycl::queue q{sycl::gpu_selector_v, q_prop};
    cs_glob_sycl_queue = q;
    if (q.get_device().is_gpu()) {
      device_id = 0;
      cs_alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;
      cs_alloc_mode_read_mostly = CS_ALLOC_HOST_DEVICE_SHARED;
    }
  }
  catch (sycl::exception const& ex) {
    sycl::queue q{sycl::cpu_selector_v, q_prop};
    cs_glob_sycl_queue = q;
  }

#endif

  /* Return default device id */

  _cs_glob_sycl_device_id = device_id;

  return device_id;
}

#endif /* defined(HAVE_SYCL) */

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
cs_omp_target_select_default_device(void)
{
  int device_id = omp_get_initial_device();

  int n_devices = omp_get_num_devices();

  if (getenv("OMP_DEFAULT_DEVICE") != nullptr) {
    device_id = atoi(getenv("OMP_DEFAULT_DEVICE"));
  }
  else if (getenv("LIBOMPTARGET_DEVICETYPE") != nullptr) {
    device_id = omp_get_default_device();
  }
  else if (n_devices > 1) {

    if (cs_glob_rank_id > -1) {
      device_id = cs_glob_node_rank_id*n_devices / cs_glob_node_n_ranks;
      if (device_id >= n_devices)
        device_id = n_devices - 1;
    }

    else
      device_id = omp_get_default_device();

    assert(device_id > -1 && device_id < n_devices);

  }

  omp_set_default_device(device_id);
  cs_mem_set_omp_target_device_id(device_id);

  if (device_id >= 0) {
    cs_alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;
    cs_alloc_mode_read_mostly = CS_ALLOC_HOST_DEVICE_SHARED;
  }

  /* Also detect whether MPI is device-aware,
     when this can be set dynamically. */

#if defined(I_MPI_VERSION)

  {
    const char *p = getenv("I_MPI_OFFLOAD");
    if (p != nullptr) {
      if (atoi(p) > 0)
        cs_mpi_device_support = 1;
    }
  }

#endif

  /* Return default device id */

  return device_id;
}

#endif /* defined(HAVE_OPENMP_TARGET) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
