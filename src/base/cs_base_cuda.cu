/*============================================================================
 * Low-level functions and global variables definition for CUDA.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

END_C_DECLS
