/*============================================================================
 * Definitions, global variables, and base functions for accelerators.
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

#include <assert.h>
#include <errno.h>

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <map>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

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

typedef struct
{
  void  *host_ptr;          //!< host pointer
  void  *device_ptr;        //!< host pointer

  size_t           size;    //! allocation size
  cs_alloc_mode_t  mode;    //!< allocation mode
} _cs_base_accel_mem_map;

/*============================================================================
 *  Global variables
 *============================================================================*/

static std::map<void *, _cs_base_accel_mem_map> _hd_alloc_map;

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
 * \brief Allocate memory on host and device for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
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

void *
cs_malloc_hd(cs_alloc_mode_t  mode,
             size_t           ni,
             size_t           size,
             const char      *var_name,
             const char      *file_name,
             int              line_num)
{
  _cs_base_accel_mem_map  me = {
    .host_ptr = NULL,
    .device_ptr = NULL,
    .size = ni * size,
    .mode = mode};

  if (mode <= CS_ALLOC_HOST_DEVICE)
    me.host_ptr = bft_mem_malloc(ni, size, var_name, file_name, line_num);

#if defined(HAVE_CUDA)

  // if (mode == CS_ALLOC_HOST_DEVICE_SHARED)
  //   delay host allocation to when we request device address

  if (mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    me.host_ptr = cs_cuda_mem_malloc_managed(me.size,
                                             var_name,
                                             file_name,
                                             line_num);
    me.device_ptr = me.host_ptr;
  }

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for managed allocation

#endif

  _hd_alloc_map[me.host_ptr] = me;

  /* Return pointer to allocated memory */

  return me.host_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given host pointer.
 *
 * \param [in]  host_ptr   host pointer to free
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 */
/*----------------------------------------------------------------------------*/

void
cs_free_hd(void        *host_ptr,
           const char  *var_name,
           const char  *file_name,
           int          line_num)
{
  if (_hd_alloc_map.count(host_ptr) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host pointer matching %p."), __func__, host_ptr);

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];

  if (me.mode <= CS_ALLOC_HOST_DEVICE)
    bft_mem_free(me.host_ptr, var_name, file_name, line_num);

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.mode >= CS_ALLOC_HOST_DEVICE)
      cs_cuda_mem_free(me.device_ptr, var_name, file_name, line_num);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif

  }

  _hd_alloc_map.erase(host_ptr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return device pointer for a given host pointer.
 *
 * \param [in]  host_ptr   host pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_get_device_ptr(void  *host_ptr)
{
  if (_hd_alloc_map.count(host_ptr) == 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host pointer matching %p."), __func__, host_ptr);
    return NULL;
  }

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];

  if (me.device_ptr == NULL) {
    if (me.mode == CS_ALLOC_HOST_DEVICE) {
#if defined(HAVE_CUDA)

      me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);
      _hd_alloc_map[host_ptr] = me;

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
    }
  }

  return me.device_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a host pointer is associated with a  device.
 *
 * \returns allocation mode associated with host pointer
 */
/*----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_check_device_ptr(void  *host_ptr)
{
  if (_hd_alloc_map.count(host_ptr) == 0)
    return CS_ALLOC_HOST;

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];
  return me.mode;
}

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

void *
cs_associate_device_ptr(void    *host_ptr,
                        size_t   ni,
                        size_t   size)
{
  if (_hd_alloc_map.count(host_ptr) == 0) {

    _cs_base_accel_mem_map  me = {
      .host_ptr = host_ptr,
      .device_ptr = NULL,
      .size = ni * size,
      .mode = CS_ALLOC_HOST_DEVICE};

    _hd_alloc_map[me.host_ptr] = me;

  }

  return cs_get_device_ptr(host_ptr);
}

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

void
cs_dissassociate_device_ptr(void  *host_ptr)
{
  if (_hd_alloc_map.count(host_ptr) > 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[me.host_ptr];

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.mode == CS_ALLOC_HOST_DEVICE)
      cs_cuda_mem_free(me.device_ptr, "me.device_ptr", __FILE__, __LINE__);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif

  }

  _hd_alloc_map.erase(host_ptr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of host-device allocations
 *
 * \returns current number of host-device allocations.
 */
/*----------------------------------------------------------------------------*/

int
cs_get_n_allocations_hd(void)
{
  return _hd_alloc_map.size();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
