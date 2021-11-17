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
#include <string.h>

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

static std::map<const void *, _cs_base_accel_mem_map> _hd_alloc_map;

static bool _initialized = false;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reallocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate reallocation function based on
 * the requested mode, and allows introspection of the allocated memory.
 *
 * \param [in]  host_ptr   host pointer
 * \param [in]  ni         number of elements
 * \param [in]  size       element size
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

static void *
_realloc_host(void            *host_ptr,
              size_t           ni,
              size_t           size,
              const char      *var_name,
              const char      *file_name,
              int              line_num)
{
  return cs_realloc_hd(host_ptr,
                       CS_ALLOC_HOST,
                       ni, size,
                       var_name, file_name, line_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate allocation function based on
 * the requested mode, and allows introspection of the allocated memory.
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

static void
_initialize(void)
{
  bft_mem_alternative_set(cs_get_allocation_hd_size,
                          _realloc_host,
                          cs_free_hd);

  _initialized = true;
}

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

#endif

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate allocation function based on
 * the requested mode, and allows introspection of the allocated memory.
 *
 * If separate pointers are used on the host and device,
 * the host pointer is returned.
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
cs_malloc_hd(cs_alloc_mode_t   mode,
             size_t            ni,
             size_t            size,
             const char       *var_name,
             const char       *file_name,
             int               line_num)
{
  if (_initialized == false)
    _initialize();

  _cs_base_accel_mem_map  me = {
    .host_ptr = NULL,
    .device_ptr = NULL,
    .size = ni * size,
    .mode = mode};

  if (mode < CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = bft_mem_malloc(ni, size, var_name, file_name, line_num);

#if defined(HAVE_CUDA)

  else if (mode == CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = cs_cuda_mem_malloc_host(me.size,
                                          var_name,
                                          file_name,
                                          line_num);

  else if (mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    me.host_ptr = cs_cuda_mem_malloc_managed(me.size,
                                             var_name,
                                             file_name,
                                             line_num);
    me.device_ptr = me.host_ptr;
  }

  else if (mode == CS_ALLOC_DEVICE)
    me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                              var_name,
                                              file_name,
                                              line_num);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for managed allocation

#endif

  if (me.host_ptr != NULL)
    _hd_alloc_map[me.host_ptr] = me;
  else if (me.device_ptr != NULL)
    _hd_alloc_map[me.device_ptr] = me;

  /* Return pointer to allocated memory */

  return me.host_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reallocate memory on host and device for ni elements of size bytes.
 *
 * This function calls the appropriate reallocation function based on
 * the requested mode, and allows introspection of the allocated memory.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \param [in]  ptr        pointer to previously allocated memory
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
cs_realloc_hd(void            *ptr,
              cs_alloc_mode_t  mode,
              size_t           ni,
              size_t           size,
              const char      *var_name,
              const char      *file_name,
              int              line_num)
{
  void *ret_ptr = ptr;
  size_t new_size = ni*size;

  if (ptr == NULL) {
    return cs_malloc_hd(mode, ni, size, var_name, file_name, line_num);
  }
  else if (new_size == 0) {
    cs_free_hd(ptr, var_name, file_name, line_num);
    return NULL;
  }

  _cs_base_accel_mem_map  me;

  if (_hd_alloc_map.count(ptr) == 0) {  /* Case where memory was allocated
                                           on host only (through BFT_MALLOC) */
    me = {.host_ptr = ptr,
          .device_ptr = NULL,
          .size = bft_mem_get_block_size(ptr),
          .mode = CS_ALLOC_HOST};
    _hd_alloc_map[me.host_ptr] = me;
  }
  else {
    me = _hd_alloc_map[ptr];
  }

  if (new_size == me.size) {
    if (me.host_ptr != NULL)
      return me.host_ptr;
    else
      return me.device_ptr;
  }

  if (   me.mode <= CS_ALLOC_HOST
      && me.mode == mode) {
    me.host_ptr = bft_mem_realloc(me.host_ptr, ni, size,
                                  var_name, file_name, line_num);
    me.size = new_size;
    _hd_alloc_map.erase(ptr);
    _hd_alloc_map[me.host_ptr] = me;
  }
  else {
    ret_ptr = cs_malloc_hd(me.mode, 1, me.size,
                           var_name, file_name, line_num);

    memcpy(ret_ptr, ptr, me.size);

    cs_free_hd(ptr, var_name, file_name, line_num);
  }

  return ret_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given host pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \param [in]  ptr        pointer to free
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 */
/*----------------------------------------------------------------------------*/

void
cs_free_hd(void        *ptr,
           const char  *var_name,
           const char  *file_name,
           int          line_num)
{
  if (ptr == NULL)
    return;

  if (_hd_alloc_map.count(ptr) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host or device pointer matching %p."),
              __func__, ptr);

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

  if (me.mode < CS_ALLOC_HOST_DEVICE_PINNED)
    bft_mem_free(me.host_ptr, var_name, file_name, line_num);

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.host_ptr != NULL)
      cs_cuda_mem_free(me.host_ptr, var_name, file_name, line_num);

    if (me.mode != CS_ALLOC_HOST_DEVICE_SHARED)
      cs_cuda_mem_free(me.device_ptr, var_name, file_name, line_num);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif

  }

  _hd_alloc_map.erase(ptr);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given pointer.
 *
 * Compared to \cs_free_hd, this function also allows freeing memory
 * allocated through BFT_MEM_MALLOC / bft_mem_malloc.
 *
 * \param [in]  ptr        pointer to free
 * \param [in]  var_name   allocated variable name string
 * \param [in]  file_name  name of calling source file
 * \param [in]  line_num   line number in calling source file
 */
/*----------------------------------------------------------------------------*/

void
cs_free(void        *ptr,
        const char  *var_name,
        const char  *file_name,
        int          line_num)
{
  if (ptr == NULL)
    return;

  else if (_hd_alloc_map.count(ptr) == 0) {
    bft_mem_free(ptr, var_name, file_name, line_num);
    return;
  }
  else
    cs_free_hd(ptr, var_name, file_name, line_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * If memory is not allocated on device yet at the call site, it will
 * be allocated automatically by this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_get_device_ptr(void  *ptr)
{
  if (ptr == NULL)
    return NULL;

  if (_hd_alloc_map.count(ptr) == 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host or device pointer matching %p."), __func__, ptr);
    return NULL;
  }

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

  if (me.device_ptr == NULL) {
    if (me.mode == CS_ALLOC_HOST_DEVICE) {
#if defined(HAVE_CUDA)

      me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);
      _hd_alloc_map[ptr] = me;

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
    }
  }

  return me.device_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given constant pointer.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

const void *
cs_get_device_ptr_const(const void  *ptr)
{
  if (ptr == NULL)
    return NULL;

  if (_hd_alloc_map.count(ptr) == 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host or device pointer matching %p."), __func__, ptr);
    return NULL;
  }

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

  return me.device_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if a pointer is associated with a device.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * \returns allocation mode associated with pointer
 */
/*----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_check_device_ptr(const void  *ptr)
{
  if (_hd_alloc_map.count(ptr) == 0)
    return CS_ALLOC_HOST;

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];
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
  if (_hd_alloc_map.count(host_ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.mode == CS_ALLOC_HOST_DEVICE)
      cs_cuda_mem_free(me.device_ptr, "me.device_ptr", __FILE__, __LINE__);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set allocation mode for an already allocated pointer.
 *
 * If the allocation mode is different from the previous one,
 * the associated memory will be reallocated with the desired mode,
 * and the previous allocation freed.
 *
 * \param [in, out]  host_ptr   pointer to host pointer to modify
 * \param [in]       mode       desired allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_set_alloc_mode(void             **host_ptr,
                  cs_alloc_mode_t    mode)
{
  if (host_ptr == NULL)
    return;

  void *ret_ptr = *host_ptr;

  void *_host_ptr = *host_ptr;

  if (_host_ptr == NULL)
    return;

  if (_hd_alloc_map.count(_host_ptr) == 0) {

    _cs_base_accel_mem_map  me = {
      .host_ptr = _host_ptr,
      .device_ptr = NULL,
      .size = bft_mem_get_block_size(_host_ptr),
      .mode = CS_ALLOC_HOST};

    _hd_alloc_map[me.host_ptr] = me;

  }

  cs_alloc_mode_t old_mode = cs_check_device_ptr(_host_ptr);

  if (mode != old_mode) {

    _cs_base_accel_mem_map  me = _hd_alloc_map[_host_ptr];

    if (old_mode == CS_ALLOC_HOST_DEVICE)
      cs_dissassociate_device_ptr(_host_ptr);

    if (   mode == CS_ALLOC_HOST_DEVICE_SHARED
        || old_mode == CS_ALLOC_HOST_DEVICE_SHARED) {

      ret_ptr = cs_malloc_hd(mode, 1, me.size,
                             "me.host_ptr", __FILE__, __LINE__);

      /* TODO: check if we have multiple OpenMP threads, in which
         case applying a "first-touch" policy might be useful here */

      memcpy(ret_ptr, _host_ptr, me.size);

      cs_free_hd(_host_ptr, "me.host_ptr", __FILE__, __LINE__);

    }

  }

  *host_ptr = ret_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize data from host to device.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 *
 * Depending on the allocation type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * This function assumes the provided pointer was allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  host pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_sync_h2d(void  *ptr)
{
  if (_hd_alloc_map.count(ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

#if defined(HAVE_CUDA)

  switch (me.mode) {
  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;
  case CS_ALLOC_HOST_DEVICE:
    cs_cuda_copy_h2d(me.device_ptr, me.host_ptr, me.size);
    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    cs_cuda_copy_h2d_async(me.device_ptr, me.host_ptr, me.size);
    break;
  case CS_ALLOC_HOST_DEVICE_SHARED:
    cs_cuda_prefetch_h2d(me.device_ptr, me.size);
    break;
  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;
  }

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize data from device to host.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be passed to this
 * function.
 *
 * Depending on the allocaton type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * This function assumes the provided pointer was allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_sync_d2h(void  *ptr)
{
  if (_hd_alloc_map.count(ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

#if defined(HAVE_CUDA)

  switch (me.mode) {
  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;
  case CS_ALLOC_HOST_DEVICE:
    cs_cuda_copy_d2h(me.host_ptr, me.device_ptr, me.size);
    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    cs_cuda_copy_d2h_async(me.host_ptr, me.device_ptr, me.size);
    break;
  case CS_ALLOC_HOST_DEVICE_SHARED:
    cs_cuda_prefetch_d2h(me.host_ptr, me.size);
    break;
  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;
  }

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prefetch data from host to device.
 *
 * This function should only be used on arrays using shared host and device
 * memory, shuch as those allocated using CS_ALLOC_HOST_DEVICE_SHARED.
 * It should be usable on a subset of such an array.
 *
 * \param [in, out]  ptr   pointer to data to prefetch
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_prefetch_h2d(void    *ptr,
                size_t   size)
{
#if defined(HAVE_CUDA)

  cs_cuda_prefetch_h2d(ptr, size);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prefetch data from device to host.
 *
 * This function should only be used on arrays using shared host and device
 * memory, shuch as those allocated using CS_ALLOC_HOST_DEVICE_SHARED.
 * It should be usable on a subset of such an array.
 *
 * \param [in, out]  ptr   pointer to data to prefetch
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_prefetch_d2h(void    *ptr,
                size_t   size)
{
#if defined(HAVE_CUDA)

  cs_cuda_prefetch_d2h(ptr, size);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for shared allocation

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from host to device.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on device
 * \param [in, out]  src   pointer to source data on host
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_h2d(void        *dest,
            const void  *src,
            size_t       size)
{
#if defined(HAVE_CUDA)

  cs_cuda_copy_h2d(dest, src, size);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for copy

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to host.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on host
 * \param [in, out]  src   pointer to source data on device
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_d2h(void        *dest,
            const void  *src,
            size_t       size)
{
#if defined(HAVE_CUDA)

  cs_cuda_copy_d2h(dest, src, size);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for copy

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy data from device to device.
 *
 * This function should be usable on subsets of arrays allocated on the host
 * and device.
 *
 * \param [out]      dest  pointer to destination data on host
 * \param [in, out]  src   pointer to source data on device
 * \param [in]       size  number of bytes to prefetch
 */
/*----------------------------------------------------------------------------*/

void
cs_copy_d2d(void        *dest,
            const void  *src,
            size_t       size)
{
#if defined(HAVE_CUDA)

  cs_cuda_copy_d2d(dest, src, size);

#elif defined(HAVE_ONEAPI)

  // TODO add OneApi wrapper for copy

#endif
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
/*!
 * \brief Check if a given host pointer is allocated with associated with
 *        cs_alloc_hd or cs_realloc_hd.
 *
 * \returns allocated memory size, or zero if not allocated with this
 *          mechanism.
 */
/*----------------------------------------------------------------------------*/

size_t
cs_get_allocation_hd_size(void  *host_ptr)
{
  if (_hd_alloc_map.count(host_ptr) == 0)
    return 0;

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];
  return me.size;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
