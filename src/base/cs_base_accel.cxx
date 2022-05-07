/*============================================================================
 * Definitions, global variables, and base functions for accelerators.
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

#define _CS_ALLOC_LOG 0

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

typedef struct
{
  void  *host_ptr;          //!< host pointer
  void  *device_ptr;        //!< host pointer

  size_t           size;    //! allocation size
  cs_alloc_mode_t  mode;    //!< allocation mode

#if _CS_ALLOC_LOG == 1
  char  var_name[32];       //!< variable name
  char  src_line[64];       //!< source line info
#endif

} _cs_base_accel_mem_map;

/*============================================================================
 *  Global variables
 *============================================================================*/

static std::map<const void *, _cs_base_accel_mem_map> _hd_alloc_map;

static bool _initialized = false;

/*! Default "host+device" allocation mode */

cs_alloc_mode_t  cs_alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;

/* Keep track of active device id using OpenMP; usually queried dynamically,
   but saving the value in this variable can be useful when debugging */

int  cs_glob_omp_target_device_id = -1;

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
 * \brief Initialize memory mapping on device.
 */
/*----------------------------------------------------------------------------*/

static void
_initialize(void)
{
  bft_mem_alternative_set(cs_get_allocation_hd_size,
                          _realloc_host,
                          cs_free_hd);

  _initialized = true;

  if (cs_get_device_id() < 0)
    cs_alloc_mode = CS_ALLOC_HOST;
}

#if defined(HAVE_OPENMP_TARGET)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of OpenMP device memory.
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

static void *
_omp_target_mem_malloc_device(size_t        n,
                              const char   *var_name,
                              const char   *file_name,
                              int           line_num)
{
  void *ptr = omp_target_alloc_device(n, cs_glob_omp_target_device_id);

  if (ptr == NULL)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes on device\n"
              "  running: %s",
              (unsigned long long)n, __func__);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of host memory using OpenMP Offload.
 *
 * No OpenMP standard way to mimick cudaMallocHost today aka Host pinned memory
 * allocation + GPU driver acceleration (DMA/zero copy).
 *
 * Closest is Intel proprietary omp_target_alloc_host (accepted in OMP 6.0) or
 * new omp allocator (pinned) + explicit data transfer
 * Note: omp_target_alloc_host supports implicit data transfert.
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

static void *
_omp_target_mem_malloc_host(size_t        n,
                            const char   *var_name,
                            const char   *file_name,
                            int           line_num)
{
  void *ptr = NULL;

#if defined(__INTEL_LLVM_COMPILER)
  ptr = omp_target_alloc_host(n, cs_glob_omp_target_device_id);
#else
  assert(0 && "Not implemented yet");
#endif

  if (ptr == NULL)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes on host\n"
              "  running: %s",
              (unsigned long long)n, __func__);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of OpenMP Offload managed memory.
 *
 * Standards define pragma unified_shared_memory to drive
 * omp_target_alloc to allocate USM
 *
 * Intel proprietary omp_target_alloc_shared (accepted in OMP 6.0) is
 * another convenient way to do so.
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

static void *
_omp_target_mem_malloc_managed(size_t        n,
                               const char   *var_name,
                               const char   *file_name,
                               int           line_num)
{
#if defined(__INTEL_LLVM_COMPILER)

  void *ptr = omp_target_alloc_shared(n, cs_glob_omp_target_device_id);

#else

#pragma omp requires unified_shared_memory
  void *ptr = omp_target_alloc(n, cs_glob_omp_target_device_id);

#endif

  if (ptr == NULL)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes\n"
              "  running: %s",
              (unsigned long long)n, __func__);

  return ptr;
}

#endif /* defined(HAVE_OPENMP_TARGET) */

#if _CS_ALLOC_LOG == 1

/*----------------------------------------------------------------------------*/
/*!
 * \brief Loc variable namev and line if enabled.
 *
 * If separate pointers are used on the host and device,
 * the host pointer is returned.
 *
 * \param [in, out]  me         memory map element
 * \param [in]       var_name   allocated variable name string
 * \param [in]       file_name  name of calling source file
 * \param [in]       line_num   line number in calling source file
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

static void
_log_call_site(_cs_base_accel_mem_map  &me,
               const char              *var_name,
               const char              *file_name,
               int                      line_num)
{
  memset(me.var_name, 0, 32);
  memset(me.src_line, 0, 64);
  strncpy(me.var_name, var_name, 63);
  const char *_file_name = file_name;
  for (size_t i = 0; file_name[i] != '\0'; i++) {
    if (file_name[i] == '/')
      _file_name = file_name + i + 1;
  }
  snprintf(me.src_line, 63, "%s:%d", _file_name, line_num);
}

#endif

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

#elif defined (HAVE_OPENMP_TARGET)

  retval = omp_get_default_device();

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
  if (_initialized == false) {
   _initialize();
   mode = cs_alloc_mode;
  }

  if (ni == 0)
    return NULL;

  _cs_base_accel_mem_map  me = {
    .host_ptr = NULL,
    .device_ptr = NULL,
    .size = ni * size,
    .mode = mode};

#if _CS_ALLOC_LOG == 1
  _log_call_site(me, var_name, file_name, line_num);
#endif

  if (mode < CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = bft_mem_malloc(ni, size, var_name, file_name, line_num);

  // Device allocation will be postponed later thru call to
  // cs_get_device_ptr. This applies for CS_ALLOC_HOST_DEVICE
  // and CS_ALLOC_HOST_DEVICE_PINNED modes

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

#elif defined(HAVE_OPENMP_TARGET)

  else if (mode == CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = _omp_target_mem_malloc_host(me.size,
                                              var_name,
                                              file_name,
                                              line_num);

  else if (mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    me.host_ptr = _omp_target_mem_malloc_managed(me.size,
                                                 var_name,
                                                 file_name,
                                                 line_num);
    me.device_ptr = me.host_ptr;
  }

  else if (mode == CS_ALLOC_DEVICE)
    me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                  var_name,
                                                  file_name,
                                                  line_num);

#endif

  if (me.host_ptr != NULL)
    _hd_alloc_map[me.host_ptr] = me;
  else if (me.device_ptr != NULL)
    _hd_alloc_map[me.device_ptr] = me;

  /* Return pointer to allocated memory */

  if (me.host_ptr != NULL)
    return me.host_ptr;
  else
    return me.device_ptr;
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
 * If the allocation parameters are unchanged, no actual reallocation
 * occurs on the host.
 *
 * If the device uses a separate allocation, it is freed, and a new
 * allocation is delayed (as per initial allocation) so as to invalidate copies
 * which will not be up to date anymore after the associated values
 * modification.
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

    if (me.device_ptr != me.host_ptr && me.device_ptr != NULL) {
#if defined(HAVE_CUDA)
      cs_cuda_mem_free(me.device_ptr, var_name, file_name, line_num);
      me.device_ptr = NULL;
#elif defined(HAVE_OPENMP_TARGET)
      omp_target_free(me.device_ptr, cs_glob_omp_target_device_id);
      me.device_ptr = NULL;
#endif
    }
  }

#if _CS_ALLOC_LOG == 1
  _log_call_site(me, var_name, file_name, line_num);
#endif

  if (new_size == me.size && mode == me.mode) {
    if (me.host_ptr != NULL)
      return me.host_ptr;
    else
      return me.device_ptr;
  }

  if (   me.mode <= CS_ALLOC_HOST_DEVICE
      && me.mode == mode) {
    me.host_ptr = bft_mem_realloc(me.host_ptr, ni, size,
                                  var_name, file_name, line_num);
    me.size = new_size;
    _hd_alloc_map.erase(ptr);
    _hd_alloc_map[me.host_ptr] = me;

    ret_ptr = me.host_ptr;
  }
  else {
    ret_ptr = cs_malloc_hd(mode, 1, me.size,
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

  if (me.mode < CS_ALLOC_HOST_DEVICE_PINNED) {
    bft_mem_free(me.host_ptr, var_name, file_name, line_num);
    me.host_ptr = NULL;
  }

  if (me.host_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED) {
      cs_cuda_mem_free(me.host_ptr, var_name, file_name, line_num);
      me.device_ptr = NULL;
    }

    else
      cs_cuda_mem_free_host(me.host_ptr, var_name, file_name, line_num);

    me.host_ptr = NULL;

#elif defined(HAVE_OPENMP_TARGET)

    omp_target_free(me.host_ptr, cs_glob_omp_target_device_id);
    if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED)
      me.device_ptr = NULL;

    me.host_ptr = NULL;

#endif

  }

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    cs_cuda_mem_free(me.device_ptr, var_name, file_name, line_num);
    me.device_ptr = NULL;

#elif defined(HAVE_OPENMP_TARGET)

    omp_target_free(me.device_ptr, cs_glob_omp_target_device_id);
    me.device_ptr = NULL;

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

  /* Allocate on device if not done yet */

  if (me.device_ptr == NULL) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {
#if defined(HAVE_CUDA)

      me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

      if (omp_target_associate_ptr(me.host_ptr, me.device_ptr, me.size, 0,
                                   cs_glob_omp_target_device_id))
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: Can't associate host pointer %p to device pointer %p."),
                  "omp_target_associate_ptr", me.host_ptr, me.device_ptr);

#endif

      _hd_alloc_map[ptr] = me;

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
 * If memory is not allocated on device yet at the call site, it will
 * be allocated automatically by this function.
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

  /* Allocate and sync on device if not done yet */

  if (me.device_ptr == NULL) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {
#if defined(HAVE_CUDA)

      me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

#endif

      _hd_alloc_map[ptr] = me;
      cs_sync_h2d(ptr);

    }
  }

  return me.device_ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return matching device pointer for a given constant pointer,
 *        prefetching if applicable.
 *
 * If separate pointers are used on the host and device, the host pointer
 * should be used with this function. In this case, it is assumed that
 * the host and device values have already been synchronized, unless
 * memory is not allocated on device yet at the call site, in which case
 * it will be allocated automatically by this function.
 *
 * \param [in]  ptr  pointer
 *
 * \returns pointer to device memory.
 */
/*----------------------------------------------------------------------------*/

const void *
cs_get_device_ptr_const_pf(const void  *ptr)
{
  if (ptr == NULL)
    return NULL;

  if (_hd_alloc_map.count(ptr) == 0) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: No host or device pointer matching %p."), __func__, ptr);
    return NULL;
  }

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

  /* Allocate and sync on device if not done yet */

  if (me.device_ptr == NULL) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {
#if defined(HAVE_CUDA)

      me.device_ptr = cs_cuda_mem_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

#endif

      _hd_alloc_map[ptr] = me;
      cs_sync_h2d(ptr);

    }
  }

  /* Prefetch if shared */

  else if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED)
    cs_prefetch_h2d(me.host_ptr, me.size);

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
cs_disassociate_device_ptr(void  *host_ptr)
{
  if (_hd_alloc_map.count(host_ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[host_ptr];

  if (me.device_ptr != NULL) {

#if defined(HAVE_CUDA)

    if (me.mode == CS_ALLOC_HOST_DEVICE)
      cs_cuda_mem_free(me.device_ptr, "me.device_ptr", __FILE__, __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

    if (me.mode == CS_ALLOC_HOST_DEVICE)
      omp_target_free(me.device_ptr, cs_glob_omp_target_device_id);

#endif

  }

  _hd_alloc_map.erase(host_ptr);
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
      cs_disassociate_device_ptr(_host_ptr);

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
cs_sync_h2d(const void  *ptr)
{
  if (_hd_alloc_map.count(ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

  if (me.device_ptr == NULL)
    me.device_ptr = const_cast<void *>(cs_get_device_ptr_const(ptr));

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

#elif defined(HAVE_OPENMP_TARGET)

  switch (me.mode) {

  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;
  case CS_ALLOC_HOST_DEVICE:
    omp_target_memcpy(me.device_ptr, me.host_ptr, me.size, 0, 0,
                      cs_glob_omp_target_device_id, omp_get_initial_device());
    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    {
    char *host_ptr = (char *)me.device_ptr;
    #pragma omp target enter data map(to:host_ptr[:me.size]) \
      nowait device(cs_glob_omp_target_device_id)
    }
    break;
  case CS_ALLOC_HOST_DEVICE_SHARED:
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target enter data map(to:host_ptr[:me.size]) \
        nowait device(cs_glob_omp_target_device_id)
    }
    break;
  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initiate synchronization of data from host to device for
 *        future access.
 *
 * If separate pointers are used on the host and device,
 * the host pointer should be used with this function.
 * In this case, synchronization is started (asynchronously
 * if the allocation mode supports it).
 *
 * In other cases, synchronization will be delayed until actual use.
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
cs_sync_h2d_future(const void  *ptr)
{
  if (_hd_alloc_map.count(ptr) == 0)
    return;

  _cs_base_accel_mem_map  me = _hd_alloc_map[ptr];

#if defined(HAVE_CUDA)

  switch (me.mode) {
  case CS_ALLOC_HOST_DEVICE:
    cs_cuda_copy_h2d(me.device_ptr, me.host_ptr, me.size);
    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    cs_cuda_copy_h2d_async(me.device_ptr, me.host_ptr, me.size);
    break;
  default:
    break;
  }

#elif defined(HAVE_OPENMP_TARGET)

  switch (me.mode) {
  case CS_ALLOC_HOST_DEVICE:
    omp_target_memcpy(me.device_ptr, me.host_ptr, me.size, 0, 0,
                      cs_glob_omp_target_device_id, omp_get_initial_device());
    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    {
      char *host_ptr = (char *)me.device_ptr;
      #pragma omp target enter data map(to:host_ptr[:me.size]) \
        nowait device(cs_glob_omp_target_device_id)
    }
    break;
  default:
    break;
  }

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

#elif defined(HAVE_OPENMP_TARGET)

  switch (me.mode) {
  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;
  case CS_ALLOC_HOST_DEVICE:
    omp_target_memcpy(me.host_ptr, me.device_ptr, me.size, 0, 0,
                      omp_get_initial_device(), cs_glob_omp_target_device_id);

    break;
  case CS_ALLOC_HOST_DEVICE_PINNED:
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target exit data map(from:host_ptr[:me.size]) \
        nowait device(cs_glob_omp_target_device_id)
    }
    break;
  case CS_ALLOC_HOST_DEVICE_SHARED:
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target exit data map(from:host_ptr[:me.size]) \
        nowait device(cs_glob_omp_target_device_id)
    }
    break;
  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;
  }

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

#elif defined(HAVE_OPENMP_TARGET)

  char *host_ptr = (char *)ptr;
  #pragma omp target enter data map(to:host_ptr[:size]) \
    nowait device(cs_glob_omp_target_device_id)

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

#elif defined(HAVE_OPENMP_TARGET)

  char *host_ptr = (char *)ptr;
  #pragma omp target exit data map(from:host_ptr[:size]) \
    nowait device(cs_glob_omp_target_device_id)

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

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    cs_glob_omp_target_device_id, omp_get_initial_device());

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

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    omp_get_initial_device(), cs_glob_omp_target_device_id);

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

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    cs_glob_omp_target_device_id, cs_glob_omp_target_device_id);

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

  if (cs_glob_rank_id > -1 && n_devices > 1) {

    device_id = cs_glob_node_rank_id*n_devices / cs_glob_node_n_ranks;

    assert(device_id > -1 && device_id < n_devices);

  }

  omp_set_default_device(device_id);

  cs_glob_omp_target_device_id = device_id;

  return device_id;
}

#endif /* defined(HAVE_OPENMP_TARGET) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
