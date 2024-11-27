#ifndef CS_MEM_H
#define CS_MEM_H

/*============================================================================
 * Base memory allocation wrappers with optional tracing
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

/*----------------------------------------------------------------------------*/

/* BFT headers */

#include "bft_error.h"

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/*!
 * Allocation modes for accelerated code.
 */

typedef enum {

  CS_ALLOC_HOST,                /*!< allocation on host only */
  CS_ALLOC_HOST_DEVICE,         /*!< allocation on host and device */
  CS_ALLOC_HOST_DEVICE_PINNED,  /*!< allocation on host and device,
                                  using page-locked memory on host
                                  if possible */
  CS_ALLOC_HOST_DEVICE_SHARED,  /*!< allocation on host and device,
                                  using mapped/shared memory */
  CS_ALLOC_DEVICE               /*!< allocation on device only */

} cs_alloc_mode_t;

/*============================================================================
 * Public macros
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls cs_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define CS_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) cs_mem_malloc(_ni, sizeof(_type), \
                               #_ptr, __FILE__, __LINE__)

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls cs_mem_malloc_hd(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer is returned.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 *   _mode <-- allocation mode.
 */

#define CS_MALLOC_HD(_ptr, _ni, _type, _mode) \
_ptr = (_type *) cs_mem_malloc_hd(_mode, _ni, sizeof(_type), \
                                  #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls cs_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define CS_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) cs_mem_realloc(_ptr, _ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls cs_mem_realloc_hd(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * If the allocation parameters are unchanged, no actual reallocation
 * occurs.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 *   _mode <-- allocation mode.
 */

#define CS_REALLOC_HD(_ptr, _ni, _type, _mode) \
_ptr = (_type *) cs_mem_realloc_hd(_ptr, _mode, _ni, sizeof(_type), \
                                   #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls cs_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be used with this
 * function.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#define CS_FREE(_ptr) \
cs_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

/*
 * Allocate aligned memory for _ni items of type _type.
 *
 * This macro calls cs_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr    --> pointer to allocated memory.
 *   _align <-- alignment.
 *   _ni    <-- number of items.
 *   _type  <-- element type.
 */

#define CS_MEMALIGN(_ptr, _align, _ni, _type) \
_ptr = (_type *) cs_mem_memalign(_align, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

#if defined(HAVE_ACCEL)

extern cs_alloc_mode_t  cs_alloc_mode;
extern cs_alloc_mode_t  cs_alloc_mode_read_mostly;

#else

#define cs_alloc_mode CS_ALLOC_HOST
#define cs_alloc_mode_read_mostly CS_ALLOC_HOST

#endif

/*============================================================================
 * Semi-private function prototypes
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Initialize memory handling.
 *
 * This function should be called before any other cs_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * If the log file name argument is non-null but is an empty string,
 * memory management be tracked, but not logged in detail, so only
 * statistics will be available.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 *----------------------------------------------------------------------------*/

void
cs_mem_init(const char  *log_file_name);

/*----------------------------------------------------------------------------
 * End memory handling.
 *
 * This function should be called after all other cs_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 *----------------------------------------------------------------------------*/

void
cs_mem_end(void);

/*----------------------------------------------------------------------------
 * Indicates if cs_mem_...() functions are initialized.
 *
 * returns:
 *   1 if cs_mem_init has been called, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_mem_initialized(void);

/*----------------------------------------------------------------------------
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the cs_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * Allocation couting and logging to trace file will be done if
 * both required by the cs_mem_init options and if file_name != nullptr.
 * If required but file_name == nullptr, it must be handled by the caller,
 * using `cs_mem_log_mem_op`.
 *
 * parameters:
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 *----------------------------------------------------------------------------*/

void *
cs_mem_malloc(size_t       ni,
              size_t       size,
              const char  *var_name,
              const char  *file_name,
              int          line_num);

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*
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
cs_mem_malloc_hd(cs_alloc_mode_t   mode,
                 size_t            ni,
                 size_t            size,
                 const char       *var_name,
                 const char       *file_name,
                 int               line_num);

#else

inline static void *
cs_mem_malloc_hd(cs_alloc_mode_t   mode,
                 size_t            ni,
                 size_t            size,
                 const char       *var_name,
                 const char       *file_name,
                 int               line_num)
{
  CS_UNUSED(mode);
  return cs_mem_malloc(ni, size, var_name, file_name, line_num);
}

#endif

/*----------------------------------------------------------------------------
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the cs_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, cs_alloc() called).
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num   -> line number in calling source file
 *
 * returns:
 *   pointer to allocated memory.
 *----------------------------------------------------------------------------*/

void *
cs_mem_realloc(void        *ptr,
               size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*----------------------------------------------------------------------------*/
/*
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

#if defined(HAVE_ACCEL)

void *
cs_mem_realloc_hd(void            *ptr,
                  cs_alloc_mode_t  mode,
                  size_t           ni,
                  size_t           size,
                  const char      *var_name,
                  const char      *file_name,
                  int              line_num);

#else

inline static void *
cs_mem_realloc_hd(void             *ptr,
                  cs_alloc_mode_t   mode,
                  size_t            ni,
                  size_t            size,
                  const char       *var_name,
                  const char       *file_name,
                  int               line_num)
{
  CS_UNUSED(mode);
  return cs_mem_realloc(ptr, ni, size, var_name, file_name, line_num);
}

#endif

/*----------------------------------------------------------------------------
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the cs_error() errorhandler if it fails to
 * free the corresponding memory. In case of a null pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, cs_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   null pointer.
 *----------------------------------------------------------------------------*/

void *
cs_mem_free(void        *ptr,
            const char  *var_name,
            const char  *file_name,
            int          line_num);

/*----------------------------------------------------------------------------
 * Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the cs_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function cs_mem_have_memalign() indicates if this
 * type of allocation may be used on this system.
 *
 * parameters:
 *   alignment <-- alignent.
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 *----------------------------------------------------------------------------*/

void *
cs_mem_memalign(size_t       alignment,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through cs_mem_...() (in kB).
 *----------------------------------------------------------------------------*/

size_t
cs_mem_size_current(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through cs_mem_...() (in kB).
 *----------------------------------------------------------------------------*/

size_t
cs_mem_size_max(void);

/*----------------------------------------------------------------------------
 * Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, cs_mem_memalign() will always fail.
 *
 * returns:
 *   1 if memory aligned allocation is possible, 0 otherwise.
 *----------------------------------------------------------------------------*/

int
cs_mem_have_memalign(void);

/*----------------------------------------------------------------------------*/
/* Returns the error handler associated with the cs_mem_...() functions.
 *
 * returns:
 *   pointer to the error handler function.
 *----------------------------------------------------------------------------*/

bft_error_handler_t *
cs_mem_error_handler_get(void);

/*----------------------------------------------------------------------------
 * Associates an error handler with the cs_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after cs_print_flush() is called), and the general error handler used
 * by cs_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameter:
 *   handler <-- pointer to the error handler function.
 *----------------------------------------------------------------------------*/

void
cs_mem_error_handler_set(bft_error_handler_t *handler);

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

#if defined(HAVE_OPENMP_TARGET)

/*----------------------------------------------------------------------------*/
/*
 * \brief Set OpenMp target device id
 *
 * \param [in]  device_id  device id to use for OpenMP device memory handling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_set_omp_target_device_id(int  device_id);

#endif

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return memory allocation stats, if available.
 *
 * Availability of statistics depends on the cs_mem_init options.
 *
 * \param [out]  alloc_cur   current allocation size, or nullptr
 * \param [out]  alloc_max   max allocation size, or nullptr
 * \param [out]  n_allocs    total number of allocations, or nullptr
 * \param [out]  n_reallocs  total number of reallocations, or nullptr
 * \param [out]  n_frees     total number of frees, or nullptr
 * \param [out]  n_current   total number of current allocations, or nullptr
 *
 * \return 1 if stats are available, O otherwise.
 *----------------------------------------------------------------------------*/

int
cs_mem_stats(uint64_t  *alloc_cur,
             uint64_t  *alloc_max,
             uint64_t  *n_allocs,
             uint64_t  *n_reallocs,
             uint64_t  *n_frees,
             uint64_t  *n_current);

#if defined(HAVE_ACCEL)

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
            size_t       size);

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
            size_t       size);

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
            size_t       size);

#endif // defined(HAVE_ACCEL)

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

#if defined(HAVE_ACCEL)

void *
cs_get_device_ptr(void  *ptr);

#else

inline static void *
cs_get_device_ptr(void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const(T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr_const(reinterpret_cast<void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

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

#if defined(HAVE_ACCEL)

const void *
cs_get_device_ptr_const(const void  *ptr);

#else

inline static const void *
cs_get_device_ptr_const(const void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const(const T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr_const(reinterpret_cast<const void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

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

#if defined(HAVE_ACCEL)

const void *
cs_get_device_ptr_const_pf(const void  *ptr);

#else

inline static const void *
cs_get_device_ptr_const_pf(const void  *ptr)
{
  return ptr;
}

#endif

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template <class T>
inline const T *
cs_get_device_ptr_const_pf(const T *ptr)
{
  const void *ptr_v
    = cs_get_device_ptr_const_pf(reinterpret_cast<const void *>(ptr));

  return (const T *)ptr_v;
}

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

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

#if defined(HAVE_ACCEL)

cs_alloc_mode_t
cs_check_device_ptr(const void  *ptr);

#else

inline static cs_alloc_mode_t
cs_check_device_ptr(const void  *ptr)
{
  CS_UNUSED(ptr);
  return CS_ALLOC_HOST;
}

#endif

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

#if defined(HAVE_ACCEL)

void *
cs_associate_device_ptr(void    *host_ptr,
                        size_t   ni,
                        size_t   size);

#else

#define cs_associate_device_ptr(_host_ptr, _ni, _size);

#endif

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

#if defined(HAVE_ACCEL)

void
cs_disassociate_device_ptr(void  *host_ptr);

#else

#define cs_disassociate_device_ptr(_host_ptr);

#endif

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

#if defined(HAVE_ACCEL)

void
cs_set_alloc_mode(void             **host_ptr,
                  cs_alloc_mode_t     mode);

#else

#define cs_set_alloc_mode(_host_ptr, mode);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set allocation mode for an already allocated pointer using
 *        pass by reference semantics for the pointer (C++ only).
 *
 * In C++, this function is preferred to direct use of cs_set_alloc_mode,
 * as it allows a more concise syntax and does not require additional casting.
 *
 * \param [in, out]  host_ptr   reference to host pointer to modify
 * \param [in]       mode       desired allocation mode
 */
/*----------------------------------------------------------------------------*/

END_C_DECLS

#if defined(__cplusplus) && defined(HAVE_ACCEL)

template<typename T>
static inline void
cs_set_alloc_mode_r(T*                &host_ptr,
                    cs_alloc_mode_t    mode)
{
  void *p = host_ptr;
  cs_set_alloc_mode(&p, mode);
  host_ptr = (T *)p;
}

#else

#define cs_set_alloc_mode_r(_host_ptr, mode);

#endif // __cplusplus && HAVE_ACCEL

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_mem_advise_set_read_mostly(void  *ptr);

#else

#define cs_mem_advise_set_read_mostly(ptr);

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_mem_advise_unset_read_mostly(void  *ptr);

#else

#define cs_mem_advise_unset_read_mostly(ptr);

#endif

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

#if defined(HAVE_ACCEL)

void
cs_sync_h2d(const void  *ptr);

#else

static inline void
cs_sync_h2d(const void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

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

#if defined(HAVE_ACCEL)

void
cs_sync_d2h(void  *ptr);

#else

static inline void
cs_sync_d2h(void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize data from device to host, only if needed.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be passed to this
 * function.
 *
 * Depending on the allocation type, this can imply a copy, data prefetch,
 * or a no-op.
 *
 * No operation occurs if the provided pointer was not allocated using
 * CS_MALLOC_HD or CS_REALLOC_HD, as it uses the associated mapping to
 * determine associated metadata.
 *
 * \param [in, out]  ptr  pointer to values to copy or prefetch
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_ACCEL)

void
cs_sync_d2h_if_needed(void  *ptr);

#else

static inline void
cs_sync_d2h_if_needed(void  *ptr)
{
  CS_UNUSED(ptr);
}

#endif

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

#if defined(HAVE_ACCEL)

void
cs_prefetch_h2d(void    *ptr,
                size_t   size);

#else

static inline void
cs_prefetch_h2d(void    *ptr,
                size_t   size)
{
  CS_UNUSED(ptr);
  CS_UNUSED(size);
}

#endif

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

#if defined(HAVE_ACCEL)

void
cs_prefetch_d2h(void    *ptr,
                size_t   size);

#else

static inline void
cs_prefetch_d2h(void    *ptr,
                size_t   size)
{
  CS_UNUSED(ptr);
  CS_UNUSED(size);
}

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_MEM_H */
