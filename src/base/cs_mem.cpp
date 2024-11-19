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

/*
  Define _GNU_SOURCE if necessary before including any headers, to ensure
  the correct feature macros are defined first.
*/

#if defined(__linux__) && !defined(_GNU_SOURCE)
#  define _GNU_SOURCE
#endif

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <map>

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#if defined(HAVE_CUDA)
#include "cs_mem_cuda_priv.h"
#endif

#include "cs_mem.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_mem.cpp
        Base memory allocation wrappers with optional tracing.

  The memory managment function provided here provide optional logging,
  and tracking of non-freed pointers.

  Since in most of the intended applications, failure to allocate memory
  is considered fatal, failed allocations from these functions are
  considedered as errors, which are fatal by default but can be handled
  differently if an appropriate error handler is defined. So additional
  checking of the return values in the calling code is not needed.

  The functions provided here are otherwise based on the matching C library
  functions.
*/

/*-------------------------------------------------------------------------------
 * Local macro documentation
 *-----------------------------------------------------------------------------*/

/*! \fn CS_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls cs_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn CS_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls cs_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn CS_MEMALIGN(_ptr, _align, _ni, _type)
 * \brief Allocate aligned memory for _ni elements of type _type.
 *
 * This macro calls cs_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr   pointer to allocated memory.
 * \param [in]  _align alignment.
 * \param [in]  _ni    number of elements.
 * \param [in]  _type  element type.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-------------------------------------------------------------------------------
 * Local macro definitions
 *-----------------------------------------------------------------------------*/

/* Directory name separator
   (historically, '/' for Unix/Linux, '\' for Windows, ':' for Mac
   but '/' should work for all on modern systems) */

#define DIR_SEPARATOR '/'

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*
 * Structure defining an allocated memory block
 */

typedef struct
{
  void  *host_ptr;          //!< host pointer
#if defined(HAVE_ACCEL)
  void  *device_ptr;        //!< device pointer
#endif

  size_t           size;    //! allocation size
#if defined(HAVE_ACCEL)
  cs_alloc_mode_t  mode;    //!< allocation mode
#endif

} cs_mem_block_t;

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default memory error handler.
 *
 * Memory status is written to stderr, and the general error handler
 * used by bft_error() is then called (which results in the
 * termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_cs_mem_error_handler_default(const char  *file_name,
                              int          line_num,
                              int          sys_error_code,
                              const char  *format,
                              va_list      arg_ptr);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static int  _cs_mem_global_init_mode = 0;

static FILE *_cs_mem_global_file = nullptr;

static std::map<const void *, cs_mem_block_t> _cs_alloc_map;

static size_t  _cs_mem_global_alloc_cur = 0;
static size_t  _cs_mem_global_alloc_max = 0;

static size_t  _cs_mem_global_n_allocs = 0;
static size_t  _cs_mem_global_n_reallocs = 0;
static size_t  _cs_mem_global_n_frees = 0;

static bft_error_handler_t  *_cs_mem_error_handler
                              = (_cs_mem_error_handler_default);

#if defined(HAVE_OPENMP)
static omp_lock_t _cs_mem_lock;
#endif

#if defined(HAVE_ACCEL)
static bool _ignore_prefetch = false;
#endif

/*! Default queue for SYCL */
#if defined(SYCL_LANGUAGE_VERSION)
extern sycl::queue  cs_glob_sycl_queue;
#endif

#if defined(HAVE_OPENMP_TARGET)
/* Active device id */
static int  _omp_target_device_id = -1;
  #if defined(HAVE_OPENMP_TARGET_USM)
  #  pragma omp requires unified_shared_memory
  #endif
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

/*! Default "host+device" allocation mode */
/*----------------------------------------*/

#if defined(HAVE_ACCEL)

cs_alloc_mode_t  cs_alloc_mode = CS_ALLOC_HOST;
cs_alloc_mode_t  cs_alloc_mode_read_mostly = CS_ALLOC_HOST;

#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Given a character string representing a file name, returns
 * pointer to that part of the string corresponding to the base name.
 *
 * parameters:
 *   file_name:      <-- full name of source file.
 *
 * return:
 *   pointer to part of file name corresponding to base name.
 *----------------------------------------------------------------------------*/

static const char *
_cs_mem_basename(const char  *file_name)
{
  int i;

  if (file_name == nullptr)
    return nullptr;

  for (i = strlen(file_name) - 1;
       i > 0 && file_name[i] != DIR_SEPARATOR;
       i--);

  if (file_name[i] == DIR_SEPARATOR)
    i++;

  return (file_name + i);
}

/*-----------------------------------------------------------------------------
 * Determines values associated with an array representing a
 * long integer
 *
 * parameters:
 *   counter: <-- counter to update.
 *   value:   --> counter values in output unit [out].
 *   unit:    --> counter value unit : ' ', 'k', 'm', 'g', 't', or 'p'
 *                for bytes, Kilobytes, Megabytes, Gigabytes, Terabytes,
 *                or Petabytes [out].
 *----------------------------------------------------------------------------*/

static void
_cs_mem_size_val(const size_t    counter,
                 unsigned long   value[2],
                 char           *unit)
{
  int i;
  size_t _counter[2];
  const char units[] = {' ', 'k', 'm', 'g', 't', 'p', 'e'};

  for (i = 0, _counter[0] = counter, _counter[1] = 0;
       _counter[0] >= 1024 && i < 6;
       i++) {
    _counter[1] = _counter[0] % 1024;
    _counter[0] /= 1024;
  }

  value[0] = _counter[0];
  value[1] = _counter[1];
  *unit = units[i];
}

/*-----------------------------------------------------------------------------
 * Memory usage summary.
 *----------------------------------------------------------------------------*/

static void
_cs_mem_summary(FILE  *f)
{
  char unit;
  unsigned long value[2];
  size_t mem_usage;

  if (f == nullptr)
    return;

  fprintf(f, "\n\n");
  fprintf(f, "Memory allocation summary\n"
          "-------------------------\n\n");

  /* Available memory usage information */

  _cs_mem_size_val(_cs_mem_global_alloc_cur, value, &unit);
  fprintf(f,
          "Theoretical current allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  _cs_mem_size_val(_cs_mem_global_alloc_max, value, &unit);
  fprintf(f,
          "Theoretical maximum allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  fprintf(f,
          "\n"
          "Number of allocations:   %lu\n"
          "          reallocations: %lu\n"
          "          frees:         %lu\n\n",
          (unsigned long)_cs_mem_global_n_allocs,
          (unsigned long)_cs_mem_global_n_reallocs,
          (unsigned long)_cs_mem_global_n_frees);

  if (bft_mem_usage_initialized() == 1) {

    /* Maximum measured memory */

    mem_usage = bft_mem_usage_max_pr_size();
    if (mem_usage > 0) {
      fprintf(f,
              "Maximum program memory measure:  %8lu kB\n",
              (unsigned long)mem_usage);
    }

    /* Current measured memory */

    mem_usage = bft_mem_usage_pr_size();
    if (mem_usage > 0)
      fprintf(f,
              "Current program memory measure:   %8lu kB\n",
              (unsigned long)mem_usage);
  }
}

/*-----------------------------------------------------------------------------
 * Default memory error handler.
 *
 * Memory status is written to stderr (after bft_print_flush() is called),
 * and the general error handler used by bft_error() is then called
 * (which results in the termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 *----------------------------------------------------------------------------*/

static void
_cs_mem_error_handler_default(const char  *file_name,
                               int          line_num,
                               int          sys_error_code,
                               const char  *format,
                               va_list      arg_ptr)
{
  bft_error_handler_t * general_err_handler;

  bft_printf_flush();

  _cs_mem_summary(stderr);

  general_err_handler = bft_error_handler_get();
  general_err_handler(file_name, line_num, sys_error_code, format, arg_ptr);
}

/*-----------------------------------------------------------------------------
 * Calls the error handler (set by cs_mem_error_handler_set() or default).
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameters:
 *   file_name      <-- name of source file from which failed cs_mem_...()
 *                      function was called.
 *   line_num       <-- line of source file from which failed cs_mem_...()
 *                      function was called.
 *   sys_error_code <-- error code if error in system or libc call,
 *                      0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   ...            <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

static void
_cs_mem_error(const char  *file_name,
               int          line_num,
               int          sys_error_code,
               const char  *format,
               ...)
{
  va_list  arg_ptr;

  if (_cs_mem_global_file != nullptr) {
    _cs_mem_summary(_cs_mem_global_file);
    fflush(_cs_mem_global_file);
  }

  va_start(arg_ptr, format);

  _cs_mem_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*-----------------------------------------------------------------------------
 * Call error function when pointer is not found in allocation info.
 *
 * parameters:
 *   p: <-- allocated block's start adress.
 *
 * returns:
 *   corresponding _cs_mem_block structure.
 *----------------------------------------------------------------------------*/

static void
_cs_mem_block_info_error(const void  *p)
{
  _cs_mem_error(__FILE__, __LINE__, 0,
                 _("Adress [%p] does not correspond to "
                   "the beginning of an allocated block."),
                 p);
}

/*-----------------------------------------------------------------------------
 * Fill a cs_mem_block_t structure for an allocated pointer.
 *----------------------------------------------------------------------------*/

static inline cs_mem_block_t
_cs_mem_block_new(void          *p_new,
                   const size_t   size_new)
{
#if defined(HAVE_ACCEL)
  cs_mem_block_t  mib = {
    .host_ptr = p_new,
    .device_ptr = nullptr,
    .size = size_new,
    .mode = CS_ALLOC_HOST};
#else
  cs_mem_block_t  mib = {
    .host_ptr = p_new,
    .size = size_new};
#endif

  return mib;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the cs_mem_block structure corresponding to a given
 * allocated block.
 *
 * \param [in]  p_get  allocated block's start adress.
 *
 * \return  corresponding cs_mem_block structure.
 */
/*----------------------------------------------------------------------------*/

static cs_mem_block_t
_get_block_info(const void  *p_get)
{
  auto it = _cs_alloc_map.find(p_get);

  if (it == _cs_alloc_map.end())
    _cs_mem_block_info_error(p_get);

  return it->second;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the cs_mem_block structure corresponding to a given
 * allocated block if available.
 *
 * If no block info is available, return block with null pointers
 * and zero size.
 *
 * \param [in]  p_get  allocated block's start adress.
 *
 * \return  corresponding cs_mem_block structure.
 */
/*----------------------------------------------------------------------------*/

static cs_mem_block_t
_get_block_info_try(const void  *p_get)
{
  cs_mem_block_t mbi;

  auto it = _cs_alloc_map.find(p_get);
  if (it != _cs_alloc_map.end())
    mbi = it->second;
  else {
    mbi.host_ptr = nullptr;
#if defined(HAVE_ACCEL)
    mbi.device_ptr = nullptr;
#endif
    mbi.size = 0;
#if defined(HAVE_ACCEL)
    mbi.mode = CS_ALLOC_HOST;
#endif
  }

  return mbi;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log matching memory operation if logging is enabled
 *
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 * \param [in] old_block pointer to old block info, if present
 * \param [in] new_block pointer to new block info, if present
 */
/*----------------------------------------------------------------------------*/

static void
_update_block_info(const char            *var_name,
                   const char            *file_name,
                   int                    line_num,
                   const cs_mem_block_t  *old_block,
                   const cs_mem_block_t  *new_block)
{
  const void *p_m_old = nullptr, *p_m_new = nullptr;
  size_t old_size = 0, new_size = 0;
  cs_alloc_mode_t old_mode = CS_ALLOC_HOST, new_mode = CS_ALLOC_HOST;

  if (old_block != nullptr) {
    p_m_old = old_block->host_ptr;
#if defined(HAVE_ACCEL)
    old_mode = old_block->mode;
    if (old_mode == CS_ALLOC_DEVICE)
      p_m_old = old_block->device_ptr;
#endif
    old_size = old_block->size;
  }

  if (new_block != nullptr) {
    p_m_new = new_block->host_ptr;
#if defined(HAVE_ACCEL)
    new_mode = new_block->mode;
    if (new_mode == CS_ALLOC_DEVICE)
      p_m_new = new_block->device_ptr;
    assert(   new_block->device_ptr == nullptr
           || new_block->mode != CS_ALLOC_HOST);
#endif
    new_size = new_block->size;
  }

  if (   _cs_mem_global_init_mode > 1
      ||  old_mode > CS_ALLOC_HOST || new_mode > CS_ALLOC_HOST) {

#if defined(HAVE_OPENMP)
    int in_parallel = omp_in_parallel();
    if (in_parallel)
      omp_set_lock(&_cs_mem_lock);
#endif

    /* Update map */

    if (old_block != nullptr && (new_block == nullptr || p_m_new != p_m_old))
      _cs_alloc_map.erase(p_m_old);
    if (new_block != nullptr && p_m_new != nullptr)
      _cs_alloc_map[p_m_new] = *new_block;

    /* Memory allocation counting */

    if (_cs_mem_global_init_mode > 1) {
      long long size_diff = new_size - old_size;

      int cat_id = 3;
      if (p_m_old == nullptr) {
        cat_id = 0;
        _cs_mem_global_n_allocs += 1;
      }
      else if (p_m_new == nullptr) {
        cat_id = 2;
        _cs_mem_global_n_frees += 1;
      }
      else if (old_size != new_size) {
        cat_id = 1;
        _cs_mem_global_n_reallocs += 1;
      }

      _cs_mem_global_alloc_cur += size_diff;

      if (_cs_mem_global_alloc_max < _cs_mem_global_alloc_cur)
        _cs_mem_global_alloc_max = _cs_mem_global_alloc_cur;

      /* Log to file */

      if (_cs_mem_global_file != nullptr) {

        static const char cat_s[4][8]
          = {"  alloc", "realloc", "   free", "mapping"};
        char c_sgn = '+';
        if (size_diff < 0) {
          c_sgn = '-';
          size_diff = - size_diff;
        }

        const void *p_old = nullptr, *p_new = nullptr;
        if (old_block != nullptr)
          p_old = old_block->host_ptr;
        if (new_block != nullptr)
          p_new = new_block->host_ptr;

        fprintf(_cs_mem_global_file, "\n%s: %-27s:%6d : %-39s: %9lu",
                cat_s[cat_id], _cs_mem_basename(file_name), line_num,
                var_name, (unsigned long)new_size);
        fprintf(_cs_mem_global_file, " : (%c%9lu) : %12lu : [%14p] : [%14p]",
                c_sgn, (unsigned long)size_diff,
                (unsigned long)_cs_mem_global_alloc_cur,
                p_old, p_new);
#if defined(HAVE_ACCEL)
        static const char *alloc_mode_s[5] = {
          "host",
          "host/device",
          "host/device pinned",
          "host/device shared",
          "device"};
        const void *p_old_d
          = (old_block != nullptr) ? old_block->device_ptr : nullptr;
        const void *p_new_d
          = (new_block != nullptr) ? new_block->device_ptr : nullptr;
        fprintf(_cs_mem_global_file, " : [%14p] : [%14p] : %s",
                p_old_d, p_new_d, alloc_mode_s[new_mode]);
#endif
        fflush(_cs_mem_global_file);

      } /* End of log to file */

    } /* End counting and logging */

#if defined(HAVE_OPENMP)
    if (in_parallel)
      omp_unset_lock(&_cs_mem_lock);
#endif

  } /* End if map needs to be updated */
}

#if defined(SYCL_LANGUAGE_VERSION)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of SYCL device memory.
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
_sycl_mem_malloc_device(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num)
{
  void *ptr = sycl::malloc_device(n, cs_glob_sycl_queue);

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[sycl::malloc_device error]: unable to allocate %llu bytes\n"
              "  on device running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of host memory using SYCL.
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
_sycl_mem_malloc_host(size_t        n,
                      const char   *var_name,
                      const char   *file_name,
                      int           line_num)
{
  void *ptr = sycl::malloc_host(n, cs_glob_sycl_queue);

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[sycl::malloc_host error]: unable to allocate %llu bytes\n"
              "  on device running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

  return ptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate n bytes of SYCL USM shared memory.
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
_sycl_mem_malloc_shared(size_t        n,
                        const char   *var_name,
                        const char   *file_name,
                        int           line_num)
{
  void *ptr = sycl::malloc_shared(n, cs_glob_sycl_queue);

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[sycl::malloc_shared error]: unable to allocate %llu bytes\n"
              "  on device running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

  return ptr;
}

#elif defined(HAVE_OPENMP_TARGET)

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
#if defined(__INTEL_LLVM_COMPILER)
  void *ptr = omp_target_alloc_device(n, _omp_target_device_id);
#else
  void *ptr = omp_target_alloc(n, _omp_target_device_id);
#endif

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes on device\n"
              "  running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

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
  void *ptr = nullptr;

#if defined(__INTEL_LLVM_COMPILER)
  ptr = omp_target_alloc_host(n, _omp_target_device_id);
#else
  bft_error(file_name, line_num, 0,
            "[OpenMP target error]: unified shared memory not supported\n"
            "  running %s for variable %s.",
            __func__, var_name);
#endif

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes on host\n"
              "  running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

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

  void *ptr = omp_target_alloc_shared(n, _omp_target_device_id);

#elif defined(HAVE_OPENMP_TARGET_USM)

  void *ptr = omp_target_alloc(n, _omp_target_device_id);

#else

  void *ptr = nullptr;
  bft_error(file_name, line_num, 0,
            "[OpenMP target error]: unified shared memory not supported\n"
            "  running %s for variable %s.",
            __func__, var_name);

#endif

  if (ptr == nullptr)
    bft_error(file_name, line_num, 0,
              "[OpenMP offload error]: unable to allocate %llu bytes\n"
              "  running %s for variable %s.",
              (unsigned long long)n, __func__, var_name);

  return ptr;
}

#endif /* defined(HAVE_OPENMP_TARGET) */

#if defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free memory on host and device for a given host pointer
 *        for allocations with mode > CS_ALLOC_HOST
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

static void
_free_hd(cs_mem_block_t  &me,
         const char      *var_name,
         const char      *file_name,
         int              line_num)
{
  if (me.host_ptr != nullptr) {

#if defined(HAVE_CUDA)

    if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED)
      cs_mem_cuda_free(me.host_ptr, var_name, file_name, line_num);
    else
      cs_mem_cuda_free_host(me.host_ptr, var_name, file_name, line_num);

#elif defined(SYCL_LANGUAGE_VERSION)

    sycl::free(me.host_ptr, cs_glob_sycl_queue);

#elif defined(HAVE_OPENMP_TARGET)

    omp_target_free(me.host_ptr, _omp_target_device_id);

#endif

  }

  if (me.device_ptr != nullptr && me.device_ptr != me.host_ptr) {

#if defined(HAVE_CUDA)

    cs_mem_cuda_free(me.device_ptr, var_name, file_name, line_num);

#elif defined(SYCL_LANGUAGE_VERSION)

    sycl::free(me.device_ptr, cs_glob_sycl_queue);

#elif defined(HAVE_OPENMP_TARGET)

    omp_target_free(me.device_ptr, _omp_target_device_id);

#endif

  }
}

#endif // defined(HAVE_ACCEL)

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

#if defined(HAVE_OPENMP_TARGET)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set OpenMp target device id
 *
 * \param [in]  device_id  device id to use for OpenMP device memory handling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_set_omp_target_device_id(int  device_id)
{
  _omp_target_device_id = device_id;
}

#endif

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize memory handling.
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
 * \param log_file_name name of optional log_file (if nullptr, no log).
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_init(const char *log_file_name)
{
#if defined(HAVE_OPENMP)
  if (omp_in_parallel()) {
    if (omp_get_thread_num() != 0)
      return;
  }
  omp_init_lock(&_cs_mem_lock);
#endif

  if (_cs_mem_global_init_mode > 0) {
    _cs_mem_error(__FILE__, __LINE__, 0,
                  _("cs_mem_init() has already been called"));
  }

#if defined(HAVE_ACCEL)
  _cs_mem_global_init_mode = 2;
  const char s[] = "CS_HD_IGNORE_PREFETCH";
  if (getenv(s) != nullptr) {
    int i = atoi(getenv(s));
    if (i > 0)
      _ignore_prefetch = true;
  }
#else
  _cs_mem_global_init_mode = 1;
#endif

  if (log_file_name != nullptr) {

    _cs_mem_global_init_mode = 2;

    if (strlen(log_file_name) > 0) {

      _cs_mem_global_file = fopen(log_file_name, "w");

      /*
        If the file could not be opened, we do not abort, as it is not
        absolutely necessary. We silently continue.
        (We could warn the user, but this would require either using
        bft_printf(), which we prefer to keep independent of the cs_mem_...()
        functions to avoid evental crossed definitions when user-defined, or
        "warning handling" similar to error handling, with a possibility
        of user-defined warning handlers, as we are not sure if the calling
        code uses stderr (especially in a distributed environment). This
        is probably not worth the bother.
      */

      if (_cs_mem_global_file == nullptr)
        fprintf(stderr,
                _("Failure to open memory log file \"%s\"\n"),
                log_file_name);

    }

  }

  /* Log file header */

  if (_cs_mem_global_file != nullptr) {

    fprintf(_cs_mem_global_file,
            "       :     FILE NAME              : LINE  :"
            "  POINTER NAME                          : N BYTES   :"
            " (+- N BYTES) : TOTAL BYTES  : [   OLD ADRESS ] : [   NEW ADRESS ]");
#if defined(HAVE_ACCEL)
    fprintf(_cs_mem_global_file,
            " : [OLD ADRESS (D)] : [NEW ADRESS (D)] : MODE");
#endif
    fprintf(_cs_mem_global_file, "\n");

    fprintf(_cs_mem_global_file,
            "-------:----------------------------:-------:"
            "----------------------------------------:-----------:"
            "-----------------------------:------------------:-----------------");
#if defined(HAVE_ACCEL)
    fprintf(_cs_mem_global_file,
            "-:------------------:------------------:-----");
#endif
    fprintf(_cs_mem_global_file, "\n");

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief End memory handling.
 *
 * This function should be called after all other cs_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_end(void)
{
  if (_cs_mem_global_init_mode == 0)
    return;

#if defined(HAVE_OPENMP)
  if (omp_in_parallel()) {
    if (omp_get_thread_num() != 0)
      return;
  }
  omp_destroy_lock(&_cs_mem_lock);
#endif

  _cs_mem_global_init_mode = 0;

  if (_cs_mem_global_file != nullptr) {

    /* Memory usage summary */

    _cs_mem_summary(_cs_mem_global_file);

    /* List of non-freed pointers */

    size_t non_free = _cs_alloc_map.size();

    if (non_free > 0) {

      fprintf(_cs_mem_global_file, "List of non freed pointers:\n");

      for (auto const& x : _cs_alloc_map) {
        const void *p = x.first;
        // const cs_mem_block_t b = x.second;

        fprintf(_cs_mem_global_file,"[%p]\n", p);
      }

    }

    fprintf(_cs_mem_global_file,
            "Number of non freed pointers remaining: %lu\n",
            non_free);

    fclose(_cs_mem_global_file);
  }

  /* Reset defaults in case of later initialization */

  _cs_alloc_map.clear();

  _cs_mem_global_n_allocs = 0;
  _cs_mem_global_n_reallocs = 0;
  _cs_mem_global_n_frees = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicates if cs_mem_...() functions are initialized.
 *
 * \returns 1 if cs_mem_init has been called, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_mem_initialized(void)
{
  return (_cs_mem_global_init_mode > 0) ? 1 : 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * Allocation couting and logging to trace file will be done if
 * both required by the cs_mem_init options and if file_name != nullptr.
 * If required but file_name == nullptr, it must be handled by the caller,
 * using \ref cs_mem_log_mem_op.
 *
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_mem_malloc(size_t       ni,
              size_t       size,
              const char  *var_name,
              const char  *file_name,
              int          line_num)
{
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return nullptr;

  /* Allocate memory and check return */

  void  *p_new = malloc(alloc_size);

  if (p_new == nullptr) {
    _cs_mem_error(file_name, line_num, errno,
                  _("Failure to allocate \"%s\" (%lu bytes)"),
                  var_name, (unsigned long)alloc_size);
    return nullptr;
  }
  else if (_cs_mem_global_init_mode < 2)
    return p_new;

  cs_mem_block_t mib = _cs_mem_block_new(p_new, alloc_size);

  if (file_name != nullptr)
    _update_block_info(var_name,
                       file_name,
                       line_num,
                       nullptr,
                       &mib);

  /* Return pointer to allocated memory */

  return p_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if nullptr, cs_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_mem_realloc(void        *ptr,
               size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num)
{
  size_t new_size = ni * size;

  /*
    Behave as cs_mem_malloc() if the previous pointer is equal to nullptr.
    Note that the operation will then appear as a first allocation
    ('alloc') in the _cs_mem_global_file trace file.
  */

  if (ptr == nullptr)
    return cs_mem_malloc(ni,
                         size,
                         var_name,
                         file_name,
                         line_num);

  /*
    Behave as cs_mem_free() if the requested size is zero.
    Note that in this case, the operation will appear as 'free'
    in the _cs_mem_global_file trace file.
  */

  else if (ni == 0)
    return cs_mem_free(ptr,
                       var_name,
                       file_name,
                       line_num);

  /* When possible, get previous allocation information. */

  cs_mem_block_t mib_old;
  if (_cs_mem_global_init_mode > 1)
    mib_old = _get_block_info(ptr);
  else
    mib_old = _get_block_info_try(ptr);

  /* If the old size is known to equal the new size,
     nothing needs to be done. */

  if (new_size == mib_old.size) {
    return ptr;
  }

  /* In the general case, we have a true reallocation. */

#if defined(HAVE_ACCEL)
  if (mib_old.mode >= CS_ALLOC_HOST_DEVICE_PINNED) {
    return cs_mem_realloc_hd(ptr, CS_ALLOC_HOST, ni, size,
                             var_name, file_name, line_num);
  }
#endif

  void *p_new = realloc(ptr, new_size);

  if (file_name != nullptr) {
    cs_mem_block_t mib_new = _cs_mem_block_new(p_new, new_size);

    _update_block_info(var_name,
                       file_name,
                       line_num,
                       &mib_old,
                       &mib_new);
  }

  return p_new;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * free the corresponding memory. In case of a null pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if nullptr, cs_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns null pointer.
 */
/*----------------------------------------------------------------------------*/

void *
cs_mem_free(void        *ptr,
            const char  *var_name,
            const char  *file_name,
            int          line_num)
{
  /* null pointer case (non-allocated location) */

  if (ptr == nullptr)
    return nullptr;

  /* General case (free allocated memory) */

  /* When possible, get previous allocation information. */

  cs_mem_block_t mib_old;
  if (_cs_mem_global_init_mode > 1)
    mib_old = _get_block_info(ptr);
  else
    mib_old = _get_block_info_try(ptr);

#if defined(HAVE_ACCEL)
  if (mib_old.mode < CS_ALLOC_HOST_DEVICE_PINNED)
    free(ptr);
  else {
    _free_hd(mib_old, var_name, file_name, line_num);
  }
#else
  free(ptr);
#endif

  if (file_name != nullptr)
    _update_block_info(var_name,
                       file_name,
                       line_num,
                       &mib_old,
                       nullptr);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the bft_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function cs_mem_have_memalign() indicates if this
 * type of allocation may be used on this system.
 *
 * \param [in]  alignment alignment.
 * \param [in]  ni        number of elements.
 * \param [in]  size      element size.
 * \param [in]  var_name  allocated variable name string.
 * \param [in]  file_name name of calling source file.
 * \param [in]  line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */
/*----------------------------------------------------------------------------*/

void *
cs_mem_memalign(size_t       alignment,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
#if defined(HAVE_POSIX_MEMALIGN)

  int         retval;
  void       *p_loc;
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return nullptr;

  /* Allocate memory and check return */

  retval = posix_memalign(&p_loc, alignment, alloc_size);

  if (retval != 0) {
    switch (retval) {
    case EINVAL:
      _cs_mem_error(file_name, line_num, 0,
                    _("Alignment %lu for \"%s\" not a power of 2\n"
                      "or a multiple of sizeof(void *) = %lu"),
                    (unsigned long)alignment, var_name,
                    (unsigned long)(sizeof(void *)));
      break;
    default:
      _cs_mem_error(file_name, line_num, 0,
                    _("Failure to allocate \"%s\" (%lu bytes)"),
                    var_name, (unsigned long)alloc_size);
    }
    return nullptr;
  }
  else if (_cs_mem_global_init_mode < 2)
    return p_loc;

  cs_mem_block_t mib = _cs_mem_block_new(p_loc, alloc_size);

  if (file_name != nullptr)
    _update_block_info(var_name,
                       file_name,
                       line_num,
                       nullptr,
                       &mib);

  /* Return pointer to allocated memory */

  return p_loc;

#else

  _cs_mem_error(file_name, line_num, errno,
                _("No aligned allocation function available on this system"));

  return nullptr;

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through cs_mem_...() (in kB).
 */
/*----------------------------------------------------------------------------*/

size_t
cs_mem_size_current(void)
{
  return (_cs_mem_global_alloc_cur / 1024);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through cs_mem_...() (in kB).
 */
/*----------------------------------------------------------------------------*/

size_t
cs_mem_size_max(void)
{
  return (_cs_mem_global_alloc_max / 1024);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return current theoretical dynamic memory allocated.
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
 */
/*----------------------------------------------------------------------------*/

int
cs_mem_stats(uint64_t  *alloc_cur,
             uint64_t  *alloc_max,
             uint64_t  *n_allocs,
             uint64_t  *n_reallocs,
             uint64_t  *n_frees,
             uint64_t  *n_current)
{
  int retval = 0;

  if (_cs_mem_global_init_mode > 1) {
    if (alloc_cur != nullptr)
      *alloc_cur = _cs_mem_global_alloc_cur;
    if (alloc_max != nullptr)
      *alloc_max = _cs_mem_global_alloc_max;
    if (n_allocs != nullptr)
      *n_allocs = _cs_mem_global_n_allocs;
    if (n_reallocs != nullptr)
      *n_reallocs = _cs_mem_global_n_reallocs;
    if (n_frees != nullptr)
      *n_frees = _cs_mem_global_n_frees;
    if (n_current != nullptr)
      *n_current = _cs_alloc_map.size();

    retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, cs_mem_memalign() will always fail.
 *
 * \returns 1 if memory aligned allocation is possible, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_mem_have_memalign(void)
{
#if defined(HAVE_POSIX_MEMALIGN)
  return 1;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Returns the error handler associated with the cs_mem_...() functions.
 *
 * \return pointer to the error handler function.
 */
/*----------------------------------------------------------------------------*/

bft_error_handler_t *
cs_mem_error_handler_get(void)
{
  return _cs_mem_error_handler;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Associates an error handler with the cs_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * \param handler pointer to the error handler function [in].
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_error_handler_set(bft_error_handler_t *handler)
{
  _cs_mem_error_handler = handler;
}

#if defined(HAVE_ACCEL)

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
cs_mem_malloc_hd(cs_alloc_mode_t   mode,
                 size_t            ni,
                 size_t            size,
                 const char       *var_name,
                 const char       *file_name,
                 int               line_num)
{
  if (ni == 0)
    return nullptr;

  cs_mem_block_t  me = {
    .host_ptr = nullptr,
    .device_ptr = nullptr,
    .size = ni * size,
    .mode = mode};

  if (mode < CS_ALLOC_HOST_DEVICE_PINNED) {
    me.host_ptr = cs_mem_malloc(ni, size, var_name, nullptr, 0);
  }

  // Device allocation will be postponed later thru call to
  // cs_get_device_ptr. This applies for CS_ALLOC_HOST_DEVICE
  // and CS_ALLOC_HOST_DEVICE_PINNED modes

#if defined(HAVE_CUDA)

  else if (mode == CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = cs_mem_cuda_malloc_host(me.size,
                                          var_name,
                                          file_name,
                                          line_num);

  else if (mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    me.host_ptr = cs_mem_cuda_malloc_managed(me.size,
                                             var_name,
                                             file_name,
                                             line_num);
    me.device_ptr = me.host_ptr;
  }

  else if (mode == CS_ALLOC_DEVICE)
    me.device_ptr = cs_mem_cuda_malloc_device(me.size,
                                              var_name,
                                              file_name,
                                              line_num);

#elif defined(SYCL_LANGUAGE_VERSION)

  else if (mode == CS_ALLOC_HOST_DEVICE_PINNED)
    me.host_ptr = _sycl_mem_malloc_host(me.size,
                                        var_name,
                                        file_name,
                                        line_num);

  else if (mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    me.host_ptr = _sycl_mem_malloc_shared(me.size,
                                          var_name,
                                          file_name,
                                          line_num);
    me.device_ptr = me.host_ptr;
  }

  else if (mode == CS_ALLOC_DEVICE)
    me.device_ptr = _sycl_mem_malloc_device(me.size,
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

  if (file_name != nullptr)
    _update_block_info(var_name, file_name, line_num,
                       nullptr, &me);

  /* Return pointer to allocated memory */

  if (me.host_ptr != nullptr)
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
cs_mem_realloc_hd(void            *ptr,
                  cs_alloc_mode_t  mode,
                  size_t           ni,
                  size_t           size,
                  const char      *var_name,
                  const char      *file_name,
                  int              line_num)
{
  void *ret_ptr = ptr;
  size_t new_size = ni*size;

  if (ptr == nullptr) {
    return cs_mem_malloc_hd(mode, ni, size, var_name, file_name, line_num);
  }
  else if (new_size == 0) {
    cs_mem_free(ptr, var_name, file_name, line_num);
    return nullptr;
  }

  cs_mem_block_t me = _get_block_info_try(ptr);

  if (new_size == me.size && mode == me.mode) {
    if (me.host_ptr != nullptr)
      return me.host_ptr;
    else
      return me.device_ptr;
  }

  cs_mem_block_t me_old = me;
  me.mode = mode;

  if (   me_old.mode <= CS_ALLOC_HOST_DEVICE
      && me.mode <= CS_ALLOC_HOST_DEVICE) {
    me.host_ptr = cs_mem_realloc(me_old.host_ptr, ni, size,
                                 var_name, nullptr, 0);
    me.size = new_size;
    ret_ptr = me.host_ptr;

    if (me.device_ptr != nullptr) {
#if defined(HAVE_CUDA)
      cs_mem_cuda_free(me.device_ptr, var_name, file_name, line_num);
#elif defined(SYCL_LANGUAGE_VERSION)
      sycl::free(me.device_ptr, cs_glob_sycl_queue);
#elif defined(HAVE_OPENMP_TARGET)
      omp_target_free(me.device_ptr, _omp_target_device_id);
#endif
      me.device_ptr = nullptr;
    }
  }

  else {
    size_t copy_size = me_old.size;
    if (new_size < copy_size)
      copy_size = new_size;
    me.size = new_size;
    me.device_ptr = nullptr;

    ret_ptr = cs_mem_malloc_hd(mode, 1, me.size,
                               var_name, nullptr, 0);

    if (me_old.mode < CS_ALLOC_DEVICE) {
      if (me.mode < CS_ALLOC_DEVICE)
        memcpy(ret_ptr, ptr, copy_size);
      else
        cs_copy_h2d(ret_ptr, ptr, copy_size);
    }
    else { /* if (me.mode == CS_ALLOC_DEVICE) */
      if (me.mode < CS_ALLOC_DEVICE)
        cs_copy_d2h(ret_ptr, ptr, copy_size);
      else
        cs_copy_d2d(ret_ptr, ptr, copy_size);
    }

    if (me.mode < CS_ALLOC_DEVICE) {
      me.host_ptr = ret_ptr;
      if (mode == CS_ALLOC_HOST_DEVICE_SHARED)
        me.device_ptr = ret_ptr;
    }
    else
      me.device_ptr = ret_ptr;

    cs_mem_free(ptr, var_name, nullptr, 0);
  }

  if (file_name != nullptr)
    _update_block_info(var_name, file_name, line_num,
                       &me_old, &me);

  return ret_ptr;
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
  if (ptr == nullptr)
    return nullptr;

  cs_mem_block_t me = _get_block_info_try(ptr);
  if (me.mode == CS_ALLOC_HOST) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p has no device association."), __func__, ptr);
    return nullptr;
  }

  /* Allocate on device if not done yet */

  if (me.device_ptr == nullptr) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {

      cs_mem_block_t me_old = me;

#if defined(HAVE_CUDA)

      me.device_ptr = cs_mem_cuda_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(SYCL_LANGUAGE_VERSION)

      me.device_ptr = _sycl_mem_malloc_device(me.size,
                                              "me.device_ptr",
                                              __FILE__,
                                              __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

      if (omp_target_associate_ptr(me.host_ptr, me.device_ptr, me.size, 0,
                                   _omp_target_device_id))
        bft_error(__FILE__, __LINE__, 0,
                  _("%s: Can't associate host pointer %p to device pointer %p."),
                  "omp_target_associate_ptr", me.host_ptr, me.device_ptr);

#endif

      _update_block_info("me.device_ptr", __FILE__, __LINE__,
                         &me_old, &me);
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
  if (ptr == nullptr)
    return nullptr;

  cs_mem_block_t me = _get_block_info_try(ptr);
  if (me.mode == CS_ALLOC_HOST) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p has no device association."), __func__, ptr);
    return nullptr;
  }

  /* Allocate and sync on device if not done yet */

  if (me.device_ptr == nullptr) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {

      cs_mem_block_t me_old = me;

#if defined(HAVE_CUDA)

      me.device_ptr = cs_mem_cuda_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(SYCL_LANGUAGE_VERSION)

      me.device_ptr = _sycl_mem_malloc_device(me.size,
                                              "me.device_ptr",
                                              __FILE__,
                                              __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

#endif

      _update_block_info("me.device_ptr", __FILE__, __LINE__,
                         &me_old, &me);
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
  if (ptr == nullptr)
    return nullptr;

  cs_mem_block_t me = _get_block_info_try(ptr);
  if (me.mode == CS_ALLOC_HOST) {
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p has no device association."), __func__, ptr);
    return nullptr;
  }

  /* Allocate and sync on device if not done yet */

  if (me.device_ptr == nullptr) {
    if (   me.mode == CS_ALLOC_HOST_DEVICE
        || me.mode == CS_ALLOC_HOST_DEVICE_PINNED) {

      cs_mem_block_t me_old = me;

#if defined(HAVE_CUDA)

      me.device_ptr = cs_mem_cuda_malloc_device(me.size,
                                                "me.device_ptr",
                                                __FILE__,
                                                __LINE__);

#elif defined(SYCL_LANGUAGE_VERSION)

      me.device_ptr = _sycl_mem_malloc_device(me.size,
                                              "me.device_ptr",
                                              __FILE__,
                                              __LINE__);

#elif defined(HAVE_OPENMP_TARGET)

      me.device_ptr = _omp_target_mem_malloc_device(me.size,
                                                    "me.device_ptr",
                                                    __FILE__,
                                                    __LINE__);

#endif

      _update_block_info("me.device_ptr", __FILE__, __LINE__,
                         &me_old, &me);
      cs_sync_h2d(ptr);

    }
  }

  /* Prefetch if shared */

  else if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED) {
    if (_ignore_prefetch == false)
      cs_prefetch_h2d(me.host_ptr, me.size);
  }

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
  cs_mem_block_t me = _get_block_info_try(ptr);

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
  cs_mem_block_t  me_old = _get_block_info_try(host_ptr);

  if (me_old.mode == CS_ALLOC_HOST) {

    cs_mem_block_t  me = {
      .host_ptr = host_ptr,
      .device_ptr = nullptr,
      .size = ni * size,
      .mode = CS_ALLOC_HOST_DEVICE};

    _update_block_info("me.host_ptr", __FILE__, __LINE__,
                       &me_old, &me);

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
  cs_mem_block_t me = _get_block_info_try(host_ptr);

  if (   me.device_ptr != nullptr
      && (   me.mode == CS_ALLOC_HOST_DEVICE
          || me.mode == CS_ALLOC_HOST_DEVICE_PINNED)) {

    cs_mem_block_t me_old = me;

#if defined(HAVE_CUDA)

    cs_mem_cuda_free(me.device_ptr, "me.device_ptr", __FILE__, __LINE__);

#elif defined(SYCL_LANGUAGE_VERSION)

    sycl::free(me.device_ptr, cs_glob_sycl_queue);

#elif defined(HAVE_OPENMP_TARGET)

    omp_target_disassociate_ptr(me.host_ptr, _omp_target_device_id);
    omp_target_free(me.device_ptr, _omp_target_device_id);

#endif

    me.device_ptr = nullptr;

    _update_block_info("me.device_ptr", __FILE__, __LINE__,
                       &me_old, &me);
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
  if (host_ptr == nullptr)
    return;

  void *ret_ptr = *host_ptr;

  void *_host_ptr = *host_ptr;

  if (_host_ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info(_host_ptr);

  if (mode != me.mode) {

    if (me.mode == CS_ALLOC_HOST_DEVICE)
      cs_disassociate_device_ptr(_host_ptr);

    if (   mode == CS_ALLOC_HOST_DEVICE_SHARED
        || me.mode == CS_ALLOC_HOST_DEVICE_SHARED) {

      ret_ptr = cs_mem_malloc_hd(mode, 1, me.size,
                                 "me.host_ptr", __FILE__, __LINE__);

      /* TODO: check if we have multiple OpenMP threads, in which
         case applying a "first-touch" policy might be useful here */

      if (me.mode < CS_ALLOC_DEVICE) {
        if (mode < CS_ALLOC_DEVICE)
          memcpy(ret_ptr, _host_ptr, me.size);
        else
          cs_copy_h2d(ret_ptr, _host_ptr, me.size);
      }
      else { /* if (old_mode == CS_ALLOC_DEVICE) */
        cs_copy_d2h(ret_ptr, _host_ptr, me.size);
      }

      cs_mem_free(_host_ptr, "me.host_ptr", __FILE__, __LINE__);

    }

    *host_ptr = ret_ptr;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_advise_set_read_mostly(void  *ptr)
{
  if (ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info_try(ptr);

  if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED) {

#if defined(HAVE_CUDA)

    cs_mem_cuda_set_advise_read_mostly(me.device_ptr, me.size);

#endif

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Advise memory system that a given allocation will be mostly read.
 *
 * \param [in]   ptr   pointer to allocation
 */
/*----------------------------------------------------------------------------*/

void
cs_mem_advise_unset_read_mostly(void  *ptr)
{
  if (ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info_try(ptr);

  if (me.mode == CS_ALLOC_HOST_DEVICE_SHARED) {

#if defined(HAVE_CUDA)

    cs_mem_cuda_unset_advise_read_mostly(me.device_ptr, me.size);

#endif

  }
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
  if (ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info(ptr);

  if (me.device_ptr == nullptr)
    me.device_ptr = const_cast<void *>(cs_get_device_ptr_const(ptr));

  switch (me.mode) {

  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;

  case CS_ALLOC_HOST_DEVICE:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_h2d(me.device_ptr, me.host_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.device_ptr, me.host_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      omp_target_memcpy(me.device_ptr, me.host_ptr, me.size, 0, 0,
                        _omp_target_device_id, omp_get_initial_device());
    }
    #endif
    break;

  case CS_ALLOC_HOST_DEVICE_PINNED:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_h2d_async(me.device_ptr, me.host_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.device_ptr, me.host_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      char *host_ptr = (char *)me.device_ptr;
      #pragma omp target enter data map(to:host_ptr[:me.size]) \
        nowait device(_omp_target_device_id)
    }
    #endif
    break;

  case CS_ALLOC_HOST_DEVICE_SHARED:
    if (_ignore_prefetch)
      return;

    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_prefetch_h2d(me.device_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.prefetch(me.device_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target enter data map(to:host_ptr[:me.size]) \
        nowait device(_omp_target_device_id)
    }
    #endif
    break;

  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize data from device to host.
 *
 * If separate allocations are used on the host and device
 * (mode == CS_ALLOC_HOST_DEVICE), the host pointer should be passed to this
 * function.
 *
 * Depending on the allocation type, this can imply a copy, data prefetch,
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
  if (ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info(ptr);

  switch (me.mode) {

  case CS_ALLOC_HOST:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on host only."),
              __func__, ptr);
    break;

  case CS_ALLOC_HOST_DEVICE:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_d2h(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      omp_target_memcpy(me.host_ptr, me.device_ptr, me.size, 0, 0,
                        omp_get_initial_device(), _omp_target_device_id);
    }
    #endif
    break;

  case CS_ALLOC_HOST_DEVICE_PINNED:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_d2h_async(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target exit data map(from:host_ptr[:me.size]) \
        nowait device(_omp_target_device_id)
    }
    #endif
    break;

  case CS_ALLOC_HOST_DEVICE_SHARED:
    if (_ignore_prefetch)
      return;

    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_prefetch_d2h(me.host_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.prefetch(me.host_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target exit data map(from:host_ptr[:me.size]) \
        nowait device(_omp_target_device_id)
    }
    #endif
    break;

  case CS_ALLOC_DEVICE:
    bft_error(__FILE__, __LINE__, 0,
              _("%s: %p allocated on device only."),
              __func__, ptr);
    break;

  }
}

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

void
cs_sync_d2h_if_needed(void  *ptr)
{
  if (ptr == nullptr)
    return;

  cs_mem_block_t me = _get_block_info_try(ptr);

  switch (me.mode) {

  case CS_ALLOC_HOST_DEVICE:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_d2h(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      omp_target_memcpy(me.host_ptr, me.device_ptr, me.size, 0, 0,
                        omp_get_initial_device(), _omp_target_device_id);
    }
    #endif
    break;

  case CS_ALLOC_HOST_DEVICE_PINNED:
    #if defined(HAVE_CUDA)
    {
      cs_mem_cuda_copy_d2h_async(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(SYCL_LANGUAGE_VERSION)
    {
      cs_glob_sycl_queue.memcpy(me.host_ptr, me.device_ptr, me.size);
    }
    #elif defined(HAVE_OPENMP_TARGET)
    {
      char *host_ptr = (char *)me.host_ptr;
      #pragma omp target exit data map(from:host_ptr[:me.size]) \
        nowait device(_omp_target_device_id)
    }
    #endif
    break;

  default:
    break;

  }
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
  if (ptr == nullptr)
    return;

#if defined(HAVE_CUDA)

  cs_mem_cuda_prefetch_h2d(ptr, size);

#elif defined(SYCL_LANGUAGE_VERSION)

  cs_glob_sycl_queue.prefetch(ptr, size);

#elif defined(HAVE_OPENMP_TARGET)

  char *host_ptr = (char *)ptr;
  #pragma omp target enter data map(to:host_ptr[:size]) \
    nowait device(_omp_target_device_id)

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
  if (ptr == nullptr)
    return;

#if defined(HAVE_CUDA)

  cs_mem_cuda_prefetch_d2h(ptr, size);

#elif defined(SYCL_LANGUAGE_VERSION)

  cs_glob_sycl_queue.prefetch(ptr, size);

#elif defined(HAVE_OPENMP_TARGET)

  char *host_ptr = (char *)ptr;
  #pragma omp target exit data map(from:host_ptr[:size]) \
    nowait device(_omp_target_device_id)

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
  if (src == nullptr)
    return;

#if defined(HAVE_CUDA)

  cs_mem_cuda_copy_h2d(dest, src, size);

#elif defined(SYCL_LANGUAGE_VERSION)

  cs_glob_sycl_queue.memcpy(dest, src, size);

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    _omp_target_device_id, omp_get_initial_device());

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
  if (src == nullptr)
    return;

#if defined(HAVE_CUDA)

  cs_mem_cuda_copy_d2h(dest, src, size);

#elif defined(SYCL_LANGUAGE_VERSION)

  cs_glob_sycl_queue.memcpy(dest, src, size);

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    omp_get_initial_device(), _omp_target_device_id);

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
  if (src == nullptr)
    return;

#if defined(HAVE_CUDA)

  cs_mem_cuda_copy_d2d(dest, src, size);

#elif defined(SYCL_LANGUAGE_VERSION)

  cs_glob_sycl_queue.memcpy(dest, src, size);

#elif defined(HAVE_OPENMP_TARGET)

  omp_target_memcpy(dest, src, size, 0, 0,
                    _omp_target_device_id, _omp_target_device_id);

#endif
}

#endif // defined(HAVE_ACCEL)

/*----------------------------------------------------------------------------*/

END_C_DECLS
