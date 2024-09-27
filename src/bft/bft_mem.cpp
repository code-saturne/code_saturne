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

/*-----------------------------------------------------------------------------*/

/*
 * Standard C and C++ library headers
 */

#include <map>

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Optional library and BFT headers
 */

#include "bft_error.h"
#include "bft_mem_usage.h"
#include "bft_mem.h"
#include "bft_printf.h"

/*-----------------------------------------------------------------------------*/

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file bft_mem.c
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

/*! \fn BFT_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls bft_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn BFT_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls bft_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn BFT_FREE(_ptr)
 * \brief Free allocated memory.
 *
 * This macro calls bft_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to nullptr to avoid accidental reuse.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 */

/*! \fn BFT_MEMALIGN(_ptr, _align, _ni, _type)
 * \brief Allocate aligned memory for _ni elements of type _type.
 *
 * This macro calls bft_mem_memalign(), automatically setting the
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
_bft_mem_error_handler_default(const char  *file_name,
                               int          line_num,
                               int          sys_error_code,
                               const char  *format,
                               va_list      arg_ptr);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static int  _bft_mem_global_init_mode = 0;

static FILE *_bft_mem_global_file = nullptr;

static std::map<const void *, cs_mem_block_t> _bft_alloc_map;

static size_t  _bft_mem_global_alloc_cur = 0;
static size_t  _bft_mem_global_alloc_max = 0;

static size_t  _bft_mem_global_n_allocs = 0;
static size_t  _bft_mem_global_n_reallocs = 0;
static size_t  _bft_mem_global_n_frees = 0;

static bft_error_handler_t  *_bft_mem_error_handler
                              = (_bft_mem_error_handler_default);

static bft_mem_realloc_t   *_bft_alt_realloc_func = nullptr;
static bft_mem_free_t      *_bft_alt_free_func = nullptr;

#if defined(HAVE_OPENMP)
static omp_lock_t _bft_mem_lock;
#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Given a character string representing a file name, returns
 * pointer to that part of the string corresponding to the base name.
 *
 * parameters:
 *   file_name:      <-- full name of source file.
 *
 * return:
 *   pointer to part of file name corresponding to base name.
 */

static const char *
_bft_mem_basename(const char  *file_name)
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

/*
 * Determines values associated with an array representing a
 * long integer
 *
 * parameters:
 *   counter: <-- counter to update.
 *   value:   --> counter values in output unit [out].
 *   unit:    --> counter value unit : ' ', 'k', 'm', 'g', 't', or 'p'
 *                for bytes, Kilobytes, Megabytes, Gigabytes, Terabytes,
 *                or Petabytes [out].
 */

static void
_bft_mem_size_val(const size_t    counter,
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

/*
 * Memory usage summary.
 */

static void _bft_mem_summary(FILE  *f)
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

  _bft_mem_size_val(_bft_mem_global_alloc_cur, value, &unit);
  fprintf(f,
          "Theoretical current allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  _bft_mem_size_val(_bft_mem_global_alloc_max, value, &unit);
  fprintf(f,
          "Theoretical maximum allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  fprintf(f,
          "\n"
          "Number of allocations:   %lu\n"
          "          reallocations: %lu\n"
          "          frees:         %lu\n\n",
          (unsigned long)_bft_mem_global_n_allocs,
          (unsigned long)_bft_mem_global_n_reallocs,
          (unsigned long)_bft_mem_global_n_frees);

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

/*
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
 */

static void
_bft_mem_error_handler_default(const char  *file_name,
                               int          line_num,
                               int          sys_error_code,
                               const char  *format,
                               va_list      arg_ptr)
{
  bft_error_handler_t * general_err_handler;

  bft_printf_flush();

  _bft_mem_summary(stderr);

  general_err_handler = bft_error_handler_get();
  general_err_handler(file_name, line_num, sys_error_code, format, arg_ptr);
}

/*
 * Calls the error handler (set by bft_mem_error_handler_set() or default).
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameters:
 *   file_name      <-- name of source file from which failed bft_mem_...()
 *                      function was called.
 *   line_num       <-- line of source file from which failed bft_mem_...()
 *                      function was called.
 *   sys_error_code <-- error code if error in system or libc call,
 *                      0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   ...            <-- variable arguments based on format string.
 */

static void
_bft_mem_error(const char  *file_name,
               int          line_num,
               int          sys_error_code,
               const char  *format,
               ...)
{
  va_list  arg_ptr;

  if (_bft_mem_global_file != nullptr) {
    _bft_mem_summary(_bft_mem_global_file);
    fflush(_bft_mem_global_file);
  }

  va_start(arg_ptr, format);

  _bft_mem_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*
 * Call error function when pointer is not found in allocation info.
 *
 * parameters:
 *   p: <-- allocated block's start adress.
 *
 * returns:
 *   corresponding _bft_mem_block structure.
 */

static void
_bft_mem_block_info_error(const void  *p)
{
  _bft_mem_error(__FILE__, __LINE__, 0,
                 _("Adress [%p] does not correspond to "
                   "the beginning of an allocated block."),
                 p);
}

/*
 * Fill a cs_mem_block_t structure for an allocated pointer.
 */

static inline cs_mem_block_t
_bft_mem_block_new(void          *p_new,
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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Semi private function definitions
 *============================================================================*/

/*!
 * \brief Return the cs_mem_block structure corresponding to a given
 * allocated block.
 *
 * \param [in]  p_get  allocated block's start adress.
 *
 * \return  corresponding cs_mem_block structure.
 */

cs_mem_block_t
bft_mem_get_block_info(const void  *p_get)
{
  auto it = _bft_alloc_map.find(p_get);

  if (it == _bft_alloc_map.end())
    _bft_mem_block_info_error(p_get);

  return it->second;
}

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

cs_mem_block_t
bft_mem_get_block_info_try(const void  *p_get)
{
  cs_mem_block_t mbi;

  auto it = _bft_alloc_map.find(p_get);
  if (it != _bft_alloc_map.end())
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

/*!
 * \brief Log matching memory operation if logging is enabled
 *
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 * \param [in] old_block pointer to old block info, if present
 * \param [in] new_block pointer to new block info, if present
 */

void
bft_mem_update_block_info(const char            *var_name,
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

  if (   _bft_mem_global_init_mode > 1
      ||  old_mode > CS_ALLOC_HOST || new_mode > CS_ALLOC_HOST) {

#if defined(HAVE_OPENMP)
    int in_parallel = omp_in_parallel();
    if (in_parallel)
      omp_set_lock(&_bft_mem_lock);
#endif

    /* Update map */

    if (old_block != nullptr && (new_block == nullptr || p_m_new != p_m_old))
      _bft_alloc_map.erase(p_m_old);
    if (new_block != nullptr && p_m_new != nullptr)
      _bft_alloc_map[p_m_new] = *new_block;

    /* Memory allocation counting */

    if (_bft_mem_global_init_mode > 1) {
      long long size_diff = new_size - old_size;

      int cat_id = 3;
      if (p_m_old == nullptr) {
        cat_id = 0;
        _bft_mem_global_n_allocs += 1;
      }
      else if (p_m_new == nullptr) {
        cat_id = 2;
        _bft_mem_global_n_frees += 1;
      }
      else if (old_size != new_size) {
        cat_id = 1;
        _bft_mem_global_n_reallocs += 1;
      }

      _bft_mem_global_alloc_cur += size_diff;

      if (_bft_mem_global_alloc_max < _bft_mem_global_alloc_cur)
        _bft_mem_global_alloc_max = _bft_mem_global_alloc_cur;

      /* Log to file */

      if (_bft_mem_global_file != nullptr) {

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

        fprintf(_bft_mem_global_file, "\n%s: %-27s:%6d : %-39s: %9lu",
                cat_s[cat_id], _bft_mem_basename(file_name), line_num,
                var_name, (unsigned long)new_size);
        fprintf(_bft_mem_global_file, " : (%c%9lu) : %12lu : [%14p] : [%14p]",
                c_sgn, (unsigned long)size_diff,
                (unsigned long)_bft_mem_global_alloc_cur,
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
        fprintf(_bft_mem_global_file, " : [%14p] : [%14p] : %s",
                p_old_d, p_new_d, alloc_mode_s[new_mode]);
#endif
        fflush(_bft_mem_global_file);

      } /* End of log to file */

    } /* End counting and logging */

#if defined(HAVE_OPENMP)
    if (in_parallel)
      omp_unset_lock(&_bft_mem_lock);
#endif

  } /* End if map needs to be updated */
}

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Initialize memory handling.
 *
 * This function should be called before any other bft_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * \param log_file_name name of optional log_file (if nullptr, no log).
 */

void
bft_mem_init(const char *log_file_name)
{
#if defined(HAVE_OPENMP)
  if (omp_in_parallel()) {
    if (omp_get_thread_num() != 0)
      return;
  }
  omp_init_lock(&_bft_mem_lock);
#endif

  if (_bft_mem_global_init_mode > 0) {
    _bft_mem_error(__FILE__, __LINE__, 0,
                   _("bft_mem_init() has already been called"));
  }

#if defined(HAVE_ACCEL)
  _bft_mem_global_init_mode = 2;
#else
  _bft_mem_global_init_mode = 1;
#endif

  if (log_file_name != nullptr) {

    _bft_mem_global_init_mode = 2;

    _bft_mem_global_file = fopen(log_file_name, "w");

    /*
      If the file could not be opened, we do not abort, as it is not
      absolutely necessary. We silently continue.
      (We could warn the user, but this would require either using
      bft_printf(), which we prefer to keep independent of the bft_mem_...()
      functions to avoid evental crossed definitions when user-defined, or
      "warning handling" similar to error handling, with a possibility
      of user-defined warning handlers, as we are not sure if the calling
      code uses stderr (especially in a distributed environment). This
      is probably not worth the bother.
    */

    if (_bft_mem_global_file == nullptr)
      fprintf(stderr,
              _("Failure to open memory log file \"%s\"\n"),
              log_file_name);

  }

  /* Log file header */

  if (_bft_mem_global_file != nullptr) {

    fprintf(_bft_mem_global_file,
            "       :     FILE NAME              : LINE  :"
            "  POINTER NAME                          : N BYTES   :"
            " (+- N BYTES) : TOTAL BYTES  : [   OLD ADRESS ] : [   NEW ADRESS ]");
#if defined(HAVE_ACCEL)
    fprintf(_bft_mem_global_file,
            " : [OLD ADRESS (D)] : [NEW ADRESS (D)] : MODE");
#endif
    fprintf(_bft_mem_global_file, "\n");

    fprintf(_bft_mem_global_file,
            "-------:----------------------------:-------:"
            "----------------------------------------:-----------:"
            "-----------------------------:------------------:-----------------");
#if defined(HAVE_ACCEL)
    fprintf(_bft_mem_global_file,
            "-:------------------:------------------:-----");
#endif
    fprintf(_bft_mem_global_file, "\n");

  }

}

/*!
 * \brief End memory handling.
 *
 * This function should be called after all other bft_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void
bft_mem_end(void)
{
  if (_bft_mem_global_init_mode == 0)
    return;

#if defined(HAVE_OPENMP)
  if (omp_in_parallel()) {
    if (omp_get_thread_num() != 0)
      return;
  }
  omp_destroy_lock(&_bft_mem_lock);
#endif

  _bft_mem_global_init_mode = 0;

  if (_bft_mem_global_file != nullptr) {

    /* Memory usage summary */

    _bft_mem_summary(_bft_mem_global_file);

    /* List of non-freed pointers */

    size_t non_free = _bft_alloc_map.size();

    if (non_free > 0) {

      fprintf(_bft_mem_global_file, "List of non freed pointers:\n");

      for (auto const& x : _bft_alloc_map) {
        const void *p = x.first;
        // const cs_mem_block_t b = x.second;

        fprintf(_bft_mem_global_file,"[%p]\n", p);
      }

    }

    fprintf(_bft_mem_global_file,
            "Number of non freed pointers remaining: %lu\n",
            non_free);

    fclose(_bft_mem_global_file);
  }

  /* Reset defaults in case of later initialization */

  _bft_alloc_map.clear();

  _bft_mem_global_n_allocs = 0;
  _bft_mem_global_n_reallocs = 0;
  _bft_mem_global_n_frees = 0;
}

/*!
 * \brief Indicates if bft_mem_...() functions are initialized.
 *
 * \returns 1 if bft_mem_init has been called, 0 otherwise.
 */

int
bft_mem_initialized(void)
{
  return (_bft_mem_global_init_mode > 0) ? 1 : 0;
}

/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * Allocation couting and logging to trace file will be done if
 * both required by the bft_mem_init options and if file_name != nullptr.
 * If required but file_name == nullptr, it must be handled by the caller,
 * using \ref bft_mem_log_mem_op.
 *
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */

void *
bft_mem_malloc(size_t       ni,
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
    _bft_mem_error(file_name, line_num, errno,
                   _("Failure to allocate \"%s\" (%lu bytes)"),
                   var_name, (unsigned long)alloc_size);
    return nullptr;
  }
  else if (_bft_mem_global_init_mode < 2)
    return p_new;

  cs_mem_block_t mib = _bft_mem_block_new(p_new, alloc_size);

  if (file_name != nullptr)
    bft_mem_update_block_info(var_name,
                              file_name,
                              line_num,
                              nullptr,
                              &mib);

  /* Return pointer to allocated memory */

  return p_new;
}

/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if nullptr, bft_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */

void *
bft_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
  size_t new_size = ni * size;

  /*
    Behave as bft_mem_malloc() if the previous pointer is equal to nullptr.
    Note that the operation will then appear as a first allocation
    ('alloc') in the _bft_mem_global_file trace file.
  */

  if (ptr == nullptr)
    return bft_mem_malloc(ni,
                          size,
                          var_name,
                          file_name,
                          line_num);

  /*
    Behave as bft_mem_free() if the requested size is zero.
    Note that in this case, the operation will appear as 'free'
    in the _bft_mem_global_file trace file.
  */

  else if (ni == 0)
    return bft_mem_free(ptr,
                        var_name,
                        file_name,
                        line_num);

  /* When possible, get previous allocation information. */

  cs_mem_block_t mib_old;
  if (_bft_mem_global_init_mode > 1)
    mib_old = bft_mem_get_block_info(ptr);
  else
    mib_old = bft_mem_get_block_info_try(ptr);

  /* If the old size is known to equal the new size,
     nothing needs to be done. */

  if (new_size == mib_old.size) {
    return ptr;
  }

  /* In the general case, we have a true reallocation. */

#if defined(HAVE_ACCEL)
  if (mib_old.mode >= CS_ALLOC_HOST_DEVICE_PINNED) {
    return _bft_alt_realloc_func(ptr, ni, size,
                                 var_name, file_name, line_num);
  }
#endif

  void *p_new = realloc(ptr, new_size);

  if (file_name != nullptr) {
    cs_mem_block_t mib_new = _bft_mem_block_new(p_new, new_size);

    bft_mem_update_block_info(var_name,
                              file_name,
                              line_num,
                              &mib_old,
                              &mib_new);
  }

  return p_new;
}

/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * free the corresponding memory. In case of a nullptr pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if nullptr, bft_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns nullptr pointer.
 */

void *
bft_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
  /* nullptr pointer case (non-allocated location) */

  if (ptr == nullptr)
    return nullptr;

  /* General case (free allocated memory) */

  /* When possible, get previous allocation information. */

  cs_mem_block_t mib_old;
  if (_bft_mem_global_init_mode > 1)
    mib_old = bft_mem_get_block_info(ptr);
  else
    mib_old = bft_mem_get_block_info_try(ptr);

#if defined(HAVE_ACCEL)
  if (mib_old.mode >= CS_ALLOC_HOST_DEVICE_PINNED) {
    _bft_alt_free_func(ptr, var_name, file_name, line_num);
    return nullptr;
  }
#endif

  free(ptr);
  if (   mib_old.host_ptr != nullptr
      && file_name != nullptr)
    bft_mem_update_block_info(var_name,
                              file_name,
                              line_num,
                              &mib_old,
                              nullptr);

  return nullptr;
}

/*!
 * \brief Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the bft_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function bft_mem_have_memalign() indicates if this
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

void *
bft_mem_memalign(size_t       alignment,
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
      _bft_mem_error(file_name, line_num, 0,
                     _("Alignment %lu for \"%s\" not a power of 2\n"
                       "or a multiple of sizeof(void *) = %lu"),
                     (unsigned long)alignment, var_name,
                     (unsigned long)(sizeof(void *)));
      break;
    default:
      _bft_mem_error(file_name, line_num, 0,
                     _("Failure to allocate \"%s\" (%lu bytes)"),
                     var_name, (unsigned long)alloc_size);
    }
    return nullptr;
  }
  else if (_bft_mem_global_init_mode < 2)
    return p_loc;

  cs_mem_block_t mib = _bft_mem_block_new(p_loc, alloc_size);

  if (file_name != nullptr)
    bft_mem_update_block_info(var_name,
                              file_name,
                              line_num,
                              nullptr,
                              &mib);

  /* Return pointer to allocated memory */

  return p_loc;

#else

  _bft_mem_error(file_name, line_num, errno,
                 _("No aligned allocation function available on this system"));

  return nullptr;

#endif
}

/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through bft_mem_...() (in kB).
 */

size_t
bft_mem_size_current(void)
{
  return (_bft_mem_global_alloc_cur / 1024);
}

/*!
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through bft_mem_...() (in kB).
 */

size_t
bft_mem_size_max(void)
{
  return (_bft_mem_global_alloc_max / 1024);
}

/*!
 * \brief Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, bft_mem_memalign() will always fail.
 *
 * \returns 1 if memory aligned allocation is possible, 0 otherwise.
 */

int
bft_mem_have_memalign(void)
{
#if defined(HAVE_POSIX_MEMALIGN)
  return 1;
#else
  return 0;
#endif
}

/*!
 * \brief Returns the error handler associated with the bft_mem_...() functions.
 *
 * \return pointer to the error handler function.
 */

bft_error_handler_t *
bft_mem_error_handler_get(void)
{
  return _bft_mem_error_handler;
}

/*!
 * \brief Associates an error handler with the bft_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * \param handler pointer to the error handler function [in].
 */

void
bft_mem_error_handler_set(bft_error_handler_t *handler)
{
  _bft_mem_error_handler = handler;
}

/*
 * Associates alternative functions with the bft_mem_...() functions.
 *
 * When memory allocated with another mechanism is reallocated or
 * freed using a bft_mem_... function, this allows trying the
 * matching alternative function rather than throwing an error.
 *
 * Though using matching methods is recommended, this allows handling
 * compatibility between methods which might be used in different parts
 * of the code.
 *
 * parameter:
 *   realloc_func <-- pointer to alternative reallocation function.
 *   free_func    <-- pointer to alternative free function.
 */

void
bft_mem_alternative_set(bft_mem_realloc_t   *realloc_func,
                        bft_mem_free_t      *free_func)
{
  _bft_alt_realloc_func = realloc_func;
  _bft_alt_free_func = free_func;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
