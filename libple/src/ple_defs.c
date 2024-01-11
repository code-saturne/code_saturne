/*============================================================================
 * Base definitions and utility functions for PLE.
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2024  EDF S.A.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*!
 * \file ple_defs.c
 *
 * \brief Configurable error and memory handling with default functions.
 */

/*-----------------------------------------------------------------------------*/

#include "ple_config_defs.h"

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

#if (__STDC_VERSION__ <202311L)
# include <stdbool.h>
#endif

#if defined (HAVE_GETTIMEOFDAY)
#include <sys/time.h>
#endif

#if defined (HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#elif defined(_POSIX_SOURCE)
#include <sys/times.h>
#include <unistd.h>
#endif

/*
 * Optional library and PLE headers
 */

#include "ple_defs.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-------------------------------------------------------------------------------
 * Local macro documentation
 *-----------------------------------------------------------------------------*/

/*! \fn PLE_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls ple_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn PLE_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls ple_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn PLE_FREE(_ptr)
 * \brief Free allocated memory.
 *
 * This macro calls ple_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 */

/*-----------------------------------------------------------------------------
 * Internationalization (future)
 *-----------------------------------------------------------------------------*/

#if defined(ENABLE_NLS)

#  include <libintl.h>
#  define _(String) dgettext(PACKAGE,String)
#  ifdef gettext_noop
#    define N_(String) gettext_noop(String)
#  else
#    define N_(String) String
#  endif /* gettext_noop */

#else

#  define _(String) (String)
#  define N_(String) String
#  define textdomain(String) (String)
#  define gettext(String) (String)
#  define dgettext(Domain,String) (String)
#  define dcgettext(Domain,String,Type) (String)
#  define bindtextdomain(Domain,Directory) (Domain)

#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Defauly error handler.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_ple_error_default(const char  *file_name,
                   int          line_num,
                   int          sys_error_code,
                   const char  *format,
                   va_list      arg_ptr)
{
  fprintf(stderr, "\n");

  if (sys_error_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_error_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  assert(0);

  exit(EXIT_FAILURE);
}

/*
 * Default memory allocation.
 *
 * This function simply wraps malloc() and calls _ple_error() if it fails
 * to allocate the required memory.
 *
 * parameters:
 *   ni        <-- number of elements.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 */

static void *
_ple_mem_malloc_default(size_t       ni,
                        size_t       size,
                        const char  *var_name,
                        const char  *file_name,
                        int          line_num)
{
  void  *p_ret;
  size_t  alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  /* Allocate memory and check return */

  p_ret = malloc(alloc_size);

  if (p_ret == NULL)
    ple_error(file_name, line_num, errno,
              _("Failure to allocate \"%s\" (%lu bytes)"),
              var_name, (unsigned long)alloc_size);

  return p_ret;
}

/*
 * Default memory reallocation.
 *
 * This function simply wraps realloc() and calls _ple_error() if it fails
 * to reallocate the required memory.
 *
 * parameters:
 *   ptr       <-- pointer to previous memory location
 *   ni        <-- number of elements.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to reallocated memory.
 */

static void *
_ple_mem_realloc_default(void        *ptr,
                         size_t       ni,
                         size_t       size,
                         const char  *var_name,
                         const char  *file_name,
                         int          line_num)
{
  void  *p_ret;

  size_t realloc_size = ni * size;

  p_ret = realloc(ptr, realloc_size);

  if (realloc_size != 0 && p_ret == NULL)
    ple_error(file_name, line_num, errno,
              _("Failure to reallocate \"%s\" (%lu bytes)"),
              var_name, (unsigned long)realloc_size);

  return p_ret;
}

/*
 * Default memory free.
 *
 * This function simply wraps free().
 *
 * parameters:
 *   ptr       <-- pointer to previous memory location
 *   var_name  <-- allocated variable name string
 *   file_name <-- name of calling source file
 *   line_num  <-- line number in calling source file
 *
 * returns:
 *   NULL pointer.
 */

static void *
_ple_mem_free_default(void        *ptr,
                      const char  *var_name,
                      const char  *file_name,
                      int          line_num)
{
  PLE_UNUSED(var_name);
  PLE_UNUSED(file_name);
  PLE_UNUSED(line_num);

  if (ptr != NULL)
    free(ptr);

  return NULL;
}

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

/* Standard output handler */

static ple_printf_t *_ple_printf = vprintf;

/* Error handler */

static ple_error_handler_t *_ple_error = _ple_error_default;

/* Memory allocation */

static ple_mem_malloc_t *_ple_mem_malloc = _ple_mem_malloc_default;
static ple_mem_realloc_t *_ple_mem_realloc = _ple_mem_realloc_default;
static ple_mem_free_t *_ple_mem_free = _ple_mem_free_default;

/* Timer */

static bool _ple_timer_initialized = false;

#if defined (HAVE_GETTIMEOFDAY)
static struct timeval  _ple_timer_wtime_tv_start;
#else
static time_t _ple_timer_wtime_start;
#endif

#if defined (HAVE_GETRUSAGE)
#elif defined(_POSIX_SOURCE)
static time_t _ple_timer_unit = 0;
#else
static clock_t _ple_timer_clock_start;
#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

static void
_ple_timer_initialize(void)
{
#if defined (HAVE_GETTIMEOFDAY)
  (void)gettimeofday(&_ple_timer_wtime_tv_start, NULL);
#else
  (void)time(&_ple_timer_wtime_start);
#endif

#if defined (HAVE_GETRUSAGE)
#elif defined(_POSIX_SOURCE)
  _ple_timer_unit = (double)sysconf(_SC_CLK_TCK);
#else
  _ple_timer_clock_start = clock();
#endif

  _ple_timer_initialized = true;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by ple_printf_function_set().
 *
 * \param [in] format format string, as printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int
ple_printf(const char *const format,
           ...)
{
  int  retval;
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  retval = _ple_printf(format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*!
 * \brief Returns function associated with the ple_printf() function.
 *
 * \return pointer to the vprintf() or replacement function.
 */

ple_printf_t *
ple_printf_function_get(void)
{
  return _ple_printf;
}

/*!
 * \brief Associates a vprintf() type function with the ple_printf() function.
 *
 * \param [in] fct pointer to a vprintf() type function.
 */

void
ple_printf_function_set(ple_printf_t *const fct)
{
  _ple_printf = fct;
}

/*!
 * \brief Calls the error handler (set by ple_error_handler_set() or default).
 *
 * With the default error handler, an error message is output to stderr,
 * and the current process exits with an EXIT_FAILURE code.
 *
 * \param [in] file_name      name of source file from which error handler
 *                            called.
 * \param [in] line_num       line of source file from which error handler
 *                            called.
 * \param [in] sys_error_code error code if error in system or libc call,
 *                            0 otherwise.
 * \param [in] format         format string, as printf() and family.
 * \param [in] ...            variable arguments based on format string.
 */

void
ple_error(const char  *const file_name,
          const int          line_num,
          const int          sys_error_code,
          const char  *const format,
          ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  _ple_error(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*!
 * \brief Returns the error handler associated with the ple_error() function.
 *
 * \return pointer to the error handler function.
 */

ple_error_handler_t *
ple_error_handler_get(void)
{
  return _ple_error;
}

/*!
 * \brief Associates an error handler with the ple_error() function.
 *
 * \param [in] handler pointer to the error handler function.
 */

void
ple_error_handler_set(ple_error_handler_t  *handler)
{
  _ple_error = handler;
}

/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * allocate the required memory.
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
ple_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num)
{
  return _ple_mem_malloc(ni, size, var_name, file_name, line_num);
}

/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, ple_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */

void *
ple_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
  return _ple_mem_realloc(ptr, ni, size, var_name, file_name, line_num);
}

/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, ple_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns NULL pointer.
 */

void *
ple_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
  return _ple_mem_free(ptr, var_name, file_name, line_num);
}

/*!
 * \brief Return the function pointers associated with PLE's memory management.
 *
 * All arguments are optional.
 *
 * \param [out] malloc_func  pointer to ple_mem_malloc function pointer
 *                           (or NULL).
 * \param [out] realloc_func pointer to ple_mem_realloc function pointer
 *                           (or NULL).
 * \param [out] free_func    pointer to ple_mem_free function pointer
 *                           (or NULL).
 */

void
ple_mem_functions_get(ple_mem_malloc_t   **malloc_func,
                      ple_mem_realloc_t  **realloc_func,
                      ple_mem_free_t     **free_func)
{
  if (malloc_func != NULL)
    *malloc_func = _ple_mem_malloc;

  if (realloc_func != NULL)
    *realloc_func = _ple_mem_realloc;

  if (free_func != NULL)
    *free_func = _ple_mem_free;
}

/*!
 * \brief Associate functions to modifiy PLE's memory management.
 *
 * All arguments are optional, so the previously set functions pointers will
 * not be modified if an argument value is NULL.
 *
 * \param [in] malloc_func  ple_mem_malloc function pointer (or NULL).
 * \param [in] realloc_func ple_mem_realloc function pointer (or NULL).
 * \param [in] free_func    ple_mem_free function pointer (or NULL).
 */

void
ple_mem_functions_set(ple_mem_malloc_t   *malloc_func,
                      ple_mem_realloc_t  *realloc_func,
                      ple_mem_free_t     *free_func)
{
  if (malloc_func != NULL)
    _ple_mem_malloc = malloc_func;

  if (realloc_func != NULL)
    _ple_mem_realloc = realloc_func;

  if (free_func != NULL)
    _ple_mem_free = free_func;
}

/*!
 * \brief Return Wall clock time
 *
 * \return elapsed time from first call of a function of the ple_timer_...()
 *         series, or -1 if unable to compute.
 */

double
ple_timer_wtime(void)
{
  double this_wtime = -1.;

  /* Ensure initialization */

  if (_ple_timer_initialized == false)
    _ple_timer_initialize();

  /* Compute elapsed time */

#if defined (HAVE_GETTIMEOFDAY)

  {
    struct timeval  wtime_tv_current;

    if (gettimeofday(&wtime_tv_current, NULL) == 0) {

      /* Perform carry for later subtraction */
      if (_ple_timer_wtime_tv_start.tv_usec > wtime_tv_current.tv_usec) {
        int nsec = (_ple_timer_wtime_tv_start.tv_usec - wtime_tv_current.tv_usec)
                   / 1000000 + 1;
        wtime_tv_current.tv_usec += 1000000 * nsec;
        wtime_tv_current.tv_sec -= nsec;
      }
      if (  wtime_tv_current.tv_usec - _ple_timer_wtime_tv_start.tv_usec
          > 1000000) {
        int nsec = (wtime_tv_current.tv_usec - _ple_timer_wtime_tv_start.tv_usec)
                   / 1000000;
        wtime_tv_current.tv_usec -= 1000000 * nsec;
        wtime_tv_current.tv_sec += nsec;
      }

      this_wtime =   (  wtime_tv_current.tv_sec
                      - _ple_timer_wtime_tv_start.tv_sec)
                   + (  wtime_tv_current.tv_usec
                      - _ple_timer_wtime_tv_start.tv_usec) * 1.e-6 ;
    }
  }

#else

  {
    time_t wtime_current;

    if (time(&wtime_current) != (time_t)-1)
      this_wtime = difftime(wtime_current, _ple_timer_wtime_start);
  }

#endif

  return this_wtime;
}

/*!
 * \brief Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see ple_timer_cpu_time_method()), at least one of
 * the ple_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * \return current CPU time usage, or -1 if unable to compute.
 */

double
ple_timer_cpu_time(void)
{
  double cpu_time = -1.;

  /* Ensure initialization */

  if (_ple_timer_initialized == 0)
    _ple_timer_initialize();

  /* Compute CPU time */

#if defined (HAVE_GETRUSAGE)

  {
    struct rusage  usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
      cpu_time  =    usage.ru_utime.tv_sec  + usage.ru_stime.tv_sec
                  + (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) * 1.e-6;
    }
  }

#elif defined(_POSIX_SOURCE)

  {
    static struct tms  ptimer;

    if (_ple_timer_unit != -1 && times(&ptimer) != -1) {
      cpu_time =   ((double)(ptimer.tms_utime + ptimer.tms_stime))
                 / _ple_timer_unit;
    }
  }

#else /* Use minimal C library function */

  if (_ple_timer_clock_start != -1) {

    static clock_t  clock_current;

    clock_current = clock();
    if (clock_current != (clock_t)-1)
      cpu_time
        = ((double)(clock_current - _ple_timer_clock_start)) / CLOCKS_PER_SEC;

  }

#endif

  return cpu_time;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
