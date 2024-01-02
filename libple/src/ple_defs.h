#ifndef __PLE_DEFS_H__
#define __PLE_DEFS_H__

/*============================================================================
 * Definitions, global variables, and base functions
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

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

/*----------------------------------------------------------------------------*/

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls ple_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define PLE_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) ple_mem_malloc(_ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls ple_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define PLE_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) ple_mem_realloc(_ptr, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls ple_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#ifdef __cplusplus /* avoid casting from void for C++ */

#define PLE_FREE(_ptr) \
ple_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

#else

#define PLE_FREE(_ptr) \
_ptr = ple_mem_free(_ptr, #_ptr, __FILE__, __LINE__)

#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * General C types such as size_t which should be known
 *----------------------------------------------------------------------------*/

/* Obtain definitions such as that of size_t */

#include <stddef.h>

/*----------------------------------------------------------------------------
 * Basic types used by PLE.
 * They may be modified here to better map to a given library, with the
 * following constraints:
 *  - ple_lnum_t must be signed
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_LONG_LNUM)
  typedef long     ple_lnum_t;     /* Local integer index or number */
#else
  typedef int      ple_lnum_t;     /* Local integer index or number */
#endif

typedef double   ple_coord_t;    /* Real number (coordinate value) */

/*----------------------------------------------------------------------------
 * MPI datatypes.
 *----------------------------------------------------------------------------*/

#if defined(PLE_HAVE_MPI)

#define PLE_MPI_TAG      (int)('P'+'L'+'E') /* MPI tag for PLE operations */

#if defined(PLE_HAVE_LONG_LNUM)
#  define PLE_MPI_LNUM   MPI_LONG        /* MPI type for ple_lnum_t type */
#else
#  define PLE_MPI_LNUM   MPI_INT         /* MPI type for ple_lnum_t type */
#endif

#define PLE_MPI_COORD    MPI_DOUBLE      /* MPI type for ple_coord_t type */

#endif

/*----------------------------------------------------------------------------
 * Macro used to silence "unused argument" warnings.
 *
 * This is useful when a function must match a given function pointer
 * type, but does not use all possible arguments.
 *----------------------------------------------------------------------------*/

#define PLE_UNUSED(x) (void)(x)

/*----------------------------------------------------------------------------
 * Macros for compilation with a C++ compiler
 *----------------------------------------------------------------------------*/

#undef PLE_BEGIN_C_DECLS
#undef   PLE_END_C_DECLS

#if defined(__cplusplus)
#  define PLE_BEGIN_C_DECLS  extern "C" {
#  define   PLE_END_C_DECLS  }
#else
#  define PLE_BEGIN_C_DECLS
#  define   PLE_END_C_DECLS
#endif

/*----------------------------------------------------------------------------
 * Macros for scoping of examples
 *----------------------------------------------------------------------------*/

#undef PLE_BEGIN_EXAMPLE_SCOPE
#undef PLE_END_EXAMPLE_SCOPE

#define PLE_BEGIN_EXAMPLE_SCOPE  {
#define   PLE_END_EXAMPLE_SCOPE  }

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int
(ple_printf_t) (const char  *const format,
                va_list            arg_ptr);

typedef void
(ple_error_handler_t) (const char  *file_name,
                       const int    line_num,
                       const int    sys_error_code,
                       const char  *format,
                       va_list      arg_ptr);

typedef void *
(ple_mem_malloc_t)(size_t       ni,
                   size_t       size,
                   const char  *var_name,
                   const char  *file_name,
                   int          line_num);

typedef void *
(ple_mem_realloc_t)(void        *ptr,
                    size_t       ni,
                    size_t       size,
                    const char  *var_name,
                    const char  *file_name,
                    int          line_num);

typedef void *
(ple_mem_free_t)(void        *ptr,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num);

/*============================================================================
 * Public function prototypes
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
           ...);

/* Returns function associated with the ple_printf() function.
 *
 * returns:
 *   pointer to the vprintf() or replacement function.
 */

ple_printf_t *
ple_printf_function_get(void);

/*
 * Associates a vprintf() type function with the ple_printf() function.
 *
 * parameters:
 *   f <-- pointer to a vprintf() type function.
 */

void
ple_printf_function_set(ple_printf_t  *f);

/*
 * Calls the error handler (set by ple_error_handler_set() or default).
 *
 * With the default error handler, an error message is output to stderr,
 * and the current process exits with an EXIT_FAILURE code.
 *
 * parameters:
 *   file_name      <-- name of source file from which error handler called.
 *   line_num       <-- line of source file from which error handler called.
 *   sys_error_code <-- error code if error in system or libc call,
 *                      0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   ...            <-- variable arguments based on format string.
 */

void
ple_error(const char  *file_name,
          const int    line_num,
          const int    sys_error_code,
          const char  *format,
          ...);

/*
 * Returns the error handler associated with the ple_error() function.
 *
 * returns:
 *   pointer to the error handler function.
 */

ple_error_handler_t *
ple_error_handler_get(void);

/*
 * Associates an error handler with the ple_error() function.
 *
 * parameters:
 *   handler <-- pointer to the error handler function.
 */

void
ple_error_handler_set(ple_error_handler_t  *handler);

/*
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * allocate the required memory.
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
 */

void *
ple_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, ple_alloc() called).
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num   -> line number in calling source file
 *
 * returns:
 *   pointer to allocated memory.
 */

void *
ple_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the ple_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, ple_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   NULL pointer.
 */

void *
ple_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num);

/* Return the function pointers associated with PLE's memory management.
 *
 * All arguments are optional.
 *
 * parameters:
 *   malloc_func  <-- pointer to ple_mem_malloc function pointer (or NULL).
 *   realloc_func <-- pointer to ple_mem_realloc function pointer (or NULL).
 *   free_func    <-- pointer to ple_mem_free function pointer (or NULL).
 */

void
ple_mem_functions_get(ple_mem_malloc_t   **malloc_func,
                      ple_mem_realloc_t  **realloc_func,
                      ple_mem_free_t     **free_func);

/* Associate functions to modifiy PLE's memory management.
 *
 * All arguments are optional, so the previously set functions pointers will
 * not be modified if an argument value is NULL.
 *
 * parameters:
 *   malloc_func  <-- ple_mem_malloc function pointer (or NULL).
 *   realloc_func <-- ple_mem_realloc function pointer (or NULL).
 *   free_func    <-- ple_mem_free function pointer (or NULL).
 */

void
ple_mem_functions_set(ple_mem_malloc_t   *malloc_func,
                      ple_mem_realloc_t  *realloc_func,
                      ple_mem_free_t     *free_func);

/*
 * Return Wall clock time
 *
 * returns:
 *   elapsed time from first call of a function of the ple_timer_...()
 *   series, or -1 if unable to compute.
 */

double
ple_timer_wtime(void);

/*
 * Return CPU time.
 *
 * Note that in the rare case that only the minimal C library clock()
 * method is available (see ple_timer_cpu_time_method()), at least one of
 * the ple_timer_...() functions (possibly this one) must be called
 * upon program start for this function to be used. In addition,
 * in this case, time may "loop" back to 0 every multiple of
 * 2^size_t / CLOCKS_PER_SEC seconds.
 *
 * returns:
 *   current CPU time usage, or -1 if unable to compute.
 */

double
ple_timer_cpu_time(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_DEFS_H__ */
