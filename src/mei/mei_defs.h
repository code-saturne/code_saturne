#ifndef __MEI_MEM_H__
#define __MEI_MEM_H__

/*============================================================================
 * Base definitions and utility functions for MEI.
 *============================================================================*/

/*
  This file is part of the "Mathematical Expression Interpreter" library.

  Copyright (C) 2009  EDF

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

#include <stdarg.h>

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if defined(__STDC_VERSION__)
#    if (__STDC_VERSION__ == 199901L)
#        include <stddef.h>
#    else
#        include <stdlib.h>
#    endif
#else
#    include <stdlib.h>
#endif

/* MEI headers */

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int (mei_printf_t) (const char     *const format,
                            va_list               arg_ptr);

typedef void
(mei_error_handler_t) (const char    *file_name,
                       const int      line_num,
                       const int      sys_error_code,
                       const char    *format,
                       va_list        arg_ptr);

typedef void *
(mei_mem_malloc_t)(size_t       ni,
                   size_t       size,
                   const char  *var_name,
                   const char  *file_name,
                   int          line_num);

typedef void *
(mei_mem_realloc_t)(void        *ptr,
                    size_t       ni,
                    size_t       size,
                    const char  *var_name,
                    const char  *file_name,
                    int          line_num);

typedef void *
(mei_mem_free_t)(void        *ptr,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num);

/*============================================================================
 * Public macros
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls mei_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define MEI_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) mei_mem_malloc(_ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls mei_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define MEI_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) mei_mem_realloc(_ptr, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls mei_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#ifdef __cplusplus /* avoid casting from void for C++ */

#define MEI_FREE(_ptr) \
mei_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

#else

#define MEI_FREE(_ptr) \
_ptr = mei_mem_free(_ptr, #_ptr, __FILE__, __LINE__)

#endif /* __cplusplus */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/* Returns function associated with the mei_printf() function.
 *
 * returns:
 *   pointer to the vprintf() or replacement function.
 */

mei_printf_t *
mei_printf_function_get(void);

/*
 * Associates a vprintf() type function with the mei_printf() function.
 *
 * parameters:
 *   f <-- pointer to a vprintf() type function.
 */

void
mei_printf_function_set(mei_printf_t  *f);

/*
 * Calls the error handler (set by mei_error_handler_set() or default).
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
mei_error(const char  *file_name,
          const int    line_num,
          const int    sys_error_code,
          const char  *format,
          ...);

/*
 * Returns the error handler associated with the mei_error() function.
 *
 * returns:
 *   pointer to the error handler function.
 */

mei_error_handler_t *
mei_error_handler_get(void);

/*
 * Associates an error handler with the mei_error() function.
 *
 * parameters:
 *   handler <-- pointer to the error handler function.
 */

void
mei_error_handler_set(mei_error_handler_t  *handler);

/*
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
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
mei_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, mei_alloc() called).
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
mei_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, mei_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   NULL pointer.
 */

void *
mei_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num);

/* Return the function pointers associated with MEI's memory management.
 *
 * All arguments are optional.
 *
 * parameters:
 *   malloc_func  <-- pointer to mei_mem_malloc function pointer (or NULL).
 *   realloc_func <-- pointer to mei_mem_realloc function pointer (or NULL).
 *   free_func    <-- pointer to mei_mem_free function pointer (or NULL).
 */

void
mei_mem_functions_get(mei_mem_malloc_t   **malloc_func,
                      mei_mem_realloc_t  **realloc_func,
                      mei_mem_free_t     **free_func);

/*Associate functions to modifiy MEI's memory management.
 *
 * All arguments are optional, so the previously set functions pointers will
 * not be modified if an argument value is NULL.
 *
 * parameters:
 *   malloc_func  <-- mei_mem_malloc function pointer (or NULL).
 *   realloc_func <-- mei_mem_realloc function pointer (or NULL).
 *   free_func    <-- mei_mem_free function pointer (or NULL).
 */

void
mei_mem_functions_set(mei_mem_malloc_t   *malloc_func,
                      mei_mem_realloc_t  *realloc_func,
                      mei_mem_free_t     *free_func);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MEI_MEM_H__ */
