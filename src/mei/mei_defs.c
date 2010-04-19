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

/*!
 * \file mei_defs.c
 *
 * \brief Configurable error and memory handling with default functions.
 */

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Optional library and MEI headers
 */

#include "mei_defs.h"

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

/*! \fn MEI_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls mei_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn MEI_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls mei_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn MEI_FREE(_ptr)
 * \brief Free allocated memory.
 *
 * This macro calls mei_mem_free(), automatically setting the
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
_mei_error_default(const char  *file_name,
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
 * This function simply wraps malloc() and calls _mei_error() if it fails
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
_mei_mem_malloc_default(size_t       ni,
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
        mei_error(file_name, line_num, errno,
                  _("Failure to allocate \"%s\" (%lu bytes)"),
                  var_name, (unsigned long)alloc_size);

    return p_ret;
}

/*
 * Default memory reallocation.
 *
 * This function simply wraps realloc() and calls _mei_error() if it fails
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
_mei_mem_realloc_default(void        *ptr,
                         size_t       ni,
                         size_t       size,
                         const char  *var_name,
                         const char  *file_name,
                         int          line_num)
{
    void  *p_ret;

    size_t realloc_size = ni * size;

    p_ret = realloc(ptr, realloc_size);

    if (size != 0 && p_ret == NULL)
        mei_error(file_name, line_num, errno,
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
_mei_mem_free_default(void        *ptr,
                      const char  *var_name,
                      const char  *file_name,
                      int          line_num)
{
    if (ptr != NULL)
        free(ptr);

    return NULL;
}

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static mei_error_handler_t *_mei_error = _mei_error_default;

static mei_mem_malloc_t *_mei_mem_malloc = _mei_mem_malloc_default;
static mei_mem_realloc_t *_mei_mem_realloc = _mei_mem_realloc_default;
static mei_mem_free_t *_mei_mem_free = _mei_mem_free_default;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Calls the error handler (set by mei_error_handler_set() or default).
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
mei_error(const char  *const file_name,
          const int          line_num,
          const int          sys_error_code,
          const char  *const format,
          ...)
{
    va_list  arg_ptr;

    va_start(arg_ptr, format);

    _mei_error(file_name, line_num, sys_error_code, format, arg_ptr);

    va_end(arg_ptr);
}

/*!
 * \brief Returns the error handler associated with the mei_error() function.
 *
 * \return pointer to the error handler function.
 */

mei_error_handler_t *
mei_error_handler_get(void)
{
    return _mei_error;
}

/*!
 * \brief Associates an error handler with the mei_error() function.
 *
 * \param [in] handler pointer to the error handler function.
 */

void
mei_error_handler_set(mei_error_handler_t  *handler)
{
    _mei_error = handler;
}

/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
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
mei_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num)
{
    return _mei_mem_malloc(ni, size, var_name, file_name, line_num);
}

/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, mei_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */

void *
mei_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
    return _mei_mem_realloc(ptr, ni, size, var_name, file_name, line_num);
}

/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the mei_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, mei_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns NULL pointer.
 */

void *
mei_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
    return _mei_mem_free(ptr, var_name, file_name, line_num);
}

/*!
 * \brief Return the function pointers associated with MEI's memory management.
 *
 * All arguments are optional.
 *
 * \param [out] malloc_func  pointer to mei_mem_malloc function pointer
 *                           (or NULL).
 * \param [out] realloc_func pointer to mei_mem_realloc function pointer
 *                           (or NULL).
 * \param [out] free_func    pointer to mei_mem_free function pointer
 *                           (or NULL).
 */

void
mei_mem_functions_get(mei_mem_malloc_t   **malloc_func,
                      mei_mem_realloc_t  **realloc_func,
                      mei_mem_free_t     **free_func)
{
    if (malloc_func != NULL)
        *malloc_func = _mei_mem_malloc;

    if (realloc_func != NULL)
        *realloc_func = _mei_mem_realloc;

    if (free_func != NULL)
        *free_func = _mei_mem_free;
}

/*!
 * \brief Associate functions to modifiy MEI's memory management.
 *
 * All arguments are optional, so the previously set functions pointers will
 * not be modified if an argument value is NULL.
 *
 * \param [in] malloc_func  mei_mem_malloc function pointer (or NULL).
 * \param [in] realloc_func mei_mem_realloc function pointer (or NULL).
 * \param [in] free_func    mei_mem_free function pointer (or NULL).
 */

void
mei_mem_functions_set(mei_mem_malloc_t   *malloc_func,
                      mei_mem_realloc_t  *realloc_func,
                      mei_mem_free_t     *free_func)
{
    if (malloc_func != NULL)
        _mei_mem_malloc = malloc_func;

    if (realloc_func != NULL)
        _mei_mem_realloc = realloc_func;

    if (free_func != NULL)
        _mei_mem_free = free_func;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
