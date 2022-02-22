#ifndef __ECS_MEM_H__
#define __ECS_MEM_H__

/*============================================================================
 * Base memory allocation wrappers with optional tracing
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

#include "ecs_def.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public macros
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls ecs_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define ECS_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) ecs_mem_malloc(_ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls ecs_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define ECS_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) ecs_mem_realloc(_ptr, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls ecs_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#ifdef __cplusplus /* avoid casting from void for C++ */

#define ECS_FREE(_ptr) \
ecs_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

#else

#define ECS_FREE(_ptr) \
_ptr = ecs_mem_free(_ptr, #_ptr, __FILE__, __LINE__)

#endif /* __cplusplus */

/*
 * Allocate aligned memory for _ni items of type _type.
 *
 * This macro calls ecs_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr    --> pointer to allocated memory.
 *   _align <-- alignment.
 *   _ni    <-- number of items.
 *   _type  <-- element type.
 */

#define ECS_MEMALIGN(_ptr, _align, _ni, _type) \
_ptr = (_type *) ecs_mem_memalign(_align, _ni, sizeof(_type), \
                                  #_ptr, __FILE__, __LINE__)

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Initialize memory handling.
 *
 * This function should be called before any other ecs_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

void
ecs_mem_init(const char  *log_file_name);

/*
 * End memory handling.
 *
 * This function should be called after all other ecs_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void
ecs_mem_end(void);

/*
 * Initialize memory handling.
 *
 * This function should be called before any other ecs_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

/*
 * Indicates if ecs_mem_...() functions are initialized.
 *
 * returns:
 *   1 if ecs_mem_init has been called, 0 otherwise.
 */

int
ecs_mem_initialized(void);

/*
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
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
ecs_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, ecs_alloc() called).
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
ecs_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, ecs_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   NULL pointer.
 */

void *
ecs_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num);

/*
 * Return current theoretical dynamic memory allocated.
 *
 * returns:
 *   current memory handled through ecs_mem_...() (in kB).
 */

size_t
ecs_mem_size_current(void);

/*
 * Return maximum theoretical dynamic memory allocated.
 *
 * returns:
 *   maximum memory handled through ecs_mem_...() (in kB).
 */

size_t
ecs_mem_size_max(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __ECS_MEM_H__ */
