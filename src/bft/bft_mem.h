#ifndef __BFT_MEM_H__
#define __BFT_MEM_H__

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

/*============================================================================
 * Public macros
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls bft_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define BFT_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) bft_mem_malloc(_ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls bft_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define BFT_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) bft_mem_realloc(_ptr, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls bft_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#ifdef __cplusplus /* avoid casting from void for C++ */

#define BFT_FREE(_ptr) \
bft_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

#else

#define BFT_FREE(_ptr) \
_ptr = bft_mem_free(_ptr, #_ptr, __FILE__, __LINE__)

#endif /* __cplusplus */

/*
 * Allocate aligned memory for _ni items of type _type.
 *
 * This macro calls bft_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr    --> pointer to allocated memory.
 *   _align <-- alignment.
 *   _ni    <-- number of items.
 *   _type  <-- element type.
 */

#define BFT_MEMALIGN(_ptr, _align, _ni, _type) \
_ptr = (_type *) bft_mem_memalign(_align, _ni, sizeof(_type), \
                                  #_ptr, __FILE__, __LINE__)

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef size_t
(bft_mem_get_size_t)(void  *ptr);

typedef void *
(bft_mem_realloc_t)(void        *ptr,
                    size_t       ni,
                    size_t       size,
                    const char  *var_name,
                    const char  *file_name,
                    int          line_num);

typedef void
(bft_mem_free_t)(void        *ptr,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num);

/*============================================================================
 * Semi private function definitions
 *============================================================================*/

END_C_DECLS

#ifdef __cplusplus /* only for C++ API */

/*----------------------------------------------------------------------------*/

/*
 * Return the cs_mem_block structure corresponding to a given
 * allocated block.
 *
 * parameters:
 *   p_get: <-- allocated block's start adress.
 *
 * returns:
 *   corresponding cs_mem_block structure.
 */

cs_mem_block_t
bft_mem_get_block_info(const void  *p_get);

/*
 * Return the cs_mem_block structure corresponding to a given
 * allocated block if available.
 *
 * If no block info is available, return block with null pointers
 * and zero size.
 *
 * parameters:
 *   p_get: <-- allocated block's start adress.
 *
 * returns:
 *   pointer tocorresponding cs_mem_block structure.
 */

cs_mem_block_t
bft_mem_get_block_info_try(const void  *p_get);

/*
 * \brief Update block information map if enabled.
 *
 * parameters
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *   old_block <-- pointer to old block info, if present
 *   new_block <-- pointer to new block info, if present
 */

void
bft_mem_update_block_info(const char            *var_name,
                          const char            *file_name,
                          int                    line_num,
                          const cs_mem_block_t  *old_block,
                          const cs_mem_block_t  *new_block);

/*----------------------------------------------------------------------------*/

#endif // __cplusplus

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Initialize memory handling.
 *
 * This function should be called before any other bft_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

void
bft_mem_init(const char  *log_file_name);

/*
 * End memory handling.
 *
 * This function should be called after all other bft_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void
bft_mem_end(void);

/*
 * Initialize memory handling.
 *
 * This function should be called before any other bft_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

/*
 * Indicates if bft_mem_...() functions are initialized.
 *
 * returns:
 *   1 if bft_mem_init has been called, 0 otherwise.
 */

int
bft_mem_initialized(void);

/*
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * Allocation couting and logging to trace file will be done if
 * both required by the bft_mem_init options and if file_name != nullptr.
 * If required but file_name == nullptr, it must be handled by the caller,
 * using `bft_mem_log_mem_op`.
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
bft_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, bft_alloc() called).
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
bft_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the bft_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, bft_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   NULL pointer.
 */

void *
bft_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num);

/*
 * Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the bft_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function bft_mem_have_memalign() indicates if this
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
 */

void *
bft_mem_memalign(size_t       alignment,
                 size_t       ni,
                 size_t       size,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num);

/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through bft_mem_...() (in kB).
 */

size_t
bft_mem_size_current(void);

/*!
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through bft_mem_...() (in kB).
 */

size_t
bft_mem_size_max(void);

/*
 * Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, bft_mem_memalign() will always fail.
 *
 * returns:
 *   1 if memory aligned allocation is possible, 0 otherwise.
 */

int
bft_mem_have_memalign(void);

/* Returns the error handler associated with the bft_mem_...() functions.
 *
 * returns:
 *   pointer to the error handler function.
 */

bft_error_handler_t *
bft_mem_error_handler_get(void);

/*
 * Associates an error handler with the bft_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bft_print_flush() is called), and the general error handler used
 * by bft_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameter:
 *   handler <-- pointer to the error handler function.
 */

void
bft_mem_error_handler_set(bft_error_handler_t *handler);

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
                        bft_mem_free_t      *free_func);

/*----------------------------------------------------------------------------*/
/*
 * Set the level of memory log output when a log file is used.
 */
/*----------------------------------------------------------------------------*/

void
bft_mem_log_level_set
(
  int level /*<[in] memory log level (int). 0 for light and 1 for full. */
);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __BFT_MEM_H__ */
