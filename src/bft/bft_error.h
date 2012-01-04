#ifndef __BFT_ERROR_H__
#define __BFT_ERROR_H__

/*============================================================================
 * Base error handling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/* Standard C library headers */

#include <stdarg.h>

/* BFT library headers */

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

typedef void (bft_error_handler_t) (const char    *const file_name,
                                    const int            line_num,
                                    const int            sys_error_code,
                                    const char    *const format,
                                    va_list              arg_ptr);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Calls the error handler (set by bft_error_handler_set() or default).
 *
 * With the default error handler, bft_print_flush() is called, an error
 * message is output to stderr, and the current process exits with an
 * EXIT_FAILURE code.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   ... :           <-- variable arguments based on format string.
 */

void
bft_error(const char  *const file_name,
          const int          line_num,
          const int          sys_error_code,
          const char  *const format,
          ...);

/*
 * Returns the error handler associated with the bft_error() function.
 *
 * returns:
 *  pointer to the error handler function.
 */

bft_error_handler_t *
bft_error_handler_get(void);

/*
 * Associates an error handler with the bft_error() function.
 *
 * parameters:
 *   handler: <-- pointer to the error handler function.
 */

void
bft_error_handler_set(bft_error_handler_t  *const handler);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFT_ERROR_H__ */
