/*============================================================================
 * Base error handling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*-----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*-----------------------------------------------------------------------------*/

/*
 * Standard C library headers
 */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Optional library and BFT headers
 */

#include "bft_printf.h"
#include "bft_error.h"

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Additional doxygen documentation
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for bft_error.h) */

/*!
 * \typedef bft_error_handler_t
 *
 * \brief Function pointer to opaque error handler.
 *
 * \param [in] file_name      name of source file from which error handler
 *                            called.
 * \param [in] line_num       line of source file from which error handler
 *                            called.
 * \param [in] sys_error_code error code if error in system or libc call,
 *                            0 otherwise.
 * \param [in] format         format string, as printf() and family.
 * \param [in, out] arg_ptr   pointer to variable argument list based on format
 *                            string.
 *
 * In an MPI environment, it is recommended to replace the default
 * error handler. This requires using the following headers:
 *
 * \snippet bft_error_example.c my_error_handler_headers
 *
 * An error handler function similar to the BFT default with MPI-awareness
 * added looks like:
 *
 * \snippet bft_error_example.c my_error_handler_body
 *
 * In a more complex environment, \c MPI_COMM_WORLD could be replaced
 * by another communicator.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default error handler.
 *
 * An error message is output to stderr (stdout is flushed first),
 * and the current process is terminated.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bft_error_handler_default(const char     *const file_name,
                           const int             line_num,
                           const int             sys_error_code,
                           const char     *const format,
                           va_list               arg_ptr);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static bft_error_handler_t  *_bft_error_handler = (_bft_error_handler_default);

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default error handler.
 *
 * An error message is output to stderr (after bft_print_flush() is called),
 * and the current process exits with an EXIT_FAILURE code.
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bft_error_handler_default(const char     *const file_name,
                           const int             line_num,
                           const int             sys_error_code,
                           const char     *const format,
                           va_list               arg_ptr)
{
  bft_printf_flush();

  fprintf(stderr, "\n");

  if (sys_error_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_error_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  assert(0);

  exit(EXIT_FAILURE);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Calls the error handler (set by bft_error_handler_set() or default).
 *
 * With the default error handler, bft_print_flush() is called, an error
 * message is output to stderr, and the current process exits with an
 * EXIT_FAILURE code.
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
bft_error(const char  *const file_name,
          const int          line_num,
          const int          sys_error_code,
          const char  *const format,
          ...)
{
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  _bft_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*!
 * \brief Returns the error handler associated with the bft_error() function.
 *
 * \return pointer to the error handler function.
 */

bft_error_handler_t *
bft_error_handler_get(void)
{
  return _bft_error_handler;
}

/*!
 * \brief Associates an error handler with the bft_error() function.
 *
 * \param [in] handler pointer to the error handler function.
 */

void
bft_error_handler_set(bft_error_handler_t  *const handler)
{
  _bft_error_handler = handler;
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
