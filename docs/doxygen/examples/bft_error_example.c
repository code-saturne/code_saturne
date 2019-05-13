/*============================================================================
 * BFT Examples
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*
 * Standard C library and BFT headers
 */

/*! [my_error_handler_headers] */
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
# include <mpi.h>
#endif

#include "bft_intl.h"
/*! [my_error_handler_headers] */
#include "bft_error.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*
 * MPI-aware error handler.
 *
 * An error message is output to stderr (after bft_print_flush() is called),
 * and the current process (or process group) is terminated.
 *
 * parameters:
 *   file_name:      --> name of source file from which error handler called.
 *   line_num:       --> line of source file from which error handler called.
 *   sys_error_code: --> error code if error in system or libc call, 0 otherwise.
 *   format:         --> format string, as printf() and family.
 *   arg_ptr:        --> variable argument list based on format string.
 */

/*! [my_error_handler_body] */
void
my_error_handler(const char     *const file_name,
                 const int             line_num,
                 const int             sys_error_code,
                 const char     *const format,
                 const va_list         arg_ptr)
{
  bft_printf_flush();

  fprintf(stderr, "\n");

  if (sys_error_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_error_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), file_name, line_num);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  assert(0);   /* Use assert to avoit exiting under debugger */

#if defined(HAVE_MPI)
  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
#endif /* HAVE_MPI */

  exit(EXIT_FAILURE);
}
/*! [my_error_handler_body] */

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
