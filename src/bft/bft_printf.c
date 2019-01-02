/*============================================================================
 * Base user-definable printf() wrapper or replacement.
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

/*----------------------------------------------------------------------------*/

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

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Additional Doxygen documentation
 *-----------------------------------------------------------------------------*/

/* Associated typedef documentation (for bft_printf.h) */

/*!
 * \typedef bft_printf_proxy_t
 *
 * \brief Function pointer for printf() type functions.
 *
 * \param [in] format       format string, as printf() and family.
 * \param [in, out] arg_ptr pointer to variable argument list based on format
 *                          string.
 */

/*!
 * \typedef bft_printf_flush_proxy_t
 *
 * \brief Function pointer for fflush(stdout) type functions.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default bft_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_bft_printf_flush_proxy_default(void);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static bft_printf_proxy_t        *_bft_printf_proxy = vprintf;
static bft_printf_flush_proxy_t  *_bft_printf_flush_proxy
                                    = _bft_printf_flush_proxy_default;

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Default bft_printf_flush() proxy.
 *
 * returns:
 *   return code of fflush(stdout).
 */

static int
_bft_printf_flush_proxy_default(void)
{
  return fflush(stdout);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Replacement for printf() with modifiable behavior.
 *
 * This function calls vprintf() by default, or a function with similar
 * arguments indicated by bft_printf_proxy_set().
 *
 * \param [in] format format string, as printf() and family.
 * \param [in] ...    variable arguments based on format string.
 *
 * \return number of characters printed, not counting the trailing '\\0'
 *         used to end output strings
 */

int
bft_printf(const char *const format,
           ...)
{
  int  retval;
  va_list  arg_ptr;

  va_start(arg_ptr, format);

  retval = _bft_printf_proxy(format, arg_ptr);

  va_end(arg_ptr);

  return retval;
}

/*!
 * \brief Flush for output of bft_printf() with modifiable behavior.
 *
 * This function calls fflush(stdout) if bft_printf()'s default behavior is
 * used. If bft_printf's behavior is modified with bft_printf_proxy_set(),
 * bft_printf_flush()'s behavior may have to be also adjusted with
 * bft_printf_flush_proxy_set().
 *
 * \return using the default behavior, the return value is that of
 *         fflush(stdout): O upon successful completion, EOF otherwise
 *         (with errno set to indicate the error).
 */

int
bft_printf_flush(void)
{
  return _bft_printf_flush_proxy();
}

/*!
 * \brief Returns function associated with the bft_printf() function.
 *
 * \return pointer to the vprintf() or replacement function.
 */

bft_printf_proxy_t *
bft_printf_proxy_get(void)
{
  return _bft_printf_proxy;
}

/*!
 * \brief Associates a vprintf() type function with the bft_printf() function.
 *
 * \param [in] fct pointer to a vprintf() type function.
 */

void
bft_printf_proxy_set(bft_printf_proxy_t *const fct)
{
  _bft_printf_proxy = fct;
}

/*!
 * \brief Returns function associated with bft_printf_flush().
 *
 * \return pointer to the bft_printf_flush() proxy.
 */

bft_printf_flush_proxy_t *
bft_printf_flush_proxy_get(void)
{
  return _bft_printf_flush_proxy;
}

/*!
 * \brief Associates a proxy function with bft_printf_flush().
 *
 * \warning
 *   bft_printf() is called by the default bft_error() error handler
 *   (so as to ensure that the error text appears at the end of the
 *   program output), so a bft_print_flush replacement must not itself
 *   call (directly or indirectly) bft_error() if the default error
 *   handler is used.
 *
 * \param [in] fct pointer to a function similar to {return fflush(stdout)}.
 */

void
bft_printf_flush_proxy_set(bft_printf_flush_proxy_t *const fct)
{
  _bft_printf_flush_proxy = fct;
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
