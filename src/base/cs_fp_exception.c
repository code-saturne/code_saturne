/*============================================================================
 * Floating-point exception handling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

/* On glibc-based systems, define _GNU_SOURCE so as to enable floating-point
   error exceptions; With Intel compilers, optimized code may raise such
   exceptions due to speculative execution, so we only enable raising of such
   exceptions for code compiled in debug mode, where reduced optimization
   should not lead to such exceptions, and locating the "true" origin of
   floating-point exceptions is helpful.
   _GNU_SOURCE must be defined before including any headers, to ensure
   the correct feature macros are defined first. */

#if defined(__linux__) || defined(__linux) || defined(linux)
#  define CS_FPE_TRAP
#  define _GNU_SOURCE
#endif

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(CS_FPE_TRAP)
#include <fenv.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_fp_exception.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Floating point exception handling */

#if defined(CS_FPE_TRAP)
static int    _fenv_set = 0;
static int    _fenv_save = 0;
static fenv_t _fenv_old;     /* Old exception mask */
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Enable floating-point exception trapping.
 *
 * Uses a counter to handle nested calls.
 */
/*----------------------------------------------------------------------------*/

void
cs_fp_exception_enable_trap(void)
{
#if defined(CS_FPE_TRAP)
  if (_fenv_set == 0) {
    if (fegetenv(&_fenv_old) == 0) {
      feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      _fenv_set = 1;
      /* To revert to initial behavior: fesetenv(&_fenv_old); */
    }
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Disable floating-point exception trapping.
 *
 * Uses a counter to handle nested calls.
 */
/*----------------------------------------------------------------------------*/

void
cs_fp_exception_disable_trap(void)
{
#if defined(CS_FPE_TRAP)
  if (_fenv_save == 0) {
    if (fegetenv(&_fenv_old) == 0) {
      _fenv_save += 1;
      fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    }
  }
  else
    _fenv_save += 1;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Restore floating-point exception trapping.
 */
/*----------------------------------------------------------------------------*/

void
cs_fp_exception_restore_trap(void)
{
#if defined(CS_FPE_TRAP)
  if (_fenv_save) {
    _fenv_save -= 1;
    if (_fenv_save == 0)
      fesetenv(&_fenv_old);
  }
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
