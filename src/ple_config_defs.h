#ifndef __PLE_CONFIG_DEFS_H__
#define __PLE_CONFIG_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
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
#include "ple_config_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Basic macros
 *============================================================================*/

/*============================================================================
 * Internationalization (future)
 *============================================================================*/

#if defined(ENABLE_NLS) && defined(HAVE_GETTEXT)
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
#endif /* ENABLE_NLS && HAVE_GETTEXT */

/*============================================================================
 * C99 Qualifiers
 *============================================================================*/

/* inline provided by ple_config.h if necessary */

/* restrict type qualifier (standard in C99) */

#if defined(__GNUC__)
#define restrict __restrict
#else
#define restrict
#endif

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_CONFIG_DEFS_H__ */
