#ifndef __PLE_CONFIG_DEFS_H__
#define __PLE_CONFIG_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2020  EDF S.A.

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

/*============================================================================
 * Definitions that may not always be provided directly by the system
 *============================================================================*/

/*
 * Usually stdint.h is included by inttypes.h, but only inttypes.h exists
 * on certain systems, such as Tru64Unix.
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* C99 _Bool type */

#if HAVE_STDBOOL_H
# include <stdbool.h>
#else
# if !HAVE__BOOL
#  ifdef __cplusplus
typedef bool _Bool;
#  else
    define _Bool signed char;
#  endif
# endif
# define bool _Bool
# define false 0
# define true 1
# define __bool_true_false_are_defined 1
#endif

/* int32_t type */

#if !defined(HAVE_INT32_T)
# if (PLE_SIZEOF_INT == 4)
typedef int int32_t;
# elif (PLE_SIZEOF_SHORT == 4)
typedef short int32_t;
# else
#  error
# endif
#endif

/* int64_t type */

#if !defined(HAVE_INT64_T)
# if (PLE_SIZEOF_INT == 8)
typedef int int64_t;
# elif (PLE_SIZEOF_LONG == 8)
typedef long int64_t;
# elif (HAVE_LONG_LONG == 8)  /* PLE_SIZEOF_LONG_LONG not generally available */
typedef long long int64_t;
# else
#  error
# endif
#endif

/* uint32_t type */

#if !defined(HAVE_UINT32_T)
# if (PLE_SIZEOF_INT == 4)
typedef unsigned uint32_t;
# elif (PLE_SIZEOF_SHORT == 4)
typedef unsigned short uint32_t;
# else
#  error
# endif
#endif

/* uint64_t type */

#if !defined(HAVE_UINT64_T)
# if (PLE_SIZEOF_INT == 8)
typedef unsigned uint64_t;
# elif (PLE_SIZEOF_LONG == 8)
typedef unsigned long uint64_t;
# elif (HAVE_LONG_LONG) /* PLE_SIZEOF_LONG_LONG not generally available */
typedef unsigned long long uint64_t;
# else
#  error
# endif
#endif

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_CONFIG_DEFS_H__ */
