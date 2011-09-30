#ifndef __CS_DEFS_H__
#define __CS_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Autoconf-defined macros
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#  include "cs_config.h"
#endif

/*============================================================================
 * Internationalization
 *============================================================================*/

#if defined(ENABLE_NLS) && defined(HAVE_GETTEXT)

#  include <libintl.h>
#  define _(String) dgettext(PACKAGE, String)
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
#  define bindtextdomain(Domain, Directory) (Domain)

#endif /* ENABLE_NLS && HAVE_GETTEXT */

/*============================================================================
 * Parallelism
 *============================================================================*/

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)
#  include <mpi.h>
#endif

#if defined(HAVE_OPENMP)
#  include <omp.h>
#endif

/*============================================================================
 * C99 Qualifiers
 *============================================================================*/

/* inline provided by cs_config.h if necessary */

#if !defined(__STDC_VERSION__)
#  define __STDC_VERSION__ 1989
#endif

/*
 * Redefinition of "inline" et "restrict" qualifiers incompatible with
 * some C89 compilers (standard in C99)
 */

#if (__STDC_VERSION__ < 199901L)

#  if defined(__GNUC__)
#    define inline __inline__
#    define restrict __restrict__
#  else
#    define inline
#    define restrict
#  endif

#endif

/*============================================================================
 * Definitions that may not always be provided directly by the system
 *============================================================================*/

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if HAVE_STDDEF_H
#  include <stddef.h>
#else
#  include <stdlib.h>
#endif

/*
 * Usually stdint.h is included by inttypes.h, but only inttypes.h exists
 * on certain systems, such as Tru64Unix.
 */

#if HAVE_STDINT_H
#  include <stdint.h>
#elif HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

/*
 * Obtain the definition of off_t.
 */

#if defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#endif

/* C99 _Bool type */

#if HAVE_STDBOOL_H
#  include <stdbool.h>
#else
#  ifndef HAVE__BOOL
#    ifdef __cplusplus
typedef bool _Bool;
#    else
#      define _Bool signed char;
#    endif
#  endif
#  define bool _Bool
#  define false 0
#  define true 1
#  define __bool_true_false_are_defined 1
#endif

/* int32_t type */

#if !defined(HAVE_INT32_T)
#  if (SIZEOF_INT == 4)
typedef int int32_t;
#  elif (SIZEOF_SHORT == 4)
typedef short int32_t;
#  else
#    error
#  endif
#endif

/* int64_t type */

#if !defined(HAVE_INT64_T)
#  if (SIZEOF_INT == 8)
typedef int int64_t;
#  elif (SIZEOF_LONG == 8)
typedef long int64_t;
#  elif (HAVE_LONG_LONG == 8)  /* SIZEOF_LONG_LONG not generally available */
typedef long long int64_t;
#  else
#    error
#  endif
#endif

/* uint32_t type */

#if !defined(HAVE_UINT32_T)
#  if (SIZEOF_INT == 4)
typedef unsigned uint32_t;
#  elif (SIZEOF_SHORT == 4)
typedef unsigned short uint32_t;
#  else
#    error
#  endif
#endif

/* uint64_t type */

#if !defined(HAVE_UINT64_T)
#  if (SIZEOF_INT == 8)
typedef unsigned uint64_t;
#  elif (SIZEOF_LONG == 8)
typedef unsigned long uint64_t;
#  elif (HAVE_LONG_LONG) /* SIZEOF_LONG_LONG not generally available */
typedef unsigned long long uint64_t;
#  else
#    error
#  endif
#endif

/*============================================================================
 * General types and macros used throughout Code_Saturne
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Basic types used by Code_Saturne
 * They may be modified here to better map to a given library, with the
 * following constraints:
 *  - cs_lnum_t must be signed
 *  - cs_gnum_t may be signed or unsigned
 *----------------------------------------------------------------------------*/

/* Global integer index or number */

#if defined(HAVE_LONG_GNUM)
  #if (SIZEOF_LONG == 8)
    typedef unsigned long       cs_gnum_t;
  #elif (SIZEOF_LONG_LONG == 8)
    typedef unsigned long long  cs_gnum_t;
  #else
    #error
  #endif
#else
  typedef unsigned  cs_gnum_t;
#endif

/* Other types */

typedef int      cs_lnum_t;     /* Local integer index or number */
typedef double   cs_coord_t;    /* Real number (coordinate value) */

typedef int      cs_int_t;      /* Fortran integer */
typedef double   cs_real_t;     /* Fortran double precision */
typedef char     cs_byte_t;     /* Byte (untyped memory unit) */

/* Mappings to MPI datatypes */
/*---------------------------*/

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

#  define CS_MPI_INT       MPI_INT         /* If cs_int_t is an int */
#  define CS_MPI_REAL      MPI_DOUBLE      /* If cs_real_t is a double */

/* MPI type for cs_gnum_t integer type (depends on configuration) */

#  if defined(HAVE_LONG_GNUM)
#    if (SIZEOF_LONG == 8)
#      define CS_MPI_GNUM     MPI_UNSIGNED_LONG
#    elif (SIZEOF_LONG_LONG == 8)
#      if defined(MPI_UNSIGNED_LONG_LONG)
#        define CS_MPI_GNUM     MPI_UNSIGNED_LONG_LONG
#      elif defined(MPI_LONG_LONG)
#        define CS_MPI_GNUM     MPI_LONG_LONG
#      endif
#    endif
#    if !defined(CS_MPI_GNUM)
#      error
#    endif
#  else
#    define CS_MPI_GNUM       MPI_UNSIGNED
#  endif

#  define CS_MPI_LNUM     MPI_INT         /* MPI type for cs_lnum_t type */
#  define CS_MPI_COORD    MPI_DOUBLE      /* MPI type for cs_coord_t type */

#endif /* defined(HAVE_MPI) && !defined(CS_IGNORE_MPI) */

/* Macros for compilation with a C++ compiler */
/*--------------------------------------------*/

#undef BEGIN_C_DECLS
#undef   END_C_DECLS

#if defined(__cplusplus)
#  define BEGIN_C_DECLS  extern "C" {
#  define   END_C_DECLS  }
#else
#  define BEGIN_C_DECLS
#  define   END_C_DECLS
#endif

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_DEFS_H__ */
