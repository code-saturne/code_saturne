#ifndef _ECS_DEF_H_
#define _ECS_DEF_H_

/*============================================================================
 * Definitions, global variables, and base functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_config.h"

/*
 * Standard C library headers
 */

/*============================================================================
 * Definition of C99 type which may not be available with older tools.
 *============================================================================*/

#if HAVE_STDDEF_H
# include <stddef.h>
#else
# include <stdio.h>
#endif

/*
 * Usually, stdint.h is included by inttypes.h, but only inttypes.h
 * may be found on some systems, such as Tru64Unix.
 */

#if HAVE_STDINT_H
# include <stdint.h>
#elif HAVE_INTTYPES_H
# include <inttypes.h>
#endif

/* _Bool */

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

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

#if defined(HAVE_LONG_GNUM)
  #if (SIZEOF_LONG == 8)
    typedef long                ecs_int_t;   /* Integer */
    typedef unsigned long       ecs_size_t;  /* Index size */
  #elif (SIZEOF_LONG_LONG == 8)
    typedef long long           ecs_int_t;   /* Integer */
    typedef unsigned long long  ecs_size_t;  /* Index size */
  #else
    #error
  #endif
#else
  typedef int           ecs_int_t;      /* Integer */
  typedef size_t        ecs_size_t;     /* Index size */
#endif

typedef double          ecs_coord_t;    /* Real (floating point) */
typedef char            ecs_byte_t;     /* Byte (untyped memory) */

/* Type enumeration */

typedef enum {
  ECS_TYPE_char,
  ECS_TYPE_ecs_coord_t,
  ECS_TYPE_ecs_int_t,
  ECS_TYPE_ecs_size_t,
  ECS_TYPE_size_t,
  ECS_TYPE_void
} ecs_type_t;

/* Element types */

typedef enum {

  ECS_ELT_TYP_NUL,         /*  No type */

  ECS_ELT_TYP_FAC_TRIA,    /*  Triangle */
  ECS_ELT_TYP_FAC_QUAD,    /*  Quadrangle */

  ECS_ELT_TYP_CEL_TETRA,   /*  Tetrahedron */
  ECS_ELT_TYP_CEL_PYRAM,   /*  Pyramid */
  ECS_ELT_TYP_CEL_PRISM,   /*  Prism */
  ECS_ELT_TYP_CEL_HEXA,    /*  Hexahedron */

  ECS_ELT_TYP_FAC_POLY,    /*  Polygon */
  ECS_ELT_TYP_CEL_POLY,    /*  Polyhedron */

  ECS_ELT_TYP_FIN

} ecs_elt_typ_t ;

/* Restrict qualifier (standard in C99) */

#if defined(__GNUC__)
#define restrict __restrict
#else
#define restrict
#endif

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Classical macros */

#define ECS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value */
#define ECS_MIN(a, b)  ((a) > (b) ?  (b) : (a))  /* Minimum */
#define ECS_MAX(a, b)  ((a) < (b) ?  (b) : (a))  /* Maximum */

/* Directory separator character */

#define ECS_PATH_SEP             '/'

#define ECS_REAL_PRECISION    1.e-13

#define ECS_STR_SIZE       80

#define ECS_PAS_NUL         0
#define ECS_PAS_UNITE       1

#define ECS_LNG_AFF_STR           43
#define ECS_LNG_AFF_ENT            8
#define ECS_LNG_AFF_REE_MANTIS    11
#define ECS_LNG_AFF_REE_PRECIS     2

#define ECS_FMT_AFF_REE_PARAM     "%.15E"

/*
 * Internationalization macros.
 */

#if defined(ENABLE_NLS)

#include <libintl.h>
#define _(String) gettext(String)
#define gettext_noop(String) String
#define N_(String) gettext_noop(String)

#else

#define _(String) String
#define N_(String) String
#define textdomain(Domain)
#define bindtextdomain(Package, Directory)

#endif

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

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

extern char      ecs_glob_build_date[];

extern int       ecs_glob_have_cgns;    /* CGNS library support*/
extern int       ecs_glob_cgns_ver_maj;
extern int       ecs_glob_cgns_ver_min;
extern int       ecs_glob_cgns_ver_rel;

extern int       ecs_glob_have_med;     /* MED library support */

/* Element type determination based on an element's number of vertices,
   for elements of dimension 2 (faces) or 3 (cells);
   Beyond 8 vertices, we always have ECS_ELT_TYP_FAC_POLY for faces
   and ECS_ELT_TYP_FAC_POLY for cells */

extern  const ecs_elt_typ_t  ecs_glob_typ_elt[2][9];

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize error handling
 *----------------------------------------------------------------------------*/

void
ecs_init_gestion_erreur(void);

/*----------------------------------------------------------------------------
 * Exit function.
 *----------------------------------------------------------------------------*/

void
ecs_exit(int  statut);

/*----------------------------------------------------------------------------
 * Print a warning.
 *----------------------------------------------------------------------------*/

void
ecs_warn(void);

/*----------------------------------------------------------------------------
 * Abort on error.
 *----------------------------------------------------------------------------*/

void
ecs_error(const char  *file_name,
          const int    line_num,
          const int    sys_error_code,
          const char  *format,
          ...);

/*----------------------------------------------------------------------------
 * Print a string with a given column width.
 *----------------------------------------------------------------------------*/

void
ecs_print_padded_str(const char  *str,
                     int          width);

/*----------------------------------------------------------------------------*/

#endif /* _ECS_DEF_H_ */
