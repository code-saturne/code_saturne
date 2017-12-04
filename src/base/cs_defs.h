#ifndef __CS_DEFS_H__
#define __CS_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*============================================================================
 * Autoconf-defined macros
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#  include "cs_config.h"
#endif

/*============================================================================
 * Internationalization
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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

#ifdef __cplusplus
}
#endif /* __cplusplus */

/*============================================================================
 * Parallelism
 *============================================================================*/

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

#  include <mpi.h>

#  if !defined(MPI_VERSION) /* Defined in up-to-date MPI versions */
#    define MPI_VERSION 1
#  endif

#  if MPI_VERSION == 1
#    define MPI_Info       int
#    define MPI_INFO_NULL  0
#  endif

#endif

#if defined(HAVE_OPENMP)
#  include <omp.h>
#endif

/*============================================================================
 * C99 Qualifiers
 *============================================================================*/

#ifndef __cplusplus /* C */

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

#else /* C++ */

#  ifndef HAVE_RESTRICT /* Must be provided by caller */
#    define restrict
#  endif

#endif /* __cplusplus */

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
#  ifndef __cplusplus
#    ifndef HAVE__BOOL
#      define _Bool signed char;
#    endif
#    define bool _Bool
#    define false 0
#    define true 1
#  else
#    define _Bool bool;
#  endif
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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * Variable value type.
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_DATATYPE_NULL,      /* empty datatype */
  CS_CHAR,               /* character values */
  CS_FLOAT,              /* 4-byte floating point values */
  CS_DOUBLE,             /* 8-byte floating point values */
  CS_UINT16,             /* 2-byte unsigned integer values */
  CS_INT32,              /* 4-byte signed integer values */
  CS_INT64,              /* 8-byte signed integer values */
  CS_UINT32,             /* 4-byte unsigned integer values */
  CS_UINT64              /* 8-byte unsigned integer values */

} cs_datatype_t;

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

typedef int                 cs_lnum_t;    /* Local integer index or number */
typedef double              cs_coord_t;   /* Real number (coordinate value) */

typedef int                 cs_int_t;    /* Fortran integer */
typedef double              cs_real_t;   /* Fortran double precision */
typedef char                cs_byte_t;   /* Byte (untyped memory unit) */
typedef unsigned short int  cs_flag_t;   /* Flag for storing metadata */

/* Vector or array block types */

typedef int        cs_lnum_2_t[2];      /* Vector of 2 local numbers */

typedef double     cs_coord_3_t[3];     /* Vector of 3 real (coordinate)
                                           values */

typedef cs_real_t  cs_real_2_t[2];      /* Vector of 2 real values */
typedef cs_real_t  cs_real_3_t[3];      /* Vector of 3 real values */
typedef cs_real_t  cs_real_4_t[4];      /* Vector of 4 real values */
typedef cs_real_t  cs_real_6_t[6];      /* Vector of 6 real values
                                           (for symmetric tensor) */
typedef cs_real_t  cs_real_9_t[9];      /* Vector of 9 real values */

typedef cs_real_t  cs_real_33_t[3][3];  /* Matrix of 3x3 real values */
typedef cs_real_t  cs_real_66_t[6][6];  /* Matrix of 6x6 real values */
typedef cs_real_t  cs_real_99_t[9][9];  /* Matrix of 9x9 real values */

typedef cs_real_t  cs_real_34_t[3][4];  /* Matrix of 3x4 real values */

typedef cs_real_t  cs_real_63_t[6][3];  /* Matrix of 6x3 real values */

typedef cs_real_33_t  cs_real_332_t[2];  /* vector of 2 3x3 matrices
                                            of real values */
typedef cs_real_66_t  cs_real_662_t[2];  /* vector of 2 6x6 matrices
                                            of real values */

typedef struct {

  double   val;  /* Value */
  int      id;   /* Id related to value */

} cs_double_int_t;

/* Vector-valued quantity stored using its measure (i.e. length) and
   its direction given by a unitary vector */
typedef struct {

  double  meas;
  double  unitv[3];

} cs_nvec3_t;

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

#  define CS_MPI_FLAG     MPI_UNSIGNED_SHORT /* MPI type for cs_flag_t type */
#  define CS_MPI_LNUM     MPI_INT            /* MPI type for cs_lnum_t type */
#  define CS_MPI_COORD    MPI_DOUBLE         /* MPI type for cs_coord_t type */

#endif /* defined(HAVE_MPI) && !defined(CS_IGNORE_MPI) */

/* Mappings to Code_Saturne datatypes */
/*------------------------------------*/

#if defined(HAVE_LONG_GNUM)
# define CS_GNUM_TYPE     CS_UINT64
#elif (SIZEOF_INT == 8)
# define CS_GNUM_TYPE     CS_UINT64
#else
# define CS_GNUM_TYPE     CS_UINT32
#endif

#if (SIZEOF_INT == 8)
# define CS_LNUM_TYPE     CS_INT64
#else
# define CS_LNUM_TYPE     CS_INT32
#endif

#if (SIZEOF_INT == 8)
# define CS_INT_TYPE      CS_INT64
#else
# define CS_INT_TYPE      CS_INT32
#endif

#define CS_FLAG_TYPE      CS_UINT16
#define CS_REAL_TYPE      CS_DOUBLE
#define CS_COORD_TYPE     CS_DOUBLE

/* Minimum size for OpenMP loops
 *  (based on initial benchmarking, which will need updates)
 *-----------------------------------------------------------*/

#if defined(__bgq__) && defined(__xlc__)
#  define CS_THR_MIN 1014
#else
#  define CS_THR_MIN 128
#endif

/* Cache line size, or multiple thereof */
/*--------------------------------------*/

#define CS_CL_SIZE 64

/*----------------------------------------------------------------------------
 * Type independent min an max (caution: the argument is evaluated)
 *----------------------------------------------------------------------------*/

#define CS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value of a */
#define CS_MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */
#define CS_MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*----------------------------------------------------------------------------
 * Variable interlace type:
 * {x1, y1, z1, x2, y2, z2, ...,xn, yn, zn} if interlaced
 * {x1, x2, ..., xn, y1, y2, ..., yn, z1, z2, ..., zn} if non interlaced
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_INTERLACE,          /* Variable is interlaced */
  CS_NO_INTERLACE        /* Variable is not interlaced */

} cs_interlace_t;

/*----------------------------------------------------------------------------
 * Macro used to silence "unused argument" warnings.
 *
 * This is useful when a function must match a given function pointer
 * type, but does not use all possible arguments.
 *----------------------------------------------------------------------------*/

#define CS_UNUSED(x) (void)(x)
#define CS_NO_WARN_IF_UNUSED(x) (void)(x)

/*----------------------------------------------------------------------------
 * Macros for compilation with a C++ compiler
 *----------------------------------------------------------------------------*/

#undef BEGIN_C_DECLS
#undef   END_C_DECLS

#if defined(__cplusplus)
#  define BEGIN_C_DECLS  extern "C" {
#  define   END_C_DECLS  }
#else
#  define BEGIN_C_DECLS
#  define   END_C_DECLS
#endif

/*----------------------------------------------------------------------------
 * Macros for Fortran interoperability
 *----------------------------------------------------------------------------*/

/*
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 */

#if !defined (__hpux)
#define CS_PROCF(x, y) x##_
#else
#define CS_PROCF(x, y) x
#endif

/*
 * Macro used to handle automatic "Fortran string length" arguments
 * (not used by Code_Saturne calls, but set by many compilers).
 * Some compilers, like the Fujitsu VPP 5000 compiler in its time, may not
 * support the variable length lists in mixed C/Fortran calls.
 */

#if defined (__uxpv__)  /* Fujitsu VPP 5000 case */
#define CS_ARGF_SUPP_CHAINE
#else
#define CS_ARGF_SUPP_CHAINE , ...
#endif

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Sizes and names associated with datatypes */

extern const size_t   cs_datatype_size[];
extern const char    *cs_datatype_name[];

/* MPI Datatypes associated with Code_Saturne datatypes */

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

extern MPI_Datatype   cs_datatype_to_mpi[];

#endif

/* Global variables indicationg task state */

extern int  cs_glob_n_threads;      /* Number of threads */

extern int  cs_glob_rank_id;        /* Rank in main MPI communicator */
extern int  cs_glob_n_ranks;        /* Size of main MPI communicator */

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

extern MPI_Comm       cs_glob_mpi_comm;      /* Main MPI intra-communicator */

#endif

/*=============================================================================
 * Public functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Given a base index i, return the next index aligned with a size m.
 *
 * \param[in]  i   base index
 * \param[in]  m   block size to align with
 *
 * \return aligned index
 */
/*----------------------------------------------------------------------------*/

inline static cs_lnum_t
cs_align(cs_lnum_t  i,
         cs_lnum_t  m)
{
  return ((i > 0) ? ((i-1)/m+1)*m : 0);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_DEFS_H__ */
