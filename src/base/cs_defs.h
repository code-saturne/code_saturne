#ifndef __CS_DEFS_H__
#define __CS_DEFS_H__

/*============================================================================
 * Base macro and typedef definitions for system portability
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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
#  include "base/cs_config.h"
#endif

#ifdef __cplusplus
#include <type_traits>
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

#  define _LIBINTL_H /* Prevent inclusion of <libintl.h> by other files
                        with incorrect or missing checks;
                        TODO locate files causing issues to avoid
                        requiring this workaround */

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

#if _OPENMP >= 201307     /* simd construct available from OpenMP 4.0 */
#undef HAVE_OPENMP_SIMD
#define HAVE_OPENMP_SIMD 1
#endif

#endif

/* Do we have accelerator support ? */

#if defined(HAVE_CUDA)
#define HAVE_ACCEL 1
#elif defined(HAVE_SYCL)
#define HAVE_ACCEL 1
#elif defined(HAVE_OPENMP_TARGET)
#define HAVE_ACCEL 1
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

#  if defined(__GNUC__) || defined(__clang__) || defined (__NVCC__)
#    define restrict __restrict__
#  elif defined(_MSC_VER)
#    define restrict __restrict
#  else
#    ifndef HAVE_RESTRICT /* Must be provided by caller */
#      define restrict
#    endif
#  endif

#endif /* __cplusplus */

/* Definition of a DEPRECATED decorator which may be called in the code */

#if defined(__GNUC__) || defined(__clang__)
#  define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#  define DEPRECATED __declspec(deprecated)
#else
#  define DEPRECATED
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

/* C++ assert necessary for template */
#if defined(__cplusplus)
#include <typeinfo>
#include "assert.h"
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
 * General types and macros used throughout code_saturne
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
 * Basic types used by code_saturne
 * They may be modified here to better map to a given library, with the
 * following constraints:
 *  - cs_lnum_t must be signed
 *  - cs_gnum_t may be signed or unsigned
 *----------------------------------------------------------------------------*/

/* Global integer index or number */

#if defined(HAVE_LONG_GNUM)
  typedef uint64_t  cs_gnum_t;
#else
  typedef unsigned  cs_gnum_t;
#endif

/* Local integer index or number */

#if defined(HAVE_LONG_LNUM)
  typedef int64_t  cs_lnum_t;
#else
  typedef int  cs_lnum_t;
#endif

/* Other types */

typedef double             cs_coord_t;   /* Real number (coordinate value) */

typedef double              cs_real_t;   /* Fortran double precision */
typedef char                cs_byte_t;   /* Byte (untyped memory unit) */
typedef unsigned short int  cs_flag_t;   /* Flag storing metadata */

typedef double              cs_nreal_t;  /* Real number (normalized value) */
typedef double              cs_dreal_t;  /* Local distance */
typedef double              cs_rreal_t;  /* Reconstruction distance */

/* Vector or array block types */

typedef cs_lnum_t  cs_lnum_2_t[2];      /* Vector of 2 local numbers */
typedef cs_lnum_t  cs_lnum_3_t[3];      /* Vector of 3 local numbers */

typedef cs_coord_t cs_coord_3_t[3];         /* Vector of 3 real (coordinate)
                                               values */

typedef cs_real_t  cs_real_2_t[2];          /* Vector of 2 real values */
typedef cs_real_t  cs_real_3_t[3];          /* Vector of 3 real values */
typedef cs_real_t  cs_real_4_t[4];          /* Vector of 4 real values */
typedef cs_real_t  cs_real_6_t[6];          /* Vector of 6 real values
                                               (for symmetric tensor) */
typedef cs_real_t  cs_real_9_t[9];          /* Vector of 9 real values */
typedef cs_real_t  cs_real_10_t[10];        /* Vector of 10 real values */

typedef cs_real_t  cs_real_23_t[2][3];      /* Matrix of 2x3 real values */

typedef cs_real_t  cs_real_33_t[3][3];      /* Matrix of 3x3 real values */
typedef cs_real_t  cs_real_66_t[6][6];      /* Matrix of 6x6 real values */
typedef cs_real_t  cs_real_99_t[9][9];      /* Matrix of 9x9 real values */

typedef cs_real_t  cs_real_333_t[3][3][3];  /* tensor of 3x3x3 real values */

typedef cs_real_t  cs_real_34_t[3][4];      /* Matrix of 3x4 real values */

typedef cs_real_t  cs_real_63_t[6][3];      /* Matrix of 6x3 real values */

typedef cs_real_t  cs_real_69_t[6][9];      /* Matrix of 6x9 real values */

typedef cs_real_33_t  cs_real_332_t[2];     /* vector of 2 3x3 matrices
                                               of real values */
typedef cs_real_66_t  cs_real_662_t[2];     /* vector of 2 6x6 matrices
                                               of real values */

typedef cs_nreal_t  cs_nreal_3_t[3];        /* Vector of normalized real values
                                               (i.e. unit vector) */
typedef cs_dreal_t  cs_dreal_3_t[3];        /* Vector joining 2 local coordinates */
typedef cs_rreal_t  cs_rreal_3_t[3];        /* Reconstruction distance Vector */

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

/* MPI type for cs_lnum_t type */

#  if defined(HAVE_LONG_LNUM)
#    define CS_MPI_LNUM     MPI_LONG
#  else
#    define CS_MPI_LNUM     MPI_INT
#  endif

#  define CS_MPI_EFLAG    MPI_UNSIGNED       /* MPI type for cs_mflag_t type */
#  define CS_MPI_FLAG     MPI_UNSIGNED_SHORT /* MPI type for cs_flag_t type */
#  define CS_MPI_COORD    MPI_DOUBLE         /* MPI type for cs_coord_t type */

#endif /* defined(HAVE_MPI) && !defined(CS_IGNORE_MPI) */

/* Mappings to code_saturne datatypes */
/*------------------------------------*/

#if defined(HAVE_LONG_GNUM)
# define CS_GNUM_TYPE     CS_UINT64
#elif (SIZEOF_INT == 8)
# define CS_GNUM_TYPE     CS_UINT64
#else
# define CS_GNUM_TYPE     CS_UINT32
#endif

#if defined(HAVE_LONG_LNUM)
# if (SIZEOF_LONG == 8)
#  define CS_LNUM_TYPE     CS_INT64
# else
#  define CS_LNUM_TYPE     CS_INT32
# endif
#else
# if (SIZEOF_INT == 8)
#  define CS_LNUM_TYPE     CS_INT64
# else
#  define CS_LNUM_TYPE     CS_INT32
# endif
#endif

#if (SIZEOF_INT == 8)
# define CS_INT_TYPE      CS_INT64
#else
# define CS_INT_TYPE      CS_INT32
#endif

#if (SIZEOF_INT == 8)
# define CS_UINT_TYPE     CS_UINT64
#else
# define CS_UINT_TYPE     CS_UINT32
#endif

#define CS_FLAG_TYPE      CS_UINT16
#define CS_EFLAG_TYPE     CS_UINT_TYPE
#define CS_REAL_TYPE      CS_DOUBLE
#define CS_COORD_TYPE     CS_DOUBLE

/* Minimum size for OpenMP loops
 *  (will need benchmarking and tuning for various systems)
 *---------------------------------------------------------*/

#define CS_THR_MIN 128

/* Cache line size, or multiple thereof */
/*--------------------------------------*/

#define CS_CL_SIZE 64

/*----------------------------------------------------------------------------
 * Type independent min an max (caution: the argument is evaluated)
 *----------------------------------------------------------------------------*/

#define CS_ABS(a)     ((a) <  0  ? -(a) : (a))  /*!< Absolute value of a */
#define CS_MIN(a,b)   ((a) < (b) ?  (a) : (b))  /*!< Minimum of a et b */
#define CS_MAX(a,b)   ((a) > (b) ?  (a) : (b))  /*!< Maximum of a et b */

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
 * Macros for function type qualifiers
 *----------------------------------------------------------------------------*/

#ifdef __CUDACC__

#define CS_F_HOST __host__
#define CS_F_DEVICE __device__
#define CS_F_HOST_DEVICE __host__ __device__
#define CS_V_CONSTANT __constant__

#else

#define CS_F_HOST
#define CS_F_DEVICE
#define CS_F_HOST_DEVICE
#define CS_V_CONSTANT static const

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

/*=============================================================================
 * Global variables
 *============================================================================*/

/* Empty but non-null string */

extern const char     cs_empty_string[];

/* Sizes and names associated with datatypes */

extern const size_t   cs_datatype_size[];
extern const char    *cs_datatype_name[];

/* MPI Datatypes associated with code_saturne datatypes */

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

extern MPI_Datatype   cs_datatype_to_mpi[];

#endif

/* Global variables indicating task state */

extern int  cs_glob_n_threads;     /* Number of threads */

extern int  cs_glob_rank_id;       /* Rank in main MPI communicator */
extern int  cs_glob_n_ranks;       /* Size of main MPI communicator */

extern int  cs_glob_node_rank_id;  /* Rank on node in main MPI communicator */
extern int  cs_glob_node_n_ranks;  /* Number of ranks on node of main
                                      MPI communicator */

#if defined(HAVE_MPI) && !defined(CS_IGNORE_MPI)

extern MPI_Comm       cs_glob_mpi_comm;      /* Main MPI intra-communicator */

#endif

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a given structure
 *
 * \param[in, out]  s  pointer to a structure cast on-the-fly
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_destructor_t)(void  *s);

/*=============================================================================
 * Public functions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
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
/*!
 * \brief Retrieve the associated thread id (0 if no OpenMP or if outside an
 *        OpenMP region)
 *
 * \return the id of the OpenMP thread
 */
/*----------------------------------------------------------------------------*/

inline static int
cs_get_thread_id(void)
{
#if defined(HAVE_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus

/*=============================================================================
 * Public C++ templates
 *============================================================================*/

namespace cs {

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return absolute value of a number.
 *
 * This function is overloaded in cs_math.h for floating-point values.
 *
 * \tparam T value type
 *
 * \param[in]  a  value
 *
 * \return absolute value of a
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline T
abs(const T  a)
{
  return ((a) <  0  ? -(a) : (a));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return minimum of two values.
 *
 * This function is overloaded in cs_math.h for floating-point values.
 *
 * \tparam T value type
 *
 * \param[in]  a   first value
 * \param[in]  b   second value
 *
 * \return minimum of a and b
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline T
min(const T  a,
    const T  b)
{
  return ((a) < (b) ?  (a) : (b));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief  Return maximum of two values.
 *
 * This function is overloaded in cs_math.h for floating-point values.
 *
 * \tparam T value type
 *
 * \param[in]  a   first value
 * \param[in]  b   second value
 *
 * \return maximum of a and b
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline T
max(const T  a,
    const T  b)
{
  return ((a) > (b) ?  (a) : (b));
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Clamp function for a given scalar value.
 *
 * This function is overloaded in cs_math.h for floating-point values.
 *
 * \tparam T value type
 *
 * \param[in] x    initial value
 * \param[in] xmin min value for clamping
 * \param[in] xmax max value for clamping
 *
 * \return clamped value which is 'x' if xmin < x < xmax or lower/upper bound
 * otherwise
 */
/*----------------------------------------------------------------------------*/

template <typename T>
CS_F_HOST_DEVICE inline T
clamp(const T x,
      const T xmin,
      const T xmax)
{
  T x_tmp = ((xmin) > (x) ?  (xmin) : (x));
  return ((xmax) < (x_tmp) ?  (xmax) : (x_tmp));
}

} // namespace cs

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the cs_datatype_t from a typename
 *
 * \tparam T : datatype
 */
/*----------------------------------------------------------------------------*/

template <typename T>
struct dependent_false : std::false_type {};

template <typename T>
static inline cs_datatype_t
cs_datatype_from_type()
{
  static_assert(dependent_false<T>::value, "Unknown datatype");
  return CS_DATATYPE_NULL;
}

/*----------------------------------------------------------------------------*/
/* Specialized versions of the templated function                             */
/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<char>()
{
  return CS_CHAR;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<float>()
{
  return CS_FLOAT;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<double>()
{
  return CS_DOUBLE;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<uint16_t>()
{
  return CS_UINT16;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<uint32_t>()
{
  return CS_UINT32;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<uint64_t>()
{
  return CS_UINT64;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<int32_t>()
{
  return CS_INT32;
}

/*----------------------------------------------------------------------------*/

template <>
constexpr inline cs_datatype_t
cs_datatype_from_type<int64_t>()
{
  return CS_INT64;
}

/*----------------------------------------------------------------------------*/

#endif /* __cplusplus */

/*----------------------------------------------------------------------------*/

#endif /* __CS_DEFS_H__ */
