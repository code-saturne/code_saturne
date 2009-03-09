/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_BASE_H__
#define __CS_BASE_H__

/*============================================================================
 * Definitions, global variables, and base functions
 *============================================================================*/


/* Definition of the C langage version used (C89 or C99) */

#if defined(__STDC_VERSION__)
#  define _CS_STDC_VERSION __STDC_VERSION__
#else
#  define _CS_STDC_VERSION 1989
#endif

/*
 * Redefinition of "inline" et "restrict" qualifiers incompatible with
 * some C89 compilers (standard in C99)
 */

#if (_CS_STDC_VERSION < 199901L)

#  if defined(__GNUC__)
#    define inline __inline__
#    define restrict __restrict__
#  else
#    define inline
#    define restrict
#  endif

#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

#if (_CS_STDC_VERSION >= 199901L)
#include <stdint.h>
#endif

#if defined(HAVE_MPI)

#include <mpi.h>

#if defined(HAVE_MPE)
#include <mpe.h>
#endif

#endif

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Application type name */

#define CS_APP_NAME     "Code_Saturne"
#define CS_APP_VERSION  PACKAGE_VERSION  /* PACKAGE_VERSION from autoconf */

/* System type name */

#if defined(__sgi__) || defined(__sgi) || defined(sgi)
#define _CS_ARCH_IRIX_64

#elif defined(__hpux__) || defined(__hpux) || defined(hpux)
#define _CS_ARCH_HP_UX

#elif defined(__blrts__) || defined(__bgp__)
#define _CS_ARCH_Blue_Gene

#elif defined(__linux__) || defined(__linux) || defined(linux)
#define _CS_ARCH_Linux

#elif defined(__sun__) || defined(__sun) || defined(sun)
#define _CS_ARCH_SunOS

#elif defined(__uxpv__) || defined(__uxpv) || defined(uxpv)
#define _CS_ARCH_UNIX_System_V

#endif

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
 * Some compilers, like the Fujitsu VPP 5000 compiler, may not
 * support the variable length lists in mixed C/Fortran calls.
 */

#if defined (__uxpv__)  /* Fujitsu VPP 5000 case */
#define CS_ARGF_SUPP_CHAINE
#else
#define CS_ARGF_SUPP_CHAINE , ...
#endif

/* On certain architectures such as IBM Blue Gene, some operations may
 * be better optimized on memory-aligned data (if 0 here, no alignment
 * is leveraged). This alignment is not exploited yet in Code_Saturne. */

#if defined(__blrts__) || defined(__bgp__)
#define CS_MEM_ALIGN 16
#else
#define CS_MEM_ALIGN 0
#endif

/* "Classical" macros */

#define CS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value of a */
#define CS_MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */
#define CS_MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*
 * Macros for internationalization via gettext() or a similar
 * function (to mark translatable character strings)
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

#undef BEGIN_C_DECLS
#undef   END_C_DECLS

#if defined(__cplusplus)
#define BEGIN_C_DECLS  extern "C" {
#define   END_C_DECLS  }
#else
#define BEGIN_C_DECLS
#define   END_C_DECLS
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef int              cs_int_t;      /* Integer */
typedef double           cs_real_t;     /* Floating-point real */
typedef char             cs_byte_t;     /* Byte (untyped memory unit) */

typedef cs_real_t        cs_point_t[3];

typedef enum {                          /* Boolean */
  CS_FALSE ,
  CS_TRUE
} cs_bool_t;

#if !defined(false)
#define false CS_FALSE
#endif

#if !defined(true)
#define true CS_TRUE
#endif

/* Mappings to MPI datatypes */

#if defined(HAVE_MPI)

#define CS_MPI_INT       MPI_INT         /* If cs_int_t is an int */
#define CS_MPI_REAL      MPI_DOUBLE      /* If cs_real_t is a double */

#endif /* defined(HAVE_MPI) */

/* Datatype enumeration to transmit a data's type to a function */

typedef enum {
  CS_TYPE_char,
  CS_TYPE_cs_int_t,
  CS_TYPE_cs_real_t,
  CS_TYPE_cs_bool_t,
  CS_TYPE_cs_point_t,
  CS_TYPE_void
} cs_type_t;

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

extern int  cs_glob_n_threads;      /* Number of threads */

extern int  cs_glob_rank_id;        /* Rank of process in group */
extern int  cs_glob_n_ranks;        /* Number of processes in group */

#if defined(HAVE_MPI)
extern MPI_Comm  cs_glob_mpi_comm;    /* Intra-communicator */
#endif

/* Global variables used for MPE instrumentation */

#if defined(HAVE_MPI) && defined(HAVE_MPE)
extern int  cs_glob_mpe_broadcast_a;
extern int  cs_glob_mpe_broadcast_b;
extern int  cs_glob_mpe_synchro_a;
extern int  cs_glob_mpe_synchro_b;
extern int  cs_glob_mpe_send_a;
extern int  cs_glob_mpe_send_b;
extern int  cs_glob_mpe_rcv_a;
extern int  cs_glob_mpe_rcv_b;
extern int  cs_glob_mpe_reduce_a;
extern int  cs_glob_mpe_reduce_b;
extern int  cs_glob_mpe_compute_a;
extern int  cs_glob_mpe_compute_b;
#endif

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Call exit routine from Fortran code
 *
 * Fortran interface:
 *
 * SUBROUTINE CSEXIT (STATUS)
 * *****************
 *
 * INTEGER          STATUS      : --> : 0 for success, 1+ for error
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *status
);

/*----------------------------------------------------------------------------
 * CPU time used since execution start
 *
 * Fortran interface:
 *
 * SUBROUTINE DMTMPS (TCPU)
 * *****************
 *
 * DOUBLE PRECISION TCPU        : --> : CPU time (user + system)
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *tcpu
);

/*----------------------------------------------------------------------------
 * Check that main integer working array memory reservation fits within
 * the allocated size of IA.
 *
 * Fortran interface:
 *
 * SUBROUTINE IASIZE (CALLER, MEMINT)
 * *****************
 *
 * CHARACTER*6      CALLER      : --> : Name of calling subroutine
 * INTEGER          MEMINT      : --> : Last required element in IA
 *----------------------------------------------------------------------------*/

void CS_PROCF (iasize, IASIZE)
(
 const char   caller[6],
 cs_int_t    *memint
);

/*----------------------------------------------------------------------------
 * Check that main floating-point working array memory reservation fits
 * within the allocated size of RA.
 *
 * Fortran interface:
 *
 * SUBROUTINE RASIZE (CALLER, MEMINT)
 * *****************
 *
 * CHARACTER*6      CALLER      : --> : Name of calling subroutine
 * INTEGER          MEMRDP      : --> : Last required element in RA
 *----------------------------------------------------------------------------*/

void CS_PROCF (rasize, RASIZE)
(
 const char   caller[6],
 cs_int_t    *memrdp
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Complete MPI initialization.
 *
 * Global variables `cs_glob_n_ranks' (number of Code_Saturne processes)
 * and `cs_glob_rank_id' (rank of local process) are set by this function.
 *
 * parameters:
 *   argc    <-- pointer to number of command line arguments.
 *   argv    <-- pointer to command line arguments array.
 *   app_num <-- -1 if MPI is not needed, or application number in
 *               MPI_COMM_WORLD of this instance of Code_Saturne.
 *----------------------------------------------------------------------------*/

void
cs_base_mpi_init(int     *argc,
                 char  ***argv,
                 int      app_num);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Exit, with handling for both normal and error cases.
 *
 * Finalize MPI if necessary.
 *
 * parameters:
 *   status <-- value to be returned to the parent:
 *              EXIT_SUCCESS / 0 for the normal case,
 *              EXIT_FAILURE or other nonzero code for error cases.
 *----------------------------------------------------------------------------*/

void
cs_exit(int  status);

/*----------------------------------------------------------------------------
 * Initialize error and signal handlers.
 *----------------------------------------------------------------------------*/

void
cs_base_error_init(void);

/*----------------------------------------------------------------------------
 * Initialize management of memory allocated through BFT.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init(void);

/*----------------------------------------------------------------------------
 * Allocate Fortran work arrays and prepare for their use.
 *
 * parameters:
 *   iasize <-- integer working array size (maximum number of values)
 *   rasize <-- floating-point working array size (maximum number of values)
 *   ia     --> pointer to integer working array
 *   ra     --> pointer to floating-point working array
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init_work(size_t       iasize,
                      size_t       rasize,
                      cs_int_t   **ia,
                      cs_real_t  **ra);

/*----------------------------------------------------------------------------
 * Finalize management of memory allocated through BFT.
 *
 * A summary of the consumed memory is given.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_fin(void);

/*----------------------------------------------------------------------------
 * Print summary of running time, including CPU and elapsed times.
 *----------------------------------------------------------------------------*/

void
cs_base_bilan_temps(void);

/*----------------------------------------------------------------------------
 * Print available system information.
 *----------------------------------------------------------------------------*/

void
cs_base_system_info(void);

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This is necessary for good consistency of messages output from C or
 * from Fortran, and to handle parallel and serial logging options.
 *----------------------------------------------------------------------------*/

void
cs_base_bft_printf_set(void);

/*----------------------------------------------------------------------------
 * Print a warning message header.
 *
 * parameters:
 *   file_name <-- name of source file
 *   line_nume <-- line number in source file
 *----------------------------------------------------------------------------*/

void
cs_base_warn(const char  *file_name,
             int          line_num);

/*----------------------------------------------------------------------------
 * Convert a character string from the Fortran API to the C API.
 *
 * Eventual leading and trailing blanks are removed.
 *
 * parameters:
 *   f_str <-- Fortran string
 *   f_len <-- Fortran string length
 *
 * returns:
 *   pointer to C string
 *----------------------------------------------------------------------------*/

char *
cs_base_string_f_to_c_create(const char  *f_str,
                             int          f_len);

/*----------------------------------------------------------------------------
 * Free a string converted from the Fortran API to the C API.
 *
 * parameters:
 *   str <-> pointer to C string
 *----------------------------------------------------------------------------*/

void
cs_base_string_f_to_c_free(char  **c_str);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BASE_H__ */
