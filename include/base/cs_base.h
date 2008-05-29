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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-----------------------------------------------------------------------------*/

#if defined(_CS_HAVE_MPI)

#include <mpi.h>

#if defined(_CS_HAVE_MPE)
#include <mpe.h>
#endif

#endif

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Nom du système */

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

#define CS_DIM_3              3                 /* Spatial dimension */

/* "Classical" macros */

#define CS_ABS(a)     ((a) <  0  ? -(a) : (a))  /* Absolute value of a */
#define CS_MIN(a,b)   ((a) > (b) ?  (b) : (a))  /* Minimum of a et b */
#define CS_MAX(a,b)   ((a) < (b) ?  (b) : (a))  /* Maximum of a et b */

/*
 * Macros for future internationalization via gettext() or a similar
 * function (to mark translatable character strings)
 */

#undef _
#define _(String) String

#undef N_
#define N_(String) String

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

#if defined(_CS_HAVE_MPI)

#define CS_MPI_INT       MPI_INT         /* If cs_int_t is an int;
                                            otherwise redefine MPI_xxx */
#define CS_MPI_REAL      MPI_DOUBLE      /* If cs_real_t is a double;
                                            otherwise redefine as MPI_REAL */
#define CS_MPI_REAL_INT  MPI_DOUBLE_INT  /* If cs_real_t est un double ;
                                            otherwise redefine as MPI_REAL_INT */

typedef struct
{
  cs_real_t val;
  cs_int_t  rang;
} cs_mpi_real_int_t;

#endif /* defined(_CS_HAVE_MPI) */

/* Datatype enumeration to transmit a data's type to a function */

typedef enum {
  CS_TYPE_char,
  CS_TYPE_cs_int_t,
  CS_TYPE_cs_real_t,
  CS_TYPE_cs_bool_t,
  CS_TYPE_cs_point_t,
  CS_TYPE_void
} cs_type_t;

/* Exit function prototype */

typedef void (cs_exit_t) (int status);

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

extern cs_int_t  cs_glob_base_rang;        /* Rank of process in group */
extern cs_int_t  cs_glob_base_nbr;         /* Number of processes in group */

#if defined(_CS_HAVE_MPI)
extern MPI_Comm  cs_glob_base_mpi_comm;    /* Intra-communicator */
#endif

/* Global variables used for MPE instrumentation */

#if defined(_CS_HAVE_MPI) && defined(_CS_HAVE_MPE)
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
 * Fonction d'arret depuis du code Fortran
 *
 * Interface Fortran :
 *
 * SUBROUTINE CSEXIT (STATUT)
 * *****************
 *
 * INTEGER          STATUT      : --> : 0 pour succès, 1 ou + pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *const statut
);


/*----------------------------------------------------------------------------
 * Temps CPU écoulé depuis le début de l'exécution
 *
 * Interface Fortran :
 *
 * SUBROUTINE DMTMPS (TCPU)
 * *****************
 *
 * DOUBLE PRECISION TCPU        : --> : temps CPU (utilisateur + système)
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *const tcpu
);


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI.
 *
 * Global variables `cs_glob_base_nbr' (number of Code_Saturne processes)
 * and `cs_glob_base_rang' (rank of local process) are set by this function.
 *
 * parameters:
 *   argc     <-- pointer to number of command line arguments.
 *   argv     <-- pointer to command line arguments array.
 *   rang_deb <-- rank of the first process of this group in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

void
cs_base_mpi_init(int     *argc,
                 char  ***argv,
                 int      rang_deb);

/*----------------------------------------------------------------------------
 * Finalize MPI.
 *----------------------------------------------------------------------------*/

void
cs_base_mpi_fin(void);

#endif /* defined(_CS_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Exit, with handling for both normal and error cases.
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
cs_base_erreur_init(void);

/*----------------------------------------------------------------------------
 * Initialize management of memory allocated through BFT.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init(void);



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
cs_base_info_systeme(void);

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
 *----------------------------------------------------------------------------*/

char *
cs_base_chaine_f_vers_c_cree(const char  *f_str,
                             int          l_en);

/*----------------------------------------------------------------------------
 * Free a string converted from the Fortran API to the C API.
 *
 * parameters:
 *   str <-- C string
 *----------------------------------------------------------------------------*/

char  *
cs_base_chaine_f_vers_c_detruit(char  *c_str);

/*----------------------------------------------------------------------------
 * Modify program exit behavior setting a specific exit function.
 *
 * parameters:
 *   exit_func <-- pointer to exit function which should be used
 *----------------------------------------------------------------------------*/

void
cs_base_exit_set(cs_exit_t *exit_func);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_BASE_H__ */
