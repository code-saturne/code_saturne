#ifndef _SYR_DEF_H_
#define _SYR_DEF_H_

/*============================================================================
 * Definitions, Global variables, and basic functions
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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
 * Type definitions
 *============================================================================*/

/* Enum for datatype description */

typedef enum {
  SYR_TYPE_char,
  SYR_TYPE_int,
  SYR_TYPE_float,
  SYR_TYPE_double,
  SYR_TYPE_void
} syr_type_t;

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* "Classical" macros */

#define SYR_MIN(a,b)   ((a) > (b) ?  (b) : (a))
#define SYR_MAX(a,b)   ((a) < (b) ?  (b) : (a))

/* Copy of SYRTHES macros for C/Fortran interoperability */

#if    defined(__sgi)     || defined(__uxpv__) || defined(__aix__) \
    || defined(__linux__)

#define name2(a,b)  a##b
#if !defined(__xlc__)
#define proc(x,y) name2(x,_)
#else
#define proc(x,y) x
#endif
#define proci(x)  x

#elif defined(sun) || defined(__alpha)

#define name2(a,b)  a/**/b
#define proc(x,y) name2(x,_)
#define proci(x)  x

#elif defined(CRAY)

#define name2(a,b)  a/**/b
#define proc(x,y) y
#define proci(x)  x

#else

#define name2(a,b)  a/**/b
#define proc(x,y)  x
#define proci(x) name2(x,_)

#endif

/*============================================================================
 * Global variables
 *============================================================================*/

extern char syr_glob_build_date[];  /* Build date */

extern MPI_Comm  syr_glob_mpi_comm;

/*============================================================================
 * Function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Error initialization and handling
 *----------------------------------------------------------------------------*/

void
syr_errhandler_initialize(void);

/*----------------------------------------------------------------------------
 * Initialize MPI communication
 *----------------------------------------------------------------------------*/

void
syr_mpi_initialize(int    *argc,
                   char  **argv[]);

/*----------------------------------------------------------------------------
 * Finalize MPI communication
 *----------------------------------------------------------------------------*/

void
syr_mpi_finalize(void);

/*----------------------------------------------------------------------------
 * Force abort of MPI communication (for atexit)
 *
 * This function only forces finalization if an MPI coupling world is defined,
 * so it should do nothing in case of a normal exit, in which the MPI coupling
 * world info structure should have been destroyed through a regular call
 * to syr_mpi_finalize.
 *----------------------------------------------------------------------------*/

void
syr_mpi_exit_force(void);

/*----------------------------------------------------------------------------
 * Recover rank information on a given application number
 *
 * parameters:
 *   app_name  <-- application name
 *   root_rank --> associated root rank
 *   n_ranks   --> number of associated ranks
 *----------------------------------------------------------------------------*/

void
syr_mpi_appinfo(const char  *app_name,
                int         *root_rank,
                int         *n_ranks);

/*----------------------------------------------------------------------------
 * Exit / Stop
 *----------------------------------------------------------------------------*/

void
syr_exit(int status);

/*----------------------------------------------------------------------------
 * Print warning
 *----------------------------------------------------------------------------*/

void
syr_warn(void);

/*----------------------------------------------------------------------------*/

#endif /* _SYR_DEF_H_ */
