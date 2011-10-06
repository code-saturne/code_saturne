#ifndef __FVM_PARALL_H__
#define __FVM_PARALL_H__

/*============================================================================
 * Base functions for parallelism
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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

#if defined(HAVE_MPI)

#define FVM_MPI_TAG      (int)('F'+'V'+'M') /* MPI tag for FVM operations */

/* MPI type for fvm_gnum_t integer type (depends on configuration) */

#if defined(HAVE_LONG_GNUM)
  #if (SIZEOF_LONG == 8)
    #define FVM_MPI_GNUM     MPI_UNSIGNED_LONG
  #elif (SIZEOF_LONG_LONG == 8)
    #if defined(MPI_UNSIGNED_LONG_LONG)
      #define FVM_MPI_GNUM     MPI_UNSIGNED_LONG_LONG
    #elif defined(MPI_LONG_LONG)
      #define FVM_MPI_GNUM     MPI_LONG_LONG
    #endif
  #endif
  #if !defined(FVM_MPI_GNUM)
    #error
  #endif
#else
  #define FVM_MPI_GNUM       MPI_UNSIGNED
#endif

#define FVM_MPI_LNUM     MPI_INT         /* MPI type for fvm_lnum_t type */
#define FVM_MPI_COORD    MPI_DOUBLE      /* MPI type for fvm_coord_t type */

#endif

/*=============================================================================
 * Static global variables
 *============================================================================*/

/* MPI Datatypes associated with fvm datatypes */

#if defined(HAVE_MPI)

extern MPI_Datatype  fvm_datatype_to_mpi[];

#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Return default MPI communicator for FVM library functions.
 *
 * returns:
 *   handle to MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
fvm_parall_get_mpi_comm(void);

/*----------------------------------------------------------------------------
 * Set default MPI communicator for FVM library functions.
 *
 * parameters:
 *   comm <-- handle to MPI communicator
 *----------------------------------------------------------------------------*/

void
fvm_parall_set_mpi_comm(const MPI_Comm  comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return rank of current process among associated program processes.
 *
 * returns:
 *   rank of current process in current communicator, or 0 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvm_parall_get_rank(void);

/*----------------------------------------------------------------------------
 * Return number of processes associated with the current program.
 *
 * returns:
 *   number of processes in current communicator, or 1 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvm_parall_get_size(void);

/*----------------------------------------------------------------------------
 * Sum counters on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void
fvm_parall_counter(fvm_gnum_t  cpt[],
                   const int   n);
#else
#define fvm_parall_counter(_cpt, _n)
#endif

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void
fvm_parall_counter_max(fvm_lnum_t  cpt[],
                       const int   n);
#else
#define fvm_parall_counter_max(_cpt, _n)
#endif

/*----------------------------------------------------------------------------
 * Return minimum recommended scatter or gather buffer size.
 *
 * This is used by FVM's internal strided and indexed array scatter/gather
 * algorithms, for non MPI-IO Input/output.
 *
 * returns:
 *   minimum recommended gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

size_t
fvm_parall_get_min_coll_buf_size(void);

/*----------------------------------------------------------------------------
 * Define minimum recommended scatter or gather buffer size.
 *
 * This is used by FVM's internal strided and indexed array scatter/gather
 * algorithms, for non MPI-IO Input/output.
 *
 * parameters:
 *   minimum recommended gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

void
fvm_parall_set_min_coll_buf_size(size_t buffer_size);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVM_PARALL_H__ */
