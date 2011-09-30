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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_parall.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

#if defined(HAVE_MPI)

/* Basic communicator info */

static MPI_Comm  _fvm_mpi_parall_comm = MPI_COMM_NULL;  /* Intra-communicator */
static int       _fvm_mpi_parall_size = 1;
static int       _fvm_mpi_parall_rank = 0;

/* MPI Datatypes associated with fvm datatypes */

MPI_Datatype  fvm_datatype_to_mpi[] = {MPI_DATATYPE_NULL,
                                       MPI_CHAR,
                                       MPI_FLOAT,
                                       MPI_DOUBLE,
                                       MPI_INT,            /* FVM_INT32 */
                                       MPI_LONG_INT,       /* FVM_INT64 */
                                       MPI_UNSIGNED,       /* FVM_UINT32 */
                                       MPI_UNSIGNED_LONG}; /* FVM_UINT64 */

/* Minimum recommended scatter/gather buffer size */

static size_t _fvm_parall_min_coll_buf_size = 1024*1024*8;

#endif

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Return default MPI communicator for FVM library functions.
 *
 * returns:
 *   handle to MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
fvm_parall_get_mpi_comm(void)
{
  return _fvm_mpi_parall_comm;
}

/*----------------------------------------------------------------------------
 * Set default MPI communicator for FVM library functions.
 *
 * parameters:
 *   comm <-- handle to MPI communicator
 *----------------------------------------------------------------------------*/

void
fvm_parall_set_mpi_comm(const MPI_Comm  comm)
{
  int mpi_flag;

  MPI_Initialized(&mpi_flag);

  /* Set communicator */

  _fvm_mpi_parall_comm = comm;

  if (mpi_flag != 0 && _fvm_mpi_parall_comm != MPI_COMM_NULL) {
    MPI_Comm_size(_fvm_mpi_parall_comm, &_fvm_mpi_parall_size);
    MPI_Comm_rank(_fvm_mpi_parall_comm, &_fvm_mpi_parall_rank);
  }
  else {
    _fvm_mpi_parall_size = 1;
    _fvm_mpi_parall_rank = 0;
  }

  /* Check (correct) fvm_datatype_to_mpi values */

  if (mpi_flag != 0)
  {
    int size_short, size_int, size_long, size_long_long;

    MPI_Type_size(MPI_SHORT, &size_short);
    MPI_Type_size(MPI_INT,   &size_int);
    MPI_Type_size(MPI_LONG,  &size_long);

#if defined(MPI_LONG_LONG)
    MPI_Type_size(MPI_LONG_LONG, &size_long_long);
#else
    size_long_long = 0;
#endif

    if (size_int == 4) {
      fvm_datatype_to_mpi[FVM_INT32] = MPI_INT;
      fvm_datatype_to_mpi[FVM_UINT32] = MPI_UNSIGNED;
    }
    else if (size_short == 4) {
      fvm_datatype_to_mpi[FVM_INT32] = MPI_SHORT;
      fvm_datatype_to_mpi[FVM_UINT32] = MPI_UNSIGNED_SHORT;
    }
    else if (size_long == 4) {
      fvm_datatype_to_mpi[FVM_INT32] = MPI_LONG;
      fvm_datatype_to_mpi[FVM_UINT32] = MPI_UNSIGNED_LONG;
    }

    if (size_int == 8) {
      fvm_datatype_to_mpi[FVM_INT64] = MPI_INT;
      fvm_datatype_to_mpi[FVM_UINT64] = MPI_UNSIGNED;
    }
    else if (size_long == 8) {
      fvm_datatype_to_mpi[FVM_INT64] = MPI_LONG;
      fvm_datatype_to_mpi[FVM_UINT64] = MPI_UNSIGNED_LONG;
    }
#if defined(MPI_LONG_LONG)
    else if (size_long_long == 8) {
      fvm_datatype_to_mpi[FVM_INT64] = MPI_LONG_LONG;
#if defined(MPI_UNSIGNED_LONG_LONG)
      fvm_datatype_to_mpi[FVM_UINT64] = MPI_UNSIGNED_LONG_LONG;
#else
      fvm_datatype_to_mpi[FVM_UINT64] = MPI_LONG_LONG;
#endif
    }
#endif
  }

}

#endif

/*----------------------------------------------------------------------------
 * Return rank of current process among associated program processes.
 *
 * returns:
 *   rank of current process in current communicator, or 0 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvm_parall_get_rank(void)
{
#if defined(HAVE_MPI)
  return _fvm_mpi_parall_rank;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------
 * Return number of processes associated with the current program.
 *
 * returns:
 *   number of processes in current communicator, or 1 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvm_parall_get_size(void)
{
#if defined(HAVE_MPI)
  return _fvm_mpi_parall_size;
#else
  return 1;
#endif
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Sum counters on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

void
fvm_parall_counter(cs_gnum_t   cpt[],
                   const int   n)
{

  if (_fvm_mpi_parall_size > 1) {

    int        i;
    cs_gnum_t *sum;
    cs_gnum_t _sum[64];

    if (n > 64)
      BFT_MALLOC(sum, n, cs_gnum_t);
    else
      sum = _sum;

    MPI_Allreduce(cpt, sum, n, CS_MPI_GNUM, MPI_SUM,
                  _fvm_mpi_parall_comm);

    for (i = 0; i < n ; i++)
      cpt[i] = sum[i];

    if (sum != _sum)
      BFT_FREE(sum);

  }

}

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

void
fvm_parall_counter_max(cs_lnum_t   cpt[],
                       const int   n)
{

  if (_fvm_mpi_parall_size > 1) {

    int        i;
    cs_lnum_t *maxval;
    cs_lnum_t _maxval[64];

    if (n > 64)
      BFT_MALLOC(maxval, n, cs_lnum_t);
    else
      maxval = _maxval;

    MPI_Allreduce(cpt, maxval, n, CS_MPI_LNUM, MPI_MAX,
                  _fvm_mpi_parall_comm);

    for (i = 0; i < n ; i++)
      cpt[i] = maxval[i];

    if (maxval != _maxval)
      BFT_FREE(maxval);

  }

}

#endif /* defined(HAVE_MPI) */

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
fvm_parall_get_min_coll_buf_size(void)
{
#if defined(HAVE_MPI)
  return _fvm_parall_min_coll_buf_size;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------
 * Define minimum recommended gather buffer size.
 *
 * This is used by FVM's internal strided and indexed array scatter/gather
 * algorithms, for non MPI-IO Input/output.
 *
 * parameters:
 *   minimum recommended gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

void
fvm_parall_set_min_coll_buf_size(size_t buffer_size)
{
#if defined(HAVE_MPI)
  _fvm_parall_min_coll_buf_size = buffer_size;
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
