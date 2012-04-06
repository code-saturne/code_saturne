/*============================================================================
 * Base macro and typedef definitions for system portability
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/* Sizes associated with datatypes */

const size_t  cs_datatype_size[] = {0,
                                    1,
                                    sizeof(float),
                                    sizeof(double),
                                    4,
                                    8,
                                    4,
                                    8};

const char   *cs_datatype_name[] = {"",
                                    "char",
                                    "float",
                                    "double",
                                    "int32",
                                    "int64",
                                    "uint32",
                                    "uint64"};

#if defined(HAVE_MPI)

/* MPI Datatypes associated with Code_Saturne datatypes */

MPI_Datatype  cs_datatype_to_mpi[] = {MPI_DATATYPE_NULL,
                                      MPI_CHAR,
                                      MPI_FLOAT,
                                      MPI_DOUBLE,
                                      MPI_INT,            /* CS_INT32 */
                                      MPI_LONG_INT,       /* CS_INT64 */
                                      MPI_UNSIGNED,       /* CS_UINT32 */
                                      MPI_UNSIGNED_LONG}; /* CS_UINT64 */


#endif

/* Global variables indicationg task state */

int  cs_glob_n_threads = 1;    /* Number of threads */

int  cs_glob_rank_id = -1;     /* Rank of process in communicator */
int  cs_glob_n_ranks =  1;     /* Number of processes in communicator */

#if defined(HAVE_MPI)

MPI_Comm  cs_glob_mpi_comm = MPI_COMM_NULL;   /* Main MPI intra-communicator */

#endif

/*----------------------------------------------------------------------------*/

END_C_DECLS
