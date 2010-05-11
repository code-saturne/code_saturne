#ifndef __PLE_COUPLING_H__
#define __PLE_COUPLING_H__

/*============================================================================
 * Set up communication with coupled codes.
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2010  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

#include "ple_config.h"

#if defined(PLE_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/* Opaque code coupling information structure */

typedef struct _ple_coupling_mpi_world_t  ple_coupling_mpi_world_t;

/* Info structure for code coupling */

typedef struct {

  int          app_num;   /* Application number, or -1 if undefined */
  int          root_rank; /* Application root rank in MPI_COMM_WORLD */
  int          n_ranks;   /* Number of ranks associated with application */
  const char  *app_type;  /* Application type name (may be empty) */
  const char  *app_name;  /* Application instance name (may be empty) */

} ple_coupling_mpi_world_info_t;

#endif /* defined(PLE_HAVE_MPI) */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI_COMM_WORLD.
 *
 * The application communicator app_comm is usually obtained from
 * MPI_COMM_WORLD using MPI_Comm_split, with app_num corresponding to
 * the "color" argument in that function.
 *
 * As this function requires communication between applications, it
 * is a collective function in MPI_COMM_WORLD.
 *
 * parameters:
 *   app_num   <-- application number in MPI_COMM_WORLD (nonnegative).
 *   app_name  <-- name of current application.
 *   case_name <-- name of current case, or NULL.
 *   app_comm  <-- communicator associated with local application.
 *
 * returns:
 *   PLE coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_world_t *
ple_coupling_mpi_world_create(int          app_num,
                              const char  *app_type,
                              const char  *app_name,
                              MPI_Comm     app_comm);

/*----------------------------------------------------------------------------
 * Free an PLE coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-> pointer to structure that should be freed.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_world_destroy(ple_coupling_mpi_world_t **w);

/*----------------------------------------------------------------------------
 * Return the number of applications in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   number of application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_world_n_apps(const ple_coupling_mpi_world_t  *w);

/*----------------------------------------------------------------------------
 * Return the id of the local application in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   id of the local application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_world_get_app_id(const ple_coupling_mpi_world_t  *w);

/*----------------------------------------------------------------------------
 * Return application information in MPI_COMM_WORLD.
 *
 * parameters:
 *   w      <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *   app_id <-- application id
 *
 * returns:
 *   application information structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_world_info_t
ple_coupling_mpi_world_get_info(const ple_coupling_mpi_world_t  *w,
                                int                              app_id);

/*----------------------------------------------------------------------------
 * Create an intracommunicator from a local and distant communicator
 * within MPI_COMM_WORLD.
 *
 * parameters:
 *   app_comm      <-- communicator associated with local application
 *   distant_root  <-- rank of distant group leader in MPI_COMM_WORLD
 *   new_comm      --> pointer to new communicator
 *   local_range   --> first and past-the last ranks of local application
 *                     in new communicator
 *   distant_range --> first and past-the last ranks of distant application
 *                     in new communicator
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_intracomm_create(MPI_Comm   app_comm,
                                  int        distant_root,
                                  MPI_Comm  *new_comm,
                                  int        local_range[2],
                                  int        distant_range[2]);

/*----------------------------------------------------------------------------
 * Dump printout of an PLE coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_world_dump(const ple_coupling_mpi_world_t  *w);

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_COUPLING_H__ */
