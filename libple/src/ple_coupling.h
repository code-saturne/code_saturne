#ifndef __PLE_COUPLING_H__
#define __PLE_COUPLING_H__

/*============================================================================
 * Set up communication with coupled codes.
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2020  EDF S.A.

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

/*
 * Mask for synchronization flag
 */

/* Command bits */

#define PLE_COUPLING_INIT             (1 << 0)  /* Not yet synchronized */

#define PLE_COUPLING_NO_SYNC          (1 << 1)  /* Not synchronized */
#define PLE_COUPLING_STOP             (1 << 2)  /* Will stop immediately */
#define PLE_COUPLING_LAST             (1 << 3)  /* Last synchronization */

/* Time stepping bits */

#define PLE_COUPLING_NEW_ITERATION    (1 << 4)
#define PLE_COUPLING_REDO_ITERATION   (1 << 5)

/* Time step value handling bits */

#define PLE_COUPLING_TS_MIN           (1 << 6)  /* Use smallest time step */
#define PLE_COUPLING_TS_LEADER        (1 << 7)  /* Prescribe time step for all
                                                   members of group (only one
                                                   member may set this flag) */

/* Calculation type or state information bits */

#define PLE_COUPLING_UNSTEADY         (1 << 8)
#define PLE_COUPLING_STEADY           (1 << 9)
#define PLE_COUPLING_CONVERGED        (1 << 10)

/* Optional user code information bits */

#define PLE_COUPLING_USER_1           (1 << 11)
#define PLE_COUPLING_USER_2           (1 << 12)
#define PLE_COUPLING_USER_3           (1 << 13)
#define PLE_COUPLING_USER_4           (1 << 14)

/*============================================================================
 * Type definitions
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/* Opaque code coupling information structure */

typedef struct _ple_coupling_mpi_set_t  ple_coupling_mpi_set_t;

/* Info structure for code coupling */

typedef struct {

  int          status;    /* Status flag for synchronization info */
  int          root_rank; /* Application root rank in MPI_COMM_WORLD */
  int          n_ranks;   /* Number of ranks associated with application */
  const char  *app_type;  /* Application type name (may be empty) */
  const char  *app_name;  /* Application instance name (may be empty) */

} ple_coupling_mpi_set_info_t;

#endif /* defined(PLE_HAVE_MPI) */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(PLE_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build a group id within a communicator based on its name.
 *
 * If multiple groups are present, ids are number from 0 to n_groups - 1,
 * based on the odering of group names. If all processes have the same
 * group name, the returned value is -1.
 *
 * The returned id may typically be used as a "color" argument for
 * MPI_Comm_split().
 *
 * As this function requires communication between applications, it
 * is a collective function in comm.
 *
 * parameters:
 *   comm       <-- MPI communicator.
 *   group_name <-- name associated with current group.
 *
 * returns:
 *   id associated with local name.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_name_to_id(MPI_Comm     comm,
                            const char  *group_name);

/*----------------------------------------------------------------------------
 * Discover other applications in a set with a common communicator.
 *
 * In most cases, the base communicator is MPI_COMM_WORLD, and the local
 * application communicator app_comm is usually obtained from it using
 * MPI_Comm_split, but other combinations may be possible using MPI-2
 * process management functions.
 *
 * As this function requires communication between applications, it
 * is a collective function in base_comm.
 *
 * parameters:
 *   sync_flag <-- 1 if application is to be synchronized at each
 *                 time step, 0 if independent from others.
 *   app_type  <-- name of current application type (software name).
 *   app_name  <-- name of current application (data/case name).
 *   base_comm <-- communicator associated with all applications.
 *   app_comm  <-- communicator associated with local application.
 *
 * returns:
 *   PLE coupling MPI set info structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_set_t *
ple_coupling_mpi_set_create(int          sync_flag,
                            const char  *app_type,
                            const char  *app_name,
                            MPI_Comm     base_comm,
                            MPI_Comm     app_comm);

/*----------------------------------------------------------------------------
 * Free an PLE coupling MPI set info structure.
 *
 * parameters:
 *   s <-> pointer to structure that should be freed.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_destroy(ple_coupling_mpi_set_t **s);

/*----------------------------------------------------------------------------
 * Return the number of applications in a coupled set.
 *
 * parameters:
 *   s <-- pointer to PLE coupling MPI set info structure.
 *
 * returns:
 *   number of application in set's common communicator.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_set_n_apps(const ple_coupling_mpi_set_t  *s);

/*----------------------------------------------------------------------------
 * Return the id of the local application in a coupled set.
 *
 * parameters:
 *   s <-- pointer to PLE coupling MPI set info structure.
 *
 * returns:
 *   id of the local application in set's common communicator.
 *----------------------------------------------------------------------------*/

int
ple_coupling_mpi_set_get_app_id(const ple_coupling_mpi_set_t  *s);

/*----------------------------------------------------------------------------
 * Return application information in set's common communicator.
 *
 * parameters:
 *   s      <-- pointer to PLE coupling MPI set info structure.
 *   app_id <-- application id
 *
 * returns:
 *   application information structure.
 *----------------------------------------------------------------------------*/

ple_coupling_mpi_set_info_t
ple_coupling_mpi_set_get_info(const ple_coupling_mpi_set_t  *s,
                              int                            app_id);

/*----------------------------------------------------------------------------
 * Synchronize applications in a set.
 *
 * Note that if a member of the set has used a PLE_COUPLING_STOP or
 * PLE_COUPLING_LAST flag when calling ple_coupling_mpi_set_create() or
 * or at the previous call to this function, it will not be synchronized
 * anymore (i.e. the PLE_COUPLING_NO_SYNC flag will be added).
 *
 * parameters:
 *   s         <-- pointer to PLE coupling MPI set info structure.
 *   sync_flag <-- synchronization info for current application.
 *   time_step <-- time step for current application.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_synchronize(ple_coupling_mpi_set_t  *s,
                                 int                      sync_flag,
                                 double                   time_step);

/*----------------------------------------------------------------------------
 * Get status of applications in a set.
 *
 * This function allows access to the status flag of each synchronized
 * application in the set. It may be used immediately after
 * ple_coupling_mpi_set_create(), and flags are updated after each
 * call to ple_coupling_mpi_set_synchronize().
 *
 * parameters:
 *   s <-- pointer to PLE coupling MPI set info structure.
 *
 * returns:
 *   a pointer to the set's array of status flags
 *----------------------------------------------------------------------------*/

const int *
ple_coupling_mpi_set_get_status(const ple_coupling_mpi_set_t  *s);

/*----------------------------------------------------------------------------
 * Get time steps in a set.
 *
 * This function may be called after ple_coupling_mpi_set_synchronize()
 * to access the time step values of each synchronized application in the set.
 *
 * parameters:
 *   s <-- pointer to PLE coupling MPI set info structure.
 *
 * returns:
 *   a pointer to the set's array of time steps
 *----------------------------------------------------------------------------*/

const double *
ple_coupling_mpi_set_get_timestep(const ple_coupling_mpi_set_t  *s);

/*----------------------------------------------------------------------------
 * Create an intracommunicator from a local and distant communicator
 * within a base communicator.
 *
 * Note that if a member of the set has used a PLE_COUPLING_STOP or
 * PLE_COUPLING_LAST flag when calling ple_coupling_mpi_set_create() or
 * or at the previous call to this function, it will not be synchronized
 * anymore (i.e. the PLE_COUPLING_NO_SYNC flag will be added).
 *
 * parameters:
 *   base_comm     <-- communicator associated with both applications
 *   app_comm      <-- communicator associated with local application
 *   distant_root  <-- rank of distant group leader in base_comm
 *   new_comm      --> pointer to new communicator
 *   local_range   --> first and past-the last ranks of local application
 *                     in new communicator
 *   distant_range --> first and past-the last ranks of distant application
 *                     in new communicator
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_intracomm_create(MPI_Comm   base_comm,
                                  MPI_Comm   app_comm,
                                  int        distant_root,
                                  MPI_Comm  *new_comm,
                                  int        local_range[2],
                                  int        distant_range[2]);

/*----------------------------------------------------------------------------
 * Dump printout of an PLE coupling MPI set info structure.
 *
 * parameters:
 *   w <-- pointer to PLE coupling MPI set info structure.
 *----------------------------------------------------------------------------*/

void
ple_coupling_mpi_set_dump(const ple_coupling_mpi_set_t  *s);

#endif /* defined(PLE_HAVE_MPI) */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PLE_COUPLING_H__ */
