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

#ifndef __CS_COUPLING_H__
#define __CS_COUPLING_H__

/*============================================================================
 * Common functionnality for various coupling types.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI root communicator.
 *
 * parameters:
 *   app_num  <-- application number for this instance of Code_Saturne (>= 0)
 *   app_name <-- optional name of this instance of Code_Saturne, or NULL.
 *----------------------------------------------------------------------------*/

void
cs_coupling_discover_mpi_apps(int          app_num,
                              const char  *app_name);

/*----------------------------------------------------------------------------
 * Finalize MPI coupling helper structures.
 *----------------------------------------------------------------------------*/

void
cs_coupling_finalize(void);

/*----------------------------------------------------------------------------
 * Return info on other applications in the same MPI root communicator.
 *
 * returns:
 *   info on other applications structure.
 *----------------------------------------------------------------------------*/

const fvm_coupling_mpi_world_t *
cs_coupling_get_mpi_apps(void);

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_COUPLING_H__ */
