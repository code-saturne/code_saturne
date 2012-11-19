/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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

/*============================================================================
 * Common functionnality for various coupling types.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/


/*============================================================================
 *  Global variables
 *============================================================================*/

#if defined(FVM_HAVE_MPI)

static int                       _cs_glob_coupling_mpi_app_num = -1;
static fvm_coupling_mpi_world_t *_cs_glob_coupling_mpi_app_world = NULL;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI root communicator.
 *
 * parameters:
 *   app_num  <-- application number for this instance of Code_Saturne (>= 0)
 *----------------------------------------------------------------------------*/

void
cs_coupling_discover_mpi_apps(int  app_num)
{
  if (app_num > -1 && cs_glob_mpi_comm != MPI_COMM_WORLD) {

    int i, n_apps, app_id;

    /* App_type contains a string such as
       "Code_Saturne 2.0.0" or "NEPTUNE_CFD 1.2.1" */

    const char app_type[] = CS_APP_NAME " " CS_APP_VERSION;

    _cs_glob_coupling_mpi_app_num = app_num;

    if (cs_glob_rank_id < 1) {
      bft_printf(_("\n"
                   "Applications accessible through MPI:\n"
                   "------------------------------------\n\n"));
      bft_printf_flush();
    }

    _cs_glob_coupling_mpi_app_world
      = fvm_coupling_mpi_world_create(app_num,
                                      app_type,
                                      NULL,
                                      cs_glob_mpi_comm);

    n_apps = fvm_coupling_mpi_world_n_apps(_cs_glob_coupling_mpi_app_world);
    app_id = fvm_coupling_mpi_world_get_app_id(_cs_glob_coupling_mpi_app_world);

    if (cs_glob_rank_id < 1) {

      const char local_add[] = " (this instance)";
      const char nolocal_add[] = "";

      for (i = 0; i< n_apps; i++) {
        const char *is_local = nolocal_add;
        fvm_coupling_mpi_world_info_t
          ai = fvm_coupling_mpi_world_get_info(_cs_glob_coupling_mpi_app_world,
                                               i);
        if (i == app_id)
          is_local = local_add;
        bft_printf(_("  %d; type:      \"%s\"%s\n"
                     "     case name: \"%s\"\n"
                     "     lead rank: %d; n_ranks: %d\n\n"),
                   ai.app_num, ai.app_type, is_local,
                   ai.app_name, ai.root_rank, ai.n_ranks);
      }

      bft_printf_flush();
    }
  }
}

/*----------------------------------------------------------------------------
 * Finalize MPI coupling helper structures.
 *----------------------------------------------------------------------------*/

void
cs_coupling_finalize(void)
{
  if (_cs_glob_coupling_mpi_app_world != NULL)
    fvm_coupling_mpi_world_destroy(&_cs_glob_coupling_mpi_app_world);
}

/*----------------------------------------------------------------------------
 * Return info on other applications in the same MPI root communicator.
 *
 * returns:
 *   info on other applications structure.
 *----------------------------------------------------------------------------*/

const fvm_coupling_mpi_world_t *
cs_coupling_get_mpi_apps(void)
{
  return _cs_glob_coupling_mpi_app_world;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/

END_C_DECLS
