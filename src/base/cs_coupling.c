/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

#include <fvm_nodal_extract.h>
#include <fvm_point_location.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

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

#if defined(PLE_HAVE_MPI)

static int                       _cs_glob_coupling_mpi_app_num = -1;
static ple_coupling_mpi_world_t *_cs_glob_coupling_mpi_app_world = NULL;

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
 *   app_name <-- optional name of this instance of Code_Saturne, or NULL.
 *----------------------------------------------------------------------------*/

void
cs_coupling_discover_mpi_apps(int          app_num,
                              const char  *app_name)
{
  if (app_num > -1 && cs_glob_mpi_comm != MPI_COMM_WORLD) {

    int i, n_apps, app_id;

    /* App_type contains a string such as
       "Code_Saturne 1.4.0" or "NEPTUNE_CFD 1.2.1" */

    const char app_type[] = CS_APP_NAME " " CS_APP_VERSION;

    _cs_glob_coupling_mpi_app_num = app_num;

    if (cs_glob_rank_id < 1) {
      bft_printf(_("\n"
                   "Applications accessible through MPI:\n"
                   "------------------------------------\n\n"));
      bft_printf_flush();
    }

    _cs_glob_coupling_mpi_app_world
      = ple_coupling_mpi_world_create(app_num,
                                      app_type,
                                      app_name,
                                      cs_glob_mpi_comm);

    n_apps = ple_coupling_mpi_world_n_apps(_cs_glob_coupling_mpi_app_world);
    app_id = ple_coupling_mpi_world_get_app_id(_cs_glob_coupling_mpi_app_world);

    if (cs_glob_rank_id < 1) {

      const char local_add[] = " (this instance)";
      const char nolocal_add[] = "";

      for (i = 0; i< n_apps; i++) {
        const char *is_local = nolocal_add;
        ple_coupling_mpi_world_info_t
          ai = ple_coupling_mpi_world_get_info(_cs_glob_coupling_mpi_app_world,
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
    ple_coupling_mpi_world_destroy(&_cs_glob_coupling_mpi_app_world);
}

/*----------------------------------------------------------------------------
 * Return info on other applications in the same MPI root communicator.
 *
 * returns:
 *   info on other applications structure.
 *----------------------------------------------------------------------------*/

const ple_coupling_mpi_world_t *
cs_coupling_get_mpi_apps(void)
{
  return _cs_glob_coupling_mpi_app_world;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Compute extents of a mesh representation
 *
 * parameters:
 *   mesh          <-- pointer to mesh representation structure
 *   n_max_extents <-- maximum number of sub-extents (such as element extents)
 *                     to compute, or -1 to query
 *   tolerance     <-- addition to local extents of each element:
 *                     extent = base_extent * (1 + tolerance)
 *   extents       <-> extents associated with mesh:
 *                     x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 * returns:
 *   the number of extents computed
 *----------------------------------------------------------------------------*/

ple_lnum_t
cs_coupling_mesh_extents(const void  *mesh,
                         ple_lnum_t   n_max_extents,
                         double       tolerance,
                         double       extents[])
{
  const fvm_nodal_t  *m = mesh;
  ple_lnum_t retval = 0;

  if (m == NULL)
    return 0;

  /* In query mode, return maximum extents available
     (currently limited to 1) */

  if (n_max_extents < 0)
    retval = 1;

  /* If n_max_extents > 0 return global mesh extents */

  else if (n_max_extents > 0) {
    fvm_nodal_extents(m, tolerance, extents);
    retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to the id of a given element + 1 in
 * concatenated sections of same element dimension.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   tolerance    <-- associated tolerance
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, 0 - 1 if inside,
 *                    and > 1 if outside a volume element, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh(const void         *mesh,
                          double              tolerance,
                          ple_lnum_t          n_points,
                          const ple_coord_t   point_coords[],
                          ple_lnum_t          location[],
                          float               distance[])
{
  fvm_point_location_nodal((const fvm_nodal_t *)mesh,
                           tolerance,
                           0, /* Do not locate on parents */
                           n_points,
                           point_coords,
                           location,
                           distance);
}

/*----------------------------------------------------------------------------
 * Find elements in a given mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * Location is relative to parent element numbers.
 *
 * parameters:
 *   mesh         <-- pointer to mesh representation structure
 *   tolerance    <-- associated tolerance
 *   n_points     <-- number of points to locate
 *   point_coords <-- point coordinates
 *   location     <-> number of element containing or closest to each
 *                    point (size: n_points)
 *   distance     <-> distance from point to element indicated by
 *                    location[]: < 0 if unlocated, 0 - 1 if inside,
 *                    and > 1 if outside a volume element, or absolute
 *                    distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
cs_coupling_point_in_mesh_p(const void         *mesh,
                            double              tolerance,
                            ple_lnum_t          n_points,
                            const ple_coord_t   point_coords[],
                            ple_lnum_t          location[],
                            float               distance[])
{
  fvm_point_location_nodal((const fvm_nodal_t *)mesh,
                           tolerance,
                           1, /* Locate on parents */
                           n_points,
                           point_coords,
                           location,
                           distance);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
