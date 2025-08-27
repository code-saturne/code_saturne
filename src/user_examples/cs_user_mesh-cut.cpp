/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "mesh/cs_mesh_cut.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Cut mesh cells by their associated planes.
 *
 * A plane is defined as a normal + an origin.
 *
 * The user must provide the plane normals and origins. In this example,
 * they are read from a checkpoint file, after having run the
 * cs_porosity_from_scan algorithm on a custom .pts file.
 *
 * After this modification, a joining operation should be called on the group,
 * "auto:transformed_internal_faces".
 *
 * The boundary polygons created by the cuts are bundled in the group,
 * "auto:closing_polygons".
 *
 * The cells lying in the positive/negative half-spaces defined by the plane
 * normals are bundled in the group, "auto::positve_cells"/
 * "auto::negative_cells".
 *
 * These group informations are useful for further mesh processing.
 */
/*----------------------------------------------------------------------------*/

static int
_read_plane_data(const cs_mesh_t *mesh,
                 cs_real_t *p_normals[],
                 cs_real_t *p_origins[],
                 const char *file_name,
                 const char *absolute_file_path,
                 const char *var_normals_name,
                 const char *var_origins_name)
{
  cs_restart_t *restart;
  int n_values_per_location;

  restart = cs_restart_create(file_name,
                              absolute_file_path,
                              CS_RESTART_MODE_READ);

  n_values_per_location = 3;
  *p_normals = NULL;
  CS_MALLOC(*p_normals, n_values_per_location*mesh->n_cells,
      cs_real_t);

  if (CS_RESTART_SUCCESS != cs_restart_read_section(restart, var_normals_name,
      1, // location at cells
      n_values_per_location,
      CS_TYPE_cs_real_t,
      *p_normals)) {
    fprintf(stderr, "Error reading plane normals!\n");
    CS_FREE(p_normals);
    cs_restart_destroy(&restart);
    return 1;
  }

  *p_origins = NULL;
  CS_MALLOC(*p_origins, n_values_per_location*mesh->n_cells,
      cs_real_t);

  if (CS_RESTART_SUCCESS != cs_restart_read_section(restart, var_origins_name,
      1, // location at cells
      n_values_per_location,
      CS_TYPE_cs_real_t,
      *p_origins)) {
    fprintf(stderr, "Error reading plane origins!\n");
    CS_FREE(p_normals);
    CS_FREE(p_origins);
    cs_restart_destroy(&restart);
    return 1;
  }

  cs_restart_destroy(&restart);
  return 0;
}

void
cs_user_mesh_modify(cs_mesh_t *mesh)
{
  cs_real_t *p_normals = nullptr, *p_origins = nullptr;

  /* Modify the input strings here. */
  const char *file_name = "file_name.csc";
  const char *absolute_file_path = "/home/...";
  const char *var_normals_name = "plane_normals";
  const char *var_origins_name = "plane_origins";

  int ret = _read_plane_data(mesh,
                             &p_normals, &p_origins,
                             file_name,
                             absolute_file_path,
                             var_normals_name, var_origins_name);

  if (ret != 0) return;

  cs_mesh_cut(mesh, (cs_real_3_t *)p_normals, (cs_real_3_t *)p_origins);

  CS_FREE(p_normals);
  CS_FREE(p_origins);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
