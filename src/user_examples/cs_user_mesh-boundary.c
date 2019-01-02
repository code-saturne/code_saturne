/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between Code_Saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundaries into a mesh.
 *
 * \param[in,out] mesh pointer to cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_boundary(cs_mesh_t  *mesh)
{
  /* Example: insert boundary along a plane */
  /*----------------------------------------*/

  /*! [mesh_thinwall] */
  {
    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;

    /* example of multi-line character string */

    const char criteria[] = "plane[0, -1, 0, 0.5, epsilon = 0.0001]"
                            " or plane[-1, 0, 0, 0.5, epsilon = 0.0001]";

    BFT_MALLOC(selected_faces, mesh->n_i_faces, cs_lnum_t);

    cs_selector_get_i_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    cs_mesh_boundary_insert(mesh,
                            n_selected_faces,
                            selected_faces);

    BFT_FREE(selected_faces);
  }
  /*! [mesh_thinwall] */

  /* Example: boundary separating cell groups */
  /*------------------------------------------*/

  /*! [mesh_boundary_cells] */
  {
    cs_lnum_t   n_selected_cells = 0;
    cs_lnum_t  *selected_cells = NULL;

    const char criteria[] = "box[0.5, 0.5, 0, 1, 1, 0.05]";

    BFT_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_cells,
                              selected_cells);

    cs_mesh_boundary_insert_separating_cells(mesh,
                                             "zone_interface",
                                             n_selected_cells,
                                             selected_cells);

    BFT_FREE(selected_cells);
  }
  /*! [mesh_boundary_cells] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
