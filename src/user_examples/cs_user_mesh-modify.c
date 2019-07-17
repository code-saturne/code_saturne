/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh modification function examples.
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

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.c
 *
 * \brief Mesh modification example.
 *
 * See \subpage cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  /* Example: modify vertex coordinates */
  /*------------------------------------*/

  /* Divide coordinates by 1000 (millimetres to metres).
   *
   * Warning:
   *
   *   This is incompatible with pre-processed periodicity,
   *   as the periodicity transformation is not updated.
   *
   *   With periodicity, using a coordinate transformation matrix
   *   in cs_user_mesh_input is preferred. */

  /*! [mesh_modify_coords] */
  {
    cs_lnum_t  vtx_id;
    const double  coo_mult = 1. / 1000.;

    for (vtx_id = 0; vtx_id < mesh->n_vertices; vtx_id++) {
      mesh->vtx_coord[vtx_id*3]     *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 1] *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 2] *= coo_mult;
    }

    /* Set mesh modification flag if it should be saved for future re-use. */

    mesh->modified = 1;
  }
  /*! [mesh_modify_coords] */

  /* Extrude mesh at boundary faces of group "outlet".
     We use a regular extrusion here */

  /*! [mesh_modify_extrude_1] */
  {
    int n_layers = 2;
    double thickness = 1.0;
    double expansion_factor = 1.5;

    const char criteria[] = "outlet";

    /* Select boudary faces */

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;

    BFT_MALLOC(selected_faces, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    /* Extrude selected boundary */

    cs_mesh_extrude_constant(mesh,
                             false, /* Maintain groups on interio faces? */
                             n_layers,
                             thickness,
                             expansion_factor,
                             n_selected_faces,
                             selected_faces);

    /* Free temporary memory */

    BFT_FREE(selected_faces);

  }
  /*! [mesh_modify_extrude_1] */

  /* Advanced mesh extrusion: impose a direction
   *
   */

  /*! [mesh_modify_extrude_2] */
  {
    /* Define extrusion parameters for each face */

    int n_zones = 2;
    const char *sel_criteria[] = {"wall_1", "wall_2"};
    const int zone_layers[] = {2, 4};
    const double zone_thickness[] = {0.1, 0.3};
    const float zone_expansion[] = {0.8, 0.7};

    cs_mesh_extrude_face_info_t *efi = cs_mesh_extrude_face_info_create(mesh);

    cs_lnum_t n_faces;
    cs_lnum_t *face_list;

    BFT_MALLOC(face_list, mesh->n_b_faces, cs_lnum_t);

    for (int z_id = 0; z_id < n_zones; z_id++) {

      cs_selector_get_b_face_list(sel_criteria[z_id], &n_faces, face_list);

      cs_mesh_extrude_set_info_by_zone(efi,
                                       zone_layers[z_id],
                                       zone_thickness[z_id],
                                       zone_expansion[z_id],
                                       n_faces,
                                       face_list);

    }

    BFT_FREE(face_list);

    /* Determine vertex values for extrusion */

    cs_mesh_extrude_vectors_t *e = cs_mesh_extrude_vectors_create(efi);

    /* Overwrite the total coord_shift */
    cs_real_3_t *vtx_coords = (cs_real_3_t *)mesh->vtx_coord;

    for (cs_lnum_t id = 0; id < e->n_vertices; id++) {

      cs_lnum_t v_id = e->vertex_ids[id];

      e->coord_shift[id][0] = 0.;
      e->coord_shift[id][1] = 0.;
      /* Percentage of the original z coordinate */
      e->coord_shift[id][2] = 0.01 * vtx_coords[v_id][2];
    }

    /* Extrude mesh with this  */

    cs_mesh_extrude_face_info_destroy(&efi);

    cs_mesh_extrude(mesh,
                    e,
                    /* Maintain group classes of  interior faces previously on boundary */
                    true);

    cs_mesh_extrude_vectors_destroy(&e);
  }
  /*! [mesh_modify_extrude_2] */

  /* Add a group to cells in a given region */

  /*! [mesh_modify_groups_1] */
  {
    cs_lnum_t   n_selected_elts = 0;
    cs_lnum_t  *selected_elts = NULL;

    const char criteria[] = "box[0.5, 0.5, 0, 1, 1, 0.05]";

    BFT_MALLOC(selected_elts, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_elts,
                              selected_elts);

    cs_mesh_group_cells_add(mesh,
                            "source_region",
                            n_selected_elts,
                            selected_elts);

    BFT_FREE(selected_elts);

    /* Mark mesh as modified to save it */

    mesh->modified = 1;
  }
  /*! [mesh_modify_groups_1] */

  /* Add a group to boundary faces in a given region */

  /*! [mesh_modify_groups_2] */
  {
    cs_lnum_t   n_selected_elts = 0;
    cs_lnum_t  *selected_elts = NULL;

    const char criteria[] = "box[0.5, 0.5, 0, 1, 1, 0.05]";

    BFT_MALLOC(selected_elts, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_elts,
                                selected_elts);

    cs_mesh_group_b_faces_add(mesh,
                              "source_region",
                              n_selected_elts,
                              selected_elts);

    BFT_FREE(selected_elts);

    /* Mark mesh as modified to save it */

    mesh->modified = 1;
  }
  /*! [mesh_modify_groups_2] */


  /* Insert boundary layers on selected zones.
   *
   * Warning:
   *
   *   This is incompatible with pre-processed periodicity,
   *   as the periodicity transformation is not updated.
   *
   *   With periodicity, using a coordinate transformation matrix
   *   in cs_user_mesh_input is preferred. */

  /*! [mesh_modify_boundary_layer] */
  {
    /* Define extrusion parameters for each face */

    int n_zones = 2;
    const char *sel_criteria[] = {"wall_1", "wall_2"};
    const int zone_layers[] = {2, 4};
    const double zone_thickness[] = {0.1, 0.3};
    const float zone_expansion[] = {0.8, 0.7};

    cs_mesh_extrude_face_info_t *efi = cs_mesh_extrude_face_info_create(mesh);

    cs_lnum_t n_faces;
    cs_lnum_t *face_list;

    BFT_MALLOC(face_list, mesh->n_b_faces, cs_lnum_t);

    for (int z_id = 0; z_id < n_zones; z_id++) {

      cs_selector_get_b_face_list(sel_criteria[z_id], &n_faces, face_list);

      cs_mesh_extrude_set_info_by_zone(efi,
                                       zone_layers[z_id],
                                       zone_thickness[z_id],
                                       zone_expansion[z_id],
                                       n_faces,
                                       face_list);

    }

    BFT_FREE(face_list);

    /* Determine vertex values for extrusion */

    cs_mesh_extrude_vectors_t *e = cs_mesh_extrude_vectors_create(efi);

    /* Insert boundary layer */

    cs_mesh_extrude_face_info_destroy(&efi);

    cs_mesh_boundary_layer_insert(mesh, e, 0.2, false, 0, NULL);

    cs_mesh_extrude_vectors_destroy(&e);
  }
  /*! [mesh_modify_boundary_layer] */

  /* Refine a selected portion of a mesh */

  /*! [mesh_modify_refine_1] */
  {
    const char criteria[] = "box[0, 0, 0, 0.5, 0.5, 0.5]";

    cs_lnum_t   n_selected_cells = 0;
    cs_lnum_t  *selected_cells = NULL;

    BFT_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_cells,
                              selected_cells);

    cs_mesh_refine_simple_selected(mesh,
                                   true,              /* conforming or not */
                                   n_selected_cells,
                                   selected_cells);

    BFT_FREE(selected_cells);
  }
  /*! [mesh_modify_refine_1] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
