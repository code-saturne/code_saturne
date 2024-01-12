/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh modification function examples.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.c
 *
 * \brief Mesh modification example.
 *
 * See \ref cs_user_mesh for examples.
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
    const double  coo_mult = 1. / 1000.;

    for (cs_lnum_t vtx_id = 0; vtx_id < mesh->n_vertices; vtx_id++) {
      mesh->vtx_coord[vtx_id*3]     *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 1] *= coo_mult;
      mesh->vtx_coord[vtx_id*3 + 2] *= coo_mult;
    }

    /* Set mesh modification flag if it should be saved for future re-use. */

    mesh->modified |= CS_MESH_MODIFIED;
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

    /* Select boundary faces */

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;

    BFT_MALLOC(selected_faces, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    /* Extrude selected boundary */

    cs_mesh_extrude_constant(mesh,
                             false, /* Maintain groups on interior faces? */
                             n_layers,
                             thickness,
                             expansion_factor,
                             n_selected_faces,
                             selected_faces);

    /* Free temporary memory */

    BFT_FREE(selected_faces);

  }
  /*! [mesh_modify_extrude_1] */

  /* Advanced mesh extrusion: impose a direction */

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

    /* Extrude mesh with this */

    cs_mesh_extrude_face_info_destroy(&efi);

    cs_mesh_extrude(mesh,
                    e,
                    true); /* Maintain group classes of interior
                              faces previously on boundary */

    cs_mesh_extrude_vectors_destroy(&e);
  }
  /*! [mesh_modify_extrude_2] */

  /* Extrude mesh at boundary faces of group "walls".
   * The resulting extruded cells are added to a new
   * group of cells called "solid" */

  /*! [mesh_modify_extrude_3] */
  {
    int n_layers = 2;
    double thickness = 1.0;
    double expansion_factor = 1.5;

    const char criteria[] = "walls";

    /* Save the initial number of cells */

    cs_lnum_t n_prev_cells = mesh->n_cells ;

    /* Select boundary faces */

    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;

    BFT_MALLOC(selected_faces, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list(criteria,
                                &n_selected_faces,
                                selected_faces);

    /* Extrude selected boundary */

    cs_mesh_extrude_constant(mesh,
                             false, /* Maintain groups on interior faces? */
                             n_layers,
                             thickness,
                             expansion_factor,
                             n_selected_faces,
                             selected_faces);

    /* Free temporary memory */

    BFT_FREE(selected_faces);

    /* Compute the number of extruded cells */

    cs_lnum_t n_selected_elts = mesh->n_cells - n_prev_cells ;

    /* Among all the cells, only select the cells above
     * the initial number of cells (before extrusion). */

    cs_lnum_t  *selected_elts = NULL;
    BFT_MALLOC(selected_elts, mesh->n_cells, cs_lnum_t);

    for(int i=0; i<n_selected_elts; i++)
      selected_elts[i] = n_prev_cells + i;

    /* Add selected cells to a new group called "solid" */
    cs_mesh_group_cells_add(mesh,
                            "solid",
                            n_selected_elts,
                            selected_elts);

    BFT_FREE(selected_elts);
  }
  /*! [mesh_modify_extrude_3] */

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

    mesh->modified |= CS_MESH_MODIFIED;
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

    mesh->modified |= CS_MESH_MODIFIED;
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

  /* Refine all cells which are intersected by a surface defined by an
   * STL file.
   */
  /*! [mesh_modify_refine_2] */
  {
    /* Create the cs_stl_mesh structure with a name "STLMESH1" */
    cs_stl_mesh_t *stl_mesh = cs_stl_mesh_add("STLMESH1");

    /* Define the stl file to read */
    cs_stl_file_read(stl_mesh,        /* pointer to cs_stl_mesh structure */
                     "cad_file.stl"); /* Name of the stl file to read. */

    /* Refine the mesh using the stl file */
    cs_stl_refine(stl_mesh, /* pointer to cs_stl_mesh_t structure */
                  3,        /* Number of refinement levels, here 3 */
                  2);       /* Propagate refinement criteria over 2 more layers */
  }
  /*! [mesh_modify_refine_2] */

  /* Remove cells from a selection
   * Note: if present, remove periodicity info first */
  /*! [mesh_modify_remove_cells] */
  {
    cs_lnum_t   n_selected_elts = 0;
    cs_lnum_t  *selected_elts = NULL;

    const char criteria[] = "box[-250, -250, 0, 250, 250, 100]";

    BFT_MALLOC(selected_elts, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
                              &n_selected_elts,
                              selected_elts);

    char *flag;
    BFT_MALLOC(flag, mesh->n_cells, char);

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      flag[i] = 0;
    }

    for (cs_lnum_t i = 0; i < n_selected_elts; i++) {
      flag[selected_elts[i]] = 1;
    }

    cs_mesh_remove_cells(mesh, flag, "[Building]");

    BFT_FREE(selected_elts);
    BFT_FREE(flag);

    /* Mark for re-partitioning */
    mesh->modified |= CS_MESH_MODIFIED_BALANCE;
  }
  /*! [mesh_modify_remove_cells] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply partial modifications to the mesh after the preprocessing
 *        stage, but before initial postprocessing mesh building.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 * \param[in,out] mesh_quantities pointer to a cs_mesh_quantities_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_mesh_modify_partial(cs_mesh_t             *mesh,
                            cs_mesh_quantities_t  *mesh_quantities)
{
  {
    /*! [mesh_modify_ignore_symmetry_faces] */
    cs_lnum_t   n_faces = 0;
    cs_lnum_t  *face_ids = NULL;

    BFT_MALLOC(face_ids, mesh->n_b_faces, cs_lnum_t);

    cs_selector_get_b_face_list("symmetry",
                                &n_faces,
                                face_ids);

    cs_preprocess_mesh_selected_b_faces_ignore(mesh,
                                               mesh_quantities,
                                               n_faces,
                                               face_ids);

    BFT_FREE(face_ids);
    /*! [mesh_modify_ignore_symmetry_faces] */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
