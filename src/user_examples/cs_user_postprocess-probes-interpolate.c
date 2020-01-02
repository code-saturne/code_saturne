/*============================================================================
 * Define postprocessing output.
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

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, NULL otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts)
{
  CS_UNUSED(mesh_name);
  CS_UNUSED(cat_id);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(n_vertices);
  CS_UNUSED(cell_list);
  CS_UNUSED(i_face_list);
  CS_UNUSED(b_face_list);
  CS_UNUSED(vertex_list);

  if (probes != NULL) {

/*! [post_probes_interpolate_var_1] */
    int n_p_fields = 2;
    const char *p_field_names[] = {"velocity", "temperature"};

    for (int i = 0; i < n_p_fields; i++) {

      cs_field_t *f = cs_field_by_name_try(p_field_names[i]);

      if (f != NULL) {

        /* use different name to avoid conflict with field name in case already
           present in probe set through default output */

        char f_name[128], p_name[128];
        snprintf(f_name, 63, "%s", f->name); f_name[127] = '\0';
        snprintf(p_name, 63, "%s_p", f->name); p_name[63] = '\0';

        cs_post_write_probe_values
          (mesh_id,
           CS_POST_WRITER_ALL_ASSOCIATED,    /* writer id filter */
           p_name,                           /* var_name */
           f->dim,                           /* var_dim */
           CS_POST_TYPE_cs_real_t,
           1,                                /* parent location id */
           cs_interpolate_from_location_p1,  /* P1 interpolation */
           f_name,                           /* interpolation input */
           f->val,
           ts);

      }

    }

  }
  /*! [post_probes_interpolate_var_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
