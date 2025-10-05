/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, null otherwise
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
 * \param[in]       ts           time step status structure, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values
(
  [[maybe_unused]] const char            *mesh_name,
  [[maybe_unused]] int                    mesh_id,
  [[maybe_unused]] int                    cat_id,
  [[maybe_unused]] cs_probe_set_t        *probes,
  [[maybe_unused]] cs_lnum_t              n_cells,
  [[maybe_unused]] cs_lnum_t              n_i_faces,
  [[maybe_unused]] cs_lnum_t              n_b_faces,
  [[maybe_unused]] cs_lnum_t              n_vertices,
  [[maybe_unused]] const cs_lnum_t        cell_list[],
  [[maybe_unused]] const cs_lnum_t        i_face_list[],
  [[maybe_unused]] const cs_lnum_t        b_face_list[],
  [[maybe_unused]] const cs_lnum_t        vertex_list[],
  [[maybe_unused]] const cs_time_step_t  *ts
)
{
  if (probes != nullptr) {

/*! [post_probes_interpolate_var_1] */
    int n_p_fields = 2;
    const char *p_field_names[] = {"velocity", "temperature"};

    for (int i = 0; i < n_p_fields; i++) {

      cs_field_t *f = cs_field_by_name_try(p_field_names[i]);

      if (f != nullptr) {

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
