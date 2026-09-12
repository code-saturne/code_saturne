/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <cstdio>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*============================================================================
 * User function definitions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*
 * User function for output of values on a post-processing mesh.
 *
 * Export the per-face IWF state and diagnostic quantities on boundary meshes.
 *
 * \param[in]  mesh_name    name of the output mesh for the current call
 * \param[in]  mesh_id      id of the output mesh for the current call
 * \param[in]  cat_id       category id of the output mesh for the current call
 * \param[in]  probes       associated probe set, or nullptr
 * \param[in]  n_cells      local number of cells of the post-processing mesh
 * \param[in]  n_i_faces    local number of interior faces
 * \param[in]  n_b_faces    local number of boundary faces
 * \param[in]  n_vertices   local number of vertices
 * \param[in]  cell_list    list of cells, or nullptr
 * \param[in]  i_face_list  list of interior faces, or nullptr
 * \param[in]  b_face_list  list of boundary faces, or nullptr
 * \param[in]  vertex_list  list of vertices, or nullptr
 * \param[in]  ts           time-step status structure, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_values(
  [[maybe_unused]] const char            *mesh_name,
  int                                    mesh_id,
  [[maybe_unused]] int                    cat_id,
  [[maybe_unused]] cs_probe_set_t        *probes,
  [[maybe_unused]] cs_lnum_t              n_cells,
  [[maybe_unused]] cs_lnum_t              n_i_faces,
  cs_lnum_t                               n_b_faces,
  [[maybe_unused]] cs_lnum_t              n_vertices,
  [[maybe_unused]] const cs_lnum_t        cell_list[],
  [[maybe_unused]] const cs_lnum_t        i_face_list[],
  const cs_lnum_t                         b_face_list[],
  [[maybe_unused]] const cs_lnum_t        vertex_list[],
  const cs_time_step_t                   *ts)
{
  if (n_b_faces == 0)
    return;

  const bool iwf_configured_locally
    = cs_wall_functions_iwf_n_faces() > 0;

  const char *names[] = {
    "IWF active",
    "IWF converged",
    "IWF h_cv",
    "IWF T_infinity",
    "IWF Ri",
    "IWF y_cross",
    "IWF u_tau",
    "IWF u_ref"
  };

  for (int var_id = 0; var_id < 8; var_id++) {
    cs_array<cs_real_t> values(n_b_faces);

    for (cs_lnum_t i = 0; i < n_b_faces; i++) {
      const cs_lnum_t face_id
        = (b_face_list != nullptr) ? b_face_list[i] : i;
      cs_wall_functions_iwf_result_t result = {};
      if (iwf_configured_locally)
        cs_wall_functions_iwf_get_result(face_id, &result);

      switch (var_id) {
      case 0:
        values[i] = result.active ? 1. : 0.;
        break;
      case 1:
        values[i] = result.converged ? 1. : 0.;
        break;
      case 2:
        values[i] = result.h_cv;
        break;
      case 3:
        values[i] = result.t_infinity;
        break;
      case 4:
        values[i] = result.ri;
        break;
      case 5:
        values[i] = result.y_cross;
        break;
      case 6:
        values[i] = result.u_tau;
        break;
      case 7:
        values[i] = result.u_ref;
        break;
      default:
        values[i] = 0.;
      }
    }

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      names[var_id],
                      1,
                      true,
                      false,
                      CS_POST_TYPE_cs_real_t,
                      nullptr,
                      nullptr,
                      values,
                      ts);
  }
}

END_C_DECLS
