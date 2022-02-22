/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_writers
void
cs_user_postprocess_writers(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_meshes
void
cs_user_postprocess_meshes(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_probes
void
cs_user_postprocess_probes(void)
{

}

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

#pragma weak cs_user_postprocess_values
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
  CS_UNUSED(mesh_id);
  CS_UNUSED(cat_id);
  CS_UNUSED(probes);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(n_vertices);
  CS_UNUSED(cell_list);
  CS_UNUSED(i_face_list);
  CS_UNUSED(b_face_list);
  CS_UNUSED(vertex_list);
  CS_UNUSED(ts);
}

/*----------------------------------------------------------------------------*/
/*!
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param[in]  nt_max_abs  maximum time step number
 * \param[in]  nt_cur_abs  current time step number
 * \param[in]  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_activate
void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{
  CS_UNUSED(nt_max_abs);
  CS_UNUSED(nt_cur_abs);
  CS_UNUSED(t_cur_abs);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
