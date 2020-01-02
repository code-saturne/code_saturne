/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_geom.h"
#include "cs_interpolate.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_probe.h"
#include "cs_time_plot.h"

#include "cs_field_pointer.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/
/*----------------------------------------------------------------------------
 * Example function for selection of cells with scalar field values above
 * a certain threshold.
 *
 * In this case, the selection is base on the value of the temperaturefield
 * being above 21 degrees.
 *
 * parameters:
 *   input    <-> pointer to input (unused here)
 *   n_cells  --> number of selected cells
 *   cell_ids --> array of selected cell ids (0 to n-1 numbering)
 *----------------------------------------------------------------------------*/

static void
_t_lt_21_select(void        *input,
                cs_lnum_t   *n_cells,
                cs_lnum_t  **cell_ids)
{
  cs_lnum_t i;
  cs_lnum_t _n_cells = 0;
  cs_lnum_t *_cell_ids = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  cs_field_t *f = cs_field_by_name("temperature"); /* Get access to field */

  if (f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "No field with name \"temperature\" defined");

  /* Before time loop, field is defined, but has no values yet,
     so ignore that case (postprocessing mesh will be initially empty) */

  if (f->val != NULL) {

    BFT_MALLOC(_cell_ids, m->n_cells, cs_lnum_t); /* Allocate selection list */

    for (i = 0; i < m->n_cells; i++) {
      if (f->val[i] < 21) {
        _cell_ids[_n_cells] = i;
        _n_cells += 1;
      }
    }

    BFT_REALLOC(_cell_ids, _n_cells, cs_lnum_t); /* Adjust size (good practice,
                                                    but not required) */

  }

  /* Set return values */

  *n_cells = _n_cells;
  *cell_ids = _cell_ids;
}
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

void
cs_user_postprocess_writers(void)
{
  /* Define additional writers */

  cs_post_define_writer(1,                            /* writer_id */
                        "user",                       /* writer name */
                        "postprocessing",             /* directory name */
                        "EnSight Gold",               /* format_name */
                        "",                           /* format_options */
                        FVM_WRITER_TRANSIENT_CONNECT, /* time dependency */
                        false,                        /* output_at_start */
                        true,                         /* output_at_end */
                        5,                            /* frequency_n */
                        -1);                          /* frequency_t */

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

void
cs_user_postprocess_meshes(void)
{
  /* Build a (time varying) volume mesh containing cells
     with values of field named "temperature" > < 21 */

  const int n_writers = 1;
  const int writer_ids[] = {1};  /* Associate to writer 1 */

  /* Define postprocessing mesh */

  cs_post_define_volume_mesh_by_func(1,               /* mesh id */
                                     "T_lt_21",
                                     _t_lt_21_select,
                                     NULL,            /* _t_lt_21_select_input */
                                     true,            /* time varying */
                                     false,           /* add_groups */
                                     false,           /* auto_variables */
                                     n_writers,
                                     writer_ids);
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
  /* Output of the temperature on output mesh 1 */

  if (mesh_id == 1) {

    cs_field_t *f = cs_field_by_name("temperature"); /* Get access to field */

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      cs_field_get_label(f),          /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      true,                           /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      f->val,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
