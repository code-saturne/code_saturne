/*============================================================================
 * Management of the post-processing
 *============================================================================*/

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

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_ctwr.h"
#include "cs_field.h"
#include "cs_function.h"
#include "cs_function_default.h"
#include "cs_lagr_tracking.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"

#include "cs_post.h"
#include "cs_post_util.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_post_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local types and structures
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Default output format and options */

static bool  _default_functions_are_registered = false;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Default additional output of mesh and time-dependent variables for the
 * call to cs_post_write_vars.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the output function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_var()
 *   or similar before the data pointed to goes out of scope.
 *
 * parameters:
 *   input       <-> pointer to optional (untyped) value or structure;
 *                   unused here.
 *   mesh_id     <-- id of the output mesh for the current call
 *   cat_id      <-- category id of the output mesh for the current call
 *   ent_flag    <-- indicate global presence of cells (ent_flag[0]), interior
 *                   faces (ent_flag[1]), boundary faces (ent_flag[2]),
 *                   particles (ent_flag[3]) or probes (ent_flag[4])
 *   n_cells     <-- local number of cells of post_mesh
 *   n_i_faces   <-- local number of interior faces of post_mesh
 *   n_b_faces   <-- local number of boundary faces of post_mesh
 *   cell_ids    <-- list of cells (0 to n-1) of post-processing mesh
 *   i_face_ids  <-- list of interior faces (0 to n-1) of post-processing mesh
 *   b_face_ids  <-- list of boundary faces (0 to n-1) of post-processing mesh
 *   ts          <-- time step status structure
 *----------------------------------------------------------------------------*/

static void
_write_additional_vars(void                  *input,
                       int                    mesh_id,
                       int                    cat_id,
                       int                    ent_flag[5],
                       cs_lnum_t              n_cells,
                       cs_lnum_t              n_i_faces,
                       cs_lnum_t              n_b_faces,
                       const cs_lnum_t        cell_ids[],
                       const cs_lnum_t        i_face_ids[],
                       const cs_lnum_t        b_face_ids[],
                       const cs_time_step_t  *ts)
{
  CS_UNUSED(input);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(cell_ids);
  CS_UNUSED(i_face_ids);

  if (cat_id == CS_POST_MESH_BOUNDARY && ent_flag[1] == 0) {

    const int n_fields = cs_field_n_fields();
    const int vis_key_id = cs_field_key_id("post_vis");

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t  *f = cs_field_by_id(f_id);

      if (f->location_id != CS_MESH_LOCATION_CELLS)
        continue;

      if (! (cs_field_get_key_int(f, vis_key_id) & CS_POST_BOUNDARY_NR))
        continue;

      cs_real_t  *b_face_val = NULL;
      BFT_MALLOC(b_face_val,
                 n_b_faces * (cs_lnum_t)(f->dim),
                 cs_real_t);

      cs_function_field_boundary_nr(f->location_id,
                                    n_b_faces,
                                    b_face_ids,
                                    f,
                                    b_face_val);

      char name[80];
      snprintf(name, 79, "bc_%s", cs_field_get_label(f)); name[78] = '\0';

      cs_post_write_var(mesh_id,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        name,
                        f->dim,
                        true,
                        false, /* use_parent */
                        CS_POST_TYPE_cs_real_t,
                        NULL,
                        NULL,
                        b_face_val,
                        ts);

      BFT_FREE(b_face_val);

    }

  }

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *----------------------------------------------------------------------------*/

void
cs_post_default_write_meshes(void)
{
  int t_top_id
    = cs_timer_stats_switch(cs_timer_stats_id_by_name("postprocessing_stage"));

  cs_post_write_meshes(cs_glob_time_step);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *----------------------------------------------------------------------------*/

void
cs_post_default_write_variables(void)
{
  /* Register function for first pass */

  if (_default_functions_are_registered == false) {
    cs_post_add_time_mesh_dep_output(_write_additional_vars, NULL);
    _default_functions_are_registered = true;
  }

  /* Call main post-processing function */

  if (cs_glob_time_step->nt_cur > -1)
    cs_post_write_vars(cs_glob_time_step);
  else
    cs_post_write_vars(NULL);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
