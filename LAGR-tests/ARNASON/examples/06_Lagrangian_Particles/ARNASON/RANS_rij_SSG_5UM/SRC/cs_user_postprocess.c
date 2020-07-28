/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* Code_Saturne version 5.1-alpha */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_probe.h"
#include "cs_time_plot.h"

#include "cs_field_pointer.h"
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a profile based on centers of cells cut by a vertical segment
 *
 * Here, the input points to string describing a segment's x coordinate.
 *
 * \param[in]   input   pointer to segment start and end:
 *                      [x0, y0, z0, x1, y1, z1]
 * \param[out]  n_elts  number of selected coordinates
 * \param[out]  coords  coordinates of selected elements.
 * \param[out]  s       curvilinear coordinates of selected elements
 *----------------------------------------------------------------------------*/

static void
_cell_x_profile_probes_define(void          *input,
                              cs_lnum_t     *n_elts,
                              cs_real_3_t  **coords,
                              cs_real_t    **s)
{
  const char *z_str = (const char *)input;

  const cs_real_t x0 = -45e-3, x1 = 45e-3;
  const cs_real_t y = 1e-6;
  const cs_real_t z = atof(z_str);

  cs_real_t seg[6] = {x0, y, z, x1, y, z};

  /* Call general cell selection function with adapted input */

  cs_cell_segment_intersect_probes_define(seg, n_elts, coords, s);
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
  /* redine profile output options */

  cs_post_define_writer(CS_POST_WRITER_PROFILES,  /* writer_id */
                        "",                       /* writer name */
                        "profiles",
                        "plot",                   /* format name */
                        "csv",                    /* format options */
                        FVM_WRITER_FIXED_MESH,
                        false,                    /* output_at_start */
                        true,                     /* output at end */
                        -1,                       /* time step frequency */
                        -1.0);                    /* time value frequency */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_probes(void)
{
  /* x vertical profiles; note we declare this array as "static" so that
     its members are garanteed to be present when the _cell_profile_select
     function is called */

  static const char *line_defs[]
    = {"-0.318", "Z318",
       "-0.502", "Z502",
       "-0.679", "Z679",
       "-1.320", "Z132"};

  for (int i = 0; i < 4; i++) {

    /* Define profiles */

    cs_probe_set_t *pset
      = cs_probe_set_create_from_local(line_defs[i*2+1],
                                       _cell_x_profile_probes_define,
                                       (void *)line_defs[i*2]);  /* input */

    /* Associate writers */

    const int writer_ids[] = {CS_POST_WRITER_PROFILES};
    cs_probe_set_associate_writers(pset, 1, writer_ids);

  }
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
  if (probes != NULL) {

    const cs_mesh_t *m = cs_glob_mesh;
    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

    const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

    int stat_type = cs_lagr_stat_type_from_attr_id(CS_LAGR_VELOCITY);

    cs_field_t *stat_vel
      = cs_lagr_stat_get_moment(stat_type,
                                CS_LAGR_MOMENT_MEAN,
                                0,
                                -1);

    cs_field_t *stat_vol_frac
      = cs_lagr_stat_get_moment(CS_LAGR_STAT_VOLUME_FRACTION,
                                CS_LAGR_MOMENT_MEAN,
                                0,
                                -1);

    const char *name = cs_probe_set_get_name(probes);

    const cs_time_step_t *ts_post = (ts->nt_cur == ts->nt_max) ? NULL : ts;

    /* For "x" profiles on fixed "z" planes
       ------------------------------------ */

    if (strncmp(name, "Z", 1) == 0) {

      cs_real_t *val;
      BFT_MALLOC(val, n_cells, cs_real_t);

      char var_name[64];

      /* Loop on columns */

      for (int col = 0; col < 6; col++) {

        switch(col) {

        case 0:
          {
            strncpy(var_name, "Z", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = cell_cen[c_id][2];
            }
          }
          break;

        case 1:
          {
            strncpy(var_name, "x/R", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              cs_real_t r = cell_cen[c_id][0];
              val[i] = r / 45e-3;
            }
          }
          break;

        case 2:
          {
            strncpy(var_name, "Vx", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = stat_vel->val[c_id*3];
            }
          }
          break;

        case 3:
        case 4:
          {
            if (col == 3)
              strncpy(var_name, "Vy/Vy_center", 64);
            else
              strncpy(var_name, "Vz/Vz_center", 64);

            int j = col - 2;

            /* get values and find point on profile closest to {0, 0, z} */

            cs_real_t min_d = HUGE_VAL;
            cs_real_t ref_v = 0;

            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = stat_vel->val[c_id*3 + j];
              cs_real_t d = sqrt(  cell_cen[c_id][0]*cell_cen[c_id][0]
                                 + cell_cen[c_id][1]*cell_cen[c_id][1]);
              if (d < min_d) {
                ref_v = val[i];
                min_d = d;
              }
            }

            cs_parall_min_loc_vals(1, &min_d, &ref_v);

            bft_printf("\nStat %s0 equal to %g for plane %s\n",
                       var_name, ref_v, name);

            /* now normalize by value at point on profile closest to {0, 0, z} */

            for (cs_lnum_t i = 0; i < n_cells; i++)
              val[i] /= ref_v;

          }
          break;

        case 5:
          {
            strncpy(var_name, "Ym/Ym_center", 64);

            /* get values and find point on profile closest to {0, 0, z} */

            cs_real_t min_d = HUGE_VAL;
            cs_real_t ref_v = 0;

            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = stat_vol_frac->val[c_id];
              cs_real_t d = sqrt(  cell_cen[c_id][0]*cell_cen[c_id][0]
                                 + cell_cen[c_id][1]*cell_cen[c_id][1]);
              if (d < min_d) {
                ref_v = val[i];
                min_d = d;
              }
            }

            cs_parall_min_loc_vals(1, &min_d, &ref_v);

            bft_printf("\nStat %s0 equal to %g for plane %s\n",
                       var_name, ref_v, name);

            /* now normalize by value at point on profile closest to {0, 0, z} */

            for (cs_lnum_t i = 0; i < n_cells; i++)
              val[i] /= ref_v;

          }
          break;

        }

        cs_post_write_probe_values
          (mesh_id,
           CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
           var_name,                       /* var_name */
           1,                              /* var_dim */
           CS_POST_TYPE_cs_real_t,
           0,                              /* parent location id */
           NULL,                           /* default interpolation */
           NULL,                           /* interpolation input */
           val,
           ts_post);

      }

      BFT_FREE(val);

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
