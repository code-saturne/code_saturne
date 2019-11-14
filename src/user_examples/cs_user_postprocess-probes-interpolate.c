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
#include "cs_field_operator.h"
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

/*! [post_interpolate_p1_func] */

/*----------------------------------------------------------------------------
 * Interpolate values defined on a mesh location at a given set of
 * points using a P0 interpolation.
 *
 * This function assumes the input is a field name. If no field matches
 * the name, a "generic" interpolation (assuming homogeneous Neumann boundary
 * conditions) is used).
 *
 * \param[in, out]  input           pointer to optional (untyped) value
 *                                  or structure.
 * \param[in]       datatype        associated datatype
 * \param[in]       val_dim         dimension of data values
 * \param[in]       n_points        number of interpolation points
 * \param[in]       point_location  location of points in mesh elements
 * \param[in]       point_coords    point coordinates
 * \param[in]       location_vals   values at mesh location
 * \param[out]      point_vals      interpolated values at points
 *----------------------------------------------------------------------------*/

static void
_cs_interpolate_from_location_p1(void                *input,
                                 cs_datatype_t        datatype,
                                 int                  val_dim,
                                 cs_lnum_t            n_points,
                                 const cs_lnum_t      point_location[],
                                 const cs_real_3_t    point_coords[],
                                 const void          *location_vals,
                                 void                *point_vals)
{
  /* If used with a non-real argument type, use P0 interpolation */

  if (   datatype != CS_REAL_TYPE
      || (val_dim != 1 && val_dim != 3 && val_dim != 6)) {
    cs_interpolate_from_location_p0(input,
                                    datatype,
                                    val_dim,
                                    n_points,
                                    point_location,
                                    point_coords,
                                    location_vals,
                                    point_vals);
    return;
  }

  /* Main usage */

  cs_mesh_t *m = cs_glob_mesh;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

  cs_gradient_type_t gradient_type = CS_GRADIENT_ITER;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_gradient_type_by_imrgra(cs_glob_space_disc->imrgra,
                             &gradient_type,
                             &halo_type);

  cs_field_t *f = NULL;
  if (input != NULL) {
    const char *name = input;
    f = cs_field_by_name_try(name);
  }

  switch(val_dim) {
  case 1:
    {
      const cs_real_t *c_vals = (const cs_real_t *)location_vals;
      cs_real_t *p_vals = (cs_real_t *)point_vals;

      cs_real_3_t *grad;
      BFT_MALLOC(grad, m->n_cells_with_ghosts, cs_real_3_t);
      if (f != NULL)
        cs_field_gradient_scalar(f,
                                 false, /* use_previous_t */
                                 1, /* inc */
                                 true, /* _recompute_cocg */
                                 grad);
      else
        cs_gradient_scalar_synced_input("user",
                                        gradient_type,
                                        halo_type,
                                        1,      /* inc */
                                        true,   /* recompute_cocg */
                                        100,    /* n_r_sweeps */
                                        0,      /* tr_dim */
                                        0,      /* hyd_p_flag */
                                        0,      /* w_stride */
                                        0,      /* verbosity */
                                        -1,     /* clip_mode */
                                        1e-5,   /* epsilon */
                                        0,      /* extrap */
                                        1.5,    /* clip_coeff */
                                        NULL,   /* f_ext[] */
                                        NULL,   /* bc_coeff_a */
                                        NULL,   /* bc_coeff_b */
                                        c_vals,
                                        NULL,   /* c_weight */
                                        NULL,   /* cpl */
                                        grad);

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          p_vals[i] = c_vals[c_id] + grad[c_id][0]*d[0]
                                   + grad[c_id][1]*d[1]
                                   + grad[c_id][2]*d[2];
        }
        else
          p_vals[i] = 0;
      }

      BFT_FREE(grad);
    }
    break;

  case 3:
    {
      const cs_real_3_t *c_vals = (const cs_real_3_t *)location_vals;
      cs_real_3_t *p_vals = (cs_real_3_t *)point_vals;

      cs_real_33_t *grad;
      BFT_MALLOC(grad, m->n_cells_with_ghosts, cs_real_33_t);
      if (f != NULL)
        cs_field_gradient_vector(f,
                                 false, /* use_previous_t */
                                 1, /* inc */
                                 grad);
      else
        cs_gradient_vector_synced_input("user",
                                        gradient_type,
                                        halo_type,
                                        1,      /* inc */
                                        100,    /* n_r_sweeps */
                                        0,      /* verbosity */
                                        -1,     /* clip_mode */
                                        1e-5,   /* epsilon */
                                        1.5,    /* clip_coeff */
                                        NULL,   /* bc_coeff_a */
                                        NULL,   /* bc_coeff_b */
                                        c_vals,
                                        NULL,   /* c_weight */
                                        NULL,   /* cpl */
                                        grad);

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          for (cs_lnum_t j = 0; j < 3; j++) {
            p_vals[i][j] = c_vals[c_id][j] + grad[c_id][j][0]*d[0]
                                           + grad[c_id][j][1]*d[1]
                                           + grad[c_id][j][2]*d[2];
          }
        }
        else {
          for (cs_lnum_t j = 0; j < 6; j++)
            p_vals[i][j] = 0;
        }
      }

      BFT_FREE(grad);
    }
    break;

  case 6:
    {
      const cs_real_6_t *c_vals = (const cs_real_6_t *)location_vals;
      cs_real_6_t *p_vals = (cs_real_6_t *)point_vals;

      cs_real_63_t *grad;
      BFT_MALLOC(grad, m->n_cells_with_ghosts, cs_real_63_t);
      if (f != NULL)
        cs_field_gradient_tensor(f,
                                 false, /* use_previous_t */
                                 1, /* inc */
                                 grad);
      else
        cs_gradient_tensor_synced_input("user",
                                        gradient_type,
                                        halo_type,
                                        1,      /* inc */
                                        100,    /* n_r_sweeps */
                                        0,      /* verbosity */
                                        -1,     /* clip_mode */
                                        1e-5,   /* epsilon */
                                        1.5,    /* clip_coeff */
                                        NULL,   /* bc_coeff_a */
                                        NULL,   /* bc_coeff_b */
                                        c_vals,
                                        grad);

      for (cs_lnum_t i = 0; i < n_points; i++) {
        cs_lnum_t c_id = point_location[i];
        if (c_id > -1) {
          cs_real_t d[3] = {point_coords[i][0] - cell_cen[c_id][0],
                            point_coords[i][1] - cell_cen[c_id][1],
                            point_coords[i][2] - cell_cen[c_id][2]};

          for (cs_lnum_t j = 0; j < 6; j++) {
            p_vals[i][j] = c_vals[c_id][j] + grad[c_id][j][0]*d[0]
                                           + grad[c_id][j][1]*d[1]
                                           + grad[c_id][j][2]*d[2];
          }
        }
        else {
          for (cs_lnum_t j = 0; j < 6; j++)
            p_vals[i][j] = 0;
        }
      }

      BFT_FREE(grad);
    }
    break;

  default:
    assert(0);
  }
}

/*! [post_interpolate_p1_func] */

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
    const char *name = cs_probe_set_get_name(probes);

    int n_p_fields = 2;
    const char *p_field_names[] = {"velocity", "temperature"};

    for (int i = 0; i < n_p_fields; i++) {

      cs_field_t *f = cs_field_by_name_try(p_field_names[i]);

      if (f != NULL) {

        /* use different name to avoid conflict with field name in case already
           present in probe set through default output */

        char p_name[64];
        snprintf(p_name, 63, "%s_p", f->name); p_name[63] = '\0';

        cs_post_write_probe_values
          (mesh_id,
           CS_POST_WRITER_ALL_ASSOCIATED,    /* writer id filter */
           p_name,                           /* var_name */
           f->dim,                           /* var_dim */
           CS_POST_TYPE_cs_real_t,
           1,                                /* parent location id */
           _cs_interpolate_from_location_p1, /* P1 interpolation */
           f->name,                          /* interpolation input */
           f->val,
           ts);

      }

    }

  }
  /*! [post_probes_interpolate_var_1] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
