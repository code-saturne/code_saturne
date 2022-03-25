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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a profile based on centers of faces defined by a given
 *        criterion
 *
 * Here, the input points to string describing a selection criterion.
 *
 * \param[in]   input   pointer to selection criterion
 * \param[out]  n_elts  number of selected coordinates
 * \param[out]  coords  coordinates of selected elements.
 * \param[out]  s       curvilinear coordinates of selected elements
 */
/*----------------------------------------------------------------------------*/

/*! [post_profile_advanced_df_2] */
static void
_b_face_criterion_probes_define(void          *input,
                                cs_lnum_t     *n_elts,
                                cs_real_3_t  **coords,
                                cs_real_t    **s)
{
  const char *criterion = (const char *)input;

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_lnum_t   n_faces;
  cs_lnum_t  *face_ids;

  BFT_MALLOC(face_ids, m->n_b_faces, cs_lnum_t);
  cs_selector_get_b_face_list(criterion, &n_faces, face_ids);

  cs_real_3_t *_coords;
  cs_real_t *_s;
  BFT_MALLOC(_coords, n_faces, cs_real_3_t);
  BFT_MALLOC(_s, n_faces, cs_real_t);

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    for (cs_lnum_t j = 0; j < 3; j++)
      _coords[i][j] = mq->b_face_cog[face_ids[i]*3 + j];
    _s[i] = _coords[i][0];
  }

  BFT_FREE(face_ids);

  /* Set return values */

  *n_elts = n_faces;
  *coords = _coords;
  *s = _s;
}
/*! [post_profile_advanced_df_2] */

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
  /* redefine profile output options */

/*! [post_writer_profile] */
  cs_post_define_writer(CS_POST_WRITER_PROFILES,  /* writer_id */
                        "",                       /* writer name */
                        "profiles",
                        "plot",                   /* format name */
                        "dat",                    /* format options */
                        FVM_WRITER_FIXED_MESH,
                        false,                    /* output_at_start */
                        true,                     /* output at end */
                        -1,                       /* time step frequency */
                        -1.0);                    /* time value frequency */
/*! [post_writer_profile] */
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

/*! [post_profile_advanced_1] */
  const char *line_defs[]
    = {"-5.87", "buicesat-6",
       "2.59", "buicesat03",
       "5.98", "buicesat06",
       "12.75", "buicesat13",
       "13.56", "buicesat14",
       "16.14", "buicesat16",
       "16.93", "buicesat17",
       "19.53", "buicesat19",
       "20.32", "buicesat20",
       "22.91", "buicesat23",
       "23.71", "buicesat24",
       "26.3", "buicesat26",
       "27.09", "buicesat27",
       "29.69", "buicesat29",
       "30.48", "buicesat30",
       "33.07", "buicesat33",
       "33.87", "buicesat34",
       "39.85", "buicesat40",
       "46.62", "buicesat47",
       "53.39", "buicesat53",
       "60.17", "buicesat60",
       "66.94", "buicesat67",
       "73.71", "buicesat74"};

  for (int i = 0; i < 23; i++) {

    /* Define profiles */

    const cs_real_t x = atof(line_defs[i*2]);

    const cs_real_t start_coords[3] = {x, 0.1, 0.05};
    const cs_real_t end_coords[3]   = {x, -5, 0.05};

    cs_probe_set_t *pset
      = cs_probe_set_create_from_segment(line_defs[i*2+1],
                                         -1,  /* Automatic sampling */
                                         start_coords,
                                         end_coords);

    /* Associate writers */

    const int writer_ids[] = {CS_POST_WRITER_PROFILES};
    cs_probe_set_associate_writers(pset, 1, writer_ids);

  }
/*! [post_profile_advanced_1] */

  /* Define top and bottom profiles */

/*! [post_profile_advanced_2] */
  static const char *wall_defs[]
    = {"UP",   "buicstr",
       "DOWN", "buicinc"};

  for (int i = 0; i < 2; i++) {

    cs_probe_set_t *pset
      = cs_probe_set_create_from_local(wall_defs[i*2+1],
                                       _b_face_criterion_probes_define,
                                       (void *)wall_defs[i*2]);  /* input */

    cs_probe_set_option(pset, "boundary", "true");

    /* Associate writers */

    const int writer_ids[] = {CS_POST_WRITER_PROFILES};
    cs_probe_set_associate_writers(pset, 1, writer_ids);

  }
/*! [post_profile_advanced_2] */
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
  CS_UNUSED(mesh_name);
  CS_UNUSED(cat_id);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_vertices);
  CS_UNUSED(i_face_list);
  CS_UNUSED(vertex_list);

  if (probes != NULL) {

/*! [post_profile_advanced_var_0] */
    const cs_mesh_t *m = cs_glob_mesh;
    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

    const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

    const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;
    const cs_real_t uref = cs_glob_turb_ref_values->uref;
    const cs_real_t uref2 = uref*uref;

    const char *name = cs_probe_set_get_name(probes);

    bool is_profile = false;
    cs_probe_set_get_post_info(probes,
                               NULL,
                               NULL,
                               &is_profile,
                               NULL,
                               NULL,
                               NULL,
                               NULL,
                               NULL);

    const cs_time_step_t *ts_post = ts;
    if (is_profile && ts->nt_cur == ts->nt_max)
      ts = NULL;

    /* Common variables */

    const cs_real_t href = 1.;
/*! [post_profile_advanced_var_0] */

    /* For "x" profiles
       ---------------- */

/*! [post_profile_advanced_var_1] */
    if (strncmp(name, "buicesat", strlen("buicesat")) == 0) {

      cs_real_t *val;
      BFT_MALLOC(val, n_cells, cs_real_t);

      char var_name[64];

      /* mean x */

      cs_real_t x_sum[] = {0, 0};
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        cs_lnum_t c_id = cell_list[i];
        x_sum[0] += cell_cen[c_id][0];
      }
      x_sum[1] = n_cells;
      cs_parall_sum(2, CS_REAL_TYPE, x_sum);

      cs_real_t xpos = x_sum[0]/x_sum[1];

      /* Reynolds stresses */

      const cs_turb_model_t *turb_mdl = cs_glob_turb_model;
      const cs_turb_rans_model_t *turb_rans_mdl = cs_glob_turb_rans_model;

      cs_real_6_t *rij = NULL;
      BFT_MALLOC(rij, n_cells, cs_real_6_t);

      if (   turb_mdl->itytur == 2
          || turb_mdl->itytur == 5
          || turb_mdl->itytur == 6) {

        cs_field_interpolate_t interpolation_type = CS_FIELD_INTERPOLATE_MEAN;
        cs_post_evm_reynolds_stresses(interpolation_type,
                                      n_cells,
                                      cell_list,
                                      NULL, /* coords */
                                      rij);

      }
      else if (turb_mdl->itytur == 3) {

        cs_real_6_t *cvar_rij = (cs_real_6_t *)CS_F_(rij)->val;
        for (cs_lnum_t i = 0; i < n_cells; i++) {
          cs_lnum_t c_id = cell_list[i];
          for (cs_lnum_t j = 0; j < 6; j++)
            rij[i][j] = cvar_rij[c_id][j];
        }

      }

      /* Reynolds stresses invariants:
       * compute xsi and eta invariant of the Lumley triangle */

      cs_real_2_t *inv = NULL;
      BFT_MALLOC(inv, n_cells, cs_real_2_t);


      cs_post_anisotropy_invariant(n_cells,
                                   cell_list,
                                   NULL, /* coords */
                                   inv);

      /* Loop on columns */

      for (int col = 0; col < 9; col++) {

        switch(col) {

        case 0:
          {
            strncpy(var_name, "U*10+x/h", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = vel[c_id][0]*10 + xpos;
            }
          }
          break;

        case 1:
          {
            strncpy(var_name, "Y/H", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = mq->cell_cen[c_id*3 + 1] / href;
            }
          }
          break;

        case 2:
          {
            strncpy(var_name, "U/Uc", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = vel[c_id][0] / uref;
            }
          }
          break;

        case 3:
          {
            strncpy(var_name, "uu/Uc^2", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              val[i] = rij[i][0] / uref2;
            }
          }
          break;

        case 4:
          {
            strncpy(var_name, "uv/Uc^2", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              val[i] = rij[i][3] / uref2;
            }
          }
          break;

        case 5:
          {
            strncpy(var_name, "vv/Uc^2", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              val[i] = rij[i][1] / uref2;
            }
          }
          break;

        case 6:
          {
            strncpy(var_name, "X", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = cell_cen[c_id][0];
            }
          }
          break;

        case 7:
          {
            strncpy(var_name, "eta", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = inv[c_id][0];
            }
          }
          break;

        case 8:
          {
            strncpy(var_name, "xsi", 64);
            for (cs_lnum_t i = 0; i < n_cells; i++) {
              cs_lnum_t c_id = cell_list[i];
              val[i] = inv[c_id][1];
            }
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

      BFT_FREE(inv);
      BFT_FREE(rij);
      BFT_FREE(val);

    }
/*! [post_profile_advanced_var_1] */

    /* Pressure and friction coefficients
       ---------------------------------- */

/*! [post_profile_advanced_var_2] */
    else if (   strcmp(name, "buicstr") == 0
             || strcmp(name, "buicinc") == 0) {

      const cs_lnum_t *b_face_cells = m->b_face_cells;
      const cs_real_3_t *face_cog = (const cs_real_3_t *)mq->b_face_cog;
      const cs_real_t *distb = mq->b_dist;
      const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
      const cs_real_t *pres = CS_F_(p)->val;

      cs_real_t div_half_ro0_uref2 = 1. / (0.5 * phys_pro->ro0 * uref2);

      cs_real_t *val;
      BFT_MALLOC(val, n_b_faces, cs_real_t);

      char var_name[64];

      /* Reference pressure */

      cs_real_t xyz_ref[3] = {-1.7, -0.5, 0.005};

      cs_real_t pref = 0;
      cs_lnum_t pref_id;
      int pref_rank;

      cs_geom_closest_point(m->n_cells,
                            (const cs_real_3_t *)(mq->cell_cen),
                            xyz_ref,
                            &pref_id,
                            &pref_rank);

      if (pref_rank == cs_glob_rank_id)
        pref = pres[pref_id];
      cs_parall_bcast(pref_rank, 1, CS_REAL_TYPE, &pref);

      /* Stresses */

      cs_real_3_t *stresses;
      BFT_MALLOC(stresses, n_b_faces, cs_real_3_t);
      cs_post_stress_tangential(n_b_faces, b_face_list, stresses);

      /* Loop on columns */

      for (int col = 0; col < 5; col++) {

        switch(col) {

        case 0:
          {
            strncpy(var_name, "X/H", 64);
            for (cs_lnum_t i = 0; i < n_b_faces; i++) {
              cs_lnum_t f_id = b_face_list[i];
              val[i] = face_cog[f_id][0] / href;
            }
          }
          break;

        case 1:
          {
            strncpy(var_name, "CP", 64);
            for (cs_lnum_t i = 0; i < n_b_faces; i++) {
              cs_lnum_t f_id = b_face_list[i];
              cs_lnum_t c_id = b_face_cells[f_id];
              val[i] = (pres[c_id] - pref) * div_half_ro0_uref2;
            }
          }
          break;

        case 2:
          {
            strncpy(var_name, "CF", 64);
            for (cs_lnum_t i = 0; i < n_b_faces; i++) {
              val[i] = cs_math_3_norm(stresses[i]) * div_half_ro0_uref2;
            }
          }
          break;

        case 3:
          {
            strncpy(var_name, "U/UREF", 64);
            for (cs_lnum_t i = 0; i < n_b_faces; i++) {
              /* previous value in val[i] from case2:
                 norm(stresses[i])/(0.5.ro0.uref^2) */
              val[i] = copysign(val[i], stresses[i][0]);
            }
          }
          break;

        case 4:
          {
            strncpy(var_name, "YPLUS", 64);
            for (cs_lnum_t i = 0; i < n_b_faces; i++) {
              cs_lnum_t f_id = b_face_list[i];
              cs_lnum_t c_id = b_face_cells[f_id];
              val[i] = sqrt(fabs(vel[c_id][0])*distb[f_id]*phys_pro->viscl0);
            }
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

      BFT_FREE(stresses);
      BFT_FREE(val);
    }
/*! [post_profile_advanced_var_2] */

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
