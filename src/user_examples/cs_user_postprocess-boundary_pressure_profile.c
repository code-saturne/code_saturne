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
  /* redefine profile output options such as format, frequency, etc... */

  cs_post_define_writer(CS_POST_WRITER_PROFILES,  /* writer_id */
                        "",                       /* writer name */
                        "profiles",
                        "plot",                   /* format name */
                        "dat",                    /* format options */
                        FVM_WRITER_FIXED_MESH,
                        false,                    /* output_at_start */
                        true,                     /* output at end */
                        1,                       /* time step frequency */
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
  /* Alias global pointers to mesh and mesh quantities structures */
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* A probe will be located at each boundary face center */

  /* Create probes set of name foil_profile from criterion FOIL_WALL */
  cs_probe_set_t *pset
    = cs_probe_set_create_from_local("foil_profile", /* probes set name */
                                     cs_b_face_criterion_probes_define,
                                                      /* probe def. function */
                                     (void *)"FOIL_WALL"); /* input */

  /* Indicate that the probes are located on the boundary */
  cs_probe_set_option(pset, "boundary", "true");

  /* Associate profile writer to this probes set */
  const int writer_ids[] = {CS_POST_WRITER_PROFILES};
  cs_probe_set_associate_writers(pset, 1, writer_ids);
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
  /* function possibly called for each postprocessing mesh, hence also for
     each probes set */

  /* check if current mesh is a probes set */
  if (probes != NULL) {

    const char *name = cs_probe_set_get_name(probes);

    /* check that current probes set is foil_profile */
    if (strncmp(name, "foil_profile", strlen("foil_profile")) == 0) {

      const cs_mesh_t *m = cs_glob_mesh;
      const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
      const cs_real_3_t *b_face_cog = (const cs_real_3_t *)mq->b_face_cog;

      cs_real_t *val;
      BFT_MALLOC(val, n_b_faces, cs_real_t);

      /* x coordinate */
      for (cs_lnum_t i = 0; i < n_b_faces; i++) {
        cs_lnum_t face_id = b_face_list[i];
        val[i] = b_face_cog[face_id][0];
      }

      /* post-process x coordinate */
      cs_post_write_probe_values
        (mesh_id,
         CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
         "X",                           /* var_name */
         1,                              /* var_dim */
         CS_POST_TYPE_cs_real_t,
         0,                              /* parent location id */
         NULL,                           /* default interpolation */
         NULL,                           /* interpolation input */
         val,
         ts);

      /* compute pressure on selected boundary faces */
      cs_post_b_pressure(n_b_faces, b_face_list, val);

      const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
      cs_real_t p0 = phys_pro->p0; /* reference pressure */
      cs_real_t ro0 = phys_pro->ro0; /* reference density */

      const cs_real_t uref = cs_glob_turb_ref_values->uref; /*ref. velocity */
      const cs_real_t uref2 = uref*uref;

      /* reference values can be set in GUI */

      /* 1/(1/2 rho U^2) */
      cs_real_t div_half_ro0_uref2 = 1. / (0.5 * ro0 * uref2);

      /* compute CP at each selected boundary face */
      for (cs_lnum_t i = 0; i < n_b_faces; i++)
        val[i] = (val[i] - p0) * div_half_ro0_uref2;

      /* post-process CP */
      cs_post_write_probe_values
        (mesh_id,
         CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
         "CP",                           /* var_name */
         1,                              /* var_dim */
         CS_POST_TYPE_cs_real_t,
         0,                              /* parent location id */
         NULL,                           /* default interpolation */
         NULL,                           /* interpolation input */
         val,
         ts);

      BFT_FREE(val);
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
