/*============================================================================
 * Define postprocessing output.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "bft_printf.h"
#include "bft_mem.h"
#include "bft_error.h"

#include "fvm_periodicity.h"

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

#include "cs_lagr_adh.h"
#include "cs_lagr_car.h"
#include "cs_lagr_clogging.h"
#include "cs_lagr_coupling.h"
#include "cs_lagr_deposition_model.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_extract.h"
#include "cs_lagr_geom.h"
#include "cs_lagr_gradients.h"
#include "cs_lagr.h"
#include "cs_lagr_head_losses.h"
#include "cs_lagr_injection.h"
#include "cs_lagr_lec.h"
#include "cs_lagr_log.h"
#include "cs_lagr_new.h"
#include "cs_lagr_options.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_poisson.h"
#include "cs_lagr_post.h"
#include "cs_lagr_precipitation_model.h"
#include "cs_lagr_print.h"
#include "cs_lagr_prototypes.h"
#include "cs_lagr_query.h"
#include "cs_lagr_restart.h"
#include "cs_lagr_resuspension.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_sde.h"
#include "cs_lagr_sde_model.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_tracking.h"




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
  /* Output de la vitesse moyenne en Z des particules dans chaque cellule
     ------------------------------------------------------ */


  if (cat_id == CS_POST_MESH_VOLUME) {
  
    cs_real_t *s_cell_Vx;
    BFT_MALLOC(s_cell_Vx, n_cells, cs_real_t);

    cs_real_t *s_cell_Vy;
    BFT_MALLOC(s_cell_Vy, n_cells, cs_real_t);

    cs_real_t *s_cell_Vz;
    BFT_MALLOC(s_cell_Vz, n_cells, cs_real_t);

    cs_real_t *s_cell_var_V_xx;
    BFT_MALLOC(s_cell_var_V_xx, n_cells, cs_real_t);

    cs_real_t *s_cell_var_V_yy;
    BFT_MALLOC(s_cell_var_V_yy, n_cells, cs_real_t);

    cs_real_t *s_cell_var_V_zz;
    BFT_MALLOC(s_cell_var_V_zz, n_cells, cs_real_t);
    
    cs_real_t *s_cell_var_V_xy;
    BFT_MALLOC(s_cell_var_V_xy, n_cells, cs_real_t);

    cs_real_t *s_cell_var_V_yz;
    BFT_MALLOC(s_cell_var_V_yz, n_cells, cs_real_t);

    cs_real_t *s_cell_var_V_xz;
    BFT_MALLOC(s_cell_var_V_xz, n_cells, cs_real_t);

    cs_real_t *s_cell_partvolfrac;
    BFT_MALLOC(s_cell_partvolfrac, n_cells, cs_real_t);

    cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
    const cs_lagr_attribute_map_t *p_am = p_set->p_am;
    
    for (cs_lnum_t i = 0; i < n_cells; i++) {

        s_cell_Vx[i] = 0;
        s_cell_Vy[i] = 0;
        s_cell_Vz[i] = 0;
        
        s_cell_var_V_xx[i] = 0;
        s_cell_var_V_yy[i] = 0;
        s_cell_var_V_zz[i] = 0;
        s_cell_var_V_xy[i] = 0;
        s_cell_var_V_yz[i] = 0;
        s_cell_var_V_xz[i] = 0;

        s_cell_partvolfrac[i] = 0;

        cs_lnum_t cell_id = cell_list[i];

        cs_lnum_t j=0;
        
        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

            unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

            cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

            cs_real_t *part_vel_seen
              = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

            cs_real_t diam = cs_lagr_particle_get_real(particle, p_am,
              CS_LAGR_DIAMETER);

            cs_real_t p_weight = cs_lagr_particle_get_real(particle, p_am,
              CS_LAGR_STAT_WEIGHT);

            cs_real_t vol = cs_glob_mesh_quantities->cell_vol[cell_id];

            
            if (iel == cell_id) {
                j=j+1;

                s_cell_Vx[i] = s_cell_Vx[i]+part_vel_seen[0];
                s_cell_Vy[i] = s_cell_Vy[i]+part_vel_seen[1];
                s_cell_Vz[i] = s_cell_Vz[i]+part_vel_seen[2];

                s_cell_partvolfrac[i] = s_cell_partvolfrac[i]+diam*diam*diam *
                  cs_math_pi / 6.0 * p_weight / vol;

            }
            
        }
        
        if (j > 0) {
          s_cell_Vx[i] = s_cell_Vx[i]/j;
          s_cell_Vy[i] = s_cell_Vy[i]/j;
          s_cell_Vz[i] = s_cell_Vz[i]/j;
        }

        cs_lnum_t k=0;
        
        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {
            unsigned char *particle = p_set->p_buffer + p_am->extents * npt;
            cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

            cs_real_t *part_vel_seen
              = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

            
            if (iel == cell_id) {
                k=k+1;
                s_cell_var_V_xx[i] = s_cell_var_V_xx[i]+(part_vel_seen[0]-s_cell_Vx[i])*(part_vel_seen[0]-s_cell_Vx[i]);
                s_cell_var_V_yy[i] = s_cell_var_V_yy[i]+(part_vel_seen[1]-s_cell_Vy[i])*(part_vel_seen[1]-s_cell_Vy[i]);
                s_cell_var_V_zz[i] = s_cell_var_V_zz[i]+(part_vel_seen[2]-s_cell_Vz[i])*(part_vel_seen[2]-s_cell_Vz[i]);
                s_cell_var_V_xy[i] = s_cell_var_V_xy[i]+(part_vel_seen[0]-s_cell_Vx[i])*(part_vel_seen[1]-s_cell_Vy[i]);
                s_cell_var_V_yz[i] = s_cell_var_V_yz[i]+(part_vel_seen[1]-s_cell_Vy[i])*(part_vel_seen[2]-s_cell_Vz[i]);
                s_cell_var_V_xz[i] = s_cell_var_V_xz[i]+(part_vel_seen[0]-s_cell_Vx[i])*(part_vel_seen[2]-s_cell_Vz[i]);
            }

        }
        
        if (k > 0) {
          s_cell_var_V_xx[i] = s_cell_var_V_xx[i]/k;
          s_cell_var_V_yy[i] = s_cell_var_V_yy[i]/k;
          s_cell_var_V_zz[i] = s_cell_var_V_zz[i]/k;
          s_cell_var_V_xy[i] = s_cell_var_V_xy[i]/k;
          s_cell_var_V_yz[i] = s_cell_var_V_yz[i]/k;
          s_cell_var_V_xz[i] = s_cell_var_V_xz[i]/k;
        }
    }
    

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "moy_Vs_x",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_Vx,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "moy_Vs_y",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_Vy,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "moy_Vs_z",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_Vz,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);



    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_xx",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_xx,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_yy",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_yy,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_zz",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_zz,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_xy",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_xy,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_yz",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_yz,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "var_Vs_xz",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_var_V_xz,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    cs_post_write_var(mesh_id,
                      CS_POST_WRITER_ALL_ASSOCIATED,  /* writer id filter */
                      "partvolfraction",                  /* var_name */
                      1,                              /* var_dim */
                      true,                           /* interlace, */
                      false,                          /* use_parent */
                      CS_POST_TYPE_cs_real_t,         /* var_type */
                      s_cell_partvolfrac,                         /* cel_vals */
                      NULL,                           /* i_face_vals */
                      NULL,                           /* b_face_vals */
                      ts);

    BFT_FREE(s_cell_Vx);
    BFT_FREE(s_cell_Vy);
    BFT_FREE(s_cell_Vz);
    
    BFT_FREE(s_cell_var_V_xx);
    BFT_FREE(s_cell_var_V_yy);
    BFT_FREE(s_cell_var_V_zz);
    BFT_FREE(s_cell_var_V_xy);
    BFT_FREE(s_cell_var_V_yz);
    BFT_FREE(s_cell_var_V_xz);
    BFT_FREE(s_cell_partvolfrac);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * Override default frequency or calculation end based output.
 *
 * This allows fine-grained control of activation or deactivation,
 *
 * \param  nt_max_abs  maximum time step number
 * \param  nt_cur_abs  current time step number
 * \param  t_cur_abs   absolute time at the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_user_postprocess_activate(int     nt_max_abs,
                             int     nt_cur_abs,
                             double  t_cur_abs)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
