/*============================================================================
 * Handle successive mesh preprocessing operations.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_boundary_zone.h"
#include "cs_ext_neighborhood.h"
#include "cs_gui.h"
#include "cs_gui_mesh.h"
#include "cs_internal_coupling.h"
#include "cs_join.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_cartesian.h"
#include "cs_mesh_from_builder.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_renumber.h"
#include "cs_mesh_save.h"
#include "cs_mesh_to_builder.h"
#include "cs_mesh_warping.h"
#include "cs_parall.h"
#include "cs_partition.h"
#include "cs_porous_model.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_preprocessor_data.h"
#include "cs_timer_stats.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_preprocess.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_preprocess.c
        Handle successive preprocessing operations.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for Fortran functions used only through this program unit.
 *============================================================================*/

/* majgeo */

extern void cs_f_majgeo(const cs_lnum_t    *ncel,
                        const cs_lnum_t    *ncelet,
                        const cs_lnum_t    *nfac,
                        const cs_lnum_t    *nfabor,
                        const cs_lnum_t    *nsom,
                        const cs_lnum_t    *lndfac,
                        const cs_lnum_t    *lndfbr,
                        const int64_t      *ncelgb,
                        const int64_t      *nfacgb,
                        const int64_t      *nfbrgb,
                        const int64_t      *nsomgb,
                        const cs_lnum_2_t   ifacel[],
                        const cs_lnum_t     ifabor[],
                        const int           ifmfbr[],
                        const int           ifmcel[],
                        const cs_lnum_t     ipnfac[],
                        const cs_lnum_t     nodfac[],
                        const cs_lnum_t     ipnfbr[],
                        const cs_lnum_t     nodfbr[],
                        const int           isympa[],
                        const int           isolid_0[],
                        const cs_real_t    *volmin,
                        const cs_real_t    *volmax,
                        const cs_real_t    *voltot,
                        const cs_real_t     xyzcen[],
                        const cs_real_t     surfac[],
                        const cs_real_t     surfbo[],
                        const cs_real_t     suffac[],
                        const cs_real_t     suffbo[],
                        const cs_real_t     cdgfac[],
                        const cs_real_t     cdgfbo[],
                        const cs_real_t     xyznod[],
                        const cs_real_t     volume[],
                        const cs_real_t     cell_f_vol[],
                        const cs_real_t     surfan[],
                        const cs_real_t     surfbn[],
                        const cs_real_t     suffan[],
                        const cs_real_t     suffbn[],
                        const cs_real_t     dist[],
                        const cs_real_t     distb[],
                        const cs_real_t     pond[],
                        const cs_real_t     dijpf[],
                        const cs_real_t     diipb[],
                        const cs_real_t     dofij[]);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine if preprocessing is needed.
 *
 * Preprocessing is ignored when a ./restart/mesh_input file is present but
 * no ./mesh_input file or directory is present. In this case, restart mesh
 * file is read, and all other preprocessing steps are skipped.
 *
 * \return  true if preprocessing is needed, false if only reading is needed.
 */
/*----------------------------------------------------------------------------*/

bool
cs_preprocess_mesh_is_needed(void)
{
  bool retval = true;

  int is_restart = (cs_preprocessor_data_is_restart()) ? 1 : 0;

#if defined(HAVE_MPI)
  if (cs_glob_rank_id >= 0)
    MPI_Bcast(&is_restart,  1, MPI_INT,  0, cs_glob_mpi_comm);
#endif

  if (is_restart)
    retval = false;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define all mesh preprocessing operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_define(void)
{
  if (cs_preprocess_mesh_is_needed() == false)
    return;

  /* Define meshes to read */

  cs_user_mesh_input();

  /* Define joining and periodicity parameters if requested
     Must be done before initi1 for the sake of verification */

  cs_gui_mesh_define_joinings();
  cs_user_join();

  cs_gui_mesh_define_periodicities();
  cs_user_periodicity();

  cs_gui_mesh_warping();
  cs_user_mesh_warping();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply all mesh preprocessing operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh(cs_halo_type_t   halo_type)
{
  double  t1, t2;

  int t_stat_id = cs_timer_stats_id_by_name("mesh_processing");

  int t_top_id = cs_timer_stats_switch(t_stat_id);

  bool allow_modify = cs_preprocess_mesh_is_needed();

  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Disable all writers until explicitely enabled for this stage */

  cs_post_disable_writer(0);

  /* Set partitioning options */

  {
    int j_id;
    bool join = false;
    bool join_periodic = false;

    for (j_id = 0; j_id < cs_glob_n_joinings; j_id++) {
      if ((cs_glob_join_array[j_id])->param.perio_type == FVM_PERIODICITY_NULL)
        join = true;
      else
        join_periodic = true;
    }

    cs_partition_set_preprocess_hints(join, join_periodic);
    cs_gui_partition();
    cs_user_partition();
  }

  /* Read Preprocessor output */

  cs_preprocessor_data_read_mesh(m,
                                 cs_glob_mesh_builder);

  if (allow_modify) {

    /* Join meshes / build periodicity links if necessary */

    cs_join_all(true);

    /* Insert boundaries if necessary */

    cs_gui_mesh_boundary(m);
    cs_user_mesh_boundary(m);

    cs_internal_coupling_preprocess(m);

  }

  /* Initialize extended connectivity, ghost cells and other remaining
     parallelism-related structures */

  cs_mesh_init_halo(m, cs_glob_mesh_builder, halo_type, m->verbosity, true);
  cs_mesh_update_auxiliary(m);

  if (allow_modify) {

    /* Possible geometry modification */

    cs_gui_mesh_extrude(m);
    cs_user_mesh_modify(m);

    /* Discard isolated faces if present */

    cs_post_add_free_faces();
    cs_mesh_discard_free_faces(m);

    /* Smoothe mesh if required */

    cs_gui_mesh_smoothe(m);
    cs_user_mesh_smoothe(m);

    /* Triangulate warped faces if necessary */

    {
      double  cwf_threshold = -1.0;
      int  cwf_post = 0;

      cs_mesh_warping_get_defaults(&cwf_threshold, &cwf_post);

      if (cwf_threshold >= 0.0) {

        t1 = cs_timer_wtime();
        cs_mesh_warping_cut_faces(m, cwf_threshold, cwf_post);
        t2 = cs_timer_wtime();

        bft_printf(_("\n Cutting warped boundary faces (%.3g s)\n"), t2-t1);

      }
    }

    /* Now that mesh modification is finished, save mesh if modified */

    cs_gui_mesh_save_if_modified(m);
    cs_user_mesh_save(m); /* Disable or force */
  }

  bool need_partition = cs_partition_get_preprocess();
  if (m->modified & CS_MESH_MODIFIED_BALANCE)
    need_partition = true;

  bool need_save = false;
  if (   (m->modified > 0 && m->save_if_modified > 0)
      || m->save_if_modified > 1)
    need_save = true;

  if (need_partition) {
    if (need_save) {
      cs_mesh_save(m, cs_glob_mesh_builder, NULL, "mesh_output.csm");
      need_save = false;
    }
    else
      cs_mesh_to_builder(m, cs_glob_mesh_builder, true, NULL);

    cs_partition(m, cs_glob_mesh_builder, CS_PARTITION_MAIN);
    cs_mesh_from_builder(m, cs_glob_mesh_builder);
    cs_mesh_init_halo(m, cs_glob_mesh_builder, halo_type, m->verbosity, true);
    cs_mesh_update_auxiliary(m);
  }

  else if (need_save)
    cs_mesh_save(m, NULL, NULL, "mesh_output.csm");

  m->n_b_faces_all = m->n_b_faces;
  m->n_g_b_faces_all = m->n_g_b_faces;

  /* Destroy the temporary structure used to build the main mesh */

  cs_mesh_builder_destroy(&cs_glob_mesh_builder);

  /* Destroy cartesian mesh builder if necessary */
  cs_mesh_cartesian_params_destroy();

  /* Renumber mesh based on code options */

  cs_user_numbering();

  cs_renumber_mesh(m);

  /* Initialize group classes */

  cs_mesh_init_group_classes(m);

  /* Print info on mesh */

  cs_mesh_print_info(m, _("Mesh"));

  /* Second pass to define internal coupling locators */
  cs_internal_coupling_map(m);

  /* Compute geometric quantities related to the mesh */

  bft_printf_flush();

  t1 = cs_timer_wtime();

  /* If fluid_solid mode is activated: disable solid cells for the dynamics */
  cs_velocity_pressure_model_t *vp_model = cs_get_glob_velocity_pressure_model();
  if (vp_model->fluid_solid)
    mq->has_disable_flag = 1;

  cs_mesh_quantities_compute(m, mq);

  cs_mesh_bad_cells_detect(m, mq);
  cs_user_mesh_bad_cells_tag(m, mq);
  t2 = cs_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_init_selectors();
  cs_mesh_location_build(m, -1);
  cs_volume_zone_build_all(true);
  cs_volume_zone_print_info();
  cs_boundary_zone_build_all(true);
  cs_boundary_zone_print_info();

  cs_ext_neighborhood_reduce(m, mq);

  /* If fluid_solid mode is activateed, disable solid cells for the dynamics */
  cs_porous_model_init_disable_flag();
  if (vp_model->fluid_solid) {
    assert(mq->has_disable_flag == 1);
    cs_volume_zone_tag_cell_type(CS_VOLUME_ZONE_SOLID, 1, mq->c_disable_flag);
  }

  /* For debugging purposes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_mesh_dump(m);
  cs_mesh_quantities_dump(m, mq);
#endif

  /* Re-enable writers disabled when entering this stage */

  cs_post_enable_writer(0);

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply numbering changes to ignore selected boundary faces.
 *
 * \param[in, out]  m         pointer to mesh structure
 * \param[in, out]  mq        pointer to mesh quantities structure
 * \param[in]       n_faces   number of selected faces
 * \param[in]       face_ids  ids of selected faces
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_selected_b_faces_ignore(cs_mesh_t             *m,
                                           cs_mesh_quantities_t  *mq,
                                           cs_lnum_t              n_faces,
                                           const cs_lnum_t        face_ids[])
{
  double  t1, t2;

  int t_stat_id = cs_timer_stats_id_by_name("mesh_processing");

  int t_top_id = cs_timer_stats_switch(t_stat_id);

  cs_mesh_update_auxiliary(m);

  /* Renumber mesh based on code options */

  cs_renumber_b_faces_select_ignore(m,
                                    n_faces,
                                    face_ids);

  /* Rebuild some structures */

  cs_mesh_update_b_cells(m);
  cs_mesh_init_group_classes(m);

  /* Print info on mesh */

  cs_mesh_print_info(m, _("Mesh"));

  /* Second pass to define internal coupling locators */

  cs_internal_coupling_map(m);

  /* Compute geometric quantities related to the mesh */

  bft_printf_flush();

  t1 = cs_timer_wtime();

  cs_mesh_quantities_compute(m, mq);

  t2 = cs_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_update_selectors(cs_glob_mesh);
  cs_mesh_location_build(cs_glob_mesh, -1);
  cs_volume_zone_build_all(true);
  cs_volume_zone_print_info();
  cs_boundary_zone_build_all(true);
  cs_boundary_zone_print_info();

  /* For debugging purposes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_mesh_dump(cs_glob_mesh);
  cs_mesh_quantities_dump(cs_glob_mesh, cs_glob_mesh_quantities);
#endif

  cs_timer_stats_switch(t_top_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update fortran arrays relative to the global mesh.
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_update_fortran(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  int64_t n_g_cells = m->n_g_cells;
  int64_t n_g_i_faces = m->n_g_i_faces;
  int64_t n_g_b_faces = m->n_g_b_faces;
  int64_t n_g_vertices = m->n_g_vertices;

  cs_f_majgeo(&(m->n_cells),
              &(m->n_cells_with_ghosts),
              &(m->n_i_faces),
              &(m->n_b_faces),
              &(m->n_vertices),
              &(m->i_face_vtx_connect_size),
              &(m->b_face_vtx_connect_size),
              &n_g_cells,
              &n_g_i_faces,
              &n_g_b_faces,
              &n_g_vertices,
              (const cs_lnum_2_t *)(m->i_face_cells),
              m->b_face_cells,
              m->b_face_family,
              m->cell_family,
              m->i_face_vtx_idx,
              m->i_face_vtx_lst,
              m->b_face_vtx_idx,
              m->b_face_vtx_lst,
              mq->b_sym_flag,
              mq->c_disable_flag,
              &(mq->min_vol),
              &(mq->max_vol),
              &(mq->tot_vol),
              mq->cell_cen,
              mq->i_face_normal,
              mq->b_face_normal,
              mq->i_f_face_normal,
              mq->b_f_face_normal,
              mq->i_face_cog,
              mq->b_face_cog,
              m->vtx_coord,
              mq->cell_vol,
              mq->cell_f_vol,
              mq->i_face_surf,
              mq->b_face_surf,
              mq->i_f_face_surf,
              mq->b_f_face_surf,
              mq->i_dist,
              mq->b_dist,
              mq->weight,
              mq->dijpf,
              mq->diipb,
              mq->dofij);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map some mesh arrays for use on device.
 *
 * More elements may be mapped dependin on which arrays are used in
 * accelerated algorithms.
 *
 * \param[in]  alloc_mode  chosen allocation mode
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_update_device(cs_alloc_mode_t  alloc_mode)
{
#if 0  /* For local testing of various modes */
  alloc_mode = CS_ALLOC_HOST_DEVICE_PINNED;
#endif

  cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_b_cells = m->n_b_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Mesh structures
     --------------- */

  {
    CS_REALLOC_HD(m->i_face_cells, n_i_faces, cs_lnum_2_t, alloc_mode);
    CS_REALLOC_HD(m->b_face_cells, n_b_faces, cs_lnum_t, alloc_mode);

    if (m->b_cells != NULL)
      CS_REALLOC_HD(m->b_cells, n_b_cells, cs_lnum_t, alloc_mode);
  }

  if (m->cell_cells_idx != NULL) {
    CS_REALLOC_HD(m->cell_cells_idx, n_cells+1, cs_lnum_t, alloc_mode);
    CS_REALLOC_HD(m->cell_cells_lst, m->cell_cells_idx[n_cells], cs_lnum_t,
                  alloc_mode);
  }

  /* Additional adjacencies
     ---------------------- */

  cs_mesh_adjacencies_update_device(alloc_mode);

  /* Update Fortran mappings as some addresses may have changed */

  cs_preprocess_mesh_update_fortran();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
