/*============================================================================
 * Handle successive mesh preprocessing operations.
 *============================================================================*/

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
#include "cs_gui.h"
#include "cs_gui_mesh.h"
#include "cs_internal_coupling.h"
#include "cs_join.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mesh.h"
#include "cs_mesh_from_builder.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_renumber.h"
#include "cs_mesh_save.h"
#include "cs_mesh_to_builder.h"
#include "cs_mesh_warping.h"
#include "cs_parall.h"
#include "cs_partition.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_preprocessor_data.h"
#include "cs_timer_stats.h"
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

extern void cs_f_majgeo(const cs_int_t     *ncel,
                        const cs_int_t     *ncelet,
                        const cs_int_t     *nfac,
                        const cs_int_t     *nfabor,
                        const cs_int_t     *nsom,
                        const cs_int_t     *lndfac,
                        const cs_int_t     *lndfbr,
                        const int64_t      *ncelgb,
                        const int64_t      *nfacgb,
                        const int64_t      *nfbrgb,
                        const int64_t      *nsomgb,
                        const cs_int_t     *nfml,
                        const cs_lnum_2_t   ifacel[],
                        const cs_int_t      ifabor[],
                        const cs_int_t      ifmfbr[],
                        const cs_int_t      ifmcel[],
                        const cs_int_t      ipnfac[],
                        const cs_int_t      nodfac[],
                        const cs_int_t      ipnfbr[],
                        const cs_int_t      nodfbr[],
                        const cs_int_t      isympa[],
                        const cs_int_t      isolid_0[],
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
 * \brief Define all mesh preprocessing operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_define(void)
{
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

  cs_preprocessor_data_read_mesh(cs_glob_mesh,
                                 cs_glob_mesh_builder);

  /* Join meshes / build periodicity links if necessary */

  cs_join_all(true);

  /* Insert boundaries if necessary */

  cs_gui_mesh_boundary(cs_glob_mesh);
  cs_user_mesh_boundary(cs_glob_mesh);

  cs_internal_coupling_preprocess(cs_glob_mesh);

  /* Initialize extended connectivity, ghost cells and other remaining
     parallelism-related structures */

  cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type);
  cs_mesh_update_auxiliary(cs_glob_mesh);

  /* Possible geometry modification */

  cs_gui_mesh_extrude(cs_glob_mesh);
  cs_user_mesh_modify(cs_glob_mesh);

  /* Discard isolated faces if present */

  cs_post_add_free_faces();
  cs_mesh_discard_free_faces(cs_glob_mesh);

  /* Smoothe mesh if required */

  cs_gui_mesh_smoothe(cs_glob_mesh);
  cs_user_mesh_smoothe(cs_glob_mesh);

  /* Triangulate warped faces if necessary */

  {
    double  cwf_threshold = -1.0;
    int  cwf_post = 0;

    cs_mesh_warping_get_defaults(&cwf_threshold, &cwf_post);

    if (cwf_threshold >= 0.0) {

      t1 = cs_timer_wtime();
      cs_mesh_warping_cut_faces(cs_glob_mesh, cwf_threshold, cwf_post);
      t2 = cs_timer_wtime();

      bft_printf(_("\n Cutting warped faces (%.3g s)\n"), t2-t1);

    }
  }

  /* Now that mesh modification is finished, save mesh if modified */

  cs_user_mesh_save(cs_glob_mesh); /* Disable or force */

  bool partition_preprocess = cs_partition_get_preprocess();
  if (cs_glob_mesh->modified > 0 || partition_preprocess) {
    if (partition_preprocess) {
      if (cs_glob_mesh->modified > 0)
        cs_mesh_save(cs_glob_mesh, cs_glob_mesh_builder, NULL, "mesh_output");
      else
        cs_mesh_to_builder(cs_glob_mesh, cs_glob_mesh_builder, true, NULL);
      cs_partition(cs_glob_mesh, cs_glob_mesh_builder, CS_PARTITION_MAIN);
      cs_mesh_from_builder(cs_glob_mesh, cs_glob_mesh_builder);
      cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type);
      cs_mesh_update_auxiliary(cs_glob_mesh);
    }
    else
      cs_mesh_save(cs_glob_mesh, NULL, NULL, "mesh_output");
  }

  /* Destroy the temporary structure used to build the main mesh */

  cs_mesh_builder_destroy(&cs_glob_mesh_builder);

  /* Renumber mesh based on code options */

  cs_user_numbering();

  cs_renumber_mesh(cs_glob_mesh);

  /* Initialize group classes */

  cs_mesh_init_group_classes(cs_glob_mesh);

  /* Print info on mesh */

  cs_mesh_print_info(cs_glob_mesh, _("Mesh"));

  /* Second pass to define internal coupling locators */
  cs_internal_coupling_map(cs_glob_mesh);

  /* Compute geometric quantities related to the mesh */

  bft_printf_flush();

  t1 = cs_timer_wtime();

  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);

  if (cs_glob_porous_model == 3) {
    cs_mesh_init_fluid_sections(cs_glob_mesh, cs_glob_mesh_quantities);
    cs_mesh_quantities_fluid_compute(cs_glob_mesh, cs_glob_mesh_quantities);
  }
  cs_mesh_bad_cells_detect(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_user_mesh_bad_cells_tag(cs_glob_mesh, cs_glob_mesh_quantities);
  t2 = cs_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  /* Initialize selectors and locations for the mesh */

  cs_mesh_init_selectors();
  cs_mesh_location_build(cs_glob_mesh, -1);
  cs_volume_zone_build_all(true);
  cs_boundary_zone_build_all(true);

  /* For debugging purposes */

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_mesh_dump(cs_glob_mesh);
  cs_mesh_quantities_dump(cs_glob_mesh, cs_glob_mesh_quantities);
#endif

  /* Re-enable writers disabled when entering this stage */

  cs_post_enable_writer(0);

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
              &(m->n_families),
              (const cs_lnum_2_t *)(m->i_face_cells),
              m->b_face_cells,
              m->b_face_family,
              m->cell_family,
              m->i_face_vtx_idx,
              m->i_face_vtx_lst,
              m->b_face_vtx_idx,
              m->b_face_vtx_lst,
              mq->b_sym_flag,
              mq->c_solid_flag,
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

END_C_DECLS
