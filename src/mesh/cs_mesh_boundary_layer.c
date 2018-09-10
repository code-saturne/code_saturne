/*============================================================================
 * Insert boundary cell layers into the mesh.
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_boundary_zone.h"
#include "cs_cdo_main.h"
#include "cs_domain.h"

#include "cs_log.h"
#include "cs_interface.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_deform.h"
#include "cs_mesh_extrude.h"
#include "cs_mesh_group.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_boundary_layer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Temporary pointers (for callback for private mesh location) */

static const cs_mesh_extrude_vectors_t  *_extrude_vectors = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of boundary faces attached which are associated
 *         to vertices with boundary layer insertion.
 *
 * If non-empty and not containing all elements, a list of elements
 * of the parent mesh belonging to the location should be allocated
 * (using BFT_MALLOC) and defined by this function when called.
 * This list's lifecycle is then managed by the mesh location object.
 *
 * \param [in]   m            pointer to associated mesh structure.
 * \param [in]   location_id  id of associated location.
 * \param [out]  n_elts       number of selected elements
 * \param [out]  elt_list     list of selected elements.
 */
/*----------------------------------------------------------------------------*/

static void
_transfer_bl_faces_selection(void              *input,
                             const cs_mesh_t   *m,
                             int                location_id,
                             cs_lnum_t         *n_elts,
                             cs_lnum_t        **elt_ids)
{
  CS_UNUSED(input);
  CS_UNUSED(m);
  CS_UNUSED(location_id);

  if (_extrude_vectors != NULL) {
    const cs_lnum_t _n_sel_faces = _extrude_vectors->n_faces;

    *n_elts = _n_sel_faces;

    BFT_MALLOC(*elt_ids, _n_sel_faces, cs_lnum_t);
    memcpy(*elt_ids,
           _extrude_vectors->face_ids,
           _n_sel_faces*sizeof(cs_lnum_t));
  }
  else {
    *n_elts = 0;
    *elt_ids = NULL;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Flag vertices for limiter.
 *
 * \param[in]   m                mesh
 * \param[in]   cell_vol_cmp     comparative cell volume (< 0 for limit)
 * \param[out]  vtx_flag         vertex flag (0 for unlimited, 1 for limited)
 */
/*----------------------------------------------------------------------------*/

static void
_flag_vertices_for_limiter(const cs_mesh_t  *m,
                           const cs_real_t  *cell_vol_cmp,
                           char             *vtx_flag)
{
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t n_vertices = m->n_vertices;

  /* Flag vertices adjacent to cells with bad volumes */

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    vtx_flag[i] = 0;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    bool flag_vertices = false;
    cs_lnum_t c_id0 = m->i_face_cells[f_id][0];
    cs_lnum_t c_id1 = m->i_face_cells[f_id][0];
    if (c_id0 > -1 && c_id0 < n_cells) {
      if (cell_vol_cmp[c_id0] <= 0)
        flag_vertices = true;
    }
    if (c_id1 > -1 && c_id1 < n_cells) {
      if (cell_vol_cmp[c_id1] <= 0)
        flag_vertices = true;
    }
    if (flag_vertices) {
      cs_lnum_t s_id = m->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->i_face_vtx_idx[f_id+1];
      for (cs_lnum_t i = s_id; i < e_id; i++)
        vtx_flag[m->i_face_vtx_lst[i]] = 1;
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    bool flag_vertices = false;
    cs_lnum_t c_id0 = m->b_face_cells[f_id];
    if (c_id0 > -1) {
      if (cell_vol_cmp[c_id0] <= 0)
        flag_vertices = true;
    }
    if (flag_vertices) {
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t i = s_id; i < e_id; i++)
        vtx_flag[m->b_face_vtx_lst[i]] = 1;
    }
  }

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_max(m->vtx_interfaces,
                         n_vertices,
                         1,    /* stride */
                         true, /* interlace */
                         CS_CHAR,
                         vtx_flag);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Limit extrusion vector definitions.
 *
 * \param[in]       vtx_flag      per vertex reduction indicator flag
 *                                (size: m->n_vertices)
 *                                interior faces previously on boundary
 * \param[in, out]  e             extrusion vector definitions
 *
 * \return  local number of vertices at which extrusion is reduced
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_extrude_vector_limit(const char                 *vtx_flag,
                      cs_mesh_extrude_vectors_t  *e)
{
  cs_lnum_t n_limited = 0;

  if (e->distribution_idx != NULL) {

    cs_lnum_t n = e->distribution_idx[0];

    for (cs_lnum_t i = 0; i < e->n_vertices; i++) {
      cs_lnum_t s_id = e->distribution_idx[i];
      cs_lnum_t e_id = e->distribution_idx[i+1];
      cs_lnum_t n_layers = e->n_layers[i];
      if (vtx_flag[i] > 0 && n_layers > 0) {
        cs_real_t r = 0;
        if (n_layers > 1) {
          r = e->distribution[e_id-2];
          for (cs_lnum_t j = s_id; j < e_id-1; j++) {
            e->distribution[j] /= r;
            if (e->distribution[j] > 1) /* in case of truncation error */
              e->distribution[j] = 1;
          }
        }
        n_layers -= 1;
        e->n_layers[i] = n_layers;
        for (cs_lnum_t j = 0; j < 3; j++)
          e->coord_shift[i][j] *= r;
        n_limited += 1;
      }
      e->distribution_idx[i] = n;
      for (cs_lnum_t j = 0; j < n_layers; j++)
        e->distribution[n++] = e->distribution[s_id+j];
    }
    e->distribution_idx[e->n_vertices] = n;
  }
  else { /* if (distribution_idx == NULL) */

    for (cs_lnum_t i = 0; i < e->n_vertices; i++) {
      if (vtx_flag[i] > 0 && e->n_layers[i] > 0) {
        cs_lnum_t n_layers = e->n_layers[i];
        double r = (double)(n_layers - 1) / (double)n_layers;
        n_layers -= 1;
        e->n_layers[i] = n_layers;
        if (n_layers == 0)
          r = 0;
        for (cs_lnum_t j = 0; j < 3; j++)
          e->coord_shift[i][0] *= r;
        n_limited += 1;
      }
    }

  }

  return n_limited;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert mesh boundary layers.
 *
 * \param[in, out]  m                  mesh
 * \param[in, out]  e                  extrusion vector definitions
 * \param[in]       min_volume_factor  cell volume multiplier threshold:
 *                                     extrusion is reduced on vertices
 *                                     adjacent to cells whose volume is
 *                                     reduced below this; < 0 to ignore
 * \param[in]       interior_gc        if true, maintain group classes of
 *                                     interior faces previously on boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_layer_insert(cs_mesh_t                  *m,
                              cs_mesh_extrude_vectors_t  *e,
                              cs_real_t                   min_volume_factor,
                              bool                        interior_gc)
{
  cs_timer_t t0 = cs_timer_time();

  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Ensure mesh quantities and locations are up to date in case
     of call during preprocessing stage */

  {
    cs_mesh_quantities_compute_preprocess(m, mq);

    cs_mesh_init_selectors();
    cs_mesh_location_build(m, -1);
  }

  /* Define associated boundary zone */

  _extrude_vectors = e;

  const char *z_name = "_boundary_layer_insert";
  int z_id[1] = {-1};

  {
    const cs_zone_t  *z
      = cs_boundary_zone_by_name_try(z_name);
    if (z != NULL) {
      z_id[0] = z->id;
      assert(z->type & CS_BOUNDARY_ZONE_PRIVATE);
    }
  }
  if (z_id[0] < 0)
    z_id[0] = cs_boundary_zone_define_by_func(z_name,
                                              _transfer_bl_faces_selection,
                                              NULL,
                                              CS_BOUNDARY_ZONE_PRIVATE);

  cs_boundary_zone_build_private(z_id[0]);

  /* Local activation of CDO module if required */

  // cs_domain_t  *domain = cs_domain_create();
  cs_domain_t  *domain = cs_glob_domain;
  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_WITH_FV);

  cs_mesh_deform_define_dirichlet_bc_zones(1, z_id);

  cs_mesh_deform_activate();

  cs_cdo_initialize_setup(domain);

  cs_cdo_initialize_structures(domain, m, mq);

  /* Deactive logging and visualization for deformation
     fields, as they are reset to 0 anyways after extrusion */

  const char *eq_name[] = {"mesh_deform_x", "mesh_deform_y", "mesh_deform_z"};
  for (int i = 0; i < 3; i++) {
    cs_field_t *f = cs_field_by_name_try(eq_name[i]);
    if (f != NULL) {
      cs_field_set_key_int(f, cs_field_key_id("log"), 0);
      cs_field_set_key_int(f, cs_field_key_id("post_vis"), 0);
    }
  }

  /* Compute or access reference volume for displacement limiter */

  const cs_lnum_t n_cells_ini = m->n_cells;

  const cs_real_t *cell_vol_ref = cs_glob_mesh_quantities->cell_vol;

  bool compute_displacement = true;

  while (compute_displacement) {

    /* Now prescribe displacement (invert extrusion direction) */

    cs_real_3_t *_c_shift;
    BFT_MALLOC(_c_shift, e->n_vertices, cs_real_3_t);
#   pragma omp parallel for if (m->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < e->n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _c_shift[i][j] = - e->coord_shift[i][j];
    }
    cs_mesh_deform_prescribe_displacement(e->n_vertices,
                                          e->vertex_ids,
                                          (const cs_real_3_t *)_c_shift);
    BFT_FREE(_c_shift);

    /* Now deform mesh */

    cs_mesh_deform_solve_displacement(domain);

    _extrude_vectors = NULL;

    const cs_real_3_t *vd = cs_mesh_deform_get_displacement();

    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      m->vtx_coord[i*3]     += vd[i][0];
      m->vtx_coord[i*3 + 1] += vd[i][1];
      m->vtx_coord[i*3 + 2] += vd[i][2];
    }

    /* Check deformation is acceptable */

    compute_displacement = false;

    if (min_volume_factor > 0 && min_volume_factor < 1) {

      cs_gnum_t  counts[3] = {0, 0, 0};

      cs_real_t *cell_vol_cmp = cs_mesh_quantities_cell_volume(m);

      for (cs_lnum_t i = 0; i < n_cells_ini; i++) {
        if (cell_vol_cmp[i] <= 0) {
          cell_vol_cmp[i] = -1;
          counts[0] += 1;
        }
        else if (cell_vol_cmp[i] < cell_vol_ref[i]*min_volume_factor) {
          cell_vol_cmp[i] = -1;
          counts[1] += 1;
        }
      }

      const cs_lnum_t n_vertices = m->n_vertices;

      char *vtx_flag;
      BFT_MALLOC(vtx_flag, n_vertices, char);

      /* Flag vertices adjacent to cells with bad volumes */

      _flag_vertices_for_limiter(m,
                                 cell_vol_cmp,
                                 vtx_flag);

      BFT_FREE(cell_vol_cmp);

      /* Now adjust extrusion vectors structure,
         removing a layer at flagged vertices */

      counts[2] = _extrude_vector_limit(vtx_flag, e);

      BFT_FREE(vtx_flag);

      cs_parall_sum(3, CS_GNUM_TYPE, counts);

      if (counts[2] > 0) {

        bft_printf
          (_("\nBoundary layer insertion:\n"
             "  %llu cells would have a negative volume\n"
             "  %llu cells would have a volume reduced by more than %g\n"
             "    (which is the user-defined threshold)\n"
             "  reducing insertion at adjacent boundary vertices.\n"),
           (unsigned long long)counts[0], (unsigned long long)counts[1],
           min_volume_factor);

        compute_displacement = true;

      }
      else if (counts[0] > 1) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf
          (_("%llu cells would have a negative volume after boundary insertion\n"
             "but none of these are adjacent to an inserted boundary.\n"
             "Unable to detemine appropriate insertion limitation."),
           (unsigned long long)counts[0]);
      }

      if (counts[0] > 1 || counts[2] > 1) {
        for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
          m->vtx_coord[i*3]     -= vd[i][0];
          m->vtx_coord[i*3 + 1] -= vd[i][1];
          m->vtx_coord[i*3 + 2] -= vd[i][2];
        }
      }

    } /* end of displacements computation and checking loop */

  }

  cs_mesh_deform_finalize();

  cell_vol_ref = NULL;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(time_count, domain->tcs, time_count);

  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO> Total runtime", time_count.wall_nsec*1e-9);
  cs_cdo_finalize(domain);

  cs_mesh_extrude(m, e, interior_gc);

  cs_mesh_quantities_free_all(mq);

  m->modified = 1;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
