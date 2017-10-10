/*============================================================================
 * Insert boundary cell layers into the mesh.
 *============================================================================*/

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

#include "fvm_io_num.h"

#include "cs_boundary_zone.h"
#include "cs_domain.h"

#include "cs_log.h"
#include "cs_mesh_builder.h"
#include "cs_mesh_deform.h"
#include "cs_mesh_extrude.h"
#include "cs_mesh_group.h"
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

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define extrusion parameters by boundary zone.
 *
 * The caller is responsible for freeing the n_layers, distribution, and
 * coord_shift arrays.
 *
 * \param[in]   n_zones         number of zones
 * \param[in]   zone_ids        ids of zones
 * \param[in]   zone_layers     number of extrusion layers per zone, or
 *                              NULL for default (1)
 * \param[in]   zone_thickness  thickness for each zone, or NULL for
 *                              default (based on neighboring cell
 *                              sizes)
 * \param[in]   zone_expansion  parameter for each zone, or NULL for
 *                              default (0.8)
 * \param[out]  n_layers        number of layers for each vertex
 * \param[out]  distribution    layer distribution for each vertex
 * \param[out]  coord_shift     coordinate shift for each vertex
 */
/*----------------------------------------------------------------------------*/

static void
_boundary_layer_define_by_zone(int              n_zones,
                               const int        zone_ids[],
                               const int        zone_layers[],
                               const double     zone_thickness[],
                               const float      zone_expansion[],
                               cs_lnum_t      **n_layers,
                               float          **distribution,
                               cs_coord_3_t   **coord_shift)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_lnum_t *_n_layers = NULL;
  float *_expansion = NULL;
  cs_coord_3_t *_coord_shift = NULL;

  /* Determine vertices to move */

  cs_real_t  *w = NULL;
  int        *c = NULL;

  BFT_MALLOC(_n_layers, m->n_vertices, cs_lnum_t);
  BFT_MALLOC(_expansion, m->n_vertices, float);
  BFT_MALLOC(_coord_shift, m->n_vertices, cs_coord_3_t);

  BFT_MALLOC(w, m->n_vertices, cs_real_t);
  BFT_MALLOC(c, m->n_vertices, int);

  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    _n_layers[i] = 0;
    _expansion[i] = 0;
    _coord_shift[i][0] = 0;
    _coord_shift[i][1] = 0;
    _coord_shift[i][2] = 0;
    w[i] = 0;
    c[i] = 0;
  }

  /* Global check for zone thickness specification */

  int z_thickness_spec = 1;
  cs_real_t *v_b_thickness = NULL;

  if (zone_thickness == NULL && coord_shift == NULL)
    z_thickness_spec = 0;
  cs_parall_min(1, CS_INT_TYPE, &z_thickness_spec);

  if (z_thickness_spec < 1) {
    BFT_MALLOC(v_b_thickness, m->n_vertices, cs_real_t);
    cs_mesh_quantities_b_thickness_v(m, mq, 4, v_b_thickness);
  }

  /* Loop on zones */

  for (int i = 0; i < n_zones; i++) {

    int z_layers = (zone_layers != NULL) ? zone_layers[i] : 1;
    float z_expansion = (zone_expansion != NULL) ? zone_expansion[i] : 0.8;
    double z_thickness = (zone_thickness != NULL) ? zone_thickness[i] : -1.;

    /* Now determine other parameters */

    cs_lnum_t z_id = zone_ids[i];
    const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);

    for (cs_lnum_t j = 0; j < z->n_faces; j++) {

      cs_lnum_t f_id = z->face_ids[j];
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      const cs_real_t *f_n = mq->b_face_normal + f_id*3;
      const cs_real_t f_s = cs_math_3_norm(f_n);

      if (v_b_thickness != NULL) {
        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          _n_layers[v_id] += z_layers;
          _expansion[v_id] += z_expansion;
          for (cs_lnum_t l = 0; l < 3; l++)
            _coord_shift[v_id][l] += v_b_thickness[v_id] * f_n[l];
          w[v_id] += f_s;
          c[v_id] += 1;
        }
      }
      else {
        for (cs_lnum_t k = s_id; k < e_id; k++) {
          cs_lnum_t v_id = m->b_face_vtx_lst[k];
          _n_layers[v_id] += z_layers;
          _expansion[v_id] += z_expansion;
          for (cs_lnum_t l = 0; l < 3; l++)
            _coord_shift[v_id][l] += z_thickness * f_n[l];
          w[v_id] += f_s;
          c[v_id] += 1;
        }
      }

    }

  }

  /* Handle parallelism */

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_LNUM_TYPE,
                         _n_layers);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_FLOAT,
                         _expansion);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         3,
                         true,
                         CS_COORD_TYPE,
                         _coord_shift);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         w);
    cs_interface_set_sum(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_INT_TYPE,
                         c);
  }

  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    if (c[i] > 0) {
      _n_layers[i] /= c[i];
      _expansion[i] /= c[i];
      for (cs_lnum_t l = 0; l < 3; l++)
        _coord_shift[i][l] /= w[i];
    }
  }

  BFT_FREE(c);
  BFT_FREE(w);

  /* Return values */

  *n_layers = _n_layers;

  float * _distribution;
  cs_lnum_t d_shift = 0;

  for (cs_lnum_t i = 0; i < m->n_vertices; i++)
    d_shift += (*n_layers)[i];

  BFT_MALLOC(_distribution, d_shift, float);

  d_shift = 0;
  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    int n_l = (*n_layers)[i];
    /* Compute distribution for each extruded vertex */
    if (n_l > 0) {
      float *_d = _distribution + d_shift;
      _d[0] = 1;
      for (cs_lnum_t l_id = 1; l_id < n_l; l_id++)
        _d[l_id] = _d[l_id-1]*_expansion[i];
      double d_tot = 0;
      for (cs_lnum_t l_id = 0; l_id < n_l; l_id++)
        d_tot += _d[l_id];
      _d[0] = 1./d_tot;
      for (cs_lnum_t l_id = 1; l_id < n_l - 1; l_id++)
        _d[l_id] = _d[l_id-1] + _d[l_id]/d_tot;
      _d[n_l-1] = 1.0;
      d_shift += n_l;
    }
  }

  *distribution = _distribution;

  BFT_FREE(_expansion);

  *coord_shift = _coord_shift;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert mesh boundary layers at end of computation.
 *
 * \param[in]  n_zones         number of zones
 * \param[in]  zone_ids        ids of zones
 * \param[in]  zone_layers     number of extrusion layers per zone, or
 *                             NULL for default (1)
 * \param[in]  zone_thickness  thickness for each zone, or NULL for
 *                             default (based on neighboring cell
 *                             sizes)
 * \param[in]  zone_expansion  parameter for each zone, or NULL for
 *                             default (0.8)
 * \param[in]  n_layers        optional specification of number of layers for
 *                             each vertex, or NULL
 * \param[in]  distribution    optional specification of layer distribution
 *                             for each vertex, or NULL
 * \param[in]  coord_shift     optional definition of coordinate shift for
 *                             each vertex, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_layer_insert(int                  n_zones,
                              const int            zone_ids[],
                              const int            zone_layers[],
                              const double         zone_thickness[],
                              const float          zone_expansion[],
                              const cs_lnum_t     *n_layers,
                              const float         *distribution,
                              const cs_coord_3_t  *coord_shift)
{
  cs_mesh_t *m = cs_glob_mesh;

  cs_timer_t t0 = cs_timer_time();

  /* Initialize some mesh structures if needed */

  cs_mesh_quantities_compute_preprocess(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_mesh_init_selectors();
  cs_mesh_location_build(cs_glob_mesh, -1);
  cs_boundary_zone_build_all(true);

  /* Determine vertices to move */

  const cs_lnum_t     *_n_layers = n_layers;
  const float         *_distribution = distribution;
  const cs_coord_3_t  *_coord_shift = coord_shift;

  cs_lnum_t     *_n_layers_tmp = NULL;
  float         *_distribution_tmp = NULL;
  cs_coord_3_t  *_coord_shift_tmp = NULL;

  int need_define_by_zone = 0;
  if (n_layers == NULL || distribution == NULL || coord_shift == NULL)
    need_define_by_zone = 1;
  cs_parall_max(1, CS_INT_TYPE, &need_define_by_zone);

  if (need_define_by_zone)
    _boundary_layer_define_by_zone(n_zones,
                                   zone_ids,
                                   zone_layers,
                                   zone_thickness,
                                   zone_expansion,
                                   &_n_layers_tmp,
                                   &_distribution_tmp,
                                   &_coord_shift_tmp);

  if (_n_layers == NULL)
    _n_layers = n_layers;

  if (_distribution == NULL)
    _distribution = distribution;

  if (_coord_shift == NULL)
    _coord_shift = coord_shift;

  if (sizeof(cs_coord_t) == sizeof(cs_real_t))
    cs_mesh_deform_prescribe_displacement(m->n_vertices,
                                          NULL,
                                          (const cs_real_3_t *)_coord_shift);
  else {
    cs_real_3_t *_c_shift;
    BFT_MALLOC(_c_shift, m->n_vertices, cs_real_3_t);
#   pragma omp parallel for if (m->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _c_shift[i][j] = _coord_shift[i][j];
    }
    cs_mesh_deform_prescribe_displacement(m->n_vertices,
                                          NULL,
                                          (const cs_real_3_t *)_c_shift);
    BFT_FREE(_c_shift);
  }

  /* Now deform mesh */

  cs_mesh_deform_solve_displacement();

  const cs_real_3_t *vd = cs_mesh_deform_get_displacement();

  for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
    m->vtx_coord[i*3]     += vd[i][0];
    m->vtx_coord[i*3 + 1] += vd[i][1];
    m->vtx_coord[i*3 + 2] += vd[i][2];
  }

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(time_count, cs_glob_domain->tcs, time_count);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO> Total runtime", time_count.wall_nsec*1e-9);

  /* Now add extrusion */

  {
    cs_lnum_t   n_selected_faces = 0;
    cs_lnum_t  *selected_faces = NULL;

    BFT_MALLOC(selected_faces, m->n_b_faces, cs_lnum_t);

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
      selected_faces[i] = -1;

    for (int i = 0; i < n_zones; i++) {
      cs_lnum_t z_id = zone_ids[i];
      const cs_boundary_zone_t  *z = cs_boundary_zone_by_id(z_id);
      for (cs_lnum_t j = 0; j < z->n_faces; j++)
        selected_faces[z->face_ids[j]] = 1;
    }

    for (cs_lnum_t i = 0; i < m->n_b_faces; i++) {
      if (selected_faces[i] > -1) {
        selected_faces[n_selected_faces] = i;
        n_selected_faces += 1;
      }
    }

    BFT_REALLOC(selected_faces, n_selected_faces, cs_lnum_t);

    /* Extrude selected boundary */

#if 0
    cs_mesh_extrude(m,
                    true, /* interior_gc */
                    n_selected_faces,
                    m->n_vertices,
                    selected_faces,
                    NULL,
                    _n_layers,
                    (const cs_coord_3_t *)_coord_shift,
                    _distribution);
#endif
    /* Free temporary memory */

    BFT_FREE(selected_faces);

  }

  BFT_FREE(_n_layers_tmp);
  BFT_FREE(_distribution_tmp);
  BFT_FREE(_coord_shift_tmp);

  m->modified = 1;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
