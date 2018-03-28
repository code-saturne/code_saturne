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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert mesh boundary layers.
 *
 * \param[in, out]  m             mesh
 * \param[in]       e             extrusion vector definitions
 * \param[in]       interior_gc   if true, maintain group classes of
 *                                interior faces previously on boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_layer_insert(cs_mesh_t                        *m,
                              const cs_mesh_extrude_vectors_t  *e,
                              bool                              interior_gc)
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

  cs_domain_t  *domain = cs_domain_create();
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

  /* Now prescribe displacement (invert extrusion direction) */

  cs_real_3_t *_c_shift;
  BFT_MALLOC(_c_shift, e->n_vertices, cs_real_3_t);
# pragma omp parallel for if (m->n_vertices > CS_THR_MIN)
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

  cs_mesh_deform_finalize();

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  CS_TIMER_COUNTER_ADD(time_count, domain->tcs, time_count);

  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO> Total runtime", time_count.wall_nsec*1e-9);
  cs_domain_free(&domain);

  cs_mesh_extrude(m, e, interior_gc);

  cs_mesh_quantities_free_all(mq);

  m->modified = 1;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
