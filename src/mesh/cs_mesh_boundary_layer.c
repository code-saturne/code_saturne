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

  if (sizeof(cs_coord_t) == sizeof(cs_real_t))
    cs_mesh_deform_prescribe_displacement(m->n_vertices,
                                          NULL,
                                          (const cs_real_3_t *)e->coord_shift);
  else {
    cs_real_3_t *_c_shift;
    BFT_MALLOC(_c_shift, m->n_vertices, cs_real_3_t);
#   pragma omp parallel for if (m->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_vertices; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        _c_shift[i][j] = e->coord_shift[i][j];
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

  cs_mesh_extrude(m, e, interior_gc);

  m->modified = 1;
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
