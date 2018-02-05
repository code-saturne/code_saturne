/*============================================================================
 * User definitions of porous media.
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

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"
#include "cs_volume_zone.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_porosity.c
 *
 * \brief This function computes the porosity (volume factor \f$ \epsilon \f$
 * when porosity module is activated
 * (iporos greater than 1 in cs_user_parameters.f90).
 *
 * See \subpage cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the beginning of the simulation only.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_porosity(void)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Get the cell porosity field value */
  cs_real_t *cpro_porosi = cs_field_by_name("porosity")->val;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)m->i_face_cells;

  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)mq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)mq->b_face_cog;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)mq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)mq->i_face_normal;
  cs_real_3_t *restrict i_f_face_normal
    = (cs_real_3_t *restrict)mq->i_f_face_normal;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)mq->b_face_normal;
  cs_real_3_t *restrict b_f_face_normal
    = (cs_real_3_t *restrict)mq->b_f_face_normal;

  const cs_real_t *i_f_face_surf = mq->i_f_face_surf;
  const cs_real_t *i_face_surf = mq->i_face_surf;


  /* First set cell porosity value, the the fluid cell volume will be
   * automatically deduced */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

    cs_real_t x = cell_cen[cell_id][0];

    if (x < 20.)
      cpro_porosi[cell_id] = 1.;
    else
      cpro_porosi[cell_id] = 0.5;
  }

  /* Second set interior face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_real_t x = i_face_cog[face_id][0];

    cs_real_t face_porosity = 1.;
    if (x > 19.9999)
      face_porosity = 0.5;

    for (int i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

    mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

  }

  /* Third set boundary face values */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_real_t x = b_face_cog[face_id][0];

    cs_real_t face_porosity = 1.;
    if (x > 19.9999)
      face_porosity = 0.5;

    for (int i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

    mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
