/*============================================================================
 * User definitions of porous media.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_porosity.c
 *
 * \brief User definitions of porous media.
 *
 * See \subpage cs_porosity for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the porosity (volume factor \f$ \epsilon \f$
 *        when the porosity model is activated.
 *        (\ref cs_glob_porous_model > 0).
 *
 * This function is called at the beginning of the simulation only.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_porosity(void)
{
  /*!< [init_poro_mq] */
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

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
  /*!< [init_poro_mq] */

  /* Get the cell porosity field value */
  /*!< [init_poro_pro] */
  cs_real_t *cpro_porosi = cs_field_by_name("porosity")->val;
  /*!< [init_poro_pro] */

  /* First, set cell porosity value; the fluid cell volume will be
   * automatically deduced */

  /*!< [set_poro_cells_1] */
  for (cs_lnum_t cell_id = 0; cell_id < m->n_cells; cell_id++) {

    cs_real_t x = cell_cen[cell_id][0];

    if (x < 20.)
      cpro_porosi[cell_id] = 1.;
    else
      cpro_porosi[cell_id] = 0.5;
  }

  /* synchronize for use in fluid face factor calculation below */
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_field_synchronize(cs_field_by_name("porosity"),
                       halo_type);

  /*!< [set_poro_cells_1] */

  /* Set interior face values */

  /*!< [set_poro_i_faces_1] */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_real_t x = i_face_cog[face_id][0];

    cs_real_t face_porosity = 1.;
    if (x > 19.9999)
      face_porosity = 0.5;

    for (int i = 0; i < 3; i++)
      i_f_face_normal[face_id][i] = face_porosity * i_face_normal[face_id][i];

    mq->i_f_face_surf[face_id] = cs_math_3_norm(i_f_face_normal[face_id]);

  }
  /*!< [set_poro_i_faces_1] */

  /* Set boundary face values */

  /*!< [set_poro_b_faces_1] */
  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_real_t x = b_face_cog[face_id][0];

    cs_real_t face_porosity = 1.;
    if (x > 19.9999)
      face_porosity = 0.5;

    for (int i = 0; i < 3; i++)
      b_f_face_normal[face_id][i] = face_porosity * b_face_normal[face_id][i];

    mq->b_f_face_surf[face_id] = cs_math_3_norm(b_f_face_normal[face_id]);

  }
  /*!< [set_poro_b_faces_1] */

  /* Four set face factor */
  if (mq->i_f_face_factor != NULL) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      cs_real_t face_porosity =
        CS_MAX(
            i_f_face_surf[face_id] / i_face_surf[face_id],
            cs_math_epzero);

      mq->i_f_face_factor[face_id][0] = cpro_porosi[ii] / face_porosity;
      mq->i_f_face_factor[face_id][1] = cpro_porosi[jj] / face_porosity;
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
