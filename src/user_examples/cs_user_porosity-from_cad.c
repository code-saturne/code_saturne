/*============================================================================
 * User definitions of porous media.
 *============================================================================*/

/* VERS */

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

/*!< [user_poro_cad_proto] */
#include "cs_cad_intersect.h"
/*!< [user_poro_cad_proto] */

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_porosity-from_cad.c
 *
 * \brief User definitions of porous media.
 *
 * See \ref cs_porosity for examples.
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
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_porosity(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

/*!< [user_poro_cad_zone] */
  const cs_zone_t  *z_poro = cs_volume_zone_by_name("Zone_1");
/*!< [user_poro_cad_zone] */

/*!< [user_poro_cad_init] */
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_i_faces = m->n_i_faces;
  const cs_lnum_t n_b_faces = m->n_b_faces;

  /* Temporary arrays for face porosity */

  cs_real_t  *cell_porosity = cs_field_by_name("porosity")->val;

  cs_real_t  *i_face_porosity, *b_face_porosity;

  BFT_MALLOC(i_face_porosity, m->n_i_faces, cs_real_t);
  BFT_MALLOC(b_face_porosity, m->n_b_faces, cs_real_t);

  for (cs_lnum_t i = 0; i < m->n_i_faces; i++)
    i_face_porosity[i] = 1;
  for (cs_lnum_t i = 0; i < m->n_b_faces; i++)
    b_face_porosity[i] = 1;
  /*!< [user_poro_cad_init] */

  /* Intersect with CAD to obtain new quantities */

  /*!< [user_poro_cad_intersect] */
  cs_cad_intersect(m,
                   "Compound_1.step",
                   CS_CAD_INTERSECT_CUT,
                   z_poro->n_elts,
                   z_poro->elt_ids,
                   cell_porosity,
                   NULL,  /* modified cell centers */
                   i_face_porosity,
                   NULL,  /* modified interior face centers */
                   b_face_porosity,
                   NULL);  /* modified boundary face centers */
  /*!< [user_poro_cad_intersect] */

  /*!< [user_poro_cad_quantities] */
  /* synchronize ghost cells */

  if (m->halo != NULL)
    cs_halo_sync_var(m->halo, CS_HALO_EXTENDED, cell_porosity);

  /* Set interior and boundary face values */

  if (mq->i_f_face_surf != mq->i_face_surf) {

    const cs_real_3_t *restrict i_face_normal
      = (const cs_real_3_t *restrict)mq->i_face_normal;
    cs_real_3_t *restrict i_f_face_normal
      = (cs_real_3_t *restrict)mq->i_f_face_normal;
    cs_real_t *i_f_face_surf = mq->i_f_face_surf;

    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_real_t f_porosity = i_face_porosity[f_id];
      for (cs_lnum_t i = 0; i < 3; i++)
        i_f_face_normal[f_id][i] = f_porosity * i_face_normal[f_id][i];
      i_f_face_surf[f_id] = cs_math_3_norm(i_f_face_normal[f_id]);
    }

  }

  if (mq->b_f_face_surf != mq->b_face_surf) {

    const cs_real_3_t *restrict b_face_normal
      = (const cs_real_3_t *restrict)mq->b_face_normal;
    cs_real_3_t *restrict b_f_face_normal
      = (cs_real_3_t *restrict)mq->b_f_face_normal;
    cs_real_t *b_f_face_surf = mq->b_f_face_surf;

    for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
      cs_real_t f_porosity = b_face_porosity[f_id];
      for (cs_lnum_t i = 0; i < 3; i++)
        b_f_face_normal[f_id][i] = f_porosity * b_face_normal[f_id][i];
      b_f_face_surf[f_id] = cs_math_3_norm(b_f_face_normal[f_id]);
    }

  }

  BFT_FREE(i_face_porosity);
  BFT_FREE(b_face_porosity);

  /* Four set face factor */

  if (mq->i_f_face_factor != NULL) {
    const cs_lnum_2_t *i_face_cells
      = (const cs_lnum_2_t *)m->i_face_cells;
    const cs_real_t *i_face_surf = mq->i_face_surf;
    const cs_real_t *i_f_face_surf = mq->i_f_face_surf;

    for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {
      cs_lnum_t ii = i_face_cells[f_id][0];
      cs_lnum_t jj = i_face_cells[f_id][1];

      cs_real_t face_porosity
        = CS_MAX(i_f_face_surf[f_id] / i_face_surf[f_id],
                 cs_math_epzero);

      mq->i_f_face_factor[f_id][0] = cell_porosity[ii] / face_porosity;
      mq->i_f_face_factor[f_id][1] = cell_porosity[jj] / face_porosity;
    }
  }
  /*!< [user_poro_cad_quantities] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
