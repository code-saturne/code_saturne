/*============================================================================
 * Handling of porosity due to particle deposition.
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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_porous_model.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_porosity.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*=============================================================================
 * Local Enumeration definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
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
 * \brief Compute interior face porosity for Lagrangian internal deposition.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_porosity(void)
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;
  const cs_lagr_model_t *lagr_model = cs_glob_lagr_model;

  cs_lagr_particle_set_t  *particles = cs_glob_lagr_particle_set;

  cs_lagr_internal_condition_t *internal_conditions
    = cs_glob_lagr_internal_conditions;

  if (internal_conditions == NULL)
    return;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;

  cs_real_t *covered_surface = NULL;
  BFT_MALLOC(covered_surface, mesh->n_cells_with_ghosts, cs_real_t);

  /* Initialization */
  for (cs_lnum_t cell_id = 0; cell_id < n_cells_ext; cell_id++)
    covered_surface[cell_id] = 0.;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces ; face_id++) {
    /* Internal face flagged as internal deposition */
    if (internal_conditions->i_face_zone_id[face_id] >= 0) {
      for (cs_lnum_t j = 0; j < 3; j++)
        fvq->i_f_face_normal[3*face_id+j] = fvq->i_face_normal[3*face_id+j];
    }
  }

  if (lagr_model->deposition == 1) {

    for (cs_lnum_t ip = 0; ip < particles->n_particles; ip++) {

      if (cs_lagr_particles_get_flag(particles, ip,
                                     CS_LAGR_PART_IMPOSED_MOTION)) {

        cs_lnum_t cell_id
          = cs_lagr_particles_get_lnum(particles, ip, CS_LAGR_CELL_ID);

        cs_real_t diam2
          = cs_math_pow2(cs_lagr_particles_get_real(particles, ip,
                                                    CS_LAGR_DIAMETER));

        covered_surface[cell_id] += cs_math_pi * 0.25 * diam2
          * cs_lagr_particles_get_real(particles, ip, CS_LAGR_FOULING_INDEX)
          * cs_lagr_particles_get_real(particles, ip, CS_LAGR_STAT_WEIGHT);

      }
    }
  }

  /* Synchronization */
  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, covered_surface);

  /* Compute fluid section and clip it to 0 if negative */
  for (cs_lnum_t face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    /* Internal face flagged as internal deposition */
    if (internal_conditions->i_face_zone_id[face_id] >= 0) {
      cs_lnum_t cell_id1 = mesh->i_face_cells[face_id][0];;
      cs_lnum_t cell_id2 = mesh->i_face_cells[face_id][1];
      /* Remove from the particle area from fluid section */
      for (cs_lnum_t id = 0; id < 3; id++)
        fvq->i_f_face_normal[3*face_id + id]
          -= (covered_surface[cell_id1] + covered_surface[cell_id2])
             * fvq->i_face_normal[3*face_id + id]
             / fvq->i_face_surf[face_id];

      /* If S_fluid . S is negative, that means we removed too much surface
       * from fluid surface */
      cs_real_t temp
        = cs_math_3_dot_product(fvq->i_f_face_normal + 3*face_id,
                                fvq->i_face_normal + 3*face_id);
      if (temp <= 0.) {
        for (cs_lnum_t j = 0; j < 3; j++)
          fvq->i_f_face_normal[3*face_id+j] = 0.;
      }
      fvq->i_f_face_surf[face_id]
        = cs_math_3_norm(fvq->i_f_face_normal + 3*face_id);
    }
  }

  BFT_FREE(covered_surface);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
