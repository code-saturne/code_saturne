/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <stdio.h>

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
 * \file cs_user_extra_operations-boundary_forces.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function).
 *
 * This is an example of cs_user_extra_operations.c which computes the total
 * force on a boundary zone.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  /*! [boundary_forces_ex1] */
  {
    const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

    cs_field_t *b_forces = cs_field_by_name_try("boundary_forces");

    if (b_forces != NULL) {
      cs_real_3_t total_b_forces = {0., 0., 0.};
      cs_real_3_t *bpro_forces = (cs_real_3_t *) b_forces->val;

      /* get zone from its name, here "selected_wall" */
      const cs_zone_t *zn = cs_boundary_zone_by_name("selected_wall");

      for (cs_lnum_t e_id = 0; e_id < zn->n_elts; e_id++) {
        cs_lnum_t face_id = zn->elt_ids[e_id];
        for (cs_lnum_t i = 0; i < 3; i++)
          total_b_forces[i] += bpro_forces[face_id][i];
      }

      /* parallel sum */
      cs_parall_sum(3, CS_REAL_TYPE, total_b_forces);
    }
  }
  /*! [boundary_forces_ex1] */

  /*! [boundary_forces_ex2] */
  {
    const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;
    const cs_real_3_t *b_f_face_normal =
      (cs_real_3_t *)domain->mesh_quantities->b_f_face_normal;

    cs_real_3_t total_b_p_forces = {0., 0., 0.};

    /* get zone from its name, here "selected_wall" */
    const cs_zone_t *zn = cs_boundary_zone_by_name("selected_wall");

    /* compute static pressure on selected boundary faces */
    cs_real_t *p_b_val;
    BFT_MALLOC(p_b_val, zn->n_elts, cs_real_t);
    cs_post_b_pressure(zn->n_elts, zn->elt_ids, p_b_val);

    for (cs_lnum_t e_id = 0; e_id < zn->n_elts; e_id++) {
      cs_lnum_t face_id = zn->elt_ids[e_id];
      for (cs_lnum_t i = 0; i < 3; i++)
        total_b_p_forces[i] += p_b_val[e_id]*b_f_face_normal[face_id][i];
    }

    BFT_FREE(p_b_val);

    /* parallel sum */
    cs_parall_sum(3, CS_REAL_TYPE, total_b_p_forces);
  }
  /*! [boundary_forces_ex2] */
}

END_C_DECLS
