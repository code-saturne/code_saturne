/*============================================================================
 * User definition of boundary conditions.
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
#include <stdio.h>
#include <string.h>


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in]     nvar          total number of variable BC's
 * \param[in]     bc_type       boundary face types
 * \param[in]     icodcl        boundary face code
 *                                - 1  -> Dirichlet
 *                                - 2  -> convective outlet
 *                                - 3  -> flux density
 *                                - 4  -> sliding wall and u.n=0 (velocity)
 *                                - 5  -> friction and u.n=0 (velocity)
 *                                - 6  -> roughness and u.n=0 (velocity)
 *                                - 9  -> free inlet/outlet (velocity)
 *                                inflowing possibly blocked
 * \param[in]     rcodcl        boundary condition values
 *                                rcodcl(3) = flux density value
 *                                (negative for gain) in W/m2
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable id "ivar" and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *  icodcl[ivar*n_b_faces + face_id]
 *
 * - Dirichlet value defined at:
 *   rcodcl[ivar*n_b_faces + face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   rcodcl[(ivar +   nvar)*n_b_faces + face_id] = interior exchange
 *
 * - Flux density defined at:
 *   rcodcl[(ivar + 2*nvar)*n_b_faces + face_id];
 *
 * For a given field f, the "ivar" variable id may be obtained as follows:
 *   int ivar = cs_field_get_key_int(f, cs_field_key_id("variable_id")) - 1;
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         bc_type[],
                            int         icodcl[],
                            cs_real_t   rcodcl[])
{
  /* Shorter aliases to global variables */

  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_zone_t *z = cs_boundary_zone_by_name("Inlet");

  /* Get access to boundary condition arrays for temperature field */

  const cs_field_t *f = cs_field_by_name("temperature");
  cs_lnum_t ivar
    = cs_field_get_key_int(f, cs_field_key_id("variable_id")) - 1;

  /* Prescribe temperature (Dirichlet) value based on time */

  if (ts->t_cur < 3.8) {
    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      cs_lnum_t face_id = z->elt_ids[i];
      rcodcl[ivar*n_b_faces + face_id] = 20. + 100.*ts->t_cur;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      cs_lnum_t face_id = z->elt_ids[i];
      rcodcl[ivar*n_b_faces + face_id] = 400.;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
