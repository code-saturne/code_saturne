/*============================================================================
 * User definition of boundary conditions.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_elec_model.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_quantities.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_time_step.h"
#include "cs_selector.h"

#if defined(HAVE_MEDCOUPLING_LOADER)
#include "cs_medcoupling_remapper.hxx"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
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
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(int         nvar,
                            int         bc_type[],
                            int         icodcl[],
                            cs_real_t   rcodcl[])
{
  CS_UNUSED(nvar);
  CS_UNUSED(bc_type);

  /* TEST TO REMOVE FOR USAGE */
  if (true)
    return;

  /* Getting the value of a user defined parameter t_inlet */
  cs_real_t t_bnd = cs_notebook_parameter_value_by_name("t_inlet");

  /* Variables needed for boundary condition sub-selection */
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const int keyvar = cs_field_key_id("variable_id");
  cs_field_t *scalar = cs_field_by_name_try("scalar1");
  int ivar = cs_field_get_key_int(scalar, keyvar) - 1;

  cs_lnum_t  nelts = 0;
  cs_lnum_t *lstelt = NULL;

  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list("inlet", &nelts, lstelt);

  for (cs_lnum_t ielt = 0; ielt < nelts; ielt++) {
    cs_lnum_t face_id = lstelt[ielt];

    icodcl[ivar*n_b_faces + face_id] = 1;
    rcodcl[ivar*n_b_faces + face_id] = t_bnd;
  }

  BFT_FREE(lstelt);

  /* Use boundary zones */
  const cs_zone_t *zone = cs_boundary_zone_by_name("grid_inlet");

  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {
    cs_lnum_t face_id = zone->elt_ids[i];

    /* Dirichlet */
    icodcl[ivar*n_b_faces + face_id] = 1;
    rcodcl[ivar*n_b_faces + face_id] = 0.;
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
