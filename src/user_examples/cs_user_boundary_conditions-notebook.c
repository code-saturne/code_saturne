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

#include "cs_headers.h"

#if defined(HAVE_MEDCOUPLING_LOADER)
#include "cs_medcoupling_remapper.hxx"
#endif

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

  /* Get a user parameter defined in the GUI notebook */
  /*! [user_defined_param] */
  cs_real_t t_bnd = cs_notebook_parameter_value_by_name("t_inlet");
  /*! [user_defined_param] */


  /* Define some variables needed for the boundary condition specification */

  /*! [bc_param] */
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;


  const int keyvar = cs_field_key_id("variable_id");
  cs_field_t *scalar = cs_field_by_name_try("scalar1");
  int iscal = cs_field_get_key_int(scalar, keyvar) - 1;

  cs_lnum_t  nelts = 0;
  cs_lnum_t *lstelt = NULL;

  /*! [bc_param] */

  /*! [apply_bc] */
  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list("inlet", &nelts, lstelt);

  for (cs_lnum_t ielt = 0; ielt < nelts; ielt++) {
    cs_lnum_t f_id = lstelt[ielt];

    icodcl[iscal*n_b_faces + f_id] = 1;
    rcodcl[iscal*n_b_faces + f_id] = t_bnd;
  }

  BFT_FREE(lstelt);
  /*! [apply_bc] */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
