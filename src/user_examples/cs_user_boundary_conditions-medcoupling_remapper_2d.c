/*============================================================================
 * User definition of boundary conditions.
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
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(bc_type);

  /*! [loc_var_def] */

  /* Variables needed for boundary condition sub-selection */

  /* MEDCoupling Remapper structure:
   * ------------------------------- */

  /* Number of fields to interpolate from the MED file */
  const int  nremapper_fields = 1;

  /* Names of the fields to read */
  const char  **field_names = NULL;
  BFT_MALLOC(field_names, nremapper_fields, const char *);
  field_names[0] = "TEMPERATURE";

  /*! [loc_var_def] */

  /*! [init] */
  /* Indexes needed to read the time step from the
   * file (-1, -1 if only one exists) */
  int it0 = -1;
  int it1 = -1;
  /*! [init] */

  /*! [remapper] */

  /* We request a remapper with a given name. If it does not exist,
   * the function returns a NULL pointer. */
  cs_medcoupling_remapper_t *r
    = cs_medcoupling_remapper_by_name_try("scalar_bc");

  /* If the returned pointer is NULL (first call), we create the
   * corresponding remapper */

  if (r == NULL) {

    /* Space dimension of the elements (2 for faces, 3 for cells) */
    int elts_dim = 2;

    /* Path to file */
    const char file_name[] = "/home/myname/study/2Dmap_Tfluid.med";

    /* The remapper is created. We retrieve its id from the function.
     * The function inputs are:
     * 1) Name of the remapper
     * 2) dimension of the mesh elements
     * 3) selection criteria for the boundary condition zone
     * 4) path to the med file
     * 5) number of fields to interpolate
     * 6) names of the fields to interpolate
     * 7 + 8) time iteration index and order */
    int r_id = cs_medcoupling_remapper_initialize("scalar_bc",
                                                  elts_dim,
                                                  "inlet",
                                                  file_name,
                                                  nremapper_fields,
                                                  field_names,
                                                  it0,
                                                  it1);

    /* Retrieve the pointer */
    r = cs_medcoupling_remapper_by_id(r_id);

    /* We create the interpolation matrix => Here it is only called once
     * since the mesh is not moving */
    cs_medcoupling_remapper_setup(r);

  }
  /*! [remapper] */

  /* If the med data needs for a translation or rotation for the geometrical
   * superposition with the target code_saturne mesh: */

  /*! [trans_rota] */
  if (false) {
    /* Translation using a tranlsation vector. Here it is (1, 0, 0) */
    cs_real_t translation_vector[3] = {1.0, 0.0, 0.0};
    cs_medcoupling_remapper_translate(r, translation_vector);

    /* Rotation using an invariant point, the rotation axis and rotation angle
       Here, center is O=(0,0,0) and z-axis (0,0,1). Angle is in radians, here
       it is ~pi/4 */
    cs_real_t rot_center[3] = {0.0, 0.0, 0.0};
    cs_real_t rot_axis[3] = {0.0, 0.0, 1.0};
    cs_real_t rot_angle = 0.7853981;

    cs_medcoupling_remapper_rotate(r, rot_center, rot_axis, rot_angle);

    /* Update of the interpolation matrix */
    cs_medcoupling_remapper_setup(r);
  }
  /*! [trans_rota] */

  /*! [copy_values] */

  /* We retrieve an array containing the interpolated values.
   * Inputs are:
   * 1) remapper
   * 2) id of the field to interpolate
   * 3) a default value (if no intersection is obtained) */
  cs_real_t *bc_scalar = cs_medcoupling_remapper_copy_values(r, 0, -5.0);

  /* We impose a dirichlet condition on all the faces of the boundary condition
   * related to the zone "inlet" */

  /*! [copy_values] */

  /*! [dirichlet_condition] */
  const cs_zone_t *z = cs_boundary_zone_by_name("inlet");

  cs_field_t *scalar = cs_field_by_name_try("scalar1");
  int       *icodcl  = scalar->bc_coeffs->icodcl;
  cs_real_t *rcodcl1 = scalar->bc_coeffs->rcodcl1;

  for (cs_lnum_t ielt = 0; ielt < z->n_elts; ielt++) {
    cs_lnum_t f_id = z->elt_ids[ielt];

    icodcl[f_id] = 1;
    rcodcl1[f_id] = bc_scalar[ielt];
  }
  /*! [dirichlet_condition] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
