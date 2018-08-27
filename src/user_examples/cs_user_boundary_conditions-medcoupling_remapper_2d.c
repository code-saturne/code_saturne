/*============================================================================
 * User definition of boundary conditions.
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

  /* TEST TO REMOVE */
  if (true) {
    return;
  }

#if defined(HAVE_MEDCOUPLING_LOADER)

  /* Variables needed for boundary condition sub-selection */
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  cs_lnum_t *lstelt = NULL;
  cs_lnum_t  nelts;

  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);


  cs_field_t *f;
  /* -----------------------------------------------------------
   * MEDCoupling Remapper structure:
   * ----------------------------------------------------------- */

  /* Number of fields to interpolate from the medfile */
  const int    nremapper_fields = 1;

  /* Names of the fields to read */
  const char **field_names = NULL;
  BFT_MALLOC(field_names, nremapper_fields, const char *);
  field_names[0] = "TEMPERATURE";

  /* Indexes needed to read the time step from the
   * file (-1, -1 if only one exists) */
  int it0 = -1;
  int it1 = -1;

  /* We ask for a remapper with a certain name. If it does not exist,
   * the function returns a NULL pointer. */
  cs_medcoupling_remapper_t *r = cs_medcoupling_remapper_by_name_try("scalar_bc");

  /* If the returned pointer is NULL (first call), we create the
   * corresponding remapper */
  if (r == NULL) {

    /* Space dimension of the elements (2 for faces, 3 for cells) */
    int elts_dim = 2;

    /* The remapper is created. We retrieve its id from the function. The function
     * inputs are:
     * 1) Name of the remapper
     * 2) dimension of the mesh elements
     * 3) selection criteria for the boundary condition zone
     * 4) path to the med file
     * 5) number of fields to interpolate
     * 6) names of the fields to interpolate
     * 7 + 8) time iteration index and order
     */
    int r_id = cs_medcoupling_remapper_initialize("scalar_bc",
                                                  elts_dim,
                                                  "inlet",
                                                  "/home/i76777/Etudes/PARAMEDMEM/BC_TEST/carte2D_Tfluid.med",
                                                  nremapper_fields,
                                                  field_names,
                                                  it0,
                                                  it1);

    /* We retrieve the pointer */
    r = cs_medcoupling_remapper_by_id(r_id);

    /* We create the interpolation matrix => Here it is only called once
     * since the mesh is not moving */
    cs_medcoupling_remapper_setup(r);

  }

  /* If the med data needs for a translation or rotation for the geometrical
   * superposition with the target Code_Saturne mesh:
   */
  if (false) {
    // Translation using a tranlsation vector. Here it is (1, 0, 0)
    cs_real_t translation_vector[3];
    translation_vector[0] = 1.0;
    translation_vector[1] = 0.0;
    translation_vector[2] = 0.0;
    cs_medcoupling_remapper_translate(r, translation_vector);

    // Rotation using an invariant point, the rotation axis and rotation angle
    // Here, center is O=(0,0,0) and z-axis (0,0,1). Angle is in radians, here
    // it is ~pi/4

    cs_real_t rot_center[3];
    rot_center[0] = 0.0;
    rot_center[1] = 0.0;
    rot_center[2] = 0.0;

    cs_real_t rot_axis[3];
    rot_axis[0] = 0.0;
    rot_axis[1] = 0.0;
    rot_axis[2] = 1.0;

    cs_real_t rot_angle = 0.7853981;

    cs_medcoupling_remapper_rotate(r, rot_center, rot_axis, rot_angle);

    // Update of the interpolation matrix
    cs_medcoupling_remapper_setup(r);
  }



  /* We retrieve an array containing the interpolated values.
   * Inputs are:
   * 1) remapper
   * 2) id of the field to interpolate
   * 3) a default value (if no intersection is obtained )
   */
  cs_real_t *bc_scalar = cs_medcoupling_remapper_copy_values(r, 0, -5.0);

  /* We impose a dirichlet condition on all the faces of the boundary condition
   * related to the zone "inlet" */
  const int keyvar = cs_field_key_id("variable_id");

  cs_field_t *scalar = cs_field_by_name_try("scalar1");
  int iscal = cs_field_get_key_int(scalar, keyvar) - 1;

  cs_selector_get_b_face_list("inlet", &nelts, lstelt);

  for (cs_lnum_t ielt = 0; ielt < nelts; ielt++) {
    cs_lnum_t f_id = lstelt[ielt];

    icodcl[iscal*n_b_faces + f_id] = 1;
    rcodcl[iscal*n_b_faces + f_id] = bc_scalar[ielt];
  }

  BFT_FREE(lstelt);

#endif

  return;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
