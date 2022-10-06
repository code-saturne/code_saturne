/*============================================================================
 * User definition of physical properties.
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  const int n_fields = cs_field_n_fields();
  const int keysca = cs_field_key_id("scalar_id");
  const int kivisl = cs_field_key_id("diffusivity_id");

  /* Index field for the capacity table (C = grad(theta)/grad(h)) */
  cs_real_t *capacity = cs_field_by_name("capacity")->val;

  /* Field for the saturation */
  cs_real_t *saturation = cs_field_by_name("saturation")->val;

  /* Darcian velocity */
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  /* Field for permeability for the flow part */
  cs_real_t *permeability = NULL;
  cs_real_6_t *tensor_permeability = NULL;
  cs_field_t *f_perm = cs_field_by_name("permeability");

  if (f_perm->dim == 1)
    permeability = f_perm->val;
  else if (f_perm->dim == 6)
    tensor_permeability = (cs_real_6_t *)f_perm->val;

  /* Example physical properties of two fully saturated soils
   * ======================================================== */

  /* Flow part
   * ========= */

  /*![richards_flow_soils]*/

  /* Set parameters for soil 1 */
  const cs_zone_t *zn = cs_volume_zone_by_name_try("SOIL1");

  /* Loop on cell of the list */
  for (cs_lnum_t ii = 0; ii < zn->n_elts; ii++) {

    /* Cell number */
    const cs_lnum_t c_id = zn->elt_ids[ii];

    /* Set saturation and capacity (0 if saturated) */
    saturation[c_id] = 0.6;
    capacity[c_id] = 0.;

    /* Set permeability */
    if (f_perm->dim == 1) {
       permeability[c_id] = 1.e-1;
    }
    else if (f_perm->dim == 6) {
      tensor_permeability[c_id][0] = 1.e-1;
      tensor_permeability[c_id][1] = 1.e-1;
      tensor_permeability[c_id][2] = 1.e-2;
      tensor_permeability[c_id][3] = 0.0;
      tensor_permeability[c_id][4] = 0.0;
      tensor_permeability[c_id][5] = 0.0;
    }
  }

  /* Set parameters for soil 2 */
  zn = cs_volume_zone_by_name_try("SOIL2");

  /* Loop on cell of the list */
  for (cs_lnum_t ii = 0; ii < zn->n_elts; ii++) {

    /* Cell number */
    const cs_lnum_t c_id = zn->elt_ids[ii];

    /* Set saturation and capacity (0 if saturated) */
    saturation[c_id] = 0.4;
    capacity[c_id] = 0.;

    /* Set permeability */
    if (f_perm->dim == 1) {
       permeability[c_id] = 5.e-1;
    }
    else if (f_perm->dim == 6) {
      tensor_permeability[c_id][0] = 5.e-1;
      tensor_permeability[c_id][1] = 5.e-1;
      tensor_permeability[c_id][2] = 5.e-2;
      tensor_permeability[c_id][3] = 0.0;
      tensor_permeability[c_id][4] = 0.0;
      tensor_permeability[c_id][5] = 0.0;
    }
  }
  /*![richards_flow_soils]*/

  /*![richards_flow_solut]*/

  /* Transport part
   * ============== */

  /* Set parameters for soil 1 */
  const cs_zone_t *z1 = cs_volume_zone_by_name_try("SOIL1");

  /* Set parameters for soil 2 */
  const cs_zone_t *z2 = cs_volume_zone_by_name_try("SOIL2");

  for (int f_id = 0; f_id < n_fields; f_id++) {

    cs_field_t *fld = cs_field_by_id(f_id);

    /* Here we only handle user or model scalar-type variables */

    int sc_id = -1;
    if (fld->type & CS_FIELD_VARIABLE)
      sc_id = cs_field_get_key_int(fld, keysca) - 1;
    if (sc_id < 0)
      continue;

    const int ifcvsl = cs_field_get_key_int(fld, kivisl);
    if (ifcvsl < 0)
      continue;

    cs_real_t *cpro_vscalt = cs_field_by_id(ifcvsl)->val;

    /* Definition of the isotropic diffusion (dispersion and moleculer
       diffusion) */

    cs_real_t molecular_diffusion = 1.e-6;
    cs_real_t darcy_isotropic_dispersion = 1.;

    for (cs_lnum_t ii = 0; ii < z1->n_elts; ii++) {

      const cs_lnum_t c_id = z1->elt_ids[ii];
      const cs_real_t velocity_norm = cs_math_3_norm(vel[c_id]);

      cpro_vscalt[c_id] =   darcy_isotropic_dispersion*velocity_norm
                          + saturation[c_id]*molecular_diffusion;
    }

    molecular_diffusion = 1.e-8;
    darcy_isotropic_dispersion = 0.2;

    for (cs_lnum_t ii = 0; ii < z2->n_elts; ii++) {

      const cs_lnum_t c_id = z2->elt_ids[ii];
      const cs_real_t velocity_norm = cs_math_3_norm(vel[c_id]);

      cpro_vscalt[c_id] =   darcy_isotropic_dispersion*velocity_norm
                          + saturation[c_id]*molecular_diffusion;
    }
  }
  /*![richards_flow_solut]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
