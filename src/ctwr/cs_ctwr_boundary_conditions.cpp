/*============================================================================
 * Cooling towers functions for boundary conditions
 *============================================================================*/

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_atmo.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_boundary_conditions.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


/*----------------------------------------------------------------------------
 * Automatic boundary condition for cooling towers
 *----------------------------------------------------------------------------*/

void
cs_ctwr_bcond(void)
{
  /* Mesh-related data */
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const int *bc_type = cs_glob_bc_type;

  /* Fluid properties and physical variables */
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;
  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t *vel_rcodcl1 = CS_F_(vel)->bc_coeffs->rcodcl1;
  cs_field_t *y_l_r= cs_field_by_name("ym_l_r");
  cs_field_t *yh_l_r= cs_field_by_name_try("ymh_l_r");
  cs_field_t *yh_l_p= cs_field_by_name("yh_l_packing");
  cs_field_t *y_l_p = cs_field_by_name("y_l_packing");
  cs_field_t *ym_w = cs_field_by_name("ym_water");
  cs_field_t *t_h = cs_field_by_name("temperature"); /* Humid air temp */
  cs_real_t tkelvin = cs_physical_constants_celsius_to_kelvin;

  const cs_real_t xhum = air_prop->humidity0;
  cs_real_t ref_temp = phys_pro->t0;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    cs_atmo_option_t *aopt = cs_glob_atmo_option;

    /* Express Dirichlet directly in term of theta and not in real
     * temperature */
    cs_real_t pref = cs_glob_atmo_constants->ps;
    cs_real_t rair = phys_pro->r_pg_cnst;
    cs_real_t cp0 = phys_pro->cp0;
    cs_real_t rscp = rair/cp0;
    cs_real_t clatev = phys_pro->clatev;

    /* Ref temperature is potential temperature */
    ref_temp = (aopt->meteo_t0 - clatev/cp0 * aopt->meteo_ql0)
                   * pow(pref/ aopt->meteo_psea, rscp);
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {

    if (bc_type[face_id] == CS_INLET || bc_type[face_id] == CS_FREE_INLET
        || bc_type[face_id] == CS_CONVECTIVE_INLET) {

      /* The turbulence BC values are calculated upstream using the base
         mechanism, so nothing specific is needed here. */

      /* Boundary conditions for the transported temperature of the humid air and
       * of the liquid water injected in the packing zones
       * --> Bulk values if not specified by the user
       * Assuming humid air is at conditions '0' */

      /* For humid air temperature */
      if (t_h->bc_coeffs->rcodcl1[face_id] > 0.5 * cs_math_infinite_r)
        t_h->bc_coeffs->rcodcl1[face_id] = ref_temp;

      /* For water mass fraction */
      if (ym_w->bc_coeffs->rcodcl1[face_id] > 0.5 * cs_math_infinite_r)
        ym_w->bc_coeffs->rcodcl1[face_id] = xhum / (1 + xhum);

      /* For injected liquid in the packing*/
      if (y_l_p->bc_coeffs->rcodcl1[face_id] > 0.5 * cs_math_infinite_r)
        y_l_p->bc_coeffs->rcodcl1[face_id] = 0.;

      /* For injected liquid enthalpy in the packing*/
      if (yh_l_p->bc_coeffs->rcodcl1[face_id] > 0.5 * cs_math_infinite_r) {
        cs_real_t t_l = phys_pro->t0 - tkelvin;
        cs_real_t b_h_l = cs_liq_t_to_h(t_l);
        /* Y_l . h_l is transported (not only h_l) */
        cs_real_t b_yh_l = b_h_l * y_l_p->bc_coeffs->rcodcl1[face_id];
        yh_l_p->bc_coeffs->rcodcl1[face_id] = b_yh_l;
      }
    }

    /* For walls -> 0 flux for previous variables
     * Dirichlet condition y_l_r = 0 to mimic water basin drain and avoid rain
     * accumulation on the floor */

    else if (   bc_type[face_id] == CS_SMOOTHWALL
             || bc_type[face_id] == CS_ROUGHWALL) {

      t_h->bc_coeffs->icodcl[face_id] = 3;
      t_h->bc_coeffs->rcodcl3[face_id] = 0.;
      ym_w->bc_coeffs->icodcl[face_id] = 3;
      ym_w->bc_coeffs->rcodcl3[face_id] = 0.;

      yh_l_p->bc_coeffs->icodcl[face_id] = 3;
      yh_l_p->bc_coeffs->rcodcl3[face_id] = 0.;
      y_l_p->bc_coeffs->icodcl[face_id] = 3;
      y_l_p->bc_coeffs->rcodcl3[face_id] = 0.;

      y_l_r->bc_coeffs->icodcl[face_id] = 1;
      y_l_r->bc_coeffs->rcodcl1[face_id] = 0.;
      if (yh_l_r != NULL) {
        yh_l_r->bc_coeffs->icodcl[face_id] = 1;
        yh_l_r->bc_coeffs->rcodcl1[face_id] = 0.;
      }

    }
  }

  /* Extra variables to load if we solve rain velocity */

  const cs_ctwr_option_t *ct_opt = cs_glob_ctwr_option;

  if (ct_opt->solve_rain_velocity) {
    char f_name[80];
    int class_id = 1;

    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *vp = cs_field_by_name(f_name);

    for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
      for (cs_lnum_t i = 0; i < 3; i++ ) {
        if (   bc_type[face_id] == CS_INLET
            || bc_type[face_id] == CS_FREE_INLET) {
          vp->bc_coeffs->icodcl[face_id] = 1;
          vp->bc_coeffs->rcodcl1[n_b_faces*i + face_id]
            = vel_rcodcl1[n_b_faces * i + face_id];
        }

        else if (   bc_type[face_id] == CS_SMOOTHWALL
                 || bc_type[face_id] == CS_ROUGHWALL) {
          vp->bc_coeffs->icodcl[face_id] = 1;
          vp->bc_coeffs->rcodcl1[n_b_faces*i + face_id] = 0.;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
