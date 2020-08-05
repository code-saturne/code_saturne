/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_gui_util.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_random.h"

#include "cs_mesh.h"

#include "cs_lagr.h"
#include "cs_lagr_new.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Experimental measurements
   ------------------------- */

/* transverse coordinate */

static const cs_real_t
  _zi[] = {0.e-3, 1.e-3, 1.5e-3, 2.0e-3, 2.5e-3,
           3.0e-3, 3.5e-3, 4.0e-3, 4.5e-3, 5.0e-3};

/* particle volume fraction */

static const cs_real_t
  _lvf[] = {0.377e-4, 2.236e-4, 3.014e-4, 4.306e-4, 5.689e-4,
            8.567e-4, 7.099e-4, 4.520e-4, 2.184e-4, 0.377e-4};

/* vertical mean velocity of the particles */

static const cs_real_t
  _ui[] = {5.544, 8.827, 9.068, 9.169, 8.923,
           8.295, 7.151, 6.048, 4.785, 5.544};

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Computation of particle injection profile.
 *
 * Note: if the input pointer is non-NULL, it must point to valid data
 * when the selection function is called, so that value or structure should
 * not be temporary (i.e. local);
 *
 * parameters:
 *   zone_id     <-- id of associated mesh zone
 *   location_id <-- id of associated mesh location
 *   input       <-- pointer to optional (untyped) value or structure.
 *   n_elts      <-- number of zone elements
 *   elt_ids     <-- ids of zone elements
 *   profile     <-- weight of a given zone element (size: n_elts)
 *----------------------------------------------------------------------------*/

static void
_injection_profile(int               zone_id,
                   int               location_id,
                   const void       *input,
                   cs_lnum_t         n_elts,
                   const cs_lnum_t   elt_ids[],
                   cs_real_t         profile[])
{
  const cs_real_3_t  *b_face_coords
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  const int itmx = 8;

  /* Loop on elements
     ---------------- */

  for (cs_lnum_t ei = 0; ei < n_elts; ei++) {

    /* Face center */

    const cs_lnum_t face_id = elt_ids[ei];

    const cs_real_t z = b_face_coords[face_id][2];

    /* Interpolation */

    int i = 0;

    if (z > _zi[0]) {
      for (i = 0; i < itmx; i++) {
        if (z >= _zi[i] && z < _zi[i+1])
          break;
      }
    }

    /* Compute volume fraction and statistical weight */

    cs_real_t up = _ui[i] +(z-_zi[i])*(_ui[i+1]-_ui[i])/(_zi[i+1]-_zi[i]);
    cs_real_t lvfp = _lvf[i] + (z-_zi[i])*(_lvf[i+1]-_lvf[i])/(_zi[i+1]-_zi[i]);

    /* number of particles in the cell */

    profile[ei] = lvfp * up;

    assert(profile[ei] >= 0.);

  }
}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle boundary conditions.
 *
 * This is used for the definition of inlet and other boundaries,
 * based on predefined boundary zones (\ref cs_boundary_zone_t).
 *
 * \param[in] bc_type    type of the boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_boundary_conditions(const int  bc_type[])
{
  cs_mesh_t *mesh = cs_glob_mesh;

  cs_lagr_zone_data_t *lagr_bcs = cs_lagr_get_boundary_conditions();

  const cs_real_3_t *b_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;
  const cs_real_3_t *b_face_normal
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_normal;

  /* Assign injection profile function */
  {
    const cs_boundary_zone_t  *z = cs_boundary_zone_by_name("inlet_p");
    int set_id = 0;

    cs_lagr_injection_set_t *zis
      = cs_lagr_get_injection_set(lagr_bcs, z->id, set_id);

    zis->injection_profile_func = _injection_profile;

    /* Compute total flow rate based on known volume fractions
       and velocity in each section (total flow rate would be simpler
       if available) */

    cs_real_t flow_rate = 0;

    cs_real_t p_rho = zis->density;

    for (int iz = 0; iz < 9; iz++) {
      /* segment length and associated surface */
      cs_real_t zl = _zi[iz+1] - _zi[iz];
      cs_real_t zsurf = zl * 0.005; /* zl * mesh thickness; TODO automate */
      /* mean flow rate (volume fraction*velocity*density) */
      cs_real_t vf = (_lvf[iz+1]*_ui[iz] + _lvf[iz+1]*_ui[iz]) * 0.5*p_rho;
      /* contribution to total flow rate */
      flow_rate += vf*zsurf;
    }

    zis->flow_rate = flow_rate;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
