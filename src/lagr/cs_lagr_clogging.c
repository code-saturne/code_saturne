/*============================================================================
 * Methods for particle clogging modeling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Functions dealing with the particle clogging modeling
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_search.h"
#include "cs_lagr_utils.h"
#include "cs_halo.h"
#include "cs_lagr_dlvo.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_clogging.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure declaration
 *============================================================================*/

static cs_lagr_clogging_param_t cs_lagr_clogging_param;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double _pi = 3.14159265358979323846;

/*============================================================================
 * Public function for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clogging initialization.
 *
 * - Retrieve various parameters for storing in global structure.
 * - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
CS_PROCF (cloginit, CLOGINIT)(const cs_real_t   *faraday_cst,
                              const cs_real_t   *free_space_permit,
                              const cs_real_t   *water_permit,
                              const cs_real_t   *ionic_strength,
                              const cs_real_t   *jamming_limit,
                              const cs_real_t   *min_porosity,
                              const cs_real_t    temperature[],
                              const cs_real_t   *phi1,
                              const cs_real_t   *phi2,
                              const cs_real_t   *cstham,
                              const cs_real_t   *dcutof,
                              const cs_real_t   *lambwl,
                              const cs_real_t   *kboltz
 )
{
#define PG_CST 8.314  /* Ideal gas constant */

  cs_lnum_t ifac;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clogging_param          */

  cs_lagr_clogging_param.faraday_cst = *faraday_cst;
  cs_lagr_clogging_param.free_space_permit = *free_space_permit;
  cs_lagr_clogging_param.water_permit = *water_permit;
  cs_lagr_clogging_param.ionic_strength = *ionic_strength;
  cs_lagr_clogging_param.jamming_limit = *jamming_limit;
  cs_lagr_clogging_param.min_porosity = *min_porosity;
  cs_lagr_clogging_param.phi1 = *phi1;
  cs_lagr_clogging_param.phi2 = *phi2;
  cs_lagr_clogging_param.cstham = *cstham;
  cs_lagr_clogging_param.dcutof = *dcutof;
  cs_lagr_clogging_param.lambwl = *lambwl;
  cs_lagr_clogging_param.kboltz = *kboltz;

  /* Allocate memory for the temperature and Debye length arrays */

  if (cs_lagr_clogging_param.temperature == NULL)
    BFT_MALLOC(cs_lagr_clogging_param.temperature, mesh->n_b_faces, cs_real_t);

  if (cs_lagr_clogging_param.debye_length == NULL)
    BFT_MALLOC(cs_lagr_clogging_param.debye_length, mesh->n_b_faces, cs_real_t);

  /* Store the temperature */

  for (ifac = 0; ifac < mesh->n_b_faces; ifac++)
    cs_lagr_clogging_param.temperature[ifac] = temperature[ifac];

  /* Computation and storage of the Debye length                */

  for (ifac = 0; ifac < mesh->n_b_faces; ifac++)

    cs_lagr_clogging_param.debye_length[ifac]
      = pow(2e3 * pow(cs_lagr_clogging_param.faraday_cst,2)
            * cs_lagr_clogging_param.ionic_strength
            /  (  cs_lagr_clogging_param.water_permit
                * cs_lagr_clogging_param.free_space_permit * PG_CST
                * cs_lagr_clogging_param.temperature[ifac]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" cstfar = %g\n", cs_lagr_clogging_param.faraday_cst);
  bft_printf(" epsvid = %g\n", cs_lagr_clogging_param.free_space_permit);
  bft_printf(" epseau = %g\n", cs_lagr_clogging_param.water_permit);
  bft_printf(" fion   = %g\n", cs_lagr_clogging_param.ionic_strength);
  bft_printf(" temp[1]   = %g\n", cs_lagr_clogging_param.temperature[0]);
  bft_printf(" debye[1]   = %g\n", cs_lagr_clogging_param.debye_length[0]);
  bft_printf(" phi1   = %g\n", cs_lagr_clogging_param.phi1);
  bft_printf(" phi2  = %g\n", cs_lagr_clogging_param.phi2);
#endif
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clogging finalization.
 *
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_clogging_finalize(void)
{
  BFT_FREE(cs_lagr_clogging_param.temperature);
  BFT_FREE(cs_lagr_clogging_param.debye_length);
}

/*----------------------------------------------------------------------------
 * Clogging:
 *
 * - Compute the number of deposited particles in contact with the depositing
 *   particle
 * - Re-compute the energy barrier if this number is greater than zero
 *
 * parameters:
 *   particle         <-- pointer to particle data
 *   attr_map         <-- pointer to attribute map
 *   face_id          <-- id of face neighboring particle
 *   face_area        <-- area of face
 *   energy_barrier   <-> energy barrier
 *   surface_coverage <-> surface coverage
 *   limit            <-> jamming limit
 *   mporos           <-> minimum porosity
 *
 * returns:
 *   number of deposited particles in contact with the depositing particle
 *----------------------------------------------------------------------------*/

int
cs_lagr_clogging_barrier(const void                     *particle,
                         const cs_lagr_attribute_map_t  *attr_map,
                         cs_lnum_t                       face_id,
                         cs_real_t                       face_area,
                         cs_real_t                      *energy_barrier,
                         cs_real_t                      *surface_coverage,
                         cs_real_t                      *limit,
                         cs_real_t                      *mporos)
{
  cs_real_t contact_area;
  cs_real_t deposited_radius;

  cs_real_t mean_nb_cont;

  cs_lnum_t  dim_aux = 1, contact_count[1], nbtemp[12000];
  cs_lnum_t  param1;
  cs_real_t  param2, value;
  cs_lnum_t  k,i;

  /* Computation of the number of particles in contact with */
  /* the depositing particle */

  /* Assuming monodispersed calculation */

  double p_stat_weight
    = cs_lagr_particle_get_real(particle, attr_map, CS_LAGR_STAT_WEIGHT);
  double p_diameter
    = cs_lagr_particle_get_real(particle, attr_map, CS_LAGR_STAT_WEIGHT);
  cs_real_t depositing_radius = p_diameter * 0.5;

  deposited_radius = depositing_radius;

  contact_area = p_stat_weight
                 * _pi * pow(2. * pow(deposited_radius * depositing_radius, 0.5)
                             + deposited_radius,2);

  mean_nb_cont =   contact_area
                 * (*surface_coverage) / (_pi * pow(deposited_radius,2));

  /* Assuming Poisson distribution */

  CS_PROCF(fische, FISCHE)(&dim_aux, &mean_nb_cont, contact_count);

  value = 700.;
  if (mean_nb_cont > value) {
    param1 = mean_nb_cont / value;
    param2 = fmod(mean_nb_cont,value);
    assert(param1 < 12000); /* TODO use dynamic allocation or set bound */

    CS_PROCF(fische, FISCHE)(&dim_aux, &param2, contact_count);
    CS_PROCF(fische, FISCHE)(&param1, &value, nbtemp);

    for (k = 0; k < param1; k++) {
      contact_count[0] = contact_count[0] + nbtemp[k];
    }
  }

  /* If the surface coverage is above the jamming limit,
     we are in multilayer deposition, so the contact number
     must be greater than zero  */

   /* The surface coverage must be greater than zero */
  if (*surface_coverage > 1e-15)

    if (((_pi * pow(depositing_radius,2) * p_stat_weight)/face_area
         + (*surface_coverage)) > cs_lagr_clogging_param.jamming_limit)
      contact_count[0] +=1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (mean_nb_cont > 0) {
    bft_printf("mean number = %g\n", mean_nb_cont);
    bft_printf("calculated number = %d\n",  contact_count[0]);
  }
#endif

  if (contact_count[0] != 0) {

    *energy_barrier = 0;

    /* Computation of the energy barrier */
    for (i = 0; i < 101; i++) {

      cs_real_t  step = 1e-10;

      cs_real_t distcc =   cs_lagr_clogging_param.dcutof + i*step
                         + depositing_radius + deposited_radius;

      cs_real_t var1
        = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                              deposited_radius,
                                              depositing_radius,
                                              cs_lagr_clogging_param.lambwl,
                                              cs_lagr_clogging_param.cstham);

      cs_real_t var2
        = cs_lagr_edl_sphere_sphere(distcc,
                                    deposited_radius,
                                    depositing_radius,
                                    cs_lagr_clogging_param.phi1,
                                    cs_lagr_clogging_param.phi2,
                                    cs_lagr_clogging_param.kboltz,
                                    cs_lagr_clogging_param.temperature[face_id],
                                    cs_lagr_clogging_param.debye_length[face_id],
                                    cs_lagr_clogging_param.free_space_permit,
                                    cs_lagr_clogging_param.water_permit);

      cs_real_t var = contact_count[0] * (var1 + var2);

      if (var > *energy_barrier)
        *energy_barrier = var;
      if (var < 0)
        *energy_barrier = 0;

    }

    *energy_barrier =  *energy_barrier / (0.5 * p_diameter);
  }

  *limit = cs_lagr_clogging_param.jamming_limit;
  *mporos = cs_lagr_clogging_param.min_porosity;

  return contact_count[0];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
