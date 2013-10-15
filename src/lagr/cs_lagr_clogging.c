/*============================================================================
 * Methods for particle clogging modeling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_clogging.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


/*============================================================================
 * Local structure declaration
 *============================================================================*/

static cs_lagr_clog_param_t cs_lagr_clog_param;


/*============================================================================
 * Public function for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Clogging initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
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

  int ifac;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clog_param          */

  cs_lagr_clog_param.faraday_cst = *faraday_cst;
  cs_lagr_clog_param.free_space_permit = *free_space_permit;
  cs_lagr_clog_param.water_permit = *water_permit;
  cs_lagr_clog_param.ionic_strength = *ionic_strength;
  cs_lagr_clog_param.jamming_limit = *jamming_limit;
  cs_lagr_clog_param.min_porosity = *min_porosity;
  cs_lagr_clog_param.phi1 = *phi1;
  cs_lagr_clog_param.phi2 = *phi2;
  cs_lagr_clog_param.cstham = *cstham;
  cs_lagr_clog_param.dcutof = *dcutof;
  cs_lagr_clog_param.lambwl = *lambwl;
  cs_lagr_clog_param.kboltz = *kboltz;

  /* Allocate memory for the temperature and Debye length arrays */

  BFT_MALLOC(cs_lagr_clog_param.temperature, mesh->n_b_faces, cs_real_t);
  BFT_MALLOC(cs_lagr_clog_param.debye_length, mesh->n_b_faces, cs_real_t);

  /* Store the temperature */

  for (ifac = 0; ifac < mesh->n_b_faces ; ifac++)
    cs_lagr_clog_param.temperature[ifac] = temperature[ifac];

  /* Computation and storage of the Debye length                */

  for (ifac = 0; ifac < mesh->n_b_faces ; ifac++)

    cs_lagr_clog_param.debye_length[ifac] =
      pow(2e3 * pow(cs_lagr_clog_param.faraday_cst,2) * cs_lagr_clog_param.ionic_strength /
      (cs_lagr_clog_param.water_permit * cs_lagr_clog_param.free_space_permit * PG_CST *
      cs_lagr_clog_param.temperature[ifac]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" cstfar = %g\n", cs_lagr_clog_param.faraday_cst);
  bft_printf(" epsvid = %g\n", cs_lagr_clog_param.free_space_permit);
  bft_printf(" epseau = %g\n", cs_lagr_clog_param.water_permit);
  bft_printf(" fion   = %g\n", cs_lagr_clog_param.ionic_strength);
  bft_printf(" temp[1]   = %g\n", cs_lagr_clog_param.temperature[0]);
  bft_printf(" debye[1]   = %g\n", cs_lagr_clog_param.debye_length[0]);
  bft_printf(" phi1   = %g\n", cs_lagr_clog_param.phi1);
  bft_printf(" phi2  = %g\n", cs_lagr_clog_param.phi2);
#endif



}

/*----------------------------------------------------------------------------
 * Clogging ending
 * Deallocate the arrays storing temperature and Debye length
 *----------------------------------------------------------------------------*/


void clogend()
{
  BFT_FREE(cs_lagr_clog_param.temperature);
  BFT_FREE(cs_lagr_clog_param.debye_length);
}


/*----------------------------------------------------------------------------
 * Clogging:
 *
 * - Calculation of the number of deposited particles in contact with the
 *   depositing particle
 * - Re-calculation of the energy barrier in this number is greater than zero
 *----------------------------------------------------------------------------*/

cs_int_t
clogging_barrier(cs_lagr_particle_t     particle,
                 cs_int_t               face_id,
                 cs_real_t              face_area,
                 cs_real_t              *energy_barrier,
                 cs_real_t              *surface_coverage,
                 cs_real_t              *limit,
                 cs_real_t              *mporos

  )
{

#define PI 3.141592653589793

  cs_real_t contact_area;
  cs_real_t depositing_radius = particle.diameter * 0.5;
  cs_real_t deposited_radius;

  cs_real_t mean_nb_cont;

  cs_int_t  dim_aux = 1, contact_number[1], nbtemp[12000];
  cs_real_t param2, value ;
  cs_int_t  param1;
  cs_lnum_t  k,i;

/* Computation of the number of particles in contact with */
/* the depositing particle */

/* Assuming monodisperse calculation */

  deposited_radius = depositing_radius;

  contact_area = particle.stat_weight *
                 PI * pow(2. * pow(deposited_radius * depositing_radius, 0.5) + deposited_radius,2);

  mean_nb_cont = contact_area * (*surface_coverage) / (PI * pow(deposited_radius,2));

/* Assuming Poisson distribution */

  CS_PROCF(fische, FISCHE)(&dim_aux, &mean_nb_cont, contact_number);

  value = 700.;
  if (mean_nb_cont > value) {
    param1 = mean_nb_cont / value;
    param2 = fmod(mean_nb_cont,value);

    CS_PROCF(fische, FISCHE)(&dim_aux, &param2, contact_number);
    CS_PROCF(fische, FISCHE)(&param1, &value ,nbtemp);

    for (k = 0; k < param1; k++) {
      contact_number[0] =  contact_number[0] + nbtemp[k];
    }
  }

/* If the surface coverage is above the jamming limit,
   we are in multilayer deposition, so the contact number
   must be greater than zero  */

  if (*surface_coverage > 1e-15)   /* The surface coverage must be greater than zero */

    if (((PI * pow(depositing_radius,2) *  particle.stat_weight)/face_area + (*surface_coverage)) > cs_lagr_clog_param.jamming_limit)

      contact_number[0] +=1;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (mean_nb_cont > 0) {
    bft_printf("mean number = %g\n", mean_nb_cont);
    bft_printf("calculated number = %d\n",  contact_number[0]);
  }
#endif


  if (contact_number[0] != 0) {

    /* Computation of the energy barrier */
    for (i = 0; i < 101; i++) {

      cs_real_t  step = 1e-10;

      cs_real_t distcc = cs_lagr_clog_param.dcutof + i * step + depositing_radius + deposited_radius;

      cs_real_t var1 = vdwsa(distcc,deposited_radius,depositing_radius);

      cs_real_t var2 = edlsa(face_id,distcc,deposited_radius,depositing_radius);

      cs_real_t var = contact_number[0] * (var1 + var2);

      if (var > 0.) {
        *energy_barrier = var / (0.5 * particle.diameter);}
      else {
        *energy_barrier = 0;}

    }
  }

  *limit = cs_lagr_clog_param.jamming_limit;
  *mporos = cs_lagr_clog_param.min_porosity;

  return contact_number[0];
}

/*----------------------------------------------------------------------------
 *    Calculation of the Van der Waals interaction between two spheres
 *    following the formula from Gregory (1981a)
 *----------------------------------------------------------------------------*/

cs_real_t
vdwsa(           cs_real_t              distcc,
                 cs_real_t              rpart1,
                 cs_real_t              rpart2
                 )
{
  cs_real_t var = - cs_lagr_clog_param.cstham * rpart1 * rpart2 / (6 * (distcc - rpart1 - rpart2)
                  * (rpart1 + rpart2)) * (1 - 5.32 * (distcc - rpart1 - rpart2)
                  / cs_lagr_clog_param.lambwl * log(1 +  cs_lagr_clog_param.lambwl / (distcc - rpart1 - rpart2) / 5.32));

  return var;
}

/*----------------------------------------------------------------------------
 *     Calculation of the EDL interaction between two spheres
 *     following the formula from Bell & al (1970)
 *     based on the McCartney & Levine method
 *----------------------------------------------------------------------------*/

cs_real_t
edlsa(            cs_int_t               ifac,
                  cs_real_t              distcc,
                  cs_real_t              rpart1,
                  cs_real_t              rpart2
                 )
{
#define PI 3.141592653589793

  cs_real_t charge = 1.6e-19;

 /* Reduced zeta potential */
  cs_real_t lphi1 =  charge * cs_lagr_clog_param.phi1 /  cs_lagr_clog_param.kboltz / cs_lagr_clog_param.temperature[ifac];
  cs_real_t lphi2 =  charge * cs_lagr_clog_param.phi1 /  cs_lagr_clog_param.kboltz / cs_lagr_clog_param.temperature[ifac];

  cs_real_t tau = rpart1 / (1. / cs_lagr_clog_param.debye_length[ifac]);

  /* Extended reduced zeta potential */
  /* (following the work from Ohshima et al, 1982, JCIS, 90, 17-26) */

  lphi1 = 8. * tanh(lphi1 / 4.) /
         ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
         * pow(tanh(lphi1 / 4.),2),0.5));

  lphi2 = 8. * tanh(lphi2 / 4.) /
          ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
          * pow(tanh(lphi2 / 4.),2),0.5));

  cs_real_t alpha = sqrt(rpart2 * (distcc - rpart2) / (rpart1 * (distcc - rpart1)))
                   + sqrt(rpart1 * (distcc - rpart1) / (rpart2 * (distcc - rpart2)));

  cs_real_t omega1 = pow(lphi1,2) + pow(lphi2,2) + alpha * lphi1 * lphi2;

  cs_real_t omega2 = pow(lphi1,2) + pow(lphi2,2) - alpha * lphi1 * lphi2;

  cs_real_t gamma = sqrt(rpart1 * rpart2 / (distcc - rpart1) / (distcc - rpart2))
                    *exp(1 / cs_lagr_clog_param.debye_length[ifac] * (rpart1 + rpart2 - distcc));

  cs_real_t var = 2 * PI * cs_lagr_clog_param.free_space_permit * cs_lagr_clog_param.water_permit
                   * pow((cs_lagr_clog_param.kboltz * cs_lagr_clog_param.temperature[ifac] / charge),2)
                   * rpart1 * rpart2 * (distcc - rpart1) * (distcc - rpart2)
                   / (distcc * ( distcc * ( rpart1  + rpart2) - pow(rpart1,2) - pow(rpart2,2)))
                   * (omega1 * log(1 + gamma) + omega2 * log(1 - gamma));

  return var;

}
/*----------------------------------------------------------------------------*/

/* Delete local macro definitions */

END_C_DECLS
