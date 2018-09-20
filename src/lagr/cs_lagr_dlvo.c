/*============================================================================
 * Calculation of DLVO forces
 *============================================================================*/

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

/*============================================================================
 * Functions dealing with the DLVO forces
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
#include "cs_halo.h"

#include "cs_lagr.h"
#include "cs_lagr_roughness.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_dlvo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local macro declarations
 *============================================================================*/

#define PG_CST 8.314  /* Ideal gas constant */

/*============================================================================
 * Local structure declarations
 *============================================================================*/

static cs_lagr_dlvo_param_t cs_lagr_dlvo_param;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double _pi = 3.14159265358979323846;

/* Cut-off distance for adhesion forces (assumed to be the Born distance) */
static const cs_real_t  _d_cut_off = 1.65e-10;

/* Boltzmann constant */
static const double _k_boltzmann = 1.38e-23;

/* Free space permittivity */
static const cs_real_t _free_space_permit = 8.854e-12;

/* Faraday constant */
static const cs_real_t _faraday_cst = 9.648e4;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * DLVO initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
cs_lagr_dlvo_init(const cs_real_t   water_permit,
                  const cs_real_t   ionic_strength,
                  const cs_real_t   temperature[],
                  const cs_real_t   valen,
                  const cs_real_t   phi_p,
                  const cs_real_t   phi_s,
                  const cs_real_t   cstham,
                  const cs_real_t   csthpp,
                  const cs_real_t   lambda_vdw)
{
  cs_lnum_t iel;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clog_param          */

  cs_lagr_dlvo_param.water_permit = water_permit;
  cs_lagr_dlvo_param.ionic_strength = ionic_strength;
  cs_lagr_dlvo_param.valen = valen;
  cs_lagr_dlvo_param.phi_p = phi_p;
  cs_lagr_dlvo_param.phi_s = phi_s;
  cs_lagr_dlvo_param.cstham = cstham;
  cs_lagr_dlvo_param.cstham = csthpp;
  cs_lagr_dlvo_param.lambda_vdw = lambda_vdw;

  /* Allocate memory for the temperature and Debye length arrays */

  if (cs_lagr_dlvo_param.temperature == NULL)
    BFT_MALLOC(cs_lagr_dlvo_param.temperature, mesh->n_cells, cs_real_t);

  if (cs_lagr_dlvo_param.debye_length == NULL)
    BFT_MALLOC(cs_lagr_dlvo_param.debye_length, mesh->n_cells, cs_real_t);

  /* Store the temperature */

  for (iel = 0; iel < mesh->n_cells; iel++)
    cs_lagr_dlvo_param.temperature[iel] = temperature[iel];

  /* Computation and storage of the Debye length */

  for (iel = 0; iel < mesh->n_cells ; iel++)

    cs_lagr_dlvo_param.debye_length[iel]
      =   pow(2e3 * pow(_faraday_cst,2)
        * cs_lagr_dlvo_param.ionic_strength /
        (cs_lagr_dlvo_param.water_permit
         * _free_space_permit * PG_CST
         * cs_lagr_dlvo_param.temperature[iel]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" epseau = %g\n", cs_lagr_dlvo_param.water_permit);
  bft_printf(" fion   = %g\n", cs_lagr_dlvo_param.ionic_strength);
  bft_printf(" temp[0]   = %g\n", cs_lagr_dlvo_param.temperature[0]);
  bft_printf(" valen   = %g\n", cs_lagr_dlvo_param.valen);
  bft_printf(" debye[0]   = %g\n", cs_lagr_dlvo_param.debye_length[0]);
  bft_printf(" phi_p   = %g\n", cs_lagr_dlvo_param.phi_p);
  bft_printf(" phi_s  = %g\n", cs_lagr_dlvo_param.phi_s);
#endif

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_dlvo_finalize()
{
  BFT_FREE(cs_lagr_dlvo_param.temperature);
  BFT_FREE(cs_lagr_dlvo_param.debye_length);
}

/*----------------------------------------------------------------------------
 * Compute the energy barrier for a smooth wall.
 *----------------------------------------------------------------------------*/

void
cs_lagr_barrier(const void                     *particle,
                const cs_lagr_attribute_map_t  *attr_map,
                cs_lnum_t                       iel,
                cs_real_t                      *energy_barrier)
{
  cs_lnum_t i;
  cs_real_t rpart = cs_lagr_particle_get_real(particle, attr_map,
                                              CS_LAGR_DIAMETER) * 0.5;

  *energy_barrier = 0.;

  cs_real_t barr = 0.;
  cs_real_t distp = 0.;

  /* Computation of the energy barrier */

  for (i = 0; i < 1001; i++) {

    cs_real_t  step = cs_lagr_dlvo_param.debye_length[iel]/30.0;

    /* Interaction between the sphere and the plate */

    distp = _d_cut_off + i * step;

    cs_real_t var1
      = cs_lagr_van_der_waals_sphere_plane(distp,
                                           rpart,
                                           cs_lagr_dlvo_param.lambda_vdw,
                                           cs_lagr_dlvo_param.cstham);

    cs_real_t var2
      = cs_lagr_edl_sphere_plane(distp,
                                 rpart,
                                 cs_lagr_dlvo_param.valen,
                                 cs_lagr_dlvo_param.phi_p,
                                 cs_lagr_dlvo_param.phi_s,
                                 cs_lagr_dlvo_param.temperature[iel],
                                 cs_lagr_dlvo_param.debye_length[iel],
                                 cs_lagr_dlvo_param.water_permit);

    barr = (var1 + var2);

    if (barr >  *energy_barrier)
      *energy_barrier = barr;
    if ( *energy_barrier < 0)
      *energy_barrier = 0;
   }

  *energy_barrier = *energy_barrier / rpart;
}

/*----------------------------------------------------------------------------
 * Compute the energy barrier for two particles.
 *----------------------------------------------------------------------------*/

void
cs_lagr_barrier_pp(cs_real_t                       dpart,
                   cs_lnum_t                       iel,
                   cs_real_t                      *energy_barrier)
{
  cs_real_t rpart = dpart * 0.5;

  *energy_barrier = 0.;

  /* Computation of the energy barrier */

  for (int i = 0; i < 1001; i++) {

    cs_real_t  step = cs_lagr_dlvo_param.debye_length[iel]/30.0;

    /* Interaction between two spheres */

    cs_real_t distcc = _d_cut_off + i * step + 2.0 * rpart;

    cs_real_t var1
      = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                            rpart,
                                            rpart,
                                            cs_lagr_dlvo_param.lambda_vdw,
                                            cs_lagr_dlvo_param.csthpp);

    cs_real_t var2
      = cs_lagr_edl_sphere_sphere(distcc,
                                  rpart,
                                  rpart,
                                  cs_lagr_dlvo_param.valen,
                                  cs_lagr_dlvo_param.phi_p,
                                  cs_lagr_dlvo_param.phi_p,
                                  cs_lagr_dlvo_param.temperature[iel],
                                  cs_lagr_dlvo_param.debye_length[iel],
                                  cs_lagr_dlvo_param.water_permit);

    cs_real_t barr = (var1 + var2);

    if (barr >  *energy_barrier)
      *energy_barrier = barr;
    if (*energy_barrier < 0)
      *energy_barrier = 0;
   }

  *energy_barrier = *energy_barrier / rpart;
}

/*----------------------------------------------------------------------------
 * Van der Waals interaction between a sphere and a plane
 * using formulas from Czarnecki (large distances)
 *                 and Gregory (small distances)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_plane (cs_real_t distp,
                                    cs_real_t rpart,
                                    cs_real_t lambda_vdw,
                                    cs_real_t cstham)
{
  cs_real_t var;

  if (distp < (lambda_vdw / 2 / _pi)) {
    var = -cstham * rpart / (6 * distp)
          * (1 / (1 + 14 * distp / lambda_vdw + 5 * _pi/4.9
          *  pow(distp,3) / lambda_vdw / pow(rpart,2)));
  }
  else {
    var = cstham
      * ((2.45 * lambda_vdw ) /(60. * _pi)
         * (  (distp - rpart) / pow(distp,2)
            - (distp + 3. * rpart) / pow(distp + 2. * rpart,2))
            - 2.17 / 720. / pow(_pi,2) * pow(lambda_vdw,2) * ((distp - 2. * rpart)
            / pow(distp,3) - (distp + 4. * rpart) / pow(distp + 2. * rpart,3))
            + 0.59 / 5040. / pow(_pi,3) * pow(lambda_vdw,3)
            * ((distp - 3. * rpart) /  pow(distp,4) - (distp + 5. * rpart)
            / pow(distp + 2. * rpart,4)));
  }

  return var;
}

/*----------------------------------------------------------------------------
 * Calculation of the Van der Waals interaction between two spheres
 * following the formula from Gregory (1981a)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_sphere(cs_real_t  distcc,
                                    cs_real_t  rpart1,
                                    cs_real_t  rpart2,
                                    cs_real_t  lambda_vdw,
                                    cs_real_t  cstham)
{
  cs_real_t var = - cstham * rpart1 * rpart2 / (6 * (distcc - rpart1 - rpart2)
               * (rpart1 + rpart2)) * (1 - 5.32 * (distcc - rpart1 - rpart2)
              / lambda_vdw * log(1 + lambda_vdw / (distcc - rpart1 - rpart2) / 5.32));

  return var;
}

/*----------------------------------------------------------------------------
 * Electric Double Layer (EDL) interaction between a sphere and a plane
 * using the formula from Bell & al (1970)
 * based on the McCartney & Levine method
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_edl_sphere_plane(cs_real_t  distp,
                         cs_real_t  rpart,
                         cs_real_t  valen,
                         cs_real_t  phi1,
                         cs_real_t  phi2,
                         cs_real_t  temp,
                         cs_real_t  debye_length,
                         cs_real_t  water_permit)
{

  cs_real_t charge = 1.6e-19;
  /* Reduced zeta potential */
  cs_real_t lphi1 =  valen * charge * phi1 /  _k_boltzmann / temp;
  cs_real_t lphi2 =  valen * charge * phi2 /  _k_boltzmann / temp;

  cs_real_t tau = rpart / debye_length;

  /* Extended reduced zeta potential */
  /* (following the work from Ohshima et al, 1982, JCIS, 90, 17-26) */

  lphi1 = 8. * tanh(lphi1 / 4.) /
         ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
         * pow(tanh(lphi1 / 4.),2),0.5));

  lphi2 = 4. * tanh(lphi2 / 4.) ;

  cs_real_t alpha =   sqrt((distp + rpart) / rpart)
                    + sqrt(rpart / (distp + rpart));
  cs_real_t omega1 = pow(lphi1,2) + pow(lphi2,2) + alpha * lphi1 * lphi2;
  cs_real_t omega2 = pow(lphi1,2) + pow(lphi2,2) - alpha * lphi1 * lphi2;
  cs_real_t gamma = sqrt(rpart / (distp + rpart)) * exp(-1./debye_length * distp);

  cs_real_t var = 2 * _pi * _free_space_permit * water_permit
    * pow((_k_boltzmann * temp / (1. * valen) / charge),2)
                  * rpart * (distp + rpart) / (distp + 2 * rpart)
                  * (omega1 * log(1 + gamma) + omega2 * log(1 - gamma));

  return var;
}

/*----------------------------------------------------------------------------
 * Calculation of the EDL interaction between two spheres
 * using the formula from Bell & al (1970)
 * based on the McCartney & Levine method
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_edl_sphere_sphere(cs_real_t  distcc,
                          cs_real_t  rpart1,
                          cs_real_t  rpart2,
                          cs_real_t  valen,
                          cs_real_t  phi1,
                          cs_real_t  phi2,
                          cs_real_t  temp,
                          cs_real_t  debye_length,
                          cs_real_t  water_permit)
{
  cs_real_t charge = 1.6e-19;

  /* Reduced zeta potential */
  cs_real_t lphi1 =  valen * charge * phi1 /  _k_boltzmann / temp;
  cs_real_t lphi2 =  valen * charge * phi2 /  _k_boltzmann / temp;


  /* Extended reduced zeta potential */
  /* (following the work from Ohshima et al, 1982, JCIS, 90, 17-26) */

  cs_real_t tau1 = rpart1 / debye_length;
  lphi1 = 8. * tanh(lphi1 / 4.) /
    ( 1. + pow(1. - (2. * tau1 + 1.) / (pow(tau1 + 1,2))
               * pow(tanh(lphi1 / 4.),2),0.5));

  cs_real_t tau2 = rpart2 / debye_length;
  lphi2 = 8. * tanh(lphi2 / 4.) /
         ( 1. + pow(1. - (2. * tau2 + 1.) / (pow(tau2 + 1,2))
          * pow(tanh(lphi2 / 4.),2),0.5));

  cs_real_t alpha =    sqrt(rpart2 * (distcc - rpart2)
                    / (rpart1 * (distcc - rpart1)))
                     + sqrt(rpart1 * (distcc - rpart1)
                    / (rpart2 * (distcc - rpart2)));

  cs_real_t omega1 = pow(lphi1,2) + pow(lphi2,2) + alpha * lphi1 * lphi2;

  cs_real_t omega2 = pow(lphi1,2) + pow(lphi2,2) - alpha * lphi1 * lphi2;

  cs_real_t gamma = sqrt(rpart1 * rpart2 / (distcc-rpart1) / (distcc-rpart2))
                    *exp(1. / debye_length * (rpart1 + rpart2 - distcc));

  cs_real_t var = 2 * _pi * _free_space_permit * water_permit
                  * pow((_k_boltzmann * temp / charge),2)
                  * rpart1 * rpart2 * (distcc - rpart1) * (distcc - rpart2)
                  / (distcc * (  distcc * (rpart1  + rpart2)
                               - pow(rpart1,2) - pow(rpart2,2)))
                  * (omega1 * log(1 + gamma) + omega2 * log(1 - gamma));

  return var;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
