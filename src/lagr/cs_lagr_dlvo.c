/*============================================================================
 * Calculation of DLVO forces
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
#include "cs_lagr_utils.h"
#include "cs_halo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_dlvo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

const double _pi = 4 * atan(1);

/*============================================================================
 * Public function for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Van der Waals interaction between a sphere and a plane
 * using formulas from Czarnecki (large distances)
 *                 and Gregory (small distances)
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_van_der_waals_sphere_plane (cs_real_t distp,
                                    cs_real_t rpart,
                                    cs_real_t lambwl,
                                    cs_real_t cstham)
{
  cs_real_t var;

  if (distp < (lambwl / 2 / _pi)) {
    var = -cstham * rpart / (6 * distp)
          * (1 / (1 + 14 * distp / lambwl + 5 * _pi/4.9
          *  pow(distp,3) / lambwl / pow(rpart,2)));
  }
  else {
    var = cstham
      * ((2.45 * lambwl ) /(60 * _pi)
         * (  (distp - rpart) / pow(distp,2)
            - (distp + 3 * rpart) / pow(distp + 2 * rpart,2)
            - 2.17 / 720 / pow(_pi,2) * pow(lambwl,2) * ((distp - 2 * rpart)
            / pow(distp,3) - (distp + 4 * rpart) / pow(distp + 2 * rpart,3))
            + 0.59 / 5040 / pow(_pi,3) * pow(lambwl,3)
            * ((distp - 3 * rpart) /  pow(distp,4) - (distp + 5 * rpart)
            / pow(distp + 2 * rpart,4))));
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
                                    cs_real_t  lambwl,
                                    cs_real_t  cstham)
{
  cs_real_t var = - cstham * rpart1 * rpart2 / (6 * (distcc - rpart1 - rpart2)
               * (rpart1 + rpart2)) * (1 - 5.32 * (distcc - rpart1 - rpart2)
              / lambwl * log(1 + lambwl / (distcc - rpart1 - rpart2) / 5.32));

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
                         cs_real_t  phi1,
                         cs_real_t  phi2,
                         cs_real_t  kboltz,
                         cs_real_t  temp,
                         cs_real_t  debye_length,
                         cs_real_t  free_space_permit,
                         cs_real_t  water_permit)
{
  cs_real_t charge = 1.6e-19;
  /* Reduced zeta potential */
  cs_real_t lphi1 =  charge * phi1 /  kboltz / temp;
  cs_real_t lphi2 =  charge * phi2 /  kboltz / temp;

  cs_real_t tau = rpart / (1. / debye_length);

  /* Extended reduced zeta potential */
  /* (following the work from Ohshima et al, 1982, JCIS, 90, 17-26) */

  lphi1 = 8. * tanh(lphi1 / 4.) /
         ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
         * pow(tanh(lphi1 / 4.),2),0.5));

  lphi2 = 8. * tanh(lphi2 / 4.) /
          ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
          * pow(tanh(lphi2 / 4.),2),0.5));

  cs_real_t alpha =   sqrt((distp + rpart) / rpart)
                    + sqrt(rpart / (distp + rpart));
  cs_real_t omega1 = pow(lphi1,2) + pow(lphi2,2) + alpha * lphi1 * lphi2;
  cs_real_t omega2 = pow(lphi1,2) + pow(lphi2,2) - alpha * lphi1 * lphi2;
  cs_real_t gamma = sqrt(rpart / (distp + rpart)) * exp(-1./debye_length * distp);

  cs_real_t var = 2 * _pi * free_space_permit * water_permit
                  * pow((kboltz * temp / charge),2)
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
                          cs_real_t  phi1,
                          cs_real_t  phi2,
                          cs_real_t  kboltz,
                          cs_real_t  temp,
                          cs_real_t  debye_length,
                          cs_real_t  free_space_permit,
                          cs_real_t  water_permit)
{
  cs_real_t charge = 1.6e-19;

  /* Reduced zeta potential */
  cs_real_t lphi1 =  charge * phi1 /  kboltz / temp;
  cs_real_t lphi2 =  charge * phi2 /  kboltz / temp;

  cs_real_t tau = rpart1 / (1. / debye_length);

  /* Extended reduced zeta potential */
  /* (following the work from Ohshima et al, 1982, JCIS, 90, 17-26) */

  lphi1 = 8. * tanh(lphi1 / 4.) /
    ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
               * pow(tanh(lphi1 / 4.),2),0.5));

  lphi2 = 8. * tanh(lphi2 / 4.) /
         ( 1. + pow(1. - (2. * tau + 1.) / (pow(tau + 1,2))
          * pow(tanh(lphi2 / 4.),2),0.5));

  cs_real_t alpha =    sqrt(rpart2 * (distcc - rpart2)
                    / (rpart1 * (distcc - rpart1)))
                     + sqrt(rpart1 * (distcc - rpart1)
                    / (rpart2 * (distcc - rpart2)));

  cs_real_t omega1 = pow(lphi1,2) + pow(lphi2,2) + alpha * lphi1 * lphi2;

  cs_real_t omega2 = pow(lphi1,2) + pow(lphi2,2) - alpha * lphi1 * lphi2;

  cs_real_t gamma = sqrt(rpart1 * rpart2 / (distcc-rpart1) / (distcc-rpart2))
                    *exp(1. / debye_length * (rpart1 + rpart2 - distcc));

  cs_real_t var = 2 * _pi * free_space_permit * water_permit
                  * pow((kboltz * temp / charge),2)
                  * rpart1 * rpart2 * (distcc - rpart1) * (distcc - rpart2)
                  / (distcc * (  distcc * (rpart1  + rpart2)
                               - pow(rpart1,2) - pow(rpart2,2)))
                  * (omega1 * log(1 + gamma) + omega2 * log(1 - gamma));

  return var;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
