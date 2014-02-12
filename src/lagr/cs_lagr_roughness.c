/*============================================================================
 * Methods for roughness surface
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
 * Functions dealing with the roughness surface modeling
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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_roughness.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro declarations
 *============================================================================*/

#define PG_CST 8.314  /* Ideal gas constant */

/*============================================================================
 * Local structure declarations
 *============================================================================*/

static cs_lagr_roughness_param_t cs_lagr_roughness_param;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double _pi = 3.14159265358979323846;

/*============================================================================
 * Public function for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Roughness initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
CS_PROCF (roughness_init, ROUGHNESS_INIT)(const cs_real_t   *faraday_cst,
                                          const cs_real_t   *free_space_permit,
                                          const cs_real_t   *water_permit,
                                          const cs_real_t   *ionic_strength,
                                          const cs_real_t    temperature[],
                                          const cs_real_t   *phi1,
                                          const cs_real_t   *phi2,
                                          const cs_real_t   *cstham,
                                          const cs_real_t   *dcutof,
                                          const cs_real_t   *lambwl,
                                          const cs_real_t   *kboltz,
                                          const cs_real_t   *espasg,
                                          const cs_real_t   *denasp,
                                          const cs_real_t   *rayasp,
                                          const cs_real_t   *rayasg)
{
  int ifac;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clog_param          */

  cs_lagr_roughness_param.faraday_cst = *faraday_cst;
  cs_lagr_roughness_param.free_space_permit = *free_space_permit;
  cs_lagr_roughness_param.water_permit = *water_permit;
  cs_lagr_roughness_param.ionic_strength = *ionic_strength;
  cs_lagr_roughness_param.phi1 = *phi1;
  cs_lagr_roughness_param.phi2 = *phi2;
  cs_lagr_roughness_param.cstham = *cstham;
  cs_lagr_roughness_param.dcutof = *dcutof;
  cs_lagr_roughness_param.lambwl = *lambwl;
  cs_lagr_roughness_param.kboltz = *kboltz;
  cs_lagr_roughness_param.espasg = *espasg;
  cs_lagr_roughness_param.denasp = *denasp;
  cs_lagr_roughness_param.rayasp = *rayasp;
  cs_lagr_roughness_param.rayasg = *rayasg;

  /* Allocate memory for the temperature and Debye length arrays */

  if (cs_lagr_roughness_param.temperature == NULL)
    BFT_MALLOC(cs_lagr_roughness_param.temperature, mesh->n_b_faces, cs_real_t);

  if (cs_lagr_roughness_param.debye_length == NULL)
    BFT_MALLOC(cs_lagr_roughness_param.debye_length, mesh->n_b_faces, cs_real_t);

  /* Store the temperature */

  for (ifac = 0; ifac < mesh->n_b_faces; ifac++)
    cs_lagr_roughness_param.temperature[ifac] = temperature[ifac];

  /* Computation and storage of the Debye length */

  for (ifac = 0; ifac < mesh->n_b_faces ; ifac++)

    cs_lagr_roughness_param.debye_length[ifac]
      =   pow(2e3 * pow(cs_lagr_roughness_param.faraday_cst,2)
        * cs_lagr_roughness_param.ionic_strength /
        (cs_lagr_roughness_param.water_permit
         * cs_lagr_roughness_param.free_space_permit * PG_CST
         * cs_lagr_roughness_param.temperature[ifac]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" cstfar = %g\n", cs_lagr_roughness_param.faraday_cst);
  bft_printf(" epsvid = %g\n", cs_lagr_roughness_param.free_space_permit);
  bft_printf(" epseau = %g\n", cs_lagr_roughness_param.water_permit);
  bft_printf(" fion   = %g\n", cs_lagr_roughness_param.ionic_strength);
  bft_printf(" temp[1]   = %g\n", cs_lagr_roughness_param.temperature[0]);
  bft_printf(" debye[1]   = %g\n", cs_lagr_roughness_param.debye_length[0]);
  bft_printf(" phi1   = %g\n", cs_lagr_roughness_param.phi1);
  bft_printf(" phi2  = %g\n", cs_lagr_roughness_param.phi2);
#endif

}

/*----------------------------------------------------------------------------
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_roughness_finalize()
{
  BFT_FREE(cs_lagr_roughness_param.temperature);
  BFT_FREE(cs_lagr_roughness_param.debye_length);
}

/*----------------------------------------------------------------------------
 * Compute the energy barrier for a rough wall.
 *----------------------------------------------------------------------------*/

cs_real_t
cs_lagr_roughness_barrier(cs_lagr_particle_t   particle,
                          cs_int_t             face_id,
                          cs_real_t           *energy_barrier)
{
  cs_int_t  dim_aux = 1;
  cs_real_t param2, value;
  cs_int_t  param1;
  cs_lnum_t k, i;
  cs_real_t dismin;
  cs_real_t nbr1, nbr2;
  cs_real_t nmoyap ;
  cs_int_t  nbrasg[1],nbrasp[1],ntmp[1],nbtemp[12000];
  cs_int_t  nbasg,nbasp;

  /* Computation of the number of particles in contact with */
  /* the depositing particle */

  /* Assuming monodisperse calculation */

  cs_real_t scovap =   cs_lagr_roughness_param.denasp
                     * _pi * pow(cs_lagr_roughness_param.rayasp, 2);
  cs_real_t scovag = _pi *  pow(cs_lagr_roughness_param.rayasg, 2)
                          / pow(cs_lagr_roughness_param.espasg, 2);

  cs_real_t scovtot = scovap + scovag;

  /* Number of large-scale asperities*/

  cs_real_t rpart = 0.5 * particle.diameter;

  cs_real_t nmoyag = (2.0 * rpart + cs_lagr_roughness_param.rayasg)
                     / cs_lagr_roughness_param.rayasg * scovag;

  CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyag, nbrasg);

  value = 700.;

  if (nmoyag > value) {
    param1 = nmoyag / value;
    param2 = fmod(nmoyag,value);

    CS_PROCF(fische, FISCHE)(&dim_aux, &param2, nbrasg);
    CS_PROCF(fische, FISCHE)(&param1, &value ,nbtemp);

    for (k = 0; k < param1; k++)
      nbrasg[0] =  nbrasg[0] + nbtemp[k];
  }

  if (nbrasg[0] > 1) {
    nmoyag =  1 + 2 * cs_lagr_roughness_param.dcutof
              * (2.0 * rpart + 2.0 * cs_lagr_roughness_param.rayasg + 4.0
                 * cs_lagr_roughness_param.dcutof)
              / pow(cs_lagr_roughness_param.rayasg,2) * scovag;

    CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyag, ntmp);
    nbasg = ntmp[0];
    if (nbasg < 1)
      nbasg = 1;
  }
  else
    nbasg = nbrasg[0];

  /* Nb of small-scale asperities */

  /* 1st case: no large-scale asperities */
  if (nbasg == 0) {
    nmoyap = (2.0 * rpart + cs_lagr_roughness_param.rayasp)
             / cs_lagr_roughness_param.rayasp * scovap;

    CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyap, nbrasp);

    value = 700.;
    if (nmoyap > value) {
      param1 = nmoyap / value;
      param2 = fmod(nmoyap,value);

      CS_PROCF(fische, FISCHE)(&dim_aux, &param2, nbrasp);
      CS_PROCF(fische, FISCHE)(&param1, &value ,nbtemp);

      for (k = 0; k < param1; k++)
        nbrasp[0] =  nbrasp[0] + nbtemp[k];
    }

    if (nbrasp[0] > 1) {
      nmoyap =  1 + 2 * cs_lagr_roughness_param.dcutof
                * (2.0 * rpart + 2.0 * cs_lagr_roughness_param.rayasp + 4.0
                   * cs_lagr_roughness_param.dcutof)
              / pow(cs_lagr_roughness_param.rayasp,2) * scovap;

      CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyap, ntmp);
      nbasp = ntmp[0];
      if (nbasp < 1)
        nbasp = 1;
    }
    else
      nbasp = nbrasp[0];

    /* Determine the minimal distance between the particle and the plate */

    if (nbrasp[0] < 1)
      nbr1 = nbrasp[0];
    else
      nbr1 = 1;
    dismin = cs_lagr_roughness_param.rayasp * nbr1;
  }
  else {
    cs_real_t paramh = 0.5 * (2.0 * rpart + cs_lagr_roughness_param.rayasp)
                       * cs_lagr_roughness_param.rayasp
                       / (rpart + cs_lagr_roughness_param.rayasg);

    nmoyap = 1 + paramh * (2 * cs_lagr_roughness_param.rayasg - paramh)
            / pow(cs_lagr_roughness_param.rayasp,2) * scovap;

    CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyap, nbrasp);

    if (nbrasp[0] > 1) {
      paramh = 0.5 * (  2. * rpart + 2 * cs_lagr_roughness_param.rayasp
                      + 4.0 * cs_lagr_roughness_param.dcutof ) * 2.0
        * cs_lagr_roughness_param.dcutof
        / (  rpart + cs_lagr_roughness_param.rayasg
           + cs_lagr_roughness_param.rayasp + cs_lagr_roughness_param.dcutof);
      nmoyap =  1 + paramh * (2 * cs_lagr_roughness_param.rayasg - paramh)
             / pow( cs_lagr_roughness_param.rayasp,2) * scovap;
      CS_PROCF(fische, FISCHE)(&dim_aux, &nmoyap, ntmp);
      nbasp = ntmp[0];
      if (nbasp < 1 )
        nbasp = 1;
    }
    else
      nbasp = nbrasp[0];

    nbasp = nbasp * nbasg;
    nbrasp[0] = nbrasp[0] * nbrasg[0];
    if (nbasp < 1)
      nbr1 = nbasp;
    else
      nbr1 = 1;

    if (nbasg < 1)
      nbr2 = nbasg;
    else
      nbr2 = 1;

    dismin =   cs_lagr_roughness_param.rayasp * nbr1
             + cs_lagr_roughness_param.rayasg * nbr2;
  }

  *energy_barrier = 0.;

  cs_real_t barr = 0.;
  cs_real_t distp = 0.;

  /* Computation of the energy barrier */

  for (i = 0; i < 101; i++) {

    cs_real_t  step = 1e-10;

    /* Interaction between the sphere and the plate */

    distp = dismin + cs_lagr_roughness_param.dcutof + i * step;

    cs_real_t var1
      = cs_lagr_van_der_waals_sphere_plane(distp,
                                           rpart,
                                           cs_lagr_roughness_param.lambwl,
                                           cs_lagr_roughness_param.cstham);

    cs_real_t var2
      = cs_lagr_edl_sphere_plane(distp,
                                 rpart,
                                 cs_lagr_roughness_param.phi1,
                                 cs_lagr_roughness_param.phi2,
                                 cs_lagr_roughness_param.kboltz,
                                 cs_lagr_roughness_param.temperature[face_id],
                                 cs_lagr_roughness_param.debye_length[face_id],
                                 cs_lagr_roughness_param.free_space_permit,
                                 cs_lagr_roughness_param.water_permit);

    barr = (var1 + var2) * (1 - scovtot);

    /* Interaction between the sphere and small-scale asperities */

    cs_real_t distcc =   cs_lagr_roughness_param.dcutof + step * i
                       + rpart + cs_lagr_roughness_param.rayasp;

    var1 = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                               rpart,
                                               cs_lagr_roughness_param.rayasp,
                                               cs_lagr_roughness_param.lambwl,
                                               cs_lagr_roughness_param.cstham);
    var2
      = cs_lagr_edl_sphere_sphere(distcc,
                                  rpart,
                                  cs_lagr_roughness_param.rayasp,
                                  cs_lagr_roughness_param.phi1,
                                  cs_lagr_roughness_param.phi2,
                                  cs_lagr_roughness_param.kboltz,
                                  cs_lagr_roughness_param.temperature[face_id],
                                  cs_lagr_roughness_param.debye_length[face_id],
                                  cs_lagr_roughness_param.free_space_permit,
                                  cs_lagr_roughness_param.water_permit);

    barr = barr + (var1 + var2) * nbasp;

    /* Interaction between the sphere and large-scale asperities */

    distcc =   cs_lagr_roughness_param.dcutof + step * i
             + rpart + cs_lagr_roughness_param.rayasg;

    var1 = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                               rpart,
                                               cs_lagr_roughness_param.rayasg,
                                               cs_lagr_roughness_param.lambwl,
                                               cs_lagr_roughness_param.cstham);
    var2
      = cs_lagr_edl_sphere_sphere(distcc,
                                  rpart,
                                  cs_lagr_roughness_param.rayasg,
                                  cs_lagr_roughness_param.phi1,
                                  cs_lagr_roughness_param.phi2,
                                  cs_lagr_roughness_param.kboltz,
                                  cs_lagr_roughness_param.temperature[face_id],
                                  cs_lagr_roughness_param.debye_length[face_id],
                                  cs_lagr_roughness_param.free_space_permit,
                                  cs_lagr_roughness_param.water_permit);

    barr = barr  + (var1 + var2) * nbasg;

    if (barr >  *energy_barrier)
      *energy_barrier = barr;
    if (barr < 0)
      *energy_barrier = 0;
   }

  *energy_barrier = *energy_barrier / (0.5 * particle.diameter);

  return  *energy_barrier;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
