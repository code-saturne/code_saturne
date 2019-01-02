/*============================================================================
 * Surface roughness
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_random.h"
#include "cs_search.h"
#include "cs_halo.h"
#include "cs_lagr_dlvo.h"

#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_roughness.h"

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

static cs_lagr_roughness_param_t _lagr_roughness_param
= {.water_permit = 0,
   .ionic_strength = 0,
   .phi_p = 0,
   .phi_s = 0,
   .temperature = NULL,
   .valen = 0,
   .debye_length = NULL,
   .cstham = 0,
   .lambda_vdw = 0,
   .espasg = 0,
   .denasp = 0,
   .rayasp = 0,
   .rayasg = 0};

cs_lagr_roughness_param_t *cs_lagr_roughness_param = &_lagr_roughness_param;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double _pi = 3.14159265358979323846;

/* Free space permittivity */
static const cs_real_t _free_space_permit = 8.854e-12;

/* Faraday constant */
static const cs_real_t _faraday_cst = 9.648e4;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Roughness initialization:
 *  - Retrieve various parameters for storing in global structure
 *  - Compute and store the Debye screening length
 *----------------------------------------------------------------------------*/

void
roughness_init (const cs_real_t   *water_permit,
                const cs_real_t   *ionic_strength,
                const cs_real_t    temperature[],
                const cs_real_t   *valen,
                const cs_real_t   *phi_p,
                const cs_real_t   *phi_s,
                const cs_real_t   *cstham,
                const cs_real_t   *lambda_vdw,
                const cs_real_t   *espasg,
                const cs_real_t   *denasp,
                const cs_real_t   *rayasp,
                const cs_real_t   *rayasg)
{
  int iel;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clog_param          */

  cs_lagr_roughness_param->water_permit = *water_permit;
  cs_lagr_roughness_param->ionic_strength = *ionic_strength;
  cs_lagr_roughness_param->valen = *valen;
  cs_lagr_roughness_param->phi_p = *phi_p;
  cs_lagr_roughness_param->phi_s = *phi_s;
  cs_lagr_roughness_param->cstham = *cstham;
  cs_lagr_roughness_param->cstham = *lambda_vdw;
  cs_lagr_roughness_param->espasg = *espasg;
  cs_lagr_roughness_param->denasp = *denasp;
  cs_lagr_roughness_param->rayasp = *rayasp;
  cs_lagr_roughness_param->rayasg = *rayasg;

  /* Allocate memory for the temperature and Debye length arrays */

  if (cs_lagr_roughness_param->temperature == NULL)
    BFT_MALLOC(cs_lagr_roughness_param->temperature, mesh->n_cells, cs_real_t);

  if (cs_lagr_roughness_param->debye_length == NULL)
    BFT_MALLOC(cs_lagr_roughness_param->debye_length, mesh->n_cells, cs_real_t);

  /* Store the temperature */

  for (iel = 0; iel < mesh->n_cells; iel++)
    cs_lagr_roughness_param->temperature[iel] = temperature[iel];

  /* Computation and storage of the Debye length */

  for (iel = 0; iel < mesh->n_cells ; iel++)

    cs_lagr_roughness_param->debye_length[iel]
      =   pow(2e3 * pow(_faraday_cst,2)
        * cs_lagr_roughness_param->ionic_strength /
        (cs_lagr_roughness_param->water_permit
         * _free_space_permit * PG_CST
         * cs_lagr_roughness_param->temperature[iel]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" epseau = %g\n", cs_lagr_roughness_param->water_permit);
  bft_printf(" fion   = %g\n", cs_lagr_roughness_param->ionic_strength);
  bft_printf(" temp[1]   = %g\n", cs_lagr_roughness_param->temperature[0]);
  bft_printf(" valen   = %g\n", cs_lagr_roughness_param->valen);
  bft_printf(" debye[1]   = %g\n", cs_lagr_roughness_param->debye_length[0]);
  bft_printf(" phi_p   = %g\n", cs_lagr_roughness_param->phi_p);
  bft_printf(" phi_s  = %g\n", cs_lagr_roughness_param->phi_s);
#endif

}

/*----------------------------------------------------------------------------
 * Deallocate the arrays storing temperature and Debye length.
 *----------------------------------------------------------------------------*/

void
cs_lagr_roughness_finalize()
{
  BFT_FREE(cs_lagr_roughness_param->temperature);
  BFT_FREE(cs_lagr_roughness_param->debye_length);
}

/*----------------------------------------------------------------------------
 * Compute the energy barrier for a rough wall.
 *
 * parameters:
 *   particle       <-- pointer to particle data
 *   attr_map       <-- pointer to attribute map
 *   iel            <-- id of cell where the particle is
 *   energy_barrier <-> energy barrier
 *----------------------------------------------------------------------------*/

void
cs_lagr_roughness_barrier(const void                     *particle,
                          const cs_lagr_attribute_map_t  *attr_map,
                          cs_lnum_t                       iel,
                          cs_real_t                      *energy_barrier)
{
  cs_int_t  i;
  cs_real_t param2, value;
  cs_lnum_t param1,contact,compt_max;
  cs_lnum_t iclas, ints, np, iasp;
  cs_real_t rpart2[2],udlvor[500];
  cs_real_t distasp, posasp1[2000], posasp2[2000],
    posasp3[2000], posasp4[2000], disminp;
  cs_real_t scov[2], seff[1];
  cs_lnum_t nbtemp[12000];
  cs_lnum_t nbasp[1], nclas, nbaspt[1], nasptot;
  cs_real_t random;

  contact = 0;
  compt_max = 5000;

  /* Computation of the surface coverage */

  nclas = 2;
  cs_real_t scovtot = 0.;
  for (iclas = 0; iclas < nclas; iclas++) {
    scov[0] = _pi *  pow(cs_lagr_roughness_param->rayasg, 2)
      / pow(cs_lagr_roughness_param->espasg, 2);

    scov[1] =  cs_lagr_roughness_param->denasp
      * _pi * pow(cs_lagr_roughness_param->rayasp, 2);

    scovtot = scovtot + scov[iclas];

    rpart2[0] = cs_lagr_roughness_param->rayasg;
    rpart2[1] =  cs_lagr_roughness_param->rayasp;
  }


  cs_real_t rpart = cs_lagr_particle_get_real(particle, attr_map,
                                              CS_LAGR_DIAMETER) * 0.5;

  /* Creation of asperities */

  for (iclas = 0; iclas < nclas; iclas++) {
    seff[iclas] = 0.;
    nbasp[iclas] = 0;
  }

  nasptot = 0;
  cs_lnum_t nasp = 0;
  cs_real_t dismin = 0.;

  for (iclas = 0; iclas < nclas; iclas++) {
      rpart2[0] = cs_lagr_roughness_param->rayasg;
      rpart2[1] =  cs_lagr_roughness_param->rayasp;

    seff[iclas]  = 2.5 * _pi * (2. * rpart + rpart2[iclas] + 10.
                                *  cs_lagr_roughness_param->debye_length[iel])
                       *  (  rpart2[iclas] + 10.
                           * cs_lagr_roughness_param->debye_length[iel]);

    value = 700.;

    cs_real_t value2 = seff[iclas] * scov[iclas] / _pi / pow(rpart2[iclas],2);

    if (value2 > 700) {
      param1 = value2 / 700;
      param2 = fmod(value2 , 700.);
      cs_random_poisson(1, param2, nbaspt);
      cs_random_poisson(param1, value , nbtemp);
      for (ints = 0; ints < param1; ints++) {
        nbaspt[0] = nbaspt[0] + nbtemp[ints];
      }
      nbasp[iclas] = nbaspt[0];
    }
    else {
      cs_random_poisson(1, value2 , nbaspt);
      nbasp[iclas] = nbaspt[0];
    }

    /* Placement of asperities */
    cs_lnum_t iboucle;

    for (i = 0; i < nbasp[iclas]; i++) {
      iboucle  = 0;
      do {
        contact = 0;
        cs_random_uniform(1, &random);

        posasp1[i + nasptot] = pow(seff[iclas] /_pi, 0.5) * random;
        posasp2[i + nasptot] = 2 * _pi * random;
        posasp3[i + nasptot] = 0.;
        posasp4[i + nasptot] = rpart2[iclas];

        /* No contact between two asperities of a same class */
        for (iasp = 0; iasp < nasp - nasptot; iasp++) {

          distasp = pow(posasp1[i + nasptot]
                        * cos(posasp2[i + nasptot])
                        - posasp1[iasp + nasptot] * cos(posasp2[iasp + nasptot]), 2)
                  + pow(posasp1[i + nasptot]
                        * sin(posasp2[i + nasptot])
                        - posasp1[iasp + nasptot] * sin(posasp2[iasp + nasptot]), 2)
                  + pow(posasp3[iasp + nasptot] - posasp3[i + nasptot], 2);

          if (distasp <  pow(posasp4[iasp + nasptot] + posasp4[i + nasptot], 2)) {
            iboucle = iboucle + 1;
            contact = contact + 1;
          }
        }
      } while (contact != 0 && iboucle < compt_max);

      if (iboucle > compt_max) {
        cs_random_uniform(1, &random);

        posasp1[i + nasptot] = pow(seff[iclas] / _pi,0.5) * 2.;
        posasp2[i + nasptot] = 2 * _pi * random;
        posasp3[i + nasptot] = 0.;
        posasp4[i + nasptot] = rpart2[iclas];
      }

   /* No contact between two asperities of various class */
      for (iasp = 0; iasp <  nasptot; iasp++) {
        distasp =  pow(  posasp1[i+nasptot] *  cos(posasp2[i+nasptot])
                       - posasp1[iasp]      * cos(posasp1[iasp]), 2)
                 + pow(  posasp1[i+nasptot] * sin(posasp2[i+nasptot])
                       - posasp1[iasp]      * sin(posasp2[iasp]), 2)
                 + pow(posasp3[iasp] - posasp3[i+nasptot], 2);

        if (  distasp < pow(posasp4[iasp], 2) && pow(posasp3[i + nasptot], 2)
            < (pow(posasp4[iasp], 2) - distasp)) {
          posasp3[i + nasptot] =   sqrt(pow(posasp4[iasp], 2) - distasp)
                                 +  posasp3[iasp];
        }
      }

      nasp = nasp + 1;
    }

   /* Number of asperities on the surface */
    nasptot = nasptot + nbasp[iclas];
    /* End of the loop on asperity size */
  }


  /* Determination of the mimnimal distance*/
  for (iasp = 0 ; iasp < nasptot;iasp++) {
    if ( posasp1[iasp] < (rpart + posasp4[iasp])) {
      disminp =   sqrt(pow(rpart + posasp4[iasp], 2) - pow(posasp1[iasp], 2))
                - rpart + posasp3[iasp];
    }
    else {
      disminp = 0.;
    }
    if (disminp > dismin) {
      dismin = disminp;
    }
  }

   /* Calculation of the energy barrier */

   /* Loop on the separation distance */
  for (np = 0; np <  500; np++) {
    udlvor[np] = 0.;
    cs_real_t distp =   dismin + (np + 1)
                      * cs_lagr_roughness_param->debye_length[iel]/30.0;

    /* DLVO between the particle and the rough plate */

    /* Sum of the interaction {particle-plate} and {particule-asperity} */

    /* Sphere-plate interaction */
    cs_real_t var1
      = cs_lagr_van_der_waals_sphere_plane(distp,
                                           rpart,
                                           cs_lagr_roughness_param->cstham,
                                           cs_lagr_roughness_param->lambda_vdw);

    cs_real_t var2
      = cs_lagr_edl_sphere_plane(distp,
                                 rpart,
                                 cs_lagr_roughness_param->valen,
                                 cs_lagr_roughness_param->phi_p,
                                 cs_lagr_roughness_param->phi_s,
                                 cs_lagr_roughness_param->temperature[iel],
                                 cs_lagr_roughness_param->debye_length[iel],
                                 cs_lagr_roughness_param->water_permit);

    udlvor[np] = (var1 + var2) * (1. - scovtot);

    /* Sphere-asperity interactions */

    for (iasp = 0; iasp <  nasptot; iasp++) {
      cs_real_t distcc = sqrt(  pow(distp + rpart- posasp3[iasp], 2)
                              + pow(posasp1[iasp],2));

      var1
        = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                              rpart,
                                              posasp4[iasp],
                                              cs_lagr_roughness_param->cstham,
                                              cs_lagr_roughness_param->lambda_vdw);

      var2 = cs_lagr_edl_sphere_sphere(distcc,
                                       rpart,
                                       posasp4[iasp],
                                       cs_lagr_roughness_param->valen,
                                       cs_lagr_roughness_param->phi_p,
                                       cs_lagr_roughness_param->phi_s,
                                       cs_lagr_roughness_param->temperature[iel],
                                       cs_lagr_roughness_param->debye_length[iel],
                                       cs_lagr_roughness_param->water_permit);

      udlvor[np] =   udlvor[np] + (var1 + var2)
                   * (distp + rpart - posasp3[iasp]) / distcc;
    }

  } /* End of the loop on the separation distance */

  /* Tracking of the energy barrier */
  cs_real_t barren = 0.;
  for (np = 0; np <  500; np++) {
    if (udlvor[np] > barren) {
      barren = udlvor[np];
    }
  }
  if (barren < 0.) barren = 0.;
  *energy_barrier = barren / rpart;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
