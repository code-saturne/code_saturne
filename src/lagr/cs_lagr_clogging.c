/*============================================================================
 * Methods for particle clogging modeling
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
#include "cs_random.h"
#include "cs_search.h"
#include "cs_halo.h"

#include "cs_lagr.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_roughness.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_clogging.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local structure declaration
 *============================================================================*/

static cs_lagr_clogging_param_t cs_lagr_clogging_param;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double _pi = 3.14159265358979323846;

/* Cut-off distance for adhesion forces (assumed to be the Born distance) */
static const cs_real_t  _d_cut_off = 1.65e-10;

/* Free space permittivity */
static const cs_real_t _free_space_permit = 8.854e-12;

/* Faraday constant */
static const cs_real_t _faraday_cst = 9.648e4;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
cloginit (const cs_real_t   *water_permit,
          const cs_real_t   *ionic_strength,
          const cs_real_t   *jamming_limit,
          const cs_real_t   *min_porosity,
          const cs_real_t   *diam_mean,
          const cs_real_t    temperature[],
          const cs_real_t   *valen,
          const cs_real_t   *phi_p,
          const cs_real_t   *phi_s,
          const cs_real_t   *cstham,
          const cs_real_t   *csthpp,
          const cs_real_t   *lambda_vdw
          )
{
#define PG_CST 8.314  /* Ideal gas constant */

  cs_lnum_t iel;

  const cs_mesh_t  *mesh = cs_glob_mesh;

  /* Retrieve physical parameters related to clogging modeling */
  /* and fill the global structure cs_lagr_clogging_param          */

  cs_lagr_clogging_param.water_permit = *water_permit;
  cs_lagr_clogging_param.ionic_strength = *ionic_strength;
  cs_lagr_clogging_param.jamming_limit = *jamming_limit;
  cs_lagr_clogging_param.min_porosity = *min_porosity;
  cs_lagr_clogging_param.diam_mean = *diam_mean;
  cs_lagr_clogging_param.valen = *valen;
  cs_lagr_clogging_param.phi_p = *phi_p;
  cs_lagr_clogging_param.phi_s = *phi_s;
  cs_lagr_clogging_param.cstham = *cstham;
  cs_lagr_clogging_param.csthpp = *csthpp;
  cs_lagr_clogging_param.lambda_vdw = *lambda_vdw;

  /* Allocate memory for the temperature and Debye length arrays */

  if (cs_lagr_clogging_param.temperature == NULL)
    BFT_MALLOC(cs_lagr_clogging_param.temperature, mesh->n_cells, cs_real_t);

  if (cs_lagr_clogging_param.debye_length == NULL)
    BFT_MALLOC(cs_lagr_clogging_param.debye_length, mesh->n_cells, cs_real_t);

  /* Store the temperature */

  for (iel = 0; iel < mesh->n_cells; iel++)
    cs_lagr_clogging_param.temperature[iel] = temperature[iel];

  /* Computation and storage of the Debye length                */

  for (iel = 0; iel < mesh->n_cells; iel++)

    cs_lagr_clogging_param.debye_length[iel]
      = pow(2e3 * pow(_faraday_cst,2)
            * cs_lagr_clogging_param.ionic_strength
            /  (  cs_lagr_clogging_param.water_permit
                * _free_space_permit * PG_CST
                * cs_lagr_clogging_param.temperature[iel]), -0.5);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  bft_printf(" epseau = %g\n", cs_lagr_clogging_param.water_permit);
  bft_printf(" valen   = %g\n", cs_lagr_clogging_param.valen);
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
 *   iel              <-- id of cell where the particle is
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
                         cs_lnum_t                       iel,
                         cs_real_t                      *energy_barrier,
                         cs_real_t                      *surface_coverage,
                         cs_real_t                      *limit,
                         cs_real_t                      *mporos)
{
  cs_real_t contact_area;
  cs_real_t deposited_radius;

  cs_real_t mean_nb_cont;

  cs_lnum_t  dim_aux = 1, contact_count[1];
  cs_real_t  value;
  cs_lnum_t  i;

  /* Computation of the number of particles in contact with */
  /* the depositing particle */

  /* Assuming monodispersed calculation */

  double p_diameter
    = cs_lagr_particle_get_real(particle, attr_map, CS_LAGR_DIAMETER);
  cs_real_t depositing_radius = p_diameter * 0.5;

  deposited_radius = depositing_radius;

  contact_area = _pi * pow(2. * pow(deposited_radius * depositing_radius, 0.5)
                           + deposited_radius,2);

  mean_nb_cont =   contact_area
                 * (*surface_coverage) / (_pi * pow(deposited_radius,2));

  /* Assuming Poisson distribution */

  value = 700.;
  if (mean_nb_cont > value) {
    cs_real_t  contact_count_r;
    cs_random_normal(1, &contact_count_r);
    contact_count[0] = (int) contact_count_r * pow(mean_nb_cont,0.5) + mean_nb_cont;
  }
  else {
    cs_random_poisson(dim_aux, mean_nb_cont, contact_count);
  }

  /* If the surface coverage is above the jamming limit,
     we are in multilayer deposition, so the contact number
     must be greater than zero  */

   /* The surface coverage must be greater than zero */
  if (*surface_coverage > cs_lagr_clogging_param.jamming_limit) {
    contact_count[0] +=1;
  }

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (mean_nb_cont > 0) {
    bft_printf("mean number = %g\n", mean_nb_cont);
    bft_printf("calculated number = %d\n",  contact_count[0]);
  }
#endif

  if (contact_count[0] == 0) {

    *energy_barrier = 0.0;

    /* Computation of the energy barrier */

    for (i = 0; i < 101; i++) {

      cs_real_t  step = cs_lagr_clogging_param.debye_length[iel]/30.0;

      cs_real_t distp = _d_cut_off + i*step;

      cs_real_t var1
        = cs_lagr_van_der_waals_sphere_plane(distp,
                                             depositing_radius,
                                             cs_lagr_clogging_param.lambda_vdw,
                                             cs_lagr_clogging_param.cstham);

      cs_real_t var2
        = cs_lagr_edl_sphere_plane(distp,
                                   depositing_radius,
                                   cs_lagr_clogging_param.valen,
                                   cs_lagr_clogging_param.phi_p,
                                   cs_lagr_clogging_param.phi_s,
                                   cs_lagr_clogging_param.temperature[iel],
                                   cs_lagr_clogging_param.debye_length[iel],
                                   cs_lagr_clogging_param.water_permit);

      cs_real_t var = var1 + var2;

      if (var > *energy_barrier)
        *energy_barrier = var;
      if (var < 0.)
        *energy_barrier = 0.;

    }

    *energy_barrier =  *energy_barrier / (0.5 * p_diameter);
  }

  else if (contact_count[0] > 0) {

    *energy_barrier = 0.0;

    /* Computation of the energy barrier */
    for (i = 0; i < 101; i++) {

      cs_real_t  step = cs_lagr_clogging_param.debye_length[iel]/30.0;

      cs_real_t distcc =   _d_cut_off + i*step
                         + depositing_radius + deposited_radius;

      cs_real_t var1
        = cs_lagr_van_der_waals_sphere_sphere(distcc,
                                              deposited_radius,
                                              depositing_radius,
                                              cs_lagr_clogging_param.lambda_vdw,
                                              cs_lagr_clogging_param.csthpp);

      cs_real_t var2
        = cs_lagr_edl_sphere_sphere(distcc,
                                    deposited_radius,
                                    depositing_radius,
                                    cs_lagr_clogging_param.valen,
                                    cs_lagr_clogging_param.phi_p,
                                    cs_lagr_clogging_param.phi_p,
                                    cs_lagr_clogging_param.temperature[iel],
                                    cs_lagr_clogging_param.debye_length[iel],
                                    cs_lagr_clogging_param.water_permit);

      cs_real_t var = contact_count[0] * (var1 + var2);

      if (var > *energy_barrier)
        *energy_barrier = var;
      if (var < 0.)
        *energy_barrier = 0.;

    }

    *energy_barrier =  *energy_barrier / (0.5 * p_diameter);
  }

  *limit = cs_lagr_clogging_param.jamming_limit;
  *mporos = cs_lagr_clogging_param.min_porosity;

  return contact_count[0];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
