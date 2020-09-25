/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* Code_Saturne version 6.0.2 */

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

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
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
#include "cs_halo.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_notebook.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"
#include "cs_thermal_model.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr.h"
#include "cs_lagr_new.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_sde.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Global variables
 *============================================================================*/

static cs_real_t _m_flow[4];

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define inlet conditions based on experimental data for a given particle
 *
 * parameters:
 *   p_set  <-> particle
 *   ip     <-- particle id
 *----------------------------------------------------------------------------*/

static void
_inlet2(cs_lagr_particle_set_t  *p_set,
        cs_lnum_t                ip)
{
  const int itmx = 4;

  /* Data initializations with experimental measurements
     --------------------------------------------------- */

  unsigned char *particle = p_set->p_buffer + p_set->p_am->extents * ip;
  const cs_real_t *part_coords = cs_lagr_particle_attr_const(particle,
                                                             p_set->p_am,
                                                             CS_LAGR_COORDS);

  cs_real_t xx = part_coords[0];
  cs_real_t yy = part_coords[1];
  cs_real_t ray = sqrt(xx*xx+yy*yy);
  cs_real_t teta = 0.0;
  
/*  printf("%f\n", xx);
  printf("%f\n", yy);
  printf("ray");
  printf("%f\n", ray);*/
  
  if(xx != 0.0 && yy != 0.0)
  {
    teta = atan(yy/xx);
  }
  else if (xx == 0.0)
  {
    teta = 0.0;
  }
  else
  {
    teta = 3.14159/2.0;
  }

  /* radial distance M */
  cs_real_t  rr[] = {0.e-3, 2.e-3 , 4.e-3, 6.0e-3, 8.e-3, 10.e-3};

  /* axial velocity M/S */
  cs_real_t  va[] = {3.989888, 3.955918, 3.818365, 3.615492, 3.235091, 2.148246};

  /* radial velocity M/S */
  cs_real_t  vr[] = {0.021573, 0.019587, 0.019053, 0.018448, 0.015862, 0.098097};

  /* tangential velocity M/S */
  cs_real_t  vt[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  /* ECART TYPE axial fluctuation velocity M/S */
  cs_real_t  vaf[] = {0.492889, 0.248873, 0.276544, 0.339334, 0.45922, 0.60341};

  /* ECART TYPE radial fluctuation velocity M/S */
  cs_real_t  vrf[] = {0.145902, 0.153009, 0.161393, 0.174045, 0.232101, 0.198226};

  /* ECART TYPE tangential fluctuation velocity M/S */
  cs_real_t  vtf[] = {0.145902, 0.153009, 0.161393, 0.174045, 0.232101, 0.198226};

  /* Point localisation
     ------------------ */

  bool TROUVE = false;
  int it = 0;
  int itt = 0;

  if (ray < rr[itmx] && ray > rr[0])
  {
    for (it = 0; it <= itmx; it++)
    {
      if(ray >= rr[it] && ray < rr[it+1])
      {
        if(TROUVE)
        {
          printf("error search");
          exit(0);
        }
        else
        {
          itt = it;
          TROUVE = true;
        }
      }
    }
  }
  else if (ray >= rr[itmx])
  {
    itt = itmx;
  }
  else if (ray <= rr[0])
  {
    itt = 0;
  }
  else
  {
    printf("error interpolation : inlet2");
    exit(0);
  }


  /* Interpolation
     ------------- */

  cs_real_t vap = 0.0;
  cs_real_t vrp = 0.0;
  cs_real_t vapf = 0.0;
  cs_real_t vrpf = 0.0;
  cs_real_t vtpf = 0.0;

  if(ray <= rr[itmx+1])
  {
    vap  = va[itt] +(ray-rr[itt])*(va[itt+1] -va[itt]) /(rr[itt+1]-rr[itt]);
    vrp  = vr[itt] +(ray-rr[itt])*(vr[itt+1] -vr[itt]) /(rr[itt+1]-rr[itt]);
    vapf = vaf[itt]+(ray-rr[itt])*(vaf[itt+1]-vaf[itt])/(rr[itt+1]-rr[itt]);
    vrpf = vrf[itt]+(ray-rr[itt])*(vrf[itt+1]-vrf[itt])/(rr[itt+1]-rr[itt]);
    vtpf = vtf[itt]+(ray-rr[itt])*(vtf[itt+1]-vtf[itt])/(rr[itt+1]-rr[itt]);
  }
  else
  {
    vap  = va[itmx+1];
    vrp  = vr[itmx+1];
    vapf = vaf[itmx+1];
    vrpf = vrf[itmx+1];
    vtpf = vtf[itmx+1];
  }

  /* Turbulence
     ---------- */

  cs_real_t vgauss[3];

  cs_random_normal(3, vgauss);


  /* Calculations of the instantaneous particle velocity */

  cs_real_t vy1 = vrp  + vgauss[1] * vrpf;
  cs_real_t vy2 = 0.0  + vgauss[2] * vtpf;

  cs_real_t *part_vel
    = cs_lagr_particle_attr(particle, p_set->p_am, CS_LAGR_VELOCITY);

  part_vel[0] = vy1*cos(teta)-vy2*sin(teta);
  part_vel[1] = vy1*sin(teta)+vy2*cos(teta);
  part_vel[2] = vap + vgauss[0] * vapf;

/*  printf("part_vel");
  printf("%f\n", part_vel[0]);
  printf("%f\n", part_vel[1]);
  printf("%f\n", part_vel[2]);
*/
  
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * User function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particles.
 *
 * This function is called after the initialization of the new particles in
 * order to modify them according to new particle profiles (injection
 * profiles, position of the injection point, statistical weights,
 * correction of the diameter if the standard-deviation option is activated).
 *
 * This function is called for each injection zone and set. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modidied by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               zone injection set data
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  NULL otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_in(cs_lagr_particle_set_t         *particles,
                const cs_lagr_injection_set_t  *zis,
                const cs_lnum_t                 particle_range[2],
                const cs_lnum_t                 particle_face_id[],
                const cs_real_t                 visc_length[])
{
  /* Loop on new particles
     --------------------- */

  /* Modifications occur after all the initializations related to
     the particle injection. */

  /* if new particles have entered the domain  */
  for (cs_lnum_t ip = particle_range[0]; ip < particle_range[1]; ip++) {
    _inlet2(particles, ip);
  }
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
