/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* code_saturne version 6.2-alpha */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

static cs_real_t _m_flow[4];

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function (non-mandatory intervention)
 *
 * User-defined modifications on the variables at the end of the
 * Lagrangian time step and calculation of user-defined
 * additional statistics on the particles.
 *
 * \param[in]  dt      time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_extra_operations(const cs_real_t  dt[])
{

  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  /* Print results
     -------------*/

  cs_time_step_t *ts = cs_get_glob_time_step();

  char file_name1[20];
  sprintf(file_name1,"results_to_plot");

  if (ts->nt_cur == 1) {
    FILE * file1 = fopen(file_name1,"w+");
    fprintf(file1,"# Time, Omega_x, Omega_y, Omega_z, Gamma_x, Gamma_y, Gamma_z, gamma_x, gamma_y, gamma_z, upus_x, upus_y, upus_z, <xp2>_x, <xp2>_y, <xp2>_z, <Up2>_x, <Up2>_y, <Up2>_z, <Us2>_x, <Us2>_y, <Us2>_z, <UpUs>_x, <UpUs>_y, <UpUs>_z\n");
    fclose(file1);
  }

  /* Write results of the variances in the file */
  FILE * file1 = fopen(file_name1,"a");
  {
    /* Computing X-components */
    
    cs_real_t Bx_p = cs_notebook_parameter_value_by_name("bx");
    cs_real_t taup = cs_notebook_parameter_value_by_name("taup");
    cs_real_t dtot = cs_glob_lagr_time_step->ttclag;

    cs_real_t tlag_x = cs_notebook_parameter_value_by_name("tlag_x");
    cs_real_t tlag_y = cs_notebook_parameter_value_by_name("tlag_y");
    cs_real_t tlag_z = cs_notebook_parameter_value_by_name("tlag_z");
    

    cs_real_t Omega_xp_x = pow(Bx_p * tlag_x/(tlag_x-taup), 2) *
      ( pow(tlag_x - taup, 2) * dtot +
        pow(tlag_x, 3) * 0.5 * (1. - exp(-2.*dtot/tlag_x) ) +
        pow(taup, 3) * 0.5 * (1. - exp(-2.*dtot/taup) ) -
        2. * pow(tlag_x, 2) * (tlag_x-taup) * (1. - exp(-dtot/tlag_x) ) +
        2. * pow(taup, 2) * (tlag_x-taup) * (1. - exp(-dtot/taup) ) -
        2. * pow(tlag_x * taup, 2) / (tlag_x+taup) * ( 1. - exp(-dtot/tlag_x)*exp(-dtot/taup) )
        );
        
    cs_real_t Omega_xp_y = pow(Bx_p * tlag_y/(tlag_y-taup), 2) *
      ( pow(tlag_y - taup, 2) * dtot +
        pow(tlag_y, 3) * 0.5 * (1. - exp(-2.*dtot/tlag_y) ) +
        pow(taup, 3) * 0.5 * (1. - exp(-2.*dtot/taup) ) -
        2. * pow(tlag_y, 2) * (tlag_y-taup) * (1. - exp(-dtot/tlag_y) ) +
        2. * pow(taup, 2) * (tlag_y-taup) * (1. - exp(-dtot/taup) ) -
        2. * pow(tlag_y * taup, 2) / (tlag_y+taup) * ( 1. - exp(-dtot/tlag_y)*exp(-dtot/taup) )
        );

    cs_real_t Omega_xp_z = pow(Bx_p * tlag_z/(tlag_z-taup), 2) *
      ( pow(tlag_z - taup, 2) * dtot +
        pow(tlag_z, 3) * 0.5 * (1. - exp(-2.*dtot/tlag_z) ) +
        pow(taup, 3) * 0.5 * (1. - exp(-2.*dtot/taup) ) -
        2. * pow(tlag_z, 2) * (tlag_z-taup) * (1. - exp(-dtot/tlag_z) ) +
        2. * pow(taup, 2) * (tlag_z-taup) * (1. - exp(-dtot/taup) ) -
        2. * pow(tlag_z * taup, 2) / (tlag_z+taup) * ( 1. - exp(-dtot/tlag_z)*exp(-dtot/taup) )
        );

    cs_real_t Gamma_up_x = pow(Bx_p * tlag_x/(tlag_x-taup), 2) *
      ( 0.5 * tlag_x * (1. - exp(-2. * dtot/tlag_x) ) +
        0.5 * taup * (1. - exp(-2. * dtot/taup) ) -
        2. * tlag_x * taup / (tlag_x+taup) * (1. - exp(-dtot/tlag_x)*exp(-dtot/taup) )
        );
        
    cs_real_t Gamma_up_y = pow(Bx_p * tlag_y/(tlag_y-taup), 2) *
      ( 0.5 * tlag_y * (1. - exp(-2. * dtot/tlag_y) ) +
        0.5 * taup * (1. - exp(-2. * dtot/taup) ) -
        2. * tlag_y * taup / (tlag_y+taup) * (1. - exp(-dtot/tlag_y)*exp(-dtot/taup) )
        );

    cs_real_t Gamma_up_z = pow(Bx_p * tlag_z/(tlag_z-taup), 2) *
      ( 0.5 * tlag_z * (1. - exp(-2. * dtot/tlag_z) ) +
        0.5 * taup * (1. - exp(-2. * dtot/taup) ) -
        2. * tlag_z * taup / (tlag_z+taup) * (1. - exp(-dtot/tlag_z)*exp(-dtot/taup) )
        );

    cs_real_t gamma_us_x =  pow(Bx_p ,2) *
      0.5 * tlag_x * (1. - exp(-2. * dtot/tlag_x) );

    cs_real_t gamma_us_y =  pow(Bx_p ,2) *
      0.5 * tlag_y * (1. - exp(-2. * dtot/tlag_y) );

    cs_real_t gamma_us_z =  pow(Bx_p ,2) *
      0.5 * tlag_z * (1. - exp(-2. * dtot/tlag_z) );

    cs_real_t gamma_usup_x =  pow(Bx_p ,2) * tlag_x / (tlag_x - taup) * tlag_x *
      ( 0.5 * (1. - exp(-2.*dtot/tlag_x))
        - taup/(taup+tlag_x) * (1. - exp(-dtot/tlag_x)*exp(-dtot/taup)) );

    cs_real_t gamma_usup_y =  pow(Bx_p ,2) * tlag_y / (tlag_y - taup) * tlag_y *
      ( 0.5 * (1. - exp(-2.*dtot/tlag_y))
        - taup/(taup+tlag_y) * (1. - exp(-dtot/tlag_y)*exp(-dtot/taup)) );

    cs_real_t gamma_usup_z =  pow(Bx_p ,2) * tlag_z / (tlag_z - taup) * tlag_z *
      ( 0.5 * (1. - exp(-2.*dtot/tlag_z))
        - taup/(taup+tlag_z) * (1. - exp(-dtot/tlag_z)*exp(-dtot/taup)) );

    cs_real_t var_xpxp_x = 0.0;
    cs_real_t var_xpxp_y = 0.0;
    cs_real_t var_xpxp_z = 0.0;
    
    cs_real_t var_UpUp_x = 0.0;
    cs_real_t var_UpUp_y = 0.0;
    cs_real_t var_UpUp_z = 0.0;
    
    cs_real_t var_UsUs_x = 0.0;
    cs_real_t var_UsUs_y = 0.0;
    cs_real_t var_UsUs_z = 0.0;
    
    cs_real_t var_UpUs_x = 0.0;
    cs_real_t var_UpUs_y = 0.0;
    cs_real_t var_UpUs_z = 0.0;

    for (cs_lnum_t p_id = 0; p_id < p_set->n_particles; p_id++) {

      cs_real_t *part_vel      = cs_lagr_particles_attr(p_set, p_id,
                                                       CS_LAGR_VELOCITY);
      cs_real_t *part_vel_seen = cs_lagr_particles_attr(p_set, p_id,
                                                       CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_coords   = cs_lagr_particles_attr(p_set, p_id,
                                                       CS_LAGR_COORDS);

      var_xpxp_x +=
        pow(part_coords[0], 2);
      var_xpxp_y +=
        pow(part_coords[1], 2);
      var_xpxp_z +=
        pow(part_coords[2], 2);

      var_UpUp_x +=
        pow(part_vel[0], 2);
      var_UpUp_y +=
        pow(part_vel[1], 2);
      var_UpUp_z +=
        pow(part_vel[2], 2);

      var_UsUs_x +=
        pow(part_vel_seen[0], 2);
      var_UsUs_y +=
        pow(part_vel_seen[1], 2);
      var_UsUs_z +=
        pow(part_vel_seen[2], 2);

      var_UpUs_x +=
        part_vel_seen[0] * part_vel[0];
      var_UpUs_y +=
        part_vel_seen[1] * part_vel[1];
      var_UpUs_z +=
        part_vel_seen[2] * part_vel[2];

    }

    fprintf(file1," %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
           dtot,
           Omega_xp_x,
           Omega_xp_y,
           Omega_xp_z,
           Gamma_up_x,
           Gamma_up_y,
           Gamma_up_z,
           gamma_us_x,
           gamma_us_y,
           gamma_us_z,
           gamma_usup_x,
           gamma_usup_y,
           gamma_usup_z,
           var_xpxp_x / p_set->n_particles,
           var_xpxp_y / p_set->n_particles,
           var_xpxp_z / p_set->n_particles,
           var_UpUp_x / p_set->n_particles,
           var_UpUp_y / p_set->n_particles,
           var_UpUp_z / p_set->n_particles,
           var_UsUs_x / p_set->n_particles,
           var_UsUs_y / p_set->n_particles,
           var_UsUs_z / p_set->n_particles,
           var_UpUs_x / p_set->n_particles,
           var_UpUs_y / p_set->n_particles,
           var_UpUs_z / p_set->n_particles);
  }
  fclose(file1);

}

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

  const int ntcabs = cs_glob_time_step->nt_cur;

  cs_lagr_zone_data_t  *lagr_bdy_conditions
    = cs_lagr_get_boundary_conditions();

  /* Simple changes to selected attributes
     ------------------------------------- */

  for (cs_lnum_t p_id = particle_range[0]; p_id < particle_range[1]; p_id++) {

    /* Particle position */
    cs_real_t *part_coord
      = cs_lagr_particles_attr(particles, p_id, CS_LAGR_COORDS);

    /* Particle velocity  */
    cs_real_t *part_vel
      = cs_lagr_particles_attr(particles, p_id, CS_LAGR_VELOCITY);

    /* Seen velocity  */
    cs_real_t *vel_seen
      = cs_lagr_particles_attr(particles, p_id, CS_LAGR_VELOCITY_SEEN);

    for (cs_lnum_t i = 0; i < 3; i++) {
      part_coord[i] = 0.;
      part_vel[i] = 0.;
      vel_seen[i] = 0.;
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the calculation of the particle relaxation time
 *  with respect to the chosen formulation for the drag coefficient
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 *
 *      \tau_c = \frac{m_p{C_p}_p}{PId_p^2h_e}
 *
 *      \tau_c  : Thermal relaxation time (value to be computed)
 *
 *      m_p    : Particle mass
 *
 *      {C_p}_p   : Particle specific heat
 *
 *      d_p    : Particle diameter
 *
 *      h_e    : Coefficient of thermal exchange
 *
 *  The coefficient of thermal exchange is calculated from a Nusselt number,
 *  itself evaluated by a correlation (Ranz-Marshall by default)
 *
 *      \nu =  \frac{h_ed_p}{\lambda} = 2 + 0.55{\Re_e}_p^{0.5}P_{rt}^{0.33}
 *
 *      \lambda : Thermal conductivity of the carrier field
 *
 *      {\Re_e}_p     : Particle Reynolds number
 *
 *      P_{rt}    : Prandtl number
 *
 * \param[in]   p_id   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[in]   cp_f   specific heat of the fluid at particle position
 * \param[in]   k_f    diffusion coefficient of the fluid at particle position
 * \param[out]  taup   thermal relaxation time
 * \param[in]   dt     time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_rt(cs_lnum_t        p_id,
                cs_real_t        re_p,
                cs_real_t        uvwr,
                cs_real_t        rho_f,
                cs_real_t        rho_p,
                cs_real_t        nu_f,
                cs_real_t        taup[],
                const cs_real_t  dt[])
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  unsigned char *particle = p_set->p_buffer + p_am->extents * p_id;

  taup[p_id] = cs_notebook_parameter_value_by_name("taup");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
