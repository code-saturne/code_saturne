/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* VERS */

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
#include "cs_lagr_geom.h"

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
  const int itmx = 8;

  /* Data initializations with experimental measurements
     --------------------------------------------------- */

  unsigned char *particle = p_set->p_buffer + p_set->p_am->extents * ip;
  const cs_real_t *part_coords = cs_lagr_particle_attr_const(particle,
                                                             p_set->p_am,
                                                             CS_LAGR_COORDS);
  cs_real_t z = part_coords[2];

  /* transverse coordinate */
  cs_real_t  zi[] = {0.e-3 , 1.e-3 , 1.5e-3, 2.0e-3, 2.5e-3,
                     3.0e-3, 3.5e-3, 4.0e-3, 4.5e-3, 5.0e-3};

  /* vertical mean velocity of the particles */
  cs_real_t  ui[] = {5.544e0, 8.827e0, 9.068e0, 9.169e0, 8.923e0,
                     8.295e0, 7.151e0, 6.048e0, 4.785e0, 5.544e0};

  /* transverse mean velocity of the particles */
  cs_real_t  wi[] = { 0.e0   , 0.179e0, 0.206e0, 0.221e0, 0.220e0,
                      0.223e0, 0.206e0, 0.190e0, 0.195e0, 0.504e0};

  /* fluctuation of the vertical velocity of the particles */
  cs_real_t  uf[] = { 0.352e0, 0.352e0, 0.275e0, 0.252e0, 0.367e0,
                      0.516e0, 0.657e0, 0.872e0, 1.080e0, 0.792e0};

  /* fluctuation of the transverse velocity of the particles */
  cs_real_t  wf[] = { 0.058e0, 0.058e0, 0.056e0, 0.056e0, 0.060e0,
                      0.063e0, 0.058e0, 0.072e0, 0.091e0, 0.232e0};

#if 0
  /* shear-stress (currently not used) of the particle velocity */
  cs_real_t  uvi[] = {0.0017e0, 0.0017e0,  0.0016e0,  0.0027e0,  0.0077e0,
                      0.0146e0, 0.0206e0,  0.0447e0,  0.0752e0,  0.1145e0};
#endif

  /* Interpolation
     ------------- */

  int it = 0;

  if (z > zi[0]) {
    for (it = 0; it < itmx; it++) {
      if (z >= zi[it] && z < zi[it+1])
        break;
    }
  }

  /* Calculation of particles velocity
     --------------------------------- */

  cs_real_t up  =   ui[it] +(z - zi[it]) * (ui[it+1] - ui[it])
                  / (zi[it+1] - zi[it]);

  /* The value of the mean transverse velocity is currently set to zero
   * due to uncertainties on this variable */

  cs_real_t wp;
  if (false)
    wp  = wi[it] + (z - zi[it]) * (wi[it+1] - wi[it]) / (zi[it+1] - zi[it]);
  else
    wp = 0.0;

  cs_real_t upp = uf[it] +   (z - zi[it]) * (uf[it+1] - uf[it])
                           / (zi[it+1] - zi[it]);
  cs_real_t wpp = wf[it] +   (z - zi[it]) * (wf[it+1] - wf[it])
                           / (zi[it+1] - zi[it]);

  /* Calculations of the instantaneous particle velocity */

  cs_real_t vgauss[2];

  cs_random_normal(2, vgauss);

  cs_real_t *part_vel
    = cs_lagr_particle_attr(particle, p_set->p_am, CS_LAGR_VELOCITY);
  part_vel[0] = up  + vgauss[0] * upp;
  part_vel[1] = 0.0;
  part_vel[2] = wp + vgauss[1] * wpp;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of an external force field acting on the particles.
 *
 * It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
 * By default gravity and drag force are the only forces acting on the particles
 * (the gravity components gx gy gz are assigned in the GUI or in usipsu)
 *
 * \param[in]     dt_p     time step (for the cell)
 * \param[in]     taup     particle relaxation time
 * \param[in]     tlag     relaxation time for the flow
 * \param[in]     piil     term in the integration of the sde
 * \param[in]     bx       characteristics of the turbulence
 * \param[in]     tsfext   infos for the return coupling
 * \param[in]     vagaus   Gaussian random variables
 * \param[in]     gradpr   pressure gradient
 * \param[in]     gradvf   gradient of the flow velocity
 * \param[in,out] rho_p     particle density
 * \param[out]    fextla   user external force field (m/s^2)$
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_ef(cs_real_t            dt_p,
                const cs_real_t      taup[],
                const cs_real_3_t    tlag[],
                const cs_real_3_t    piil[],
                const cs_real_33_t   bx[],
                const cs_real_t      tsfext[],
                const cs_real_33_t   vagaus[],
                const cs_real_3_t    gradpr[],
                const cs_real_33_t   gradvf[],
                cs_real_t            rho_p[],
                cs_real_3_t          fextla[])
{
  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++){
    fextla[ip][0] = 0;
    fextla[ip][1] = 0;
    fextla[ip][2] = 0;
  }
}

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

  /* Example: computation of the particle mass flow rate on 4 planes
     --------------------------------------------------------------- */

  {
    cs_real_t zz[4] = {0.1e0, 0.15e0, 0.20e0, 0.25e0};

    /* If we are in an unsteady case, or if the beginning of the steady stats
     * is not reached yet, all statistics are reset to zero at each time
     step before entering this function.*/

    if(   cs_glob_lagr_time_scheme->isttio == 0
       || cs_glob_time_step->nt_cur <= cs_glob_lagr_stat_options->nstist) {
      for (cs_lnum_t iplan = 0; iplan < 4; iplan++)
        _m_flow[iplan] = 0.0;

    }

    for (cs_lnum_t iplan = 0; iplan < 4; iplan++) {

      for (cs_lnum_t npt = 0; p_set->n_particles; npt++) {

        unsigned char *part = p_set->p_buffer + p_am->extents * npt;

        cs_lnum_t iel = cs_lagr_particle_get_cell_id(part, p_am);

        if( iel >= 0 ) {

          const cs_real_t *part_coords
            = cs_lagr_particle_attr_const(part, p_am, CS_LAGR_COORDS);
          const cs_real_t *prev_part_coords
            = cs_lagr_particle_attr_n_const(part, p_am, 1, CS_LAGR_COORDS);

          if (    part_coords[0] > zz[iplan]
              && prev_part_coords[0] <= zz[iplan])
            _m_flow[iplan] +=  cs_lagr_particle_get_real(part, p_am,
                                                        CS_LAGR_STAT_WEIGHT)
                             * cs_lagr_particle_get_real(part, p_am,
                                                         CS_LAGR_MASS);

        }

      }

    }

    cs_real_t stat_age = cs_lagr_stat_get_age();

    for (cs_lnum_t iplan = 0; iplan < 4; iplan++)
      bft_printf(" Particle mass flow at Z(%d): %e14.5)",
                 iplan,
                 _m_flow[iplan]/stat_age);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Impose the motion of a particle flagged CS_LAGR_PART_IMPOSED_MOTION.
 *
 * User-defined modifications on the particle position and its
 * velocity.
 *
 * \param[in]   coords    old particle coordinates
 * \param[in]   dt        time step (per particle)
 * \param[out]  disp      particle dispacement
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_imposed_motion(const cs_real_t  coords[3],
                            cs_real_t        dt,
                            cs_real_t        disp[3])
{
  /* Angular velocity */
  cs_real_t omega = 1.0;

  /* Here we impose the particle to move arround a cylinder with
   * the axis  is (*, 0, 1) */
  cs_real_t rcost = (coords[1] - 0.0);
  cs_real_t rsint = (coords[2] - 1.0);

  /* Imposed displacement */
  disp[0] = 0.;
  disp[1] = rcost * (cos(omega*dt) - 1.0 ) - rsint * sin(omega*dt);
  disp[2] = rsint * (cos(omega*dt) - 1.0 ) + rcost * sin(omega*dt);
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

    /* velocity */

    cs_real_t *part_vel
      = cs_lagr_particles_attr(particles, p_id, CS_LAGR_VELOCITY);
    part_vel[0] = 1.0;
    part_vel[1] = 0.0;
    part_vel[2] = 0.0;

    /* diameter */

    cs_lagr_particles_set_real(particles, p_id, CS_LAGR_DIAMETER, 5e-05);

    /* Temperature profile */

    cs_lagr_particles_set_real(particles, p_id, CS_LAGR_TEMPERATURE, 20.0);

    /* Statistical weight profile */

    cs_lagr_particles_set_real(particles, p_id, CS_LAGR_STAT_WEIGHT, 0.01);

    /* User variables */
    for (int attr_id = CS_LAGR_USER;
         attr_id < CS_LAGR_USER + cs_glob_lagr_model->n_user_variables;
         attr_id++) {
      cs_real_t *user_var = cs_lagr_particle_attr(particles, p_id, attr_id);
      *user_var = 0.;
    }

  }

  if (particle_range[1] <= particle_range[0])
    return;

  /* Modifications occur after all the initializations related to
     the particle injection. */

  /* if new particles have entered the domain  */
  for (cs_lnum_t ip = particle_range[0]; ip < particle_range[1]; ip++) {
    _inlet2(particles, ip);
  }

  /*
   * Trick to average the statistics at iteration nstist
   * starting from an unsteady two-coupling calculation
   *                                                      */
  if (cs_glob_time_step->nt_cur > cs_glob_lagr_stat_options->nstist) {

    cs_glob_lagr_source_terms->nstits = cs_glob_lagr_stat_options->nstist;
    cs_glob_lagr_time_scheme->isttio = 1;

  }

  /* Simulation of the instantaneous turbulent fluid flow velocities seen
     by the solid particles along their trajectories
     -------------------------------------------------------------------- */

  /* In the previous operations, the particle data has been set with the
   * components of the instantaneous velocity (fluctuation + mean value) seen
   * by the particles.
   *
   * When the velocity of the flow is modified as above, most of the time
   * the user knows only the mean value. In some flow configurations and some
   * injection conditions, it may be necessary to reconstruct the
   * fluctuating part.
   * That is why the following function may be called.
   * Caution:
   *   - this turbulent component must be reconstructed only on the modified
   *     velocities of the flow seen.
   *   - the reconstruction must be adapted to the case. */

  if (false) {
    int time_id = 1;
    if (cs_glob_time_step->nt_cur == 1)
      time_id = 0;
    cs_lagr_new_particle_init(particle_range,
                              time_id,
                              visc_length);
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
 *                m   Cp
 *                 p    p
 *       Tau = ---------------
 *          c          2
 *                PI d    h
 *                    p    e
 *
 *      Tau  : Thermal relaxation time (value to be computed)
 *         c
 *
 *      m    : Particle mass
 *       p
 *
 *      Cp   : Particle specific heat
 *        p
 *
 *      d    : Particle diameter
 *       p
 *
 *      h    : Coefficient of thermal exchange
 *       e
 *
 *  The coefficient of thermal exchange is calculated from a Nusselt number,
 *  itself evaluated by a correlation (Ranz-Marshall by default)
 *
 *             h  d
 *              e  p
 *      Nu = --------  = 2 + 0.55 Re **(0.5) Prt**(0.33)
 *            Lambda                p
 *
 *      Lambda : Thermal conductivity of the carrier field
 *
 *      Re     : Particle Reynolds number
 *        p
 *
 *      Prt    : Prandtl number
 *
 * \param[in]   id_p   particle id
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
cs_user_lagr_rt(cs_lnum_t        id_p,
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

  unsigned char *particle = p_set->p_buffer + p_am->extents * id_p;
  cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);

  /*===============================================================================
   * Relaxation time with the standard (Wen-Yu) formulation of the drag coefficient
   *===============================================================================*/

  /* This example gives the standard relaxation time as an indication:*/

  cs_real_t fdr;

  cs_real_t cd1  = 0.15;
  cs_real_t cd2  = 0.687;

  if (re_p <= 1000)
    fdr = 18.0 * nu_f * (1.0 + cd1 * pow(re_p, cd2)) / (p_diam * p_diam);

  else
    fdr = (0.44 * 3.0 / 4.0) * uvwr / p_diam;

  taup[id_p] = rho_p / rho_f / fdr;

  /*===============================================================================
   * Computation of the relaxation time with the drag coefficient of
   * S.A. Morsi and A.J. Alexander, J. of Fluid Mech., Vol.55, pp 193-208 (1972)
   *===============================================================================*/

  cs_real_t rec1 =  0.1;
  cs_real_t rec2 =  1.0;
  cs_real_t rec3 =  10.0;
  cs_real_t rec4 = 200.0;

  cs_real_t dd2 = p_diam * p_diam;

  if (re_p <= rec1)
    fdr = 18.0 * nu_f / dd2;

  else if (re_p <= rec2)
    fdr = 3.0/4.0 * nu_f / dd2 * (22.73 + 0.0903 / re_p + 3.69 * re_p);

  else if (re_p <= rec3)
    fdr = 3.0/4.0 * nu_f / dd2 * (29.1667 - 3.8889 / re_p + 1.222 * re_p);

  else if (re_p <=rec4)
    fdr = 18.0 * nu_f / dd2 *(1.0 + 0.15 * pow(re_p, 0.687));

  else
    fdr = (0.44 * 3.0 / 4.0) * uvwr / p_diam;

  taup[id_p] = rho_p / rho_f / fdr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the computation of the thermal relaxation time
 *        of the particles with respect to the chosen formulation of
 *        the Nusselt number.
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 * \param[in]   id_p   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[in]   cp_f   specific heat of the fluid at particle position
 * \param[in]   k_f    diffusion coefficient of the fluid at particle position
 * \param[out]  tauc   thermal relaxation time
 * \param[in]   dt     time step (per cell)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_rt_t(cs_lnum_t        id_p,
                  cs_real_t        re_p,
                  cs_real_t        uvwr,
                  cs_real_t        rho_f,
                  cs_real_t        rho_p,
                  cs_real_t        nu_f,
                  cs_real_t        cp_f,
                  cs_real_t        k_f,
                  cs_real_t        tauc[],
                  const cs_real_t  dt[])
{
  /* 1. Initializations: Particles management */
  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  unsigned char *particle = p_set->p_buffer + p_am->extents * id_p;

  /* 2. Standard thermal relaxation time */

  /* This example gives the standard thermal relaxation time
   * as an indication.*/

  cs_real_t prt = nu_f / k_f;

  cs_real_t fnus = 2.0 + 0.55 * sqrt(re_p) * pow(prt, 1./3.);

  cs_real_t diam = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER);
  cs_real_t cp_p = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CP);

  tauc[id_p]= diam * diam * rho_p * cp_p  / ( fnus * 6.0 * rho_f * cp_f * k_f);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User integration of the SDE for the user-defined variables.
 *
 * The variables are constant by default. The SDE must be of the form:
 *
 * \f[
 *    \frac{dT}{dt}=\frac{T - PIP}{Tca}
 * \f]
 *
 * T:   particle attribute representing the variable
 * Tca: characteristic time for the sde
 *      to be prescribed in the array auxl1
 * PIP: coefficient of the SDE (pseudo RHS)
 *      to be prescribed in the array auxl2.
 *      If the chosen scheme is first order (nordre=1) then, at the first
 *      and only call pip is expressed as a function of the quantities of
 *      the previous time step (contained in the particle data).
 *      If the chosen scheme is second order (nordre=2)
 *      then, at the first call (nor=1) pip is expressed as a function of
 *      the quantities of the previous time step, and at the second passage
 *      (nor=2) pip is expressed as a function of the quantities of the
 *      current time step.
 *
 * \param[in]  dt      time step (per cell)
 * \param[in]  taup    particle relaxation time
 * \param[in]  tlag    relaxation time for the flow
 * \param[in]  tempct  characteristic thermal time and implicit source
 *                     term of return coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_sde(const cs_real_t  dt[],
                 cs_real_t        taup[],
                 cs_real_3_t      tlag[],
                 cs_real_t        tempct[])
{
  /* Initializations
     --------------- */

  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  cs_real_t *tcarac, *pip;

  BFT_MALLOC(tcarac, p_set->n_particles, cs_real_t);
  BFT_MALLOC(pip   , p_set->n_particles, cs_real_t);

  /* Characteristic time of the current SDE
     -------------------------------------- */

  /* Loop on the additional variables */

  for (int i = 0;
       i < cs_glob_lagr_model->n_user_variables;
       i++) {

    for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

      unsigned char *part = p_set->p_buffer + p_am->extents * npt;
      cs_lnum_t iel = cs_lagr_particle_get_cell_id(part, p_am);

      cs_real_t *usr_var
        = cs_lagr_particle_attr_n(part, p_am, 0, CS_LAGR_USER);
      cs_real_t *prev_usr_var
        = cs_lagr_particle_attr_n(part, p_am, 1, CS_LAGR_USER);

      if (iel >= 0) {

        /* Characteristic time tca of the differential equation,
           This example must be adapted to the case */
        tcarac[npt] = 1.0;

        /* Prediction at the first substep;
           This example must be adapted to the case */
        if (cs_glob_lagr_time_step->nor == 1)
          pip[npt] = prev_usr_var[i];

        /* Correction at the second substep;
           This example must be adapted to the case */
        else
          pip[npt] = usr_var[i];

      }

    }

    /* Integration of the variable ipl
       ------------------------------- */

    cs_lagr_sde_attr(CS_LAGR_USER, tcarac, pip);

  }

  BFT_FREE(tcarac);
  BFT_FREE(pip);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
