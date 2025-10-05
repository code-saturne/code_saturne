/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User definition of an external force field acting on the particles.
 *
 * It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
 * By default gravity and drag force are the only forces acting on the particles
 * (the gravity components gx gy gz are assigned in the GUI or in
 * cs_user_parameters)
 * Note that taup, tlag, piil and bx are associated to the particle and
 * their value for all phases is provided
 *
 * \param[in]      dt_p     time step (for the cell)
 * \param[in]      p_id     particle id
 * \param[in]      taup     particle relaxation time
 * \param[in]      tlag     relaxation time for the flow
 * \param[in]      piil     term in the integration of the sde
 * \param[in]      bx       characteristics of the turbulence
 * \param[in]      tsfext   infos for the return coupling
 * \param[in]      vagaus   Gaussian random variables
 * \param[in,out]  rho_p     particle density
 * \param[out]     fextla   user external force field (m/s^2)$
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_ef
void
cs_user_lagr_ef([[maybe_unused]] cs_real_t            dt_p,
                [[maybe_unused]] const cs_lnum_t      p_id,
                [[maybe_unused]] const cs_real_t     *taup,
                [[maybe_unused]] const cs_real_3_t   *tlag,
                [[maybe_unused]] const cs_real_3_t   *piil,
                [[maybe_unused]] const cs_real_33_t  *bx,
                [[maybe_unused]] const cs_real_t      tsfext,
                [[maybe_unused]] const cs_real_3_t   *vagaus,
                [[maybe_unused]] const cs_real_3_t    gradpr,
                [[maybe_unused]] const cs_real_33_t   gradvf,
                [[maybe_unused]] cs_real_t            rho_p,
                [[maybe_unused]] cs_real_3_t          fextla)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User function (non-mandatory intervention)
 *
 * User-defined modifications on the variables at the end of the
 * Lagrangian time step and calculation of user-defined
 * additional statistics on the particles.
 *
 * \param[in]  dt      time step (per cell)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_extra_operations
void
cs_user_lagr_extra_operations([[maybe_unused]] const cs_real_t  dt[])
{
}

/*----------------------------------------------------------------------------*/
/*
 * Impose the motion of a particle flagged CS_LAGR_PART_IMPOSED_MOTION.
 *
 * User-defined modifications on the particle position and its
 * velocity.
 *
 * \param[in]   particles       pointer to particle set
 * \param[in]   p_id            particle id
 * \param[in]   coords          old particle coordinates
 * \param[in]   dt              time step (per particle)
 * \param[out]  disp            particle dispacement
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_imposed_motion
void
cs_user_lagr_imposed_motion
(
  [[maybe_unused]] const cs_lagr_particle_set_t *particles,
  [[maybe_unused]] cs_lnum_t                     p_id,
  [[maybe_unused]] const cs_real_t               coords[3],
  [[maybe_unused]] const cs_real_t               dt,
  [[maybe_unused]] cs_real_t                     disp[3]
)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User modification of newly injected particle location.
 *
 * This function aims at modifying injection coordinates, particle properties
 * and cell_id depending on the position are updated based on the modified
 * position after this function and before cs_user_lagr_in.
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modified by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               zone injection set data
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  null otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_in_force_coords
void
cs_user_lagr_in_force_coords
(
  [[maybe_unused]] cs_lagr_particle_set_t         *particles,
  [[maybe_unused]] const cs_lagr_injection_set_t  *zis,
  [[maybe_unused]] const cs_lnum_t                 particle_range[2],
  [[maybe_unused]] const cs_lnum_t                 particle_face_id[],
  [[maybe_unused]] const cs_real_t                 visc_length[]
)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User modification of newly injected particles.
 *
 * This function is called after the initialization of the new particles in
 * order to modify them according to new particle profiles (injection
 * profiles, statistical weights, correction of the diameter if the
 * standard-deviation option is activated); the modification of particles
 * position should preferentially be made in cs_user_lagr_in_force_coords to
 * get an initialization of particle properties coherent with the local fields.
 *
 * This function is called for each injection zone and set. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modified by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               zone injection set data
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  null otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_in
void
cs_user_lagr_in
(
  [[maybe_unused]] cs_lagr_particle_set_t         *particles,
  [[maybe_unused]] const cs_lagr_injection_set_t  *zis,
  [[maybe_unused]] const cs_lnum_t                 particle_range[2],
  [[maybe_unused]] const cs_lnum_t                 particle_face_id[],
  [[maybe_unused]] const cs_real_t                 visc_length[]
)
{
}

/*----------------------------------------------------------------------------*/
/*
 *  Modification of the calculation of the particle relaxation time
 *  with respect to the chosen formulation for the drag coefficient
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
 *
 *      \f$\tau_c = \frac{m_p{C_p}_p}{PId_p^2h_e}\f$
 *
 *      \f$\tau_c\f$  : Thermal relaxation time (value to be computed)
 *
 *      \f$m_p\f$     : Particle mass
 *
 *      \f${C_p}_p\f$ : Particle specific heat
 *
 *      \f$d_p\f$     : Particle diameter
 *
 *      \f$h_e\f$     : Coefficient of thermal exchange
 *
 *  The coefficient of thermal exchange is calculated from a Nusselt number,
 *  itself evaluated by a correlation (Ranz-Marshall by default)
 *
 *      \f$\nu = \frac{h_ed_p}{\lambda} = 2 + 0.55{\Re_e}_p^{0.5}P_{rt}^{0.33}\f$
 *
 *      \f$\lambda\f$ : Thermal conductivity of the carrier field
 *
 *      \f${\Re_e}_p\f$     : Particle Reynolds number
 *
 *     \f$ P_{rt}\f$    : Prandtl number
 *
 * \param[in]   phase_id   carrier phase_id
 * \param[in]   p_id       particle id
 * \param[in]   re_p       particle Reynolds number
 * \param[in]   uvwr       relative velocity of the particle
 *                         (flow-seen velocity - part. velocity)
 * \param[in]   rho_f      fluid density at  particle position
 * \param[in]   rho_p      particle density
 * \param[in]   nu_f       kinematic viscosity of the fluid at particle position
 * \param[out]  taup       thermal relaxation time
 * \param[in]   dt         time step associated to the particle
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_rt
void
cs_user_lagr_rt
(
  [[maybe_unused]] int              phase_id,
  [[maybe_unused]] cs_lnum_t        p_id,
  [[maybe_unused]] cs_real_t        re_p,
  [[maybe_unused]] cs_real_t        uvwr,
  [[maybe_unused]] cs_real_t        rho_f,
  [[maybe_unused]] cs_real_t        rho_p,
  [[maybe_unused]] cs_real_t        nu_f,
  [[maybe_unused]] cs_real_t       *taup,
  [[maybe_unused]] const cs_real_t  dt
)
{
}

/*----------------------------------------------------------------------------*/
/*
 * Modification of the computation of the thermal relaxation time
 * of the particles with respect to the chosen formulation of
 * the Nusselt number.
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
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
 * \param[out]  tauc   thermal relaxation time
 * \param[in]   dt     time step associated to the particle
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_rt_t
void
cs_user_lagr_rt_t
(
  [[maybe_unused]] cs_lnum_t        p_id,
  [[maybe_unused]] cs_real_t        re_p,
  [[maybe_unused]] cs_real_t        uvwr,
  [[maybe_unused]] cs_real_t        rho_f,
  [[maybe_unused]] cs_real_t        rho_p,
  [[maybe_unused]] cs_real_t        nu_f,
  [[maybe_unused]] cs_real_t        cp_f,
  [[maybe_unused]] cs_real_t        k_f,
  [[maybe_unused]] cs_real_2_t      tempct,
  [[maybe_unused]] const cs_real_t  dt
)
{
}

/*----------------------------------------------------------------------------*/
/*
 * User integration of the SDE for the user-defined variables.
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
 * Note that taup, tlag are associated to the particle and their value for
 * all phases is provided
 *
 * \param[in]  dt      time step (per cell)
 * \param[in]  p_id    particle id
 * \param[in]  taup    particle relaxation time
 * \param[in]  tlag    relaxation time for the flow
 * \param[in]  tempct  characteristic thermal time and implicit source
 *                     term of return coupling
 * \param[in]  nor     current step id (for 2nd order scheme)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_sde
void
cs_user_lagr_sde
(
  [[maybe_unused]] const cs_real_t         dt,
  [[maybe_unused]] const cs_lnum_t         p_id,
  [[maybe_unused]] const cs_real_t        *taup,
  [[maybe_unused]] const cs_real_3_t      *tlag,
  [[maybe_unused]] const cs_real_2_t       tempct,
  [[maybe_unused]] const int               nor
)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
