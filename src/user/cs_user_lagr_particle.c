/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#pragma weak cs_user_lagr_ef
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
  CS_UNUSED(dt_p);
  CS_UNUSED(taup);
  CS_UNUSED(tlag);
  CS_UNUSED(piil);
  CS_UNUSED(bx);
  CS_UNUSED(tsfext);
  CS_UNUSED(vagaus);
  CS_UNUSED(gradpr);
  CS_UNUSED(gradvf);
  CS_UNUSED(rho_p);
  CS_UNUSED(fextla);
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

#pragma weak cs_user_lagr_extra_operations
void
cs_user_lagr_extra_operations(const cs_real_t  dt[])
{
  CS_UNUSED(dt);
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

#pragma weak cs_user_lagr_imposed_motion
void
cs_user_lagr_imposed_motion(const cs_real_t  coords[3],
                            cs_real_t        dt,
                            cs_real_t        disp[3])
{
  CS_UNUSED(coords);
  CS_UNUSED(dt);
  CS_UNUSED(disp);
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
 * but may be modified by this function.
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

#pragma weak cs_user_lagr_in
void
cs_user_lagr_in(cs_lagr_particle_set_t         *particles,
                const cs_lagr_injection_set_t  *zis,
                const cs_lnum_t                 particle_range[2],
                const cs_lnum_t                 particle_face_id[],
                const cs_real_t                 visc_length[])
{
  CS_UNUSED(particles);
  CS_UNUSED(zis);
  CS_UNUSED(particle_range);
  CS_UNUSED(particle_face_id);
  CS_UNUSED(visc_length);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the calculation of the particle relaxation time
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
 * \param[in]   id_p   particle id
 * \param[in]   re_p   particle Reynolds number
 * \param[in]   uvwr   relative velocity of the particle
 *                     (flow-seen velocity - part. velocity)
 * \param[in]   rho_f  fluid density at  particle position
 * \param[in]   rho_p  particle density
 * \param[in]   nu_f   kinematic viscosity of the fluid at particle position
 * \param[out]  taup   thermal relaxation time
 * \param[in]   dt     time step (per cell)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_lagr_rt
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
  CS_UNUSED(id_p);
  CS_UNUSED(re_p);
  CS_UNUSED(uvwr);
  CS_UNUSED(rho_p);
  CS_UNUSED(rho_f);
  CS_UNUSED(nu_f);
  CS_UNUSED(taup);
  CS_UNUSED(dt);
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

#pragma weak cs_user_lagr_rt_t
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
  CS_UNUSED(id_p);
  CS_UNUSED(re_p);
  CS_UNUSED(uvwr);
  CS_UNUSED(rho_p);
  CS_UNUSED(rho_f);
  CS_UNUSED(nu_f);
  CS_UNUSED(cp_f);
  CS_UNUSED(k_f);
  CS_UNUSED(tauc);
  CS_UNUSED(dt);
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

#pragma weak cs_user_lagr_sde
void
cs_user_lagr_sde(const cs_real_t  dt[],
                 cs_real_t        taup[],
                 cs_real_3_t      tlag[],
                 cs_real_t        tempct[])
{
  CS_UNUSED(dt);
  CS_UNUSED(taup);
  CS_UNUSED(tlag);
  CS_UNUSED(tempct);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
