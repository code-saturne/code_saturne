#ifndef __CS_LAGR_PROTOTYPES_H__
#define __CS_LAGR_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "mesh/cs_mesh_bad_cells.h"

#include "cdo/cs_domain.h"

#include "lagr/cs_lagr.h"
#include "lagr/cs_lagr_tracking.h"
#include "lagr/cs_lagr_stat.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Lagrangian User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_lagr_boundary_conditions.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle boundary conditions.
 *
 * This is used for the definition of inlet and other boundaries,
 * based on predefined boundary zones (\ref cs_zone_t).
 *
 * \param[in] bc_type    type of the boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_boundary_conditions(const int  itypfb[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of a particle interaction with a boundary of type
 *        \ref CS_LAGR_BC_USER.
 *
 * \param[in, out]  particles       pointer to particle set
 * \param[in]       p_id            particle id
 * \param[in]       face_id         boundary face id
 * \param[in]       face_norm       unit face (or face subdivision) normal
 * \param[in]       c_intersect     coordinates of intersection with the face
 * \param[in]       t_intersect     relative distance (in [0, 1]) of the
 *                                  intersection point with the face relative
 *                                  to the initial trajectory segment
 * \param[in]       b_zone_id       boundary zone id of the matching face
 * \param[in, out]  event_flag      event flag in case events are available
 * \param[in, out]  tracking_state  particle tracking state
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_user_boundary_interaction(cs_lagr_particle_set_t    *particles,
                                  cs_lnum_t                  p_id,
                                  cs_lnum_t                  face_id,
                                  const cs_real_t            face_norm[3],
                                  const cs_real_t            c_intersect[3],
                                  cs_real_t                  t_intersect,
                                  int                        b_zone_id,
                                  int                       *event_flag,
                                  cs_lagr_tracking_state_t  *tracking_state);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_lagr_model.cpp
 */
/*----------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*
 * \brief User function of the Lagrangian particle-tracking module
 *
 *  User input of physical, numerical and post-processing options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_lagr_particle.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of an external force field acting on the particles.
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
 * \param[in,out]  rho_p    particle density
 * \param[out]     fextla   user external force field (m/s^2)$
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_ef(cs_real_t            dt_p,
                const cs_lnum_t      p_id,
                const cs_real_t     *taup,
                const cs_real_3_t   *tlag,
                const cs_real_3_t   *piil,
                const cs_real_33_t  *bx,
                const cs_real_t      tsfext,
                const cs_real_3_t   *vagaus,
                const cs_real_3_t    gradpr,
                const cs_real_33_t   gradvf,
                cs_real_t            rho_p,
                cs_real_3_t          fextla);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function (non-mandatory intervention)
 *
 * User-defined modifications on the variables at the end of the
 * Lagrangian time step and calculation of user-defined
 * additional statistics on the particles.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_extra_operations(const cs_real_t  dt[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Impose the motion of a particle flagged CS_LAGR_PART_IMPOSED_MOTION.
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

void
cs_user_lagr_imposed_motion(const cs_lagr_particle_set_t *particles,
                            cs_lnum_t                     p_id,
                            const cs_real_t               coords[3],
                            const cs_real_t               dt,
                            cs_real_t                     disp[3]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particle location and cell_id.
 *
 * This function is called before the initialization of the new particles in
 * order to initialize them according to new particle positions, to this aim it
 * is necessary to modify the cell_id associated to the particle accordingly.
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modified by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               injection data for this set
 * \param[in]      particle_range    start and past-the-end ids of new particles
 *                                   for this zone and class
 * \param[in]      particle_face_id  face ids of new particles if zone is
 *                                   a boundary,  NULL otherwise
 * \param[in]      visc_length       viscous layer thickness
 *                                   (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_in_force_coords(cs_lagr_particle_set_t         *particles,
                             const cs_lagr_injection_set_t  *zis,
                             const cs_lnum_t                 particle_range[2],
                             const cs_lnum_t                 particle_face_id[],
                             const cs_real_t                 visc_length[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of newly injected particles.
 *
 * This function is called after the initialization of the new particles in
 * order to modify them according to new particle profiles (injection
 * profiles, statistical weights, correction of the diameter if the
 * standard-deviation option is activated); the modification of particles
 * position should preferentially be made in cs_user_lagr_in_force_coords to
 * get an initialization of particle properties coherent with the local fields.
 *
 * This function is called for each injection zone and class. Particles
 * with ids between \c pset->n_particles and \c n_elts are initialized
 * but may be modified by this function.
 *
 * \param[in,out]  particles         particle set
 * \param[in]      zis               injection data for this set
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
                const cs_real_t                 visc_length[]);

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
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_rt(int              phase_id,
                cs_lnum_t        id_p,
                cs_real_t        re_p,
                cs_real_t        uvwr,
                cs_real_t        rho_f,
                cs_real_t        rho_p,
                cs_real_t        nu_f,
                cs_real_t       *taup,
                const cs_real_t  dt);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the calculation of the thermal relaxation time of the
 *   particles with respect to the chosen formulation of the Nusselt number.
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
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
                  cs_real_2_t      tempct,
                  const cs_real_t  dt);

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

void
cs_user_lagr_sde(const cs_real_t         dt,
                 const cs_lnum_t         p_id,
                 const cs_real_t        *taup,
                 const cs_real_3_t      *tlag,
                 const cs_real_2_t       tempct,
                 const int               nor);

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_lagr_volume_conditions.cpp
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle volume conditions.
 *
 * This is used definition of for injection
 * based on predefined volume zones (\ref cs_zone_t).
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_volume_conditions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of a particle interaction with a interior face of type
 *        \ref CS_LAGR_BC_USER.
 *
 * \param[in, out]  particles       pointer to particle set
 * \param[in]       p_id            particle id
 * \param[in]       face_id         interior face id
 * \param[in]       face_u_norm     unit face (or face subdivision) normal
 * \param[in]       c_intersect     coordinates of intersection with the face
 * \param[in]       t_intersect     relative distance (in [0, 1]) of the
 *                                  intersection point with the face relative
 *                                  to the initial trajectory segment
 * \param[in, out]  tracking_state  particle tracking state
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_user_internal_interaction(cs_lagr_particle_set_t    *particles,
                                  cs_lnum_t                  p_id,
                                  cs_lnum_t                  face_id,
                                  const cs_nreal_t           face_u_norm[3],
                                  const cs_real_t            c_intersect[3],
                                  cs_real_t                  t_intersect,
                                  cs_lagr_tracking_state_t  *tracking_state);

/*----------------------------------------------------------------------------*/

#endif /* __CS_LAGR_PROTOTYPES_H__ */

END_C_DECLS
