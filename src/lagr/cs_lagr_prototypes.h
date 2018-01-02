#ifndef __CS_LAGR_PROTOTYPES_H__
#define __CS_LAGR_PROTOTYPES_H__

/*============================================================================
 * Prototypes for Fortran functions and subroutines callable from C
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"

#include "cs_domain.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_stat.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Lagrangian User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of an external force field acting on the particles.
 *
 * It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
 * By default gravity and drag force are the only forces acting on the particles
 * (the gravity components gx gy gz are assigned in the GUI or in usipsu)
 *
 * \param[in]      dt_p     time step (for the cell)
 * \param[in]      taup     particle relaxation time
 * \param[in]      tlag     relaxation time for the flow
 * \param[in]      piil     term in the integration of the sde
 * \param[in]      bx       characteristics of the turbulence
 * \param[in]      tsfext   infos for the return coupling
 * \param[in]      vagaus   Gaussian random variables
 * \param[in]      gradpr   pressure gradient
 * \param[in]      gradvf   gradient of the flow velocity
 * \param[in,out]  rho_p     particle density
 * \param[out]     fextla   user external force field (m/s^2)$
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
                cs_real_3_t          fextla[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for the boundary conditions for the particles (inlet
 *   and treatment for the other boundaries)
 *  This routine is called after the initialization of the new particles in order
 *  to modify them according to new particle profiles.
 *
 * \param[in] time_id         time step indicator for fields
 *                            0: use fields at current time step
 *                            1: use fields at previous time step
 * \param[in] injfac          array of injection face id for every particles
 * \param[in] local_userdata  local_userdata pointer to zone/cluster specific
 *                            boundary conditions (number of injected
 *                            particles, velocity profile...)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_in(int                         time_id,
                int                        *injfac,
                cs_lagr_zone_class_data_t  *local_userdata,
                cs_real_t                   vislen[]  );

/*---------------------------------------------------------------------------------*/
/* \brief User subroutine of the Lagrangian particle-tracking module
 *
 *  User subroutine for input of calculation parameters.
 *  This parameters concerns physical, numerical and post-processing options.
 */
/*---------------------------------------------------------------------------------*/

void
cs_user_lagr_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prescribe some attributes for newly injected particles.
 *
 * This function is called at different points, at which different attributes
 * may be modified.
 *
 * \param[in,out] particle  particle structure
 * \param[in]     p_am      particle attributes map
 * \param[in]     face_id   id of particle injection face
 * \param[in]     attr_id   id of variable modifiable by this call. called for
                            CS_LAGR_PRED_VELOCITY, CS_LAGR_DIAMETER,
                            CS_LAGR_TEMPERATURE, CS_LAGR_STAT_WEIGHT
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_new_p_attr(unsigned char                  *particle,
                        const cs_lagr_attribute_map_t  *p_am,
                        cs_lnum_t                       face_id,
                        cs_lagr_attribute_t             attr_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modification of the calculation of the particle relaxation time
 *  with respect to the chosen formulation for the drag coefficient
 *
 * This function is called in a loop on the particles, so be careful
 * to avoid too costly operations.
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
                const cs_real_t  dt[]);

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
                  cs_real_t        tauc[],
                  const cs_real_t  dt[]);

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
                            const cs_real_t  dt,
                            cs_real_t        disp[3]);

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
 * \brief User function (non-mandatory intervention)
 *     Integration of the sde for the user-defined variables.
 *     The variables are constant by default.
 *     The sde must be of the form:
 * \f[
 *    \frac{dT}{dt}=\frac{T - PIP}{Tca}
 * \f]
 *     T : IIIIeme user-defined variable, given for the ip particle by
 *            T = EPTP(JVLS(IIII),IP)
 *            T = EPTPA(JVLS(IIII),IP)
 *     Tca : Characteristic time for the sde
 *           to be prescribed in the array auxl1
 *     PIP : Coefficient of the sde (pseudo right member)
 *           to be prescribed in the array auxl2
 *           If the chosen scheme is first order (nordre=1)
 *            then, at the first and only passage pip is expressed
 *            as a function of the quantities of the previous time step
 *            contained in eptpa
 *           If the chosen scheme is second order (nordre=2)
 *            then, at the first passage (nor=1) pip is expressed as
 *            a function of the quantities of the previous time step contained
 *            in eptpa, and at the second passage (nor=2) pip is expressed as
 *            a function of the quantities of the current time step
 *
 * \param[in] dt        time step (per cell)
 * \param[in] taup      particle relaxation time
 * \param[in] tlag      relaxation time for the flow
 * \param[in] tempct    characteristic thermal time and implicit source
 *                      term of return coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_sde(const cs_real_t  dt[],
                 cs_real_t        taup[],
                 cs_real_3_t      tlag[],
                 cs_real_t        tempct[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function of the Lagrangian particle-tracking module:
 *        User function for input of calculation parameters.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_model(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define particle boundary conditions.
 *
 *   This is used definition of for inlet and of the other boundaries
 *
 *   Boundary faces may be selected using the
 *    \ref cs_selector_get_b_face_num_list function.
 *
 * parameters:
 *
 * \param[in] itypfb    type of the boundary faces
 */
/* --------------------------------------------------------------------------- */

void
cs_user_lagr_boundary_conditions(const int  itypfb[]);

/*----------------------------------------------------------------------------*/

#endif /* __CS_LAGR_PROTOTYPES_H__ */

END_C_DECLS
