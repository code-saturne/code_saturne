#ifndef __CS_PRESSURE_CORRECTION_H__
#define __CS_PRESSURE_CORRECTION_H__

/*============================================================================
 * Pressure correction.
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  /*!< Pressure correction step related to the mass
   * balance equation (scalar-valued) */

  cs_equation_t  *pressure_incr;

  /*! \var pressure_incr_gradient
   * Gradient of pressure increment. Used to store the gradient
   * pressure increment. */

  cs_field_t    *pressure_incr_gradient;

  /*! \var pressure_gradient
   * Gradient of pressure. Used to store the gradient of pressure */

  cs_field_t    *pressure_gradient;

  /*! \var div_st
   * Source term on the correction step stemming from the divergence of the
   * predicted velocity */

  cs_real_t      *div_st;

  /*! \var inner_potential_flux
   * Potential flux at interior faces. Used for Rhie & Chow */

  cs_real_t      *inner_potential_flux;

  /*! \var bdy_potential_flux
   * Potential flux at boundary faces. Used for Rhie & Chow */

  cs_real_t      *bdy_potential_flux;

  /*! \var bdy_pressure_incr
   * Pressure increment at the boundary. Used as an array to set the boundary
   * condition arising from a Dirichlet on the pressure. */

  cs_real_t      *bdy_pressure_incr;

  cs_flag_t      post_flag;

} cs_pressure_correction_cdo_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Activate the pressure increment solving with Legacy FV
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_fv_activate(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the pressure correction
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_destroy_all(void);


/*----------------------------------------------------------------------------*/
/*
 * \brief Activate the pressure increment solving with CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_activate(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Activate the pressure increment, either FV or CDO
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_model_activate(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Test if pressure solving with CDO is activated
 *
 * \return true if solving with CDO is requested, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_pressure_correction_cdo_is_activated(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Start setting-up the pressure increment equation
 *         At this stage, numerical settings should be completely determined
 *         but connectivity and geometrical information is not yet available.
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_init_setup(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Finalize setting-up the pressure increment equation
 *         At this stage, numerical settings should be completely determined
 *
 * \param[in] domain     pointer to a cs_domaint_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction_cdo_finalize_setup(const cs_domain_t   *domain);

/*----------------------------------------------------------------------------*/
/*
 * \brief Perform the pressure correction step of the Navier-Stokes equations
 *        for incompressible or slightly compressible flows.
 *
 * This function solves the following Poisson equation on the pressure:
 * \f[
 *     D \left( \Delta t, \delta p \right) =
 * \divs \left( \rho \vect{\widetilde{u}}\right)
 *     - \Gamma^n
 *     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
 * \f]
 *
 * Either the legacy FV method or a CDO face-based scheme is used.
 *
 * For the legacy case, the mass flux is  updated as follows:
 * \f[
 *  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
 *                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
 * \f]
 *
 * \Remark:
 * - an iterative process is used to solve the Poisson equation.
 * - if the arak coefficient is set to 1, the the Rhie & Chow filter is
 *   activated.
 *
 * Please refer to the
 * <a href="../../theory.pdf#resopv"><b>resopv</b></a>
 * section of the theory guide for more information.
 *
 * \param[in]       iterns        Navier-Stokes iteration number
 * \param[in]       nfbpcd        number of faces with condensation source term
 * \param[in]       ncmast        number of cells with condensation source terms
 * \param[in]       ifbpcd        index of faces with condensation source term
 * \param[in]       ltmast        list of cells with condensation source terms
 *                                (1 to n numbering)
 * \param[in]       isostd        indicator of standard outlet and index
 *                                of the reference outlet face
 * \param[in]       vel           velocity
 * \param[in, out]  da_uu         velocity matrix
 * \param[in]       bc_coeffs_v   boundary condition structure for the variable
 * \param[in]       bc_coeffs_dp  boundary conditions structure for the
 *                                pressure increment
 * \param[in]       spcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, spcond is the
 *                                flow rate
 *                                \f$ \Gamma_{s,cond}^n \f$)
 * \param[in]       svcond        variable value associated to the condensation
 *                                source term (for ivar=ipr, svcond is the flow rate
 *                                \f$ \Gamma_{v, cond}^n \f$)
 * \param[in]       frcxt         external forces making hydrostatic pressure
 * \param[in]       dfrcxt        variation of the external forces
 *                                composing the hydrostatic pressure
 * \param[in]       i_visc        visc*surface/dist aux faces internes
 * \param[in]       b_visc        visc*surface/dist aux faces de bord
 */
/*----------------------------------------------------------------------------*/

void
cs_pressure_correction(int                   iterns,
                       cs_lnum_t             nfbpcd,
                       cs_lnum_t             ncmast,
                       cs_lnum_t             ifbpcd[],
                       cs_lnum_t             ltmast[],
                       const int             isostd[],
                       cs_real_t             vel[][3],
                       cs_real_t             da_uu[][6],
                       cs_field_bc_coeffs_t *bc_coeffs_v,
                       cs_field_bc_coeffs_t *bc_coeffs_dp,
                       cs_real_t             spcond[],
                       cs_real_t             svcond[],
                       cs_real_t             frcxt[][3],
                       cs_real_t             dfrcxt[][3],
                       cs_real_t             i_visc[],
                       cs_real_t             b_visc[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PRESSURE_CORRECTION_H__ */
