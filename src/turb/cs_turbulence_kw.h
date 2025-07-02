#ifndef __CS_TURBULENCE_KW_H__
#define __CS_TURBULENCE_KW_H__

/*============================================================================
 * k-w turbulence model.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip k and omega.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 * \param[in]     n_cells       number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_kw_clip(int       phase_id,
                      cs_lnum_t n_cells);


/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the k-omega equations.
 *
 * Solve the \f$ k - \omega \f$ SST for incompressible flows
 * or slightly compressible flows for one time step.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_kw(int phase_id);

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *         the \f$ k - \omega \f$ SST model.
 *
 * \f[ \mu_T = \rho A1 \dfrac{k}{\max(A1 \omega; \; S f_2)} \f]
 * with
 * \f[ S = \sqrt{  2 S_{ij} S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * and \f$ f_2 = \tanh(arg2^2) \f$
 * \f[ arg2^2 = \max(2 \dfrac{\sqrt{k}}{C_\mu \omega y}; \;
 *                   500 \dfrac{\nu}{\omega y^2}) \f]
 * where \f$ y \f$ is the distance to the wall.
 *
 * \f$ \divs{\vect{u}} \f$ is calculated at the same time than \f$ S \f$
 * for use in cs_turbulence_kw.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_kw_mu_t(int phase_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_KW_H__ */
