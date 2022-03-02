#ifndef __CS_TURBULENCE_VIS_LES_H__
#define __CS_TURBULENCE_VIS_LES_H__

/*============================================================================
 * Turbulent viscosity for LES models.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *        a dynamic Smagorinsky LES model
 *
 * \f[ smago = \dfrac{L_{ij}M_{ij}}{M_{ij}M_{ij}} \f]
 *
 * \f[ \mu_T = \rho smago L^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 *
 * \param[out]  gradv the computed velocity gradients
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_dyn(cs_real_33_t *gradv);

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *        a Smagorinsky LES model
 *
 * \f[ \mu_T = \rho (C_{S} l)^2  \sqrt{2 S_{ij}S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 *
 * \param[out]  gradv the computed velocity gradients
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_smago_const(cs_real_33_t *gradv);

/*----------------------------------------------------------------------------*/
/*! \brief Compute the turbulent viscosity for the WALE LES model
 *
 * The turbulent viscosity is:
 * \f$ \mu_T = \rho (C_{wale} L)^2 * \dfrac{(\tens{S}:\tens{Sd})^{3/2}}
 *                                         {(\tens{S} :\tens{S})^(5/2)
 *                                         +(\tens{Sd}:\tens{Sd})^(5/4)} \f$
 * with \f$ \tens{S}  = \frac{1}{2}(\gradt \vect{u} + \transpose{\gradt \vect{u}})\f$
 * and  \f$ \tens{Sd} = \deviator{(\symmetric{(\tens{S}^2)})}\f$
 *
 * \param[out]  gradv the computed velocity gradients
 !*/
/*----------------------------------------------------------------------------*/

void
cs_les_mu_t_wale(cs_real_33_t *restrict gradv);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LES_MU_T_H__ */
