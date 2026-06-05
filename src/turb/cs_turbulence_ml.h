#ifndef CS_TURBULENCE_ML_H
#define CS_TURBULENCE_ML_H

/*============================================================================
 * Mixture-length turbulence model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Calculation of turbulent viscosity for
 *         a model of length of simple mixture
 *
 * \f[ \mu_T = \rho (\kappa L)^2 \cdot \sqrt{2 S_{ij} S_{ij}} \f]
 * \f[ S_{ij} = \dfrac{\der{u_i}{x_j} + \der{u_j}{x_i}}{2}\f]
 !*/
/*----------------------------------------------------------------------------*/

extern "C" void
cs_turbulence_ml_mu_t(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_ML_H */