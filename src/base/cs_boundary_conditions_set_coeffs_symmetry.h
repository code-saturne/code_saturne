#ifndef __CS_BOUNDARY_CONDITIONS_SET_COEFFS_SYMMETRY_H__
#define __CS_BOUNDARY_CONDITIONS_SET_COEFFS_SYMMETRY_H__

/*============================================================================
 * Boundary condition management.
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
/*! \file cs_boundary_conditions_set_coeffs_symmetry.c
 *
 * \brief Compute the symmetry boundary condition coefficients.
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Symmetry boundary conditions for vectors and tensors.
 *
 * Corresponds to the code icodcl(ivar) = 4.
 *
 * Please refer to the
 * <a href="../../theory.pdf#clsyvt"><b>clsyvt</b></a> section of the
 * theory guide for more informations.
 *
 * \param[in]  velipb  value of the velocity at \f$ \centip \f$
 *                     of boundary cells
 * \param[in]  rijipb  value of \f$ R_{ij} \f$ at \f$ \centip \f$
 *                     of boundary cells
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_coeffs_symmetry(cs_real_t  velipb[][3],
                                           cs_real_t  rijipb[][6]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_SET_COEFFS_SYMMETRY_H__ */
