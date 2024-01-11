#ifndef __CS_SOLVE_EQUATION_H__
#define __CS_SOLVE_EQUATION_H__

/*============================================================================
 * Solve the convection diffusion equation with additional terms.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the convection/diffusion equation (with optional source
 *        terms and/or drift) for a scalar quantity over a time step.
 *
 * \param[in]     f          pointer to field structure
 * \param[in]     ncesmp     number of cells with mass source term
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param[in]     itspdv     indicator to compute production/dissipation
 *                           terms for a variance:
 *                           - 0: no
 *                           - 1: yes
 * \param[in]     icetsm     index of cells with mass source term
 * \param[in]     itypsm     type of mass source term for the variables
 * \param[in]     smacel     variable value associated to the mass source
 *                           term (for ivar=ipr, smacel is the mass flux
 *                           \f$ \Gamma^n \f$)
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_scalar(cs_field_t        *f,
                         cs_lnum_t          ncesmp,
                         int                iterns,
                         int                itspdv,
                         const cs_lnum_t    icetsm[],
                         int                itypsm[],
                         cs_real_t          smacel[],
                         cs_real_t          viscf[],
                         cs_real_t          viscb[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the convection/diffusion equation (with optional source
 *        terms and/or drift) for a vectorial quantity over a time step..
 *
 * \param[in]     f          pointer to field structure
 * \param[in]     ncesmp     number of cells with mass source term
 * \param[in]     iterns     Navier-Stokes iteration number
 * \param[in]     icetsm     index of cells with mass source term
 * \param[in]     itypsm     type of mass source term for the variables
 * \param[in]     smacel     variable value associated to the mass source
 *                           term (for ivar=ipr, smacel is the mass flux
 *                           \f$ \Gamma^n \f$)
 * \param         viscf      visc*surface/dist at internal faces (work array)
 * \param         viscb      visc*surface/dist at boundary faces (work array)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_equation_vector(cs_field_t       *f,
                         const cs_lnum_t   ncesmp,
                         int               iterns,
                         const cs_lnum_t   icetsm[],
                         int               itypsm[],
                         cs_real_t         smacel[],
                         cs_real_t         viscf[],
                         cs_real_t         viscb[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLVE_EQUATION_H__ */
