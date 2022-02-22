#ifndef __CS_ATPRKE_H__
#define __CS_ATPRKE_H__

/*============================================================================
 * Modify the k-epsilon turbulence model for the atmospheric module.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify the \f$k-\varepsilon\f$ turbulence model
 *        formulation for the atmospheric module
 *
 * Adjunction of a production term for buyancy in the \f$k-\varepsilon\f$
 * model in the context of the atmospheric module
 * g = g*grad(theta)/prdtur/theta
 *
 * \param[in, out]  tinstk    Implicit part of the buoyancy term (for k)
 * \param[in, out]  smbrk     Explicit part of the buoyancy term (for k)
 * \param[in, out]  smbre     Explicit part of the buoyancy term (for eps)
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_buoyancy_ke_prod(cs_real_t  *tinstk,
                         cs_real_t  *smbrk,
                         cs_real_t  *smbre);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATPRKE_H__ */
