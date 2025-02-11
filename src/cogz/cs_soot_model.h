#ifndef __CS_SOOT_MODEL_H__
#define __CS_SOOT_MODEL_H__

/*============================================================================
 * Soot production models.
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
#include "base/cs_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Specific physic subroutine : soot production model.
 *
 * Define the source terms for the soot mass fraction
 * and the precursor number for soot model.
 * The equations read: \f$ rovsdt \delta a = smbrs \f$
 *
 * \f$ rovsdt \f$ et \f$ smbrs \f$ may already contain source term
 * so must not be overwritten, but incremented.
 *
 * For stability, only positive terms should be add in \f$ rovsdt \f$.
 * There is no constraint for \f$ smbrs \f$.
 * For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
 *           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
 *           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
 *
 * Here we set \f$ rovsdt \f$ and \f$ smbrs \f$ containing \f$ \rho \Omega \f$
 *   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
 *     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.C.s^{-1} \f$,
 *     for enthalpy: \f$ J.s^{-1} \f$)
 *   - \f$ rovsdt \f$ in \f$ kg.s^{-1} \f$
 *
 * \param[in]     f_id           scalar field id
 * \param[out]    smbrs          explicit right hand side
 * \param[out]    rovsdt         implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_soot_production(int        f_id,
                   cs_real_t  smbrs[],
                   cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOOT_MODEL__ */
