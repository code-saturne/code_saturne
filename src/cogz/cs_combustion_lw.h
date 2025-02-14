#ifndef CS_COMBUSTION_LW_H
#define CS_COMBUSTION_LW_H

/*============================================================================
 * Libby-Williams gas combustion model.
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
 * Standard library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_field.h"
#include "pprt/cs_physical_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Initialize specific fields for Libby-Williams gas combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_fields_init(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties for Libby-Williams combustion model.
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_physical_prop(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute physical properties for premixed flame Libby-Williams
 *        combustion model.
 *
 * Define the source terms for a given scalar over one time step.
 *
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
 * \param[in]      f_sc          pointer to scalar field
 * \param[in,out]  smbrs         explicit right hand side
 * \param[in,out]  rovsdt        implicit terms
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_lw_source_terms(cs_field_t  *f_sc,
                              cs_real_t    smbrs[],
                              cs_real_t    rovsdt[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COMBUSTION_LW_H */
