#ifndef CS_COAL_SOURCE_TERMS_H
#define CS_COAL_SOURCE_TERMS_H

/*============================================================================
 * Coal combustion model: source term computation
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute scalar source terms for pulverized coal flame.
 *
 * \warning  the treatement is different from that of cs_user_source_terms.
 *
 * We solve: \f[ rovsdt D(var) = smbrs \f]
 *
 * rovsdt and smbrs already contain eventual user source terms.
 * So they have to be incremented and not erased.
 *
 * For stability reasons, only positive terms can be added in rovsdt.
 * There is no contraint for smbrs.
 *
 * In the case of a source term in \f$ cexp + cimp var \f$, it has to be written:
 *        - \f$ smbrs  = smbrs  + cexp + cimp var \f$
 *        - \f$ rovsdt = rovsdt + \max(-cimp,0) \f$
 *
 * Here are \f$ rovsdt \f$ and \f$ smbrs \f$ (they contain \f$ \rho volume\f$)
 *    smbrs in kg variable/s:
 *     \c i.e.: - for velocity            \f$ kg . m . s^{-2} \f$
 *              - for temperature         \f$ kg . [degres] . s^{-1} \f$
 *              - for enthalpy            \f$ J . s^{-1} \f$
 *              - rovsdt                  \f$ kg . s^{-1} \f$
 *
 * \param[in]      fld_id  scalar field id
 * \param[in,out]  smbrs   explicit second member
 * \param[in,out]  rovsdt  implicit diagonal part
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_source_terms_scalar(int        fld_id,
                            cs_real_t  smbrs[],
                            cs_real_t  rovsdt[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Pulverized coal flame: production and dissipation source terms
 *        for the variance.
 *
 * \param[in]      iscal   scalar id
 * \param[in,out]  smbrs   explicit second member
 * \param[in,out]  rovsdt  implicit diagonal part
 */
/*----------------------------------------------------------------------------*/

void
cs_coal_fp2st(int         iscal,
              cs_real_t  *smbrs,
              cs_real_t  *rovsdt);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* CS_COAL_SOURCE_TERMS_H */
