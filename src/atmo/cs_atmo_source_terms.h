#ifndef __CS_ATMO_SOURCE_TERMS_H__
#define __CS_ATMO_SOURCE_TERMS_H__

/*============================================================================
 * Main for atmospheric source terms
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

#include "base/cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Phase change source terms - Exchange terms between the injected
 *        liquid and the water vapor phase in the bulk, humid air
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 * \param[in,out] imp_st        Implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_source_term(int              f_id,
                    cs_real_t        exp_st[],
                    cs_real_t        imp_st[]);

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Additional right-hand side source terms for scalar equations
 *        taking into account dry and humid atmospheric variables.
 *        If 1D atmospheric radiative module is used additional source terms for
 *        the thermal scalar equation to take into account the radiative forcing.
 *
 * \param[in]     f_id          field id
 * \param[in,out] exp_st        Explicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_scalar_source_term(int              f_id,
                           cs_real_t        exp_st[]);


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ATMO_SOURCE_TERMS_H__ */
