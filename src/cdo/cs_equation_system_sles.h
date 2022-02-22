#ifndef __CS_EQUATION_SYSTEM_SLES_H__
#define __CS_EQUATION_SYSTEM_SLES_H__

/*============================================================================
 * Functions to handle the solving step of systems of equations hinging on the
 * cs_equation_system_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_equation.h"
#include "cs_equation_system_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the SLES associated to a system of equation
 *
 * \param[in]  n_eqs     number of equations in the system to solve
 * \param[in]  sysp      set of paremeters for the system of equations
 * \param[in]  blocks    array of the core members for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_sles_init(int                            n_eqs,
                             cs_equation_system_param_t    *sysp,
                             cs_equation_core_t           **blocks);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_SYSTEM_SLES_H__ */
