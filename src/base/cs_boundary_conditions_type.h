#ifndef __CS_BOUNDARY_CONDITIONS_TYPE_H__
#define __CS_BOUNDARY_CONDITIONS_TYPE_H__

/*============================================================================
 * Handle boundary condition type codes.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \file cs_boundary_conditions_type.c
 *
 * \brief Handle boundary condition type codes (\ref itypfb).
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*
 * \brief  Handle boundary condition type code (\ref bc_type).
 *
 * \param[in]      init      partial treatment (before time loop) if true
 * \param[in,out]  bc_type   type per boundary face
 * \param[out]     itrifb    indirection for faces ordering
 * \param[out]     isostd    standard output indicator + reference face number
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_type(bool  init,
                            int   bc_type[],
                            int   itrifb[],
                            int   isostd[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_BOUNDARY_CONDITIONS_TYPE_H__ */
