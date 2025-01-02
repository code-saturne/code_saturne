#ifndef __CS_1D_WALL_THERMAL_CHECK_H__
#define __CS_1D_WALL_THERMAL_CHECK_H__

/*============================================================================
 * Data checking for the 1D thermal wall module.
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Data checking for the 1D thermal wall module.
 *
 * \param[in]   iappel   Call number:
 *                       - 1: first call during the initialization (called once)
 *                       Checking the number of cells nfpt1d.
 *                       - 2: second call during the initialization (called once)
 *                       Checking ifpt1d, nppt1d, eppt1d and rgpt1d arrays.
 *                       - 3: called at each time step
 *                       Checking iclt1d, xlmbt1, rcpt1d and dtpt1d arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_1d_wall_thermal_check(int  iappel);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_1D_WALL_THERMAL_CHECK_H__ */
