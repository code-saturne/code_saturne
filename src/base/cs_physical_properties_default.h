#ifndef __CS_PHYSICAL_PROPERTIES_DEFAULT_H__
#define __CS_PHYSICAL_PROPERTIES_DEFAULT_H__

/*============================================================================
 *  physical properties variable in time.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fills physical properties which are variable in time
 *
 *\param[in]     iterns      Navier-Stokes sub-iterations indicator:
 *                           - if strictly negative, indicate that this
 *                                function is called outside Navier-Stokes loop
 *                           - if positive, Navier-Stokes iteration number.
 */
/*----------------------------------------------------------------------------*/

void
cs_physical_properties_update(int   iterns);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PHYSICAL_PROPERTIES_DEFAULT__H__ */
