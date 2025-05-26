#ifndef __CS_ATMO_SOLIVA_H__
#define __CS_ATMO_SOLIVA_H__

/*============================================================================
 * Atmospheric soil module - soil variables initialization
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

/*----------------------------------------------------------------------------
 * Standard headers
 *----------------------------------------------------------------------------*/

#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"
#include "base/cs_physical_constants.h"
#include "atmo/cs_atmo.h"
#include "base/cs_field.h"
#include "mesh/cs_mesh.h"
#include "atmo/cs_air_props.h"
#include "base/cs_math.h"
#include "alge/cs_divergence.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function declarations
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize soil model variables
 */
/*----------------------------------------------------------------------------*/

void cs_soliva(void);

/*----------------------------------------------------------------------------*/

#endif /* __CS_ATMO_SOLIVA_H__ */
