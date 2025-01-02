#ifndef __CS_SEARCH_H__
#define __CS_SEARCH_H__

/*============================================================================
 * Management of the mass flux, the viscosity, the density, the specific
 * heat and the tsnsa array in case of a theta-scheme.
 *===========================================================================*/

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "base/cs_base.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Management of the mass flux, the viscosity, the density, the specific
 * heat  and the tsnsa array in case of a theta-scheme.
 *
 * Please refer to the
 * <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
 * of the theory guide for more informations.
 *
 * parameters:
 *   iappel  <--  call id (before or after cs_physical_properties_update)
 *
 * returns:
 *---------------------------------------------------------------------------*/

void
cs_theta_scheme_update_var(const cs_lnum_t  iappel);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SEARCH_H__ */
