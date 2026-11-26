#ifndef CS_TURBULENCE_RIJ_PROJECTION_H
#define CS_TURBULENCE_RIJ_PROJECTION_H

/*============================================================================
 * Godunov-type projection of the convective subsystem of second-moment
 * turbulence closures.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Project the mean velocity, the Reynolds stresses and the
 *        dissipation rate over one time step of the convective subsystem.
 *
 * The projected states are stored in the cell fields
 * "algo:proj_velocity", "algo:proj_rij" and "algo:proj_epsilon"
 * (created in cs_setup.cpp when rij_discretization_scheme ==
 * CS_RIJ_SCHEME_GODUNOV_PROJECTION). They are computed once per time step,
 * from the fields at the beginning of the time step; further calls in the
 * same time step return immediately.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_projection(void);

/*----------------------------------------------------------------------------*/

#endif /* CS_TURBULENCE_RIJ_PROJECTION_H */
