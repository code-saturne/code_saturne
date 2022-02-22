#ifndef __CS_ATMO_AEROSOL_H__
#define __CS_ATMO_AEROSOL_H__

/*============================================================================
 * Main for atmospheric aerosols related functions
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

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
 * \brief This function initializes the external aerosol code
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function finalizes the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with aerosol concentrations
 *        and numbers from the external aerosol code..
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_aero(cs_real_t*);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function fills the given array with gas concentrations from
 *        the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_get_gas(cs_real_t*);

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function computes a time step of gaseous chemistry and aerosols
 *        dynamic using the external aerosol code.
 */
/*----------------------------------------------------------------------------*/

void
cs_atmo_aerosol_time_advance(void);

END_C_DECLS

#endif /* __CS_ATMO_AEROSOL_H__ */
