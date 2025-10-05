/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* VERS */

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * User definition of radiative transfer boundary conditions.
 *
 * The temperature unit here is the Kelvin
 *
 * The following face characteristics must be set:
 *  - isothp(face_id) boundary face type
 *
 * Depending on the value of isothp, other values may also need to be set:
 *  - rcodcl = conduction flux
 *  - epsp   = emissivity
 *  - xlamp  = conductivity (W/m/K)
 *  - epap   = thickness (m)
 *  - textp  = outside temperature (K)
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       bc_type       boundary face types
 * \param[in]       isothp        boundary face type for radiative transfer
 * \param[out]      tmin          min allowed value of the wall temperature
 * \param[out]      tmax          max allowed value of the wall temperature
 * \param[in]       tx            relaxation coefficient (0 < tx < 1)
 * \param[in]       dt            time step (per cell)
 * \param[in]       thwall        inside current wall temperature (K)
 * \param[in]       qincid        radiative incident flux  (W/m2)
 * \param[in]       hfcnvp        convective exchange coefficient (W/m2/K)
 * \param[in]       flcnvp        convective flux (W/m2)
 * \param[out]      xlamp         conductivity (W/m/K)
 * \param[out]      epap          thickness (m)
 * \param[out]      epsp          emissivity (>0)
 * \param[out]      textp         outside temperature (K)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_radiative_transfer_bcs
void
cs_user_radiative_transfer_bcs
(
  [[maybe_unused]] cs_domain_t      *domain,
  [[maybe_unused]] const int         bc_type[],
  [[maybe_unused]] int               isothp[],
  [[maybe_unused]] cs_real_t        *tmin,
  [[maybe_unused]] cs_real_t        *tmax,
  [[maybe_unused]] cs_real_t        *tx,
  [[maybe_unused]] const cs_real_t   dt[],
  [[maybe_unused]] const cs_real_t   thwall[],
  [[maybe_unused]] const cs_real_t   qincid[],
  [[maybe_unused]] cs_real_t         hfcnvp[],
  [[maybe_unused]] cs_real_t         flcnvp[],
  [[maybe_unused]] cs_real_t         xlamp[],
  [[maybe_unused]] cs_real_t         epap[],
  [[maybe_unused]] cs_real_t         epsp[],
  [[maybe_unused]] cs_real_t         textp[])
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
