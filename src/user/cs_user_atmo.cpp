/*============================================================================
 * User-defined functions specific to amospheric flow models.
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_atmo.cpp
 *
 * \brief User-defined functions specific to amospheric flow models.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill in vertical profiles of atmospheric properties prior to solving
 *        1D radiative transfers.
 *
 * \param[in, out] preray        pressure vertical profile
 * \param[in, out] temray        real temperature vertical profile
 * \param[in, out] romray        density vertical profile
 * \param[in, out] qvray         water vapor content vertical profile
 * \param[in, out] qlray         water liquid content vertical profile
 * \param[in, out] ncray         droplets density vertical profile
 * \param[in, out] aeroso        aerosol concentration vertical profile
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_atmo_1d_rad_prf
void
cs_user_atmo_1d_rad_prf(cs_real_t   preray[],
                        cs_real_t   temray[],
                        cs_real_t   romray[],
                        cs_real_t   qvray[],
                        cs_real_t   qlray[],
                        cs_real_t   ncray[],
                        cs_real_t   aeroso[])
{
  CS_UNUSED(qvray);
  CS_UNUSED(qlray);
  CS_UNUSED(ncray);
  CS_UNUSED(preray);
  CS_UNUSED(temray);
  CS_UNUSED(romray);
  CS_UNUSED(aeroso);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
