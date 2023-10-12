/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <string.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_user_radiative_transfer.c
 *
 * \brief User function for input of radiative transfer parameters:
 *        absorption coefficient and net radiation flux.
 *
 *  See \ref cs_user_radiative_transfer for examples.
 */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for input of radiative transfer module options.
 *
 * \deprecated Use cs_user_model instead.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_radiative_transfer_parameters
void
cs_user_radiative_transfer_parameters(void)
{

}

/*-----------------------------------------------------------------------------*/
/*!
 * \brief Absorption coefficient for radiative module
 *
 * It is necessary to define the value of the fluid's absorption coefficient Ck.
 *
 * This value is defined automatically for specific physical models, such
 * as gas and coal combustion, so this function is not used by these models.
 *
 * For a transparent medium, the coefficient should be set to 0.
 *
 * In the case of the P-1 model, we check that the optical length is at
 * least of the order of 1.
 *
 * \param[in]   bc_type  boundary face types
 * \param[out]  ck       medium's absorption coefficient (zero if transparent)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_rad_transfer_absorption
void
cs_user_rad_transfer_absorption(const int  bc_type[],
                                cs_real_t  ck[])
{
  CS_UNUSED(bc_type);
  CS_UNUSED(ck);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the net radiation flux.
 *
 * The density of net radiation flux must be calculated
 * consistently with the boundary conditions of the intensity.
 * The density of net flux is the balance between the radiative
 * emiting part of a boundary face (and not the reflecting one)
 * and the radiative absorbing part.
 *
 * \param[in]   bc_type   boundary face types
 * \param[in]   twall     inside current wall temperature (K)
 * \param[in]   qincid    radiative incident flux  (W/m2)
 * \param[in]   xlam      conductivity (W/m/K)
 * \param[in]   epa       thickness (m)
 * \param[in]   eps       emissivity (>0)
 * \param[in]   ck        absorption coefficient
 * \param[out]  net_flux  net flux (W/m2)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_rad_transfer_net_flux
void
cs_user_rad_transfer_net_flux(const int        bc_type[],
                              const cs_real_t  twall[],
                              const cs_real_t  qincid[],
                              const cs_real_t  xlam[],
                              const cs_real_t  epa[],
                              const cs_real_t  eps[],
                              const cs_real_t  ck[],
                              cs_real_t        net_flux[])
{
  CS_UNUSED(bc_type);
  CS_UNUSED(twall);
  CS_UNUSED(qincid);
  CS_UNUSED(xlam);
  CS_UNUSED(epa);
  CS_UNUSED(eps);
  CS_UNUSED(ck);
  CS_UNUSED(net_flux);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
