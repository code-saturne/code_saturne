/*============================================================================
! User synthetic turbulence inlet definition.
!
! 1) Global characteristics of synthetic turbulence inlets
! 2) Caracteristics of one specific inlet
! 3) Accurate specification of target statistics at inlet
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

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
 * \file cs_user_les_inflow.c
 *
 * \brief Generation of synthetic turbulence at LES inlets initialization.
 *
 * See \ref les_inflow for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define parameters of synthetic turbulence at LES inflow.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_les_inflow_define
void
cs_user_les_inflow_define(void)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update of the characteristics of a given synthetic turbulence inlet.
 *
 * \param[in]   zone       pointer to associated boundary zone
 * \param[out]  vel_r      reference mean velocity
 * \param[out]  k_r        reference turbulent kinetic energy
 * \param[out]  eps_r      reference turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_les_inflow_update
void
cs_user_les_inflow_update(const cs_zone_t  *zone,
                          cs_real_t         vel_r[3],
                          cs_real_t        *k_r,
                          cs_real_t        *eps_r)
{

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Definition of mean velocity, Reynolds stresses and dissipation rate
 *        for each boundary face of the given synthetic turbulence inlet.
 *
 * Accurate definition of mean velocity, Reynolds stresses and dissipation
 * rate for each boundary face of the given synthetic turbulence inlet
 *
 * Rij components are ordered as usual: XX, YY, ZZ, XY, YZ, XZ
 *
 * Arrays are pre-initialized before this function is called
 * (see \ref cs_user_les_inflow_define).
 *
 * vel_l[face_id][coo_id] = vel_r[coo_id]
 *
 * rij_l[face_id][0] = 2./3. * k_l
 * rij_l[face_id][1] = 2./3. * k_l
 * rij_l[face_id][2] = 2./3. * k_l
 * rij_l[face_id][3] = 0
 * rij_l[face_id][4] = 0
 * rij_l[face_id][5] = 0
 *
 * eps_l[face_id] = eps_r
 *
 * \param[in]       zone    pointer to associated boundary zone
 * \param[in, out]  vel_l   velocity a zone faces
 * \param[in, out]  rij_l   reynods stresses at zone faces
 * \param[in, out]  eps_l   reference turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_les_inflow_advanced
void
cs_user_les_inflow_advanced(const cs_zone_t  *zone,
                            cs_real_3_t       vel_l[],
                            cs_real_6_t       rij_l[],
                            cs_real_t         eps_l[])
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
