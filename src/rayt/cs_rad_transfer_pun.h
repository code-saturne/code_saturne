#ifndef __CS_RAD_TRANSFER_PUN_H__
#define __CS_RAD_TRANSFER_PUN_H__

/*============================================================================
 * Radiation solver boundary conditions treatment.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Radiative flux and source term computation
 *
 * \param[in]       iband     number of the i-th gray gas
 * \param[in]       bc_type   boundary face types
 * \param[in, out]  coefap    boundary condition work array for the luminance
 *                             (explicit part)
 * \param[in, out]  coefbp    boundary condition work array for the luminance
 *                             (implicit part)
 * \param[in, out]  cofafp    boundary condition work array for the diffusion
 *                             of the luminance (explicit part)
 * \param[in, out]  cofbfp    boundary condition work array for the diffusion
 *                             of the luminance (implicit part)
 * \param[in, out]  flurds    pseudo mass flux work array (interior faces)
 * \param[in, out]  flurdb    pseudo mass flux work array (boundary faces)
 * \param[in, out]  viscf     visc*surface/dist work array at interior faces
 * \param[in, out]  viscb     visc*surface/dist work array at boundary faces
 * \param[in, out]  smbrs     work array for RHS
 * \param[in, out]  rovsdt    work array for unsteady term
 * \param[in]       twall     wall temperature in Kelvin
 * \param[in, out]  ckmel     absorption coefficient for gas-particles mix
 * \param[out]      q         explicit flux density vector
 * \param[in]       abo       weights of the i-th gray gas at boundaries
 * \param[out]      int_rad_domega  integral of I dOmega
 * \param[out]      theta4    bulk absorption
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_pun(int              iband,
                    int              bc_type[],
                    cs_real_t        coefap[],
                    cs_real_t        coefbp[],
                    cs_real_t        cofafp[],
                    cs_real_t        cofbfp[],
                    cs_real_t        flurds[],
                    cs_real_t        flurdb[],
                    cs_real_t        viscf[],
                    cs_real_t        viscb[],
                    cs_real_t        smbrs[],
                    cs_real_t        rovsdt[],
                    cs_real_t        twall[],
                    cs_real_t        ckmel[],
                    cs_real_3_t      q[],
                    const cs_real_t  abo[],
                    cs_real_t        int_rad_domega[],
                    cs_real_t        theta4[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_PUN_H__ */
