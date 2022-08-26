#ifndef __CS_TURBULENCE_RIJ_H__
#define __CS_TURBULENCE_RIJ_H__

/*============================================================================
 * Rij-epsilon turbulence model.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve of epsilon for \f$ R_{ij} - \varepsilon \f$ RANS
 *        turbulence model.
 *
 * \param[in]     ncesmp      number of cells with mass source term
 * \param[in]     icetsm      index of cells with mass source term
 * \param[in]     itypsm      type of mass source term for each variable
 * \param[in]     gradv       work array for the term grad
 *                            of velocity only for iturb=31
 * \param[in]     produc      work array for production (without
 *                            rho volume) only for iturb=30
 * \param[in]     gradro      work array for \f$ \grad{rom} \f$
 * \param[in]     smacel      value associated to each variable in the mass
 *                            source terms or mass rate
 * \param[in]     viscf       visc*surface/dist at internal faces
 * \param[in]     viscb       visc*surface/dist at edge faces
 * \param[in]     rhs         working array
 * \param[in]     rovsdt      working array
 !*/
/*-------------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_eps(cs_lnum_t            ncesmp,
                            cs_lnum_t            icetsm[],
                            int                  itypsm[],
                            const cs_real_33_t   gradv[],
                            const cs_real_6_t    produc[],
                            const cs_real_3_t    gradro[],
                            cs_real_t            smacel[],
                            cs_real_t            viscf[],
                            cs_real_t            viscb[],
                            cs_real_t            rhs[],
                            cs_real_t            rovsdt[]);

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_alpha(int        f_id,
                              cs_real_t  c_durbin_l);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gravity terms for terms
 *        For \f$R_{ij}\f$
 *
 * \param[in]   gradro    work array for \f$ \grad{\rho} \f$
 * \param[out]  buoyancy  buoyancy term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_grav_st(const cs_real_t  gradro[][3],
                          cs_real_t        buoyancy[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Terms of wall echo for \f$ R_{ij} \f$
 *        \f$var  = R_{11} \: R_{22} \: R_{33} \: R_{12} \: R_{13} \: R_{23}\f$
 *        \f$comp =  1 \:  2 \:  3 \:  4 \:  5 \:  6\f$
 *
 * \param[in]     produc  production
 * \param[in,out] rhs     work array for right-hand-side
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_echo(const cs_real_t  produc[][6],
                       cs_real_t        rhs[][6]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (coupled components version).
 *
 * \param[in]  n_cells  number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip(cs_lnum_t  n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (segregated version)
 *
 * \param[in]  n_cells  number of cells
 * \param[in]  iclip    if 0, viscl0 is used; otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip_sg(cs_lnum_t  n_cells,
                          int        iclip);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_RIJ_H__ */
