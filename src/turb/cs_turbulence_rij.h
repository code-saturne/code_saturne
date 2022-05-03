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
 * \brief Solve the coupled Reynolds stress components in the
 *        \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
 *
 * \param[in]     field_id      index of current field
 * \param[in]     gradv         work array for the velocity grad term
 *                                 only for iturb=31
 * \param[in]     produc        work array for production
 * \param[in]     gradro        work array for grad rom
 *                              (without rho volume) only for iturb=30
 * \param[out]    viscf         visc*surface/dist at internal faces
 * \param[out]    viscb         visc*surface/dist at edge faces
 * \param[out]    viscce        Daly Harlow diffusion term
 * \param[out]    rhs           working array
 * \param[out]    rovsdt        working array
 * \param[out]    weighf        working array
 * \param[out]    weighb        working array
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_lrr(int              field_id,
                            const cs_real_t  gradv[][3][3],
                            const cs_real_t  produc[][6],
                            const cs_real_t  gradro[][3],
                            cs_real_t        viscf[],
                            cs_real_t        viscb[],
                            cs_real_t        viscce[][6],
                            cs_real_t        rhs[][6],
                            cs_real_t        rovsdt[][6][6],
                            cs_real_t        weighf[][2],
                            cs_real_t        weighb[]);

/*----------------------------------------------------------------------------*/
/*!/
 * \brief Solve the segregated Reynolds stress components in the
 *        \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
 *
 * \param[in]     field_id      index of current field
 *
 * \param[in]     produc        work array for production
 * \param[in]     gradro        work array for grad rom
 *                              (without rho volume) only for iturb=30
 * \param[out]    viscf         visc*surface/dist at internal faces
 * \param[out]    viscb         visc*surface/dist at edge faces
 * \param[out]    viscce        Daly Harlow diffusion term
 * \param[out]    rhs           working array
 * \param[out]    rovsdt        working array
 * \param[out]    weighf        working array
 * \param[out]    weighb        working array
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_lrr_sg(int              field_id,
                               const cs_real_t  produc[][6],
                               const cs_real_t  gradro[][3],
                               cs_real_t        viscf[],
                               cs_real_t        viscb[],
                               cs_real_t        viscce[][6],
                               cs_real_t        rhs[][6],
                               cs_real_t        rovsdt[][6][6],
                               cs_real_t        weighf[][2],
                               cs_real_t        weighb[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve the coupled Reynolds stress components in the
 *        \f$ R_{ij} - \varepsilon \f$ RANS (SSG) turbulence model.
 *
 * \param[in]     field_id      index of current field
 * \param[in]     gradv         work array for the velocity grad term
 *                                 only for iturb=31
 * \param[in]     produc        work array for production
 * \param[in]     gradro        work array for grad rom
 *                              (without rho volume) only for iturb=30
 * \param[out]    viscf         visc*surface/dist at internal faces
 * \param[out]    viscb         visc*surface/dist at edge faces
 * \param[out]    viscce        Daly Harlow diffusion term
 * \param[out]    rhs           working array
 * \param[out]    rovsdt        working array
 * \param[out]    weighf        working array
 * \param[out]    weighb        working array
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_ssg(int               field_id,
                            const cs_real_t   gradv[][3][3],
                            const cs_real_t   produc[][6],
                            const cs_real_t   gradro[][3],
                            cs_real_t         viscf[],
                            cs_real_t         viscb[],
                            cs_real_t         viscce[][6],
                            cs_real_t         rhs[][6],
                            cs_real_t         rovsdt[][6][6],
                            cs_real_t         weighf[][2],
                            cs_real_t         weighb[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Solve epsilon for \f$ R_{ij} - \varepsilon \f$ RANS
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
cs_turbulence_rij_solve_eps(cs_lnum_t        ncesmp,
                            cs_lnum_t        icetsm[],
                            int              itypsm[],
                            const cs_real_t  gradv[][3][3],
                            const cs_real_t  produc[][6],
                            const cs_real_t  gradro[][3],
                            cs_real_t        smacel[],
                            cs_real_t        viscf[],
                            cs_real_t        viscb[],
                            cs_real_t        rhs[],
                            cs_real_t        rovsdt[]);

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
