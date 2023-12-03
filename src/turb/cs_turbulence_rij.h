#ifndef __CS_TURBULENCE_RIJ_H__
#define __CS_TURBULENCE_RIJ_H__

/*============================================================================
 * Rij-epsilon turbulence model.
 *============================================================================*/

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
/*! \brief Solve the \f$ R_{ij} - \epsilon \f$ for incompressible flows or
 *         slightly compressible flows for one time step.
 *
 * Please refer to the
 * <a href="../../theory.pdf#rijeps"><b>\f$ R_{ij} - \epsilon \f$ model</b></a>
 * section of the theory guide for more informations, as well as the
 * <a href="../../theory.pdf#turrij"><b>turrij</b></a> section.
 *
 * \param[in]     phase_id     turbulent phase id (-1 for single phase flow)
 * \param[in]     ncesmp       number of cells with mass source term
 * \param[in]     icetsm       index of cells with mass source term
 * \param[in]     itypsm       mass source type for the variables
 * \param[in]     smacel       values of the variables associated to the
 *                             mass source
 *                             (for ivar=ipr, smacel is the mass flux)
 !*/
/*-----------------------------------------------------------------------------*/

void
cs_turbulence_rij(int          phase_id,
                  cs_lnum_t    ncesmp,
                  cs_lnum_t    icetsm[],
                  int          itypsm[],
                  cs_real_t    smacel[]);

/*----------------------------------------------------------------------------*/
/*! \brief Solve the equation on alpha in the framework of the Rij-EBRSM model.
 *
 * Also called for alpha of scalars for EB-DFM.
 *
 * \param[in]  f_id          field id of alpha variable
 * \param[in]  phase_id      turbulent phase id (-1 for single phase flow)
 * \param[in]  c_durbin_l    constant for the Durbin length
 !*/
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_solve_alpha(int        f_id,
                              int        phase_id,
                              cs_real_t  c_durbin_l);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize Rij-epsilon variables based on reference quantities.
 *
 * If uref is not provided (0 or negative), values are set at a large
 * negative value (-cs_math_big_r) to allow for later checks.
 *
 * \param[in]  uref    characteristic flow velocity
 * \param[in]  almax   characteristic macroscopic length of the domain
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_init_by_ref_quantities(cs_real_t  uref,
                                         cs_real_t  almax);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clip the turbulent Reynods stress tensor and the turbulent
 *        dissipation (coupled components version).
 *
 * \param[in]  phase_id   turbulent phase id (-1 for single phase flow)
 * \param[in]  n_cells    number of cells
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_clip(int        phase_id,
                       cs_lnum_t  n_cells);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the turbulent viscosity for the Reynolds Stress model.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_mu_t(int  phase_id);

/*----------------------------------------------------------------------------*/
/*! \brief Compute Rusanov equivalent diffusivity of the model.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_rij_compute_rusanov(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_RIJ_H__ */
