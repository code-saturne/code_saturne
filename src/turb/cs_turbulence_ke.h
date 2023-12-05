#ifndef __CS_TURBULENCE_KE_H__
#define __CS_TURBULENCE_KE_H__

/*============================================================================
 * k-epsilon turbulence model.
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
/*!
 * \brief Solve the k-epsilon equations.
 *
 * Solve the \f$ k - \varepsilon \f$  for incompressible flows
 * or slightly compressible flows for one time step.
 *
 * \param[in]     phase_id      turbulent phase id (-1 for single phase flow)
 * \param[in]     ncesmp        number of cells with mass source term
 * \param[in]     icetsm        index of cells with mass source term
 * \param[in]     itypsm        mass source type for the variables
 *                              size: [nvar][ncesmp]
 * \param[in]     dt            time step (per cell)
 * \param[in]     smacel        values of the variables associated to the
 *                              mass source (for the pressure variable,
 *                              smacel is the mass flux)
 *                              size: [nvar][ncesmp]
 * \param[out]    prdv2f        v2f production term
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke(int              phase_id,
                 cs_lnum_t        ncesmp,
                 cs_lnum_t        icetsm[],
                 int              itypsm[],
                 const cs_real_t *dt,
                 cs_real_t        smacel[],
                 cs_real_t       *prdv2f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clipping of the turbulent kinetic energy and turbulent dissipation.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 * \param[in]     n_cells  number of cells
 * \param[in]     iclip    indicator = 0 if viscl0 is used
 *                         otherwise viscl is used.
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_clip(int        phase_id,
                      cs_lnum_t  n_cells,
                      int        iclip);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for the K-epsilon model.
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_mu_t(int  phase_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of turbulent viscosity for
 *        the non-linear quadratic K-epsilon from
 *        Baglietto et al. (2005)
 *
 * \param[in]     phase_id   turbulent phase id (-1 for single phase flow)
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_q_mu_t(int  phase_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Calculation of non linear terms of the quadratic k-epsilon model
 *        (Baglietto et al.)
 *
 * \param[in]   phase_id  turbulent phase id (-1 for single phase flow)
 * \param[out]  rij        non linear terms of quadratic Boussinesq approximation
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_ke_q(int          phase_id,
                   cs_real_6_t  rij[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TURBULENCE_KE_H__ */
