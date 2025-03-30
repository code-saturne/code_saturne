#ifndef __CS_SOLVE_NAVIER_STOKES_H__
#define __CS_SOLVE_NAVIER_STOKES_H__

/*============================================================================
 * Solve the Navier-Stokes equations.
 *============================================================================*/

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Update total pressure (defined as a post-processed property).
 *
 * For the compressible module, the solved pressure is already
 * the total pressure.
 *
 * Note: for Eddy Viscosity Models, the TKE may be included in the
 * solved pressure.
 *
 * \param[in]   m       pointer to mesh structure
 * \param[in]   mq      pointer to mesh quantities structure
 * \param[in]   fp      pointer to fluid properties structure
 * \param[in]   gxyz_h  gravity (on host)
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_navier_stokes_update_total_pressure
  (const cs_mesh_t              *m,
   const cs_mesh_quantities_t   *mq,
   const cs_fluid_properties_t  *fp,
   const cs_real_t               gxyz_h[3]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Solve Navier-Stokes equations for incompressible or slightly
 *        compressible flows for one time step. Both convection-diffusion
 *        and continuity steps are performed.
 *
 * \param[in]       iterns        index of the iteration on Navier-Stokes
 * \param[in]       icvrge        convergence indicator
 * \param[in]       itrale        number of the current ALE iteration
 * \param[in]       isostd        indicator of standard outlet
 *                              + index of the reference face
 * \param[in]       ckupdc        head loss coefficients, if present
 * \param[in, out]  trava         working array for velocity-pressure coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_solve_navier_stokes(const int        iterns,
                       int             *icvrge,
                       const int        itrale,
                       const int        isostd[],
                       const cs_real_t  ckupdc[][6],
                       cs_real_3_t     *trava);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOLVE_NAVIER_STOKES_H__ */
