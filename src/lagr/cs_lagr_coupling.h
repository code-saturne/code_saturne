#ifndef __CS_LAGR_COUPLING_H__
#define __CS_LAGR_COUPLING_H__

/*============================================================================
 * Functions and types for the Lagrangian module
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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Prepare source terms for Lagrangian 2-way coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Increment the source terms for Lagrangian 2-way coupling with
 *         quantites attached to a given particle.
 *
 * \remark  Source terms are computed for the starting cell of a particle
 *          during a given iteration. Even if particle exits the domain,
 *          it s necessary to compute a source term matching the exchange
 *          between the carrier fluid and the particle at the beginning
 *          of the time step. If nor == 2 and the particle interacts with a
 *          boundary, then the source terms are computed as if nor == 1.
 *
 * \param[in]   p_set   pointer to particle set
 * \param[in]   p_id    particle id
 * \param[in]   dt_part remaining time step associated to the particle
 * \param[in]   rebound true if a rebound occured over last trajectory step
 * \param[in]   taup    dynamic characteristic time
 * \param[in]   force_p forces per mass unit on particles (m/s^2)
 * \param[in]   tempct  thermal characteristic time
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_increment_part_contrib(cs_lagr_particle_set_t       *p_set,
                                        const cs_lnum_t               p_id,
                                        const cs_real_t               dt_part,
                                        const bool                    rebound,
                                        const cs_real_t               taup,
                                        const cs_real_3_t             force_p,
                                        const cs_real_2_t             tempct);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize source terms for Lagrangian 2-way coupling by treating
 *         cell-attached properties.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_COUPLING_H__ */
