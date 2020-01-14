#ifndef __CS_LAGR_COUPLING_H__
#define __CS_LAGR_COUPLING_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * \brief  Compute source terms for Lagrangian 2-way coupling.
 *
 * \remark  Source terms are computed for the starting cell of a particle
 *          during a given iteration. Even if particle exits the domain,
 *          it s necessary to compute a source term matching the exchange
 *          between the carrier fluid and the particle at the beginning
 *          of the time step. If cs_glob_lagr_time_step->nor == 2 and the
 *          particle interacts with a boundary, then the source terms
 *          are computed as if nor == 1.
 *
 * \param[in]   taup    dynamic characteristic time
 * \param[in]   tempct  thermal charactersitic time
 * \param[out]  tsfext  external forces
 * \param[in]   cpgd1   devolatization term 1 for heterogeneous coal
 * \param[in]   cpgd2   devolatization term 2 for heterogeneous coal
 * \param[in]   cpght   combustion term for heterogeneous coal
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling(const cs_real_t  taup[],
                 const cs_real_t  tempct[],
                 cs_real_t        tsfext[],
                 const cs_real_t  cpgd1[],
                 const cs_real_t  cpgd2[],
                 const cs_real_t  cpght[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_COUPLING_H__ */
