#ifndef __CS_LAGR_FRAGMENTATION_H__
#define __CS_LAGR_FRAGMENTATION_H__

/*============================================================================
 * Fragmentation modeling.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Particle fragmentation algorithm.
 *
 * Parcels represent physical particles with similar properties (size).
 * The number of physical particles in a parcel is represented by the
 * statistical weight.
 *
 * For each parcel, the number of fragmentation events is generated with
 * a Poisson distribution that depends on the fragmentation kernel.
 *
 * Working hypotheses:
 *  1) Discrete diameters
 *     - Minimal particle size , called a monomere (unbreakable particle)
 *  2) Fragmentation is binary
 *     - Equally sized fragments (N0/2, or (N0+-1)/2)
 *  3) Fragmentation is treated after agglomeration
 *
 * Warning:  particle indices are not necessarily contiguous
 *           (due to agglomeration occuring before).
 *
 * Two contiguous groups: particles present before agglomeration
 *                        newly created particles (due to agglomeration)
 *
 * \param[in]  cell_id                current cell id
 * \param[in]  dt                     time step
 * \param[in]  minimum_particle_diam  minumum diameter (monomere diameter)
 * \param[in]  rho                    particles density
 * \param[in]  main_start             index of the first particle in cell
 *                                    present before the agglomeration
 * \param[in]  main_end               index after the last particle in cell,
 *                                    present before the agglomeration
 * \param[in]  agglo_start            index of the first particle in cell,
 *                                    created by the agglomeration
 * \param[in]  agglo_end              index after the last particle in cell,
 *                                    created by the agglomeration
 *
 * \returns  modified list of particles containing newly created parcels
 *           at the end
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_fragmentation(cs_real_t  dt,
                      cs_real_t  minimum_particle_diam,
                      cs_real_t  rho,
                      cs_lnum_t  main_start,
                      cs_lnum_t  main_end,
                      cs_lnum_t  agglo_start,
                      cs_lnum_t  agglo_end);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_FRAGMENTATION_H__*/
