#ifndef __CS_LAGR_AGGLO_H__
#define __CS_LAGR_AGGLO_H__

/*============================================================================
 * Agglomeration modeling.
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
 * \brief Merge two sorted arrays in a third sorted array
 *
 * \param[in]       arr1   first sorted array
 * \param[in]       arr2   second sorted array
 * \param[in]       n1     size of first sorted array
 * \param[in]       n2     size of second sorted array
 * \param[in, out]  arr3   preallocated array that will contain the sorted
 *                         merge of the two previous arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_agglo_merge_arrays(cs_lnum_2_t  arr1[],
                           cs_lnum_2_t  arr2[],
                           cs_lnum_t    n1,
                           cs_lnum_t    n2,
                           cs_lnum_2_t  arr3[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Agglomeration algorithm based on algorithms used in
 *        rare gas modelling.
 *
 * Parcels represent physical particles with similar properties (size).
 * The number of physical particles in a parcel is represented by
 * the statistical weight.
 *
 * - We select randomly a number of parcels within which agglomeration
 *   is treated.
 * - We selected randomly pairs of particles and, for each pair,
 *   the number of agglomeration events is generated with
 *   a Poisson distribution that depends on the agglomeration kernel and
 *   number of each particle (i.e. statistical weight).
 *
 * Working hypotheses:
 *
 * 1) Discrete diameters
 *    - Minimal particle size , called a monomere (unbreakable particle)
 *    - Aggregates correponds to group of monomeres sticking together.
 *    - Each class of the parcel represent one size (stored in CS_LAGR_AGGREGATE)
 * 2) Agglomeration happens between two parcels
 * 3) Particles in the same cell are contiguous in the particle list
 *
 * \param[in]  cell_id                current cell id
 * \param[in]  dt                     time step
 * \param[in]  minimum_particle_diam  minumum diameter (monomere diameter)
 * \param[in]  rho                    particles density
 * \param[in]  start_particle         index of the first particle
 * \param[in]  end_particle           index after the last particle
 *
 * \returns a modified list of particles, containing newly
 *          created parcels at the end
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_lagr_agglomeration(cs_lnum_t  cell_id,
                      cs_real_t  dt,
                      cs_real_t  minimum_particle_diam,
                      cs_real_t  rho,
                      cs_lnum_t  start_particle,
                      cs_lnum_t  end_particle);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_AGGLO_H__*/
