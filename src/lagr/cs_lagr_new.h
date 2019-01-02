#ifndef __CS_LAGR_NEW_H__
#define __CS_LAGR_NEW_H__

/*============================================================================
 * Handling of new particles.
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_lagr_particle.h"

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
 * \brief Inject a series of particles at random positions on given faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[in,out]  particles          pointer to particle set
 * \param[in]      n_faces            number of faces in zone
 * \param[in]      face_ids           ids of faces in zone
 * \param[in]      face_particle_idx  starting index of added particles
 *                                    for each face in zone
 *----------------------------------------------------------------------------*/

void
cs_lagr_new(cs_lagr_particle_set_t  *particles,
            cs_lnum_t                n_faces,
            const cs_lnum_t          face_ids[],
            const cs_lnum_t          face_particle_idx[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject a series of particles at random positions on given cells.
 *
 * \warning Currently works only for tri and quadrangular faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[in,out]  particles          pointer to particle set
 * \param[in]      n_cells            number of cells in zone
 * \param[in]      cell_ids           ids of cells in zone
 * \param[in]      cell_particle_idx  starting index of added particles
 *                                    for each cell in zone
 *----------------------------------------------------------------------------*/

void
cs_lagr_new_v(cs_lagr_particle_set_t  *particles,
              cs_lnum_t                n_cells,
              const cs_lnum_t          cell_ids[],
              const cs_lnum_t          cell_particle_idx[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization for new particles.
 *
 * The fluid velocity seen is computed here.
 *
 * \param[in]  particle_range  start and past-the-end ids of new particles
 *                             for this zone and class
 * \param[in]  time_id         associated time id (0: current, 1: previous)
 * \param[in]  visc_length     viscous layer thickness
 *                             (size: number of mesh boundary faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_particle_init(const cs_lnum_t  particle_range[2],
                          int              time_id,
                          const cs_real_t  visc_length[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_NEW_H__ */
