#ifndef CS_RENUMBER_UPDATE_H
#define CS_RENUMBER_UPDATE_H

/*============================================================================
 * Repartitioning and redistribution mesh data and fields.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Renumber fields and bc_coeffs based on new-to-old element map.
 *
 * \param[in]  cell_n2o        cells new-to-old mapping
 * \param[in]  i_face_n2o      internal faces new-to-old mapping
 * \param[in]  b_face_n2o      boundary faces new-to-old mapping
 * \param[in]  vtx_n2o         vertices new-to-old mapping
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_update_fields(const cs_lnum_t   cell_n2o[],
                          const cs_lnum_t   i_face_n2o[],
                          const cs_lnum_t   b_face_n2o[],
                          const cs_lnum_t   vtx_n2o[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Renumber bc_types based on new-to-old element map.
 *
 * \param[in]  b_face_n2o      boundary faces new-to-old mapping
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_update_bc_types(const cs_lnum_t   b_face_n2o[]);


/*----------------------------------------------------------------------------*/
/*
 * \brief Renumber mesh and update all associated data.
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_update(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Update some mesh structures.
 *
 * This function regroups all the data that needs to be updated after some
 * mesh operations, such as renumbering or redistribution.
 *
 * \param[in]       mesh            pointer to mesh
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_update_partial(void);


/*----------------------------------------------------------------------------*/

#endif /* CS_RENUMBER_UPDATE_H */
