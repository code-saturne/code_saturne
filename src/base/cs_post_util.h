#ifndef __CS_POST_UTIL_H__
#define __CS_POST_UTIL_H__

/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type Definitions
 *============================================================================*/

typedef enum {

  CS_POST_UTIL_Q_CRITERION,        /*!< Q-criterion post-treatment */

  CS_POST_UTIL_N_TYPES             /*!< Number of post utility types */

} cs_post_util_type_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Names of post utilities */

extern int cs_glob_post_util_flag[];

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the head of a turbomachinery (total pressure increase)
 *
 * \param[in]   criteria_in   selection criteria of turbomachinery suction
 * \param[in]   location_in   mesh location of turbomachinery suction
 * \param[in]   criteria_out  selection criteria of turbomachinery discharge
 * \param[in]   location_out  mesh location of turbomachinery discharge
 *
 * \return turbomachinery head
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_turbomachinery_head(const char               *criteria_in,
                            cs_mesh_location_type_t   location_in,
                            const char               *criteria_out,
                            cs_mesh_location_type_t   location_out);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the magnitude of a moment of force torque) given an
 *         axis and the stress on a specific boundary.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_list  face list
 * \param[in]   axis         axis
 *
 * \return couple about the axis
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_post_moment_of_force(cs_lnum_t     n_b_faces,
                        cs_lnum_t    *b_face_list,
                        cs_real_3_t   axis);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute tangential stress on a specific boundary.
 *
 * \param[in]   n_b_faces    number of faces
 * \param[in]   b_face_list  face list
 * \param[out]  stress       tangential stress on the specific
 *                           boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_post_stress_tangential(cs_lnum_t n_b_faces,
                          cs_lnum_t  b_face_list[],
                          cs_real_3_t stress[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute Reynolds stresses in case of Eddy Viscosity Models
 *
 * \param[in]  n_loc_cells   number of cells
 * \param[in]  cells_list    cells list
 * \param[out] rst           Reynolds stresses stored as vector
 *                           [r11,r22,r33,r12,r23,r13]
 */
/*----------------------------------------------------------------------------*/

void
cs_post_evm_reynolds_stresses(cs_lnum_t   n_loc_cells,
                              cs_lnum_t   cells_list[],
                              cs_real_6_t rst[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the Q criterion from Hunt et. al over a specified
 *        volumic region.
 *
 * \param[in]  n_loc_cells   number of cells
 * \param[in]  cells_list    cells list
 * \param[out] q_crit        Q-criterion over the specified volumic region.
 */
/*----------------------------------------------------------------------------*/

void
cs_post_q_criterion(const cs_lnum_t   n_loc_cells,
                    const cs_lnum_t   cells_list[],
                          cs_real_t   q_crit[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_UTIL_H__ */
