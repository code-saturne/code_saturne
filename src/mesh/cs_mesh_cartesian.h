#ifndef __CS_MESH_CARTESIAN_H__
#define __CS_MESH_CARTESIAN_H__

/*============================================================================
 * Cartesian mesh generation
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*============================================================================
 * Structure definitions
 *============================================================================*/

typedef enum {

  CS_MESH_CARTESIAN_CONSTANT_LAW,
  CS_MESH_CARTESIAN_GEOMETRIC_LAW,
  CS_MESH_CARTESIAN_PARABOLIC_LAW,
  CS_MESH_CARTESIAN_USER_LAW,
  CS_N_MESH_CARTESIAN_LAW_TYPE

} cs_mesh_cartesian_law_t;

typedef struct _cs_mesh_cartesian_params_t cs_mesh_cartesian_params_t;
/*============================================================================
 * Public C function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to cartesian mesh parameters structure
 *
 * \return pointer to cs_mesh_cartesian_params_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_get_params(void);

/*----------------------------------------------------------------------------*/
/*! \brief Create cartesian mesh strucutre
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_create(void);

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh with a constant step in all directions
 *
 * \param[in] ncells  Array of size 3 containing number of cells in each direction
 * \param[in] xyz     Array of size 6 containing min values, followed by
 *                    max values for the three directions.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_simple(int        ncells[3],
                                cs_real_t  xyz[6]);

/*----------------------------------------------------------------------------*/
/*! \brief Define parameters for a given direction.
 *
 * \param[in] idim         Geometrical direction: 0->X, 1->Y, 2->Z
 * \param[in] law          1D discreization law : constant, geometric or parabolic
 * \param[in] ncells       Number of cells for this direction
 * \param[in] smin         Min coordinate value for this direction
 * \param[in] smax         Max coordinate value for this direction
 * \param[in] progression  Progression value, only used for geometric or
 *                         parabolic laws.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_params(int                     idim,
                                    cs_mesh_cartesian_law_t law,
                                    int                     ncells,
                                    cs_real_t               smin,
                                    cs_real_t               smax,
                                    cs_real_t               progression);

/*----------------------------------------------------------------------------*/
/*! \brief Indicate if a cartesian mesh is to be built.
 *
 * \return 1 if mesh needs to be built, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_need_build(void);

/*----------------------------------------------------------------------------*/
/*! \brief Get number of cells in a given direction.
 *
 * \param[in] idim  Index of direction: 0->X, 1->Y, 2->Z
 *
 * \return Number of cells in corresponding direction (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_ncells(int idim);

/*----------------------------------------------------------------------------*/
/*! \brief Build unstructured connectivity needed for partitionning.
 *
 * \param[in] mb    pointer to cs_mesh_builder_t structure
 * \param[in] echo  verbosity flag
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_connectivity(cs_mesh_builder_t *mb,
                               long               echo);

/*----------------------------------------------------------------------------*/
/*! \brief Destroy cartesian mesh parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_params_destroy(void);

/*----------------------------------------------------------------------------*/

#endif /* __CS_MESH_CARTESIAN_H__ */
