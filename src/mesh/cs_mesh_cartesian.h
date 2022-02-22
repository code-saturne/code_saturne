#ifndef __CS_MESH_CARTESIAN_H__
#define __CS_MESH_CARTESIAN_H__

/*============================================================================
 * Cartesian mesh generation
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_MESH_CARTESIAN_CONSTANT_LAW,
  CS_MESH_CARTESIAN_GEOMETRIC_LAW,
  CS_MESH_CARTESIAN_PARABOLIC_LAW,
  CS_MESH_CARTESIAN_USER_LAW,
  CS_MESH_CARTESIAN_N_LAW_TYPES

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
/*! \brief Create cartesian mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_create(void);

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh with a constant step in all
 *         directions
 *
 * \param[in] ncells  Array of size 3 containing number of cells in each
 *                    direction
 * \param[in] xyz     Array of size 6 containing min values, followed by
 *                    max values for the three directions.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_simple(int        ncells[3],
                                cs_real_t  xyz[6]);

/*----------------------------------------------------------------------------*/
/*! \brief Define directions parameters based on a user input
 *
 * \param[in] idir       Direction index. 0->X, 1->Y, 2->Z
 * \param[in] ncells     Number of cells for the direction
 * \param[in] vtx_coord  Array of size ncells+1 containing 1D coordinate values
 *                       for vertices on the given direction
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_user(int       idir,
                                  int       ncells,
                                  cs_real_t vtx_coord[]);

/*----------------------------------------------------------------------------*/
/*! \brief Define direction parameters based on a piecewise definition. Each
 *         part follows a geometric (or uniform) sequence. To get the uniform
 *         sequence, set the amplification factor to 1 in the wanted part.
 *
 *         A direction is split in several parts. Each part contains a number
 *         of cells, its starting and ending position (stored in a compact way)
 *         inside part_coords, the amplification factor (f) between the first
 *         and last cell size of each part. Notice that if f = 1, this leads to
 *         a uniform refinement. If f > 1, (resp f < 1) this leads to a growing
 *         (resp. decreasing) geometric progression of the cell size when
 *         moving along the direction of increasing coordinates.
 *
 * \param[in] idir          Direction index. 0->X, 1->Y, 2->Z
 * \param[in] n_parts       Number of parts to define the direction
 * \param[in] part_coords   Position delimiting each part (size = n_parts + 1)
 * \param[in] n_part_cells  Number of cells in each part (size = n_parts)
 * \param[in] amp_factors   Amplification factor in each part (size = n_parts)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_geom_by_part(int                idir,
                                          int                n_parts,
                                          const cs_real_t    part_coords[],
                                          const cs_lnum_t    n_part_cells[],
                                          const cs_real_t    amp_factors[]);

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh based on a CSV file.
 *         CSV file needs to contain :
 *         (1) First line which is empty or contains a header
 *         (2) Second line containing number of vertices per direction:
 *             format is 'nx;ny;nz'
 *         (3) Third line is empty or contains a header
 *         (4) Fourth line and onwards contains vertices coordinates for each
 *             direction. Format is "X1[i];X2[i];X3[i]" for index i.
 *             If current vertex index is beyond max value for a given
 *             direction, an empty value is expected.
 *             For example, if for index 'j' the first direction is empty,
 *             format is : ';X2[j];X3[j]'
 *
 * \param[in] csv_file_name  name of CSV file containing mesh information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_from_csv(const char *csv_file_name);

/*----------------------------------------------------------------------------*/
/*! \brief Define parameters for a given direction.
 *
 * \param[in] idim         Geometrical direction: 0->X, 1->Y, 2->Z
 * \param[in] law          1D discretization law: constant, geometric or
 *                         parabolic
 * \param[in] ncells       Number of cells for this direction
 * \param[in] smin         Min coordinate value for this direction
 * \param[in] smax         Max coordinate value for this direction
 * \param[in] progression  Progression value, only used for geometric or
 *                         parabolic laws.
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
 * \param[in] m     pointer to cs_mesh_t structure
 * \param[in] mb    pointer to cs_mesh_builder_t structure
 * \param[in] echo  verbosity flag
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_connectivity(cs_mesh_t          *m,
                               cs_mesh_builder_t  *mb,
                               long                echo);

/*----------------------------------------------------------------------------*/
/*! \brief Destroy cartesian mesh parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_params_destroy(void);

/*----------------------------------------------------------------------------*/

#endif /* __CS_MESH_CARTESIAN_H__ */
