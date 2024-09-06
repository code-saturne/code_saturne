#ifndef __CS_MESH_CARTESIAN_H__
#define __CS_MESH_CARTESIAN_H__

/*============================================================================
 * Cartesian mesh generation
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_MESH_CARTESIAN_CONSTANT_LAW,  /* Constant step law */
  CS_MESH_CARTESIAN_GEOMETRIC_LAW, /* Geometric step law */
  CS_MESH_CARTESIAN_PARABOLIC_LAW, /* Parabolic step law */
  CS_MESH_CARTESIAN_USER_LAW,      /* User defined discretization */
  CS_MESH_CARTESIAN_N_LAW_TYPES    /* Number of step discretization laws */

} cs_mesh_cartesian_law_t;

typedef struct _cs_mesh_cartesian_params_t cs_mesh_cartesian_params_t;

/*============================================================================
 * Public C function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return number of structured meshes to build.
 *
 * \returns number of structured meshes to build.
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_number_of_meshes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to cartesian mesh parameters structure
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \return pointer to cs_mesh_cartesian_params_t structure
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_id(const int id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get function for structured mesh based on its name.
 *
 * \param[in] name  Name of mesh
 *
 * \returns pointer to corresponding mesh parameters, or NULL if mesh
 *          does not exist.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_name_try(const char *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get function for structured mesh based on its name.
 *
 * \param[in] name  Name of mesh
 *
 * \returns pointer to corresponding mesh parameters. Raises error if mesh does
 * not exist.
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_by_name(const char *name);

/*----------------------------------------------------------------------------*/
/*! \brief Create cartesian mesh structure
 *
 * \param[in] name  Name of mesh to create
 *
 * \returns pointer to newly created mesh parameters
 */
/*----------------------------------------------------------------------------*/

cs_mesh_cartesian_params_t *
cs_mesh_cartesian_create(const char *name);

/*----------------------------------------------------------------------------*/
/*! \brief Define a simple cartesian mesh with a constant step in all
 *         directions
 *
 * \param[in] name    Name of mesh to create
 * \param[in] ncells  Array of size 3 containing number of cells in each
 *                    direction
 * \param[in] xyz     Array of size 6 containing min values, followed by
 *                    max values for the three directions.
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_define_simple(const char *name,
                                int         ncells[3],
                                cs_real_t   xyz[6]);

/*----------------------------------------------------------------------------*/
/*! \brief Define directions parameters based on a user input
 *
 * \param[in] mp         Pointer to mesh parameters
 * \param[in] idir       Direction index. 0->X, 1->Y, 2->Z
 * \param[in] ncells     Number of cells for the direction
 * \param[in] vtx_coord  Array of size ncells+1 containing 1D coordinate values
 *                       for vertices on the given direction
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_user(cs_mesh_cartesian_params_t *mp,
                                  int                         idir,
                                  int                         ncells,
                                  cs_real_t                   vtx_coord[]);

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
 * \param[in] mp            Pointer to mesh parameters
 * \param[in] idir          Direction index. 0->X, 1->Y, 2->Z
 * \param[in] n_parts       Number of parts to define the direction
 * \param[in] part_coords   Position delimiting each part (size = n_parts + 1)
 * \param[in] n_part_cells  Number of cells in each part (size = n_parts)
 * \param[in] amp_factors   Amplification factor in each part (size = n_parts)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_dir_geom_by_part(cs_mesh_cartesian_params_t *mp,
                                          int                         idir,
                                          int                         n_parts,
                                          const cs_real_t  part_coords[],
                                          const cs_lnum_t  n_part_cells[],
                                          const cs_real_t  amp_factors[]);

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
 * \param[in] name           Name of new mesh
 * \param[in] csv_file_name  name of CSV file containing mesh information.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_define_from_csv(const char *name,
                                  const char *csv_file_name);

/*----------------------------------------------------------------------------*/
/*! \brief Define parameters for a given direction.
 *
 * \param[in] mp           Pointer to mesh parameters
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
cs_mesh_cartesian_define_dir_params(cs_mesh_cartesian_params_t *mp,
                                    int                         idim,
                                    cs_mesh_cartesian_law_t     law,
                                    int                         ncells,
                                    cs_real_t                   smin,
                                    cs_real_t                   smax,
                                    cs_real_t                   progression);

/*----------------------------------------------------------------------------*/
/*! \brief Indicate if a cartesian mesh is to be built.
 *
 * \return 1 if mesh needs to be built, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_need_build(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get name of structured mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns Name of the mesh
 */
/*----------------------------------------------------------------------------*/

const char *
cs_mesh_cartesian_get_name(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get group class id shift of cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns shift value
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_gc_id_shift(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set group class id shift of cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] shift Value of shift index
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_set_gc_id_shift(int  id,
                                  int  shift);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of cells of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of cells
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_cells(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of faces of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of faces
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_faces(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get global number of vertices of a cartesian mesh
 *
 * \param[in] id    Id of the cartesian mesh
 *
 * \returns number of vertices
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_cartesian_get_n_g_vtx(int  id);

/*----------------------------------------------------------------------------*/
/*! \brief Get number of cells in a given direction.
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] idim  Index of direction: 0->X, 1->Y, 2->Z
 *
 * \return Number of cells in corresponding direction (int)
 */
/*----------------------------------------------------------------------------*/

int
cs_mesh_cartesian_get_ncells(int  id,
                             int  idim);

/*----------------------------------------------------------------------------*/
/*! \brief Build unstructured connectivity needed for partitionning.
 *
 * \param[in] id    Id of the cartesian mesh
 * \param[in] m     pointer to cs_mesh_t structure
 * \param[in] mb    pointer to cs_mesh_builder_t structure
 * \param[in] echo  verbosity flag
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_block_connectivity(int                 id,
                                     cs_mesh_t          *m,
                                     cs_mesh_builder_t  *mb,
                                     long                echo);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute all global values for meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_finalize_definition(void);

/*----------------------------------------------------------------------------*/
/*! \brief Destroy cartesian mesh parameters
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_params_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set maximum number of cartesian blocks (by default is set to None)
 *
 * \param[in] n_blocks  maximum number of cartesian blocks which can be created
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_cartesian_set_max_number_of_blocks(int  n_blocks);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_CARTESIAN_H__ */
