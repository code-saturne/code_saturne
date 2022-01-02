#ifndef __CS_MESH_INTERSECT_H__
#define __CS_MESH_INTERSECT_H__

/*============================================================================
 * Postprocessing utility functions.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Type Definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a given segment
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \param[in]   input     pointer to segment start and end:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_intersect_segment_cell_select(void        *input,
                                      cs_lnum_t   *n_cells,
                                      cs_lnum_t  **cell_ids);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select cells cut by a line composed of segments
 *
 * This selection function may be used as an elements selection function
 * for postprocessing.
 *
 * In this case, the input points to a real array containing the segment's
 * start and end coordinates.
 *
 * Note: the input pointer must point to valid data when this selection
 * function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - post-processing output must be ensured using cs_post_write_meshes()
 *   with a fixed-mesh writer before the data pointed to goes out of scope;
 *
 * The caller is responsible for freeing the returned cell_ids array.
 * When passed to postprocessing mesh or probe set definition functions,
 * this is handled automatically.
 *
 * \param[in]   input     pointer to segments starts and ends:
 *                        [x0, y0, z0, x1, y1, z1]
 * \param[in]   n_points  number of vertices in the polyline
 * \param[out]  n_cells   number of selected cells
 * \param[out]  cell_ids  array of selected cell ids (0 to n-1 numbering)
 * \param[out]  seg_c_len array of length of the segment in the selected cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_intersect_polyline_cell_select(void        *input,
                                       cs_lnum_t    n_points,
                                       cs_lnum_t   *n_cells,
                                       cs_lnum_t  **cell_ids,
                                       cs_real_t  **seg_c_len);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_INTERSECT_H__ */
