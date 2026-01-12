#ifndef __CS_MESH_COARSEN_H__
#define __CS_MESH_COARSEN_H__

/*============================================================================
 * Mesh coarsening.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_base.h"
#include "mesh/cs_mesh.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine rank which should be assigned to given cells to allow
 *        mesh coarsening.
 *
 * All sub-cells of a given parent cell should be on the same rank to allow
 * coarsening, so this updates a destination rank array to ensure this.
 *
 * The dest_rank array should be initialized using a prior partitioning
 * step, or set to the current local rank id.
 *
 * \param[in]       m          mesh
 * \param[in]       cell_flag  coarsening flag for each cell (0: no 1: yes),
 *                             or null (in which case 1 is assumed uniformly)
 * \param[in, out]  dest_rank  cell initial refinement level if coarsened,
 *                             0 otherwise.
 *
 * \return  number of new cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_cells_dest_rank(cs_mesh_t   *m,
                                const int   *cell_flag,
                                int          dest_rank[]);

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \brief Coarsen flagged mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       cell_flag   coarsening flag for each cell
 *                              (0: do not coarsen; 1: coarsen)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_simple(cs_mesh_t  *m,
                       const int   cell_flag[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Coarsen selected mesh cells.
 *
 * \param[in, out]  m           mesh
 * \param[in]       n_cells     number of selected cells
 * \param[in]       cells       list of selected cells (0 to n-1)
 *                              or NULL if no indirection is needed
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_coarsen_simple_selected(cs_mesh_t        *m,
                                cs_lnum_t         n_cells,
                                const cs_lnum_t   cells[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_COARSEN_H__ */
