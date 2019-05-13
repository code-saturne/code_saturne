#ifndef __CS_MESH_BAD_CELLS_H__
#define __CS_MESH_BAD_CELLS_H__

/*============================================================================
 * Detect bad cells within meshes.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * Detection criteria type
 */

#define CS_BAD_CELL_ORTHO_NORM           (1 << 0)
#define CS_BAD_CELL_OFFSET               (1 << 1)
#define CS_BAD_CELL_LSQ_GRAD             (1 << 2)
#define CS_BAD_CELL_RATIO                (1 << 3)
#define CS_BAD_CELL_GUILT                (1 << 4)
#define CS_BAD_CELL_USER                 (1 << 5)
#define CS_BAD_CELL_TO_REGULARIZE        (1 << 6)

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define which cell quality indicators are used and when.
 *
 * Note: we assume that if a given criterion is computed at each time
 * step, it is also computed at initialization, but for visualization,
 * it is either one or the other, as visualization formats and tools
 * may not always accept both a fixed and time-varying instance of a
 * given variable.
 *
 * parameters:
 *   type_flag_mask <-- criterion type mask (0 for all)
 *   compute        <-- 0: never compute;
 *                      1: compute at initialization;
 *                      2: compute at each time step
 *   visualize      <-- 0: never visualize
 *                      1: visualize at initialization;
 *                      2: visualize at each time step
 *----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_set_options(int  type_flag_mask,
                              int  compute,
                              int  visualize);

/*----------------------------------------------------------------------------
 * Indicate which cell quality indicators are used and when.
 *
 * Each array is optional, and returns 2 flags; the first flag is used at
 * initialization, the second one at each time step.
 *
 * A flag is a mask to be compared using an "and" (&) operation with a given
 * criteria type mask (CS_BAD_CELL_ORTHO_NORM, CS_BAD_CELL_OFFSET, ...).
 *
 * parameters:
 *   compute   --> computation mask (initialization, per time step), or NULL
 *   visualize --> visualization mask (initialization, per time step), or NULL
 *----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_get_options(int  compute[2],
                              int  visualize[2]);

/*----------------------------------------------------------------------------
 * Compute bad cell quality indicators.
 *
 * parameters:
 *   mesh            <-- pointer to a mesh structure.
 *   mesh_quantities <-> pointer to a mesh quantities structures.
 *----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_detect(const cs_mesh_t       *mesh,
                         cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------
 * Post-process bad cell quality indicators.
 *
 * parameters:
 *   mesh            <-- pointer to a mesh structure.
 *   mesh_quantities <-- pointer to a mesh quantities structures.
 *----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_postprocess(const cs_mesh_t             *mesh,
                              const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BAD_CELLS_H__ */
