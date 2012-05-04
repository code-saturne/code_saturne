#ifndef __CS_MESH_BAD_CELLS_DETECTION_H__
#define __CS_MESH_BAD_CELLS_DETECTION_H__

/*============================================================================
 *  Detect and post-process bad cells within meshes.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute and post-process mesh quality indicators.
 *
 * parameters:
 *   mesh             --> pointer to a mesh structure.
 *   mesh_quantities  <-> pointer to a mesh quantities structures.
 *----------------------------------------------------------------------------*/

void
cs_mesh_bad_cells_detection(const cs_mesh_t       *mesh,
                            cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BAD_CELLS_DETECTION_H__ */
