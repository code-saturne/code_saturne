#ifndef __CS_MESH_BOUNDARY_LAYER_H__
#define __CS_MESH_BOUNDARY_LAYER_H__

/*============================================================================
 * Insert boundary cell layers into the mesh.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_defs.h"

#include "cs_mesh.h"
#include "cs_mesh_extrude.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert mesh boundary layers.
 *
 * \param[in, out]  m                  mesh
 * \param[in, out]  e                  extrusion vector definitions
 * \param[in]       min_volume_factor  cell volume multiplier threshold:
 *                                     extrusion is reduced on vertices
 *                                     adjacent to cells whose volume is
 *                                     reduced below this; < 0 to ignore
 * \param[in]       interior_gc        if true, maintain group classes of
 *                                     interior faces previously on boundary
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_layer_insert(cs_mesh_t                  *m,
                              cs_mesh_extrude_vectors_t  *e,
                              cs_real_t                   min_volume_factor,
                              bool                        interior_gc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BOUNDARY_LAYER_H__ */
