#ifndef __CS_MESH_BOUNDARY_LAYER_H__
#define __CS_MESH_BOUNDARY_LAYER_H__

/*============================================================================
 * Insert boundary cell layers into the mesh.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 * \param[in]  n_zones         number of zones
 * \param[in]  zone_ids        ids of zones
 * \param[in]  zone_layers     number of extrusion layers per zone, or
 *                             NULL for default (1)
 * \param[in]  zone_thickness  thickness for each zone, or NULL for
 *                             default (based on neighboring cell
 *                             sizes)
 * \param[in]  zone_expansion  parameter for each zone, or NULL for
 *                             default (0.8)
 * \param[in]  n_layers        optional specification of number of layers for
 *                             each vertex, or NULL
 * \param[in]  distribution    optional specification of layer distribution
 *                             for each vertex, or NULL
 * \param[in]  coord_shift     optional definition of coordinate shift for
 *                             each vertex, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_layer_insert(int                  n_zones,
                              const int            zone_ids[],
                              const int            zone_layers[],
                              const double         zone_thickness[],
                              const float          zone_expansion[],
                              const cs_lnum_t     *n_layers,
                              const float         *distribution,
                              const cs_coord_3_t  *coord_shift);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BOUNDARY_LAYER_H__ */
