#ifndef __CS_MESH_QUALITY_H__
#define __CS_MESH_QUALITY_H__

/*============================================================================
 * Compute several mesh quality criteria.
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Evaluate face warping angle for internal and border faces..
 *
 * parameters:
 *   mesh             --> pointer to a cs_mesh_t structure
 *   i_face_normal    --> internal face normal
 *   b_face_normal    --> border face normal
 *   i_face_warping   <-- face warping angle for internal faces
 *   b_face_warping   <-- face warping angle for border faces
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
cs_mesh_quality_compute_warping(const cs_mesh_t    *mesh,
                                const cs_real_t     i_face_normal[],
                                const cs_real_t     b_face_normal[],
                                cs_real_t           i_face_warping[],
                                cs_real_t           b_face_warping[]);

/*----------------------------------------------------------------------------
 * Compute mesh quality indicators
 *
 * parameters:
 *   mesh             --> pointer to a mesh structure.
 *   mesh_quantities  --> pointer to a mesh quantities structures.
 *----------------------------------------------------------------------------*/

void
cs_mesh_quality(const cs_mesh_t             *mesh,
                const cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_QUALITY_H__ */
