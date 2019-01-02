#ifndef __CS_MESH_SMOOTHER_H__
#define __CS_MESH_SMOOTHER_H__

/*============================================================================
 * Mesh smoothing.
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set fixed vertices flag based on feature angle criterion.
 *
 * parameters:
 *   mesh           <--  pointer to a cs_mesh_t structure
 *   feature_angle  <--  feature angle (bounded between 0 and 90 degrees)
 *   vtx_is_fixed   -->  array to define vertices mobility (1: fixed, 0: free)
 *----------------------------------------------------------------------------*/

void
cs_mesh_smoother_fix_by_feature(cs_mesh_t   *mesh,
                                cs_real_t    feature_angle,
                                int          vtx_is_fixed[]);

/*----------------------------------------------------------------------------
 * Unwarping smoother.
 *
 * parameters:
 *   mesh         <-- pointer to a cs_mesh_t structure
 *   vtx_is_fixed --> array to define vertices mobility (1 : fixed, 0 : free)
 *----------------------------------------------------------------------------*/

void
cs_mesh_smoother_unwarp(cs_mesh_t  *mesh,
                        const int   vtx_is_fixed[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_SMOOTHER_H__ */
