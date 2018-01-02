#ifndef __CS_JOIN_SPLIT_H__
#define __CS_JOIN_SPLIT_H__

/*============================================================================
 * Set of subroutines for cutting faces during face joining.
 *===========================================================================*/

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
 * Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_mesh.h"
#include "cs_join_set.h"
#include "cs_join_util.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *===========================================================================*/

/*============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Build new faces after the vertex fusion operation. Split initial faces into
 * subfaces and keep the historic between initial/final faces.
 *
 * parameters:
 *   param           <-- set of user-defined parameters
 *   face_normal     <-- array of normal vector on each face
 *   edges           <-- list of edges
 *   work            <-> pointer to a cs_join_mesh_t structure
 *   old2new_history <-> relation between new faces and old one:
 *                       old global face -> new local face
 *---------------------------------------------------------------------------*/

void
cs_join_split_faces(cs_join_param_t          param,
                    const cs_real_t          face_normal[],
                    const cs_join_edges_t   *edges,
                    cs_join_mesh_t         **work,
                    cs_join_gset_t         **old2new_history);

/*----------------------------------------------------------------------------
 * Update after face splitting of the local join mesh structure.
 * Send back to the original rank the new face description.
 *
 * parameters:
 *   param           <-- set of user-defined parameters
 *   work_mesh       <-- distributed mesh on faces to join
 *   gnum_rank_index <-- index on ranks for the old global face numbering
 *   o2n_hist        <-> old global face -> new local face numbering
 *   local_mesh      <-> mesh on local selected faces to be joined
 *---------------------------------------------------------------------------*/

void
cs_join_split_update_struct(const cs_join_param_t   param,
                            const cs_join_mesh_t   *work_mesh,
                            const cs_gnum_t         gnum_rank_index[],
                            cs_join_gset_t        **o2n_hist,
                            cs_join_mesh_t        **local_mesh);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_SPLIT_H__ */
