#ifndef __CS_JOIN_MERGE_H__
#define __CS_JOIN_MERGE_H__

/*============================================================================
 * Set of subroutines for:
 *  - fusing equivalent vertices,
 *  - managing tolerance reduction,
 *===========================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

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
 * Standard C library headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_intersect.h"
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
 * Creation of new vertices.
 *
 * Update list of equivalent vertices, and assign a vertex (existing or
 * newly created) to each intersection.
 *
 * parameters:
 *   verbosity          <-- verbosity level
 *   edges              <-- list of edges
 *   work               <-> joining mesh maintaining initial vertex data
 *   inter_set          <-> cs_join_inter_set_t structure including
 *                          data on edge-edge  intersections
 *   init_max_vtx_gnum  <-- initial max. global number for vertices
 *   p_n_g_new_vertices <-> pointer to the global number of new vertices
 *   p_vtx_eset         <-> pointer to a structure dealing with vertex
 *                          equivalences
 *---------------------------------------------------------------------------*/

void
cs_join_create_new_vertices(int                     verbosity,
                            const cs_join_edges_t  *edges,
                            cs_join_mesh_t         *work,
                            cs_join_inter_set_t    *inter_set,
                            cs_gnum_t               init_max_vtx_gnum,
                            cs_gnum_t              *p_n_g_new_vertices,
                            cs_join_eset_t        **p_vtx_eset);

/*----------------------------------------------------------------------------
 * Merge of equivalent vertices (and tolerance reduction if necessary)
 *
 * Define a new cs_join_vertex_t structure (stored in "work" structure).
 * Returns an updated cs_join_mesh_t and cs_join_edges_t structures.
 *
 * parameters:
 *   param            <-- set of user-defined parameters for the joining
 *   n_g_vertices_tot <-- global number of vertices (initial parent mesh)
 *   work             <-> pointer to a cs_join_mesh_t structure
 *   vtx_eset         <-- structure storing equivalences between vertices
 *                        (two vertices are equivalent if they are within
 *                        each other's tolerance)
 *---------------------------------------------------------------------------*/

void
cs_join_merge_vertices(cs_join_param_t        param,
                       cs_gnum_t              n_g_vertices_tot,
                       cs_join_mesh_t        *work,
                       const cs_join_eset_t  *vtx_eset);

/*----------------------------------------------------------------------------
 * Merge of equivalent vertices (and reduction of tolerance if necessary)
 *
 * Define a new cs_join_vertex_t structure (stored in "work" structure)
 * Returns an updated cs_join_mesh_t and cs_join_edges_t structures.
 *
 * parameters:
 *   param                <-- set of user-defined parameters for the joining
 *   n_iwm_vertices       <-- initial number of vertices (work mesh struct.)
 *   iwm_vtx_gnum         <-- initial global vertex num. (work mesh struct)
 *   init_max_vtx_gnum    <-- initial max. global number for vertices
 *   rank_face_gnum_index <-- index on face global numbering to determine
 *                            the related rank
 *   p_mesh               <-> pointer to cs_join_mesh_t structure
 *   p_edges              <-> pointer to cs_join_edges_t structure
 *   p_inter_edges        <-> pointer to a cs_join_inter_edges_t struct.
 *   p_local_mesh         <-> pointer to a cs_join_mesh_t structure
 *   p_o2n_vtx_gnum       --> array on blocks on the new global vertex
 *                            numbering for the init. vertices (before inter.)
 *---------------------------------------------------------------------------*/

void
cs_join_merge_update_struct(cs_join_param_t          param,
                            cs_lnum_t                n_iwm_vertices,
                            const cs_gnum_t          iwm_vtx_gnum[],
                            cs_gnum_t                init_max_vtx_gnum,
                            const cs_gnum_t          rank_face_gnum_index[],
                            cs_join_mesh_t         **p_mesh,
                            cs_join_edges_t        **p_edges,
                            cs_join_inter_edges_t  **p_inter_edges,
                            cs_join_mesh_t         **p_local_mesh,
                            cs_gnum_t               *p_o2n_vtx_gnum[]);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_MERGE_H__ */
