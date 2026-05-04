#ifndef __CS_MESH_ALGORITHM_H__
#define __CS_MESH_ALGORITHM_H__

/*============================================================================
 * Mesh refinement.
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_new      <-- new number of elements
 *   n2o        <-- new to old array (same as old element ids list)
 *   global_num <-> global numbering (allocated if initially nullptr)
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_algorithm_n2o_update_global_num(cs_lnum_t          n_new,
                                        const cs_lnum_t    n2o[],
                                        cs_gnum_t        **global_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Merge cells based on renumbering array.
 *
 * Interior faces separating merged cells are removed.
 *
 * \param[in, out]  m              mesh
 * \param[in]       n_new          new number of cells
 * \param[in]       c_o2n          cell old to new renumbering
 * \param[out]      i_f_n2o_pre    new-to-old map of interior faces induced by
 *                                 cell merge
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_algorithm_merge_cells(cs_mesh_t       *m,
                              cs_lnum_t        n_new,
                              const cs_lnum_t  c_o2n[],
                              cs_lnum_t       *i_f_n2o[]);

/*----------------------------------------------------------------------------
 *
 * \brief Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_old      <-- old number of elements
 *   n_g_old    <-- old global number of elements
 *   o2n_idx    <-- old to new index
 *   global_num <-> global numbering (allocated if initially nullptr)
 *
 * returns:
 *   new global number of elements
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_algorithm_o2n_idx_update_global_num(cs_lnum_t          n_old,
                                            cs_gnum_t          n_g_old,
                                            const cs_lnum_t    o2n_idx[],
                                            cs_gnum_t        **global_num);

/*----------------------------------------------------------------------------*/
/*
 * \brief Build global numbers for new vertices on edges, faces, or cells.
 *
 * These vertices are appended at the end of the initial vertex definitions.
 * The numbering arrays should be resized before calling this function
 * (to allow for vertices inserted on edges, faces, and possibly
 * cells with a single resize).
 *
 * Each call of this function updates the global number of vertices
 * member of the mesh structure.
 *
 * \param[in, out]  m          mesh
 * \param[in]       n_elts     number of parent elements
 * \param[in]       n_g_elts   global number of parent elements
 * \param[in]       elt_v_idx  for each element, start index of added vertices
 * \param[in]       g_elt_num  global number of each element
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_algorithm_build_add_vertices_gnum(cs_mesh_t       *m,
                                          cs_lnum_t        n_elts,
                                          cs_gnum_t        n_g_elts,
                                          const cs_lnum_t  elt_v_idx[],
                                          const cs_gnum_t  g_elt_num[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Sync edges flag for parallelism and determine associated
 *        added vertices global numbers.
 *
 * \param[in]       m            pointer to mesh structure
 * \param[in]       v2v          vertex adjacency
 * \param[in, out]  e_v_flag     for each edge, flag (count) for added vertices
 * \param[out]      g_edges_num  global edges number, or nullptr
 *
 * \return: global number of edges
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_mesh_algorithm_sync_edges_flag(const cs_mesh_t        *m,
                                  const cs_adjacency_t   *v2v,
                                  cs_lnum_t               e_v_flag[],
                                  cs_gnum_t              *g_edges_num);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the face -> edges connectivity which is stored in a
 *        cs_adjacency_t structure
 *
 * \param[in] m    pointer to a cs_mesh_t structure
 * \param[in] v2v  pointer to the cs_adjacency_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adjacency_t *
cs_mesh_algorithm_build_f2e_connect(const cs_mesh_t      *m,
                                    const cs_adjacency_t *v2v);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_ALGORITHM_H__ */
