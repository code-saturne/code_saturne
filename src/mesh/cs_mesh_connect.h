#ifndef __CS_MESH_CONNECT_H__
#define __CS_MESH_CONNECT_H__

/*============================================================================
 * Extract nodal connectivity mesh representations from a native mesh.
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

#include "fvm_nodal.h"

#include "cs_base.h"
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
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extract a mesh's "cells -> faces" connectivity.
 *
 * We consider a common numbering for internal and boundary faces, in which
 * boundary faces are defined first. The common id for the i-th boundary
 * face is thus i, and that of the j-th interior face is n_b_faces + j.
 *
 * If ind_cel_extr != NULL, then:
 * --- ind_cel_extr[cell_id] = id in the list to extract (0 to n-1)
 *     if cell cell_id should be extracted
 * --- ind_cel_extr[cell_id] = -1 if cells cell_id should be ignored
 *
 * parameters:
 *   mesh             <-- pointer to mesh structure
 *   extr_cell_size   <-- size of extr_cell_id[] array
 *   extr_cell_id     <-- extr_cell_id = ids of extracted cells, or -1
 *   p_cell_faces_idx --> cells -> faces index
 *   p_cell_faces_val --> cells -> faces connectivity
 *----------------------------------------------------------------------------*/

void
cs_mesh_connect_get_cell_faces(const cs_mesh_t         *mesh,
                               cs_lnum_t                extr_cell_size,
                               const cs_lnum_t          extr_cell_id[],
                               cs_lnum_t        **const p_cell_faces_idx,
                               cs_lnum_t        **const p_cell_faces_val);

/*----------------------------------------------------------------------------
 * Build a nodal connectivity structure from a subset of a mesh's cells.
 *
 * The list of cells to extract is optional (if none is given, all cells
 * faces are extracted by default); it does not need to be ordered on input,
 * but is always ordered on exit (as cells are extracted by increasing number
 * traversal, the list is reordered to ensure the coherency of the extracted
 * mesh's link to its parent cells, built using this list).
 *
 * parameters:
 *   mesh             <-- base mesh
 *   name             <-- extracted mesh name
 *   include_families <-- include family info if true
 *   cell_list_size   <-- size of cell_list[] array
 *   cell_list        <-> list of cells (1 to n), or NULL
 *
 * returns:
 *   pointer to extracted nodal mesh
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
cs_mesh_connect_cells_to_nodal(const cs_mesh_t  *mesh,
                               const char       *name,
                               bool              include_families,
                               cs_lnum_t         cell_list_size,
                               cs_lnum_t         cell_list[]);

/*----------------------------------------------------------------------------
 * Build a nodal connectivity structure from a subset of a mesh's faces.
 *
 * The lists of faces to extract are optional (if none is given, boundary
 * faces are extracted by default); they do not need to be ordered on input,
 * but they are always ordered on exit (as faces are extracted by increasing
 * number traversal, the lists are reordered to ensure the coherency of
 * the extracted mesh's link to its parent faces, built using these lists).
 *
 * parameters:
 *   mesh             <-- base mesh
 *   name             <-- extracted mesh name
 *   include_families <-- include family info if true
 *   i_face_list_size <-- size of i_face_list[] array
 *   b_face_list_size <-- size of b_face_list[] array
 *   i_face_list      <-> list of interior faces (1 to n), or NULL
 *   b_face_list      <-> list of boundary faces (1 to n), or NULL
 *
 * returns:
 *   pointer to extracted nodal mesh
 *----------------------------------------------------------------------------*/

fvm_nodal_t *
cs_mesh_connect_faces_to_nodal(const cs_mesh_t  *mesh,
                               const char       *name,
                               bool              include_families,
                               cs_lnum_t         i_face_list_size,
                               cs_lnum_t         b_face_list_size,
                               cs_lnum_t         i_face_list[],
                               cs_lnum_t         b_face_list[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a vertex to cell connectivity for marked vertices only.
 *
 * It is the caller's responsibility to free the v2c_idx and v2c arrays,
 * which are allocated by this function.
 *
 * \param[in]    mesh      pointer to mesh structure
 * \param[in]    v_flag    vertex selection flag (0: not selected, 1: selected)
 * \param[out]   v2c_idx   vertex to cells index (size: mesh->n_vertices +1)
 * \param[out]   v2c       vertex to cells
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_connect_vertices_to_cells(cs_mesh_t    *mesh,
                                  const char    v_flag[],
                                  cs_lnum_t   **v2c_idx,
                                  cs_lnum_t   **v2c);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_CONNECT_H__ */
