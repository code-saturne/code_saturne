#ifndef __FVM_NODAL_FROM_DESC_H__
#define __FVM_NODAL_FROM_DESC_H__

/*============================================================================
 * Initialization of a nodal connectivity definition based upon
 * a (possibly partial) descending connectivity
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Convert and add cells from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_cells[] argument is non-NULL, cells
 * {extr_cells[0], extr_cells[1], extr_cells[n_extr_cells - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, cells
 * {1, 2, ..., n_extr_cells} are considered.
 *
 * In addition, an optional parent_cell_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_cell_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_cells    <-- count of cells to add
 *   extr_cells      <-- optional filter list of cells to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   cell_gc_id      <-- cell -> group class ids, or NULL
 *   parent_cell_num <-- cell -> parent cell number (1 to n) if non-trivial
 *                       (i.e. if cell definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *   cell_face_list  --> numbers of faces defining polyhedra
 *----------------------------------------------------------------------------*/

void
fvm_nodal_from_desc_add_cells(fvm_nodal_t        *this_nodal,
                              const cs_lnum_t     n_extr_cells,
                              const cs_lnum_t     extr_cells[],
                              const int           n_face_lists,
                              const cs_lnum_t     face_list_shift[],
                              const cs_lnum_t    *face_vertex_idx[],
                              const cs_lnum_t    *face_vertex_num[],
                              const cs_lnum_t     cell_face_idx[],
                              const cs_lnum_t     cell_face_num[],
                              const int           cell_gc_id[],
                              const cs_lnum_t     parent_cell_num[],
                              cs_lnum_t          *cell_face_list[]);

/*----------------------------------------------------------------------------
 * Convert and add faces from an descending connectivity mesh to a nodal mesh.
 *
 * If the optional filter list extr_faces[] argument is non-NULL, faces
 * {extr_faces[0], extr_faces[1], extr_faces[n_extr_faces - 1]} are converted
 * and added to the nodal mesh. If this filter is set to NULL, faces
 * {1, 2, ..., n_extr_faces} are considered.
 *
 * In addition, an optional parent_face_num[] array may also be given, in
 * case the descending connectivity mesh definition is based on a temporary
 * subset of a parent mesh, (corresponding to the parent_face_num[] list,
 * using 1 to n numbering), and the final nodal mesh element parent numbering
 * should correspond to that parent mesh and not the temporary subset.
 *
 * parameters:
 *   this_nodal      <-> nodal mesh structure
 *   n_extr_faces    <-- count of faces to add
 *   extr_faces      <-- optional filter list of faces to extract (1 to n)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   face_gc_id      <-- face -> group class ids, or NULL (per face list)
 *   parent_face_num <-- face -> parent face number (1 to n) if non-trivial
 *                       (i.e. if face definitions correspond to a subset
 *                       of the parent mesh), NULL otherwise.
 *----------------------------------------------------------------------------*/

void
fvm_nodal_from_desc_add_faces(fvm_nodal_t        *this_nodal,
                              const cs_lnum_t     n_extr_faces,
                              const cs_lnum_t     extr_faces[],
                              const int           n_face_lists,
                              const cs_lnum_t     face_list_shift[],
                              const cs_lnum_t    *face_vertex_idx[],
                              const cs_lnum_t    *face_vertex_num[],
                              const int          *face_gc_id[],
                              const cs_lnum_t     parent_face_num[]);

/*----------------------------------------------------------------------------
 * Determination of a given cell's type.
 *
 * If the optional cell_vtx[8] array is given, it is filled with the vertex
 * indexes of cell's vertices, unless the cell is a general polyhedron.
 *
 * parameters:
 *   cell_id         <-- cell id (0 to n-1)
 *   n_face_lists    <-- number of face lists
 *   face_list_shift <-- face list to common number index shifts;
 *                       size: n_face_lists
 *   face_vertex_idx <-- face -> vertex indexes (per face list)
 *   face_vertex     <-- face -> vertex ids (per face list)
 *   cell_face_idx   <-- cell -> face indexes (1 to n)
 *   cell_face_num   <-- cell -> face numbers (1 to n)
 *   vertex_num      --> nodal connectivity of cell, if not a general
 *                       polyhedron
 *
 * returns:
 *   type of cell defined by cell_id
 *----------------------------------------------------------------------------*/

fvm_element_t
fvm_nodal_from_desc_cell(const cs_lnum_t    cell_id,
                         const int          n_face_lists,
                         const cs_lnum_t    face_list_shift[],
                         const cs_lnum_t   *face_vertex_idx[],
                         const cs_lnum_t   *face_vertex_num[],
                         const cs_lnum_t    cell_face_idx[],
                         const cs_lnum_t    cell_face_num[],
                         cs_lnum_t          vertex_num[8]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NODAL_FROM_DESC_H__ */
