/*============================================================================
 * Extract nodal connectivity mesh representations from a native mesh.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_sort.h"

#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_nodal_from_desc.h"
#include "fvm_nodal_order.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Add a subset of a mesh's faces to an external mesh.
 *
 * The lists of faces to extract are optional (if none is given, boundary
 * faces are extracted by default); they do not need to be ordered on input,
 * but they are always ordered on exit (as faces are extracted by increasing
 * number traversal, the lists are reordered to ensure the coherency of
 * the extracted mesh's link to its parent faces, built using these lists).
 *
 * parameters:
 *   mesh             <-- base mesh
 *   extr_mesh        <-- extracted mesh name
 *   boundary_flag    <-- -1 if unspecified, 0 if faces are not on boundary,
 *                        1 if faces are on boundary
 *   include_families <-- include family info if true
 *   i_face_list_size <-- size of i_face_list[] array
 *   b_face_list_size <-- size of b_face_list[] array
 *   i_face_list      <-- list of interior faces (0 to n-1), or NULL
 *   b_face_list      <-- list of boundary faces (0 to n-1), or NULL
 *----------------------------------------------------------------------------*/

static void
_add_faces_to_nodal(const cs_mesh_t  *mesh,
                    fvm_nodal_t      *extr_mesh,
                    int               boundary_flag,
                    bool              include_families,
                    cs_lnum_t         i_face_list_size,
                    cs_lnum_t         b_face_list_size,
                    const cs_lnum_t   i_face_list[],
                    const cs_lnum_t   b_face_list[])
{
  cs_lnum_t   face_id;

  cs_lnum_t  *extr_face_list = NULL;

  cs_lnum_t   face_num_shift[3];
  cs_lnum_t  *face_vertices_idx[2];
  cs_lnum_t  *face_vertices_num[2];
  const int   *_face_families[2];
  const int   **face_families = NULL;

  /* Count the number of faces to convert */

  cs_lnum_t n_max_faces = mesh->n_i_faces + mesh->n_b_faces;
  BFT_MALLOC(extr_face_list, n_max_faces, cs_lnum_t);

  /* Initialize list as marker */

  for (face_id = 0; face_id < n_max_faces; face_id++)
    extr_face_list[face_id] = -1;

  if (b_face_list != NULL) {
    for (face_id = 0; face_id < b_face_list_size; face_id++)
      extr_face_list[b_face_list[face_id]] = 0;
  }
  else {
    for (face_id = 0; face_id < b_face_list_size; face_id++)
      extr_face_list[face_id] = 0;
  }

  if (i_face_list != NULL) {
    for (face_id = 0; face_id < i_face_list_size; face_id++)
      extr_face_list[i_face_list[face_id] + mesh->n_b_faces] = 0;
  }
  else {
    for (face_id = 0; face_id < i_face_list_size; face_id++)
      extr_face_list[face_id + mesh->n_b_faces] = 0;
  }

  /* Convert marked ids to contiguous list (0 to n-1). */

  cs_lnum_t extr_face_count = 0;

  for (face_id = 0; face_id < n_max_faces; face_id++) {
    if (extr_face_list[face_id] == 0) {
      extr_face_list[extr_face_count] = face_id;
      extr_face_count++;
    }
  }

  BFT_REALLOC(extr_face_list, extr_face_count, cs_lnum_t);

  if (include_families) {
    _face_families[0] = mesh->b_face_family;
    _face_families[1] = mesh->i_face_family;
    face_families = _face_families;
  }

  /* Build the nodal connectivity */

  face_num_shift[0] = 0;
  face_num_shift[1] = mesh->n_b_faces + face_num_shift[0];
  face_num_shift[2] = mesh->n_i_faces + face_num_shift[1];

  face_vertices_idx[0] = mesh->b_face_vtx_idx;
  face_vertices_idx[1] = mesh->i_face_vtx_idx;
  face_vertices_num[0] = mesh->b_face_vtx_lst;
  face_vertices_num[1] = mesh->i_face_vtx_lst;

  fvm_nodal_from_desc_add_faces(extr_mesh,
                                boundary_flag,
                                extr_face_count,
                                extr_face_list,
                                2,
                                face_num_shift,
                                (const cs_lnum_t **)face_vertices_idx,
                                (const cs_lnum_t **)face_vertices_num,
                                face_families,
                                NULL);

  BFT_FREE(extr_face_list);
}

/*----------------------------------------------------------------------------
 * Order added faces in an external mesh.
 *
 * parameters:
 *   mesh             <-- base mesh
 *   extr_mesh        <-- extracted mesh name
 *----------------------------------------------------------------------------*/

static void
_order_nodal_faces(const cs_mesh_t  *mesh,
                   fvm_nodal_t      *extr_mesh)
{
  cs_lnum_t   face_id, i;

  cs_lnum_t   n_max_faces = 0;

  cs_gnum_t  *num_glob_fac = NULL;

  /* Count the number of faces to convert */

  n_max_faces = mesh->n_i_faces + mesh->n_b_faces;

  /* In case of parallelism or face renumbering, sort faces by
     increasing global number */

  if (mesh->global_i_face_num != NULL || mesh->global_b_face_num != NULL) {

    BFT_MALLOC(num_glob_fac, n_max_faces, cs_gnum_t);

    if (mesh->global_b_face_num == NULL) {
      for (face_id = 0; face_id < mesh->n_b_faces; face_id++)
        num_glob_fac[face_id] = face_id + 1;
    }
    else {
      for (face_id = 0; face_id < mesh->n_b_faces; face_id++)
        num_glob_fac[face_id] = mesh->global_b_face_num[face_id];
    }

    assert(mesh->n_g_b_faces + mesh->n_g_i_faces > 0);

    if (mesh->global_i_face_num == NULL) {
      for (face_id = 0, i = mesh->n_b_faces;
           face_id < mesh->n_i_faces;
           face_id++, i++)
        num_glob_fac[i] = face_id + 1 + mesh->n_g_b_faces;
    }
    else {
      for (face_id = 0, i = mesh->n_b_faces;
           face_id < mesh->n_i_faces;
           face_id++, i++)
        num_glob_fac[i] = mesh->global_i_face_num[face_id] + mesh->n_g_b_faces;
    }

  }

  /* Sort faces by increasing global number */

  fvm_nodal_order_faces(extr_mesh, num_glob_fac);
  fvm_nodal_init_io_num(extr_mesh, num_glob_fac, 2);

  if (num_glob_fac != NULL)
    BFT_FREE(num_glob_fac);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
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
                               cs_lnum_t        **const p_cell_faces_val)
{
  cs_lnum_t  cell_id, c_id1, c_id2, face_id, n_loc_cells;

  cs_lnum_t  *cell_face_count = NULL;
  cs_lnum_t  *cell_faces_idx = NULL;
  cs_lnum_t  *cell_faces_val = NULL;

  const cs_lnum_t n_b_faces = CS_MAX(mesh->n_b_faces_all, mesh->n_b_faces);

  /* Allocate and initialize cell ->faces index */

  n_loc_cells = mesh->n_cells;
  if (extr_cell_id != NULL)
    n_loc_cells = extr_cell_size;

  BFT_MALLOC(cell_faces_idx, n_loc_cells + 1, cs_lnum_t);

  for (cell_id = 0; cell_id < n_loc_cells + 1; cell_id++)
    cell_faces_idx[cell_id] = 0;

  /* Count number of faces per cell (we assign the temporary counter
     corresponding to cell_id to cell_faces_idx[cell_id + 1] rather than
     cell_faces_idx[cell_id] to simplify the next step) */

  /* Remark: test if cell_id < mesh->n_cells on internal faces so
     as to ignore ghost cells */

  for (face_id = 0; face_id < n_b_faces; face_id++) {
    cell_id = mesh->b_face_cells[face_id];
    if (extr_cell_id != NULL)
      cell_id = extr_cell_id[cell_id];
    if (cell_id > -1)
      cell_faces_idx[cell_id + 1] += 1;
  }

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    c_id1 = mesh->i_face_cells[face_id][0];
    c_id2 = mesh->i_face_cells[face_id][1];
    if (extr_cell_id != NULL) {
      if (c_id1 < mesh->n_cells)
        c_id1 = extr_cell_id[c_id1];
      else
        c_id1 = -1;
      if (c_id2 < mesh->n_cells)
        c_id2 = extr_cell_id[c_id2];
      else
        c_id2 = -1;
    }
    if (c_id1 > -1 && c_id1 < mesh->n_cells)
      cell_faces_idx[c_id1 + 1] += 1;
    if (c_id2 > -1 && c_id2 < mesh->n_cells)
      cell_faces_idx[c_id2 + 1] += 1;
  }

  /* Build cell -> faces index */

  cell_faces_idx[0] = 1;
  for (cell_id = 0; cell_id < n_loc_cells; cell_id++)
    cell_faces_idx[cell_id + 1] += cell_faces_idx[cell_id];

  /* Build array of values */

  BFT_MALLOC(cell_faces_val, cell_faces_idx[n_loc_cells] - 1, cs_lnum_t);
  BFT_MALLOC(cell_face_count, n_loc_cells, cs_lnum_t);

  for (cell_id = 0; cell_id < n_loc_cells; cell_id++)
    cell_face_count[cell_id] = 0;

  for (face_id = 0; face_id < n_b_faces; face_id++) {
    cell_id = mesh->b_face_cells[face_id];
    if (extr_cell_id != NULL)
      cell_id = extr_cell_id[cell_id];
    if (cell_id > -1) {
      cell_faces_val[cell_faces_idx[cell_id] + cell_face_count[cell_id] - 1]
        = face_id + 1;
      cell_face_count[cell_id] += 1;
    }
  }

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    c_id1 = mesh->i_face_cells[face_id][0];
    c_id2 = mesh->i_face_cells[face_id][1];
    if (extr_cell_id != NULL) {
      if (c_id1 < mesh->n_cells)
        c_id1 = extr_cell_id[c_id1];
      else
        c_id1 = -1;
      if (c_id2 < mesh->n_cells)
        c_id2 = extr_cell_id[c_id2];
      else
        c_id2 = -1;
    }
    if (c_id1 > -1 && c_id1 < mesh->n_cells) {
      cell_faces_val[cell_faces_idx[c_id1] + cell_face_count[c_id1] - 1]
        =   face_id + n_b_faces + 1;
      cell_face_count[c_id1] += 1;
    }
    if (c_id2 > -1 && c_id2 < mesh->n_cells) {
      cell_faces_val[cell_faces_idx[c_id2] + cell_face_count[c_id2] - 1]
        = -(face_id + n_b_faces + 1);
      cell_face_count[c_id2] += 1;
    }
  }

  BFT_FREE(cell_face_count);

  /* Return values */

  *p_cell_faces_idx = cell_faces_idx;
  *p_cell_faces_val = cell_faces_val;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
 {
   cs_lnum_t ipos, ival;
   /* Print arrays */
   bft_printf("dbg : cs_mesh_ret_cel_fac\n"
              "nombre de cellules extraites = %d\n", extr_cell_size);
   for (ipos = 0; ipos < extr_cell_size; ipos++) {
     bft_printf("  cellule %d\n", ipos);
     bft_printf("    cell_faces_idx[%d] = %d\n", ipos, cell_faces_idx[ipos]);
     for (ival = cell_faces_idx[ipos]     - 1;
          ival < cell_faces_idx[ipos + 1] - 1;
          ival++)
       bft_printf("      cell_faces_val[%d] = %d\n",
                  ival, cell_faces_val[ival]);
   }
   bft_printf("  cell_faces_idx[%d] = %d\n", ipos, cell_faces_idx[ipos]);
 }
#endif

}

/*----------------------------------------------------------------------------
 * Build a nodal connectivity structure from a subset of a mesh's cells.
 *
 * The list of cells to extract is optional (if none is given, all cells
 * faces are extracted by default).
 *
 * parameters:
 *   mesh             <-- base mesh
 *   name             <-- extracted mesh name
 *   include_families <-- include family info if true
 *   cell_list_size   <-- size of cell_list[] array
 *   cell_list        <-- list of cells (0 to n-1), or NULL
 *
 * returns:
 *   pointer to extracted nodal mesh
 *----------------------------------------------------------------------------*/

fvm_nodal_t  *
cs_mesh_connect_cells_to_nodal(const cs_mesh_t  *mesh,
                               const char       *name,
                               bool              include_families,
                               cs_lnum_t         cell_list_size,
                               const cs_lnum_t   cell_list[])
{
  cs_lnum_t   face_id, cell_id;

  int         null_family = 0;
  cs_lnum_t   extr_cell_count = 0, i_face_count = 0, b_face_count = 0;
  cs_lnum_t  *extr_cell_idx = NULL;
  cs_lnum_t  *extr_cell_ids = NULL;

  cs_lnum_t  *cell_face_idx = NULL;
  cs_lnum_t  *cell_face_num = NULL;

  cs_lnum_t  *i_face_list= NULL, *b_face_list = NULL;

  cs_lnum_t  face_num_shift[3];
  cs_lnum_t  *face_vertices_idx[2];
  cs_lnum_t  *face_vertices_num[2];
  cs_lnum_t  *polyhedra_faces = NULL;
  const int  *cell_family = NULL;

  fvm_nodal_t  *extr_mesh;

  const cs_lnum_t n_b_faces = CS_MAX(mesh->n_b_faces_all, mesh->n_b_faces);

  /* If a family has no attributes, it must be 1st by construction
     (as families are sorted when merging duplicates) */

  if (mesh->n_families > 0) {
    if (mesh->family_item[0] == 0)
      null_family = 1;
  }

  /* Check that the mesh contains face -> vertices connectivity */

  if (mesh->b_face_vtx_idx == NULL || mesh->i_face_vtx_idx == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The main mesh does not contain any face -> vertices\n"
                "connectivity, necessary for the nodal connectivity\n"
                "reconstruction (cs_mesh_connect_cells_to_nodal)."));

  if (include_families) {
    BFT_MALLOC(i_face_list, mesh->n_i_faces, cs_lnum_t);
    BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);
  }

  /* Count the number of cells to convert */

  if (cell_list != NULL) {

    BFT_MALLOC(extr_cell_ids, cell_list_size, cs_lnum_t);
    BFT_MALLOC(extr_cell_idx, mesh->n_cells_with_ghosts, cs_lnum_t);

    /* Initialize index as marker */

    for (cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
      extr_cell_idx[cell_id] = -1;
    for (cell_id = 0; cell_id < cell_list_size; cell_id++) {
      if (cell_list[cell_id] <= mesh->n_cells)
        extr_cell_idx[cell_list[cell_id]] = 1;
    }

    /* Also mark faces bearing families */

    if (include_families) {

      for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
        cs_lnum_t c_id_0 = mesh->i_face_cells[face_id][0];
        cs_lnum_t c_id_1 = mesh->i_face_cells[face_id][1];
        if (   (extr_cell_idx[c_id_0] == 1 || extr_cell_idx[c_id_1] == 1)
            && (mesh->i_face_family[face_id] != null_family)) {
          i_face_list[i_face_count++] = face_id;
        }
      }
      BFT_REALLOC(i_face_list, i_face_count, cs_lnum_t);

      for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
        cs_lnum_t c_id = mesh->b_face_cells[face_id];
        if (   (extr_cell_idx[c_id] == 1)
            && (mesh->b_face_family[face_id] != null_family)) {
          b_face_list[b_face_count++] = face_id;
        }
      }
      BFT_REALLOC(b_face_list, b_face_count, cs_lnum_t);

    }

    /* Convert marked ids to indexes (1 to n) and reconstruct values of
       cell_list[] to ensure that it is ordered. */

    extr_cell_count = 0;
    for (cell_id = 0; cell_id < mesh->n_cells; cell_id++) {
      if (extr_cell_idx[cell_id] == 1) {
        extr_cell_ids[extr_cell_count] = cell_id;
        extr_cell_idx[cell_id] = extr_cell_count++;
      }
    }

    assert(extr_cell_count <= cell_list_size);

  }
  else {
    extr_cell_count = CS_MIN(mesh->n_cells, cell_list_size);
    extr_cell_idx = NULL;

    if (include_families && extr_cell_count > 0) {

      for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
        cs_lnum_t c_id_0 = mesh->i_face_cells[face_id][0];
        cs_lnum_t c_id_1 = mesh->i_face_cells[face_id][1];
        if (   (c_id_0 < extr_cell_count || c_id_1 < extr_cell_count)
            && (mesh->i_face_family[face_id] != null_family)) {
          i_face_list[i_face_count++] = face_id;
        }
      }
      BFT_REALLOC(i_face_list, i_face_count, cs_lnum_t);

      for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
        cs_lnum_t c_id = mesh->b_face_cells[face_id];
        if (   (c_id < extr_cell_count)
            && (mesh->b_face_family[face_id] != null_family)) {
          b_face_list[b_face_count++] = face_id;
        }
      }
      BFT_REALLOC(b_face_list, b_face_count, cs_lnum_t);

    }
  }

  /* Extract "cells -> faces" connectivity */

  cs_mesh_connect_get_cell_faces(mesh,
                                 extr_cell_count,
                                 extr_cell_idx,
                                 &cell_face_idx,
                                 &cell_face_num);

  if (extr_cell_idx != NULL)
    BFT_FREE(extr_cell_idx);

  /* Build nodal connectivity */

  face_num_shift[0] = 0;
  face_num_shift[1] = n_b_faces + face_num_shift[0];
  face_num_shift[2] = mesh->n_i_faces + face_num_shift[1];

  face_vertices_idx[0] = mesh->b_face_vtx_idx;
  face_vertices_idx[1] = mesh->i_face_vtx_idx;
  face_vertices_num[0] = mesh->b_face_vtx_lst;
  face_vertices_num[1] = mesh->i_face_vtx_lst;

  extr_mesh = fvm_nodal_create(name, 3);

  fvm_nodal_set_parent(extr_mesh, mesh);

  if (include_families)
    cell_family = mesh->cell_family;

  fvm_nodal_from_desc_add_cells(extr_mesh,
                                extr_cell_count,
                                2,
                                face_num_shift,
                                (const cs_lnum_t **)face_vertices_idx,
                                (const cs_lnum_t **)face_vertices_num,
                                cell_face_idx,
                                cell_face_num,
                                cell_family,
                                extr_cell_ids,
                                &polyhedra_faces);

  BFT_FREE(extr_cell_ids);

  /* Also add faces bearing families */

  if (include_families) {

    _add_faces_to_nodal(mesh,
                        extr_mesh,
                        1,
                        include_families,
                        0,
                        b_face_count,
                        NULL,
                        b_face_list);

    _add_faces_to_nodal(mesh,
                        extr_mesh,
                        0,
                        include_families,
                        i_face_count,
                        0,
                        i_face_list,
                        NULL);

    _order_nodal_faces(mesh, extr_mesh);

    BFT_FREE(i_face_list);
    BFT_FREE(b_face_list);
  }

  fvm_nodal_set_shared_vertices(extr_mesh, mesh->vtx_coord);
  fvm_nodal_set_group_class_set(extr_mesh, mesh->class_defs);

  /* Free memory */

  BFT_FREE(polyhedra_faces);

  BFT_FREE(cell_face_idx);
  BFT_FREE(cell_face_num);

  /* Sort cells by increasing global number */

  fvm_nodal_order_cells(extr_mesh, mesh->global_cell_num);
  fvm_nodal_init_io_num(extr_mesh, mesh->global_cell_num, 3);

  /* Sort vertices by increasing global number */

  fvm_nodal_order_vertices(extr_mesh, mesh->global_vtx_num);
  fvm_nodal_init_io_num(extr_mesh, mesh->global_vtx_num, 0);

  /* We are done */

  return extr_mesh;
}

/*----------------------------------------------------------------------------
 * Build a nodal connectivity structure from a subset of a mesh's faces.
 *
 * The lists of faces to extract are optional (if none is given, boundary
 * faces are extracted by default).
 *
 * parameters:
 *   mesh             <-- base mesh
 *   name             <-- extracted mesh name
 *   include_families <-- include family info if true
 *   i_face_list_size <-- size of i_face_list[] array
 *   b_face_list_size <-- size of b_face_list[] array
 *   i_face_list      <-- list of interior faces (0 to n-1), or NULL
 *   b_face_list      <-- list of boundary faces (0 to n-1), or NULL
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
                               const cs_lnum_t   i_face_list[],
                               const cs_lnum_t   b_face_list[])
{
  fvm_nodal_t  *extr_mesh = NULL;

  /* Check that the mesh contains face -> vertices connectivity */

  if (mesh->b_face_vtx_idx == NULL || mesh->i_face_vtx_idx == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("The main mesh does not contain any face -> vertices\n"
                "connectivity, necessary for the nodal connectivity\n"
                "reconstruction (cs_mesh_connect_faces_to_nodal)."));

  extr_mesh = fvm_nodal_create(name, 3);

  fvm_nodal_set_parent(extr_mesh, mesh);

  _add_faces_to_nodal(mesh,
                      extr_mesh,
                      -1,
                      include_families,
                      i_face_list_size,
                      b_face_list_size,
                      i_face_list,
                      b_face_list);

  _order_nodal_faces(mesh, extr_mesh);

  fvm_nodal_set_shared_vertices(extr_mesh,
                                mesh->vtx_coord);

  /* Sort vertices by increasing global number */

  fvm_nodal_order_vertices(extr_mesh, mesh->global_vtx_num);
  fvm_nodal_init_io_num(extr_mesh, mesh->global_vtx_num, 0);

  /* Define group classes */

  if (include_families)
    fvm_nodal_set_group_class_set(extr_mesh, mesh->class_defs);

  /* We are done */

  return extr_mesh;
}

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
                                  cs_lnum_t   **v2c)
{
  const cs_lnum_t n_vertices = mesh->n_vertices;
  const cs_lnum_t n_b_faces = CS_MAX(mesh->n_b_faces_all, mesh->n_b_faces);

  /* Mark vertices which may be split (vertices lying on new boundary faces) */

  cs_lnum_t  *_v2c_idx;
  BFT_MALLOC(_v2c_idx, n_vertices+1, cs_lnum_t);

  _v2c_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    _v2c_idx[i+1] = 0;

  /* Now build vertex -> cells index
     (which will contain duplicate entries at first) */

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->i_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0) {
        if (mesh->i_face_cells[f_id][0] > -1)
          _v2c_idx[vtx_id + 1] += 1;
        if (mesh->i_face_cells[f_id][1] > -1)
          _v2c_idx[vtx_id + 1] += 1;
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t s_id = mesh->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->b_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->b_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0)
        _v2c_idx[vtx_id + 1] += 1;
    }
  }

  /* Transform count to index */

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    _v2c_idx[i+1] += _v2c_idx[i];

  /* Now define selected vertex->cell adjacency */

  cs_lnum_t *_v2c;
  BFT_MALLOC(_v2c, _v2c_idx[n_vertices], cs_lnum_t);

  cs_lnum_t *v2c_count;
  BFT_MALLOC(v2c_count, n_vertices, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v2c_count[i] = 0;

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->i_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0) {
        cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
        cs_lnum_t j = _v2c_idx[vtx_id] + v2c_count[vtx_id];
        if (c_id_0 > -1) {
          _v2c[j++] = c_id_0;
          v2c_count[vtx_id] += 1;
        }
        if (c_id_1 > -1) {
          _v2c[j++] = c_id_1;
          v2c_count[vtx_id] += 1;
        }
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    cs_lnum_t s_id = mesh->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->b_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->b_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0) {
        cs_lnum_t c_id_0 = mesh->b_face_cells[f_id];
        cs_lnum_t j = _v2c_idx[vtx_id] + v2c_count[vtx_id];
        _v2c[j] = c_id_0;
        v2c_count[vtx_id] += 1;
      }
    }
  }

  BFT_FREE(v2c_count);

  /* Order and compact adjacency array */

  cs_sort_indexed(n_vertices, _v2c_idx, _v2c);

  cs_lnum_t *tmp_v2c_idx = NULL;

  BFT_MALLOC(tmp_v2c_idx, n_vertices+1, cs_lnum_t);
  memcpy(tmp_v2c_idx, _v2c_idx, (n_vertices+1)*sizeof(cs_lnum_t));

  cs_lnum_t k = 0;

  for (cs_lnum_t i = 0; i < n_vertices; i++) {
    cs_lnum_t js = tmp_v2c_idx[i];
    cs_lnum_t je = tmp_v2c_idx[i+1];
    cs_lnum_t v2c_prev = -1;
    _v2c_idx[i] = k;
    for (cs_lnum_t j = js; j < je; j++) {
      if (v2c_prev != _v2c[j]) {
        _v2c[k++] = _v2c[j];
        v2c_prev = _v2c[j];
      }
    }
  }
  _v2c_idx[n_vertices] = k;

  assert(_v2c_idx[n_vertices] <= tmp_v2c_idx[n_vertices]);

  BFT_FREE(tmp_v2c_idx);
  BFT_REALLOC(_v2c, _v2c_idx[n_vertices], cs_lnum_t);

  *v2c_idx = _v2c_idx;
  *v2c = _v2c;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
