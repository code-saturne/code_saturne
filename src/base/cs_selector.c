/*============================================================================
 * Selection criteria for cells, boundary and interior faces.
 *============================================================================*/

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "bft_mem_usage.h"
#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

#include "cs_selector.h"

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fill a list of boundary face numbers verifying a given selection criteria.
 *
 * parameters:
 *   criteria        <-- selection criteria string
 *   n_b_faces       --> number of selected interior faces
 *   b_face_num_list --> list of selected boundary face numbers
 *                       (1 to n, preallocated to cs_glob_mesh->n_b_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_num_list(const char  *criteria,
                                cs_lnum_t   *n_b_faces,
                                cs_lnum_t    b_face_num_list[])
{
  int c_id;

  *n_b_faces = 0;

  if (cs_glob_mesh->select_b_faces == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%sd: %s is not defined at this stage."),
                __func__, "cs_glob_mesh->select_b_faces");

  c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                               criteria,
                               1,
                               n_b_faces,
                               b_face_num_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group \"%s\" in the selection criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any boundary face.\n"),
               missing, criteria);
  }
}

/*----------------------------------------------------------------------------
 * Fill a list of interior faces verifying a given selection criteria.
 *
 * parameters:
 *   criteria        <-- selection criteria string
 *   n_i_faces       --> number of selected interior faces
 *   i_face_num_list --> list of selected interior face numbers
 *                       (1 to n, preallocated to cs_glob_mesh->n_i_faces)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_i_face_num_list(const char  *criteria,
                                cs_lnum_t   *n_i_faces,
                                cs_lnum_t    i_face_num_list[])
{
  int c_id;

  *n_i_faces = 0;

  if (cs_glob_mesh->select_i_faces == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%sd: %s is not defined at this stage."),
                __func__, "cs_glob_mesh->select_i_faces");

  c_id = fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                               criteria,
                               1,
                               n_i_faces,
                               i_face_num_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_i_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_i_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group \"%s\" in the selection criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any interior face.\n"),
               missing, criteria);
  }
}

/*----------------------------------------------------------------------------
 * Fill a list of cells verifying a given selection criteria.
 *
 * parameters:
 *   criteria      <-- selection criteria string
 *   n_cells       --> number of selected cells
 *   cell_num_list --> list of selected cell numbers
 *                     (1 to n, preallocated to cs_glob_mesh->n_cells)
 *----------------------------------------------------------------------------*/

void
cs_selector_get_cell_num_list(const char  *criteria,
                              cs_lnum_t   *n_cells,
                              cs_lnum_t    cell_num_list[])
{
  int c_id;

  *n_cells = 0;

  if (cs_glob_mesh->select_b_faces == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%sd: %s is not defined at this stage."),
                __func__, "cs_glob_mesh->select_b_faces");

  c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                               criteria,
                               1,
                               n_cells,
                               cell_num_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group \"%s\" in the selection criteria:\n"
                 "\"%s\"\n"
                 " does not correspond to any cell.\n"),
               missing, criteria);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of boundary faces verifying a given selection criteria.
 *
 * \param[in]   criteria     selection criteria string
 * \param[out]  n_b_faces    number of selected boundary faces
 * \param[out]  b_face_list  list of selected boundary faces
 *                           (0 to n-1, preallocated to cs_glob_mesh->n_b_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_list(const char  *criteria,
                            cs_lnum_t   *n_b_faces,
                            cs_lnum_t    b_face_list[])
{
  int c_id;

  *n_b_faces = 0;

  if (cs_glob_mesh->select_b_faces != NULL) {

    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 criteria,
                                 0,
                                 n_b_faces,
                                 b_face_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The group \"%s\" in the selection criteria:\n"
                   "\"%s\"\n"
                   " does not correspond to any boundary face.\n"),
                 missing, criteria);

    }

  }

  else {

    cs_mesh_t *mesh = cs_glob_mesh;

    bool del_class_defs = (mesh->class_defs == NULL) ? true : false;

    cs_mesh_init_group_classes(mesh);

    cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;

    cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

    fvm_selector_t *sel_b_faces = fvm_selector_create(mesh->dim,
                                                      mesh->n_b_faces,
                                                      mesh->class_defs,
                                                      mesh->b_face_family,
                                                      1,
                                                      b_face_cog,
                                                      b_face_normal);

    c_id = fvm_selector_get_list(sel_b_faces,
                                 criteria,
                                 0,
                                 n_b_faces,
                                 b_face_list);

    BFT_FREE(b_face_cog);
    BFT_FREE(b_face_normal);

    if (del_class_defs)
      mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

    sel_b_faces = fvm_selector_destroy(sel_b_faces);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of interior faces verifying a given selection criteria.
 *
 * \param[in]   criteria     selection criteria string
 * \param[out]  n_i_faces    number of selected interior faces
 * \param[out]  i_face_list  list of selected interior faces
 *                           (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_i_face_list(const char  *criteria,
                            cs_lnum_t   *n_i_faces,
                            cs_lnum_t    i_face_list[])
{
  int c_id;

  *n_i_faces = 0;

  if (cs_glob_mesh->select_i_faces != NULL) {

    c_id = fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                                 criteria,
                                 0,
                                 n_i_faces,
                                 i_face_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_i_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_i_faces, c_id, 0);
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The group \"%s\" in the selection criteria:\n"
                   "\"%s\"\n"
                   " does not correspond to any interior face.\n"),
                 missing, criteria);
    }

  }

  else {

    cs_mesh_t *mesh = cs_glob_mesh;

    bool del_class_defs = (mesh->class_defs == NULL) ? true : false;

    cs_mesh_init_group_classes(mesh);

    cs_real_t  *i_face_cog = NULL, *i_face_normal = NULL;

    cs_mesh_quantities_i_faces(mesh, &i_face_cog, &i_face_normal);

    fvm_selector_t *sel_i_faces = fvm_selector_create(mesh->dim,
                                                      mesh->n_i_faces,
                                                      mesh->class_defs,
                                                      mesh->i_face_family,
                                                      1,
                                                      i_face_cog,
                                                      i_face_normal);

    c_id = fvm_selector_get_list(sel_i_faces,
                                 criteria,
                                 0,
                                 n_i_faces,
                                 i_face_list);

    BFT_FREE(i_face_cog);
    BFT_FREE(i_face_normal);

    if (del_class_defs)
      mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

    sel_i_faces = fvm_selector_destroy(sel_i_faces);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of cells verifying a given selection criteria.
 *
 * \param[in]   criteria   selection criteria string
 * \param[out]  n_cells    number of selected cells
 * \param[out]  cell_list  list of selected cells
 *                         (0 to n-1, preallocated to cs_glob_mesh->n_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cell_list(const char  *criteria,
                          cs_lnum_t   *n_cells,
                          cs_lnum_t    cell_list[])
{
  int c_id;

  *n_cells = 0;

  if (cs_glob_mesh->select_cells != NULL) {

    c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                                 criteria,
                                 0,
                                 n_cells,
                                 cell_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("The group \"%s\" in the selection criteria:\n"
                   "\"%s\"\n"
                   " does not correspond to any cell.\n"),
                 missing, criteria);
    }

  }

  else {

    cs_mesh_t *mesh = cs_glob_mesh;

    bool del_class_defs = (mesh->class_defs == NULL) ? true : false;

    cs_mesh_init_group_classes(mesh);

    cs_real_t  *i_face_cog = NULL, *i_face_normal = NULL;
    cs_real_t  *b_face_cog = NULL, *b_face_normal = NULL;
    cs_real_t  *cell_cen = NULL;
    BFT_MALLOC(cell_cen, mesh->n_cells_with_ghosts*3, cs_real_t);

    cs_mesh_quantities_i_faces(mesh, &i_face_cog, &i_face_normal);
    cs_mesh_quantities_b_faces(mesh, &b_face_cog, &b_face_normal);

    cs_mesh_quantities_cell_faces_cog(mesh,
                                      i_face_normal,
                                      i_face_cog,
                                      b_face_normal,
                                      b_face_cog,
                                      cell_cen);

    BFT_FREE(b_face_normal);
    BFT_FREE(b_face_cog);
    BFT_FREE(i_face_normal);
    BFT_FREE(i_face_cog);

    fvm_selector_t *sel_cells = fvm_selector_create(mesh->dim,
                                                    mesh->n_cells,
                                                    mesh->class_defs,
                                                    mesh->cell_family,
                                                    1,
                                                    cell_cen,
                                                    NULL);

    c_id = fvm_selector_get_list(sel_cells,
                                 criteria,
                                 0,
                                 n_cells,
                                 cell_list);

    BFT_FREE(cell_cen);

    if (del_class_defs)
      mesh->class_defs = fvm_group_class_set_destroy(mesh->class_defs);

    sel_cells = fvm_selector_destroy(sel_cells);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of cells verifying a given selection criteria.
 *
 * \param[in]   criteria    selection criteria string
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cell_vertices_list(const char  *criteria,
                                   cs_lnum_t   *n_vertices,
                                   cs_lnum_t    vtx_ids[])
{
  cs_lnum_t  n_cells = 0;
  cs_lnum_t  *cell_ids = NULL;

  BFT_MALLOC(cell_ids, cs_glob_mesh->n_cells, cs_lnum_t);

  cs_selector_get_cell_list(criteria, &n_cells, cell_ids);
  cs_selector_get_cell_vertices_list_by_ids(n_cells,
                                            cell_ids,
                                            n_vertices,
                                            vtx_ids);

  BFT_FREE(cell_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of vertices belonging to a given list of cells.
 *
 * \param[in]   n_cells     number of selected cells
 * \param[in]   cell_ids    ids of selected cells
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cell_vertices_list_by_ids(cs_lnum_t         n_cells,
                                          const cs_lnum_t   cell_ids[],
                                          cs_lnum_t        *n_vertices,
                                          cs_lnum_t         vtx_ids[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t _n_vertices = m->n_vertices;

  char *cell_flag;
  BFT_MALLOC(cell_flag, m->n_cells, char);

  for (cs_lnum_t i = 0; i < m->n_cells; i++)
    cell_flag[i] = 0;

  if (cell_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[cell_ids[i]] = 1;
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells; i++)
      cell_flag[i] = 1;
  }

  for (cs_lnum_t i = 0; i < _n_vertices; i++)
    vtx_ids[i] = -1;

  /* Now mark associated vertices using main connectivty
     (could be faster when some adjacencies are available) */

  for (cs_lnum_t i = 0; i < m->n_i_faces; i++) {
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t c_id = m->i_face_cells[i][j];
      if (c_id < m->n_cells) {
        if (cell_flag[c_id] != 0) {
          cs_lnum_t s_id = m->i_face_vtx_idx[i];
          cs_lnum_t e_id = m->i_face_vtx_idx[i+1];
          for (cs_lnum_t k = s_id; k < e_id; k++)
            vtx_ids[m->i_face_vtx_lst[k]] = 1;
        }
      }
    }
  }

  for (cs_lnum_t i = 0; i < m->n_b_faces; i++) {
    cs_lnum_t c_id = m->b_face_cells[i];
    if (cell_flag[c_id] != 0) {
      cs_lnum_t s_id = m->b_face_vtx_idx[i];
      cs_lnum_t e_id = m->b_face_vtx_idx[i+1];
      for (cs_lnum_t k = s_id; k < e_id; k++)
        vtx_ids[m->b_face_vtx_lst[k]] = 1;
    }
  }

  BFT_FREE(cell_flag);

  /* Now compact list */

  cs_lnum_t n = 0;
  for (cs_lnum_t i = 0; i < _n_vertices; i++) {
    if (vtx_ids[i] != -1) {
      vtx_ids[n] = i;
      n++;
    }
  }

  *n_vertices = n;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of vertices verifying a given boundary selection criteria.
 *
 * \param[in]   criteria    selection criteria string
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_vertices_list(const char *criteria,
                                     cs_lnum_t  *n_vertices,
                                     cs_lnum_t   vtx_ids[])
{

  cs_lnum_t  n_faces = 0;
  cs_lnum_t  *face_ids = NULL;

  BFT_MALLOC(face_ids, cs_glob_mesh->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list(criteria, &n_faces, face_ids);
  cs_selector_get_b_face_vertices_list_by_ids(n_faces,
                                              face_ids,
                                              n_vertices,
                                              vtx_ids);

  BFT_FREE(face_ids);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of vertices belonging to a given list of boundary faces.
 *
 * \param[in]   n_cells     number of selected cells
 * \param[in]   cell_ids    ids of selected cells
 * \param[out]  n_vertices  number of selected vertices
 * \param[out]  vtx_ids     list of selected vertices
 *                          (0 to n-1, preallocated to cs_glob_mesh->n_vertices)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_vertices_list_by_ids(cs_lnum_t         n_faces,
                                            const cs_lnum_t   face_ids[],
                                            cs_lnum_t        *n_vertices,
                                            cs_lnum_t         vtx_ids[])
{

  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t _n_vertices = m->n_vertices;

  for (cs_lnum_t i = 0; i < _n_vertices; i++)
    vtx_ids[i] = -1;

  /* Mark vertices */
  if (face_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t f_id = face_ids[i];
      cs_lnum_t s_id = m->b_face_vtx_idx[f_id];
      cs_lnum_t e_id = m->b_face_vtx_idx[f_id+1];
      for (cs_lnum_t k = s_id; k < e_id; k++)
        vtx_ids[m->b_face_vtx_lst[k]] = 1;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      cs_lnum_t s_id = m->b_face_vtx_idx[i];
      cs_lnum_t e_id = m->b_face_vtx_idx[i+1];
      for (cs_lnum_t k = s_id; k < e_id; k++)
        vtx_ids[m->b_face_vtx_lst[k]] = 1;
    }
  }

  /* Now compact list */

  cs_lnum_t n = 0;
  for (cs_lnum_t i = 0; i < _n_vertices; i++) {
    if (vtx_ids[i] != -1) {
      vtx_ids[n] = i;
      n++;
    }
  }

  *n_vertices = n;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill lists of faces at the boundary of a set of cells verifying
 *        a given selection criteria.
 *
 * \param[in]   criteria   selection criteria string
 * \param[out]  n_i_faces  number of selected interior faces
 * \param[out]  n_b_faces  number of selected boundary faces
 * \param[out]  i_face_id  list of selected interior faces
 *                         (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 * \param[out]  b_face_id  list of selected boundary faces
 *                         (0 to n-1, preallocated to cs_glob_mesh->n_b_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_cells_boundary(const char  *criteria,
                               cs_lnum_t   *n_i_faces,
                               cs_lnum_t   *n_b_faces,
                               cs_lnum_t    i_face_id[],
                               cs_lnum_t    b_face_id[])
{
  cs_lnum_t ii, n_cells;
  cs_lnum_t *cell_list, *cell_flag;

  const cs_mesh_t *mesh = cs_glob_mesh;

  /* Mark cells inside zone selection */

  BFT_MALLOC(cell_list, mesh->n_cells, cs_lnum_t);
  BFT_MALLOC(cell_flag, mesh->n_cells_with_ghosts, cs_lnum_t);

  for (ii = 0; ii < mesh->n_cells; ii++)
    cell_flag[ii] = 0;

  n_cells = 0;

  cs_selector_get_cell_list(criteria, &n_cells, cell_list);

  for (ii = 0; ii < n_cells; ii++)
    cell_flag[cell_list[ii]] = 1;

  BFT_FREE(cell_list);

  if (mesh->halo != NULL)
    cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, cell_flag);

  /* Now build lists of faces on cell boundaries */

  for (ii = 0; ii < mesh->n_i_faces; ii++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[ii][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[ii][1];
    if (cell_flag[c_id_0] != cell_flag[c_id_1]) {
      i_face_id[*n_i_faces] = ii;
      *n_i_faces += 1;
    }
  }

  for (ii = 0; ii < mesh->n_b_faces; ii++) {
    cs_lnum_t c_id = mesh->b_face_cells[ii];
    if (cell_flag[c_id] == 1) {
      b_face_id[*n_b_faces] = ii;
      *n_b_faces += 1;
    }
  }

  BFT_FREE(cell_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of cells attached to a given boundary selection criteria.
 *
 * Only cells sharing a face (not just a vertex) with the boundary are
 * selected.
 *
 * \param[in]   criteria     selection criteria string
 * \param[out]  n_b_cells    number of selected cells
 * \param[out]  b_cell_list  list of selected cells
 *                           (0 to n-1, preallocated to cs_glob_mesh->n_b_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_b_face_cells_list(const char  *criteria,
                                  cs_lnum_t   *n_b_cells,
                                  cs_lnum_t    b_cell_list[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;

  /* Get list of faces */
  cs_lnum_t  n_b_faces = 0;
  cs_lnum_t *b_face_list = NULL;

  BFT_MALLOC(b_face_list, mesh->n_b_faces, cs_lnum_t);

  cs_selector_get_b_face_list(criteria, &n_b_faces, b_face_list);

  /* Flag array initialization */
  int *cell_flag = NULL;

  BFT_MALLOC(cell_flag, mesh->n_cells, int);

  for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
    cell_flag[i] = 0;

  /* Loop on all found boundary faces.
   *  -> Check if cell was allready flagged (if attached to 2 or more faces)
   *  -> If not flagged, add cell id to list array, flag the cell and
   *     increment cell count.
   */
  *n_b_cells = 0;

  for (cs_lnum_t e_id = 0; e_id < n_b_faces; e_id++) {
    cs_lnum_t f_id = b_face_list[e_id];

    cs_lnum_t c_id = mesh->b_face_cells[f_id];

    if (cell_flag[c_id] == 0) {
      b_cell_list[*n_b_cells] = c_id;
      cell_flag[c_id] = 1;
      *n_b_cells += 1;
    }
  }

  BFT_FREE(b_face_list);
  BFT_FREE(cell_flag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of interior faces belonging to a given periodicity.
 *
 * \param[in]   perio_num  periodicity number
 * \param[out]  n_i_faces  number of selected interior faces
 * \param[out]  i_face_id  list of selected interior faces
 *                         (0 to n-1, preallocated to cs_glob_mesh->n_i_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_perio_face_list(int         perio_num,
                                cs_lnum_t  *n_i_faces,
                                cs_lnum_t   i_face_id[])
{
  int ii;
  int *face_perio_num = NULL;

  BFT_MALLOC(face_perio_num, cs_glob_mesh->n_i_faces, int);

  cs_mesh_get_face_perio_num(cs_glob_mesh, face_perio_num);

  *n_i_faces = 0;
  for (ii = 0; ii < cs_glob_mesh->n_i_faces; ii++) {
    if (CS_ABS(face_perio_num[ii]) == perio_num) {
      i_face_id[*n_i_faces] = ii;
      *n_i_faces += 1;
    }
  }

  BFT_FREE(face_perio_num);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Fill a list of families verifying a given selection criteria.
 *
 * \param[in]   criteria     selection criteria string
 * \param[out]  n_families   number of selected families
 * \param[out]  family_list  list of selected family ids (preallocated to
 *                           cs_glob_mesh->n_families + 1)
 */
/*----------------------------------------------------------------------------*/

void
cs_selector_get_family_list(const char  *criteria,
                            int         *n_families,
                            int          family_list[])
{
  int c_id;

  *n_families = 0;

  /* As all selectors were built with the same group class definitions,
     any selector may be used here. */
  c_id = fvm_selector_get_gc_list(cs_glob_mesh->select_cells,
                                  criteria,
                                  n_families,
                                  family_list);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("The group \"%s\" in the selection criteria:\n"
                 "\"%s\"\n"
                 " is not present in the mesh.\n"),
               missing, criteria);
  }

  /* Families seen in mesh are 1-based, while selector is 0-based */

  for (int i = 0; i < *n_families; i++)
    family_list[i] += 1;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
