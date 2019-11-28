/*============================================================================
 * Functions to remove mesh elements.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_interface.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_boundary.h"
#include "cs_mesh_halo.h"
#include "cs_mesh_group.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_remove.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove flagged cells.
 *
 * \param[in, out]  m            mesh
 * \param[in]       flag        cell flag (!= 0 to remove)
 * \param[in]       group_name  name of group to assign to new boundary faces,
 *                              or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_remove_cells(cs_mesh_t    *m,
                     char          flag[],
                     const char   *group_name)
{
  cs_lnum_t n_cells = m->n_cells;

  cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  cs_lnum_t n_i_faces = m->n_i_faces;
  cs_lnum_t n_b_faces = m->n_b_faces;
  cs_lnum_2_t *i_face_cells = (cs_lnum_2_t *)(m->i_face_cells);
  cs_lnum_t *b_face_cells = m->b_face_cells;

  cs_lnum_t  n_sel_faces = 0;
  cs_lnum_t *sel_faces;
  BFT_MALLOC(sel_faces, n_i_faces, cs_lnum_t);

  cs_lnum_t *c_o2n;
  BFT_MALLOC(c_o2n, n_cells_ext, cs_lnum_t);

  cs_lnum_t n_cells_new = 0;

  const int default_family_id = 1;

  /* Build old to new cells renumbering based on flags */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    if (flag[i])
      c_o2n[i] = -1;
    else
      c_o2n[i] = n_cells_new++;
  }

  cs_gnum_t n_g_cells_new = n_cells_new;
  cs_parall_counter(&n_g_cells_new, 1);

  if (n_g_cells_new == m->n_g_cells) {
    BFT_FREE(c_o2n);
    return;
  }

  if (m->verbosity > 0) {
    cs_log_printf(CS_LOG_DEFAULT, "\n");
    cs_log_separator(CS_LOG_DEFAULT);
    cs_log_printf(CS_LOG_DEFAULT,
                  _("Removing %llu cells from mesh\n\n"),
                  (unsigned long long)(m->n_g_cells - n_g_cells_new));
  }

  if (m->halo != NULL)
    cs_halo_sync_untyped(m->halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_lnum_t),
                         c_o2n);

  cs_mesh_free_rebuildable(m, true);

  /* Prepare for propagation of boundary face groups to interior faces
     of removed boundary cells; if the cell has multiple boundary faces
     with different groups, the first group wins */

  int *b_gc_id;
  BFT_MALLOC(b_gc_id, n_cells, cs_lnum_t);
  for (cs_lnum_t i = 0; i < n_cells; i++)
    b_gc_id[i] = default_family_id;

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t k = b_face_cells[i];
    if (flag[k])
      b_gc_id[k] = m->b_face_family[i];
  }

  /* Transform interior to boundary faces */

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    bool select = false;;
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t k = i_face_cells[i][j];
      if (k > -1) {
        if (c_o2n[k] == -1) {
          select = true;
          if (b_gc_id[k] != 0 && m->i_face_family[i] == default_family_id)
            m->i_face_family[i] = b_gc_id[k];
        }
      }
    }
    if (select)
      sel_faces[n_sel_faces++] = i;
  }

  BFT_FREE(b_gc_id);

  const cs_lnum_t n_b_faces_ini = m->n_b_faces;

  cs_mesh_boundary_insert_with_shared_vertices(m,
                                               n_sel_faces,
                                               sel_faces);

  BFT_FREE(sel_faces);

  i_face_cells = (cs_lnum_2_t *)(m->i_face_cells);
  b_face_cells = m->b_face_cells;
  n_i_faces = m->n_i_faces;
  n_b_faces = m->n_b_faces;

  /* Add group to new boundary faces */

  if (group_name != NULL) {

    cs_lnum_t n_b_add = m->n_b_faces - n_b_faces_ini;
    BFT_MALLOC(sel_faces, n_b_add, cs_lnum_t);
    cs_lnum_t k = 0;
    for (cs_lnum_t i = 0; i < n_b_add; i++) {
      cs_lnum_t j = n_b_add + i;
      if (m->b_face_family[j] == default_family_id)
        sel_faces[k++] = j;
    }

    cs_mesh_group_b_faces_set(m,
                              group_name,
                              k,
                              sel_faces);


    BFT_FREE(sel_faces);

  }

  /* Remove isolated boundary faces */

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t k = b_face_cells[i];
    if (k > -1) {
      if (c_o2n[k] < 0)
        b_face_cells[i] = -1;
    }
  }

  cs_gnum_t n_g_free_faces = 0;
  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t k = b_face_cells[i];
    if (k < 0)
      n_g_free_faces += 1;
  }
  cs_parall_counter(&n_g_free_faces, 1);
  m->n_g_free_faces = n_g_free_faces;

  cs_mesh_discard_free_faces(m);

  b_face_cells = m->b_face_cells;
  n_b_faces = m->n_b_faces;

  /* Update matching cell ids (i.e. discard removed cells) */

  for (cs_lnum_t i = 0; i < n_i_faces; i++) {
    for (cs_lnum_t j = 0; j < 2; j++) {
      cs_lnum_t k = i_face_cells[i][j];
      if (k > -1 && k < n_cells)
        i_face_cells[i][j] = c_o2n[k];
      else
        i_face_cells[i][j] = -1;
    }
  }

  for (cs_lnum_t i = 0; i < n_b_faces; i++) {
    cs_lnum_t k = b_face_cells[i];
    if (k > -1)
      b_face_cells[i] = c_o2n[k];
  }

  /* Now renumber cells (note that we avoid owerwriting useful values
     here only because renumbering was built in the same order) */

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    cs_lnum_t j = c_o2n[i];
    assert(i >= j);
    if (j > -1)
      m->cell_family[j] = m->cell_family[i];
  }

  if (m->global_cell_num != NULL) {
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      cs_lnum_t j = c_o2n[i];
      assert(i >= j);
      if (j > -1)
        m->global_cell_num[j] = m->global_cell_num[i];
    }
  }

  /* Update global numbering */

  if (m->global_cell_num != NULL || cs_glob_n_ranks > 1) {

    fvm_io_num_t *n_io_num
      = fvm_io_num_create_from_select(NULL, m->global_cell_num, n_cells_new, 0);

    BFT_FREE(m->global_cell_num);

    m->global_cell_num = fvm_io_num_transfer_global_num(n_io_num);
    m->n_g_cells = fvm_io_num_get_global_count(n_io_num);

    n_io_num = fvm_io_num_destroy(n_io_num);

    assert(m->n_g_cells == n_g_cells_new);

  }
  else
    m->n_g_cells = n_g_cells_new;

  m->n_cells = n_cells_new;
  m->n_cells_with_ghosts = n_cells_new;

  BFT_FREE(c_o2n);

  m->modified = CS_MAX(m->modified, 1);

  /* Rebuild ghosts */

  int mv_save = m->verbosity;
  m->verbosity = -1;

  if (   m->n_domains > 1 || m->n_init_perio > 0
      || m->halo_type == CS_HALO_EXTENDED) {
    cs_halo_type_t halo_type = m->halo_type;
    assert(m == cs_glob_mesh);
    cs_mesh_builder_t *mb = (m == cs_glob_mesh) ? cs_glob_mesh_builder : NULL;
    cs_mesh_init_halo(m, mb, halo_type);
  }

  cs_mesh_update_auxiliary(cs_glob_mesh);

  m->verbosity = mv_save;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Remove cells with negative volumes
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_remove_cells_negative_volume(cs_mesh_t  *m)
{
  cs_lnum_t n_cells = m->n_cells;

  cs_real_t *cell_vol = cs_mesh_quantities_cell_volume(m);

  cs_gnum_t n_neg = 0;

  for (cs_lnum_t i = 0; i < n_cells; i++) {
    if (cell_vol[i] <= 0)
      n_neg += 1;
  }

  cs_parall_counter(&n_neg, 1);

  if (n_neg > 0) {

    bft_printf(_("\n Will remove %llu cells with negative volume\n"),
               (unsigned long long)n_neg);

    char *flag;
    BFT_MALLOC(flag, m->n_cells, char);

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      if (cell_vol[i] <= 0)
        flag[i] = 1;
      else
        flag[i] = 0;
    }

    cs_mesh_remove_cells(m,
                         flag,
                         "[join_neg_volume]");

    BFT_FREE(flag);

  }

  BFT_FREE(cell_vol);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
