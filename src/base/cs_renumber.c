/*============================================================================
 * Optional mesh renumbering
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_defs.h"
#include "cs_prototypes.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_renumber.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_renumber.c
        Optional mesh renumbering.

  \enum cs_renumber_i_faces_type_t

  \brief Interior faces renumbering algorithm types

  \var CS_RENUMBER_I_FACES_BLOCK
       No shared cell in block.
       This should produce blocks of similar (prescribed) size across thread
       groups.
  \var CS_RENUMBER_I_FACES_MULTIPASS
       Use multipass face numbering.
       This should produce a smaller number of blocks, with a diminishing
       number of faces per thread group.

  \var CS_RENUMBER_I_FACES_NONE
       No interior face numbering.
*/

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_RENUMBER_N_SUBS  5  /* Number of categories for histograms */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* CSR (Compressed Sparse Row) graph representation */
/*--------------------------------------------------*/

/* Note that mesh cells correspond to graph vertices,
   and mesh faces to graph edges */

typedef struct {

  cs_lnum_t         n_rows;           /* Number of rows in CSR structure */
  cs_lnum_t         n_cols_max;       /* Maximum number of nonzero values
                                         on a given row */

  /* Pointers to structure arrays and info (row_index, col_id) */

  cs_lnum_t        *row_index;        /* Row index (0 to n-1) */
  cs_lnum_t        *col_id;           /* Column id (0 to n-1) */

} _csr_graph_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

int _cs_renumber_n_threads = 0;

cs_lnum_t  _min_i_subset_size = 64;
cs_lnum_t  _min_b_subset_size = 64;

cs_renumber_i_faces_type_t _i_faces_algorithm = CS_RENUMBER_I_FACES_MULTIPASS;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Redistribute family (group class) ids in case of renubering
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   n_elts   <--  Number of elements
 *   renum    <--  Pointer to renumbering array (new -> old, 1 to n)
 *   family   <->  Pointer to array of family ids (or NULL)
 *----------------------------------------------------------------------------*/

static void
_update_family(cs_int_t         n_elts,
               const cs_int_t  *renum,
               cs_int_t        *family)
{
  cs_int_t ii;
  cs_int_t *old_family;

  if (family == NULL)
    return;

  BFT_MALLOC(old_family, n_elts, cs_int_t);

  memcpy(old_family, family, n_elts*sizeof(cs_int_t));

  for (ii = 0; ii < n_elts; ii++)
    family[ii] = old_family[renum[ii] - 1];

  BFT_FREE(old_family);
}

/*----------------------------------------------------------------------------
 * Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_elts      --> number of elements in array
 *   init_num    --> initial local number of renumbered elements (1 to n)
 *   global_num  <-> global numbering (allocated if initially NULL)
 *----------------------------------------------------------------------------*/

static void
_update_global_num(size_t             n_elts,
                   const cs_lnum_t    init_num[],
                   cs_gnum_t        **global_num)
{
  size_t i;
  cs_gnum_t *_global_num = *global_num;

  if (_global_num == NULL) {

    BFT_MALLOC(_global_num, n_elts, cs_gnum_t);

    for (i = 0; i < n_elts; i++)
      _global_num[i] = init_num[i];

    *global_num = _global_num;
  }

  else {

    cs_gnum_t *tmp_global;

    BFT_MALLOC(tmp_global, n_elts, cs_gnum_t);
    memcpy(tmp_global, _global_num, n_elts*sizeof(cs_gnum_t));

    for (i = 0; i < n_elts; i++)
      _global_num[i] = tmp_global[init_num[i] - 1];

    BFT_FREE(tmp_global);
  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering of cells.
 *
 * parameters:
 *   mesh            <-> Pointer to global mesh structure
 *   renum           <-- Cells renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_cells(cs_mesh_t             *mesh,
                          const cs_int_t        *renum)
{
  cs_int_t  ii, jj, kk, face_id, n_vis, start_id, start_id_old;

  cs_int_t  *face_cells_tmp = NULL;
  cs_int_t  *new_cell_id = NULL;

  cs_int_t  face_cells_max_size = CS_MAX(mesh->n_i_faces*2, mesh->n_b_faces);
  const cs_int_t  n_cells = mesh->n_cells;

  /* If no renumbering is present, return */

  if (renum == NULL)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(face_cells_tmp, face_cells_max_size, cs_int_t);
  BFT_MALLOC(new_cell_id, mesh->n_cells_with_ghosts, cs_int_t);

  /* Build old -> new renumbering (1 to n) */

  for (ii = 0; ii < n_cells; ii++)
    new_cell_id[renum[ii] - 1] = ii;

  for (ii = n_cells; ii < mesh->n_cells_with_ghosts; ii++)
    new_cell_id[ii] = ii;

  /* Update halo connectivity */

  if (mesh->halo != NULL)
    cs_halo_renumber_cells(mesh->halo, new_cell_id);

  /* Update faces -> cells connectivity */

  memcpy(face_cells_tmp,
         mesh->i_face_cells,
         mesh->n_i_faces * 2 * sizeof(cs_int_t));

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    ii = face_cells_tmp[face_id*2] - 1;
    jj = face_cells_tmp[face_id*2 + 1] - 1;
    mesh->i_face_cells[face_id*2] = new_cell_id[ii] + 1;
    mesh->i_face_cells[face_id*2 + 1] = new_cell_id[jj] + 1;
  }

  if (mesh->n_b_faces > 0) {

    memcpy(face_cells_tmp,
           mesh->b_face_cells,
           mesh->n_b_faces * sizeof(cs_int_t));

    for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
      ii = face_cells_tmp[face_id] - 1;
      mesh->b_face_cells[face_id] = new_cell_id[ii] + 1;
    }
  }

  /* Update cell -> cells connectivity for extended neighborhood */

  if (mesh->cell_cells_lst != NULL) {

    cs_int_t *cell_cells_idx_old, *cell_cells_lst_old;
    const cs_int_t cell_cells_lst_size = mesh->cell_cells_idx[n_cells] - 1;

    BFT_MALLOC(cell_cells_idx_old, n_cells + 1, cs_int_t);
    BFT_MALLOC(cell_cells_lst_old, cell_cells_lst_size, cs_int_t);

    memcpy(cell_cells_idx_old,
           mesh->cell_cells_idx,
           (n_cells + 1)*sizeof(cs_int_t));
    memcpy(cell_cells_lst_old,
           mesh->cell_cells_lst,
           cell_cells_lst_size*sizeof(cs_int_t));

    mesh->cell_cells_idx[0] = 1;
    start_id = 0;

    for (ii = 0; ii < n_cells; ii++) {

      jj = renum[ii] - 1;
      n_vis = cell_cells_idx_old[jj+1] - cell_cells_idx_old[jj];
      start_id_old = cell_cells_idx_old[jj] - 1;

      for (kk = 0; kk < n_vis; kk++)
        mesh->cell_cells_lst[start_id + kk]
          = new_cell_id[cell_cells_lst_old[start_id_old + kk] - 1] + 1;

      start_id += n_vis;
      mesh->cell_cells_idx[ii + 1] = start_id + 1;
    }
  }

  /* Free work arrays */

  BFT_FREE(new_cell_id);
  BFT_FREE(face_cells_tmp);

  /* Update cell families and global numbering */

  _update_family(n_cells, renum, mesh->cell_family);

  _update_global_num(n_cells, renum, &(mesh->global_cell_num));

  /* Update parent cell numbers for post-processing meshes
     that may already have been built; Post-processing meshes
     built after renumbering will have correct parent numbers */

  cs_post_renum_cells(renum);
}

/*----------------------------------------------------------------------------
 * Apply renumbering to a face -> vertices connectivity.
 *
 * parameters:
 *   n_faces         <-- Number of faces
 *   face_vtx_idx    <-> Face -> vertices index (1 to n)
 *   face_vtx        <-- Face vertices
 *   renum           <-- Faces renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_update_face_vertices(cs_int_t         n_faces,
                      cs_int_t        *face_vtx_idx,
                      cs_int_t        *face_vtx,
                      const cs_int_t  *renum)
{
  if (renum != NULL && face_vtx != NULL) {

    cs_int_t ii, jj, kk, n_vtx, start_id, start_id_old;
    cs_int_t *face_vtx_idx_old, *face_vtx_old;

    const cs_int_t connect_size = face_vtx_idx[n_faces] - 1;

    BFT_MALLOC(face_vtx_idx_old, n_faces + 1, cs_int_t);
    BFT_MALLOC(face_vtx_old, connect_size, cs_int_t);

    memcpy(face_vtx_idx_old, face_vtx_idx, (n_faces+1)*sizeof(int));
    memcpy(face_vtx_old, face_vtx, connect_size*sizeof(int));

    face_vtx_idx[0] = 1;
    start_id = 0;

    for (ii = 0; ii < n_faces; ii++) {

      jj = renum[ii] - 1;
      n_vtx = face_vtx_idx_old[jj+1] - face_vtx_idx_old[jj];
      start_id_old = face_vtx_idx_old[jj] - 1;

      for (kk = 0; kk < n_vtx; kk++)
        face_vtx[start_id + kk] = face_vtx_old[start_id_old + kk];

      start_id += n_vtx;
      face_vtx_idx[ii + 1] = start_id + 1;
    }

    BFT_FREE(face_vtx_idx_old);
    BFT_FREE(face_vtx_old);
  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering of faces.
 *
 * parameters:
 *   mesh            <-> Pointer to global mesh structure
 *   renum_i         <-- Interior faces renumbering array (new -> old)
 *   renum_b         <-- Boundary faces renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_faces(cs_mesh_t             *mesh,
                          const cs_lnum_t       *renum_i,
                          const cs_lnum_t       *renum_b)
{
  cs_int_t  face_id, face_id_old;

  cs_int_t  *face_cells_old = NULL;

  const cs_int_t  n_i_faces = mesh->n_i_faces;
  const cs_int_t  n_b_faces = mesh->n_b_faces;

  /* Interior faces */

  if (renum_i != NULL) {

    /* Allocate Work array */

    BFT_MALLOC(face_cells_old, n_i_faces*2, cs_int_t);

    /* Update faces -> cells connectivity */

    memcpy(face_cells_old, mesh->i_face_cells, n_i_faces*2*sizeof(cs_int_t));

    for (face_id = 0; face_id < n_i_faces; face_id++) {
      face_id_old = renum_i[face_id] - 1;
      mesh->i_face_cells[face_id*2] = face_cells_old[face_id_old*2];
      mesh->i_face_cells[face_id*2 + 1] = face_cells_old[face_id_old*2 + 1];
    }

    BFT_FREE(face_cells_old);

    /* Update faces -> vertices connectivity */

    _update_face_vertices(n_i_faces,
                          mesh->i_face_vtx_idx,
                          mesh->i_face_vtx_lst,
                          renum_i);

    /* Update face families and global numbering */

    _update_family(n_i_faces, renum_i, mesh->i_face_family);

    _update_global_num(n_i_faces, renum_i, &(mesh->global_i_face_num));
  }

  /* Boundary faces */

  if (renum_b != NULL) {

    /* Allocate Work array */

    BFT_MALLOC(face_cells_old, n_b_faces, cs_int_t);

    /* Update faces -> cells connectivity */

    memcpy(face_cells_old, mesh->b_face_cells, n_b_faces*sizeof(cs_int_t));

    for (face_id = 0; face_id < n_b_faces; face_id++) {
      face_id_old = renum_b[face_id] - 1;
      mesh->b_face_cells[face_id] = face_cells_old[face_id_old];
    }

    BFT_FREE(face_cells_old);

    /* Update faces -> vertices connectivity */

    _update_face_vertices(n_b_faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          renum_b);

    /* Update face families and global numbering */

    _update_family(n_b_faces, renum_b, mesh->b_face_family);

    _update_global_num(n_b_faces, renum_b, &(mesh->global_b_face_num));
  }

  /* Update parent face numbers for post-processing meshes
     that may already have been built; Post-processing meshes
     built after renumbering will have correct parent numbers */

  cs_post_renum_faces(renum_i, renum_b);
}

/*----------------------------------------------------------------------------
 * Compute the minimum and the maximum of a vector (locally).
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *   min       --> minimum
 *   max       --> maximum
 *----------------------------------------------------------------------------*/

static void
_compute_local_minmax_gnum(cs_lnum_t        n_vals,
                           const cs_gnum_t  var[],
                           cs_gnum_t       *min,
                           cs_gnum_t       *max)
{
  cs_lnum_t  i;
  cs_gnum_t  _min = var[0], _max = var[0];

  for (i = 1; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

  if (min != NULL)  *min = _min;
  if (max != NULL)  *max = _max;
}

static void
_compute_local_minmax_double(cs_lnum_t        n_vals,
                             const double     var[],
                             double          *min,
                             double          *max)
{
  cs_lnum_t  i;
  double  _min = var[0], _max = var[0];

  for (i = 1; i < n_vals; i++) {
    _min = CS_MIN(_min, var[i]);
    _max = CS_MAX(_max, var[i]);
  }

  if (min != NULL)  *min = _min;
  if (max != NULL)  *max = _max;
}

/*----------------------------------------------------------------------------
 * Display the distribution of values of a vector.
 *
 * parameters:
 *   n_vals    <-- local number of elements
 *   var       <-- pointer to vector
 *----------------------------------------------------------------------------*/

static void
_display_histograms_gnum(int               n_vals,
                         const cs_gnum_t   var[])
{
  cs_gnum_t i, j, k;
  cs_gnum_t val_max, val_min;
  double step;

  cs_gnum_t count[CS_RENUMBER_N_SUBS];
  cs_gnum_t n_steps = CS_RENUMBER_N_SUBS;

  /* Compute local min and max */

  if (n_vals == 0) {
    bft_printf(_("    no value\n"));
    return;
  }

  val_max = var[0];
  val_min = var[0];
  _compute_local_minmax_gnum(n_vals, var, &val_min, &val_max);

  bft_printf(_("    minimum value =         %10llu\n"),
             (unsigned long long)val_min);
  bft_printf(_("    maximum value =         %10llu\n\n"),
             (unsigned long long)val_max);

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (val_max - val_min > 0) {

    if (val_max-val_min < n_steps)
      n_steps = CS_MAX(1, floor(val_max-val_min));

    step = (double)(val_max - val_min) / n_steps;

    /* Loop on values */

    for (i = 0; i < (cs_gnum_t)n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < val_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3llu : [ %10llu ; %10llu [ = %10llu\n",
                 (unsigned long long)(i+1),
                 (unsigned long long)(val_min + i*step),
                 (unsigned long long)(val_min + j*step),
                 (unsigned long long)(count[i]));

    bft_printf("    %3llu : [ %10llu ; %10llu ] = %10llu\n",
               (unsigned long long)n_steps,
               (unsigned long long)(val_min + (n_steps - 1)*step),
               (unsigned long long)val_max,
               (unsigned long long)(count[n_steps - 1]));

  }

  else { /* if (val_max == val_min) */
    bft_printf("    %3d : [ %10llu ; %10llu ] = %10llu\n",
               1, (unsigned long long)(val_min),
               (unsigned long long)val_max,
               (unsigned long long)n_vals);
  }

}

static void
_display_histograms_double(int           n_vals,
                           const double  var[])
{
  cs_gnum_t i, j, k;
  double val_max, val_min;
  double step;

  cs_gnum_t count[CS_RENUMBER_N_SUBS];
  cs_gnum_t n_steps = CS_RENUMBER_N_SUBS;

  /* Compute local min and max */

  if (n_vals == 0) {
    bft_printf(_("    no value\n"));
    return;
  }

  val_max = var[0];
  val_min = var[0];
  _compute_local_minmax_double(n_vals, var, &val_min, &val_max);

  bft_printf(_("    minimum value =         %10.5e\n"), val_min);
  bft_printf(_("    maximum value =         %10.5e\n\n"), val_max);

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    count[j] = 0;

  if (val_max - val_min > 0) {

    if (val_max-val_min < n_steps)
      n_steps = CS_MAX(1, floor(val_max-val_min));

    step = (double)(val_max - val_min) / n_steps;

    /* Loop on values */

    for (i = 0; i < (cs_gnum_t)n_vals; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (var[i] < val_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    %3d : [ %10.5e ; %10.5e [ = %10llu\n",
                 (int)(i+1),
                 val_min + i*step, val_min + j*step,
                 (unsigned long long)(count[i]));

    bft_printf("    %3d : [ %10.5e ; %10.5e ] = %10llu\n",
               (int)n_steps,
               val_min + (n_steps - 1)*step,
               val_max,
               (unsigned long long)(count[n_steps - 1]));

  }

  else { /* if (val_max == val_min) */
    bft_printf("    %3d : [ %10.5e ; %10.5e ] = %10llu\n",
               1, val_min, val_max,
               (unsigned long long)n_vals);
  }

}

#if defined(HAVE_IBM_RENUMBERING_LIB)

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces and cells for multiple threads.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_for_threads_ibm(cs_mesh_t  *mesh)
{
}

#endif /* defined(HAVE_IBM_RENUMBERING_LIB) */

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of a cs_lnum_t (integer) array.
 *
 * parameters:
 *   number    <-> pointer to elements that should be ordered
 *   level     <-- level of the binary tree to descend
 *   n_elts    <-- number of elements in the binary tree to descend
 *----------------------------------------------------------------------------*/

inline static void
_sort_descend_tree(cs_lnum_t  number[],
                   size_t     level,
                   size_t     n_elts)
{
  size_t lv_cur;
  cs_lnum_t num_save;

  num_save = number[level];

  while (level <= (n_elts/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_elts - 1)
      if (number[lv_cur+1] > number[lv_cur]) lv_cur++;

    if (lv_cur >= n_elts) break;

    if (num_save >= number[lv_cur]) break;

    number[level] = number[lv_cur];
    level = lv_cur;

  }

  number[level] = num_save;
}

/*----------------------------------------------------------------------------
 * Order an array of global numbers.
 *
 * parameters:
 *   number   <-> number of arrays to sort
 *   n_elts   <-- number of elements considered
 *----------------------------------------------------------------------------*/

static void
_sort_local(cs_lnum_t  number[],
            size_t     n_elts)
{
  size_t i, j, inc;
  cs_lnum_t num_save;

  if (n_elts < 2)
    return;

  /* Use shell sort for short arrays */

  if (n_elts < 20) {

    /* Compute increment */
    for (inc = 1; inc <= n_elts/9; inc = 3*inc+1);

    /* Sort array */
    while (inc > 0) {
      for (i = inc; i < n_elts; i++) {
        num_save = number[i];
        j = i;
        while (j >= inc && number[j-inc] > num_save) {
          number[j] = number[j-inc];
          j -= inc;
        }
        number[j] = num_save;
      }
      inc = inc / 3;
    }

  }

  else {

    /* Create binary tree */

    i = (n_elts / 2);
    do {
      i--;
      _sort_descend_tree(number, i, n_elts);
    } while (i > 0);

    /* Sort binary tree */

    for (i = n_elts - 1 ; i > 0 ; i--) {
      num_save   = number[0];
      number[0] = number[i];
      number[i] = num_save;
      _sort_descend_tree(number, 0, i);
    }
  }
}

/*----------------------------------------------------------------------------
 * Create a CSR graph structure from a native face-based connectivity.
 *
 * parameters:
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity (1 to n)
 *
 * returns:
 *   pointer to allocated CSR graph structure.
 *----------------------------------------------------------------------------*/

static _csr_graph_t *
_csr_graph_create(cs_lnum_t         n_cells_ext,
                  cs_lnum_t         n_faces,
                  const cs_lnum_t  *face_cell)
{
  int n_cols_max;
  cs_lnum_t ii, jj, f_id;

  cs_lnum_t  *ccount = NULL;
  bool unique_faces = true;

  _csr_graph_t  *g;

  /* Allocate and map */

  BFT_MALLOC(g, 1, _csr_graph_t);

  g->n_rows = n_cells_ext;

  BFT_MALLOC(g->row_index, g->n_rows + 1, cs_lnum_t);

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, g->n_rows, cs_lnum_t);

  for (ii = 0; ii < g->n_rows; ii++)
    ccount[ii] = 0;

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2] - 1;
    jj = face_cell[f_id*2 + 1] - 1;
    ccount[ii] += 1;
    ccount[jj] += 1;
  }

  n_cols_max = 0;

  g->row_index[0] = 0;
  for (ii = 0; ii < g->n_rows; ii++) {
    g->row_index[ii+1] = g->row_index[ii] + ccount[ii];
    if (ccount[ii] > n_cols_max)
      n_cols_max = ccount[ii];
    ccount[ii] = 0;
  }

  g->n_cols_max = n_cols_max;

  /* Build structure */

  BFT_MALLOC(g->col_id, (g->row_index[g->n_rows]), cs_lnum_t);

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2] - 1;
    jj = face_cell[f_id*2 + 1] - 1;
    g->col_id[g->row_index[ii] + ccount[ii]] = jj;
    ccount[ii] += 1;
    g->col_id[g->row_index[jj] + ccount[jj]] = ii;
    ccount[jj] += 1;
  }

  BFT_FREE(ccount);

  /* Sort line elements by column id (for better access patterns) */

  if (n_cols_max > 1) {

    for (ii = 0; ii < g->n_rows; ii++) {
      cs_lnum_t *col_id = g->col_id + g->row_index[ii];
      cs_lnum_t n_cols = g->row_index[ii+1] - g->row_index[ii];
      cs_lnum_t col_id_prev = -1;
      _sort_local(col_id, g->row_index[ii+1] - g->row_index[ii]);
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id[jj] == col_id_prev)
          unique_faces = false;
        col_id_prev = col_id[jj];
      }
    }

  }

  /* Compact elements if necessary */

  if (unique_faces == false) {

    cs_lnum_t *tmp_row_index = NULL;
    cs_lnum_t  kk = 0;

    BFT_MALLOC(tmp_row_index, g->n_rows+1, cs_lnum_t);
    memcpy(tmp_row_index, g->row_index, (g->n_rows+1)*sizeof(cs_lnum_t));

    kk = 0;

    for (ii = 0; ii < g->n_rows; ii++) {
      cs_lnum_t *col_id = g->col_id + g->row_index[ii];
      cs_lnum_t n_cols = g->row_index[ii+1] - g->row_index[ii];
      cs_lnum_t col_id_prev = -1;
      g->row_index[ii] = kk;
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          g->col_id[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    g->row_index[g->n_rows] = kk;

    assert(g->row_index[g->n_rows] < tmp_row_index[g->n_rows]);

    BFT_FREE(tmp_row_index);
    BFT_REALLOC(g->col_id, (g->row_index[g->n_rows]), cs_lnum_t);

  }

  return g;
}

/*----------------------------------------------------------------------------
 * Create a CSR graph structure from a native face-based conectivity.
 *
 * parameters:
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity (1 to n)
 *
 * returns:
 *   pointer to allocated CSR graph structure.
 *----------------------------------------------------------------------------*/

static _csr_graph_t *
_csr_graph_create_cell_face(cs_lnum_t         n_cells_ext,
                            cs_lnum_t         n_faces,
                            const cs_lnum_t  *face_cell)
{
  int n_cols_max;
  cs_lnum_t ii, jj, f_id;

  cs_lnum_t  *ccount = NULL;

  _csr_graph_t  *g;

  /* Allocate and map */

  BFT_MALLOC(g, 1, _csr_graph_t);

  g->n_rows = n_cells_ext;

  BFT_MALLOC(g->row_index, g->n_rows + 1, cs_lnum_t);

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, g->n_rows, cs_lnum_t);

  for (ii = 0; ii < g->n_rows; ii++)
    ccount[ii] = 0;

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2] - 1;
    jj = face_cell[f_id*2 + 1] - 1;
    ccount[ii] += 1;
    ccount[jj] += 1;
  }

  n_cols_max = 0;

  g->row_index[0] = 0;
  for (ii = 0; ii < g->n_rows; ii++) {
    g->row_index[ii+1] = g->row_index[ii] + ccount[ii];
    if (ccount[ii] > n_cols_max)
      n_cols_max = ccount[ii];
    ccount[ii] = 0;
  }

  g->n_cols_max = n_cols_max;

  /* Build structure */

  BFT_MALLOC(g->col_id, (g->row_index[g->n_rows]), cs_lnum_t);

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2] - 1;
    jj = face_cell[f_id*2 + 1] - 1;
    g->col_id[g->row_index[ii] + ccount[ii]] = f_id;
    ccount[ii] += 1;
    g->col_id[g->row_index[jj] + ccount[jj]] = f_id;
    ccount[jj] += 1;
  }

  BFT_FREE(ccount);

  return g;
}

/*----------------------------------------------------------------------------
 * Destroy CSR graph structure.
 *
 * parameters:
 *   g  <->  Pointer to CSR graph structure pointer
 *----------------------------------------------------------------------------*/

static void
_csr_graph_destroy(_csr_graph_t  **graph)
{
  if (graph != NULL && *graph !=NULL) {

    _csr_graph_t  *g = *graph;

    if (g->row_index != NULL)
      BFT_FREE(g->row_index);

    if (g->col_id != NULL)
      BFT_FREE(g->col_id);

    BFT_FREE(g);

    *graph = g;

  }
}

/*----------------------------------------------------------------------------
 * Build groups including independent faces.
 *
 * parameters:
 *   max_group_size  <-- max group size
 *   n_faces         <-- number of faces
 *   n_cells_ext     <-- local number of cells + ghost cells sharing a face
 *   n_faces         <-- local number of faces
 *   face_cell       <-- face -> cells connectivity (1 to n)
 *   new_to_old      --> new -> old face renumbering (1-based)
 *   n_groups        --> number of groups
 *   group_size      --> array containing the sizes of groups
 *----------------------------------------------------------------------------*/

static void
_independent_face_groups(cs_lnum_t          max_group_size,
                         cs_lnum_t          n_cells_ext,
                         cs_lnum_t          n_faces,
                         const cs_lnum_t   *face_cell,
                         cs_lnum_t         *new_to_old,
                         cs_lnum_t         *n_groups,
                         cs_lnum_t        **group_size)
{
  cs_lnum_t f_id, i, j, k;
  cs_lnum_t *group_face_ids = NULL, *face_marker = NULL;
  cs_lnum_t *old_to_new = NULL;
  _csr_graph_t *cell_faces = NULL;

  cs_lnum_t first_unmarked_face_id = 0;
  cs_lnum_t _n_groups_max = 4;
  cs_lnum_t n_marked_faces = 0;
  cs_lnum_t group_id = 0;

  cs_lnum_t *_group_size = NULL;

  BFT_MALLOC(_group_size, _n_groups_max, cs_lnum_t);

  BFT_MALLOC(old_to_new, n_faces, cs_lnum_t);
  BFT_MALLOC(face_marker, n_faces, cs_lnum_t);
  BFT_MALLOC(group_face_ids, max_group_size, cs_lnum_t);

  /* Create CSR cells -> faces graph */

  cell_faces = _csr_graph_create_cell_face(n_cells_ext,
                                           n_faces,
                                           face_cell);

  /* mark cell in a group */

  for (f_id = 0; f_id < n_faces; f_id++) {
    face_marker[f_id] = -1;
    old_to_new[f_id] = f_id;
  }

  while (n_marked_faces != n_faces) {

    cs_lnum_t  g_size = 0;

    /* Start a new group */

    for (f_id = 0; f_id < max_group_size; f_id++)
      group_face_ids[f_id] = -1;

    for (f_id = first_unmarked_face_id; f_id < n_faces; f_id++) {

      /* Search for next free face and check if it can be added
         in the current group */

      if (face_marker[f_id] == -1) {

        bool f_ok = true;

        for (i = 0; i < g_size && f_ok; i++) {

          cs_lnum_t f_cmp = group_face_ids[i];
          cs_lnum_t c_id[2] = {face_cell[f_cmp*2] - 1,
                               face_cell[f_cmp*2 + 1] - 1};

          for (j = 0; j < 2; j++) {
            cs_lnum_t start_id = cell_faces->row_index[c_id[j]];
            cs_lnum_t end_id = cell_faces->row_index[c_id[j] + 1];
            for (k = start_id; k < end_id; k++) {
              if (cell_faces->col_id[k] == f_id) {
                f_ok = false;
                break;
              }
            }
          }

        }

        /* Add the face to the group */

        if (f_ok == 1) {
          if (first_unmarked_face_id == f_id)
            first_unmarked_face_id = f_id + 1;
          face_marker[f_id] = group_id;
          group_face_ids[g_size++] = f_id;
          old_to_new[f_id] = n_marked_faces++;
        }

        /* Prepare to start new group if complete */

        if (g_size == max_group_size)
          break;

      } /* End of test on face_marker */

    } /* End of loop on faces */

    if (group_id + 1 >= _n_groups_max) {
      _n_groups_max *= 2;
      BFT_REALLOC(_group_size, _n_groups_max, cs_lnum_t);
    }
    _group_size[group_id++] = g_size;

  }

  _csr_graph_destroy(&cell_faces);

  BFT_FREE(face_marker);
  BFT_FREE(group_face_ids);

  BFT_REALLOC(_group_size, group_id, cs_lnum_t);

  /* Set return values */

  for (f_id = 0; f_id < n_faces; f_id++)
    new_to_old[old_to_new[f_id]] = f_id + 1;

  BFT_FREE(old_to_new);

  *n_groups = group_id;
  *group_size = _group_size;
}

/*----------------------------------------------------------------------------
 * Compute bounds for groups threads using only group sizes and
 * face renumbering
 *
 * parameters:
 *   n_faces         <-- local number of faces
 *   n_groups        <-- number of groups
 *   group_size      <-- array containing the sizes of groups
 *   group_index     --> index for groups
 *
 * returns:
 *   0 on success, -1 otherwise
 *----------------------------------------------------------------------------*/

static int
_thread_bounds_by_group_size(cs_lnum_t   n_faces,
                             int         n_groups,
                             int         n_threads,
                             cs_lnum_t  *group_size,
                             cs_lnum_t  *group_index)
{
  cs_lnum_t  group_id;
  cs_lnum_t  j, jr, k;

  cs_lnum_t ip = 0;
  cs_lnum_t stride = 2*n_groups;

  for (group_id = 0; group_id < n_groups; group_id++) {

    cs_lnum_t _group_size = group_size[group_id];

    j  = _group_size / n_threads;
    jr = _group_size % n_threads;

    if (j > 4) {
      for (k=0; k < n_threads; k++) {
        group_index[k*stride + group_id*2] = ip;
        ip += j;
        if (k < jr)
          ip++;
        group_index[k*stride + group_id*2+1] = ip;
      }
    }
    else {
      /* only thread 0 has elements */
      k = 0;
      group_index[k*stride + group_id*2] = ip;
      ip += _group_size;
      group_index[k*stride + group_id*2+1] = ip;
      for (k = 1; k < n_threads; k++) {
        group_index[k*stride + group_id*2]  = 0;
        group_index[k*stride + group_id*2+1] = 0;
      }
    }

  }

  if (ip != n_faces)
    return -1;

  return 0;
}

/*----------------------------------------------------------------------------
 * Pre-assign faces to threads of a given group for the multipass
 * algorithm, so as to improve load balance.
 *
 * parameters:
 *   n_i_threads     <-- number of threads required for interior faces
 *   n_g_i_threads   <-- number of threads active for interior faces for
 *                       this group
 *   g_id            <-- id of current threads group
 *   faces_list_size <-- size of list of faces to handle
 *   faces_list      <-- list of faces to handle, in lexicographical order
 *   l_face_cells    <-- face->cells connectivity,
 *                       with l_face_cells[i][0] < l_face_cells[i][0]
 *   f_t_id          <-> thread ids associated with interior faces
 *                       (local thread_id + g_id*n_ithreads, or -1 if
 *                       not determined yet)
 *   n_t_faces       --> number of faces associated with a given thread
 *   t_face_last     --> last face list if for a given thread
 *   t_cell_index    <-- index of starting and past-the end cell ids for
 *                       a given thread
 *----------------------------------------------------------------------------*/

static void
_renum_face_multipass_assign(int                         n_i_threads,
                             int                         n_g_i_threads,
                             int                         g_id,
                             cs_lnum_t                   faces_list_size,
                             const cs_lnum_t   *restrict faces_list,
                             const cs_lnum_2_t *restrict l_face_cells,
                             int               *restrict f_t_id,
                             cs_lnum_t         *restrict n_t_faces,
                             cs_lnum_t         *restrict t_face_last,
                             const cs_lnum_t   *restrict t_cell_index)
{
  int t_id;
  cs_lnum_t fl_id, f_id, c_id_0, c_id_1;

  for (t_id = 0; t_id < n_g_i_threads; t_id++) {
    n_t_faces[t_id] = 0;
    t_face_last[t_id] = faces_list_size;
  }

  t_id = 0;

  for (fl_id = 0; fl_id < faces_list_size; fl_id++) {

    f_id = faces_list[fl_id];

    c_id_0 = l_face_cells[f_id][0];
    c_id_1 = l_face_cells[f_id][1];

    /* determine thread possibly associated to this face */

    while (c_id_0 >= t_cell_index[t_id+1])
      t_id += 1;

    assert(t_id <= n_g_i_threads);

    if (   c_id_0 >= t_cell_index[t_id]
        && c_id_1 < t_cell_index[t_id+1]) {
      f_t_id[f_id] = t_id + g_id*n_i_threads;
      n_t_faces[t_id] += 1;
      t_face_last[t_id] = fl_id;
    }
    else
      f_t_id[f_id] = -1;

  }
}

/*----------------------------------------------------------------------------
 * Estimate unbalance between threads of a given group for the multipass
 * algorithm.
 *
 * Unbalance is considered to be: (max/mean - 1)
 *
 * parameters:
 *   n_i_threads     <-- number of threads required for interior faces
 *   n_t_faces       <-- number of faces associated with a given thread
 *
 * returns:
 *   estimated unbalance for this group
 *----------------------------------------------------------------------------*/

static double
_renum_face_multipass_g_unbalance(int               n_i_threads,
                                  const cs_lnum_t  *n_t_faces)
{
  int t_id;
  double n_t_faces_mean, imbalance;

  cs_lnum_t n_t_faces_sum = 0;
  cs_lnum_t n_t_faces_max = 0;

  for (t_id = 0; t_id < n_i_threads; t_id++) {
    n_t_faces_sum += n_t_faces[t_id];
    if (n_t_faces[t_id] > n_t_faces_max)
      n_t_faces_max = n_t_faces[t_id];
  }

  n_t_faces_mean = (double)n_t_faces_sum / n_i_threads;

  imbalance = (n_t_faces_max / n_t_faces_mean) - 1.0;

  return imbalance;
}

/*----------------------------------------------------------------------------
 * Redistribute faces between threads of a given group for the multipass
 * algorithm, so as to improve load balance.
 *
 * parameters:
 *   n_i_threads     <-- number of threads required for interior faces
 *   n_g_i_threads   <-- number of threads active for interior faces for
 *                       this group
 *   g_id            <-- id of current threads group
 *   relax           <-- relaxation factor
 *   faces_list_size <-- size of list of faces to handle
 *   faces_list      <-- list of faces to handle, in lexicographical order
 *   l_face_cells    <-- face->cells connectivity,
 *                       with l_face_cells[i][0] < l_face_cells[i][0]
 *   f_t_id          <-> thread ids associated with interior faces
 *                       (local thread_id + g_id*n_ithreads, or -1 if
 *                       not determined yet)
 *   n_t_faces       <-> number of faces associated with a given thread
 *   t_face_last     <-- last face list if for a given thread
 *   t_cell_index    <-> index of starting and past-the end cell ids for
 *                       a given thread
 *----------------------------------------------------------------------------*/

static void
_renum_face_multipass_redistribute(int                         n_i_threads,
                                   int                         n_g_i_threads,
                                   int                         g_id,
                                   double                      relax,
                                   cs_lnum_t                   faces_list_size,
                                   const cs_lnum_t   *restrict faces_list,
                                   const cs_lnum_2_t *restrict l_face_cells,
                                   int               *restrict f_t_id,
                                   cs_lnum_t         *restrict n_t_faces,
                                   cs_lnum_t         *restrict t_face_last,
                                   cs_lnum_t         *restrict t_cell_index)
{
  int t_id, t_id1;
  double unbalance[2];
  double n_t_faces_mean = 0.0;

  cs_lnum_t *t_cell_index_prev = NULL;

  if (n_g_i_threads < 2)
    return;

  /* Save previous cell index to allow reversal */

  BFT_MALLOC(t_cell_index_prev, n_g_i_threads+1, cs_lnum_t);
  memcpy(t_cell_index_prev,
         t_cell_index,
         (n_g_i_threads+1)*sizeof(cs_lnum_t));

  /* Estimate initial imbalance */

  unbalance[0] = _renum_face_multipass_g_unbalance(n_g_i_threads,
                                                   n_t_faces);

  /* Now ty to improve balancing */

  for (t_id = 0; t_id < n_g_i_threads; t_id++)
    n_t_faces_mean += n_t_faces[t_id];

  n_t_faces_mean /= n_g_i_threads;

  for (t_id = 0; t_id < n_g_i_threads - 1; t_id++) {

    t_id1 = t_id+1;

    cs_lnum_t t0_c_start = t_cell_index[t_id];
    cs_lnum_t t1_c_start = t_cell_index[t_id1];
    cs_lnum_t t1_c_end = t_cell_index[t_id1+1];

    cs_lnum_t n_t_faces_target = n_t_faces_mean; /* double to int */
    cs_lnum_t n_t_faces_move = n_t_faces[t_id] - n_t_faces_target;

    cs_lnum_t fl_id_end = t_face_last[t_id];

    n_t_faces_move *= relax;

    /* If t_id has too many edges, try to shift thread boundary back */

    if (n_t_faces_move > 0) {

      int f_t_id0 = t_id + g_id*n_i_threads;

      for (fl_id_end = t_face_last[t_id] - 1;
           (    fl_id_end > -1
             && l_face_cells[faces_list[fl_id_end]][0] >= t0_c_start
             && n_t_faces_move > 0);
           fl_id_end--) {
        if (f_t_id[faces_list[fl_id_end]] == f_t_id0)
          n_t_faces_move -= 1;
      }

      while (   fl_id_end < t_face_last[t_id]
             && (   l_face_cells[faces_list[fl_id_end+1]][0]
                 == l_face_cells[faces_list[fl_id_end]][0]))
        fl_id_end++;

      t_cell_index[t_id1] = l_face_cells[faces_list[fl_id_end]][0] + 1;
      if (t_cell_index[t_id1] > t1_c_start)
        t_cell_index[t_id1] = t1_c_start;

      t1_c_start = t_cell_index[t_id1];

    }

    /* If t_id has too few edges, try to shift thread boundary forward. */

    else if (n_t_faces_move < 0) {

      /* We assume the number of faces "removed" from the following
         thread is close to the number that will be gained by the
         current thread. */

      int f_t_id1 = t_id1 + g_id*n_i_threads;
      cs_lnum_t fl_id_max = CS_MIN(t_face_last[t_id1], faces_list_size - 1);

      for (fl_id_end = t_face_last[t_id];
            (   fl_id_end <= fl_id_max
             && l_face_cells[faces_list[fl_id_end]][0] <= t1_c_end
             && n_t_faces_move < 0);
           fl_id_end++) {
        if (f_t_id[faces_list[fl_id_end]] == f_t_id1)
          n_t_faces_move += 1;
      }

      fl_id_end = CS_MIN(fl_id_end, faces_list_size - 1);

      while (   fl_id_end >= t_face_last[t_id]
             && fl_id_end > 0
             && (   l_face_cells[faces_list[fl_id_end]][0]
                 == l_face_cells[faces_list[fl_id_end-1]][0]))
        fl_id_end--;

      t_cell_index[t_id1] = l_face_cells[faces_list[fl_id_end]][0];
      if (t_cell_index[t_id1] < t0_c_start)
        t_cell_index[t_id1] = t0_c_start;

      t1_c_start = t_cell_index[t_id1];

    }

  }

  /* Now reassign threads to faces */

  _renum_face_multipass_assign(n_i_threads,
                               n_g_i_threads,
                               g_id,
                               faces_list_size,
                               faces_list,
                               l_face_cells,
                               f_t_id,
                               n_t_faces,
                               t_face_last,
                               t_cell_index);

  unbalance[1] = _renum_face_multipass_g_unbalance(n_g_i_threads,
                                                   n_t_faces);

  /* If redistribution has degraded balancing (probably due to a too
     high relaxation factor value) , revert to initial distribution. */

  if (unbalance[1] > unbalance[0]) {

    memcpy(t_cell_index,
           t_cell_index_prev,
           (n_g_i_threads+1)*sizeof(cs_lnum_t));

    _renum_face_multipass_assign(n_i_threads,
                                 n_g_i_threads,
                                 g_id,
                                 faces_list_size,
                                 faces_list,
                                 l_face_cells,
                                 f_t_id,
                                 n_t_faces,
                                 t_face_last,
                                 t_cell_index);

  }

  BFT_FREE(t_cell_index_prev);
}

/*----------------------------------------------------------------------------
 * Update local face->cells connnectivity for multiple pass algorithm.
 *
 * Cells are marked and locally renumbered, so that only cells adjacent
 * to faces not yet assigned to a thread group are considered.
 *
 * parameters:
 *   n_f_cells_prev  <-- input number of cells adjacent to remaining faces
 *   faces_list_size <-- size of remaining faces to handle
 *   faces_list      <-- list of faces to handle, in lexicographical order
 *   l_face_cells    <-> face->cells connectivity,
 *                       with l_face_cells[i][0] < l_face_cells[i][0]
 *
 * returns:
 *   new number of cells adjacent to remaining faces
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_renum_face_multipass_remaining(cs_lnum_t                   n_f_cells_prev,
                                cs_lnum_t                   faces_list_size,
                                const cs_lnum_t   *restrict faces_list,
                                cs_lnum_2_t       *restrict l_face_cells)
{
  cs_lnum_t  fl_id, f_id, c_id_0, c_id_1;
  cs_lnum_t  n_f_cells_new = 0;
  cs_lnum_t *new_cell_id;

  BFT_MALLOC(new_cell_id, n_f_cells_prev, cs_lnum_t);

  for (c_id_0 = 0; c_id_0 < n_f_cells_prev; c_id_0++)
    new_cell_id[c_id_0] = -1;

  for (fl_id = 0; fl_id < faces_list_size; fl_id++) {

    f_id = faces_list[fl_id];

    c_id_0 = l_face_cells[f_id][0];
    c_id_1 = l_face_cells[f_id][1];

    if (new_cell_id[c_id_0] < 0)
      new_cell_id[c_id_0] = n_f_cells_new++;
    if (new_cell_id[c_id_1] < 0)
      new_cell_id[c_id_1] = n_f_cells_new++;

    if (new_cell_id[c_id_0] < new_cell_id[c_id_1]) {
      l_face_cells[f_id][0] = new_cell_id[c_id_0];
      l_face_cells[f_id][1] = new_cell_id[c_id_1];
    }
    else {
      l_face_cells[f_id][0] = new_cell_id[c_id_1];
      l_face_cells[f_id][1] = new_cell_id[c_id_0];
    }

  }

  BFT_FREE(new_cell_id);

  return n_f_cells_new;
}

/*----------------------------------------------------------------------------
 * Build groups including independent faces, using multiple pass algorithm
 *
 * Note: this function tries to optimize load balance between threads of
 *       a same group. It may be better to ensure that cells adjacent to
 *       faces of a same thread for a given group do not belong to a same
 *       cache line. This is not easy, so simply enforcing a minimum
 *       subset size for threads may be the simples approach.
 *
 * parameters:
 *   mesh           <-> pointer to global mesh structure
 *   n_i_threads    <-- number of threads required for interior faces
 *   group_size     <-- target group size
 *   renum_i        --> interior faces renumbering array (new -> old, 1-based)
 *   n_groups       --> number of groups of graph edges (interior faces)
 *   group_index    --> group/thread index
 *
 * returns:
 *   0 on success, -1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_face_multipass(cs_mesh_t    *mesh,
                      int           n_i_threads,
                      cs_lnum_t     renum_i[],
                      cs_lnum_t    *n_groups,
                      cs_lnum_t   **group_index)
{
  int g_id, t_id;
  cs_lnum_t fl_id, f_id, c_id_0, c_id_1;

  cs_lnum_t n_f_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_faces = mesh->n_i_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;

  double redistribute_relaxation_factor = 0.5;

  cs_lnum_t faces_list_size = n_faces, faces_list_size_new = 0;

  cs_lnum_t   _n_groups = 0;
  cs_lnum_t  *_group_index = NULL;

  cs_lnum_t *restrict faces_keys = NULL;
  cs_lnum_t *restrict faces_list = NULL;
  cs_lnum_2_t *restrict l_face_cells = NULL;
  cs_lnum_t *n_t_faces = NULL;
  cs_lnum_t *t_face_last = NULL;
  cs_lnum_t *t_cell_index = NULL;
  int *f_t_id = NULL;

  if (faces_list_size <= _min_i_subset_size)
    return -1;

  /* Initialization */

  BFT_MALLOC(faces_list, n_faces, cs_lnum_t);
  BFT_MALLOC(l_face_cells, n_faces, cs_lnum_2_t);
  BFT_MALLOC(n_t_faces, n_i_threads, cs_lnum_t);
  BFT_MALLOC(t_face_last, n_i_threads, cs_lnum_t);
  BFT_MALLOC(t_cell_index, n_i_threads + 1, cs_lnum_t);
  BFT_MALLOC(f_t_id, n_faces, int);

  /* Build lexical ordering of faces */

# pragma omp parallel for private(c_id_0, c_id_1)
  for (f_id = 0; f_id < n_faces; f_id++) {
    c_id_0 = i_face_cells[f_id][0] - 1;
    c_id_1 = i_face_cells[f_id][1] - 1;
    if (c_id_0 < c_id_1) {
      l_face_cells[f_id][0] = c_id_0;
      l_face_cells[f_id][1] = c_id_1;
    }
    else {
      l_face_cells[f_id][0] = c_id_1;
      l_face_cells[f_id][1] = c_id_0;
    }
    f_t_id[f_id] = -1;
  }

  cs_order_lnum_allocated_s(NULL,
                            (cs_lnum_t *)l_face_cells,
                            2,
                            faces_list,
                            n_faces);

  /* Add groups as required */

  for (g_id = 0; faces_list_size > _min_i_subset_size; g_id++) {

    int group_size = n_f_cells / n_i_threads;
    int n_g_i_threads = n_i_threads;

    /* Reduce number of threads for this level if required to
       ensure sufficient work per thread */

    if (faces_list_size / _min_i_subset_size  < n_g_i_threads) {
      n_g_i_threads = faces_list_size / _min_i_subset_size;
      if (! (faces_list_size % _min_i_subset_size))
        n_g_i_threads += 1;
    }

    /* Get an initial edge distribution */

    t_cell_index[0] = 0;
    for (t_id = 1; t_id < n_g_i_threads; t_id++) {
      t_cell_index[t_id] = t_cell_index[t_id-1] + group_size;
      if (t_cell_index[t_id] > n_f_cells)
        t_cell_index[t_id] = n_f_cells;
    }
    t_cell_index[n_g_i_threads] = n_f_cells;

    /* Pre-assign threads to faces (initial distribution) */

    _renum_face_multipass_assign(n_i_threads,
                                 n_g_i_threads,
                                 g_id,
                                 faces_list_size,
                                 faces_list,
                                 (const cs_lnum_2_t *restrict)l_face_cells,
                                 f_t_id,
                                 n_t_faces,
                                 t_face_last,
                                 t_cell_index);

    /* Try to redistribute the load */

    _renum_face_multipass_redistribute(n_i_threads,
                                       n_g_i_threads,
                                       g_id,
                                       redistribute_relaxation_factor,
                                       faces_list_size,
                                       faces_list,
                                       (const cs_lnum_2_t *restrict)l_face_cells,
                                       f_t_id,
                                       n_t_faces,
                                       t_face_last,
                                       t_cell_index);

    /* Update list of remaining faces */

    for (fl_id = 0; fl_id < faces_list_size; fl_id++) {

      f_id = faces_list[fl_id];

      if (f_t_id[f_id] < 0)
        faces_list[faces_list_size_new++] = f_id;

    }

    faces_list_size = faces_list_size_new;
    faces_list_size_new = 0;

    if (faces_list_size > 0)
      n_f_cells = _renum_face_multipass_remaining(n_f_cells,
                                                  faces_list_size,
                                                  faces_list,
                                                  l_face_cells);

  }

  /* Handle last group of faces */

  if (faces_list_size > 0) {

    for (fl_id = 0; fl_id < faces_list_size; fl_id++) {
      f_id = faces_list[fl_id];
      f_t_id[f_id] = g_id*n_i_threads;
    }

    g_id += 1;

    n_t_faces[0] = faces_list_size;
    for (t_id = 1; t_id < n_i_threads; t_id++)
      n_t_faces[t_id] = 0;

  }

  /* Free memory */

  BFT_FREE(l_face_cells);
  BFT_FREE(n_t_faces);
  BFT_FREE(t_face_last);
  BFT_FREE(t_cell_index);

  /* Now build final numbering and index */

  /* Build lexical ordering of faces */

  BFT_MALLOC(faces_keys, n_faces*3, cs_lnum_t);

# pragma omp parallel for private(c_id_0, c_id_1)
  for (f_id = 0; f_id < n_faces; f_id++) {
    faces_keys[f_id*3] = f_t_id[f_id];
    c_id_0 = i_face_cells[f_id][0];
    c_id_1 = i_face_cells[f_id][1];
    if (c_id_0 < c_id_1) {
      faces_keys[f_id*3 + 1] = c_id_0 - 1;
      faces_keys[f_id*3 + 2] = c_id_1 - 1;
    }
    else {
      faces_keys[f_id*3 + 1] = c_id_1 - 1;
      faces_keys[f_id*3 + 2] = c_id_0 - 1;
    }
  }

  cs_order_lnum_allocated_s(NULL,
                            faces_keys,
                            3,
                            faces_list,
                            n_faces);

  BFT_FREE(faces_keys);

  _n_groups = g_id;
  BFT_MALLOC(_group_index, _n_groups*n_i_threads*2, cs_lnum_t);

  _group_index[0] = 0;

  for (g_id=0; g_id < _n_groups; g_id++) {
    for (t_id = 0; t_id < n_i_threads; t_id++) {
      _group_index[(t_id*_n_groups + g_id)*2] = -1;
      _group_index[(t_id*_n_groups + g_id)*2 + 1] = -1;
    }
  }

  for (fl_id = 0; fl_id < n_faces; fl_id++) {

    f_id = faces_list[fl_id];
    renum_i[fl_id] = f_id + 1;

    assert(f_t_id[f_id] > -1);

    t_id = f_t_id[f_id]%n_i_threads;
    g_id = (f_t_id[f_id] - t_id) / n_i_threads;

    /* Update group index to mark maximum face id */
    _group_index[(t_id*_n_groups + g_id)*2 + 1] = fl_id + 1;

  }

  BFT_FREE(f_t_id);
  BFT_FREE(faces_list);

  /* Finalize group index */

  f_id = 0;
  for (g_id=0; g_id < _n_groups; g_id++) {
    for (t_id = 0; t_id < n_i_threads; t_id++) {
      _group_index[(t_id*_n_groups + g_id)*2] = f_id;
      f_id = CS_MAX(_group_index[(t_id*_n_groups + g_id)*2+1],
                    f_id);
    }
  }

  for (g_id=0; g_id < _n_groups; g_id++) {
    for (t_id = 0; t_id < n_i_threads; t_id++) {
      if (_group_index[(t_id*_n_groups + g_id)*2 + 1] < 0)
        _group_index[(t_id*_n_groups + g_id)*2] = -1;
    }
  }

  *n_groups = _n_groups;
  *group_index = _group_index;

  return 0;
}

/*----------------------------------------------------------------------------
 * Compute renumbering of faces using groups in which no two faces share
 * a cell.
 *
 * parameters:
 *   mesh           <-> pointer to global mesh structure
 *   n_i_threads    <-- number of threads required for interior faces
 *   max_group_size <-- target size for groups
 *   group_size     <-- target group size
 *   renum_i        --> interior faces renumbering array (new -> old, 1-based)
 *   n_i_groups     --> number of groups of interior faces
 *   i_group_index  --> group/thread index
 *
 * returns:
 *   0 on success, -1 otherwise
  *----------------------------------------------------------------------------*/

static int
_renum_i_faces_no_share_cell_in_block(cs_mesh_t    *mesh,
                                      int           n_i_threads,
                                      int           max_group_size,
                                      cs_lnum_t     renum_i[],
                                      cs_lnum_t    *n_i_groups,
                                      cs_lnum_t   **i_group_index)
{
  cs_lnum_t  *i_group_size = NULL;

  int retval = 0;

  while (   mesh->n_i_faces/max_group_size < 2*n_i_threads
         && max_group_size > _min_i_subset_size)
    max_group_size -= 64;

  if (max_group_size < _min_i_subset_size)
    max_group_size = _min_i_subset_size;
  if (max_group_size < n_i_threads*2)
    max_group_size = n_i_threads*2;

  _independent_face_groups(max_group_size,
                           mesh->n_cells_with_ghosts,
                           mesh->n_i_faces,
                           mesh->i_face_cells,
                           renum_i,
                           n_i_groups,
                           &i_group_size);

  BFT_MALLOC(*i_group_index, n_i_threads*(*n_i_groups)*2, cs_lnum_t);

  retval = _thread_bounds_by_group_size(mesh->n_i_faces,
                                        *n_i_groups,
                                        n_i_threads,
                                        i_group_size,
                                        *i_group_index);

  BFT_FREE(i_group_size);

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute renumbering of boundary faces for threads.
 *
 * As boundary faces belong to a single cell, boundary faces are
 * lexicographically ordered by their matching cell id, and subsets
 * of "almost" equal size are built, adjusted so that all boundary faces
 * sharing a cell are in the same subset.
 *
 * Usign this algorithm, a single group of subsets is required.
 *
 * parameters:
 *   mesh            <-> pointer to global mesh structure
 *   n_b_threads     <-- number of threads required for boundary faces
 *   min_subset_size <-- minimum size of subset associated to a thread
 *   renum_b         <-- interior faces renumbering array (new -> old)
 *   n_b_groups      --> number of groups of boundary faces
 *   b_group_index   --> group/thread index
 *
 * returns:
 *   0 on success, -1 otherwise
  *----------------------------------------------------------------------------*/

static int
_renum_b_faces_no_share_cell_across_thread(cs_mesh_t   *mesh,
                                           int          n_b_threads,
                                           cs_lnum_t    min_subset_size,
                                           cs_lnum_t    renum_b[],
                                           cs_lnum_t   *n_b_groups,
                                           cs_lnum_t  **b_group_index)
{
  int t_id;
  cs_lnum_t ii, subset_size, start_id, end_id;
  cs_lnum_t *order = NULL, *fc_num = NULL;

  int retval = 0;

  /* Initialization */

  *n_b_groups = 1;

  BFT_MALLOC(*b_group_index, n_b_threads*2, cs_lnum_t);

  /* Order faces lexicographically */

  BFT_MALLOC(order, mesh->n_b_faces, cs_lnum_t);
  BFT_MALLOC(fc_num, mesh->n_b_faces*2, cs_lnum_t);

  for (ii = 0; ii < mesh->n_b_faces; ii++) {
    fc_num[ii*2] = mesh->b_face_cells[ii];
    fc_num[ii*2+1] = ii;
  }

  cs_order_lnum_allocated_s(NULL, fc_num, 2, order, mesh->n_b_faces);

  BFT_FREE(fc_num);

  /* Build new numbering index */

  for (ii = 0; ii < mesh->n_b_faces; ii++)
    renum_b[ii] = order[ii] + 1;

  BFT_FREE(order);

  /* Compute target subset size */

  subset_size = mesh->n_b_faces / n_b_threads;
  if (mesh->n_b_faces % n_b_threads > 0)
    subset_size++;
  subset_size = CS_MAX(subset_size, min_subset_size);

  /* Build then adjust group / thread index */

  for (t_id = 0, end_id = 0; t_id < n_b_threads; t_id++) {

    start_id = end_id;
    end_id = (t_id+1)*subset_size;

    if (end_id < start_id)
      end_id = start_id;

    if (end_id > mesh->n_b_faces)
      end_id = mesh->n_b_faces;
    else if (end_id > 0 && end_id < mesh->n_b_faces) {
      cs_lnum_t f_id = renum_b[end_id - 1] - 1;
      cs_lnum_t c_num = mesh->b_face_cells[f_id];
      f_id = renum_b[end_id] - 1;
      while (mesh->b_face_cells[f_id] == c_num) {
        end_id += 1;
        if (end_id < mesh->n_b_faces)
          f_id = renum_b[end_id] - 1;
        else
          break;
      }
    }

    (*b_group_index)[t_id*2] = start_id;
    (*b_group_index)[t_id*2+1] = end_id;
  }

  if (mesh->n_b_faces < 1)
    retval = -1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Log statistics for bandwidth and profile.
 *
 * Bandwidth ist the maximum distance between two adjacent vertices (cells),
 * with distance measured by the difference of vertex (cell) ids.
 *
 * Profile is the sum of all the maximum distances between the i-th vertex
 * and any of its neighbors with an index j > i (as the matrix structure
 * is symmetric, this simplifies to the sum of the maximum distances between
 * a vertex and any of its neighbors).
 *
 * parameters:
 *   mesh      <-- associated mesh
 *   title     <-- title or name of mesh or matrix
 *----------------------------------------------------------------------------*/

static void
_log_bandwidth_info(const cs_mesh_t  *mesh,
                    const char       *title)
{
  cs_lnum_t cell_id, face_id;

  cs_lnum_t bandwidth = 0;
  cs_gnum_t profile = 0;
  cs_lnum_t *max_distance = NULL;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;

  BFT_MALLOC(max_distance, mesh->n_cells_with_ghosts, cs_lnum_t);

  for (cell_id = 0; cell_id < mesh->n_cells_with_ghosts; cell_id++)
    max_distance[cell_id] = 0;

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {

    cs_lnum_t cid0 = i_face_cells[face_id][0] - 1;
    cs_lnum_t cid1 = i_face_cells[face_id][1] - 1;

    cs_lnum_t distance = CS_ABS(cid1 - cid0);

    if (distance > bandwidth)
      bandwidth = distance;

    if (distance > max_distance[cid0])
      max_distance[cid0] = distance;

    if (distance > max_distance[cid1])
      max_distance[cid1] = distance;
  }

  for (cell_id = 0; cell_id < mesh->n_cells; cell_id++)
    profile += max_distance[cell_id];

  profile /= mesh->n_cells;

  BFT_FREE(max_distance);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t loc_buffer;
    cs_gnum_t *rank_buffer = NULL;
    BFT_MALLOC(rank_buffer, cs_glob_n_ranks, cs_gnum_t);

    loc_buffer = bandwidth;
    MPI_Allgather(&loc_buffer, 1, CS_MPI_GNUM,
                  rank_buffer, 1, CS_MPI_GNUM, cs_glob_mpi_comm);
    bft_printf
      (_("\n Histogram of %s matrix bandwidth per rank:\n\n"), title);
    _display_histograms_gnum(cs_glob_n_ranks, rank_buffer);

    loc_buffer = profile;
    MPI_Allgather(&loc_buffer, 1, CS_MPI_GNUM,
                  rank_buffer, 1, CS_MPI_GNUM, cs_glob_mpi_comm);

    bft_printf
      (_("\n Histogram of %s matrix profile/lines per rank:\n\n"), title);
    _display_histograms_gnum(cs_glob_n_ranks, rank_buffer);

    BFT_FREE(rank_buffer);

  } /* End if cs_glob_n_ranks > 1 */

#endif

  if (cs_glob_n_ranks == 1) {
    bft_printf
      (_("\n Matrix bandwidth for %s :          %llu\n"
         " Matrix profile/lines for %s :      %llu\n"),
       title, (unsigned long long)bandwidth,
       title, (unsigned long long)profile);
  }
}

/*----------------------------------------------------------------------------
 * Estimate unbalance between threads of a given group.
 *
 * Test local operations related to renumbering.
 *
 * Unbalance is considered to be: (max/mean - 1)
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 * returns:
 *   estimated unbalance for this group
 *----------------------------------------------------------------------------*/

static double
_estimate_imbalance(const cs_numbering_t  *face_numbering)
{
  double t_imbalance_tot = 0.0;

  if (face_numbering == NULL)
    return 0;

  if (face_numbering->type == CS_NUMBERING_THREADS) {

    int g_id;

    cs_lnum_t n_faces = 0;

    const int n_threads = face_numbering->n_threads;
    const int n_groups = face_numbering->n_groups;
    const cs_lnum_t *group_index = face_numbering->group_index;

    for (g_id = 0; g_id < n_groups; g_id++) {

      int t_id;
      double n_t_faces_mean, imbalance;

      cs_lnum_t n_t_faces_sum = 0;
      cs_lnum_t n_t_faces_max = 0;

      for (t_id = 0; t_id < n_threads; t_id++) {
        cs_lnum_t n_t_faces =   group_index[(t_id*n_groups + g_id)*2 + 1]
                              - group_index[(t_id*n_groups + g_id)*2];
        n_t_faces = CS_MAX(n_t_faces, 0);
        n_t_faces_sum += n_t_faces;
        n_t_faces_max = CS_MAX(n_t_faces, n_t_faces_max);
      }

      n_faces += n_t_faces_sum;

      n_t_faces_mean = (double)n_t_faces_sum / n_threads;

      imbalance = (n_t_faces_max / n_t_faces_mean) - 1.0;
      t_imbalance_tot += imbalance*n_t_faces_sum;

    }

    t_imbalance_tot /= n_faces;

  }

  return t_imbalance_tot;
}

/*----------------------------------------------------------------------------
 * Log statistics for threads and groups.
 *
 * parameters:
 *   elt_type_name <-- name of element type (interior of boundary face)
 *   n_domains     <-- number of MPI domains
 *   n_threads     <-- local number of threads
 *   n_groups      <-- local number of groups
 *   imbalance     <-- estimation of imbalance
 *----------------------------------------------------------------------------*/

static void
_log_threading_info(const char  *elt_type_name,
                    int          n_domains,
                    int          n_threads,
                    int          n_groups,
                    double       imbalance)
{
  /* Build histograms for number of threads, number for groups,
     and group size */

#if defined(HAVE_MPI)

  if (n_domains > 1) {

    cs_gnum_t loc_buffer;
    double d_loc_buffer;
    cs_gnum_t *rank_buffer = NULL;
    double *d_rank_buffer = NULL;

    BFT_MALLOC(rank_buffer, n_domains, cs_gnum_t);

    loc_buffer = n_threads;
    MPI_Allgather(&loc_buffer, 1, CS_MPI_GNUM,
                  rank_buffer, 1, CS_MPI_GNUM, cs_glob_mpi_comm);
    bft_printf
      (_("\n Histogram of thread pools size for %s per rank:\n\n"),
       elt_type_name);
    _display_histograms_gnum(n_domains, rank_buffer);

    loc_buffer = n_groups;
    MPI_Allgather(&loc_buffer, 1, CS_MPI_GNUM,
                  rank_buffer, 1, CS_MPI_GNUM, cs_glob_mpi_comm);
    bft_printf
      (_("\n Histogram of threading groups count for %s per rank:\n\n"),
       elt_type_name);
    _display_histograms_gnum(n_domains, rank_buffer);

    BFT_FREE(rank_buffer);

    BFT_MALLOC(d_rank_buffer, n_domains, double);

    d_loc_buffer = imbalance;
    MPI_Allgather(&d_loc_buffer, 1, MPI_DOUBLE,
                  d_rank_buffer, 1, MPI_DOUBLE, cs_glob_mpi_comm);
    bft_printf
      (_("\n Histogram of thread imbalance for %s per rank:\n\n"),
       elt_type_name);
    _display_histograms_double(n_domains, d_rank_buffer);

    BFT_FREE(rank_buffer);

  } /* End if n_domains > 1 */

#endif

  if (n_domains == 1) {
    bft_printf
      (_("\n Number of thread pools for %s :          %d\n"
         " Number of threading groups for %s :      %d\n"
         " Estimated thread imbalance for %s :      %10.5e\n"),
       elt_type_name, n_threads,
       elt_type_name, n_groups,
       elt_type_name, imbalance);
  }
}

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces and cells for multiple threads.
 *
 * Relation to graph edge coloring:
 * No graph vertex (cell) is incident to 2 edges (faces) of the same color.
 * A thread pool may thus be built, with 1 thread per color.
 * Groups may then be built, containing only cells of a given color.
 *
 * parameters:
 *   mesh  <->  Pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_for_threads(cs_mesh_t  *mesh)
{
  int  update_c = 0, update_fi = 0, update_fb = 0;
  int  n_i_groups = 1, n_b_groups = 1;
  cs_lnum_t  max_group_size = 1014;       /* Default */
  cs_lnum_t  ii;
  cs_lnum_t  *renum_c = NULL, *renum_i = NULL, *renum_b = NULL;
  cs_lnum_t  *i_group_index = NULL, *b_group_index = NULL;

  int  n_i_threads = _cs_renumber_n_threads;
  int  n_b_threads = _cs_renumber_n_threads;

  int retval = 0;

  /* Note: group indexes for n_threads and n_groups are defined as follows:
   *  group_index <-- group_index[thread_id*group_id*2 + 2*group_id] and
   *                  group_index[thread_id*group_id*2 + 2*group_id +1]
   *                  define the tart and end ids (+1) for entities in a
   *                  given group and thread (size: n_groups *2 * n_threads) */

  if (_cs_renumber_n_threads < 2)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(renum_c, mesh->n_cells_with_ghosts, cs_lnum_t);
  BFT_MALLOC(renum_i, mesh->n_i_faces, cs_lnum_t);
  BFT_MALLOC(renum_b, mesh->n_b_faces, cs_lnum_t);

  /* Initialize renumbering arrays */

  {
    for (ii = 0; ii < mesh->n_cells_with_ghosts; ii++)
      renum_c[ii] = ii+1;

    for (ii = 0; ii < mesh->n_i_faces; ii++)
      renum_i[ii] = ii+1;

    for (ii = 0; ii < mesh->n_b_faces; ii++)
      renum_b[ii] = ii+1;
  }

  /* Interior faces renumbering */
  /*----------------------------*/

  /* Adjust block size depending on the number of faces and threads */

  switch (_i_faces_algorithm) {
  case CS_RENUMBER_I_FACES_BLOCK:
    retval = _renum_i_faces_no_share_cell_in_block(mesh,
                                                   n_i_threads,
                                                   max_group_size,
                                                   renum_i,
                                                   &n_i_groups,
                                                   &i_group_index);
    break;

  case CS_RENUMBER_I_FACES_MULTIPASS:
    retval = _renum_face_multipass(mesh,
                                 n_i_threads,
                                 renum_i,
                                 &n_i_groups,
                                 &i_group_index);
    break;

  case CS_RENUMBER_I_FACES_NONE:
  default:
    retval = -1;
    break;
  }

  if (retval != 0) {
    n_i_groups = 1;
    n_i_threads = 1;
    update_fi = 0;
    BFT_FREE(i_group_index);
  }
  else
    update_fi = 1;

  /* Transfer interior face numbering information to mesh */

  if (n_i_groups *n_i_threads > 1)
    mesh->i_face_numbering = cs_numbering_create_threaded(n_i_threads,
                                                          n_i_groups,
                                                          i_group_index);
  BFT_FREE(i_group_index);

  _log_threading_info(_("interior faces"),
                      mesh->n_domains,
                      n_i_threads,
                      n_i_groups,
                      _estimate_imbalance(mesh->i_face_numbering));

  /* Boundary faces renumbering */
  /*----------------------------*/

  retval = _renum_b_faces_no_share_cell_across_thread(mesh,
                                                      n_b_threads,
                                                      _min_b_subset_size,
                                                      renum_b,
                                                      &n_b_groups,
                                                      &b_group_index);

  if (retval != 0) {
    n_b_groups = 1;
    n_b_threads = 1;
    update_fb = 0;
    BFT_FREE(b_group_index);
  }
  else
    update_fb = 1;

  /* Transfer boundary face numbering information to mesh */

  if (n_b_groups *n_b_threads > 1)
    mesh->b_face_numbering = cs_numbering_create_threaded(n_b_threads,
                                                          n_b_groups,
                                                          b_group_index);
  BFT_FREE(b_group_index);

  _log_threading_info(_("boundary faces"),
                      mesh->n_domains,
                      n_b_threads,
                      n_b_groups,
                      _estimate_imbalance(mesh->b_face_numbering));

  bft_printf("\n ----------------------------------------------------------\n");

  /* Free memory */

  if (update_c == 0)
    BFT_FREE(renum_c);

  if (update_fi == 0)
    BFT_FREE(renum_i);

  if (update_fb == 0)
    BFT_FREE(renum_b);

  /* Now update mesh connectivity */
  /*------------------------------*/

  if (renum_i != NULL || renum_b != NULL)
    _cs_renumber_update_faces(mesh,
                              renum_i,
                              renum_b);

  if (renum_c != NULL)
    _cs_renumber_update_cells(mesh,
                              renum_c);

  /* Now free remaining arrays */

  BFT_FREE(renum_i);
  BFT_FREE(renum_b);
  BFT_FREE(renum_c);

}

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces for vector machines.
 *
 * Renumbering can be cancelled using the IVECTI and IVECTB values in
 * Fortan common IVECTO: -1 indicates we should try to renumber,
 * 0 means we should not renumber. On exit, 0 means we have not found an
 * adequate renumbering, 1 means we have (and it was applied).
 *
 * If the target architecture does not enable vectorization, do as if no
 * adequate renumbering was found.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *
 * returns:
 *   1 if renumbering was tried, 0 otherwise.
 *----------------------------------------------------------------------------*/

static int
_renumber_for_vectorizing(cs_mesh_t  *mesh)
{
  int _ivect[2] = {0, 0};
  cs_int_t   ivecti = 0, ivectb = 0;
  cs_int_t  *renum_i = NULL, *renum_b = NULL;
  cs_int_t  *iworkf = NULL, *ismbs = NULL;

  cs_int_t  n_faces_max = CS_MAX(mesh->n_i_faces, mesh->n_b_faces);
  cs_int_t  n_cells_wghosts = mesh->n_cells_with_ghosts;


#if defined(__uxpvp__) /* For Fujitsu VPP5000 (or possibly successors) */

  /* Vector register numbers and lengths:
   *   4       4096 ;
   *  16       1024
   *  32        512
   *  64        256
   * 128        128
   * 256         64 */

  const int vector_size = 1024; /* Use register 16 */

#elif defined(SX) && defined(_SX) /* For NEC SX series */

  const int vector_size = 256; /* At least for NEC SX-9 */

#else

  const int vector_size = 1; /* Non-vector machines */

#endif

  /* Nothing to do if vector size = 1 */

  if (vector_size == 1)
    return 0;

  /* Allocate Work arrays */

  BFT_MALLOC(renum_i, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(renum_b, mesh->n_b_faces, cs_int_t);
  BFT_MALLOC(iworkf, n_faces_max, cs_int_t);
  BFT_MALLOC(ismbs, n_cells_wghosts, cs_int_t);

  /* Try renumbering */

  CS_PROCF(numvec, NUMVEC)(&(mesh->n_cells_with_ghosts),
                           &(mesh->n_cells),
                           &(mesh->n_i_faces),
                           &(mesh->n_b_faces),
                           &vector_size,
                           &ivecti,
                           &ivectb,
                           mesh->i_face_cells,
                           mesh->b_face_cells,
                           renum_i,
                           renum_b,
                           iworkf,
                           ismbs);

  /* Free Work arrays */

  BFT_FREE(ismbs);
  BFT_FREE(iworkf);

  /* Update mesh */

  if (ivecti > 0 || ivectb > 0) {

    cs_int_t   *_renum_i = NULL;
    cs_int_t   *_renum_b = NULL;

    if (ivecti > 0)
      _renum_i = renum_i;
    if (ivectb > 0)
      _renum_b = renum_b;

    _cs_renumber_update_faces(mesh,
                              _renum_i,
                              _renum_b);

  }

  /* Free final work arrays */

  BFT_FREE(renum_b);
  BFT_FREE(renum_i);

  /* Check renumbering (sanity check) */

  if (ivecti > 0 || ivectb > 0) {

    cs_int_t  *ismbv = NULL;
    cs_real_t  *rworkf = NULL, *rsmbs = NULL, *rsmbv = NULL;

    BFT_MALLOC(iworkf, n_faces_max, cs_int_t);
    BFT_MALLOC(ismbs, n_cells_wghosts, cs_int_t);
    BFT_MALLOC(ismbv, n_cells_wghosts, cs_int_t);
    BFT_MALLOC(rworkf, n_faces_max, cs_real_t);
    BFT_MALLOC(rsmbs, n_cells_wghosts, cs_real_t);
    BFT_MALLOC(rsmbv, n_cells_wghosts, cs_real_t);

    CS_PROCF(tstvec, TSTVEC)(&(mesh->n_cells_with_ghosts),
                             &(mesh->n_cells),
                             &(mesh->n_i_faces),
                             &(mesh->n_b_faces),
                             mesh->i_face_cells,
                             mesh->b_face_cells,
                             iworkf,
                             ismbs,
                             ismbv,
                             rworkf,
                             rsmbs,
                             rsmbv);

    BFT_FREE(rsmbv);
    BFT_FREE(rsmbs);
    BFT_FREE(rworkf);
    BFT_FREE(ismbv);
    BFT_FREE(ismbs);
    BFT_FREE(iworkf);
  }

  /* Update mesh */

  if (ivecti > 0)
    mesh->i_face_numbering
      = cs_numbering_create_vectorized(mesh->n_i_faces, vector_size);
  if (ivectb > 0)
    mesh->b_face_numbering
      = cs_numbering_create_vectorized(mesh->n_b_faces, vector_size);

  /* Output info */

  _ivect[0] = ivecti; _ivect[1] = ivectb;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    int ivect_tot[2];
    MPI_Allreduce(_ivect, ivect_tot, 2, MPI_INT, MPI_SUM,
                  cs_glob_mpi_comm);
    _ivect[0] = ivect_tot[0]; _ivect[1] = ivect_tot[1];
  }
#endif

  bft_printf(_("\n"
               " Vectorization:\n"
               " --------------\n"
               "   interior faces: %d ranks (of %d)\n"
               "   boundary faces: %d ranks\n\n"),
             _ivect[0], cs_glob_n_ranks, _ivect[1]);

  return 1;
}

/*----------------------------------------------------------------------------
 * Test local operations related to renumbering.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_test(cs_mesh_t  *mesh)
{
  cs_gnum_t face_errors[2] = {0, 0};
  cs_lnum_t *accumulator = NULL;

  if (mesh == NULL)
    return;

  if (mesh->verbosity > 0)
    bft_printf
      (_("\n"
         "Checking mesh renumbering for threads:\n"
         "-------------------------------------\n\n"));

  /* Check for interior faces */
  /*--------------------------*/

  if (mesh->i_face_numbering != NULL) {

    if (mesh->i_face_numbering->type == CS_NUMBERING_THREADS) {

      int g_id, t_id;
      cs_lnum_t f_id, c_id_0, c_id_1;

      cs_lnum_t counter = 0;

      const int n_threads = mesh->i_face_numbering->n_threads;
      const int n_groups = mesh->i_face_numbering->n_groups;
      const cs_lnum_t *group_index = mesh->i_face_numbering->group_index;

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        accumulator[c_id_0] = 0;

      for (g_id=0; g_id < n_groups; g_id++) {

#       pragma omp parallel for private(f_id, c_id_0, c_id_1)
        for (t_id=0; t_id < n_threads; t_id++) {
          for (f_id = group_index[(t_id*n_groups + g_id)*2];
               f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               f_id++) {
            c_id_0 = mesh->i_face_cells[f_id*2] - 1;
            c_id_1 = mesh->i_face_cells[f_id*2 + 1] - 1;
            accumulator[c_id_0] += 1;
            accumulator[c_id_1] += 1;
          }
        }

      }

      for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        counter += accumulator[c_id_0];

      face_errors[0] = mesh->n_i_faces*2 - counter;

      /* Additional serial test */

      if (face_errors[0] == 0) {

        for (g_id=0; g_id < n_groups; g_id++) {

          for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
            accumulator[c_id_0] = -1;

          for (t_id=0; t_id < n_threads; t_id++) {
            for (f_id = group_index[(t_id*n_groups + g_id)*2];
                 f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
                 f_id++) {
              c_id_0 = mesh->i_face_cells[f_id*2] - 1;
              c_id_1 = mesh->i_face_cells[f_id*2 + 1] - 1;
              if (   (accumulator[c_id_0] > -1 && accumulator[c_id_0] != t_id)
                  || (accumulator[c_id_1] > -1 && accumulator[c_id_1] != t_id))
                face_errors[0] += 1;
              if ( (accumulator[c_id_0] > -1 && accumulator[c_id_0] != t_id)
                || (accumulator[c_id_1] > -1 && accumulator[c_id_1] != t_id)) {
                if (mesh->verbosity > 0)
                  bft_printf("f_id %d (%d %d) g %d t %d\n",
                             f_id, c_id_0, c_id_1, g_id, t_id);
                accumulator[c_id_0] = t_id;
                accumulator[c_id_1] = t_id;
              }
            }
          }

        }

      }

      BFT_FREE(accumulator);
    }

  }

  /* Check for boundary faces */
  /*--------------------------*/

  if (mesh->b_face_numbering != NULL) {

    if (mesh->b_face_numbering->type == CS_NUMBERING_THREADS) {

      int g_id, t_id;
      cs_lnum_t f_id, c_id;

      cs_lnum_t counter = 0;

      const int n_threads = mesh->b_face_numbering->n_threads;
      const int n_groups = mesh->b_face_numbering->n_groups;
      const cs_lnum_t *group_index = mesh->b_face_numbering->group_index;

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
        accumulator[c_id] = 0;

      for (g_id=0; g_id < n_groups; g_id++) {

#       pragma omp parallel for private(f_id, c_id)
        for (t_id=0; t_id < n_threads; t_id++) {
          for (f_id = group_index[(t_id*n_groups + g_id)*2];
               f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               f_id++) {
            c_id = mesh->b_face_cells[f_id] - 1;
            accumulator[c_id] += 1;
          }
        }

      }

      for (c_id = 0; c_id < mesh->n_cells; c_id++)
        counter += accumulator[c_id];

      face_errors[1] = mesh->n_b_faces - counter;

      /* Additional serial test */

      if (face_errors[1] == 0) {

        for (g_id=0; g_id < n_groups; g_id++) {

          for (c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
            accumulator[c_id] = -1;

          for (t_id=0; t_id < n_threads; t_id++) {
            for (f_id = group_index[(t_id*n_groups + g_id)*2];
                 f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
                 f_id++) {
              c_id = mesh->b_face_cells[f_id] - 1;
              if (accumulator[c_id] > -1 && accumulator[c_id] != t_id)
                face_errors[1] += 1;
              accumulator[c_id] = t_id;
            }
          }

        }

      }

      BFT_FREE(accumulator);
    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t  g_face_errors[2];
    MPI_Allreduce(face_errors, g_face_errors, 2, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    face_errors[0] = g_face_errors[0];
    face_errors[1] = g_face_errors[1];
  }
#endif

  if (face_errors[0] != 0 || face_errors[1] != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Conflicts detected using mesh renumbering:\n"
                "  for interior faces: %llu\n"
                "  for boundary faces: %llu"),
              (unsigned long long)(face_errors[0]),
              (unsigned long long)(face_errors[1]));
}

/*----------------------------------------------------------------------------
 * Renumber mesh elements for vectorization or OpenMP depending on code
 * options and target machine.
 *
 * Currently, only the legacy vectorizing renumbering is handled.
 *
 * parameters:
 *   mesh  <->  Pointer to global mesh structure
 *
 *----------------------------------------------------------------------------*/

static void
_renumber_mesh(cs_mesh_t  *mesh)
{
  int retval = 0;
  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    _cs_renumber_n_threads = cs_glob_n_threads;

  p = getenv("CS_RENUMBER");

  if (p != NULL) {

    if (strcmp(p, "off") == 0) {
      bft_printf(_("\n Mesh renumbering off.\n\n"));
      return;
    }

#if defined(HAVE_IBM_RENUMBERING_LIB)
    if (strcmp(p, "IBM") == 0) {
      bft_printf("\n Use IBM Mesh renumbering.\n\n");
      _renumber_for_threads_ibm(mesh);
      _renumber_test(mesh);
      return;
    }
#endif

  }

  /* Try vectorizing first, then renumber for Cache / OpenMP */

  retval = _renumber_for_vectorizing(mesh);

  if (retval == 0)
    _renumber_for_threads(mesh);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the target number of threads for mesh renumbering.
 *
 * By default, the target number of threads is set to cs_glob_n_threads,
 * but the value may be forced using this function. This is mainly useful
 * for testing purposes.
 *
 * \param[in]  n_threads  target number of threads for mesh numbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_set_n_threads(int  n_threads)
{
  _cs_renumber_n_threads = n_threads;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the target number of threads for mesh renumbering.
 *
 * \return  the target number of threads for mesh numbering
 */
/*----------------------------------------------------------------------------*/

int
cs_renumber_get_n_threads(void)
{
  return _cs_renumber_n_threads;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the minimum sunset sizes when renumbering for threads.
 *
 * \param[in]  min_i_subset_size  minimum number of interior faces per
 *                                thread per group
 * \param[in]  min_b_subset_size  minimum number of boundary faces per
 *                                thread per group
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_set_min_subset_size(cs_lnum_t  min_i_subset_size,
                                cs_lnum_t  min_b_subset_size)
{
  _min_i_subset_size = min_i_subset_size;
  _min_b_subset_size = min_b_subset_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the minimum sunset sizes when renumbering for threads.
 *
 * \param[out]  min_i_subset_size  minimum number of interior faces per
 *                                 thread per group, or NULL
 * \param[out]  min_b_subset_size  minimum number of boundary faces per
 *                                 thread per group, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_get_min_subset_size(cs_lnum_t  *min_i_subset_size,
                                cs_lnum_t  *min_b_subset_size)
{
  if (min_i_subset_size != NULL)
    *min_i_subset_size = _min_i_subset_size;
  if (min_b_subset_size != NULL)
    *min_b_subset_size = _min_b_subset_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select the algorithm for interior faces renumbering.
 *
 * \param[in]  algorithm  algorithm type for interior faces renumbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_set_i_face_algorithm(cs_renumber_i_faces_type_t  algorithm)
{
  _i_faces_algorithm = algorithm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the algorithm for interior faces renumbering.
 *
 * \return  algorithm type for interior faces renumbering
 */
/*----------------------------------------------------------------------------*/

cs_renumber_i_faces_type_t
cs_renumber_get_i_face_algorithm(void)
{
  return _i_faces_algorithm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber mesh elements for vectorization or OpenMP depending on code
 * options and target machine.
 *
 * \param[in, out]  mesh             Pointer to global mesh structure
 * \param[in, out]  mesh_quantities  Pointer to global mesh quantities
 *                                   structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t             *mesh,
                 cs_mesh_quantities_t  *mesh_quantities)
{
  bool quantities_computed = false;

  if (mesh_quantities != NULL) {
    if (mesh_quantities->cell_cen != NULL)
      quantities_computed = true;
  }

  _renumber_mesh(mesh);

  if (mesh->i_face_numbering == NULL)
    mesh->i_face_numbering = cs_numbering_create_default(mesh->n_i_faces);
  if (mesh->b_face_numbering == NULL)
    mesh->b_face_numbering = cs_numbering_create_default(mesh->n_b_faces);

  _renumber_test(mesh);

  if (mesh->verbosity > 0)
    _log_bandwidth_info(mesh, _("volume mesh"));

  if (quantities_computed)
    cs_mesh_quantities_compute(mesh, mesh_quantities);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
