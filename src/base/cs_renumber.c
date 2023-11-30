/*============================================================================
 * Optional mesh renumbering
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

#include <float.h>
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
 * METIS library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HAVE_PARMETIS)
#include <parmetis.h>
#endif

#if defined(HAVE_METIS)
#include <metis.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

/*----------------------------------------------------------------------------
 * SCOTCH library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_SCOTCH)
#include <stdint.h>
#include <scotch.h>
#elif defined(HAVE_PTSCOTCH)
#include <stdint.h>
#include <ptscotch.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_morton.h"
#include "fvm_hilbert.h"

#include "cs_array.h"
#include "cs_defs.h"
#include "cs_halo.h"
#include "cs_join.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sort.h"
#include "cs_prototypes.h"

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

  \enum cs_renumber_cells_type_t

  \brief Cell renumbering algorithm types

  \var CS_RENUMBER_CELLS_SCOTCH_PART
       Subpartition for thread blocks based SCOTCH library.
  \var CS_RENUMBER_CELLS_SCOTCH_ORDER
       Fill-reducing ordering based on SCOTCH library.
  \var CS_RENUMBER_CELLS_METIS_PART
       Subpartition for thread blocks based METIS library.
  \var CS_RENUMBER_CELLS_METIS_ORDER
       Fill-reducing ordering based on METIS library.
  \var CS_RENUMBER_CELLS_MORTON
       Order cells using domain-local Morton space-filling curve.
  \var CS_RENUMBER_CELLS_HILBERT
       Order cells using domain-local Hilbert space-filling curve.
  \var CS_RENUMBER_CELLS_RCM
       Order cells using domain-local reverse Cuthill-McKee algorithm.
  \var CS_RENUMBER_CELLS_NONE
       No cells renumbering.

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
  \var CS_RENUMBER_I_FACES_SIMD
       Renumber to allow SIMD operations in interior face->cell gather
       operations (such as SpMV products with native matrix representation).
  \var CS_RENUMBER_I_FACES_NONE
       No interior face renumbering.

  \enum cs_renumber_b_faces_type_t

  \brief Boundary faces renumbering algorithm types

  \var CS_RENUMBER_B_FACES_THREAD
       Renumber for threads, with one block per thread, and no cell
       referenced by faces in different threads blocks.
  \var CS_RENUMBER_B_FACES_SIMD
       Renumber to allow SIMD operations in boundary face->cell gather
       operations.
  \var CS_RENUMBER_B_FACES_NONE
       No interior face renumbering.

  \enum cs_renumber_vertices_type_t

  \brief Vertices renumbering algorithm types

  \var CS_RENUMBER_VERTICES_BY_CELL_ADJ
       Renumbering based on the cell adjacency.
  \var CS_RENUMBER_VERTICES_BY_FACE_ADJ
       Renumbering based on the face adjacency.
  \var CS_RENUMBER_VERTICES_NONE
       No vertex renumbering.

  \enum cs_renumber_ordering_t

  \brief Ordering options for adjacency arrays

  \var CS_RENUMBER_ADJACENT_LOW
       Lexicographical ordering with lowest adjacent id first

  \var CS_RENUMBER_ADJACENT_HIGH
       Lexicographical ordering with highest adjacent id first
} ;

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_RENUMBER_N_SUBS  5  /* Number of categories for histograms */

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static int _cs_renumber_vector_size = CS_NUMBERING_SIMD_SIZE;

static int _cs_renumber_n_threads = 0;

static cs_lnum_t  _min_i_subset_size = 256;
static cs_lnum_t  _min_b_subset_size = 256;

static bool _renumber_ghost_cells = true;
static bool _cells_adjacent_to_halo_last = false;
static bool _i_faces_adjacent_to_halo_last = false;
static cs_renumber_ordering_t _i_faces_base_ordering = CS_RENUMBER_ADJACENT_LOW;

static cs_renumber_cells_type_t _cells_algorithm[] = {CS_RENUMBER_CELLS_NONE,
                                                      CS_RENUMBER_CELLS_NONE};

#if defined(HAVE_OPENMP)
static cs_renumber_i_faces_type_t _i_faces_algorithm
  = CS_RENUMBER_I_FACES_MULTIPASS;
static cs_renumber_b_faces_type_t _b_faces_algorithm
  = CS_RENUMBER_B_FACES_THREAD;
#elif defined(__NEC__)
static cs_renumber_i_faces_type_t _i_faces_algorithm
  = CS_RENUMBER_I_FACES_SIMD;
static cs_renumber_b_faces_type_t _b_faces_algorithm
   = CS_RENUMBER_B_FACES_SIMD;
#else
static cs_renumber_i_faces_type_t _i_faces_algorithm
  = CS_RENUMBER_I_FACES_NONE;
static cs_renumber_b_faces_type_t _b_faces_algorithm
   = CS_RENUMBER_B_FACES_NONE;
#endif

static cs_renumber_vertices_type_t _vertices_algorithm
  = CS_RENUMBER_VERTICES_NONE;

static const char *_cell_renum_name[]
  = {N_("sub-partitioning with LibScotch"),
     N_("fill-reducing ordering with LibScotch"),
     N_("sub-partitioning with METIS"),
     N_("fill-reducing ordering with METIS"),
     N_("Morton curve in local bounding box"),
     N_("Hilbert curve in local bounding box"),
     N_("Reverse Cuthill-McKee"),
     N_("no renumbering")};

static const char *_i_face_renum_name[]
  = {N_("coloring, no shared cell in block"),
     N_("multipass"),
     N_("vectorizing"),
     N_("adjacent cells")};

static const char *_b_face_renum_name[]
  = {N_("no shared cell across threads"),
     N_("vectorizing"),
     N_("adjacent cells")};

static const char *_vertices_renum_name[]
= {N_("adjacent cells"),
   N_("adjacent faces"),
   N_("no renumbering")};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Redistribute family (group class) ids in case of renubering
 *
 * parameters:
 *   n_elts     <--  Number of elements
 *   new_to_old <--  Pointer to renumbering array
 *   family     <->  Pointer to array of family ids (or NULL)
 *----------------------------------------------------------------------------*/

static void
_update_family(cs_lnum_t         n_elts,
               const cs_lnum_t  *new_to_old,
               int              *family)
{
  int *old_family;

  if (family == NULL)
    return;

  BFT_MALLOC(old_family, n_elts, int);

  memcpy(old_family, family, n_elts*sizeof(int));

  for (cs_lnum_t ii = 0; ii < n_elts; ii++)
    family[ii] = old_family[new_to_old[ii]];

  BFT_FREE(old_family);
}

/*----------------------------------------------------------------------------
 * Redistribute face generation level in case of renubering
 *
 * parameters:
 *   n_elts     <--  Number of elements
 *   new_to_old <--  Pointer to renumbering array
 *   family     <->  Pointer to array of family ids (or NULL)
 *----------------------------------------------------------------------------*/

static void
_update_r_gen(cs_lnum_t         n_elts,
              const cs_lnum_t  *new_to_old,
              char             *r_gen)
{
  char *old_r_gen;

  if (r_gen == NULL)
    return;

  BFT_MALLOC(old_r_gen, n_elts, char);

  memcpy(old_r_gen, r_gen, n_elts);

  for (cs_lnum_t ii = 0; ii < n_elts; ii++)
    r_gen[ii] = old_r_gen[new_to_old[ii]];

  BFT_FREE(old_r_gen);
}

/*----------------------------------------------------------------------------
 * Update a global numbering array in case of entity renumbering
 *
 * parameters:
 *   n_elts      --> number of elements in array
 *   new_to_old  --> renumbering array
 *   global_num  <-> global numbering (allocated if initially NULL)
 *----------------------------------------------------------------------------*/

static void
_update_global_num(size_t             n_elts,
                   const cs_lnum_t    new_to_old[],
                   cs_gnum_t        **global_num)
{
  cs_gnum_t *_global_num = *global_num;

  if (_global_num == NULL) {

    BFT_MALLOC(_global_num, n_elts, cs_gnum_t);

    for (size_t i = 0; i < n_elts; i++)
      _global_num[i] = new_to_old[i] + 1;

    *global_num = _global_num;
  }

  else {

    cs_gnum_t *tmp_global;

    BFT_MALLOC(tmp_global, n_elts, cs_gnum_t);
    memcpy(tmp_global, _global_num, n_elts*sizeof(cs_gnum_t));

    for (size_t i = 0; i < n_elts; i++)
      _global_num[i] = tmp_global[new_to_old[i]];

    BFT_FREE(tmp_global);
  }
}

/*----------------------------------------------------------------------------
 * Reorder ghost cells according to their position in the halo and adjacency
 *
 * Reordering keeps cells in a section based halo class (separated based on
 * domain first, transform id second), but reorders cells within a section
 * based on their adjacency.
 *
 * parameters:
 *   mesh       <-> pointer to mesh structure
 *   old_to_new <-> old to new cells numbering (non-ghost in, add ghosts out)
 *   new_to_old <-> new to old cells numbering (non-ghost in, add ghosts out)
 *----------------------------------------------------------------------------*/

static void
_reorder_halo_cells(const cs_mesh_t  *mesh,
                    cs_lnum_t         old_to_new[],
                    cs_lnum_t         new_to_old[])
{
  assert(mesh->halo != NULL);

  int n_c_domains = mesh->halo->n_c_domains;
  int n_transforms = mesh->halo->n_transforms;
  int stride = 2*n_transforms + 2;
  cs_lnum_t *index = mesh->halo->index;
  cs_lnum_t *perio_lst = mesh->halo->perio_lst;
  const cs_lnum_t n_cells = mesh->n_cells;

  cs_lnum_t *keys;

  assert(mesh->n_ghost_cells == mesh->halo->n_elts[1]);

  BFT_MALLOC(keys, mesh->n_ghost_cells*2, cs_lnum_t);

  for (int i = 0; i < n_c_domains; i++) {

    for (cs_lnum_t j = index[2*i]; j < index[2*i+1]; j++) {
      keys[j*2] = i*stride;
      keys[j*2 + 1] = -1;
    }

    for (int t_id = 0; t_id < n_transforms; t_id++) {
      int shift = 4 * n_c_domains * t_id;
      cs_lnum_t s = perio_lst[shift + 4*i];
      cs_lnum_t e = s + perio_lst[shift + 4*i + 1];
      for (cs_lnum_t j = s; j < e; j++)
        keys[j*2] = i*stride + t_id + 1;
    }

    for (cs_lnum_t j = index[2*i+1]; j < index[2*i+2]; j++) {
      keys[j*2] = i*stride + n_transforms + 1;
      keys[j*2 + 1] = -1;
    }

    for (int t_id = 0; t_id < n_transforms; t_id++) {
      int shift = 4 * n_c_domains * t_id;
      cs_lnum_t s = perio_lst[shift + 4*i + 2];
      cs_lnum_t e = s + perio_lst[shift + 4*i + 3];
      for (cs_lnum_t j = s; j < e; j++)
        keys[j*2] = i*stride + n_transforms + t_id + 2;
    }

  } /* End of loop on involved ranks */

  /* Order standard halo by adjacency */

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t gc_id_0 = mesh->i_face_cells[f_id][0] - n_cells;
    cs_lnum_t gc_id_1 = mesh->i_face_cells[f_id][1] - n_cells;
    if (gc_id_0 > -1) {
      if (keys[gc_id_0*2 + 1] < 0)
        keys[gc_id_0*2 + 1] = old_to_new[mesh->i_face_cells[f_id][1]];
    }
    if (gc_id_1 > -1) {
      if (keys[gc_id_1*2 + 1] < 0)
        keys[gc_id_1*2 + 1] = old_to_new[mesh->i_face_cells[f_id][0]];
    }
  }

  /* Order extended halo by adjacency */

  if (mesh->cell_cells_lst != NULL) {
    for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
      for (cs_lnum_t kk = mesh->cell_cells_idx[ii];
           kk < mesh->cell_cells_idx[ii];
           kk++) {
        cs_lnum_t gc_id = mesh->cell_cells_lst[kk] - n_cells;
        if (gc_id > -1) {
          if (keys[gc_id*2 + 1] < 0)
            keys[gc_id*2 + 1] = old_to_new[ii];
        }
      }
    }
  }

  /* Ensure all keys are initialized */

  for (cs_lnum_t j = 0; j < mesh->n_ghost_cells; j++) {
    if (keys[j*2 + 1] < 0)
      keys[j*2 + 1] = j;
  }

  cs_order_lnum_allocated_s(NULL,
                            keys,
                            2,
                            new_to_old + n_cells,
                            mesh->n_ghost_cells);


  BFT_FREE(keys);

  for (cs_lnum_t j = 0; j < mesh->n_ghost_cells; j++) {
    new_to_old[n_cells + j] += n_cells;
    old_to_new[new_to_old[n_cells + j]] = n_cells + j;
  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering of cells.
 *
 * parameters:
 *   mesh       <-> Pointer to global mesh structure
 *   new_to_old <-> Cells renumbering array (size: mesh->n_cells_with_ghosts,
 *                  values at indexes mesh->n_cells and above may be modified)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_cells(cs_mesh_t  *mesh,
                          cs_lnum_t  *new_to_old)
{
  cs_lnum_t  ii, jj, kk, face_id, n_vis, start_id, start_id_old;

  cs_lnum_t  *face_cells_tmp = NULL;
  cs_lnum_t  *new_cell_id = NULL;

  cs_lnum_t  face_cells_max_size = CS_MAX(mesh->n_i_faces*2, mesh->n_b_faces);
  const cs_lnum_t  n_cells = mesh->n_cells;

  /* If no renumbering is present, return */

  if (new_to_old == NULL)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(face_cells_tmp, face_cells_max_size, cs_lnum_t);
  BFT_MALLOC(new_cell_id, mesh->n_cells_with_ghosts, cs_lnum_t);

  /* Build old -> new renumbering */

  for (ii = 0; ii < n_cells; ii++)
    new_cell_id[new_to_old[ii]] = ii;

  for (ii = n_cells; ii < mesh->n_cells_with_ghosts; ii++) {
    new_to_old[ii] = ii;
    new_cell_id[ii] = ii;
  }

  /* Update halo connectivity */

  if (mesh->halo != NULL)
    cs_halo_renumber_cells(mesh->halo, new_cell_id);

  /* Try to optimize ghost cells, based on adjacency */

  if (_renumber_ghost_cells && mesh->halo != NULL) {
    _reorder_halo_cells(mesh, new_cell_id, new_to_old);
    cs_halo_renumber_ghost_cells(mesh->halo, new_to_old);
  }

  /* Update faces -> cells connectivity */

  memcpy(face_cells_tmp,
         mesh->i_face_cells,
         mesh->n_i_faces * 2 * sizeof(cs_lnum_t));

  for (face_id = 0; face_id < mesh->n_i_faces; face_id++) {
    ii = face_cells_tmp[face_id*2];
    jj = face_cells_tmp[face_id*2 + 1];
    mesh->i_face_cells[face_id][0] = new_cell_id[ii];
    mesh->i_face_cells[face_id][1] = new_cell_id[jj];
  }

  if (mesh->n_b_faces > 0) {

    memcpy(face_cells_tmp,
           mesh->b_face_cells,
           mesh->n_b_faces * sizeof(cs_lnum_t));

    for (face_id = 0; face_id < mesh->n_b_faces; face_id++) {
      ii = face_cells_tmp[face_id];
      mesh->b_face_cells[face_id] = new_cell_id[ii];
    }
  }

  /* Update cell -> cells connectivity for extended neighborhood */

  if (mesh->cell_cells_lst != NULL) {

    cs_lnum_t *cell_cells_idx_old, *cell_cells_lst_old;
    const cs_lnum_t cell_cells_lst_size = mesh->cell_cells_idx[n_cells];

    BFT_MALLOC(cell_cells_idx_old, n_cells + 1, cs_lnum_t);
    BFT_MALLOC(cell_cells_lst_old, cell_cells_lst_size, cs_lnum_t);

    memcpy(cell_cells_idx_old,
           mesh->cell_cells_idx,
           (n_cells + 1)*sizeof(cs_lnum_t));
    memcpy(cell_cells_lst_old,
           mesh->cell_cells_lst,
           cell_cells_lst_size*sizeof(cs_lnum_t));

    mesh->cell_cells_idx[0] = 0;
    start_id = 0;

    for (ii = 0; ii < n_cells; ii++) {

      jj = new_to_old[ii];
      n_vis = cell_cells_idx_old[jj+1] - cell_cells_idx_old[jj];
      start_id_old = cell_cells_idx_old[jj];

      for (kk = 0; kk < n_vis; kk++)
        mesh->cell_cells_lst[start_id + kk]
          = new_cell_id[cell_cells_lst_old[start_id_old + kk]];

      start_id += n_vis;
      mesh->cell_cells_idx[ii + 1] = start_id;
    }

    BFT_FREE(cell_cells_lst_old);
    BFT_FREE(cell_cells_idx_old);
  }

  /* Update ghost cell -> vertices connectivity for extended neighborhood */

  if (mesh->gcell_vtx_lst != NULL) {

    cs_lnum_t *gcell_vtx_idx_old, *gcell_vtx_lst_old;
    const cs_lnum_t  n_ghost_cells = mesh->n_ghost_cells;
    const cs_lnum_t gcell_vtx_lst_size
      = mesh->gcell_vtx_idx[mesh->n_ghost_cells];

    BFT_MALLOC(gcell_vtx_idx_old, n_ghost_cells + 1, cs_lnum_t);
    BFT_MALLOC(gcell_vtx_lst_old, gcell_vtx_lst_size, cs_lnum_t);

    memcpy(gcell_vtx_idx_old,
           mesh->gcell_vtx_idx,
           (n_ghost_cells + 1)*sizeof(cs_lnum_t));
    memcpy(gcell_vtx_lst_old,
           mesh->gcell_vtx_lst,
           gcell_vtx_lst_size*sizeof(cs_lnum_t));

    mesh->gcell_vtx_idx[0] = 0;
    start_id = 0;

    for (ii = 0; ii < n_ghost_cells; ii++) {

      jj = new_to_old[n_cells + ii] - n_cells;
      n_vis = gcell_vtx_idx_old[jj+1] - gcell_vtx_idx_old[jj];
      start_id_old = gcell_vtx_idx_old[jj];

      for (kk = 0; kk < n_vis; kk++)
        mesh->gcell_vtx_lst[start_id + kk]
          = gcell_vtx_lst_old[start_id_old + kk];

      start_id += n_vis;
      mesh->gcell_vtx_idx[ii + 1] = start_id;
    }

    BFT_FREE(gcell_vtx_lst_old);
    BFT_FREE(gcell_vtx_idx_old);
  }

  /* Free work arrays */

  BFT_FREE(new_cell_id);
  BFT_FREE(face_cells_tmp);

  /* Update list of boundary cells */

  cs_mesh_update_b_cells(mesh);

  /* Update cell families and global numbering */

  _update_family(n_cells, new_to_old, mesh->cell_family);

  _update_global_num(n_cells, new_to_old, &(mesh->global_cell_num));

  /* Update parent cell numbers for post-processing meshes
     that may already have been built; Post-processing meshes
     built after renumbering will have correct parent numbers */

  cs_post_renum_cells(new_to_old);
}

/*----------------------------------------------------------------------------
 * Apply renumbering to a face -> vertices connectivity.
 *
 * parameters:
 *   n_faces         <-- Number of faces
 *   face_vtx_idx    <-> Face -> vertices index
 *   face_vtx        <-- Face vertices
 *   new_to_old      <-- Faces renumbering array
 *----------------------------------------------------------------------------*/

static void
_update_face_vertices(cs_lnum_t         n_faces,
                      cs_lnum_t        *face_vtx_idx,
                      cs_lnum_t        *face_vtx,
                      const cs_lnum_t  *new_to_old)
{
  if (new_to_old != NULL && face_vtx != NULL) {

    cs_lnum_t ii, jj, kk, n_vtx, start_id, start_id_old;
    cs_lnum_t *face_vtx_idx_old, *face_vtx_old;

    const cs_lnum_t connect_size = face_vtx_idx[n_faces];

    BFT_MALLOC(face_vtx_idx_old, n_faces + 1, cs_lnum_t);
    BFT_MALLOC(face_vtx_old, connect_size, cs_lnum_t);

    memcpy(face_vtx_idx_old, face_vtx_idx, (n_faces+1)*sizeof(cs_lnum_t));
    memcpy(face_vtx_old, face_vtx, connect_size*sizeof(cs_lnum_t));

    face_vtx_idx[0] = 0;
    start_id = 0;

    for (ii = 0; ii < n_faces; ii++) {

      jj = new_to_old[ii];
      n_vtx = face_vtx_idx_old[jj+1] - face_vtx_idx_old[jj];
      start_id_old = face_vtx_idx_old[jj];

      for (kk = 0; kk < n_vtx; kk++)
        face_vtx[start_id + kk] = face_vtx_old[start_id_old + kk];

      start_id += n_vtx;
      face_vtx_idx[ii + 1] = start_id;
    }

    BFT_FREE(face_vtx_idx_old);
    BFT_FREE(face_vtx_old);
  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering of interior faces.
 *
 * parameters:
 *   mesh          <-> Pointer to global mesh structure
 *   new_to_old_i  <-- Interior faces renumbering array
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_i_faces(cs_mesh_t        *mesh,
                            const cs_lnum_t  *new_to_old_i)
{
  if (new_to_old_i != NULL) {

    cs_lnum_t  face_id, face_id_old;

    cs_lnum_2_t  *i_face_cells_old;

    const cs_lnum_t  n_i_faces = mesh->n_i_faces;

   /* Allocate Work array */

    BFT_MALLOC(i_face_cells_old, n_i_faces, cs_lnum_2_t);

    /* Update faces -> cells connectivity */

    memcpy(i_face_cells_old, mesh->i_face_cells, n_i_faces*sizeof(cs_lnum_2_t));

    for (face_id = 0; face_id < n_i_faces; face_id++) {
      face_id_old = new_to_old_i[face_id];
      mesh->i_face_cells[face_id][0] = i_face_cells_old[face_id_old][0];
      mesh->i_face_cells[face_id][1] = i_face_cells_old[face_id_old][1];
    }

    BFT_FREE(i_face_cells_old);

    /* Update faces -> vertices connectivity */

    _update_face_vertices(n_i_faces,
                          mesh->i_face_vtx_idx,
                          mesh->i_face_vtx_lst,
                          new_to_old_i);

    /* Update face families and global numbering */

    _update_family(n_i_faces, new_to_old_i, mesh->i_face_family);

    if (mesh->i_face_r_gen != NULL)
      _update_r_gen(n_i_faces, new_to_old_i, mesh->i_face_r_gen);

    _update_global_num(n_i_faces, new_to_old_i, &(mesh->global_i_face_num));

    /* Update parent face numbers for post-processing meshes
       that may already have been built; Post-processing meshes
       built after renumbering will have correct parent numbers */

    cs_post_renum_faces(new_to_old_i, NULL);

  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering of boundary faces.
 *
 * parameters:
 *   mesh         <-> Pointer to global mesh structure
 *   new_to_old_b <-- Boundary faces renumbering array
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_b_faces(cs_mesh_t        *mesh,
                            const cs_lnum_t  *new_to_old_b)
{
  if (new_to_old_b != NULL) {

    cs_lnum_t  face_id, face_id_old;

    cs_lnum_t  *b_face_cells_old = NULL;

    const cs_lnum_t  n_b_faces = mesh->n_b_faces;

    /* Allocate Work array */

    BFT_MALLOC(b_face_cells_old, n_b_faces, cs_lnum_t);

    /* Update faces -> cells connectivity */

    memcpy(b_face_cells_old, mesh->b_face_cells, n_b_faces*sizeof(cs_lnum_t));

    for (face_id = 0; face_id < n_b_faces; face_id++) {
      face_id_old = new_to_old_b[face_id];
      mesh->b_face_cells[face_id] = b_face_cells_old[face_id_old];
    }

    BFT_FREE(b_face_cells_old);

    /* Update faces -> vertices connectivity */

    _update_face_vertices(n_b_faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          new_to_old_b);

    /* Update face families and global numbering */

    _update_family(n_b_faces, new_to_old_b, mesh->b_face_family);

    _update_global_num(n_b_faces, new_to_old_b, &(mesh->global_b_face_num));

    /* Update parent face numbers for post-processing meshes
       that may already have been built; Post-processing meshes
       built after renumbering will have correct parent numbers */

    cs_post_renum_faces(NULL, new_to_old_b);

  }
}

/*----------------------------------------------------------------------------
 * Apply renumbering to a face -> vertices connectivity when the vertices
 * have been renumbered.
 *
 * parameters:
 *   n_faces         <-- Number of faces
 *   face_vtx_idx    <-- Face -> vertices index
 *   face_vtx        <-> Face vertices
 *   old_to_new      <-- Vertices renumbering array
 *----------------------------------------------------------------------------*/

static void
_update_face_vertices_v(cs_lnum_t         n_faces,
                        const cs_lnum_t  *face_vtx_idx,
                        cs_lnum_t        *face_vtx,
                        const cs_lnum_t  *old_to_new)
{
  if (old_to_new == NULL || face_vtx == NULL)
    return;

  const cs_lnum_t connect_size = face_vtx_idx[n_faces];

  /* face_vtx_idx is not modified by the vertices renumbering */

  cs_lnum_t *face_vtx_old = NULL;
  BFT_MALLOC(face_vtx_old, connect_size, cs_lnum_t);

  memcpy(face_vtx_old, face_vtx, connect_size*sizeof(cs_lnum_t));

  for (cs_lnum_t jj = 0; jj < connect_size; jj++)
    face_vtx[jj] = old_to_new[face_vtx_old[jj]];

  BFT_FREE(face_vtx_old);
}

/*----------------------------------------------------------------------------
 * Apply renumbering of vertices.
 *
 * parameters:
 *   mesh         <-> Pointer to global mesh structure
 *   n2o_v        <-- Vertices renumbering array
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_vertices(cs_mesh_t        *mesh,
                             const cs_lnum_t  *n2o_v)
{
  if (n2o_v == NULL)
    return;

  const cs_lnum_t  n_vertices = mesh->n_vertices;

  /* Allocate Work array */

  cs_real_t  *vtx_coord_old = NULL;

  BFT_MALLOC(vtx_coord_old, 3*n_vertices, cs_real_t);

  /* Update the vertex coordinates */

  cs_array_real_copy(3*n_vertices, mesh->vtx_coord, vtx_coord_old);

  cs_array_real_copy_subset(n_vertices, 3, n2o_v,
                            CS_ARRAY_SUBSET_IN, /* apply n2o_v to ref */
                            vtx_coord_old,      /* ref (in) */
                            mesh->vtx_coord);   /* dest (out) */

  BFT_FREE(vtx_coord_old);

  /* Update refinement level */

  if (mesh->vtx_r_gen != NULL)
    _update_r_gen(n_vertices, n2o_v, mesh->vtx_r_gen);

  /* Build the old --> new numbering */

  cs_lnum_t  *o2n_v = NULL;

  BFT_MALLOC(o2n_v, n_vertices, cs_lnum_t);
  for (cs_lnum_t new_id = 0; new_id < n_vertices; new_id++)
    o2n_v[n2o_v[new_id]] = new_id; /* old_id = n2o_v[new_id] */

  /* Update faces -> vertices connectivity for boundary faces */

  _update_face_vertices_v(mesh->n_b_faces,
                          mesh->b_face_vtx_idx,
                          mesh->b_face_vtx_lst,
                          o2n_v);

  /* Update faces -> vertices connectivity for interior faces */

  _update_face_vertices_v(mesh->n_i_faces,
                          mesh->i_face_vtx_idx,
                          mesh->i_face_vtx_lst,
                          o2n_v);

  /* Update global numbering */

  _update_global_num(n_vertices, n2o_v, &(mesh->global_vtx_num));

  /* Interfaces for vertices  */

  if (mesh->vtx_interfaces != NULL)
    cs_interface_set_renumber(mesh->vtx_interfaces, o2n_v);

  BFT_FREE(o2n_v);

  /* Update the numbering of parent vertices for post-processing meshes
     that may already have been built; Post-processing meshes
     built after renumbering will have correct parent numbers

     This should not be the standard behavior. Otherwise, this is
     something to do.
  */

  /* Sanity check: if not true, one should update this members of the mesh
   structure. This should not be the case since the renumbering should be called
   at the first steps of the computation */
  assert(mesh->gcell_vtx_idx == NULL);
  assert(mesh->gcell_vtx_lst == NULL);
  assert(mesh->vtx_range_set == NULL);
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

#if defined(HAVE_IBM_RENUMBERING_LIB)

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces and cells for multiple threads.
 *
 * parameters:
 *   mesh  <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_for_threads_ibm(cs_mesh_t  *mesh)
{
}

#endif /* defined(HAVE_IBM_RENUMBERING_LIB) */

/*----------------------------------------------------------------------------
 * Create a CSR graph structure from a native face-based connectivity.
 *
 * parameters:
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity
 *
 * returns:
 *   pointer to allocated adjacency structure.
 *----------------------------------------------------------------------------*/

static cs_adjacency_t *
_c2c_from_face_cell(cs_lnum_t         n_cells_ext,
                    cs_lnum_t         n_faces,
                    const cs_lnum_t  *face_cell)
{
  cs_lnum_t ii, jj, f_id;

  cs_lnum_t  *ccount = NULL;
  bool unique_faces = true;

  cs_adjacency_t  *a = cs_adjacency_create(0, 0, n_cells_ext);

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, a->n_elts, cs_lnum_t);

  for (ii = 0; ii < a->n_elts; ii++)
    ccount[ii] = 0;

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2];
    jj = face_cell[f_id*2 + 1];
    ccount[ii] += 1;
    ccount[jj] += 1;
  }

  a->idx[0] = 0;
  for (ii = 0; ii < a->n_elts; ii++) {
    a->idx[ii+1] = a->idx[ii] + ccount[ii];
    ccount[ii] = 0;
  }

  /* Build structure */

  BFT_MALLOC(a->ids, (a->idx[a->n_elts]), cs_lnum_t);

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id*2];
    jj = face_cell[f_id*2 + 1];
    a->ids[a->idx[ii] + ccount[ii]] = jj;
    ccount[ii] += 1;
    a->ids[a->idx[jj] + ccount[jj]] = ii;
    ccount[jj] += 1;
  }

  BFT_FREE(ccount);

  /* Sort line elements by column id (for better access patterns) */

  unique_faces = cs_sort_indexed(a->n_elts, a->idx, a->ids);

  /* Compact elements if necessary */

  if (unique_faces == false) {

    cs_lnum_t *tmp_idx = NULL;
    cs_lnum_t  kk = 0;

    BFT_MALLOC(tmp_idx, a->n_elts+1, cs_lnum_t);
    memcpy(tmp_idx, a->idx, (a->n_elts+1)*sizeof(cs_lnum_t));

    kk = 0;

    for (ii = 0; ii < a->n_elts; ii++) {
      cs_lnum_t *col_id = a->ids + a->idx[ii];
      cs_lnum_t n_cols = a->idx[ii+1] - a->idx[ii];
      cs_lnum_t col_id_prev = -1;
      a->idx[ii] = kk;
      for (jj = 0; jj < n_cols; jj++) {
        if (col_id_prev != col_id[jj]) {
          a->ids[kk++] = col_id[jj];
          col_id_prev = col_id[jj];
        }
      }
    }
    a->idx[a->n_elts] = kk;

    assert(a->idx[a->n_elts] < tmp_idx[a->n_elts]);

    BFT_FREE(tmp_idx);
    BFT_REALLOC(a->ids, (a->idx[a->n_elts]), cs_lnum_t);

  }

  return a;
}

/*----------------------------------------------------------------------------
 * Create a CSR graph structure from a native face-based conectivity.
 *
 * parameters:
 *   n_cells_ext <-- Local number of cells + ghost cells sharing a face
 *   n_faces     <-- Local number of faces
 *   face_cell   <-- Face -> cells connectivity
 *
 * returns:
 *   pointer to adjacency structure
 *----------------------------------------------------------------------------*/

static cs_adjacency_t *
_c2f_from_face_cell(cs_lnum_t           n_cells_ext,
                    cs_lnum_t           n_faces,
                    const cs_lnum_2_t  *face_cell)
{
  cs_lnum_t ii, jj, f_id;

  cs_lnum_t  *ccount = NULL;

  cs_adjacency_t  *a = cs_adjacency_create(0, 0, n_cells_ext);

  /* Count number of nonzero elements per row */

  BFT_MALLOC(ccount, a->n_elts, cs_lnum_t);

  for (ii = 0; ii < a->n_elts; ii++)
    ccount[ii] = 0;

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id][0];
    jj = face_cell[f_id][1];
    ccount[ii] += 1;
    ccount[jj] += 1;
  }

  a->idx[0] = 0;
  for (ii = 0; ii < a->n_elts; ii++) {
    a->idx[ii+1] = a->idx[ii] + ccount[ii];
    ccount[ii] = 0;
  }

  /* Build structure */

  BFT_MALLOC(a->ids, (a->idx[a->n_elts]), cs_lnum_t);

  for (f_id = 0; f_id < n_faces; f_id++) {
    ii = face_cell[f_id][0];
    jj = face_cell[f_id][1];
    a->ids[a->idx[ii] + ccount[ii]] = f_id;
    ccount[ii] += 1;
    a->ids[a->idx[jj] + ccount[jj]] = f_id;
    ccount[jj] += 1;
  }

  BFT_FREE(ccount);

  return a;
}

/*----------------------------------------------------------------------------
 * Classify ghost cells according to their position in the halo
 *
 * The cell_class value will based on the halo section (separated based on
 * domain first, transform id second; extended ghost cells do not intervene
 * here), with a base class id of 2.
 *
 * parameters:
 *   mesh       <-> pointer to mesh structure
 *   halo_class --> resulting cell classification
 *----------------------------------------------------------------------------*/

static void
_classify_halo_cells(const cs_mesh_t  *mesh,
                     int               halo_class[])
{
  assert(mesh->halo != NULL);

  int n_c_domains = mesh->halo->n_c_domains;
  int n_transforms = mesh->halo->n_transforms;
  int stride = n_transforms + 2;
  cs_lnum_t *index = mesh->halo->index;
  cs_lnum_t *perio_lst = mesh->halo->perio_lst;

  for (int i = 0; i < n_c_domains; i++) {

    for (cs_lnum_t j = index[2*i]; j < index[2*i+1]; j++)
      halo_class[j] = i*stride + 2;

    for (int t_id = 0; t_id < n_transforms; t_id++) {
      int shift = 4 * n_c_domains * t_id;
      cs_lnum_t s = perio_lst[shift + 4*i];
      cs_lnum_t e = s + perio_lst[shift + 4*i + 1];
      for (cs_lnum_t k = s; k < e; k++)
        halo_class[k] = i*stride + t_id + 3;
    }

    /* Extended ghost cells should not intervene here;
       we initialize for safety, but will not bother
       with transforms */

    for (cs_lnum_t j = index[2*i+1]; j < index[2*i+2]; j++)
      halo_class[j] = i*stride + (stride - 1);

  } /* End of loop on involved ranks */

}

/*----------------------------------------------------------------------------
 * Classify (non-ghost) cells according to their neighbor domain
 *
 * For cells connected to one or more ghost cells, the cell_class value
 * will be based on the halo section (separated based on domain first,
 * transform id second; extended ghost cells do not intervene here),
 * with a base class id of 2. For cells connected to a boundary face
 * that will be later joined, the cell class will be 1. For other cells,
 * the cell class will be 0.
 *
 * parameters:
 *   mesh       <->- pointer to mesh structure
 *   cell_class --> resulting cell classification
 *----------------------------------------------------------------------------*/

static void
_classify_cells_by_neighbor(const cs_mesh_t  *mesh,
                            int               cell_class[])
{
  for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
    cell_class[i] = 0;

  {
    bool *b_select_flag;
    BFT_MALLOC(b_select_flag, mesh->n_b_faces, bool);

    cs_join_mark_selected_faces(mesh, false, b_select_flag);

    for (cs_lnum_t face_id = 0; face_id < mesh->n_b_faces; face_id++) {
      if (b_select_flag[face_id] == true)
        cell_class[mesh->b_face_cells[face_id]] = 1;
    }

    BFT_FREE(b_select_flag);
  }

  if (mesh->halo != NULL) {

    int *halo_class;
    BFT_MALLOC(halo_class, mesh->n_ghost_cells, int);

    _classify_halo_cells(mesh, halo_class);

    /* Now apply this to cells */

    for (cs_lnum_t face_id = 0; face_id < mesh->n_i_faces; face_id++) {
      for (int s_id = 0; s_id < 2; s_id++) {
        cs_lnum_t c_id = mesh->i_face_cells[face_id][s_id];
        if (c_id >= mesh->n_cells) {
          int min_class = halo_class[c_id - mesh->n_cells];
          int adj_c_id = mesh->i_face_cells[face_id][(s_id+1)%2];
          cell_class[adj_c_id] = CS_MAX(cell_class[adj_c_id], min_class);
        }
      }
    }

    BFT_FREE(halo_class);
  }
}

/*----------------------------------------------------------------------------
 * Initial numbering of interior faces, based on cell adjacency.
 *
 * parameters:
 *   mesh          <-- pointer to global mesh structure
 *   base_ordering <-- pre-ordering of interior faces by
 *                     lowest or highest adjacent cell id
 *   order         --> ordering of interior faces based on cell adjacency
 *
 * returns:
 *   number of faces in ordering before first face adjacent to halo
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_order_i_faces_by_cell_adjacency(const cs_mesh_t         *mesh,
                                 cs_renumber_ordering_t   base_ordering,
                                 cs_lnum_t               *order)
{
  cs_lnum_t  *faces_keys;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;

  cs_lnum_t n_no_adj_halo = 0;

  if (base_ordering == CS_RENUMBER_ADJACENT_LOW) {

    /* If cells adjacent to halo cells are last, adjacent interior faces
       will also be last; otherwise, this may be forced */

    if (   mesh->halo != NULL
        && _i_faces_adjacent_to_halo_last
        && (! _cells_adjacent_to_halo_last)) {

      int  *halo_class;

      BFT_MALLOC(faces_keys, mesh->n_i_faces * 3, cs_lnum_t);
      BFT_MALLOC(halo_class, mesh->n_ghost_cells, int);

      _classify_halo_cells(mesh, halo_class);

      /* Build lexical ordering of faces */

#     pragma omp parallel for reduction(+:n_no_adj_halo)
      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        if (c_id_0 >= n_cells)
          faces_keys[f_id*3] = halo_class[c_id_0 - n_cells];
        else if (c_id_1 > n_cells)
          faces_keys[f_id*3] = halo_class[c_id_1 - n_cells];
        else {
          faces_keys[f_id*3] = 0;
          n_no_adj_halo += 1;
        }
        if (c_id_0 < c_id_1) {
          faces_keys[f_id*3 + 1] = c_id_0;
          faces_keys[f_id*3 + 2] = c_id_1;
        }
        else {
          faces_keys[f_id*3 + 1] = c_id_1;
          faces_keys[f_id*3 + 2] = c_id_0;
        }
      }

      BFT_FREE(halo_class);

      cs_order_lnum_allocated_s(NULL,
                                faces_keys,
                                3,
                                order,
                                n_i_faces);

    }
    else {

      BFT_MALLOC(faces_keys, mesh->n_i_faces*2, cs_lnum_t);

      /* Build lexical ordering of faces */

#     pragma omp parallel for  if (n_i_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
        cs_lnum_t c_id_0 = i_face_cells[f_id][0];
        cs_lnum_t c_id_1 = i_face_cells[f_id][1];
        if (c_id_0 < c_id_1) {
          faces_keys[f_id*2]     = c_id_0;
          faces_keys[f_id*2 + 1] = c_id_1;
        }
        else {
          faces_keys[f_id*2]     = c_id_1;
          faces_keys[f_id*2 + 1] = c_id_0;
        }
      }

      cs_order_lnum_allocated_s(NULL,
                                faces_keys,
                                2,
                                order,
                                n_i_faces);

      if (_i_faces_adjacent_to_halo_last) {

        for (cs_lnum_t i = 0; i < n_i_faces; i++) {
          cs_lnum_t f_id = order[i];
          if (faces_keys[f_id*2 + 1] > n_cells)
            break;
          else
            n_no_adj_halo += 1;
        }

      }

    }

  }

  else { /* if (i_faces_base_ordering == CS_RENUMBER_ADJACENT_HIGH) */

    BFT_MALLOC(faces_keys, mesh->n_i_faces*2, cs_lnum_t);

    /* Build lexical ordering of faces */

#   pragma omp parallel for  if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      if (c_id_0 < c_id_1) {
        faces_keys[f_id*2]     = c_id_1;
        faces_keys[f_id*2 + 1] = c_id_0;
      }
      else {
        faces_keys[f_id*2]     = c_id_0;
        faces_keys[f_id*2 + 1] = c_id_1;
      }
    }

    cs_order_lnum_allocated_s(NULL,
                              faces_keys,
                              2,
                              order,
                              n_i_faces);

    if (_i_faces_adjacent_to_halo_last) {

      for (cs_lnum_t i = 0; i < n_i_faces; i++) {
        cs_lnum_t f_id = order[i];
        if (faces_keys[f_id*2] > n_cells)
          break;
        else
          n_no_adj_halo += 1;
      }

    }

  } /* end of test on i_faces_base_ordering */

  BFT_FREE(faces_keys);

  if (mesh->halo == NULL)
    n_no_adj_halo = n_i_faces;

  return n_no_adj_halo;
}

/*----------------------------------------------------------------------------
 * Initial numbering of interior faces, based on cell adjacency.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *
 * returns:
 *   number of faces in ordering before first face adjacent to halo
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_renumber_i_faces_by_cell_adjacency(cs_mesh_t  *mesh)
{
  cs_lnum_t  *new_to_old_i;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;

  BFT_MALLOC(new_to_old_i, n_i_faces, cs_lnum_t);

  /* If cells adjacent to halo cells are last, adjacent interior faces
     will also be last; otherwise, this may be forced */

  cs_lnum_t n_no_adj_halo
    = _order_i_faces_by_cell_adjacency(mesh,
                                       _i_faces_base_ordering,
                                       new_to_old_i);

  /* Check numbering is non trivial */

  {
    cs_lnum_t f_id = 0;
    while (f_id < n_i_faces) {
      if (new_to_old_i[f_id] != f_id)
        break;
      else
        f_id++;
    }
    if (f_id == n_i_faces)
      BFT_FREE(new_to_old_i);
  }

  /* Update connectivity */

  if (new_to_old_i != NULL)
    _cs_renumber_update_i_faces(mesh, new_to_old_i);

  BFT_FREE(new_to_old_i);

  return n_no_adj_halo;
}

/*----------------------------------------------------------------------------
 * Initial numbering of boundary faces, based on cell adjacency.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_b_faces_by_cell_adjacency(cs_mesh_t  *mesh)
{
  cs_lnum_t  *new_to_old_b;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  BFT_MALLOC(new_to_old_b, n_b_faces, cs_lnum_t);

  cs_lnum_t *fc_num;
  BFT_MALLOC(fc_num, mesh->n_b_faces*2, cs_lnum_t);

  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
    fc_num[ii*2] = mesh->b_face_cells[ii];
    fc_num[ii*2+1] = ii;
  }

  cs_order_lnum_allocated_s(NULL, fc_num, 2, new_to_old_b, n_b_faces);

  BFT_FREE(fc_num);

  /* Check numbering is non trivial */

  {
    cs_lnum_t f_id = 0;
    while (f_id < n_b_faces) {
      if (new_to_old_b[f_id] != f_id)
        break;
      else
        f_id++;
    }
    if (f_id == n_b_faces)
      BFT_FREE(new_to_old_b);
  }

  /* Update connectivity */

  if (new_to_old_b != NULL)
    _cs_renumber_update_b_faces(mesh, new_to_old_b);

  BFT_FREE(new_to_old_b);
}

/*----------------------------------------------------------------------------
 * Build groups including independent faces.
 *
 * parameters:
 *   max_group_size  <-- max group size
 *   n_faces         <-- number of faces
 *   n_cells_ext     <-- local number of cells + ghost cells sharing a face
 *   n_faces         <-- local number of faces
 *   face_cell       <-- face -> cells connectivity
 *   new_to_old      --> new -> old face renumbering
 *   n_groups        --> number of groups
 *   group_size      --> array containing the sizes of groups
 *----------------------------------------------------------------------------*/

static void
_independent_face_groups(cs_lnum_t            max_group_size,
                         cs_lnum_t            n_cells_ext,
                         cs_lnum_t            n_faces,
                         const cs_lnum_2_t   *face_cell,
                         cs_lnum_t           *new_to_old,
                         int                 *n_groups,
                         cs_lnum_t          **group_size)
{
  cs_lnum_t f_id, i, j, k;
  cs_lnum_t *group_face_ids = NULL, *face_marker = NULL;
  cs_lnum_t *old_to_new = NULL;

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

  cs_adjacency_t *cell_faces = _c2f_from_face_cell(n_cells_ext,
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
          cs_lnum_t c_id[2] = {face_cell[f_cmp][0], face_cell[f_cmp][1]};

          for (j = 0; j < 2; j++) {
            cs_lnum_t start_id = cell_faces->idx[c_id[j]];
            cs_lnum_t end_id = cell_faces->idx[c_id[j] + 1];
            for (k = start_id; k < end_id; k++) {
              if (cell_faces->ids[k] == f_id) {
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

  cs_adjacency_destroy(&cell_faces);

  BFT_FREE(face_marker);
  BFT_FREE(group_face_ids);

  BFT_REALLOC(_group_size, group_id, cs_lnum_t);

  /* Set return values */

  for (f_id = 0; f_id < n_faces; f_id++)
    new_to_old[old_to_new[f_id]] = f_id;

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
 * Pre-assign faces to threads of a given group for the multipass algorithm
 * to improve load balance.
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

    if (_i_faces_base_ordering == CS_RENUMBER_ADJACENT_LOW) {
      if (   c_id_0 >= t_cell_index[t_id]
          && c_id_1 < t_cell_index[t_id+1]) {
        f_t_id[f_id] = t_id + g_id*n_i_threads;
        n_t_faces[t_id] += 1;
        t_face_last[t_id] = fl_id;
      }
      else
        f_t_id[f_id] = -1;
    }
    else { /* _i_faces_base_ordering == CS_RENUMBER_ADJACENT_HIGH */
      if (   c_id_0 < t_cell_index[t_id+1]
          && c_id_1 >= t_cell_index[t_id]) {
        assert(c_id_0 >= t_cell_index[t_id]);
        f_t_id[f_id] = t_id + g_id*n_i_threads;
        n_t_faces[t_id] += 1;
        t_face_last[t_id] = fl_id;
      }
      else
        f_t_id[f_id] = -1;
    }

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
_renum_face_multipass_redistribute(int                          n_i_threads,
                                   int                          n_g_i_threads,
                                   int                          g_id,
                                   double                       relax,
                                   cs_lnum_t                    faces_list_size,
                                   const cs_lnum_t    *restrict faces_list,
                                   const cs_lnum_2_t  *restrict l_face_cells,
                                   int                *restrict f_t_id,
                                   cs_lnum_t          *restrict n_t_faces,
                                   cs_lnum_t          *restrict t_face_last,
                                   cs_lnum_t          *restrict t_cell_index)
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

  if (_i_faces_base_ordering == CS_RENUMBER_ADJACENT_LOW) {

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

  }
  else { /* _i_faces_base_ordering == CS_RENUMBER_ADJACENT_HIGH */

    for (fl_id = 0; fl_id < faces_list_size; fl_id++) {

      f_id = faces_list[fl_id];

      c_id_0 = l_face_cells[f_id][0];
      c_id_1 = l_face_cells[f_id][1];

      if (new_cell_id[c_id_0] < 0)
        new_cell_id[c_id_0] = n_f_cells_new++;
      if (new_cell_id[c_id_1] < 0)
        new_cell_id[c_id_1] = n_f_cells_new++;

      if (new_cell_id[c_id_0] < new_cell_id[c_id_1]) {
        l_face_cells[f_id][0] = new_cell_id[c_id_1];
        l_face_cells[f_id][1] = new_cell_id[c_id_0];
      }
      else {
        l_face_cells[f_id][0] = new_cell_id[c_id_0];
        l_face_cells[f_id][1] = new_cell_id[c_id_1];
      }

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
 *   mesh                 <-> pointer to global mesh structure
 *   n_i_threads          <-- number of threads required for interior faces
 *   group_size           <-- target group size
 *   new_to_old_i         --> interior faces renumbering array
 *   n_groups             --> number of groups of graph edges (interior faces)
 *   n_no_adj_halo_groups --> number of groups with faces not adjacent to
 *                            halo cells
 *   group_index          --> group/thread index
 *
 * returns:
 *   0 on success, -1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_face_multipass(cs_mesh_t    *mesh,
                      int           n_i_threads,
                      cs_lnum_t     new_to_old_i[],
                      int          *n_groups,
                      int          *n_no_adj_halo_groups,
                      cs_lnum_t   **group_index)
{
  int g_id;
  cs_lnum_t f_id;

  cs_lnum_t n_f_cells = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_faces = mesh->n_i_faces;
  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)mesh->i_face_cells;

  double redistribute_relaxation_factor = 0.5;

  cs_lnum_t new_count = 0;
  cs_lnum_t faces_list_size = n_faces, faces_list_size_new = 0;
  cs_lnum_t faces_list_assign_size = faces_list_size;

  cs_lnum_t   _n_groups = 0;
  cs_lnum_t   _n_no_adj_halo_groups = 0;
  cs_lnum_t  *_group_index = NULL;

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

  /* Build lexical ordering of faces
     (possibly forcing faces ajacent to ghost cells last) */

  cs_lnum_t n_no_adj_halo
    = _order_i_faces_by_cell_adjacency(mesh,
                                       _i_faces_base_ordering,
                                       faces_list);

  /* Initialize final sorting keys
     (second key here based on initial ordering, first key
     will be defined later based on group and thread ids) */

  if (_i_faces_base_ordering == CS_RENUMBER_ADJACENT_LOW) {
#   pragma omp parallel for  if (n_faces > CS_THR_MIN)
    for (f_id = 0; f_id < n_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
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
  }
  else { /* _i_faces_base_ordering == CS_RENUMBER_ADJACENT_HIGH */
#   pragma omp parallel for  if (n_faces > CS_THR_MIN)
    for (f_id = 0; f_id < n_faces; f_id++) {
      cs_lnum_t c_id_0 = i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = i_face_cells[f_id][1];
      if (c_id_0 < c_id_1) {
        l_face_cells[f_id][0] = c_id_1;
        l_face_cells[f_id][1] = c_id_0;
      }
      else {
        l_face_cells[f_id][0] = c_id_0;
        l_face_cells[f_id][1] = c_id_1;
      }
      f_t_id[f_id] = -1;
    }
  }

  /* Add groups as required */

  for (g_id = 0; faces_list_size > _min_i_subset_size; g_id++) {

    int group_size = n_f_cells / n_i_threads;
    int n_g_i_threads = n_i_threads;

    /* If the number of faces not adjacent to ghost cells represents
       a significant enough portion of faces, place faces adjacent
       to ghost cells in separate (later) groups */

    if (n_no_adj_halo > faces_list_size / 2)
      faces_list_assign_size = n_no_adj_halo;
    else {
      faces_list_assign_size = faces_list_size;
      n_no_adj_halo = 0;
    }

    /* Reduce number of threads for this level if required to
       ensure sufficient work per thread */

    if (faces_list_assign_size / _min_i_subset_size  < n_g_i_threads) {
      n_g_i_threads = faces_list_assign_size / _min_i_subset_size;
      if (! (faces_list_assign_size % _min_i_subset_size))
        n_g_i_threads += 1;
    }

    /* Get an initial edge distribution */

    t_cell_index[0] = 0;
    for (int t_id = 1; t_id < n_g_i_threads; t_id++) {
      t_cell_index[t_id] = t_cell_index[t_id-1] + group_size;
      if (t_cell_index[t_id] > n_f_cells)
        t_cell_index[t_id] = n_f_cells;
    }
    t_cell_index[n_g_i_threads] = n_f_cells;

    /* Pre-assign threads to faces (initial distribution) */

    _renum_face_multipass_assign(n_i_threads,
                                 n_g_i_threads,
                                 g_id,
                                 faces_list_assign_size,
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
                                       faces_list_assign_size,
                                       faces_list,
                                       (const cs_lnum_2_t *restrict)l_face_cells,
                                       f_t_id,
                                       n_t_faces,
                                       t_face_last,
                                       t_cell_index);

    if (faces_list_assign_size < faces_list_size) {

      _n_no_adj_halo_groups = g_id+1;

      for (cs_lnum_t fl_id = faces_list_assign_size;
           fl_id < faces_list_size;
           fl_id++) {
        f_id = faces_list[fl_id];
        f_t_id[f_id] = - 1;
      }

    }

    /* Update list of remaining faces */

    for (cs_lnum_t fl_id = 0; fl_id < faces_list_size; fl_id++) {

      f_id = faces_list[fl_id];

      if (f_t_id[f_id] < 0)
        faces_list[faces_list_size_new++] = f_id;
      else
        new_to_old_i[new_count++] = f_id;

    }

    if (n_no_adj_halo > 0)
      n_no_adj_halo -= (faces_list_size - faces_list_size_new);

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

    for (cs_lnum_t fl_id = 0; fl_id < faces_list_size; fl_id++) {
      f_id = faces_list[fl_id];
      f_t_id[f_id] = g_id*n_i_threads;
      new_to_old_i[new_count++] = f_id;
    }

    g_id += 1;

    n_t_faces[0] = faces_list_size;
    for (int t_id = 1; t_id < n_i_threads; t_id++)
      n_t_faces[t_id] = 0;

  }

  assert(new_count == n_faces);

  /* Free memory */

  BFT_FREE(l_face_cells);
  BFT_FREE(n_t_faces);
  BFT_FREE(t_face_last);
  BFT_FREE(t_cell_index);

  /* Now build final numbering and index */

  /* Build lexical ordering of faces */

  _n_groups = g_id;
  BFT_MALLOC(_group_index, _n_groups*n_i_threads*2, cs_lnum_t);

  _group_index[0] = 0;

  for (g_id = 0; g_id < _n_groups; g_id++) {
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      _group_index[(t_id*_n_groups + g_id)*2] = -1;
      _group_index[(t_id*_n_groups + g_id)*2 + 1] = -1;
    }
  }

#if defined(DEBUG) && !defined(NDEBUG)
  int f_t_id_prev = - 1;
#endif

  for (cs_lnum_t fl_id = 0; fl_id < n_faces; fl_id++) {

    f_id = new_to_old_i[fl_id];

    assert(f_t_id[f_id] > -1);

#if defined(DEBUG) && !defined(NDEBUG)
    assert(f_t_id[f_id] >= f_t_id_prev);
    f_t_id_prev = f_t_id[f_id];
#endif

    int t_id = f_t_id[f_id]%n_i_threads;
    g_id = (f_t_id[f_id] - t_id) / n_i_threads;

    /* Update group index to mark maximum face id */
    _group_index[(t_id*_n_groups + g_id)*2 + 1] = fl_id + 1;

  }

  BFT_FREE(f_t_id);
  BFT_FREE(faces_list);

  /* Finalize group index */

  f_id = 0;
  for (g_id = 0; g_id < _n_groups; g_id++) {
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      _group_index[(t_id*_n_groups + g_id)*2] = f_id;
      f_id = CS_MAX(_group_index[(t_id*_n_groups + g_id)*2+1],
                    f_id);
    }
  }

  for (g_id = 0; g_id < _n_groups; g_id++) {
    for (int t_id = 0; t_id < n_i_threads; t_id++) {
      if (_group_index[(t_id*_n_groups + g_id)*2 + 1] < 0)
        _group_index[(t_id*_n_groups + g_id)*2] = -1;
    }
  }

  *n_groups = _n_groups;
  *n_no_adj_halo_groups = _n_no_adj_halo_groups;
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
 *   new_to_old_i   --> interior faces renumbering array
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
                                      cs_lnum_t     new_to_old_i[],
                                      int          *n_i_groups,
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
                           (const cs_lnum_2_t *)(mesh->i_face_cells),
                           new_to_old_i,
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
 *   new_to_old_b    <-- interior faces renumbering array
 *   b_group_index   --> group/thread index (1 group only)
 *
 * returns:
 *   0 on success, -1 otherwise
  *----------------------------------------------------------------------------*/

static int
_renum_b_faces_no_share_cell_across_thread(cs_mesh_t   *mesh,
                                           int          n_b_threads,
                                           cs_lnum_t    min_subset_size,
                                           cs_lnum_t    new_to_old_b[],
                                           cs_lnum_t  **b_group_index)
{
  int t_id;
  cs_lnum_t ii, subset_size, start_id, end_id;

  int retval = 0;

  /* Initialization */

  BFT_MALLOC(*b_group_index, n_b_threads*2, cs_lnum_t);

  /* Order faces lexicographically (with stable sort) */

  cs_lnum_t *fc_num;
  BFT_MALLOC(fc_num, mesh->n_b_faces*2, cs_lnum_t);

  for (ii = 0; ii < mesh->n_b_faces; ii++) {
    fc_num[ii*2] = mesh->b_face_cells[ii];
    fc_num[ii*2+1] = ii;
  }

  cs_order_lnum_allocated_s(NULL, fc_num, 2, new_to_old_b, mesh->n_b_faces);

  BFT_FREE(fc_num);

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
      cs_lnum_t f_id = new_to_old_b[end_id - 1];
      cs_lnum_t c_id = mesh->b_face_cells[f_id];
      f_id = new_to_old_b[end_id];
      while (mesh->b_face_cells[f_id] == c_id) {
        end_id += 1;
        if (end_id < mesh->n_b_faces)
          f_id = new_to_old_b[end_id];
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
 * Compute renumbering of interior faces for vectorizing.
 *
 * parameters:
 *   mesh         <-> pointer to global mesh structure
 *   vector_size  <-- target size for groups
 *   group_size   <-- target group size
 *   new_to_old_i --> interior faces renumbering array
 *
 * returns:
 *   0 on success, -1 otherwise
  *----------------------------------------------------------------------------*/

static int
_renum_i_faces_for_vectorizing(cs_mesh_t  *mesh,
                               int         vector_size,
                               cs_lnum_t   new_to_old_i[])
{
  int retval = -1;

  const cs_lnum_t n_i_faces = mesh->n_i_faces;
  const cs_lnum_2_t *i_face_cells = (const cs_lnum_2_t *)(mesh->i_face_cells);

  /* Initialize variables to avoid compiler warnings */

  cs_lnum_t swap_id = -1;

  /* Initialization */

  for (cs_lnum_t face_id = 0; face_id < mesh->n_i_faces; face_id++)
    new_to_old_i[face_id] = face_id;

  /* Order interior faces (we choose to place the "remainder" at the end) */

  /* determine remainder and number of complete registers */

  cs_lnum_t irelii = n_i_faces % vector_size;
  cs_lnum_t nregii = n_i_faces / vector_size;

  /* External loop */

  for (int loop_id = 0; loop_id < 100; loop_id++) {

    int mod_prev = 0; /* indicates if elements were exchanged in array new_to_old_i */

    cs_lnum_t iregic = 0; /* Previous register */

    cs_lnum_t block_id = 0;  /* Counter to avoid exchanging all elements
                                of new_to_old_i more than n times */

    /* Loop on elements of new_to_old_i */

    for (cs_lnum_t jj = 0;
         jj < mesh->n_i_faces && block_id > -1;
         jj++) {

      cs_lnum_t last_id, inext;

      /* Current register and position inside it */

      cs_lnum_t iregip = iregic;
      cs_lnum_t jregic = (jj % vector_size) + 1;
      iregic = jj / vector_size + 1;

      /* Test between last_id, start of register, and current position;
         take the worst case between remainder at beginning and end:
         remaninder at beginning */

      if (iregic == 1)
        last_id = 0;
      else if (jregic <= irelii)
        last_id = (iregic-2)*vector_size+irelii;
      else
        last_id = (iregic-1)*vector_size;

      /* Swap starting from inext, start of next register */

      if ((iregic == nregii && jregic > irelii) || (iregic == nregii+1))
        inext = 0;
      else if (jregic > irelii)
        inext = iregic*vector_size+irelii;
      else
        inext = iregic*vector_size;

      if (iregic != iregip) swap_id = inext - 1;

      block_id = 0;

      /* Test with all preceding elements since last_id:
       * swap_id indicates with which element of new_to_old_i we swap
       * mod_prev indicates we modify an already seen element
       * block_id indicates we have seen all elements and we must mix
       * (there is no solution) */

      bool test_all_since_last = true;

      while (test_all_since_last) {

        test_all_since_last = false;
        cs_lnum_t face_id = new_to_old_i[jj];

        for (cs_lnum_t ii = last_id; ii < jj; ii++) {

          cs_lnum_t cn0 = i_face_cells[new_to_old_i[ii]][0];
          cs_lnum_t cn1 = i_face_cells[new_to_old_i[ii]][1];
          cs_lnum_t cr0 = i_face_cells[face_id][0];
          cs_lnum_t cr1 = i_face_cells[face_id][1];

          if (cn0 == cr0 || cn1 == cr1 || cn0 == cr1 || cn1 == cr0) {

            swap_id += 1;

            if (swap_id >= n_i_faces) {
              swap_id = 0;
              block_id = block_id + 1;
            }
            if (swap_id < jj) mod_prev = 1;
            if (block_id >= 2) {
              block_id = -1;
              break;
            }

            cs_lnum_t itmp = new_to_old_i[swap_id];
            new_to_old_i[swap_id] = new_to_old_i[jj];
            new_to_old_i[jj] = itmp;

            test_all_since_last = true;
            break;

          }

        }

      } /* test_all_since_last */;

    } /* loop on jj (faces) */

    /* If we did not touch elements preceding the current one,
       the algorithm has succeeded */

    if (mod_prev == 0 && block_id > -1) {
      retval = 0;
      break;
    }

    /* Shuffle if there is no solution or we looped 10 times */

    if (loop_id < 100 && (((loop_id+1)%10 == 0) || block_id == -1)) {
      for (cs_lnum_t ii = 0; ii < (n_i_faces-4)/2; ii += 2) {
        cs_lnum_t jj = n_i_faces-ii-1;
        cs_lnum_t itmp = new_to_old_i[ii];
        new_to_old_i[ii] = new_to_old_i[jj];
        new_to_old_i[jj] = itmp;
      }
    }

  } /* Loop on loop_id */

  /* Checks */

  if (retval == 0) {

    cs_lnum_t iok = 0;

    cs_lnum_t *order;
    BFT_MALLOC(order, n_i_faces, cs_lnum_t);
    cs_order_lnum_allocated(NULL, new_to_old_i, order, n_i_faces);

    for (cs_lnum_t ii = 0; ii < n_i_faces; ii++) {
      if (new_to_old_i[order[ii]] !=  n_i_faces-ii-1)
        iok -= 1;
    }

    BFT_FREE(order);

    /* Classical test looping on previous faces */

    if (iok == 0) {

      for (cs_lnum_t jj = 0; jj < mesh->n_i_faces; jj++) {

        /* Current register and position inside it */

        cs_lnum_t iregic = jj / vector_size + 1;
        cs_lnum_t jregic = (jj % vector_size) + 1;

        /* Test between last_id, start of register, and current position;
           take the worst case between remainder at beginning and end:
           remaninder at beginning */

        cs_lnum_t last_id;

        if (iregic == 1)
          last_id = 0;
        else if (jregic < irelii)
          last_id = (iregic-2)*vector_size+irelii;
        else
          last_id = (iregic-1)*vector_size;

        /* Test with all preceding elements since last_id */

        for (cs_lnum_t ii = last_id; ii < jj; ii++) {

          cs_lnum_t face_id = new_to_old_i[jj];

          cs_lnum_t cn0 = i_face_cells[new_to_old_i[ii]][0];
          cs_lnum_t cn1 = i_face_cells[new_to_old_i[ii]][1];
          cs_lnum_t cr0 = i_face_cells[face_id][0];
          cs_lnum_t cr1 = i_face_cells[face_id][1];

          if (cn0 == cr0 || cn1 == cr1 || cn0 == cr1 || cn1 == cr0)
            iok -= 1;

        }
      }

    }

    if (iok != 0 && mesh->verbosity > 2) {
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(_("Faces renumbering for vectorization:\n"
                   "====================================\n\n"
                   "%llu errors in interior face renumbering array.\n\n"
                   "Faces are not renumbered, and vectorization of face loops\n"
                   "will not be forced.\n"), (unsigned long long)iok);
      retval = -1;
    }

  }

  /* Return value */

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute renumbering of boundary faces for vectorizing.
 *
 * parameters:
 *   mesh         <-> pointer to global mesh structure
 *   vector_size  <-- target size for groups
 *   group_size   <-- target group size
 *   new_to_old_b --> interior faces renumbering array
 *
 * returns:
 *   0 on success, -1 otherwise
  *----------------------------------------------------------------------------*/

static int
_renum_b_faces_for_vectorizing(cs_mesh_t  *mesh,
                               int         vector_size,
                               cs_lnum_t   new_to_old_b[])
{
  int retval = -1;

  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  cs_lnum_t *b_face_cells = mesh->b_face_cells;

  /* Initialization */

  for (cs_lnum_t face_id = 0; face_id < mesh->n_b_faces; face_id++)
    new_to_old_b[face_id] = face_id;

  /* Order boundary faces */

  /* determine remainder and number of complete registers */

  cs_lnum_t irelib = n_b_faces % vector_size;
  cs_lnum_t nregib = n_b_faces / vector_size;

  /* Maximum number of boundary faces; if < nregib, there is no solution */

  cs_lnum_t *irhss;
  BFT_MALLOC(irhss, n_cells, cs_lnum_t);

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    irhss[cell_id] = 0;

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t ii = b_face_cells[face_id];
    irhss[ii] += 1;
  }

  cs_lnum_t nfamax = 0;
  cs_lnum_t nfanp1 = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    nfamax = CS_MAX(nfamax, irhss[cell_id]);
    if (irhss[cell_id] == nregib+1)
      nfanp1 += 1;
  }

  if (nfamax > nregib+1 || (nfamax == nregib+1 && nfanp1 > irelib)) {
    BFT_FREE(irhss);
    return retval;
  }

  /* Order by decreasing number of cell boundary faces */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t cell_id = b_face_cells[face_id];
    b_face_cells[face_id] += n_cells*irhss[cell_id];
  }

  cs_lnum_t *order;

  BFT_MALLOC(order, n_b_faces, cs_lnum_t);
  cs_order_lnum_allocated(NULL, b_face_cells, order, n_b_faces);

  /* Restore connectivity */
  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++)
    b_face_cells[face_id] = b_face_cells[face_id] % n_cells;

  /* Distribute faces in registers */

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    cs_lnum_t ireg, ilig, ii;
    if (face_id <= irelib*(nregib+1)) {
      ireg = face_id % (nregib+1);
      ilig = face_id / (nregib+1);
      ii = ireg*vector_size+ilig;
    }
    else {
      cs_lnum_t face_id1 = face_id-irelib*(nregib+1);
      ireg = face_id1 % nregib;
      ilig = face_id1 / nregib + irelib;
      ii = ireg*vector_size+ilig;
    }
    new_to_old_b[ii] = order[face_id];
  }

  retval = 0;

  /* Checks */

  cs_lnum_t iok = 0;

  cs_order_lnum_allocated(NULL, new_to_old_b, order, n_b_faces);

  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
  if (new_to_old_b[order[ii]] !=  n_b_faces-ii-1)
    iok -= 1;
  }

  BFT_FREE(order);

  /* Classical test looping on previous faces */

  if (iok == 0) {

    for (cs_lnum_t jj = 0; jj < mesh->n_b_faces; jj++) {

      /* Current register and position inside it */

      cs_lnum_t iregic = jj / vector_size + 1;
      cs_lnum_t jregic = (jj % vector_size) + 1;

      /* Test between last_id, start of register, and current position;
         take the worst case between remainder at beginning and end:
         remaninder at beginning */

      cs_lnum_t last_id;

      if (iregic == 1)
        last_id = 0;
      else if (jregic < irelib)
        last_id = (iregic-2)*vector_size+irelib;
      else
        last_id = (iregic-1)*vector_size;

      /* Test with all preceding elements since last_id */

      for (cs_lnum_t ii = last_id; ii < jj; ii++) {
        cs_lnum_t face_id = new_to_old_b[jj];
        if (b_face_cells[new_to_old_b[ii]] == b_face_cells[face_id])
          iok -= 1;
      }

    }

  }

  if (iok != 0 && mesh->verbosity > 2) {
    /* TODO: add global logging info for rank 0) */
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_("Faces renumbering for vectorization:\n"
                 "====================================\n\n"
                 "%llu errors in boundary face renumbering array.\n\n"
                 "Faces are not renumbered, and vectorization of face loops\n"
                 "will not be forced.\n"), (unsigned long long)iok);
    retval = -1;
  }

  /* Return value */

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

    cs_lnum_t cid0 = i_face_cells[face_id][0];
    cs_lnum_t cid1 = i_face_cells[face_id][1];

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
 * Compute local cell centers.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   cell_center --> cell centers array
 *----------------------------------------------------------------------------*/

static void
_precompute_cell_center(const cs_mesh_t  *mesh,
                        cs_coord_t        cell_center[])
{
  int ft;
  cs_lnum_t i, j;
  cs_lnum_t vtx_id, face_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3], face_center[3];

  cs_lnum_t n_cells = mesh->n_cells;

  const cs_lnum_t  *face_cells = NULL;
  const cs_lnum_t  *face_vtx_idx = NULL;
  const cs_lnum_t  *face_vtx = NULL;
  const cs_real_t  *vtx_coord = mesh->vtx_coord;

  cs_lnum_t n_max_face_vertices = 0;

  cs_real_3_t  *face_vtx_coord = NULL;
  cs_coord_t  *weight = NULL;

  const double surf_epsilon = 1e-24;

  BFT_MALLOC(weight, n_cells, cs_coord_t);

  for (i = 0; i < n_cells; i++) {
    weight[i] = 0.0;
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] = 0.0;
  }

  for (ft = 0; ft < 2; ft++) {

    cs_lnum_t n_faces;

    if (ft == 0) {
      n_faces = mesh->n_i_faces;
      face_cells = (cs_lnum_t *)(mesh->i_face_cells);
      face_vtx_idx = mesh->i_face_vtx_idx;
      face_vtx = mesh->i_face_vtx_lst;
    }
    else {
      n_faces = mesh->n_b_faces;
      face_cells = mesh->b_face_cells;
      face_vtx_idx = mesh->b_face_vtx_idx;
      face_vtx = mesh->b_face_vtx_lst;
    }

    /* Counting and allocation */

    n_max_face_vertices = 0;

    for (face_id = 0; face_id < n_faces; face_id++) {
      n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
      if (n_max_face_vertices <= n_face_vertices)
        n_max_face_vertices = n_face_vertices;
    }

    BFT_MALLOC(face_vtx_coord, n_max_face_vertices, cs_real_3_t);

    /* Loop on each face */

    for (face_id = 0; face_id < n_faces; face_id++) {

      /* Initialization */

      cs_lnum_t tri_id;

      cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
      cs_coord_t face_surface = 0.0;

      n_face_vertices = 0;

      start_id = face_vtx_idx[face_id];
      end_id = face_vtx_idx[face_id + 1];

      /* Define the polygon (P) according to the vertices (Pi) of the face */

      for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

        cs_lnum_t shift = 3 * face_vtx[vtx_id];
        for (i = 0; i < 3; i++)
          face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
        n_face_vertices++;

      }

      /* Compute the center of gravity of the face vertices */

      for (i = 0; i < 3; i++) {
        vtx_cog[i] = 0.0;
        for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
          vtx_cog[i] += face_vtx_coord[vtx_id][i];
        vtx_cog[i] /= n_face_vertices;
      }

      /* Loop on the triangles of the face (defined by an edge of the face
         and its center of gravity) */

      for (i = 0; i < 3; i++) {
        ref_normal[i] = 0.;
        face_center[i] = 0.0;
      }

      for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

        cs_coord_t tri_surface;
        cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

        cs_lnum_t id0 = tri_id;
        cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

        /* Normal for each triangle */

        for (i = 0; i < 3; i++) {
          vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
          vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
        }

        tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
        tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
        tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

        if (tri_id == 0) {
          for (i = 0; i < 3; i++)
            ref_normal[i] = tri_normal[i];
        }

        /* Center of gravity for a triangle */

        for (i = 0; i < 3; i++) {
          tri_center[i] = (  vtx_cog[i]
                           + face_vtx_coord[id0][i]
                           + face_vtx_coord[id1][i]) / 3.0;
        }

        tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                           + tri_normal[1]*tri_normal[1]
                           + tri_normal[2]*tri_normal[2]) * 0.5;

        if ((  tri_normal[0]*ref_normal[0]
             + tri_normal[1]*ref_normal[1]
             + tri_normal[2]*ref_normal[2]) < 0.0)
          tri_surface *= -1.0;

        /* Now compute contribution to face center and surface */

        face_surface += tri_surface;

        for (i = 0; i < 3; i++) {
          face_center[i] += tri_surface * tri_center[i];
          unweighted_center[i] = tri_center[i];
        }

      } /* End of loop  on triangles of the face */

      if (face_surface > surf_epsilon) {
        for (i = 0; i < 3; i++)
          face_center[i] /= face_surface;
      }
      else {
        face_surface = surf_epsilon;
        for (i = 0; i < 3; i++)
          face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
      }

     /* Now contribute to cell centers */

      if (ft == 0) {

        cs_lnum_t cell_id_0 = face_cells[face_id*2];
        cs_lnum_t cell_id_1 = face_cells[face_id*2 + 1];

        if (cell_id_0 < n_cells) {
          for (i = 0; i < 3; i++)
            cell_center[cell_id_0*3 + i] += face_center[i]*face_surface;
          weight[cell_id_0] += face_surface;
        }

        if (cell_id_1 < n_cells) {
          for (i = 0; i < 3; i++)
            cell_center[cell_id_1*3 + i] += face_center[i]*face_surface;
          weight[cell_id_1] += face_surface;
        }

      }
      else {
        cs_lnum_t cell_id = face_cells[face_id];
        for (i = 0; i < 3; i++)
          cell_center[cell_id*3 + i] += face_center[i]*face_surface;
        weight[cell_id] += face_surface;
      }

    } /* End of loop on faces */

    BFT_FREE(face_vtx_coord);

  }

  for (i = 0; i < n_cells; i++) {
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] /= weight[i];
  }

  BFT_FREE(weight);
}

/*----------------------------------------------------------------------------
 * Determine the local extents associated with a set of coordinates
 *
 * parameters:
 *   dim      <-- spatial dimension
 *   n_coords <-- local number of coordinates
 *   coords   <-- entity coordinates; size: n_entities*dim (interlaced)
 *   extents  --> global extents (size: dim*2)
 *---------------------------------------------------------------------------*/

static void
_get_coord_extents(int               dim,
                   size_t            n_coords,
                   const cs_coord_t  coords[],
                   cs_coord_t        extents[])
{
  size_t  i, j;

  /* Get local min/max coordinates */

  for (j = 0; j < (size_t)dim; j++) {
    extents[j]       = DBL_MAX;
    extents[j + dim] = -DBL_MAX;
  }

  for (i = 0; i < n_coords; i++) {
    for (j = 0; j < (size_t)dim; j++) {
      if (coords[i*dim + j] < extents[j])
        extents[j] = coords[i*dim + j];
      if (coords[i*dim + j] > extents[j + dim])
        extents[j + dim] = coords[i*dim + j];
    }
  }
}

/*----------------------------------------------------------------------------
 * Renumber cells based on local Morton encoding.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *----------------------------------------------------------------------------*/

static void
_renum_cells_morton(const cs_mesh_t  *mesh,
                    cs_lnum_t         new_to_old[])
{
  cs_coord_t extents[6];
  fvm_morton_code_t *m_code = NULL;
  cs_lnum_t n_cells = mesh->n_cells;

  if (mesh->cell_numbering->n_no_adj_halo_elts > 0)
    n_cells = mesh->cell_numbering->n_no_adj_halo_elts;

  const int level = sizeof(fvm_morton_int_t)*8 - 1;

  /* Build Morton encoding and order it */

  cs_coord_t *cell_center;

  BFT_MALLOC(cell_center, mesh->n_cells*3, cs_coord_t);

  _precompute_cell_center(mesh, cell_center);

  _get_coord_extents(mesh->dim, mesh->n_cells, cell_center, extents);

  BFT_MALLOC(m_code, mesh->n_cells, fvm_morton_code_t);

  fvm_morton_encode_coords(mesh->dim,
                           level,
                           extents,
                           mesh->n_cells,
                           cell_center,
                           m_code);

  fvm_morton_local_order(n_cells, m_code, new_to_old);

  BFT_FREE(m_code);

  BFT_FREE(cell_center);
}

/*----------------------------------------------------------------------------
 * Renumber cells based on local Hilbert encoding.
 *
 * In the case that 2 entities have a same Hilbert code, their order
 * will be determined by lexicographical ordering of coordinates.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *----------------------------------------------------------------------------*/

static void
_renum_cells_hilbert(const cs_mesh_t  *mesh,
                     cs_lnum_t         new_to_old[])
{
  cs_coord_t extents[6];

  cs_coord_t *cell_center;

  BFT_MALLOC(cell_center, mesh->n_cells*3, cs_coord_t);

  _precompute_cell_center(mesh, cell_center);

  _get_coord_extents(mesh->dim, mesh->n_cells, cell_center, extents);

  fvm_hilbert_local_order_coords(mesh->dim,
                                 extents,
                                 mesh->n_cells,
                                 cell_center,
                                 new_to_old);

  BFT_FREE(cell_center);
}

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

/*----------------------------------------------------------------------------
 * Build Metis graph (cell -> cell connectivity, ignoring ghosts)
 *
 * parameters:
 *   mesh           <-- pointer to mesh structure
 *   n_cells        <-- number of cells considered in graph
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *----------------------------------------------------------------------------*/

static void
_metis_graph(const cs_mesh_t   *mesh,
             idx_t              n_cells,
             idx_t            **cell_idx,
             idx_t            **cell_neighbors)
{
  idx_t i;

  idx_t  *n_neighbors;
  idx_t  *_cell_idx;
  idx_t  *_cell_neighbors;

  const idx_t n_faces = mesh->n_i_faces;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, idx_t);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t id_0 = mesh->i_face_cells[i][0];
    cs_lnum_t id_1 = mesh->i_face_cells[i][1];
    if (id_0 < n_cells && id_1 < n_cells) {
      n_neighbors[id_0] += 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, idx_t);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++)
    _cell_idx[i + 1] = _cell_idx[i] + n_neighbors[i];

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_cells], idx_t);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    cs_lnum_t id_0 = mesh->i_face_cells[i][0];
    cs_lnum_t id_1 = mesh->i_face_cells[i][1];

    if (id_0 < n_cells && id_1 < n_cells) {
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = id_1;
      n_neighbors[id_0] += 1;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = id_0;
      n_neighbors[id_1] += 1;
    }

  }

  BFT_FREE(n_neighbors);

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
}

/*----------------------------------------------------------------------------
 * Compute local partition using METIS
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_part_metis(const cs_mesh_t  *mesh,
            idx_t             n_parts,
            idx_t            *cell_part)
{
  idx_t i;

  idx_t   n_constraints = 1;

  idx_t   edgecut = 0; /* <-- Number of faces on partition */

  idx_t   n_cells = mesh->n_cells;

  idx_t  *cell_idx = NULL, *cell_neighbors = NULL;

  int retval = 0, retcode = METIS_OK;

  /* If we have separated cells with no halo adjacency, only partition
     on those cells */

  if (mesh->cell_numbering->n_no_adj_halo_elts > 0) {
    for (i = mesh->cell_numbering->n_no_adj_halo_elts; i < n_cells; i++)
      cell_part[i] = n_parts - 1;
    n_cells = mesh->cell_numbering->n_no_adj_halo_elts;
  }

  _metis_graph(mesh, n_cells, &cell_idx, &cell_neighbors);

  /* Subpartition */

  if (n_parts < 8) {

    bft_printf(_("\n"
                 " Sub-partitioning cells to %d domains per rank\n"
                 "   (%s).\n"), (int)n_parts, "METIS_PartGraphRecursive");

    retcode
      = METIS_PartGraphRecursive(&n_cells,
                                 &n_constraints,
                                 cell_idx,
                                 cell_neighbors,
                                 NULL,     /* vwgt:   cell weights */
                                 NULL,     /* vsize:  size of the vertices */
                                 NULL,     /* adjwgt: face weights */
                                 &n_parts,
                                 NULL,     /* tpwgts */
                                 NULL,     /* ubvec: load imbalance tolerance */
                                 NULL,     /* options */
                                 &edgecut,
                                 cell_part);
  }

  else {

    bft_printf(_("\n"
                 " Sub-partitioning cells to %d domains per rank\n"
                 "   (%s).\n"), (int)n_parts, "METIS_PartGraphKway");

    retcode
      = METIS_PartGraphKway(&n_cells,
                            &n_constraints,
                            cell_idx,
                            cell_neighbors,
                            NULL,     /* vwgt:   cell weights */
                            NULL,     /* vsize:  size of the vertices */
                            NULL,     /* adjwgt: face weights */
                            &n_parts,
                            NULL,     /* tpwgts */
                            NULL,     /* ubvec: load imbalance tolerance */
                            NULL,     /* options */
                            &edgecut,
                            cell_part);

  }

  BFT_FREE(cell_idx);
  BFT_FREE(cell_neighbors);

  if (retcode != METIS_OK)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Renumber cells based on local Metis partitioning.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_cells_metis_part(const cs_mesh_t  *mesh,
                        cs_lnum_t         new_to_old[])
{
  idx_t *cell_part;
  cs_lnum_t *number;

  int retval = 0;

  BFT_MALLOC(cell_part, mesh->n_cells, idx_t);
  BFT_MALLOC(number, mesh->n_cells * 2, cs_lnum_t);

  retval = _part_metis(mesh, _cs_renumber_n_threads, cell_part);

  for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
    number[i*2] = cell_part[i];
    number[i*2 + 1] = i;
  }

  BFT_FREE(cell_part);

  cs_order_lnum_allocated_s(NULL, number, 2, new_to_old, mesh->n_cells);

  BFT_FREE(number);

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute local ordering using METIS
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_cells_metis_order(const cs_mesh_t  *mesh,
                         cs_lnum_t         new_to_old[])
{
  idx_t   n_cells = mesh->n_cells;
  idx_t  *perm = NULL, *iperm = NULL;

  idx_t  *cell_idx = NULL, *cell_neighbors = NULL;

  int retval = 0, retcode = METIS_OK;

  if (sizeof(idx_t) == sizeof(cs_lnum_t))
    perm = (idx_t *)new_to_old;
  else
    BFT_MALLOC(perm, n_cells, idx_t);

  /* If we have separated cells with no halo adjacency, only order
     those cells */

  if (mesh->cell_numbering->n_no_adj_halo_elts > 0) {
    for (idx_t i = mesh->cell_numbering->n_no_adj_halo_elts; i < n_cells; i++)
      new_to_old[i] = i;
    n_cells = mesh->cell_numbering->n_no_adj_halo_elts;
  }

  BFT_MALLOC(iperm, n_cells, idx_t);

  _metis_graph(mesh, n_cells, &cell_idx, &cell_neighbors);

  bft_printf(_("\n"
               " Ordering cells for fill reduction (METIS_NodeND).\n"));

  retcode = METIS_NodeND(&n_cells,
                         cell_idx,
                         cell_neighbors,
                         NULL,       /* vwgt:   cell weights */
                         NULL,       /* options */
                         perm,
                         iperm);

  BFT_FREE(iperm);

  if (sizeof(idx_t) != sizeof(cs_lnum_t)) {
    for (idx_t i = 0; i < n_cells; i++)
      new_to_old[i] = perm[i];
    BFT_FREE(perm);
  }

  BFT_FREE(cell_idx);
  BFT_FREE(cell_neighbors);

  if (retcode != METIS_OK)
    retval = 1;

  return retval;
}

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

static void
_scotch_sort_shell(SCOTCH_Num  l,
                   SCOTCH_Num  r,
                   SCOTCH_Num  a[])
{
  int i, j, h;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1) ;

  /* Sort array */
  for (; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      SCOTCH_Num v = a[i];

      j = i;
      while ((j >= l+h) && (v < a[j-h])) {
        a[j] = a[j-h];
        j -= h;
      }
      a[j] = v;

    } /* Loop on array elements */

  } /* End of loop on stride */

}

/*----------------------------------------------------------------------------
 * Build SCOTCH graph (cell -> cell connectivity)
 *
 * parameters:
 *   mesh           <-- pointer to mesh structure
 *   n_cells        <-- number of cells considered in graph
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *   mesh_to_graph  --> mesh to graph indirection
 *
 * returns:
 *   number of vertices in graph (cells non in extended halo)
 *----------------------------------------------------------------------------*/

static SCOTCH_Num
_scotch_graph(const cs_mesh_t   *mesh,
              SCOTCH_Num         n_cells,
              SCOTCH_Num       **cell_idx,
              SCOTCH_Num       **cell_neighbors,
              cs_lnum_t        **mesh_to_graph)
{
  SCOTCH_Num  i;
  SCOTCH_Num  start_id, end_id, c_id;

  SCOTCH_Num  *n_neighbors;
  SCOTCH_Num  *_cell_idx;
  SCOTCH_Num  *_cell_neighbors;
  cs_lnum_t   *_mesh_to_graph;

  SCOTCH_Num  n_graph_cells = 0;

  const SCOTCH_Num n_faces = mesh->n_i_faces;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, SCOTCH_Num);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {
    cs_lnum_t id_0 = mesh->i_face_cells[i][0];
    cs_lnum_t id_1 = mesh->i_face_cells[i][1];
    if (id_0 < n_cells && id_1 < n_cells) {
      n_neighbors[id_0] += 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, SCOTCH_Num);
  BFT_MALLOC(_mesh_to_graph, n_cells, cs_lnum_t);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++) {
    if (n_neighbors[i] > 0) {
      _cell_idx[n_graph_cells + 1] = _cell_idx[n_graph_cells] + n_neighbors[i];
      _mesh_to_graph[i] = n_graph_cells;
      n_graph_cells++;
    }
    else
      _mesh_to_graph[i] = -1;
  }

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_graph_cells], SCOTCH_Num);

  for (i = 0; i < n_graph_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    cs_lnum_t id_0 = _mesh_to_graph[mesh->i_face_cells[i][0]];
    cs_lnum_t id_1 = _mesh_to_graph[mesh->i_face_cells[i][1]];

    if (id_0 > -1 && id_1 > -1) {
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = id_1;
      n_neighbors[id_0] += 1;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = id_0;
      n_neighbors[id_1] += 1;
    }

  }

  BFT_FREE(n_neighbors);

  /* Clean graph */

  c_id = 0;
  start_id = _cell_idx[0]; /* also = 0 */
  end_id = 0;

  for (i = 0; i < n_graph_cells; i++) {

    SCOTCH_Num j, n_prev;

    end_id = _cell_idx[i+1];

    if (end_id > start_id) {

      _scotch_sort_shell(start_id, end_id, _cell_neighbors);

      n_prev = _cell_neighbors[start_id];
      _cell_neighbors[c_id] = n_prev;
      c_id += 1;

      for (j = start_id + 1; j < end_id; j++) {
        if (_cell_neighbors[j] != n_prev) {
          n_prev = _cell_neighbors[j];
          _cell_neighbors[c_id] = n_prev;
          c_id += 1;
        }
      }

    }

    start_id = end_id;
    _cell_idx[i+1] = c_id;

  }

  if (c_id < end_id)
    BFT_REALLOC(_cell_neighbors, c_id, SCOTCH_Num);

  /* Set return values */

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
  *mesh_to_graph = _mesh_to_graph;

  return n_graph_cells;
}

/*----------------------------------------------------------------------------
 * Compute local partition using SCOTCH
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_part_scotch(const cs_mesh_t  *mesh,
             SCOTCH_Num        n_parts,
             int              *cell_part)
{
  SCOTCH_Num i;

  SCOTCH_Num   n_cells = mesh->n_cells;
  SCOTCH_Num   n_ext_cells = mesh->n_cells_with_ghosts;
  SCOTCH_Num   n_graph_cells = 0;

  SCOTCH_Graph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;

  SCOTCH_Num  *cell_idx = NULL, *cell_neighbors = NULL;
  cs_lnum_t   *mesh_to_graph = NULL;

  int  retval = 0;

  /* Initialization */

  for (i = 0; i < n_ext_cells; i++)
    cell_part[i] = -1; /* 0 to n for constrained sub-partition */

  n_graph_cells = _scotch_graph(mesh,
                                n_ext_cells,
                                &cell_idx,
                                &cell_neighbors,
                                &mesh_to_graph);

  /* If we have separated cells with no halo adjacency, only partition
     on those cells; distribute others evenly */

  if (mesh->cell_numbering->n_no_adj_halo_elts > 0)
    n_cells = mesh->cell_numbering->n_no_adj_halo_elts;

  for (i = n_cells; i < n_graph_cells; i++) {
    int part_id = ceil((double)i / (double)n_cells) * n_parts - 1;
    cell_part[i] = CS_MIN(part_id, n_parts - 1);
  }

#if SCOTCH_VERSION >= 6
  bft_printf(_("\n"
               " Sub-partitioning cells to %d domains per rank\n"
               "   (%s).\n"), (int)n_parts, "SCOTCH_graphPartFixed");
#else
  bft_printf(_("\n"
               " Sub-partitioning cells to %d domains per rank\n"
               "   (%s).\n"), (int)n_parts, "SCOTCH_graphPart");
#endif

  /* Partition using libScotch */

  SCOTCH_graphInit(&grafdat);

  retval
    = SCOTCH_graphBuild(&grafdat,
                        0,                  /* baseval; 0 to n -1 numbering */
                        n_graph_cells,      /* vertnbr */
                        cell_idx,           /* verttab */
                        NULL,               /* vendtab: verttab + 1 or NULL */
                        NULL,               /* velotab: vertex weights */
                        NULL,               /* vlbltab; vertex labels */
                        cell_idx[n_graph_cells],  /* edgenbr */
                        cell_neighbors,     /* edgetab */
                        NULL);              /* edlotab */

  if (retval == 0) {

    SCOTCH_Num  *graph_part = NULL;

    BFT_MALLOC(graph_part, n_graph_cells, SCOTCH_Num);

    for (i = 0; i < n_ext_cells; i++) {
      if (mesh_to_graph[i] > -1)
        graph_part[mesh_to_graph[i]] = cell_part[i];
    }

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_graphCheck(&grafdat) == 0) {
#if SCOTCH_VERSION >= 6
      retval = SCOTCH_graphPartFixed(&grafdat, n_parts, &stradat, graph_part);
#else
      retval = SCOTCH_graphPart(&grafdat, n_parts, &stradat, graph_part);
#endif
    }

    SCOTCH_stratExit(&stradat);

    if (retval == 0) {
      for (i = 0; i < n_ext_cells; i++) {
        if (mesh_to_graph[i] > -1)
          cell_part[i] = graph_part[mesh_to_graph[i]];
      }
    }

    BFT_FREE(graph_part);

  }

  SCOTCH_graphExit(&grafdat);

  BFT_FREE(mesh_to_graph);

  BFT_FREE(cell_idx);
  BFT_FREE(cell_neighbors);

  return retval;
}

/*----------------------------------------------------------------------------
 * Renumber cells based on local SCOTCH partitioning.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_cells_scotch_part(const cs_mesh_t  *mesh,
                         cs_lnum_t         new_to_old[])
{
  int *cell_part;
  cs_lnum_t *number;

  int retval = 0;

  BFT_MALLOC(cell_part, mesh->n_cells_with_ghosts, int);
  BFT_MALLOC(number, mesh->n_cells * 2, cs_lnum_t);

  retval = _part_scotch(mesh, _cs_renumber_n_threads, cell_part);

  for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
    number[i*2] = cell_part[i];
    number[i*2 + 1] = i;
  }

  BFT_FREE(cell_part);

  cs_order_lnum_allocated_s(NULL, number, 2, new_to_old, mesh->n_cells);

  BFT_FREE(number);

  return retval;
}

/*----------------------------------------------------------------------------
 * Compute local ordering using SCOTCH
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *
 * returns:
 *   0 in case of success, 1 otherwise
 *----------------------------------------------------------------------------*/

static int
_renum_cells_scotch_order(const cs_mesh_t  *mesh,
                          cs_lnum_t         new_to_old[])
{
  SCOTCH_Num   n_cells = mesh->n_cells;
  SCOTCH_Num   n_ext_cells = mesh->n_cells_with_ghosts;
  SCOTCH_Num   n_graph_cells = 0;

  SCOTCH_Graph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;
  SCOTCH_Ordering  order;

  SCOTCH_Num  *peritab = NULL, *listtab = NULL;
  SCOTCH_Num  *cell_idx = NULL, *cell_neighbors = NULL;
  cs_lnum_t   *mesh_to_graph = NULL;

  int  retval = 0;

  /* If we have separated cells with no halo adjacency, only order
     those cells */

  if (mesh->cell_numbering->n_no_adj_halo_elts > 0) {
    for (SCOTCH_Num i = mesh->cell_numbering->n_no_adj_halo_elts;
         i < n_cells;
         i++)
      new_to_old[i] = i;
    n_cells = mesh->cell_numbering->n_no_adj_halo_elts;
  }

  n_graph_cells = _scotch_graph(mesh,
                                n_ext_cells,
                                &cell_idx,
                                &cell_neighbors,
                                &mesh_to_graph);

  BFT_MALLOC(peritab, n_graph_cells, SCOTCH_Num);
  BFT_MALLOC(listtab, n_graph_cells, SCOTCH_Num);

  for (SCOTCH_Num i = 0; i < n_cells; i++)
    listtab[i] = i;

  for (SCOTCH_Num i = n_cells; i < n_graph_cells; i++)
    peritab[i] = i; /* simple precaution, probably not required */

  bft_printf
    (_("\n"
       " Ordering cells for fill reduction (SCOTCH_graphOrderComputeList).\n"));

  /* Order using libScotch */

  SCOTCH_graphInit(&grafdat);

  retval
    = SCOTCH_graphBuild(&grafdat,
                        0,                  /* baseval; 0 to n -1 numbering */
                        n_graph_cells,      /* vertnbr */
                        cell_idx,           /* verttab */
                        NULL,               /* vendtab: verttab + 1 or NULL */
                        NULL,               /* velotab: vertex weights */
                        NULL,               /* vlbltab; vertex labels */
                        cell_idx[n_graph_cells],  /* edgenbr */
                        cell_neighbors,     /* edgetab */
                        NULL);              /* edlotab */

  if (retval == 0) {

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_graphCheck(&grafdat) == 0) {

      retval = SCOTCH_graphOrderInit(&grafdat,
                                     &order,
                                     NULL,   /* permtab */
                                     peritab,
                                     NULL,   /* cblkprt */
                                     NULL,   /* rangtab */
                                     NULL);  /* treetab */

      if (retval == 0) {

        retval = SCOTCH_graphOrderComputeList(&grafdat,
                                              &order,
                                              n_graph_cells,
                                              listtab,
                                              &stradat);  /* treetab */

        if (retval != 0) {
          for (SCOTCH_Num i = 0; i < n_cells; i++)
            peritab[i] = i;
        }

        SCOTCH_graphOrderExit(&grafdat, &order);

      }

    }

  }

  SCOTCH_graphExit(&grafdat);

  /* Free local arrays */

  BFT_FREE(listtab);

  SCOTCH_Num j = 0;
  for (SCOTCH_Num i = 0; i < n_graph_cells; i++) {
    if (peritab[i] < n_cells)
      new_to_old[j++] = peritab[i];
  }
  BFT_FREE(peritab);

  assert(j == n_cells);

  BFT_FREE(mesh_to_graph);

  BFT_FREE(cell_idx);
  BFT_FREE(cell_neighbors);

  return retval;
}

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

/*----------------------------------------------------------------------------
 * Compute local ordering using reverse Cuthill-McKee
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  --> new to old cell renumbering
 *----------------------------------------------------------------------------*/

static void
_renum_cells_rcm(const cs_mesh_t  *mesh,
                 cs_lnum_t         new_to_old[])
{
  if (mesh->n_cells < 1)
    return;

  cs_lnum_t *keys, *order;
  int *cell_class;

  BFT_MALLOC(keys, mesh->n_cells*3, cs_lnum_t);
  BFT_MALLOC(order, mesh->n_cells, cs_lnum_t);

  BFT_MALLOC(cell_class, mesh->n_cells_with_ghosts, int);

  cs_adjacency_t *a
    = _c2c_from_face_cell(mesh->n_cells_with_ghosts,
                          mesh->n_i_faces,
                          (const cs_lnum_t  *)(mesh->i_face_cells));

  cs_lnum_t l_s = 0, l_e = 0;
  bool boot = false, over_constrained = false;

  if (_cells_adjacent_to_halo_last) {

    _classify_cells_by_neighbor(mesh, cell_class);

    int cell_class_max = 2;
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      if (cell_class[i] > cell_class_max)
        cell_class_max = cell_class[i];
    }
    cell_class_max += 1;

    for (cs_lnum_t i = mesh->n_cells; i < mesh->n_cells_with_ghosts; i++)
      cell_class[i] = cell_class_max + 1;

    /* Order halo-adjacent cells in order of matching halo cells */

    for (cs_lnum_t i = mesh->n_cells_with_ghosts-1; i >= mesh->n_cells; i--) {
      for (cs_lnum_t j = a->idx[i+1]-1; j >= a->idx[i]; j--) {
        cs_lnum_t k = a->ids[j];
        if (cell_class[k] < cell_class_max) {
          cell_class[k] = cell_class_max;
          keys[l_e*3] = a->idx[k+1]- a->idx[k];
          keys[l_e*3+1] = -l_e;
          keys[l_e*3+2] = k;
          l_e++;
        }
      }
    }

    for (cs_lnum_t i = mesh->n_cells-1; i >= 0; i--) {
      if (cell_class[i] > 0 && cell_class[i] < cell_class_max) {
         /* halo-adjacent cells marked cell_class_max */
        assert(cell_class[i] == 1);
        cell_class[i] = 2;
        keys[l_e*3] = a->idx[i+1]- a->idx[i];
        keys[l_e*3+1] = -l_e;
        keys[l_e*3+2] = i;
        l_e++;
      }
    }

    if (l_e == 0)
      boot = true;

  }
  else {

    boot = true;
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {

      int  count = 0;
      for (cs_lnum_t j = a->idx[i]; j < a->idx[i+1]; j++)
        if (a->ids[j] >= mesh->n_cells)
          count++;

      cell_class[i] = 0;
      if (count == (a->idx[i+1]-a->idx[i])) {
        cell_class[i] = 2;
        keys[l_e*3  ] = a->idx[i+1]- a->idx[i];
        keys[l_e*3+1] = i;
        keys[l_e*3+2] = i;
        l_e += 1;
      }

    }

    for (cs_lnum_t i = mesh->n_cells; i < mesh->n_cells_with_ghosts; i++)
      cell_class[i] = 2;

  }

  /* If there are starting cells, search for cell of highest degree */

  if (boot) {

    cs_lnum_t  id_min = mesh->n_cells, nn_min = mesh->n_cells;

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      cs_lnum_t nn = a->idx[i+1] - a->idx[i];
      if (nn <= nn_min && cell_class[i] == 0) {
        id_min = i;
        nn_min = nn;
      }
    }

    if (id_min < mesh->n_cells) {
      cell_class[id_min] = 2;
      keys[l_e*3  ] = a->idx[id_min+1]- a->idx[id_min];
      keys[l_e*3+1] = id_min;
      keys[l_e*3+2] = id_min;
      l_e += 1;
    }

  }

  /* Now generate sets */

  cs_lnum_t *rl;
  BFT_MALLOC(rl, mesh->n_cells, cs_lnum_t);

#if defined(DEBUG) && !defined(NDEBUG)
  for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
    rl[i] = -1;
#endif

  int level = -1;

  while (true) {

    cs_order_lnum_allocated_s(NULL,
                              keys,
                              3,
                              order,
                              l_e - l_s);

    for (cs_lnum_t i = l_s; i < l_e; i++) {
      cs_lnum_t j = order[i - l_s];
      rl[i] = keys[j*3+2];
    }

    /* Generate next set */
    if (l_e >= mesh->n_cells)
      break;
    else if (l_e == l_s) {
      /* Disjoint set ? Find next starting cell, similar to boot
         (could be improved, but avoids failure) */
      cs_lnum_t  id_min = mesh->n_cells, nn_min = mesh->n_cells;
      for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
        cs_lnum_t nn = a->idx[i+1] - a->idx[i];
        if (nn <= nn_min && cell_class[i] == 0) {
          id_min = i;
          nn_min = nn;
        }
      }
      if (id_min >= mesh->n_cells) { /* no cell with class 0 */
        over_constrained = true;
        break;
      }
      cell_class[id_min] = level;
      keys[l_e*3  ] = a->idx[id_min+1]- a->idx[id_min];
      keys[l_e*3+1] = id_min;
      keys[l_e*3+2] = id_min;
      rl[l_e] = id_min;
      l_e += 1;
    }

    cs_lnum_t n = 0;
    for (cs_lnum_t l_id = l_s; l_id < l_e; l_id++) {
      cs_lnum_t i = rl[l_id];
      assert(i >= 0 && i < mesh->n_cells);
      for (cs_lnum_t j = a->idx[i+1]-1; j >= a->idx[i]; j--) {
        cs_lnum_t k = a->ids[j];
        if (cell_class[k] == 0) {
          cell_class[k] = level;
          keys[n*3] = a->idx[k+1] - a->idx[k];
          keys[n*3+1] = -n;
          keys[n*3+2] = k;
          n++;
        }
      }
    }

    l_s = l_e;
    l_e = l_s + n;

    level--;

  }

  cs_adjacency_destroy(&a);

  BFT_FREE(keys);
  BFT_FREE(order);
  BFT_FREE(cell_class);

#if defined(DEBUG) && !defined(NDEBUG)
  for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
    new_to_old[i] = -1;
#endif

  if (over_constrained == false) {
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      new_to_old[mesh->n_cells - 1 - i] = rl[i];

#if defined(DEBUG) && !defined(NDEBUG)
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      assert(new_to_old[i] != -1);
    }
#endif
  }
  else { /* if overconstrained, do not renumber */
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      new_to_old[i] = i;
  }

  BFT_FREE(rl);
}

/*----------------------------------------------------------------------------
 * Renumber cells for locality.
 *
 * parameters:
 *   mesh         <-> pointer to global mesh structure
 *   algorithm    <-- algorithm used for renumbering
 *   new_to_old_c <-- cell rnumbering array
 *
 * returns:
 *   0 if renumbering was successful, -1 if failed, 1 if not required
 *----------------------------------------------------------------------------*/

static int
_cells_locality_renumbering(cs_mesh_t                 *mesh,
                            cs_renumber_cells_type_t   algorithm,
                            cs_lnum_t                 *new_to_old_c)
{
  int retval = 0;

  /* Cells renumbering */
  /*-------------------*/

  switch (algorithm) {

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

  case CS_RENUMBER_CELLS_METIS_PART:
    retval = _renum_cells_metis_part(mesh, new_to_old_c);
    break;

  case CS_RENUMBER_CELLS_METIS_ORDER:
    retval = _renum_cells_metis_order(mesh, new_to_old_c);
    break;

#endif

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

  case CS_RENUMBER_CELLS_SCOTCH_PART:
    retval = _renum_cells_scotch_part(mesh, new_to_old_c);
    break;

  case CS_RENUMBER_CELLS_SCOTCH_ORDER:
    retval = _renum_cells_scotch_order(mesh, new_to_old_c);
    break;

#endif

  case CS_RENUMBER_CELLS_MORTON:
    _renum_cells_morton(mesh, new_to_old_c);
    break;

  case CS_RENUMBER_CELLS_HILBERT:
    _renum_cells_hilbert(mesh, new_to_old_c);
    break;

  case CS_RENUMBER_CELLS_RCM:
    _renum_cells_rcm(mesh, new_to_old_c);
    break;

  case CS_RENUMBER_CELLS_NONE:
    retval = 1;
    break;

  default:
    if (algorithm <= CS_RENUMBER_CELLS_NONE)
      bft_printf
        (_("\n"
           " Cell prenumbering of type: %s\n"
           "   not supported in this build.\n"),
         _(_cell_renum_name[algorithm]));
    else
      bft_printf
        (_("\n"
           " Cell prenumbering of type: %d\n"
           "   not supported in this build.\n"),
         (int)algorithm);

    retval = -1;
    break;
  }

  if (retval != 0) {
    for (cs_lnum_t ii = 0; ii < mesh->n_cells; ii++)
      new_to_old_c[ii] = ii;
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Ensure only cells not neighboring ghost cells are renumbered.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  <-> new to old cell renumbering
 *----------------------------------------------------------------------------*/

static void
_renum_only_no_adj_halo_cells(const cs_mesh_t  *mesh,
                              cs_lnum_t         new_to_old[])
{
  if (mesh->cell_numbering->n_no_adj_halo_elts > 0) {
    cs_lnum_t j = 0;
    cs_lnum_t n_cells = mesh->cell_numbering->n_no_adj_halo_elts;
    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      cs_lnum_t old_id = new_to_old[i];
      if (old_id < n_cells)
        new_to_old[j++] = old_id;
    }
    assert(j == n_cells);
    for (cs_lnum_t i = n_cells; i < mesh->n_cells; i++)
      new_to_old[i] = i;
  }
}

/*----------------------------------------------------------------------------
 * Ensure only cells not neighboring ghost cells are renumbered.
 *
 * parameters:
 *   mesh        <-- pointer to mesh structure
 *   new_to_old  <-> new to old cell renumbering
 *----------------------------------------------------------------------------*/

static void
_renum_adj_halo_cells_last(const cs_mesh_t  *mesh,
                           cs_lnum_t         new_to_old[])
{
  cs_lnum_t n_i_cells = 0;

  assert(mesh->cell_numbering != NULL);

  if (_cells_adjacent_to_halo_last) {

    cs_lnum_t *number;
    int *cell_class;

    BFT_MALLOC(number, mesh->n_cells*2, cs_lnum_t);

    BFT_MALLOC(cell_class, mesh->n_cells, int);

    _classify_cells_by_neighbor(mesh, cell_class);

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++) {
      number[i*2] = cell_class[i];
      if (cell_class[i] == 0)
        n_i_cells++;
    }

    BFT_FREE(cell_class);

    for (cs_lnum_t i = 0; i < mesh->n_cells; i++)
      number[new_to_old[i]*2 + 1] = i;

    cs_order_lnum_allocated_s(NULL, number, 2, new_to_old, mesh->n_cells);

    BFT_FREE(number);
  }

  if (n_i_cells > 0) {
    cs_numbering_t  *numbering = mesh->cell_numbering;
    numbering->n_no_adj_halo_elts = n_i_cells;
    numbering->n_threads = 1;
    numbering->n_groups = 2;
    BFT_REALLOC(numbering->group_index, 4, cs_lnum_t);
    numbering->group_index[0] = 0;
    numbering->group_index[1] = n_i_cells;
    numbering->group_index[2] = n_i_cells;
    numbering->group_index[3] = mesh->n_cells;
  }
}

/*----------------------------------------------------------------------------
 * Renumbering of vertices based on cell adjacency.
 *
 * parameters:
 *   mesh   <-> pointer to global mesh structure
 *   n2o_v  <-> new to old numbering for vertices
 *----------------------------------------------------------------------------*/

static void
_renumber_vertices_by_cell_adjacency(cs_mesh_t  *mesh,
                                     cs_lnum_t  *n2o_v)
{
  const cs_lnum_t n_vertices = mesh->n_vertices;

  /* Order vertices such that the couple (f_id, v_id) is scanned in an
     increasing way */

  cs_lnum_t  *v_couples = NULL;
  BFT_MALLOC(v_couples, 2*n_vertices, cs_lnum_t);

  /* Set the initial values */

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_couples[2*i] = -1, v_couples[2*i+1] = n2o_v[i];

  /* Build the cell --> faces connectivity */

  cs_adjacency_t  *c2f = cs_mesh_adjacency_c2f(mesh, 1);

  /* Map the boundary face --> vertices connectivity */

  cs_adjacency_t  *bf2v
    = cs_adjacency_create_from_i_arrays(mesh->n_b_faces,
                                        mesh->b_face_vtx_idx,
                                        mesh->b_face_vtx_lst,
                                        NULL);

  /* Map the interior face --> vertices connectivity */

  cs_adjacency_t  *if2v
    = cs_adjacency_create_from_i_arrays(mesh->n_i_faces,
                                        mesh->i_face_vtx_idx,
                                        mesh->i_face_vtx_lst,
                                        NULL);

  /* Loop on cells to build the v_couples */

  for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++) {

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->ids[j];

      if (f_id < mesh->n_i_faces) { /* Interior face */

        for (cs_lnum_t jj = if2v->idx[f_id]; jj < if2v->idx[f_id+1]; jj++) {
          if (v_couples[2*if2v->ids[jj]] == -1)
            v_couples[2*if2v->ids[jj]] = c_id;
        }

      }
      else { /* Border face */

        const cs_lnum_t  bf_id = f_id - mesh->n_i_faces;
        for (cs_lnum_t jj = bf2v->idx[bf_id]; jj < bf2v->idx[bf_id+1]; jj++) {
          if (v_couples[2*bf2v->ids[jj]] == -1)
            v_couples[2*bf2v->ids[jj]] = c_id;
        }

      }

    } /* Loop on cell faces */

  } /* Loop on cells */

  cs_order_lnum_allocated_s(NULL,
                            v_couples,
                            2,
                            n2o_v,
                            n_vertices);

  /* Free temporary array */

  BFT_FREE(v_couples);
  cs_adjacency_destroy(&c2f);
  cs_adjacency_destroy(&if2v);
  cs_adjacency_destroy(&bf2v);
}

/*----------------------------------------------------------------------------
 * Renumbering of vertices based on face adjacency.
 *
 * parameters:
 *   mesh   <-> pointer to global mesh structure
 *   n2o_v  <-> new to old numbering for vertices
 *----------------------------------------------------------------------------*/

static void
_renumber_vertices_by_face_adjacency(cs_mesh_t  *mesh,
                                     cs_lnum_t  *n2o_v)
{
  const cs_lnum_t n_vertices = mesh->n_vertices;

  /* Order vertices such that the couple (f_id, v_id) is scanned in an
     increasing way */

  cs_lnum_t  *v_couples = NULL;
  BFT_MALLOC(v_couples, 2*n_vertices, cs_lnum_t);

  /* Set the initial values */

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_couples[2*i] = -1, v_couples[2*i+1] = n2o_v[i];

  /* Loop on interior faces */

  for (cs_lnum_t i = 0; i < mesh->n_i_faces; i++) {

    for (cs_lnum_t j = mesh->i_face_vtx_idx[i]; j < mesh->i_face_vtx_idx[i+1];
         j++) {

      const cs_lnum_t  v_id = mesh->i_face_vtx_lst[j];

      if (v_couples[2*v_id] < 0)
        v_couples[2*v_id] = i;
      else {
        if (v_couples[2*v_id] > i)
          v_couples[2*v_id] = i;
      }

    }

  } /* Loop on interior faces */

  /* Loop on border faces */

  for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {

    const cs_lnum_t  bf_id = i + mesh->n_i_faces;

    for (cs_lnum_t j = mesh->b_face_vtx_idx[i]; j < mesh->b_face_vtx_idx[i+1];
         j++) {

      const cs_lnum_t  v_id = mesh->b_face_vtx_lst[j];

      if (v_couples[2*v_id] < 0)
        v_couples[2*v_id] = bf_id;
      else {
        if (v_couples[2*v_id] > bf_id)
          v_couples[2*v_id] = bf_id;
      }

    }

  } /* Loop on border faces */

  cs_order_lnum_allocated_s(NULL,
                            v_couples,
                            2,
                            n2o_v,
                            n_vertices);

  /* Free temporary array */

  BFT_FREE(v_couples);

}

/*----------------------------------------------------------------------------
 * Renumber cells for locality and possible computation/communication
 * overlap.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_cells(cs_mesh_t  *mesh)
{
  cs_lnum_t  *new_to_old_c = NULL;
  int retval = 0;
  int halo_order_stage = 0;

  if (mesh->cell_numbering != NULL)
    cs_numbering_destroy(&(mesh->cell_numbering));

  mesh->cell_numbering = cs_numbering_create_default(mesh->n_cells);

  BFT_MALLOC(new_to_old_c, mesh->n_cells_with_ghosts, cs_lnum_t);

  /* When do we reorder cells by adjacent halo ?
     0: never, 1: after pre-ordering, 2: after ordering;
     Scotch can be made aware of graph vertices with fixed
     partition or ordered as "ghosts"; Metis and space-filling
     curves can not. */

  if (_cells_adjacent_to_halo_last) {
    switch (_cells_algorithm[1]) {
    case CS_RENUMBER_CELLS_SCOTCH_PART:
    case CS_RENUMBER_CELLS_SCOTCH_ORDER:
      halo_order_stage = 1;
      break;
    default:
      halo_order_stage = 2;
    }
  }

  /* Initial numbering, and optional classification by neighbor */

  if (_cells_algorithm[0] != CS_RENUMBER_CELLS_NONE) {

    retval = _cells_locality_renumbering(mesh,
                                         _cells_algorithm[0],
                                         new_to_old_c);

    if (retval != 0 && _cells_algorithm[0] != CS_RENUMBER_CELLS_NONE)
      bft_printf
        (_("\n Cell prenumbering (%s) failed.\n"),
         _cell_renum_name[_cells_algorithm[0]]);

    if (halo_order_stage == 1)
      _renum_adj_halo_cells_last(mesh, new_to_old_c);

    if (retval == 0 || halo_order_stage == 1)
      _cs_renumber_update_cells(mesh, new_to_old_c);
  }

  /* Last stage: numbering for locality */

  retval = _cells_locality_renumbering(mesh,
                                       _cells_algorithm[1],
                                       new_to_old_c);

  if (halo_order_stage == 2)
    _renum_adj_halo_cells_last(mesh, new_to_old_c);
  else if (halo_order_stage == 1)
    _renum_only_no_adj_halo_cells(mesh, new_to_old_c);

  /* Now update mesh connectivity */
  /*------------------------------*/

  if (retval == 0 || halo_order_stage > 0)
    _cs_renumber_update_cells(mesh, new_to_old_c);

  else if (retval < 0)
    bft_printf(_("\n Cell renumbering (%s) failed.\n"),
               _cell_renum_name[_cells_algorithm[1]]);

  if (_cells_algorithm[1] != CS_RENUMBER_CELLS_NONE)
    bft_printf
      ("\n ----------------------------------------------------------\n");

  if (mesh->verbosity > 0)
    cs_numbering_log_info(CS_LOG_DEFAULT,
                          _("cells"),
                          mesh->cell_numbering);

  /* Now free remaining array */

  BFT_FREE(new_to_old_c);
}

/*----------------------------------------------------------------------------
 * Try to apply renumbering of interior faces for multiple threads.
 *
 * Relation to graph edge coloring:
 * No graph vertex (cell) is incident to 2 edges (faces) of the same color.
 * A thread pool may thus be built, with 1 thread per color.
 * Groups may then be built, containing only cells of a given color.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_i_faces(cs_mesh_t  *mesh)
{
  int  n_i_groups = 1, n_i_no_adj_halo_groups = 0;
  cs_lnum_t  max_group_size = 1014;       /* Default */
  cs_lnum_t  ii;
  cs_lnum_t  *new_to_old_i = NULL;
  cs_lnum_t  *i_group_index = NULL;

  int  n_i_threads = _cs_renumber_n_threads;

  cs_numbering_type_t numbering_type = CS_NUMBERING_DEFAULT;

  int retval = 0;

  /* Note: group indexes for n_threads and n_groups are defined as follows:
   *  group_index <-- group_index[thread_id*group_id*2 + 2*group_id] and
   *                  group_index[thread_id*group_id*2 + 2*group_id +1]
   *                  define the start and end ids (+1) for entities in a
   *                  given group and thread (size: n_groups *2 * n_threads) */

  /* Allocate Work array */

  BFT_MALLOC(new_to_old_i, mesh->n_i_faces, cs_lnum_t);

  /* Initialize renumbering array */

  for (ii = 0; ii < mesh->n_i_faces; ii++)
    new_to_old_i[ii] = ii;

  /* Interior faces renumbering */
  /*----------------------------*/

  /* We always apply renumbering by cell adjacency using a lexicographical
     ordering; when ordering by lowest id first, the resulting order should be
     similar to the "natural" ordering of faces generated during mesh
     conversion, from a nodal to a "native" (face->cells) representation.

     For some algorithms (such as multipass), this lexicographical ordering is
     included, so for better efficiency, it is not called twice. */

  /* Adjust block size depending on the number of faces and threads */

  switch (_i_faces_algorithm) {
  case CS_RENUMBER_I_FACES_BLOCK:
    numbering_type = CS_NUMBERING_THREADS;
    _renumber_i_faces_by_cell_adjacency(mesh);
    retval = _renum_i_faces_no_share_cell_in_block(mesh,
                                                   n_i_threads,
                                                   max_group_size,
                                                   new_to_old_i,
                                                   &n_i_groups,
                                                   &i_group_index);
    break;

  case CS_RENUMBER_I_FACES_MULTIPASS:
    numbering_type = CS_NUMBERING_THREADS;
    retval = _renum_face_multipass(mesh,
                                   n_i_threads,
                                   new_to_old_i,
                                   &n_i_groups,
                                   &n_i_no_adj_halo_groups,
                                   &i_group_index);
    break;

  case CS_RENUMBER_I_FACES_SIMD:
    numbering_type = CS_NUMBERING_VECTORIZE;
    _renumber_i_faces_by_cell_adjacency(mesh);
    retval = _renum_i_faces_for_vectorizing(mesh,
                                            _cs_renumber_vector_size,
                                            new_to_old_i);
    break;

  case CS_RENUMBER_I_FACES_NONE:
  default:
    _renumber_i_faces_by_cell_adjacency(mesh);
    retval = -1;
    break;
  }

  /* Update mesh if needed */

  if (retval != 0) {
    n_i_groups = 1;
    n_i_threads = 1;
  }
  else
    _cs_renumber_update_i_faces(mesh, new_to_old_i);

  /* Transfer interior face numbering information to mesh */

  if (numbering_type == CS_NUMBERING_THREADS) {
    if (retval != 0) {
      BFT_REALLOC(i_group_index, 2, cs_lnum_t);
      i_group_index[0] = 0;
      i_group_index[1] = mesh->n_i_faces;
    }
    mesh->i_face_numbering = cs_numbering_create_threaded(n_i_threads,
                                                          n_i_groups,
                                                          i_group_index);
    mesh->i_face_numbering->n_no_adj_halo_groups = n_i_no_adj_halo_groups;
    if (n_i_threads == 1)
      mesh->i_face_numbering->type = CS_NUMBERING_DEFAULT;
  }
  else if (numbering_type == CS_NUMBERING_VECTORIZE && retval == 0) {
    mesh->i_face_numbering
      = cs_numbering_create_vectorized(mesh->n_i_faces,
                                       _cs_renumber_vector_size);
  }
  else
    mesh->i_face_numbering
      = cs_numbering_create_default(mesh->n_i_faces);

  if (mesh->verbosity > 0)
    cs_numbering_log_info(CS_LOG_DEFAULT,
                          _("interior faces"),
                          mesh->i_face_numbering);

  /* Free memory */

  BFT_FREE(i_group_index);
  BFT_FREE(new_to_old_i);
}

/*----------------------------------------------------------------------------
 * Try to apply renumbering of boundary faces.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_b_faces(cs_mesh_t  *mesh)
{
  cs_lnum_t  ii;
  cs_lnum_t  *new_to_old_b = NULL;
  cs_lnum_t  *b_group_index = NULL;

  int  n_b_threads = _cs_renumber_n_threads;

  cs_numbering_type_t numbering_type = CS_NUMBERING_DEFAULT;

  int retval = 0;

  /* Note: group indexes for n_threads and 1 group are defined as follows:
   *  group_index <-- group_index[thread_id] and
   *                  group_index[thread_id +1]
   *                  define the start and end ids (+1) for entities in a
   *                  thread (size: 2 * n_threads) */

  /* Allocate Work array */

  BFT_MALLOC(new_to_old_b, mesh->n_b_faces, cs_lnum_t);

  /* Initialize renumbering array */

  for (ii = 0; ii < mesh->n_b_faces; ii++)
    new_to_old_b[ii] = ii;

  /* Boundary faces renumbering */
  /*----------------------------*/

  /* We always apply renumbering by cell adjacency; the resulting order should
     be similar to the "natural" ordering of faces generated during mesh
     conversion, from a nodal to a "native" (face->cells) representation.

     For some algorithms (such as no_share_cell_across_thread),
     this lexicographical ordering is included, so for better efficiency,
     it is not called twice. */

  switch (_b_faces_algorithm) {
  case CS_RENUMBER_B_FACES_THREAD:
    numbering_type = CS_NUMBERING_THREADS;
    retval = _renum_b_faces_no_share_cell_across_thread(mesh,
                                                        n_b_threads,
                                                        _min_b_subset_size,
                                                        new_to_old_b,
                                                        &b_group_index);
    break;

  case CS_RENUMBER_B_FACES_SIMD:
    numbering_type = CS_NUMBERING_VECTORIZE;
    _renumber_b_faces_by_cell_adjacency(mesh);
    retval = _renum_b_faces_for_vectorizing(mesh,
                                            _cs_renumber_vector_size,
                                            new_to_old_b);
    break;

  case CS_RENUMBER_B_FACES_NONE:
  default:
    _renumber_b_faces_by_cell_adjacency(mesh);
    retval = -1;
    break;
  }

  /* Update mesh if needed */
  /*-----------------------*/

  if (retval != 0) {
    n_b_threads = 1;
  }
  else
    _cs_renumber_update_b_faces(mesh, new_to_old_b);

  /* Transfer boundary face numbering information to mesh */

  if (numbering_type == CS_NUMBERING_THREADS) {
    if (retval != 0) {
      BFT_REALLOC(b_group_index, 2, cs_lnum_t);
      b_group_index[0] = 0;
      b_group_index[1] = mesh->n_b_faces;
    }
    mesh->b_face_numbering = cs_numbering_create_threaded(n_b_threads,
                                                          1,
                                                          b_group_index);
    if (n_b_threads == 1)
      mesh->b_face_numbering->type = CS_NUMBERING_DEFAULT;
  }
  else if (numbering_type == CS_NUMBERING_VECTORIZE && retval == 0) {
    mesh->b_face_numbering
      = cs_numbering_create_vectorized(mesh->n_b_faces,
                                       _cs_renumber_vector_size);
  }
  else
    mesh->b_face_numbering
      = cs_numbering_create_default(mesh->n_b_faces);

  mesh->b_face_numbering->n_no_adj_halo_groups = 0;

  if (mesh->verbosity > 0)
    cs_numbering_log_info(CS_LOG_DEFAULT,
                          _("boundary faces"),
                          mesh->b_face_numbering);

  /* Free memory */

  BFT_FREE(b_group_index);
  BFT_FREE(new_to_old_b);
}

/*----------------------------------------------------------------------------
 * Renumber vertices for locality and possible computation/communication
 * overlap.
 *
 * parameters:
 *   mesh <-> pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_vertices(cs_mesh_t  *mesh)
{
  if (_vertices_algorithm == CS_RENUMBER_VERTICES_NONE)
    return;

  cs_lnum_t  *n2o_v = NULL;

  if (mesh->vtx_numbering != NULL)
    cs_numbering_destroy(&(mesh->vtx_numbering));

  mesh->vtx_numbering = cs_numbering_create_default(mesh->n_vertices);

  BFT_MALLOC(n2o_v, mesh->n_vertices, cs_lnum_t);

  for (cs_lnum_t ii = 0; ii < mesh->n_vertices; ii++) n2o_v[ii] = ii;

  /* Vertices renumbering */
  /*----------------------*/

  switch(_vertices_algorithm) {

  case CS_RENUMBER_VERTICES_BY_CELL_ADJ:
    _renumber_vertices_by_cell_adjacency(mesh, n2o_v);
    break;

  case CS_RENUMBER_VERTICES_BY_FACE_ADJ:
    _renumber_vertices_by_face_adjacency(mesh, n2o_v);
    break;

  default:
    break; /* Nothing to do */
  }

  /* Update mesh if needed */
  /*-----------------------*/

  {
    /* Check numbering is non trivial */

    cs_lnum_t v_id = 0;
    while (v_id < mesh->n_vertices) {
      if (n2o_v[v_id] != v_id)
        break;
      else
        v_id++;
    }

    /* Update connectivity */
    if (v_id < mesh->n_vertices)
      _cs_renumber_update_vertices(mesh, n2o_v);

  }

  if (mesh->verbosity > 0)
    cs_numbering_log_info(CS_LOG_DEFAULT,
                          _("vertices"),
                          mesh->vtx_numbering);

  /* Now free remaining array */

  BFT_FREE(n2o_v);
}

/*----------------------------------------------------------------------------
 * Test local operations related to renumbering for interior faces.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_i_test(cs_mesh_t  *mesh)
{
  if (mesh == NULL)
    return;

  if (mesh->i_face_numbering != NULL) {

    cs_gnum_t face_errors = 0;
    cs_lnum_t *accumulator = NULL;

    if (mesh->i_face_numbering->type == CS_NUMBERING_THREADS) {

      cs_lnum_t counter = 0;

      const int n_threads = mesh->i_face_numbering->n_threads;
      const int n_groups = mesh->i_face_numbering->n_groups;
      const cs_lnum_t *group_index = mesh->i_face_numbering->group_index;

      if (mesh->verbosity > 1)
        bft_printf
          (_("\n"
             "Checking interior faces renumbering...\n"));

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (cs_lnum_t c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        accumulator[c_id_0] = 0;

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {
          for (cs_lnum_t f_id = group_index[(t_id*n_groups + g_id)*2];
               f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               f_id++) {
            cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
            cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
            accumulator[c_id_0] += 1;
            accumulator[c_id_1] += 1;
          }
        }

      }

      for (cs_lnum_t c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        counter += accumulator[c_id_0];

      face_errors = mesh->n_i_faces*2 - counter;

      /* Additional serial test */

      if (face_errors == 0) {

        for (int g_id = 0; g_id < n_groups; g_id++) {

          bool adj_halo = false;

          for (cs_lnum_t c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
            accumulator[c_id_0] = -1;

          for (int t_id = 0; t_id < n_threads; t_id++) {
            for (cs_lnum_t f_id = group_index[(t_id*n_groups + g_id)*2];
                 f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
                 f_id++) {
              cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
              cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
              if (   (accumulator[c_id_0] > -1 && accumulator[c_id_0] != t_id)
                  || (accumulator[c_id_1] > -1 && accumulator[c_id_1] != t_id)) {
                face_errors += 1;
                if (mesh->verbosity > 3)
                  bft_printf("f_id %ld (%ld %ld) g %d t %d\n",
                             (long)f_id, (long)c_id_0, (long)c_id_1, g_id, t_id);
              }
              accumulator[c_id_0] = t_id;
              accumulator[c_id_1] = t_id;
              if (c_id_0 >= mesh->n_cells || c_id_1 >= mesh->n_cells)
                adj_halo = true;
            }
          }

          if (adj_halo) {
            mesh->i_face_numbering->n_no_adj_halo_groups
              = CS_MIN(mesh->i_face_numbering->n_no_adj_halo_groups, g_id+1);
          }

        }

      }

      BFT_FREE(accumulator);
    }

    else if (mesh->i_face_numbering->type == CS_NUMBERING_VECTORIZE) {

      cs_lnum_t f_id, c_id_0, c_id_1;

      cs_lnum_t counter = 0;

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        accumulator[c_id_0] = 0;

#     if defined(HAVE_OPENMP_SIMD)
#       pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#     else
#       pragma dir nodep
#       pragma GCC ivdep
#     endif
      for (f_id = 0; f_id < mesh->n_i_faces; f_id++) {
        c_id_0 = mesh->i_face_cells[f_id][0];
        c_id_1 = mesh->i_face_cells[f_id][1];
        accumulator[c_id_0] += 1;
        accumulator[c_id_1] += 1;
      }

      for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
        counter += accumulator[c_id_0];

      face_errors = mesh->n_i_faces*2 - counter;

      /* Additional serial test */

      if (face_errors == 0) {

        const cs_lnum_t vector_size = mesh->i_face_numbering->vector_size;

        for (c_id_0 = 0; c_id_0 < mesh->n_cells_with_ghosts; c_id_0++)
          accumulator[c_id_0] = -1;

        for (f_id = 0; f_id < mesh->n_i_faces; f_id++) {
          cs_lnum_t block_id = f_id / vector_size;
          c_id_0 = mesh->i_face_cells[f_id][0];
          c_id_1 = mesh->i_face_cells[f_id][1];
          if (   accumulator[c_id_0] == block_id
              || accumulator[c_id_1] == block_id) {
            face_errors += 1;
            if (mesh->verbosity > 3)
              bft_printf("f_id %ld (%ld %ld) b %d\n",
                         (long)f_id, (long)c_id_0, (long)c_id_1, (int)block_id);
          }
          accumulator[c_id_0] = block_id;
          accumulator[c_id_1] = block_id;
        }

      }

      BFT_FREE(accumulator);
    }

    if (mesh->verbosity > 0) {

      cs_parall_counter(&face_errors, 1);

      if (face_errors != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("%llu conflicts detected using interior faces renumbering."),
                  (unsigned long long)face_errors);

    }
  }
}

/*----------------------------------------------------------------------------
 * Test local operations related to renumbering for boundary faces.
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_b_test(cs_mesh_t  *mesh)
{
  if (mesh == NULL)
    return;

  /* Check for boundary faces */
  /*--------------------------*/

  if (mesh->b_face_numbering != NULL) {

    cs_gnum_t face_errors = 0;
    cs_lnum_t *accumulator = NULL;

    if (mesh->verbosity > 1)
      bft_printf
        (_("\n"
           "Checking boundary faces renumbering...\n"));

    if (mesh->b_face_numbering->type == CS_NUMBERING_THREADS) {

      cs_lnum_t counter = 0;

      const int n_threads = mesh->b_face_numbering->n_threads;
      const int n_groups = mesh->b_face_numbering->n_groups;
      const cs_lnum_t *group_index = mesh->b_face_numbering->group_index;

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (cs_lnum_t c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
        accumulator[c_id] = 0;

      for (int g_id = 0; g_id < n_groups; g_id++) {

#       pragma omp parallel for
        for (int t_id = 0; t_id < n_threads; t_id++) {
          for (cs_lnum_t f_id = group_index[(t_id*n_groups + g_id)*2];
               f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
               f_id++) {
            cs_lnum_t c_id = mesh->b_face_cells[f_id];
            accumulator[c_id] += 1;
          }
        }

      }

      for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
        counter += accumulator[c_id];

      face_errors = mesh->n_b_faces - counter;

      /* Additional serial test */

      if (face_errors == 0) {

        for (int g_id = 0; g_id < n_groups; g_id++) {

          for (cs_lnum_t c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
            accumulator[c_id] = -1;

          for (int t_id = 0; t_id < n_threads; t_id++) {
            for (cs_lnum_t f_id = group_index[(t_id*n_groups + g_id)*2];
                 f_id < group_index[(t_id*n_groups + g_id)*2 + 1];
                 f_id++) {
              cs_lnum_t c_id = mesh->b_face_cells[f_id];
              if (accumulator[c_id] > -1 && accumulator[c_id] != t_id)
                face_errors += 1;
              accumulator[c_id] = t_id;
            }
          }

        }

      }

      BFT_FREE(accumulator);
    }

    if (mesh->b_face_numbering->type == CS_NUMBERING_VECTORIZE) {

      cs_lnum_t counter = 0;

      BFT_MALLOC(accumulator, mesh->n_cells_with_ghosts, cs_lnum_t);

      for (cs_lnum_t c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
        accumulator[c_id] = 0;

#       if defined(HAVE_OPENMP_SIMD)
#         pragma omp simd safelen(CS_NUMBERING_SIMD_SIZE)
#       else
#         pragma dir nodep
#         pragma GCC ivdep
#       endif
        for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
          cs_lnum_t c_id = mesh->b_face_cells[f_id];
          accumulator[c_id] += 1;
        }

      for (cs_lnum_t c_id = 0; c_id < mesh->n_cells; c_id++)
        counter += accumulator[c_id];

      face_errors = mesh->n_b_faces - counter;

      /* Additional serial test */

      if (face_errors == 0) {

        const cs_lnum_t vector_size = mesh->b_face_numbering->vector_size;

        for (cs_lnum_t c_id = 0; c_id < mesh->n_cells_with_ghosts; c_id++)
          accumulator[c_id] = -1;

        for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
          cs_lnum_t block_id = f_id / vector_size;
          cs_lnum_t c_id = mesh->b_face_cells[f_id];
          if (accumulator[c_id] == block_id)
            face_errors += 1;
          if (mesh->verbosity > 3)
            bft_printf("f_id %ld (%ld) b %d\n",
                       (long)f_id, (long)c_id, (int)block_id);
          accumulator[c_id] = block_id;
        }

      }

      BFT_FREE(accumulator);
    }

    cs_parall_counter(&face_errors, 1);

    if (face_errors != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("%llu conflicts detected using boundary faces renumbering."),
                (unsigned long long)face_errors);
  }
}

/*----------------------------------------------------------------------------
 * Renumber mesh elements for vectorization or OpenMP depending on code
 * options and target machine.
 *
 * parameters:
 *   mesh  <->  Pointer to global mesh structure
 *----------------------------------------------------------------------------*/

static void
_renumber_mesh(cs_mesh_t  *mesh)
{
  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

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
      _renumber_i_test(mesh);
      _renumber_b_test(mesh);
      return;
    }
#endif

  }

  /* Cell pre-numbering may be ignored if not useful for the
     chosen renumbering algorithm; scotch algorithms are
     made aware of the halo and halo-adjacent cells, and
     graph-based partitioning algorithms in general may use
     pre-numbering to order cells inside sub-partitions. */

  if (_cells_algorithm[0] != CS_RENUMBER_CELLS_NONE) {

    switch (_cells_algorithm[1]) {
    case CS_RENUMBER_CELLS_METIS_PART:
    case CS_RENUMBER_CELLS_SCOTCH_PART:
    case CS_RENUMBER_CELLS_RCM:
      break;
    case CS_RENUMBER_CELLS_SCOTCH_ORDER:
      if (!_cells_adjacent_to_halo_last)
        _cells_algorithm[0] = CS_RENUMBER_CELLS_NONE;
      break;
    default:
      _cells_algorithm[0] = CS_RENUMBER_CELLS_NONE;
    }

    if (   _cells_algorithm[0] == CS_RENUMBER_CELLS_NONE
        && mesh->verbosity > 0)
      bft_printf
        (_("\n"
           "   Cells pre-renumbering deactivated, as it is not useful\n"
           "   for the current numbering algorithm.\n"));
  }

  if (mesh->verbosity > 0) {

    int c_halo_adj_last = (_cells_adjacent_to_halo_last) ? 1 : 0;
    int i_halo_adj_last = (_i_faces_adjacent_to_halo_last) ? 1 : 0;
    int hi = (_i_faces_base_ordering == CS_RENUMBER_ADJACENT_LOW) ? 0 : 1;
    const char *no_yes[] = {N_("no"), N_("yes")};
    const char *low_high[] = {N_("lowest id first"), N_("highest id first")};

    bft_printf
      (_("\n"
         "   renumbering for cells:\n"
         "     pre-numbering:                       %s\n"
         "     cells adjacent to ghost cells last:  %s\n"
         "     numbering:                           %s\n"),
       _(_cell_renum_name[_cells_algorithm[0]]),
       _(no_yes[c_halo_adj_last]),
       _(_cell_renum_name[_cells_algorithm[1]]));

    bft_printf
      (_("\n"
         "   renumbering for interior faces:\n"
         "     cell adjacency pre-ordering:         %s\n"
         "     faces adjacent to ghost cells last:  %s\n"
         "     numbering:                           %s\n"),
        _(low_high[hi]),_(no_yes[i_halo_adj_last]),
       _(_i_face_renum_name[_i_faces_algorithm]));

    bft_printf
      (_("\n"
         "   renumbering for boundary faces:\n"
         "     numbering:                           %s\n"),
       _(_b_face_renum_name[_b_faces_algorithm]));

    bft_printf
      (_("\n"
         "   renumbering for vertices:\n"
         "     numbering:                           %s\n"),
       _(_vertices_renum_name[_vertices_algorithm]));

  }

  /* Renumber cells first */

  _renumber_cells(mesh);

  /* Renumber faces afterwards */

  _renumber_i_faces(mesh);
  _renumber_b_faces(mesh);

  /* Renumber vertices afterwards */

  _renumber_vertices(mesh);

  if (mesh->verbosity > 0)
    bft_printf
      ("\n ----------------------------------------------------------\n");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
  if (_cs_renumber_n_threads < 1) {
    if (n_threads > 1) {
      _i_faces_algorithm = CS_RENUMBER_I_FACES_MULTIPASS;
      _b_faces_algorithm = CS_RENUMBER_B_FACES_THREAD;
    }
  }

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
  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

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
 * \brief Select the algorithm for mesh renumbering.
 *
 * \param[in]  halo_adjacent_cells_last  if true, cells adjacent to ghost cells
 *                                       will be placed last
 *                                       (after pre-numbering)
 * \param[in]  halo_adjacent_faces_last  if true, interior faces adjacent to
 *                                       ghost cells will be placed last
 *                                       (after pre-numbering)
 * \param[in]  i_faces_base_ordering     pre-ordering of interior faces by
 *                                       lowest or highest adjacent cell id
 * \param[in]  cells_pre_numbering       algorithm for cells pre-numbering
 * \param[in]  cells_numbering           algorithm for cells numbering
 * \param[in]  i_faces_numbering         algorithm for interior faces numbering
 * \param[in]  b_faces_numbering         algorithm for boundary faces numbering
 * \param[in]  vertices_numbering        algorithm for vertices numbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_set_algorithm(bool                         halo_adjacent_cells_last,
                          bool                         halo_adjacent_faces_last,
                          cs_renumber_ordering_t       i_faces_base_ordering,
                          cs_renumber_cells_type_t     cells_pre_numbering,
                          cs_renumber_cells_type_t     cells_numbering,
                          cs_renumber_i_faces_type_t   i_faces_numbering,
                          cs_renumber_b_faces_type_t   b_faces_numbering,
                          cs_renumber_vertices_type_t  vertices_numbering)
{
  _cells_adjacent_to_halo_last = halo_adjacent_cells_last;
  _i_faces_adjacent_to_halo_last = halo_adjacent_faces_last;
  _i_faces_base_ordering = i_faces_base_ordering;

  _cells_algorithm[0] = cells_pre_numbering;
  _cells_algorithm[1] = cells_numbering;
  _i_faces_algorithm = i_faces_numbering;
  _b_faces_algorithm = b_faces_numbering;
  _vertices_algorithm = vertices_numbering;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the algorithms for mesh renumbering.
 *
 * Any argument may be passed NULL if this option is not queried.
 *
 * \param[out]  halo_adjacent_cells_last  if true, cells adjacent to ghost cells
 *                                        will be placed last
 *                                        (after pre-numbering)
 * \param[out]  halo_adjacent_faces_last  if true, interior faces adjacent to
 *                                        ghost cells will be placed last
 *                                        (after pre-numbering)
 * \param[out]  i_faces_base_ordering     pre-ordering of interior faces by
 *                                        lowest or highest adjacent cell id
 * \param[out]  cells_pre_numbering       algorithm for cells pre-numbering
 * \param[out]  cells_numbering           algorithm for cells numbering
 * \param[out]  i_faces_numbering         algorithm for interior faces numbering
 * \param[out]  b_faces_numbering         algorithm for boundary faces numbering
 * \param[out]  vertices_numbering        algorithm for vertices numbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_get_algorithm(bool                        *halo_adjacent_cells_last,
                          bool                        *halo_adjacent_faces_last,
                          cs_renumber_ordering_t      *i_faces_base_ordering,
                          cs_renumber_cells_type_t    *cells_pre_numbering,
                          cs_renumber_cells_type_t    *cells_numbering,
                          cs_renumber_i_faces_type_t  *i_faces_numbering,
                          cs_renumber_b_faces_type_t  *b_faces_numbering,
                          cs_renumber_vertices_type_t *vertices_numbering)
{
  if (halo_adjacent_cells_last != NULL)
    *halo_adjacent_cells_last = _cells_adjacent_to_halo_last;
  if (halo_adjacent_faces_last != NULL)
    *halo_adjacent_faces_last = _i_faces_adjacent_to_halo_last;
  if (i_faces_base_ordering != NULL)
    *i_faces_base_ordering = _i_faces_base_ordering;

  if (cells_pre_numbering != NULL)
    *cells_pre_numbering = _cells_algorithm[0];
  if (cells_numbering != NULL)
    *cells_numbering = _cells_algorithm[1];
  if (i_faces_numbering != NULL)
    *i_faces_numbering = _i_faces_algorithm;
  if (b_faces_numbering != NULL)
    *b_faces_numbering = _b_faces_algorithm;
  if (vertices_numbering != NULL)
    *vertices_numbering = _vertices_algorithm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber mesh elements for vectorization or threading depending on
 * code options and target machine.
 *
 * Renumbering cells may also allow improving locality (and favor faces
 * renumbering).
 * It is also possible to place cells connected to ghost cells last,
 * which may be useful to enable computation/communication overlap.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t  *mesh)
{
  bft_printf(_("\n Renumbering mesh:\n"));
  bft_printf_flush();

  _renumber_mesh(mesh);

  if (mesh->cell_numbering == NULL)
    mesh->cell_numbering = cs_numbering_create_default(mesh->n_cells);
  if (mesh->i_face_numbering == NULL)
    mesh->i_face_numbering = cs_numbering_create_default(mesh->n_i_faces);
  if (mesh->b_face_numbering == NULL)
    mesh->b_face_numbering = cs_numbering_create_default(mesh->n_b_faces);
  if (mesh->vtx_numbering == NULL)
    mesh->vtx_numbering = cs_numbering_create_default(mesh->n_vertices);

  _renumber_i_test(mesh);
  _renumber_b_test(mesh);

  if (mesh->verbosity > 0)
    _log_bandwidth_info(mesh, _("volume mesh"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber cells depending on code options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_cells(cs_mesh_t  *mesh)
{
  if (mesh->cell_numbering != NULL)
    cs_numbering_destroy(&(mesh->cell_numbering));

  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

  p = getenv("CS_RENUMBER");

  if (p != NULL) {
    if (strcmp(p, "off") == 0 || strcmp(p, "IBM") == 0) {
      if (mesh->cell_numbering == NULL)
        mesh->cell_numbering = cs_numbering_create_default(mesh->n_cells);
      return;
    }
  }

  /* Apply renumbering */

  _renumber_cells(mesh);

  if (mesh->verbosity > 0)
    bft_printf
      ("\n ----------------------------------------------------------\n");

  if (mesh->cell_numbering == NULL)
    mesh->cell_numbering = cs_numbering_create_default(mesh->n_cells);

  if (mesh->verbosity > 0)
    _log_bandwidth_info(mesh, _("volume mesh"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber interior faces for vectorization or threading depending on
 * code options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_i_faces(cs_mesh_t  *mesh)
{
  if (mesh->i_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->i_face_numbering));

  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

  p = getenv("CS_RENUMBER");

  if (p != NULL) {
    if (strcmp(p, "off") == 0 || strcmp(p, "IBM") == 0) {
      if (mesh->i_face_numbering == NULL)
        mesh->i_face_numbering = cs_numbering_create_default(mesh->n_i_faces);
      return;
    }
  }

  /* Apply renumbering */

  _renumber_i_faces(mesh);

  if (mesh->verbosity > 0)
    bft_printf
      ("\n ----------------------------------------------------------\n");

  if (mesh->i_face_numbering == NULL)
    mesh->i_face_numbering = cs_numbering_create_default(mesh->n_i_faces);

  _renumber_i_test(mesh);
}

/*----------------------------------------------------------------------------
 * Renumber interior faces by global number.
 *
 * This effectively resets the interior faces to their initial numbering.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_i_faces_by_gnum(cs_mesh_t  *mesh)
{
  if (mesh->i_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->i_face_numbering));

  if (mesh->global_i_face_num != NULL) {

    cs_lnum_t *new_to_old_i = cs_order_gnum(NULL,
                                            mesh->global_i_face_num,
                                            mesh->n_i_faces);

    _cs_renumber_update_i_faces(mesh, new_to_old_i);

    mesh->i_face_numbering
      = cs_numbering_create_default(mesh->n_i_faces);

    BFT_FREE(new_to_old_i);

    if (mesh->n_domains < 2)
      BFT_FREE(mesh->global_i_face_num);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber boundary faces for vectorization or threading depending on
 * code options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_b_faces(cs_mesh_t  *mesh)
{
  if (mesh->b_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->b_face_numbering));

  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

  p = getenv("CS_RENUMBER");

  if (p != NULL) {
    if (strcmp(p, "off") == 0 || strcmp(p, "IBM") == 0) {
      if (mesh->b_face_numbering == NULL)
        mesh->b_face_numbering = cs_numbering_create_default(mesh->n_b_faces);
      return;
    }
  }

  /* Apply renumbering */

  _renumber_b_faces(mesh);

  if (mesh->verbosity > 0)
    bft_printf
      ("\n ----------------------------------------------------------\n");

  if (mesh->b_face_numbering == NULL)
    mesh->b_face_numbering = cs_numbering_create_default(mesh->n_b_faces);

  _renumber_b_test(mesh);
}

/*----------------------------------------------------------------------------
 * Renumber boundary faces by global number.
 *
 * This effectively resets the boundary faces to their initial numbering.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_b_faces_by_gnum(cs_mesh_t  *mesh)
{
  if (mesh->b_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->b_face_numbering));

  if (mesh->global_b_face_num != NULL) {

    cs_lnum_t *new_to_old_b = cs_order_gnum(NULL,
                                            mesh->global_b_face_num,
                                            mesh->n_b_faces);

    _cs_renumber_update_b_faces(mesh, new_to_old_b);

    mesh->b_face_numbering
      = cs_numbering_create_default(mesh->n_b_faces);

    BFT_FREE(new_to_old_b);

    if (mesh->n_domains < 2)
      BFT_FREE(mesh->global_b_face_num);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber boundary faces such that selected faces appear last
 *        and will be ignored.
 *
 * Those faces will appear last, and the local number of boundary faces set
 * to the number of remaining faces; The mesh's n_b_faces_all and
 * n_g_b_faces_all allows accessing the full boundary faces list.
 *
 * \param[in, out]  mesh      pointer to global mesh structure
 * \param[in]       n_faces   number of selected faces
 * \param[in]       face_ids  number of selected faces
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_b_faces_select_ignore(cs_mesh_t        *mesh,
                                  cs_lnum_t         n_faces,
                                  const cs_lnum_t   face_ids[])
{
  cs_lnum_t  *_face_ids = NULL;

  /* Try to ensure number of total boundary faces with symmetry is updated */

  if (mesh->n_b_faces_all < mesh->n_b_faces) {
    mesh->n_g_b_faces_all = mesh->n_g_b_faces;
    mesh->n_b_faces_all = mesh->n_b_faces;
  }

  mesh->n_g_b_faces = mesh->n_g_b_faces_all;
  mesh->n_b_faces = mesh->n_b_faces_all;

  if (n_faces < 1)
    return;

  if (mesh->b_face_numbering != NULL)
    cs_numbering_destroy(&(mesh->b_face_numbering));

  /* Use selection flag to avoir list reordering issues */

  char *sel_flag = NULL;

  /* Revert to initial numbering */

  if (mesh->global_b_face_num != NULL) {
    const cs_lnum_t n = mesh->n_b_faces;
    const cs_lnum_t n_s = n_faces;
    const cs_lnum_t *s_id = face_ids;

    cs_lnum_t *new_to_old_b = cs_order_gnum(NULL,
                                            mesh->global_b_face_num,
                                            mesh->n_b_faces);
    _cs_renumber_update_b_faces(mesh, new_to_old_b);

    BFT_MALLOC(sel_flag, mesh->n_b_faces, char);
    for (cs_lnum_t i = 0; i < n; i++)
      sel_flag[i] = 0;
    for (cs_lnum_t i = 0; i < n_s; i++)
      sel_flag[s_id[i]] = 1;

    BFT_MALLOC(_face_ids, n_s, cs_lnum_t);
    cs_lnum_t j = 0;
    for (cs_lnum_t i = 0; i < mesh->n_b_faces; i++) {
      if (sel_flag[new_to_old_b[i]] != 0) {
        _face_ids[j] = i;
        j++;
      }
    }
    assert(j == n_s);

    BFT_FREE(new_to_old_b);

    if (mesh->n_domains < 2)
      BFT_FREE(mesh->global_b_face_num);
  }

  cs_lnum_t *new_to_old = NULL;
  {
    const cs_lnum_t n = mesh->n_b_faces;
    const cs_lnum_t n_s = n_faces;
    const cs_lnum_t *s_id = (_face_ids != NULL) ? _face_ids : face_ids;

    BFT_MALLOC(new_to_old, n, cs_lnum_t);
    if (sel_flag == NULL)
      BFT_MALLOC(sel_flag, n, char);

    for (cs_lnum_t i = 0; i < n; i++)
      sel_flag[i] = 0;

    for (cs_lnum_t i = 0; i < n_s; i++)
      sel_flag[s_id[i]] = 1;

    cs_lnum_t k = 0, l = n - n_s;

    for (cs_lnum_t i = 0; i < n; i++) {
      if (sel_flag[i] == 0) {
        new_to_old[k] = i;
        k++;
      }
      else {
        new_to_old[l] = i;
        l++;
      }
    }
  }

  BFT_FREE(_face_ids);
  BFT_FREE(sel_flag);

  _cs_renumber_update_b_faces(mesh, new_to_old);
  BFT_FREE(new_to_old);

  /* Also modify global boundary faces so that the remaining faces can
     be used in parallel operators requiring them in a consistent manner */

  mesh->n_b_faces = mesh->n_b_faces_all - n_faces;

  if (mesh->n_domains > 1 || mesh->global_b_face_num != NULL) {

    cs_lnum_t n_b_faces = mesh->n_b_faces;
    cs_lnum_t n_b_faces_ext = mesh->n_b_faces_all - mesh->n_b_faces;

    fvm_io_num_t *n_io_num
      = fvm_io_num_create_from_select(NULL,
                                      mesh->global_b_face_num,
                                      n_b_faces,
                                      0);
    fvm_io_num_t *n_io_num_end
      = fvm_io_num_create_from_select(NULL,
                                      mesh->global_b_face_num + mesh->n_b_faces,
                                      n_b_faces_ext,
                                      0);

    const cs_gnum_t *b_gnum = fvm_io_num_get_global_num(n_io_num);
    const cs_gnum_t *b_gnum_end = fvm_io_num_get_global_num(n_io_num_end);

    cs_gnum_t n_g_b_faces = fvm_io_num_get_global_count(n_io_num);

    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      mesh->global_b_face_num[i] = b_gnum[i];

    for (cs_lnum_t i = 0; i < n_b_faces_ext; i++)
      mesh->global_b_face_num[n_b_faces + i] = b_gnum_end[i] + n_g_b_faces;

    n_io_num = fvm_io_num_destroy(n_io_num);
    n_io_num_end = fvm_io_num_destroy(n_io_num_end);

    mesh->n_g_b_faces = n_g_b_faces;

  }
  else
    mesh->n_g_b_faces = mesh->n_b_faces;

  mesh->b_face_numbering
    = cs_numbering_create_default(mesh->n_b_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber vertices depending on code options and target machine.
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_vertices(cs_mesh_t  *mesh)
{
  if (mesh->vtx_numbering != NULL)
    cs_numbering_destroy(&(mesh->vtx_numbering));

  const char *p = NULL;

  /* Initialization */

  if (_cs_renumber_n_threads < 1)
    cs_renumber_set_n_threads(cs_glob_n_threads);

  p = getenv("CS_RENUMBER");

  if (p != NULL) {
    if (strcmp(p, "off") == 0) {
      if (mesh->vtx_numbering == NULL)
        mesh->vtx_numbering = cs_numbering_create_default(mesh->n_vertices);
      return;
    }
  }

  /* Apply renumbering */

  _renumber_vertices(mesh);

  if (mesh->verbosity > 0)
    bft_printf
      ("\n ----------------------------------------------------------\n");

  if (mesh->vtx_numbering == NULL)
    mesh->vtx_numbering = cs_numbering_create_default(mesh->n_vertices);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
