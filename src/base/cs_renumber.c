/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Optional mesh renumbering
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_prototypes.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_renumber.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Redistribute scalar values in case of renubering
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   n_elts   <--  Number of elements
 *   renum    <--  Pointer to renumbering array (new -> old, 1 to n)
 *   val      <->  Pointer to array of vector values
 *   tmp_val  <->  Working array (size n_elts)
 *----------------------------------------------------------------------------*/

static void
_update_elt_scalar(cs_int_t         n_elts,
                   const cs_int_t  *renum,
                   cs_real_t       *val,
                   cs_real_t       *tmp_val)
{
  cs_int_t  elt_id, tmp_elt_id;

  for (elt_id = 0; elt_id < n_elts; elt_id++) {
    tmp_elt_id = renum[elt_id] - 1;
    tmp_val[elt_id] = val[tmp_elt_id];
  }

  memcpy(val, tmp_val, n_elts*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------
 * Redistribute vector values in case of renubering
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   n_elts   <--  Number of elements
 *   renum    <--  Pointer to renumbering array (new -> old, 1 to n)
 *   val      <->  Pointer to array of vector values
 *   tmp_val  <->  Working array (size n_elts)
 *----------------------------------------------------------------------------*/

static void
_update_elt_vector(cs_int_t         n_elts,
                   const cs_int_t  *renum,
                   cs_real_t       *val,
                   cs_real_t       *tmp_val)
{
  int  dim_id;
  cs_int_t  elt_id, tmp_elt_id;

  for (elt_id = 0; elt_id < n_elts; elt_id++) {
    for (dim_id = 0; dim_id < 3; dim_id++) {
      tmp_elt_id = renum[elt_id] - 1;
      tmp_val[elt_id*3 + dim_id] = val[tmp_elt_id*3 + dim_id];
    }
  }

  memcpy(val, tmp_val, n_elts*3*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------
 * Update cell quantities in case they were built before renumbering.
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   mesh            <-> Pointer to global mesh structure
 *   mesh_quantities <-> Pointer to global mesh quantities structure
 *   renum           <-- Cells renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_update_cell_quantities(cs_mesh_t             *mesh,
                        cs_mesh_quantities_t  *mesh_quantities,
                        const cs_int_t        *renum)
{
  cs_real_t  *tmp_val = NULL;

  if (mesh == NULL && mesh_quantities == NULL)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(tmp_val, mesh->n_cells*3, cs_real_t);

  /* Cells */

  if (renum != NULL) {

    if (mesh_quantities->cell_cen != NULL)
      _update_elt_vector(mesh->n_cells,
                         renum,
                         mesh_quantities->cell_cen,
                         tmp_val);

    if (mesh_quantities->cell_vol != NULL)
      _update_elt_scalar(mesh->n_cells,
                         renum,
                         mesh_quantities->cell_vol,
                         tmp_val);

  }

  /* Free Work array */

  BFT_FREE(tmp_val);
}

/*----------------------------------------------------------------------------
 * Update face quantities in case they were build before renumbering.
 *
 * This is the case when the mesh is read in the obsolete 'slc' format.
 *
 * parameters:
 *   mesh            <-> Pointer to global mesh structure
 *   mesh_quantities <-> Pointer to global mesh quantities structure
 *   renum_i         <-- Interior faces renumbering array (new -> old)
 *   renum_b         <-- Boundary faces renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_update_face_quantities(cs_mesh_t             *mesh,
                        cs_mesh_quantities_t  *mesh_quantities,
                        const cs_int_t        *renum_i,
                        const cs_int_t        *renum_b)
{
  cs_real_t  *tmp_val = NULL;
  cs_int_t  n_val_max = CS_MAX(mesh->n_i_faces, mesh->n_b_faces) * 3;

  if (mesh == NULL && mesh_quantities == NULL)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(tmp_val, n_val_max, cs_real_t);

  /* Interior faces */

  if (renum_i != NULL) {

    if (mesh_quantities->i_face_normal != NULL)
      _update_elt_vector(mesh->n_i_faces,
                         renum_i,
                         mesh_quantities->i_face_normal,
                         tmp_val);

    if (mesh_quantities->i_face_cog != NULL)
      _update_elt_vector(mesh->n_i_faces,
                         renum_i,
                         mesh_quantities->i_face_cog,
                         tmp_val);

  }

  /* Boundary Faces */

  if (renum_b != NULL) {

    if (mesh_quantities->b_face_normal != NULL)
      _update_elt_vector(mesh->n_b_faces,
                         renum_b,
                         mesh_quantities->b_face_normal,
                         tmp_val);

    if (mesh_quantities->b_face_cog != NULL)
      _update_elt_vector(mesh->n_b_faces,
                         renum_b,
                         mesh_quantities->b_face_cog,
                         tmp_val);

  }

  /* Free Work arrays */

  BFT_FREE(tmp_val);
}

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
_update_global_num(size_t              n_elts,
                   const fvm_lnum_t    init_num[],
                   fvm_gnum_t        **global_num)
{
  size_t i;
  fvm_gnum_t *_global_num = *global_num;

  if (_global_num == NULL) {

      BFT_MALLOC(_global_num, n_elts, fvm_gnum_t);

      for (i = 0; i < n_elts; i++)
        _global_num[i] = init_num[i];

      *global_num = _global_num;
  }

  else {

    fvm_gnum_t *tmp_global;

    BFT_MALLOC(tmp_global, n_elts, fvm_gnum_t);
    memcpy(tmp_global, _global_num, n_elts*sizeof(fvm_gnum_t));

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
 *   mesh_quantities <-> Pointer to global mesh quantities structure
 *   renum           <-- Cells renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_cells(cs_mesh_t             *mesh,
                          cs_mesh_quantities_t  *mesh_quantities,
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

  /* Update cell quantities if present */

  _update_cell_quantities(mesh, mesh_quantities, renum);

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
 *   mesh_quantities <-> Pointer to global mesh quantities structure
 *   renum_i         <-- Interior faces renumbering array (new -> old)
 *   renum_b         <-- Boundary faces renumbering array (new -> old)
 *----------------------------------------------------------------------------*/

static void
_cs_renumber_update_faces(cs_mesh_t             *mesh,
                          cs_mesh_quantities_t  *mesh_quantities,
                          const cs_int_t        *renum_i,
                          const cs_int_t        *renum_b)
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

  /* Update associated face quantities if present */

  _update_face_quantities(mesh,
                          mesh_quantities,
                          renum_i,
                          renum_b);

  /* Update parent face numbers for post-processing meshes
     that may already have been built; Post-processing meshes
     built after renumbering will have correct parent numbers */

  cs_post_renum_faces(renum_i, renum_b);
}

/*----------------------------------------------------------------------------
 * Try to apply renumbering of faces and cells for multiple threads.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

static void
_renumber_for_threads(cs_mesh_t             *mesh,
                      cs_mesh_quantities_t  *mesh_quantities)
{
  int  update_c = 0, update_fi = 0, update_fb = 0;
  int  n_i_groups = 1, n_b_groups = 1;
  int  n_i_threads = 1, n_b_threads = 1;
  cs_int_t  *inumc = NULL, *inumfi = NULL, *inumfb = NULL;
  fvm_lnum_t *i_group_index = NULL, *b_group_index = NULL;

  /* Allocate Work arrays */

  BFT_MALLOC(inumc, mesh->n_cells, cs_int_t);
  BFT_MALLOC(inumfi, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(inumfb, mesh->n_b_faces, cs_int_t);

  /* Try renumbering */

  {
    cs_int_t ii;

    for (ii = 0; ii < mesh->n_cells; ii++)
      inumc[ii] = ii+1;

    for (ii = 0; ii < mesh->n_i_faces; ii++)
      inumfi[ii] = ii+1;

    for (ii = 0; ii < mesh->n_b_faces; ii++)
      inumfb[ii] = ii+1;
  }

  /* Only keep non-trivial renumbering arrays */

  if (update_c == 0)
    BFT_FREE(inumc);

  if (update_fi == 0)
    BFT_FREE(inumfi);

  if (update_fb == 0)
    BFT_FREE(inumfb);

  /* Now update mesh */
  /*-----------------*/

  if (inumfi != NULL || inumfb != NULL)
    _cs_renumber_update_faces(mesh,
                              mesh_quantities,
                              inumfi,
                              inumfb);

  if (inumc != NULL)
    _cs_renumber_update_cells(mesh,
                              mesh_quantities,
                              inumc);

  /* Add numbering info to mesh */
  /*----------------------------*/

  /*
   *  n_threads   <-- number of threads
   *  n_groups    <-- number of groups
   *  group_index <-- group_index[thread_id*group_id*2 + 2*group_id] and
   *                  group_index[thread_id*group_id*2 + 2*group_id +1]
   *                  define the tart and end ids (+1) for entities in a
   *                  given group and thread (size: n_groups *2 * n_threads)
   */

  if (n_i_groups >= 1) {
    mesh->i_face_numbering = cs_numbering_create_threaded(n_i_threads,
                                                          n_i_groups,
                                                          i_group_index);
    BFT_FREE(i_group_index);
  }

  if (n_b_groups >= 1) {
    mesh->b_face_numbering = cs_numbering_create_threaded(n_b_threads,
                                                          n_b_groups,
                                                          b_group_index);
    BFT_FREE(b_group_index);
  }

  /* Now free remaining arrays */

  BFT_FREE(inumfi);
  BFT_FREE(inumfb);
  BFT_FREE(inumc);
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
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

static void
_renumber_for_vectorizing(cs_mesh_t             *mesh,
                          cs_mesh_quantities_t  *mesh_quantities)
{
  int _ivect[2] = {0, 0};
  cs_int_t   ivecti = 0, ivectb = 0;
  cs_int_t  *inumfi = NULL, *inumfb = NULL;
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

  const cs_int_t vector_size = 1024; /* Use register 16 */

#elif defined(SX) && defined(_SX) /* For NEC SX series */

  const cs_int_t vector_size = 256; /* At least for NEC SX-9 */

#else

  const cs_int_t vector_size = 1; /* Non-vector machines */

#endif

  /* Nothing to do if vector size = 1 */

  if (vector_size == 1)
    return;

  /* Allocate Work arrays */

  BFT_MALLOC(inumfi, mesh->n_i_faces, cs_int_t);
  BFT_MALLOC(inumfb, mesh->n_b_faces, cs_int_t);
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
                           inumfi,
                           inumfb,
                           iworkf,
                           ismbs);

  /* Free Work arrays */

  BFT_FREE(ismbs);
  BFT_FREE(iworkf);

  /* Update mesh */

  if (ivecti > 0 || ivectb > 0) {

    cs_int_t   *_inumfi = NULL;
    cs_int_t   *_inumfb = NULL;

    if (ivecti > 0)
      _inumfi = inumfi;
    if (ivectb > 0)
      _inumfb = inumfb;

    _cs_renumber_update_faces(mesh,
                              mesh_quantities,
                              _inumfi,
                              _inumfb);

  }

  /* Free final work arrays */

  BFT_FREE(inumfb);
  BFT_FREE(inumfi);

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
    mesh->i_face_numbering = cs_numbering_create_vectorized(vector_size);
  if (ivectb > 0)
    mesh->b_face_numbering = cs_numbering_create_vectorized(vector_size);

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
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Renumber mesh elements for vectorization or OpenMP depending on code
 * options and target machine.
 *
 * Currently, only the legacy vectorizing renumbering is handled.
 *
 * parameters:
 *   mesh            <->  Pointer to global mesh structure
 *   mesh_quantities <->  Pointer to global mesh quantities structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t             *mesh,
                 cs_mesh_quantities_t  *mesh_quantities)
{
  /* Try vectorizing first */

  _renumber_for_vectorizing(mesh, mesh_quantities);

  /* Then try OpenMP if available */

  if (cs_glob_n_threads > 1)
    _renumber_for_threads(mesh, mesh_quantities);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
