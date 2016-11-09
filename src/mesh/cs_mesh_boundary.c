/*============================================================================
 * Insert thin walls into the mesh.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_mesh_builder.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_boundary.h"

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

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compare local numbers (qsort function).
 *
 * parameters:
 *   x <-> pointer to first value
 *   y <-> pointer to second value
 *
 * returns:
 *   1 if x > y, -1 if x < y, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_compare_nums(const void *x,
              const void *y)
{
  int retval = 0;
  cs_lnum_t v_diff = (*((const cs_lnum_t *)x) - *((const cs_lnum_t *)y));
  if (v_diff > 0)
    retval = 1;
  else if (v_diff < 0)
    retval = -1;
  return retval;
}

/*----------------------------------------------------------------------------
 * Build a renumbering for interior faces
 *
 * parameters:
 *   n_i_faces <-- number of interior faces
 *   list_size <-- size of selected (to be removed) list
 *   list      <-- ordered ids of of faces to remove  (0 to n-1)
 *
 * returns:
 *   interior face renumbering array (to by freed by caller)
 *----------------------------------------------------------------------------*/

static cs_lnum_t *
_i_list_renum(cs_lnum_t         n_i_faces,
              cs_lnum_t         list_size,
              const cs_lnum_t  *list)
{
  cs_lnum_t  *renum;
  BFT_MALLOC(renum, n_i_faces, cs_lnum_t);

  cs_lnum_t i_empty = 0;
  cs_lnum_t cur_id = 0;

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    if (i_empty < list_size) {

      if (face_id != list[i_empty])
        renum[face_id] = cur_id++;
      else {
        renum[face_id] = -1;
        i_empty++;
      }

    }
    else
      renum[face_id] = cur_id++;
  }

  return renum;
}

/*----------------------------------------------------------------------------
 * Remove internal faces of i_face_vtx_idx and i_face_vtx_lst from one index
 * (realloc has to be done after).
 *
 * parameters:
 *   i_face_vtx_idx  <-> interior faces -> vertices index
 *   i_face_vtx_lst  <-> interior faces -> vertices connectivity
 *   n_i_faces       <-- number of internal faces
 *   clean_list      <-- sorted list of faces to remove (1 to n)
 *   clean_list_size <-- size of clean_list
 *
 * returns:
 *   number of element removed from i_face_lst
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_clean_i_faces(cs_lnum_t         *i_face_vtx_idx,
               cs_lnum_t         *i_face_vtx_lst,
               cs_lnum_t          n_i_faces,
               const cs_lnum_t   *clean_list,
               cs_lnum_t          clean_list_size)
{
  cs_lnum_t face_id, i;
  cs_lnum_t ind_empty = 0;
  cs_lnum_t ind_full = 0;

  cs_lnum_t l_shift = 0;

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    bool remove_face = false;

    cs_lnum_t start_id = i_face_vtx_idx[face_id];
    cs_lnum_t end_id = i_face_vtx_idx[face_id + 1];
    cs_lnum_t d = end_id - start_id;

    if (ind_empty < clean_list_size) {
      if (face_id == clean_list[ind_empty]) {
        remove_face = true;
        ind_empty++;
        l_shift += d;
      }
    }

    if (remove_face == false) {
      if (face_id != ind_full) {
        for (i = i_face_vtx_idx[ind_full];
             i < i_face_vtx_idx[ind_full] + d;
             i++)
          i_face_vtx_lst[i] = i_face_vtx_lst[i + l_shift];
      }
      i_face_vtx_idx[ind_full + 1] = i_face_vtx_idx[ind_full] + d;
      ind_full++;
    }

  }

  return l_shift;
}

/*----------------------------------------------------------------------------
 * Remove internal faces of i_face_cells from a index
 * (realloc has to be done after).
 *
 * parameters:
 *   i_face_idx      <-> interior faces -> cells connectivity
 *   n_i_faces       <-- number of internal faces
 *   clean_list      <-- sorted index of faces to remove (1 to n)
 *   clean_list_size <-- size of clean_list
 *----------------------------------------------------------------------------*/

static void
_clean_i_face_cells(cs_lnum_2_t      *i_face_cells,
                    cs_lnum_t         n_i_faces,
                    const cs_lnum_t  *clean_list,
                    cs_lnum_t         clean_list_size)
{
  cs_lnum_t face_id;
  cs_lnum_t ind_empty = 0;
  cs_lnum_t ind_full = 0;

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    bool remove_face = false;

    if (ind_empty < clean_list_size) {
      if (face_id == clean_list[ind_empty]) {
        remove_face = true;
        ind_empty++;
      }
    }

    if (remove_face == false) {
      if (face_id != ind_full) {
        i_face_cells[ind_full][0] = i_face_cells[face_id][0];
        i_face_cells[ind_full][1] = i_face_cells[face_id][1];
      }
      ind_full++;
    }

  }
}

/*----------------------------------------------------------------------------
 * Remove internal faces of i_face_family from a index
 * (realloc has to be done after).
 *
 * parameters:
 *   i_face_family   <-> interior faces family
 *   n_i_faces       <-- number of internal faces
 *   clean_list      <-- sorted index of faces to remove
 *   clean_list_size <-- size of clean_list
 *----------------------------------------------------------------------------*/

static void
_clean_i_family(cs_lnum_t        *i_face_family,
                cs_lnum_t         n_i_faces,
                const cs_lnum_t  *clean_list,
                cs_lnum_t         clean_list_size)
{
  cs_lnum_t face_id;
  cs_lnum_t ind_empty = 0;
  cs_lnum_t ind_full = 0;

  for (face_id = 0; face_id < n_i_faces; face_id++) {

    bool remove_face = false;

    if (ind_empty < clean_list_size) {
      if (face_id == clean_list[ind_empty]) {
        remove_face = true;
        ind_empty++;
      }
    }

    if (remove_face == false) {
      if (face_id != ind_full)
        i_face_family[ind_full] = i_face_family[face_id];
      ind_full++;
    }

  }
}

/*----------------------------------------------------------------------------
 * Get the complement of one list.
 *
 * parameters:
 *   list_c      --> the complement of list
 *                   (preallocated array of size list_c_size)
 *   list_c_size <-- size of list_c
 *   list        <-- sorted list of faces to remove  (0 to n-1)
 *   list_size   <-- size of list
 *----------------------------------------------------------------------------*/

static void
_get_list_c(cs_lnum_t        *list_c,
            cs_lnum_t         list_c_size,
            const cs_lnum_t  *list,
            cs_lnum_t         list_size)
{
  cs_lnum_t face_id;
  cs_lnum_t i_empty = 0;
  cs_lnum_t i_full = 0;

  for (face_id = 0; face_id < list_c_size; face_id++) {

    if (i_empty < list_size) {

      if (face_id != list[i_empty])
        list_c[i_full++] = face_id;
      else
        i_empty++;

    }
    else
      list_c[i_full++] = face_id;

  }
}

/*----------------------------------------------------------------------------
 * Count number of boundary faces to add in local rank
 * (realloc has to be done after).
 *
 * parameters:
 *   i_face_cells   <-- interior faces -> cells connectivity
 *   i_face_vtx_idx <-- interior faces -> vertices list
 *   list           <-- interior faces list to check (1 to n)
 *   list_size      <-- size of list
 *   count          --> preallocated array of size 2 :
 *                       count[0] = number of boundary faces to add
 *                       count[1] = number of i_face_vtx_idx corresponding
 *----------------------------------------------------------------------------*/

static void
_count_b_faces_to_add(const cs_lnum_2_t  *i_face_cells,
                      const cs_lnum_t    *i_face_vtx_idx,
                      const cs_lnum_t    *list,
                      cs_lnum_t           list_size,
                      cs_lnum_t          *count)
{
  cs_lnum_t ii;

  cs_lnum_t n_bf = 0;
  cs_lnum_t n_bf_lst = 0;

  for (ii = 0; ii < list_size; ii++) {
    if (i_face_cells[list[ii]][0] > -1) {
      n_bf++;
      n_bf_lst +=   i_face_vtx_idx[list[ii] + 1]
                  - i_face_vtx_idx[list[ii]];
    }
    if (i_face_cells[list[ii]][1] > -1) {
      n_bf++;
      n_bf_lst +=   i_face_vtx_idx[list[ii] + 1]
                  - i_face_vtx_idx[list[ii]];
    }
  }

  count[0] = n_bf;
  count[1] = n_bf_lst;
}

/*----------------------------------------------------------------------------
 * Add boundary faces to b_face_vtx_idx, b_face_vtx_lst, b_face_cells
 * and b_face_family from one list.
 *
 * Face connectivity and family arrays should be preallocated.
 *
 * parameters:
 *   mesh                    <-> pointer to mesh structure
 *   b_face_vtx_idx          <-> boundary faces -> vertices list
 *   b_face_vtx_lst          <-> boundary faces -> vertices connectivity
 *   b_face_cells            <-> boundary faces -> cells connectivity
 *   b_face_family           <-> boundary faces family list
 *   b_face_vtx_connect_size <-- new b_face_vtx_connect_size
 *   list                    <-- list of boundary faces to add
 *   list_size               <-- size of list
 *----------------------------------------------------------------------------*/

static void
_add_b_faces(cs_mesh_t        *mesh,
             cs_lnum_t        *b_face_vtx_idx,
             cs_lnum_t        *b_face_vtx_lst,
             cs_lnum_t        *b_face_cells,
             cs_lnum_t        *b_face_family,
             cs_lnum_t         b_face_vtx_connect_size,
             const cs_lnum_t  *list,
             cs_lnum_t         list_size)
 {
  cs_lnum_t ii, jj;
  cs_lnum_t n_face_vertices;
  cs_lnum_t inc = 0;
  cs_lnum_t *i_face_vtx_idx = mesh->i_face_vtx_idx;
  cs_lnum_t *i_face_vtx_lst = mesh->i_face_vtx_lst;
  cs_lnum_2_t *i_face_cells = mesh->i_face_cells;
  int        *i_face_family = mesh->i_face_family;
  cs_lnum_t n_b_faces = mesh->n_b_faces;

  for (ii = 0; ii < list_size; ii++) {

    n_face_vertices = i_face_vtx_idx[list[ii] + 1] - i_face_vtx_idx[list[ii]];

    if (i_face_cells[list[ii]][0] > -1) {

      b_face_cells[n_b_faces + inc]
        = i_face_cells[list[ii]][0];
      b_face_vtx_idx[n_b_faces + inc + 1]
        = b_face_vtx_idx[n_b_faces + inc] + n_face_vertices;

      b_face_family[n_b_faces + inc] = i_face_family[list[ii]];

      for (jj = 0; jj < n_face_vertices; jj++) {
        b_face_vtx_lst[b_face_vtx_connect_size + jj]
          = i_face_vtx_lst[i_face_vtx_idx[list[ii]] + jj];
      }

      inc++;
      b_face_vtx_connect_size += n_face_vertices;

    }

    if (i_face_cells[list[ii]][1] > -1) {

      b_face_cells[n_b_faces + inc]
        = i_face_cells[list[ii]][1];
      b_face_vtx_idx[n_b_faces + inc + 1 ]
        = b_face_vtx_idx[n_b_faces + inc] + n_face_vertices;

      b_face_family[n_b_faces + inc] = i_face_family[list[ii]];
      for (jj = 0; jj < n_face_vertices; jj++) {
        b_face_vtx_lst[b_face_vtx_connect_size + jj]
          = i_face_vtx_lst[   i_face_vtx_idx[list[ii]]
                           +  n_face_vertices - jj - 1];
      }

      inc++;
      b_face_vtx_connect_size += n_face_vertices;
    }

  }
}

/*----------------------------------------------------------------------------
 * Refresh global numering for new added boundary faces
 *
 * parameters:
 *   i_face_cells  <-- interior faces -> cells connectivity
 *   n_g_faces     <-- old number of global boundary faces
 *   n_b_faces     <-- old number of local boundary faces
 *   list          <-- list of boundary faces to add (with gost)
 *   list_size     <-- size of list
 *   list_glob_num <-- global numbering of elements in list
 *   glob_num      --> new global numbering of boundary faces
 *----------------------------------------------------------------------------*/

static void
_refresh_b_glob_num(const cs_lnum_2_t  *i_face_cells,
                    cs_gnum_t           n_g_faces,
                    cs_lnum_t           n_b_faces,
                    const cs_lnum_t    *list,
                    cs_lnum_t           list_size,
                    const cs_gnum_t    *list_glob_num,
                    cs_gnum_t          *b_faces_glob_num)
{
  cs_lnum_t ii;
  cs_lnum_t inc = 0;

  for (ii = 0; ii < list_size; ii++) {
    if (i_face_cells[list[ii]][0] > -1) {
      b_faces_glob_num[n_b_faces + inc]
        = n_g_faces + 2*(list_glob_num[ii] - 1) + 1;
      inc++;
    }
    if (i_face_cells[list[ii]][1] > -1) {
      b_faces_glob_num[n_b_faces + inc]
        = n_g_faces + 2*(list_glob_num[ii] - 1) + 2;
      inc++;
    }
  }
}

/*----------------------------------------------------------------------------
 * Print information on a mesh structure.
 *
 * parameters:
 *   mesh  <--  pointer to mesh structure.
 *   name  <--  associated name.
 *----------------------------------------------------------------------------*/

static void
_print_mesh_info(const cs_mesh_t  *mesh,
                 const char       *name)
{
  bft_printf(_(" %s\n"
               "     Number of cells:          %llu\n"
               "     Number of interior faces: %llu\n"
               "     Number of boundary faces: %llu\n"
               "     Number of vertices:       %llu\n"),
             name,
             (unsigned long long)(mesh->n_g_cells),
             (unsigned long long)(mesh->n_g_i_faces),
             (unsigned long long)(mesh->n_g_b_faces),
             (unsigned long long)(mesh->n_g_vertices));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a mesh face interface set in case of parallelism.
 *
 * \param[in]  mesh      pointer to mesh structure to modify
 *
 *
 * \return  pointer to interface set for faces, or NULL
 */
/*----------------------------------------------------------------------------*/

static cs_interface_set_t *
_build_faces_interface_set_if_needed(const cs_mesh_t  *mesh)
{
  cs_interface_set_t  *face_ifs = NULL;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_mesh_builder_t  *mb = cs_glob_mesh_builder;

    int *periodicity_num;
    BFT_MALLOC(periodicity_num, mb->n_perio, int);

    for (int i = 0; i < mb->n_perio; i++)
      periodicity_num[i] = i+1;

    face_ifs = cs_interface_set_create
                 (mesh->n_i_faces,
                  NULL,
                  mesh->global_i_face_num,
                  mesh->periodicity,
                  mb->n_perio,
                  periodicity_num,
                  mb->n_per_face_couples,
                  (const cs_gnum_t *const *)mb->per_face_couples);

    BFT_FREE(periodicity_num);

  }

#endif /* HAVE_MPI */

  return face_ifs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update periodicity info in builder and destroy temporary
 *  faces interface set in parallel mode.
 *
 * \param[in]       mesh      pointer to mesh structure to modify
 * \param[in, out]  face_ifs  pointer to pointer to face interface set
 */
/*----------------------------------------------------------------------------*/

static void
_destroy_faces_interface_set(const cs_mesh_t      *mesh,
                             cs_interface_set_t  **face_ifs)
{
#if defined(HAVE_MPI)

  cs_mesh_builder_t  *mb = cs_glob_mesh_builder;

  if (cs_glob_n_ranks > 1 && *face_ifs != NULL) {
    cs_mesh_builder_extract_periodic_faces_g(mesh->n_init_perio,
                                             mb,
                                             mesh->periodicity,
                                             mesh->global_i_face_num,
                                             *face_ifs);
    cs_interface_set_destroy(face_ifs);
  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 * The created faces share vertices.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \param[in, out]  mesh      pointer to mesh structure to modify
 * \param[in, out]  face_ifs  pointer to interface set for faces
 * \param[in]       n_faces   number of selected (interior) faces
 * \param[in, out]  ids       list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

static void
_boundary_insert(cs_mesh_t           *mesh,
                 cs_interface_set_t  *face_ifs,
                 cs_lnum_t            n_faces,
                 cs_lnum_t            face_id[])
{
  cs_lnum_t ii;
  cs_gnum_t _n_g_b_faces, _n_g_i_faces;
  cs_lnum_t i_face_vtx_cleaned;
  cs_lnum_t count[2];

  cs_lnum_t n_i_faces = mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  cs_lnum_t new_n_i_faces = n_i_faces - n_faces;

  cs_lnum_t i_face_vtx_connect_size = mesh->i_face_vtx_connect_size;
  cs_lnum_t b_face_vtx_connect_size = mesh->b_face_vtx_connect_size;

  cs_gnum_t *b_face_gnum = NULL;

  cs_mesh_builder_t *mb = cs_glob_mesh_builder;

  cs_lnum_t *face_id_c = NULL;

  qsort(face_id, n_faces, sizeof(cs_lnum_t), &_compare_nums);

  assert(mb != NULL);

  /* In case of periodicity, we also need to renumber periodic couples in the
     builder; in parallel,we transform info into an interface, which will be
     used at the end of the mesh update; in serial, the update is immediate */

  if (mb->n_perio > 0) {

    cs_lnum_t *renum = _i_list_renum(mesh->n_i_faces,
                                     n_faces,
                                     face_id);

    if (face_ifs != NULL)
      cs_interface_set_renumber(face_ifs, renum);

    if (cs_glob_n_ranks == 1) {

      for (int i = 0; i < mb->n_perio; i++) {
        cs_lnum_t  n_c_ini = mb->n_per_face_couples[i];
        cs_lnum_t  n_c_new = 0;
        cs_gnum_t *f_c = mb->per_face_couples[i];
        for (cs_lnum_t j = 0; j < n_c_ini; j++) {
          cs_lnum_t f0 = renum[f_c[j*2] - 1];
          cs_lnum_t f1 = renum[f_c[j*2 + 1] - 1];
          if (f0 > -1 && f1 > -1) {
            f_c[n_c_new*2]   = f0+1;
            f_c[n_c_new*2+1] = f1+1;
            n_c_new++;
          }
        }
        if (n_c_new != n_c_ini) {
          BFT_REALLOC(mb->per_face_couples[i], n_c_new*2, cs_gnum_t);
          mb->n_per_face_couples[i] = n_c_new;
          mb->n_g_per_face_couples[i] = n_c_new;
        }
      }
    }

    BFT_FREE(renum);
  }

  if (mesh->global_i_face_num != NULL || cs_glob_n_ranks > 1) {

    BFT_MALLOC(face_id_c, (n_i_faces - n_faces), cs_lnum_t);
    _get_list_c(face_id_c,
                n_i_faces,
                face_id,
                n_faces);

    fvm_io_num_t *global_number_i_faces
      = fvm_io_num_create_from_select(face_id_c,
                                      mesh->global_i_face_num,
                                      new_n_i_faces,
                                      0);

    fvm_io_num_t *global_number_b_faces
      = fvm_io_num_create_from_select(face_id,
                                      mesh->global_i_face_num,
                                      n_faces,
                                      0);

    b_face_gnum
      = fvm_io_num_transfer_global_num(global_number_b_faces);

    BFT_FREE(mesh->global_i_face_num);
    mesh->global_i_face_num
      = fvm_io_num_transfer_global_num(global_number_i_faces);

    global_number_i_faces = fvm_io_num_destroy(global_number_i_faces);
    global_number_b_faces = fvm_io_num_destroy(global_number_b_faces);
  }

  _count_b_faces_to_add((const cs_lnum_2_t  *)mesh->i_face_cells,
                        mesh->i_face_vtx_idx,
                        face_id,
                        n_faces,
                        count);

  BFT_REALLOC(mesh->b_face_vtx_idx, n_b_faces + count[0]  + 1, cs_lnum_t);
  BFT_REALLOC(mesh->b_face_cells, n_b_faces + count[0], cs_lnum_t);
  BFT_REALLOC(mesh->b_face_vtx_lst,
              b_face_vtx_connect_size + count[1],
              cs_lnum_t);
  BFT_REALLOC(mesh->b_face_family, n_b_faces + count[0], cs_lnum_t);

  _add_b_faces(mesh,
               mesh->b_face_vtx_idx,
               mesh->b_face_vtx_lst,
               mesh->b_face_cells,
               mesh->b_face_family,
               b_face_vtx_connect_size,
               face_id,
               n_faces);

  mesh->n_b_faces = n_b_faces + count[0];
  mesh->b_face_vtx_connect_size =  b_face_vtx_connect_size + count[1];

  _n_g_b_faces = mesh->n_b_faces;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t _n_l_b_faces = mesh->n_b_faces;
    MPI_Allreduce(&_n_l_b_faces, &_n_g_b_faces, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    BFT_REALLOC(mesh->global_b_face_num, mesh->n_b_faces, cs_gnum_t);

    _refresh_b_glob_num((const cs_lnum_2_t *)mesh->i_face_cells,
                        mesh->n_g_b_faces,
                        n_b_faces,
                        face_id,
                        n_faces,
                        b_face_gnum,
                        mesh->global_b_face_num);

    BFT_FREE(b_face_gnum);

  }

#endif /* HAVE_MPI */

  i_face_vtx_cleaned = _clean_i_faces(mesh->i_face_vtx_idx,
                                      mesh->i_face_vtx_lst,
                                      n_i_faces,
                                      face_id,
                                      n_faces);

  _clean_i_face_cells(mesh->i_face_cells,
                      mesh->n_i_faces,
                      face_id,
                      n_faces);

  _clean_i_family(mesh->i_face_family,
                  n_i_faces,
                  face_id,
                  n_faces);

  mesh->n_i_faces = n_i_faces - n_faces;
  mesh->i_face_vtx_connect_size = i_face_vtx_connect_size - i_face_vtx_cleaned;

  _n_g_i_faces = mesh->n_i_faces;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t n_bf = 0;
    for (ii = 0; ii < mesh->n_i_faces; ii++) {
      if (mesh->i_face_cells[ii][1] > -1)
        n_bf++;
    }
    MPI_Allreduce(&n_bf, &_n_g_i_faces, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  BFT_REALLOC(mesh->i_face_vtx_idx, mesh->n_i_faces + 1, cs_lnum_t);
  BFT_REALLOC(mesh->i_face_vtx_lst, mesh->i_face_vtx_connect_size, cs_lnum_t);
  BFT_REALLOC(mesh->i_face_cells, mesh->n_i_faces, cs_lnum_2_t);
  BFT_REALLOC(mesh->i_face_family, mesh->n_i_faces, cs_lnum_t);

  mesh->n_g_b_faces = _n_g_b_faces;
  mesh->n_g_i_faces = _n_g_i_faces;

  if (cs_glob_n_ranks > 1)
    BFT_FREE(face_id_c);

  _print_mesh_info(mesh, " After addition of user-defined thin walls");
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 * The created faces share vertices.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \param[in]       mesh     pointer to mesh structure to modify
 * \param[in]       n_faces  number of selected (interior) faces
 * \param[in, out]  ids      list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert(cs_mesh_t  *mesh,
                        cs_lnum_t   n_faces,
                        cs_lnum_t   face_id[])
{
  cs_interface_set_t *face_ifs = _build_faces_interface_set_if_needed(mesh);

  _boundary_insert(mesh, NULL, n_faces, face_id);

  _destroy_faces_interface_set(mesh, &face_ifs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert a boundary into the mesh between a given set of cells
 *        and the the others.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * An optional group name may be assigned to newly created boundary faces.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * \param[in]  mesh        pointer to mesh structure to modify
 * \param[in]  group_name  group name to assign to newly created boundary
 *                         faces, or NULL
 * \param[in]  n_cells     number of separated cells
 * \param[in]  cell_id     separated cell ids
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert_separating_cells(cs_mesh_t        *mesh,
                                         const char       *group_name,
                                         cs_lnum_t         n_cells,
                                         const cs_lnum_t   cell_id[])
{
  cs_lnum_t  n_m_cells = mesh->n_cells;
  cs_lnum_t  n_i_faces = mesh->n_i_faces;

  /* Mark cells and matching faces */

  int32_t    *cell_tag;
  cs_lnum_t  *face_tag;

  BFT_MALLOC(face_tag, n_i_faces, cs_lnum_t);
  BFT_MALLOC(cell_tag, n_m_cells, int32_t);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    face_tag[f_id] = 0;

  for (cs_lnum_t c_id = 0; c_id < n_m_cells; c_id++)
    cell_tag[c_id] = -1;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    cell_tag[cell_id[c_id]] = 1;

  if (mesh->halo != NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("%s should be called before halo creation."),
              __func__);

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
    cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
    if (c_id_0 > -1 && c_id_0 < n_m_cells)
      face_tag[f_id] += cell_tag[c_id_0];
    if (c_id_1 > -1 && c_id_1 < n_m_cells)
      face_tag[f_id] += cell_tag[c_id_1];
  }

  /* Now synchronize tags in case of periodicity or parallelism */

  cs_interface_set_t  *face_ifs
    = _build_faces_interface_set_if_needed(mesh);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    cs_interface_set_sum(face_ifs,
                         n_i_faces,
                         1,
                         true,
                         CS_LNUM_TYPE,
                         face_tag);

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1) {

    cs_mesh_builder_t *mb = cs_glob_mesh_builder;

    for (int i = 0; i < mb->n_perio; i++) {
      cs_lnum_t  n_c_ini = mb->n_per_face_couples[i];
      cs_gnum_t *f_c = mb->per_face_couples[i];
      for (cs_lnum_t j = 0; j < n_c_ini; j++) {
        cs_lnum_t f0 = f_c[j*2] - 1;
        cs_lnum_t f1 = f_c[j*2 + 1] - 1;
        int f_tag = face_tag[f0] + face_tag[f1];
        face_tag[f0] = f_tag;
        face_tag[f1] = f_tag;
      }
    }
  }

  BFT_FREE(cell_tag);

  /* Transform face tag into list */

  cs_lnum_t b_count = 0;

  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {
    if (face_tag[f_id] == 0)
      face_tag[b_count++] = f_id;
  }

  BFT_REALLOC(face_tag, b_count, cs_lnum_t);

  _boundary_insert(mesh, face_ifs, b_count, face_tag);

  _destroy_faces_interface_set(mesh, &face_ifs);

  BFT_FREE(face_tag);
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
