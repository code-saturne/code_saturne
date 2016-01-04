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

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_mesh_thinwall.h"

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
  int ii, jj;
  cs_lnum_t   n_face_vertices;
  cs_lnum_t   inc = 0;
  cs_lnum_t  *i_face_vtx_idx = mesh->i_face_vtx_idx;
  cs_lnum_t  *i_face_vtx_lst = mesh->i_face_vtx_lst;
  cs_lnum_2_t  *i_face_cells = mesh->i_face_cells;
  int        *i_face_family = mesh->i_face_family;
  cs_lnum_t   n_b_faces = mesh->n_b_faces;

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Insert thin walls into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the list of selected faces is sorted by this function.
 *
 * parameters:
 *   mesh           <->  pointer to mesh structure to modify
 *   face_list      <-> list of selected (interior) faces (0 to n-1)
 *   face_list_size <-> number of selected (interior) faces
 *----------------------------------------------------------------------------*/

void
cs_create_thinwall(cs_mesh_t  *mesh,
                   cs_lnum_t  *face_list,
                   cs_lnum_t   face_list_size)
 {
  cs_lnum_t ii;
  cs_gnum_t _n_g_b_faces, _n_g_i_faces;
  cs_lnum_t i_face_vtx_cleaned;
  cs_lnum_t count[2];

  cs_lnum_t n_i_faces =  mesh->n_i_faces;
  cs_lnum_t n_b_faces = mesh->n_b_faces;
  cs_lnum_t new_n_i_faces = n_i_faces - face_list_size;

  cs_lnum_t i_face_vtx_connect_size = mesh->i_face_vtx_connect_size;
  cs_lnum_t b_face_vtx_connect_size = mesh->b_face_vtx_connect_size;

  fvm_io_num_t *global_number_i_faces = NULL;
  fvm_io_num_t *global_number_b_faces = NULL;
  const cs_gnum_t *global_order_i_faces = NULL;
  const cs_gnum_t *global_order_b_faces = NULL;

  cs_lnum_t *face_list_c = NULL;

  qsort(face_list, face_list_size, sizeof(cs_lnum_t), &_compare_nums);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    BFT_MALLOC(face_list_c, (n_i_faces - face_list_size), cs_lnum_t);
    _get_list_c(face_list_c,
                n_i_faces,
                face_list,
                face_list_size);

    global_number_i_faces
      = fvm_io_num_create_from_select(face_list_c,
                                      mesh->global_i_face_num,
                                      new_n_i_faces,
                                      0);

    global_order_i_faces = fvm_io_num_get_global_num(global_number_i_faces);

    global_number_b_faces
      = fvm_io_num_create_from_select(face_list,
                                      mesh->global_i_face_num,
                                      face_list_size,
                                      0);
    global_order_b_faces = fvm_io_num_get_global_num(global_number_b_faces);

    BFT_REALLOC(mesh->global_i_face_num, new_n_i_faces, cs_gnum_t);

    memcpy(mesh->global_i_face_num,
           global_order_i_faces,
           new_n_i_faces*sizeof(cs_gnum_t));

  }

#endif /* HAVE_MPI */

  _count_b_faces_to_add(mesh->i_face_cells,
                        mesh->i_face_vtx_idx,
                        face_list,
                        face_list_size,
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
               face_list,
               face_list_size);

  mesh->n_b_faces = n_b_faces + count[0];
  mesh->b_face_vtx_connect_size =  b_face_vtx_connect_size + count[1];

  _n_g_b_faces = mesh->n_b_faces;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_gnum_t _n_l_b_faces = mesh->n_b_faces;
    MPI_Allreduce(&_n_l_b_faces, &_n_g_b_faces, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);

    BFT_REALLOC(mesh->global_b_face_num, mesh->n_b_faces, cs_gnum_t);

    _refresh_b_glob_num(mesh->i_face_cells,
                        mesh->n_g_b_faces,
                        n_b_faces,
                        face_list,
                        face_list_size,
                        global_order_b_faces,
                        mesh->global_b_face_num);

  }

#endif /* HAVE_MPI */

  i_face_vtx_cleaned = _clean_i_faces(mesh->i_face_vtx_idx,
                                      mesh->i_face_vtx_lst,
                                      n_i_faces,
                                      face_list,
                                      face_list_size);

  _clean_i_face_cells(mesh->i_face_cells,
                      mesh->n_i_faces,
                      face_list,
                      face_list_size);

  _clean_i_family(mesh->i_face_family,
                  n_i_faces,
                  face_list,
                  face_list_size);

  mesh->n_i_faces = n_i_faces - face_list_size;
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

  _print_mesh_info(mesh, " After addition of user-defined thin walls");

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    BFT_FREE(face_list_c);
#endif
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
