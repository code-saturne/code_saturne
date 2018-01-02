/*============================================================================
 * Insert boundaries into the mesh.
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
#include "cs_mesh_connect.h"
#include "cs_mesh_group.h"
#include "cs_parall.h"
#include "cs_sort.h"

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
 * Remove internal faces of i_face_cells from an index
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
 *   n_g_b_faces   <-- old number of global boundary faces
 *   n_b_faces     <-- old number of local boundary faces
 *   list          <-- list of boundary faces to add (with gost)
 *   list_size     <-- size of list
 *   list_glob_num <-- global numbering of elements in list
 *   glob_num      --> new global numbering of boundary faces
 *----------------------------------------------------------------------------*/

static void
_refresh_b_glob_num(const cs_lnum_2_t  *i_face_cells,
                    cs_gnum_t           n_g_b_faces,
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
        = n_g_b_faces + 2*(list_glob_num[ii] - 1) + 1;
      inc++;
    }
    if (i_face_cells[list[ii]][1] > -1) {
      b_faces_glob_num[n_b_faces + inc]
        = n_g_b_faces + 2*(list_glob_num[ii] - 1) + 2;
      inc++;
    }
  }
}

/*----------------------------------------------------------------------------
 * Substep of vertex separaion at inserted boundary faces based on their
 * cell adjacency:
 * - create tuples of global vertex and lowest adjacent cell numbers
 *   for a partial indexed vertex to cell adjacency, keeping only tuples
 *   where the lowest adjacent cell number for a given vertex/cell couple
 *   is higher than the lowest adjacent cell number of that vertex.
 *
 * It is the caller's responsibility to free the returned tuples array.
 *
 * parameters:
 *   mesh          <-- pointer to mesh structure
 *   v2c_idx       <-- index of partial vertex to cell adjacency
 *   v_g_min_c_num <-- minimum global cell number associated with vertex
 *                     to cell adjacency
 *   c_num_min     <-- minimum cell global number adjacent to a given vertex
 *   n_tuples      --> number of returned tuples
 *   tuple_vtx_id  --> vertex id associated with each tuple
 *   tuples        --> tuple values
 *----------------------------------------------------------------------------*/

static void
_split_vertices_vertex_cell_tuples(cs_mesh_t         *mesh,
                                   const cs_lnum_t    v2c_idx[],
                                   const cs_gnum_t    v_g_min_c_num[],
                                   const cs_gnum_t    c_num_min[],
                                   cs_lnum_t         *n_tuples,
                                   cs_lnum_t        **tuple_vtx_id,
                                   cs_gnum_t        **tuples)
{
  cs_lnum_t v_tmp_size = 0;
  const cs_lnum_t n_vertices = mesh->n_vertices;

  cs_lnum_t *_tuple_vtx_id = NULL;
  cs_gnum_t *v_g_min_c_add = NULL, *v_tmp = NULL;

  BFT_MALLOC(_tuple_vtx_id, v2c_idx[n_vertices], cs_lnum_t);
  BFT_MALLOC(v_g_min_c_add, v2c_idx[n_vertices]*2, cs_gnum_t);

  cs_lnum_t k = 0;

  for (cs_lnum_t i = 0; i < n_vertices; i++) {
    if (v2c_idx[i+1] > v2c_idx[i]) {
      cs_lnum_t nc = v2c_idx[i+1] - v2c_idx[i];
      if (nc > v_tmp_size) {
        v_tmp_size = nc * 2;
        BFT_REALLOC(v_tmp, v_tmp_size, cs_gnum_t);
      }
      cs_lnum_t s_id = v2c_idx[i];
      cs_lnum_t e_id = v2c_idx[i+1];
      cs_lnum_t l = 0;
      for (cs_lnum_t j = s_id; j < e_id; j++)
        v_tmp[l++] = v_g_min_c_num[j];
      cs_sort_gnum_shell(0, nc, v_tmp);
      cs_gnum_t g_num_prv = c_num_min[i];
      for (l = 0; l < nc; l++) {
        if (v_tmp[l] != g_num_prv) {
          _tuple_vtx_id[k] = i;
          v_g_min_c_add[k*2+1] = v_tmp[l];
          g_num_prv = v_tmp[l];
          k++;
        }
      }
    }
  }

  BFT_FREE(v_tmp);

  /* Resize arrays */

  BFT_REALLOC(_tuple_vtx_id, k, cs_lnum_t);
  BFT_REALLOC(v_g_min_c_add, k*2, cs_gnum_t);

  /* Transform vertex index to global number */

  if (mesh->global_vtx_num != NULL) {
    for (cs_lnum_t j = 0; j < k; j++)
      v_g_min_c_add[j*2] = mesh->global_vtx_num[_tuple_vtx_id[j]];
  }
  else {
    for (cs_lnum_t j = 0; j < k; j++)
      v_g_min_c_add[j*2] = _tuple_vtx_id[j] + 1;
  }

  /* return tuples */

  *n_tuples = k;
  *tuple_vtx_id = _tuple_vtx_id;
  *tuples = v_g_min_c_add;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a mesh face interface set in case of parallelism for
 *        use by the split vertices algorithm.
 *
 * The construction is restricted to faces with a missing adjacency
 * and flagged vertices.
 *
 * The caller is responsible for freeing the ifs_face_id array and
 * ifs structure.
 *
 * \param[in]   mesh         pointer to mesh structure to modify
 * \param[in]   v_flag       flag for vertices lying on inserted boundary
 * \param[out]  n_ifs_faces  number of selected faces on parallel interface
 * \param[out]  ifs_face_id  ids of selected faces on parallel interface
 * \param[out]  ifs          matching face interface set
 *
 * \return  pointer to interface set for selected faces, or NULL
 */
/*----------------------------------------------------------------------------*/

static void
_split_vertices_build_faces_interface(const cs_mesh_t      *mesh,
                                      const char            v_flag[],
                                      cs_lnum_t            *n_ifs_faces,
                                      cs_lnum_t           **ifs_face_id,
                                      cs_interface_set_t  **ifs)
{
  cs_lnum_t  _n_ifs_faces = 0;
  cs_lnum_t  *_ifs_face_id = NULL;

  cs_interface_set_t  *_face_ifs = NULL;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    BFT_MALLOC(_ifs_face_id, mesh->n_i_faces, cs_lnum_t);

    for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {

      cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
      cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];

      if (c_id_0 < 0 || c_id_1 < 0) {
        cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
        cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
        for (cs_lnum_t i = s_id; i < e_id; i++) {
          cs_lnum_t vtx_id = mesh->i_face_vtx_lst[i];
          if (v_flag[vtx_id] != 0) {
            _ifs_face_id[_n_ifs_faces] = f_id;
            _n_ifs_faces += 1;
            break;
          }
        }
      }

    }

    BFT_REALLOC(_ifs_face_id, _n_ifs_faces, cs_lnum_t);

    _face_ifs = cs_interface_set_create(_n_ifs_faces,
                                        _ifs_face_id,
                                        mesh->global_i_face_num,
                                        NULL,
                                        0,
                                        NULL,
                                        NULL,
                                        NULL);

  }

#endif /* HAVE_MPI */

  *n_ifs_faces = _n_ifs_faces;
  *ifs_face_id = _ifs_face_id;
  *ifs = _face_ifs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build indexes for parallel update for use by the
 *        \ref _split_vertices_vtx_min_cell_num function.
 *
 * The caller is responsible for freeing the returned send_idx and
 * recv_idx arrays.
 *
 * \param[in]   mesh         pointer to mesh structure to modify
 * \param[in]   ifs          matching face interface set
 * \param[in]   ifs_face_id  ids of selected faces on interface
 * \param[in]   v2c_idx      vertex to cell index
 * \param[out]  send_idx     interface set send index
 * \param[out]  recv_idx     interface set receive index
 */
/*----------------------------------------------------------------------------*/

static void
_split_vertices_vtx_min_cell_num_ifs_index(cs_mesh_t            *mesh,
                                           cs_interface_set_t   *ifs,
                                           const cs_lnum_t       ifs_face_id[],
                                           const cs_lnum_t       v2c_idx[],
                                           cs_lnum_t           **send_idx,
                                           cs_lnum_t           **recv_idx)
{
  /* Exchange v2c info */

  int n_interfaces = cs_interface_set_size(ifs);

  cs_lnum_t ifs_tot_size = cs_interface_set_n_elts(ifs);

  cs_lnum_t *_send_idx = NULL, *_recv_idx = NULL;

  BFT_MALLOC(_send_idx, ifs_tot_size + 1, cs_lnum_t);
  BFT_MALLOC(_recv_idx, ifs_tot_size + 1, cs_lnum_t);

  /* Counting pass for v2v_info exchange
     (send/receive indexes are prepared as counts) */

  cs_lnum_t if_shift = 0;

  _send_idx[0] = 0;
  _recv_idx[0] = 0;

  for (int i = 0; i < n_interfaces; i++) {

    const cs_interface_t *interface = cs_interface_set_get(ifs, i);
    const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);
    const cs_lnum_t if_size = cs_interface_size(interface);

    for (cs_lnum_t j = 0; j < if_size; j++) {

      cs_lnum_t f_id = ifs_face_id[local_id[j]];
      cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
      cs_lnum_t lc = 0;

      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t vtx_id = mesh->i_face_vtx_lst[k];
        if (v2c_idx[vtx_id+1] > v2c_idx[vtx_id]) lc += 1;
      }

      _send_idx[if_shift + j + 1] = lc;

    }

    if_shift += if_size;

  }

  cs_interface_set_copy_array(ifs,
                              CS_LNUM_TYPE,
                              1,
                              false,
                              _send_idx + 1,
                              _recv_idx + 1);

  /* Transform count to index */

  for (cs_lnum_t j = 0; j < ifs_tot_size; j++) {
    _send_idx[j+1] += _send_idx[j];
    _recv_idx[j+1] += _recv_idx[j];
  }

  *send_idx = _send_idx;
  *recv_idx = _recv_idx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Parallel update for determination of  minimum global cell number
 *        adjacent to each possibly split vertex for a given initial
 *        cell number.
 *
 * \param[in]       mesh           pointer to mesh structure to modify
 * \param[in]       ifs            matching face interface set
 * \param[in]       ifs_face_id    ids of selected faces on interface
 * \param[in]       send_idx       send index for interface set
 * \param[in]       recv_idx       receive index for interface set
 * \param[in]       v2c_idx        vertex to cell index
 * \param[in]       v2c            vertex to cell adjacency
 * \param[out]      v2c_send       send buffer
 * \param[out]      v2c_recv       receive buffer
 * \param[in, out]  v_g_min_c_num  vertex to minimum adjacent global cell
 *                                 adjacency (based on the v2c_idx index)
 *
 * \return number of local changes
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t
_split_vertices_vtx_min_cell_num_exchange(cs_mesh_t           *mesh,
                                          cs_interface_set_t  *ifs,
                                          const cs_lnum_t      ifs_face_id[],
                                          const cs_lnum_t      send_idx[],
                                          const cs_lnum_t      recv_idx[],
                                          const cs_lnum_t      v2c_idx[],
                                          const cs_lnum_t      v2c[],
                                          cs_gnum_t            v2c_send[],
                                          cs_gnum_t            v2c_recv[],
                                          cs_gnum_t            v_g_min_c_num[])
{
  /* Exchange v2c info */

  int n_interfaces = cs_interface_set_size(ifs);

  cs_lnum_t if_shift = 0;

  for (int i = 0; i < n_interfaces; i++) {

    const cs_interface_t *interface = cs_interface_set_get(ifs, i);
    const cs_lnum_t if_size = cs_interface_size(interface);

    const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);

    for (cs_lnum_t j = 0; j < if_size; j++) {

      cs_lnum_t f_id = ifs_face_id[local_id[j]];
      cs_lnum_t c_id = mesh->i_face_cells[f_id][0];
      if (c_id < 0)
        c_id = mesh->i_face_cells[f_id][1];

      assert(   (   mesh->i_face_cells[f_id][0] < 0
                 || mesh->i_face_cells[f_id][1] < 0)
             && c_id > -1);

      cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
      cs_lnum_t lc = 0;

      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t vtx_id = mesh->i_face_vtx_lst[k];
        cs_lnum_t v_s_id = v2c_idx[vtx_id];
        cs_lnum_t v_e_id = v2c_idx[vtx_id + 1];
        for (cs_lnum_t l = v_s_id; l < v_e_id; l++) {
          if (v2c[l] == c_id) {
            v2c_send[send_idx[if_shift + j] + lc] = v_g_min_c_num[l];
            lc++;
            break;
          }
        }
      }

    }

    if_shift += if_size;

  }

  cs_interface_set_copy_indexed(ifs,
                                CS_GNUM_TYPE,
                                false,
                                send_idx,
                                recv_idx,
                                v2c_send,
                                v2c_recv);

  cs_gnum_t n_changes = 0;

  if_shift = 0;

  for (int i = 0; i < n_interfaces; i++) {

    const cs_interface_t *interface = cs_interface_set_get(ifs, i);
    const cs_lnum_t if_size = cs_interface_size(interface);

    const cs_lnum_t *local_id = cs_interface_get_elt_ids(interface);

    for (cs_lnum_t j = 0; j < if_size; j++) {

      cs_lnum_t f_id = ifs_face_id[local_id[j]];
      cs_lnum_t c_id = mesh->i_face_cells[f_id][0];
      if (c_id < 0)
        c_id = mesh->i_face_cells[f_id][1];

      cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
      cs_lnum_t lc = 0;

      for (cs_lnum_t k = s_id; k < e_id; k++) {
        cs_lnum_t vtx_id = mesh->i_face_vtx_lst[k];
        cs_lnum_t v_s_id = v2c_idx[vtx_id];
        cs_lnum_t v_e_id = v2c_idx[vtx_id + 1];
        for (cs_lnum_t l = v_s_id; l < v_e_id; l++) {
          if (v2c[l] == c_id) {
            cs_gnum_t g_num_cmp = v2c_recv[recv_idx[if_shift + j] + lc];
            if (v_g_min_c_num[l] > g_num_cmp) {
              v_g_min_c_num[l] = g_num_cmp;
              n_changes += 1;
            }
            lc++;
            break;
          }
        }
      }

    }

    if_shift += if_size;

  }

  return n_changes;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine minimum global cell number adjacent to each possibly
 *        split vertex for a given initial cell number.
 *
 * This is done iteratively by propagating interior face -> cell adjacencies.
 *
 * The caller is responsible for freeing the returned array
 *
 * \param[in]   mesh     pointer to mesh structure to modify
 * \param[in]   v_flag   flag for vertices lying on inserted boundary
 * \param[in]   v2c_idx  vertex to cell index
 * \param[in]   v2c      vertex to cell adjacency (non-empty only for
 *                       flagged vertices)
 *
 * \return  vertex to minimum adjacent global cell adjacency (based on the
 *          v2c_idx index)
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t *
_split_vertices_vtx_min_cell_num(cs_mesh_t       *mesh,
                                 const char       v_flag[],
                                 const cs_lnum_t  v2c_idx[],
                                 const cs_lnum_t  v2c[])
{
  cs_lnum_t n_vertices = mesh->n_vertices;

  cs_gnum_t *v_g_min_c_num;

  BFT_MALLOC(v_g_min_c_num, v2c_idx[n_vertices], cs_gnum_t);

  if (mesh->global_cell_num != NULL) {
    for (cs_lnum_t j = 0; j < v2c_idx[n_vertices]; j++)
      v_g_min_c_num[j] = mesh->global_cell_num[v2c[j]];
  }
  else {
    for (cs_lnum_t j = 0; j < v2c_idx[n_vertices]; j++)
      v_g_min_c_num[j] = v2c[j] + 1;
  }

  /* In parallel cases, we will need to communicate across matching
     interior faces containing flagged vertices */

#if defined(HAVE_MPI)

  cs_lnum_t  n_ifs_faces = 0;
  cs_lnum_t  *ifs_face_id = NULL;
  cs_interface_set_t  *face_ifs = NULL;
  cs_lnum_t *send_idx = NULL, *recv_idx = NULL;
  cs_gnum_t *v2c_send = NULL, *v2c_recv = NULL;

  _split_vertices_build_faces_interface(mesh,
                                        v_flag,
                                        &n_ifs_faces,
                                        &ifs_face_id,
                                        &face_ifs);

  if (face_ifs != NULL) {
    cs_lnum_t face_ifs_tot_size = cs_interface_set_n_elts(face_ifs);
    _split_vertices_vtx_min_cell_num_ifs_index(mesh,
                                               face_ifs,
                                               ifs_face_id,
                                               v2c_idx,
                                               &send_idx,
                                               &recv_idx);
    BFT_MALLOC(v2c_send, send_idx[face_ifs_tot_size], cs_gnum_t);
    BFT_MALLOC(v2c_recv, recv_idx[face_ifs_tot_size], cs_gnum_t);
  }

#endif /* HAVE_MPI */

  /* Now run multiple passes, marking equivalences as we go */

  cs_gnum_t n_changes = 0;

  do {

    n_changes = 0;

    /* loop on local faces */

    for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
      cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
      cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
      for (cs_lnum_t i = s_id; i < e_id; i++) {
        cs_lnum_t vtx_id = mesh->i_face_vtx_lst[i];
        if (v_flag[vtx_id] != 0) {
          cs_lnum_t c_id_0 = mesh->i_face_cells[f_id][0];
          cs_lnum_t c_id_1 = mesh->i_face_cells[f_id][1];
          cs_lnum_t v_s_id = v2c_idx[vtx_id];
          cs_lnum_t v_e_id = v2c_idx[vtx_id + 1];
          cs_lnum_t j0 = -1, j1 = -1;
          for (cs_lnum_t j = v_s_id; j < v_e_id; j++) {
            if (v2c[j] == c_id_0)
              j0 = j;
            else if (v2c[j] == c_id_1)
              j1 = j;
          }
          if (j0 > -1 && j1 > -1) {
            if (v_g_min_c_num[j0] < v_g_min_c_num[j1]) {
              v_g_min_c_num[j1] = v_g_min_c_num[j0];
              n_changes += 1;
            }
            else if (v_g_min_c_num[j0] > v_g_min_c_num[j1]) {
              v_g_min_c_num[j0] = v_g_min_c_num[j1];
              n_changes += 1;
            }
          }
        }
      }
    }

    /* parallel update */

#if defined(HAVE_MPI)
    if (face_ifs != NULL)
      n_changes += _split_vertices_vtx_min_cell_num_exchange(mesh,
                                                             face_ifs,
                                                             ifs_face_id,
                                                             send_idx,
                                                             recv_idx,
                                                             v2c_idx,
                                                             v2c,
                                                             v2c_send,
                                                             v2c_recv,
                                                             v_g_min_c_num);
#endif /* HAVE_MPI */

    cs_parall_counter(&n_changes, 1);

  } while (n_changes > 0);

  /* Free memory */

#if defined(HAVE_MPI)

  BFT_FREE(v2c_recv);
  BFT_FREE(v2c_send);

  BFT_FREE(recv_idx);
  BFT_FREE(send_idx);

  if (face_ifs != NULL)
    cs_interface_set_destroy(&face_ifs);

  BFT_FREE(ifs_face_id);

#endif

  return v_g_min_c_num;
}

/*----------------------------------------------------------------------------
 * Seperate vertices at inserted boundary faces based on their cell adjacency.
 *
 * parameters:
 *   mesh         <-> pointer to mesh structure
 *   b_face_start <-- start id for added boundary faces
 *   b_face_end   <-- past-the-end id for added boundary faces
 *----------------------------------------------------------------------------*/

static void
_split_vertices(cs_mesh_t  *mesh,
                cs_lnum_t   b_face_start,
                cs_lnum_t   b_face_end)
{
  cs_lnum_t n_vertices = mesh->n_vertices;

  /* We should not have a mesh vertices interface if this function
     is called before building ghost cells, but plan for the case where
     it is called after that (except where we also have periodicity */

  bool rebuild_vtx_interfaces = false;
  if (mesh->vtx_interfaces != NULL) {
    rebuild_vtx_interfaces = true;
    cs_interface_set_destroy(&(mesh->vtx_interfaces));
  }

  /* Use a local vertices interface set, ignoring periodicity */

#if defined(HAVE_MPI)
  cs_interface_set_t *ifs = NULL;

  if (cs_glob_n_ranks > 1)
    ifs = cs_interface_set_create(n_vertices, NULL, mesh->global_vtx_num,
                                  NULL, 0, NULL, NULL, NULL);
#endif /* defined(HAVE_MPI) */

  /* Mark vertices which may be split (vertices lying on new boundary faces) */

  char *v_flag;
  BFT_MALLOC(v_flag, n_vertices, char);

  for (cs_lnum_t i = 0; i < n_vertices; i++)
    v_flag[i] = 0;

  for (cs_lnum_t f_id = b_face_start; f_id < b_face_end; f_id++) {
    cs_lnum_t s_id = mesh->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->b_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++)
      v_flag[mesh->b_face_vtx_lst[i]] = 1;
  }

#if defined(HAVE_MPI)

  if (ifs != NULL)
    cs_interface_set_max(ifs, n_vertices, 1, true, CS_CHAR, v_flag);

#endif /* defined(HAVE_MPI) */

  /* Now build vertex -> cells index */

  cs_lnum_t  *v2c_idx = NULL, *v2c = NULL;

  cs_mesh_connect_vertices_to_cells(mesh,
                                    v_flag,
                                    &v2c_idx,
                                    &v2c);

  cs_gnum_t *v_g_min_c_num
    = _split_vertices_vtx_min_cell_num(mesh, v_flag, v2c_idx, v2c);

  /* Vertices with minimum adjacent cell global id will not be modified,
     so we determine the lowest minimum ajacent global id */

  cs_gnum_t *c_num_min;
  BFT_MALLOC(c_num_min, n_vertices, cs_gnum_t);

  for (cs_lnum_t i = 0; i < n_vertices; i++) {
    if (v2c_idx[i+1] == v2c_idx[i])
      c_num_min[i] = 0;
    else {
      cs_lnum_t s_id = v2c_idx[i];
      cs_lnum_t e_id = v2c_idx[i+1];
      c_num_min[i] = v_g_min_c_num[s_id];
      for (cs_lnum_t j = s_id+1; j < e_id; j++) {
        if (c_num_min[i] > v_g_min_c_num[j])
          c_num_min[i] = v_g_min_c_num[j];
      }
    }
  }

  /* Last synchronization using vertices interface */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_interface_set_min(ifs, n_vertices, 1, true, CS_GNUM_TYPE, c_num_min);
    cs_interface_set_destroy(&ifs);
  }
#endif

  /* We now build a set of referenced global vertex and cell number tuples,
     filtering out the minimum */

  cs_lnum_t  n_add_vertices = 0;
  cs_lnum_t  *v_min_c_add = NULL;
  cs_gnum_t  *v_g_min_c_add = NULL;

  _split_vertices_vertex_cell_tuples(mesh,
                                     v2c_idx,
                                     v_g_min_c_num,
                                     c_num_min,
                                     &n_add_vertices,
                                     &v_min_c_add,
                                     &v_g_min_c_add);

  /* Update local size */

  mesh->n_vertices += n_add_vertices;

  /* Update coordinates for added vertices */

  if (n_add_vertices > 0) {

    BFT_REALLOC(mesh->vtx_coord, mesh->n_vertices*3, cs_real_t);

    for (cs_lnum_t i = 0; i < n_add_vertices; i++) {
      cs_lnum_t j = v_min_c_add[i];
      cs_lnum_t k = n_vertices + i;
      for (cs_lnum_t l = 0; l < 3; l++)
        mesh->vtx_coord[k*3 + l] = mesh->vtx_coord[j*3 + l];
    }

  }

  BFT_FREE(v_min_c_add);

  /* Also update global numbering if required */

  if (cs_glob_n_ranks > 1 || mesh->global_vtx_num != NULL) {

    fvm_io_num_t *v_add_io_num = fvm_io_num_create_from_adj_s(NULL,
                                                              v_g_min_c_add,
                                                              n_add_vertices,
                                                              2);

    assert(fvm_io_num_get_local_count(v_add_io_num) == n_add_vertices);

    cs_gnum_t n_g_add_vertices = fvm_io_num_get_global_count(v_add_io_num);
    const cs_gnum_t *v_add_gnum = fvm_io_num_get_global_num(v_add_io_num);

    if (n_add_vertices > 0) {
      BFT_REALLOC(mesh->global_vtx_num, n_vertices + n_add_vertices, cs_gnum_t);
      for (cs_lnum_t i = 0; i < n_add_vertices; i++)
        mesh->global_vtx_num[n_vertices + i] = mesh->n_g_vertices + v_add_gnum[i];
    }

    v_add_io_num = fvm_io_num_destroy(v_add_io_num);

    mesh->n_g_vertices += n_g_add_vertices;
  }
  else
    mesh->n_g_vertices += (cs_gnum_t)n_add_vertices;

  /* Rebuild local vertices interface set if needed */

  if (rebuild_vtx_interfaces) {
    if (mesh->n_init_perio > 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s should be called before halo creation in case of periodicity."),
         __func__);
    mesh->vtx_interfaces = cs_interface_set_create(n_vertices + n_add_vertices,
                                                   NULL,
                                                   mesh->global_vtx_num,
                                                   NULL,
                                                   0,
                                                   NULL,
                                                   NULL,
                                                   NULL);
  }

  /* Now build indexed sets of matching local vertices */

  cs_lnum_t *v_new_vtx;
  BFT_MALLOC(v_new_vtx, v2c_idx[n_vertices], cs_lnum_t);

  {
    cs_lnum_t k = 0;

    for (cs_lnum_t i = 0; i < n_vertices; i++) {
      if (v2c_idx[i+1] > v2c_idx[i]) {
        cs_gnum_t g_vtx_num = (mesh->global_vtx_num != NULL) ?
          mesh->global_vtx_num[i] : (cs_gnum_t)i + 1;
        cs_lnum_t s_id = v2c_idx[i];
        cs_lnum_t e_id = v2c_idx[i+1];
        for (cs_lnum_t j = s_id; j < e_id; j++) {
          cs_gnum_t v_g_min_cmp = v_g_min_c_num[j];
          if (c_num_min[i] == v_g_min_cmp)
            v_new_vtx[j] = i;
          else {
            while (v_g_min_c_add[2*k] != g_vtx_num && k < n_add_vertices)
              k++;
            assert(k < n_add_vertices);
            cs_lnum_t ks = k;
            while (v_g_min_c_add[2*ks+1] != v_g_min_cmp && ks < n_add_vertices)
              ks++;
            assert(ks < n_add_vertices);
            assert(v_g_min_c_add[2*ks] == g_vtx_num);
            v_new_vtx[j] = n_vertices + ks;
          }
        }
      }
    }
  }

  BFT_FREE(v_g_min_c_num);
  BFT_FREE(v_g_min_c_add);
  BFT_FREE(c_num_min);

  /* Update face -> vertex ids based on face adjacencies */

  for (cs_lnum_t f_id = 0; f_id < mesh->n_i_faces; f_id++) {
    cs_lnum_t s_id = mesh->i_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->i_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->i_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0) {
        cs_lnum_t c_id = mesh->i_face_cells[f_id][0];
        if (c_id < 0)
          c_id = mesh->i_face_cells[f_id][1];
        cs_lnum_t v_s_id = v2c_idx[vtx_id];
        cs_lnum_t v_e_id = v2c_idx[vtx_id+1];
        cs_lnum_t j;
        for (j = v_s_id; j < v_e_id; j++) {
          if (v2c[j] == c_id) {
            mesh->i_face_vtx_lst[i] = v_new_vtx[j];
            break;
          }
        }
        assert(j < v_e_id);
      }
    }
  }

  for (cs_lnum_t f_id = 0; f_id < mesh->n_b_faces; f_id++) {
    cs_lnum_t s_id = mesh->b_face_vtx_idx[f_id];
    cs_lnum_t e_id = mesh->b_face_vtx_idx[f_id+1];
    for (cs_lnum_t i = s_id; i < e_id; i++) {
      cs_lnum_t vtx_id = mesh->b_face_vtx_lst[i];
      if (v_flag[vtx_id] != 0) {
        cs_lnum_t c_id = mesh->b_face_cells[f_id];
        cs_lnum_t v_s_id = v2c_idx[vtx_id];
        cs_lnum_t v_e_id = v2c_idx[vtx_id+1];
        for (cs_lnum_t j = v_s_id; j < v_e_id; j++) {
          if (v2c[j] == c_id) {
            mesh->b_face_vtx_lst[i] = v_new_vtx[j];
            break;
          }
        }
      }
    }
  }

  BFT_FREE(v_new_vtx);

  BFT_FREE(v2c);
  BFT_FREE(v2c_idx);

  BFT_FREE(v_flag);
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
               "   Number of cells:          %llu\n"
               "   Number of interior faces: %llu\n"
               "   Number of boundary faces: %llu\n"
               "   Number of vertices:       %llu\n\n"),
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

  if (*face_ifs != NULL) {
    if (mesh->periodicity != NULL)
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
 * \brief Insert boundary into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \param[in, out]  mesh            pointer to mesh structure to modify
 * \param[in, out]  face_ifs        pointer to interface set for faces
 * \param[in]       split_vertices  should we also split vertices ?
 * \param[in]       n_faces         number of selected (interior) faces
 * \param[in, out]  ids             list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

static void
_boundary_insert(cs_mesh_t           *mesh,
                 cs_interface_set_t  *face_ifs,
                 bool                 split_vertices,
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

  _print_mesh_info(mesh, _("Before insertion of user-defined boundary"));

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

  if (mesh->n_g_b_faces != _n_g_b_faces)
    mesh->modified = 1;

  mesh->n_g_b_faces = _n_g_b_faces;
  mesh->n_g_i_faces = _n_g_i_faces;

  if (cs_glob_n_ranks > 1)
    BFT_FREE(face_id_c);

  /* Separate vertices on each side of wall */

  if (split_vertices)
    _split_vertices(mesh, n_b_faces, mesh->n_b_faces);

  /* Now log new mesh counts */

  _print_mesh_info(mesh, _("After insertion of user-defined boundary"));
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundary into the mesh.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \param[in]       mesh     pointer to mesh structure to modify
 * \param[in]       n_faces  number of selected (interior) faces
 * \param[in, out]  face_id  list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert(cs_mesh_t  *mesh,
                        cs_lnum_t   n_faces,
                        cs_lnum_t   face_id[])
{
  cs_interface_set_t *face_ifs = _build_faces_interface_set_if_needed(mesh);

  _boundary_insert(mesh, NULL, true, n_faces, face_id);

  _destroy_faces_interface_set(mesh, &face_ifs);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert boundary into the mesh, sharing vertices on both sides.
 *
 * This is done by transforming interior faces into boundary faces.
 *
 * Note that the ids of selected faces are sorted by this function.
 *
 * \deprecated Use of this function is not recommended, as sharing vertices
 *             may cause issues with vertex-based values, gradient extended
 *             neighborhoods, and some visualization operations.

 * \param[in]       mesh     pointer to mesh structure to modify
 * \param[in]       n_faces  number of selected (interior) faces
 * \param[in, out]  face_id  list of selected (interior) faces (0 to n-1)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_boundary_insert_with_shared_vertices(cs_mesh_t  *mesh,
                                             cs_lnum_t   n_faces,
                                             cs_lnum_t   face_id[])
{
  cs_interface_set_t *face_ifs = _build_faces_interface_set_if_needed(mesh);

  _boundary_insert(mesh, NULL, false, n_faces, face_id);

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

  cs_lnum_t  n_b_faces_old = mesh->n_b_faces;

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

  _boundary_insert(mesh, face_ifs, true, b_count, face_tag);

  _destroy_faces_interface_set(mesh, &face_ifs);

  BFT_FREE(face_tag);

  /* Add group name if requested */

  if (group_name != NULL) {

    cs_lnum_t *sel_faces;
    cs_lnum_t  n_sel = mesh->n_b_faces - n_b_faces_old;
    BFT_MALLOC(sel_faces, n_sel, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_sel; i++)
      sel_faces[i] = n_b_faces_old + i;

    cs_mesh_group_b_faces_add(mesh,
                              group_name,
                              n_sel,
                              sel_faces);

    BFT_FREE(sel_faces);

  }
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
