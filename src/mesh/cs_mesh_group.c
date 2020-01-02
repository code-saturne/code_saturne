/*============================================================================
 * Management of mesh groups.
 *===========================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_periodicity.h"

#include "cs_order.h"
#include "cs_search.h"
#include "cs_sort.h"
#include "cs_join_perio.h"
#include "cs_join_post.h"
#include "cs_join_util.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *---------------------------------------------------------------------------*/

#include "cs_mesh_group.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *===========================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *===========================================================================*/

/*============================================================================
 * Private function definitions
 *===========================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compare 2 family definitions
 *
 * \param[in]  n0    number of values in family 0
 * \param[in]  n1    number of values in family 1
 * \param[in]  val0  family 0 values
 * \param[in]  val1  family 1 values
 *
 * \return  -1 if family 0 is smaller than family 1, 0 if families are equal,
 *   1 otherwise
 */
/*----------------------------------------------------------------------------*/

inline static int
_sync_compare_families(int        n0,
                       int        n1,
                       const int  val0[],
                       const int  val1[])
{
  int i;
  int retval = 0;

  for (i = 0; i < n0 && i < n1; i++) {
    if (val0[i] < val1[i]) {
      retval = -1;
      break;
    }
    else if (val0[i] > val1[i]) {
      retval = 1;
      break;
    }
  }

  if (retval == 0) {
    if (n0 < n1)
      retval = -1;
    else if (n0 > n1)
      retval = 1;
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize family combinations globally.
 *
 * \param[in, out]  n_fam       number of families
 * \param[in, out]  family_idx  index of element families
 * \param[in, out]  family      element family ancestors
 *
 * \return  renumbering from previous to synchronized family combinations
 */
/*----------------------------------------------------------------------------*/

static int *
_sync_family_combinations(int   *n_fam,
                          int  **family_idx,
                          int  **family)
{
  int i;
  int sizes[2];
  int rank_mult, rank_id, n_ranks;

  int _n_fam = *n_fam;
  int *_family_idx = NULL;
  int *_family = NULL;

  int *renum = NULL;

  BFT_MALLOC(_family_idx, _n_fam+1, int);
  memcpy(_family_idx, *family_idx, (_n_fam+1)*sizeof(int));

  BFT_MALLOC(_family, _family_idx[_n_fam], int);
  if (_family != NULL)
    memcpy(_family, *family, _family_idx[_n_fam]*sizeof(int));

  /* Merge all info by stages up to rank 0 */

  for (rank_id = cs_glob_rank_id, n_ranks = cs_glob_n_ranks, rank_mult = 1;
       n_ranks > 1;
       rank_id /= 2, n_ranks = (n_ranks+1)/2, rank_mult *= 2) {

    /* Even ranks receive and merge */

    if (rank_id %2 == 0 && rank_id + 1 < n_ranks) {

      MPI_Status status;
      int recv_rank = (rank_id + 1)*rank_mult;

      MPI_Recv(sizes, 2, MPI_INT, recv_rank, 1, cs_glob_mpi_comm, &status);

      if (sizes[0] > 0) {

        int j, k;
        int n_fam_new = 0;
        int *idx_recv = NULL, *fam_recv = NULL;
        int *idx_new = NULL, *fam_new = NULL;

        BFT_MALLOC(idx_recv, sizes[0] + 1, int);
        MPI_Recv(idx_recv, sizes[0] + 1, MPI_INT, recv_rank, 2,
                 cs_glob_mpi_comm, &status);
        BFT_MALLOC(fam_recv, sizes[1], int);
        MPI_Recv(fam_recv, sizes[1], MPI_INT, recv_rank, 3,
                 cs_glob_mpi_comm, &status);

        /* Merge data */

        BFT_MALLOC(idx_new, _n_fam + sizes[0] + 1, int);
        BFT_MALLOC(fam_new, _family_idx[_n_fam] + sizes[1], int);

        idx_new[0] = 0;
        i = 0; j = 0;

        while (i < _n_fam && j < sizes[0]) {
          int cmp;
          int n0 = _family_idx[i+1] - _family_idx[i];
          int n1 = idx_recv[j+1] - idx_recv[j];
          cmp = _sync_compare_families(n0,
                                       n1,
                                       _family + _family_idx[i],
                                       fam_recv + idx_recv[j]);
          if (cmp <= 0) {
            for (k = 0; k < n0; k++)
              fam_new[idx_new[n_fam_new] + k] = _family[_family_idx[i] + k];
            idx_new[n_fam_new + 1] = idx_new[n_fam_new] + n0;
            n_fam_new += 1;
            i += 1;
            if (cmp == 0)
              j += 1;
          }
          else if (cmp > 0) {
            for (k = 0; k < n1; k++)
              fam_new[idx_new[n_fam_new] + k] = fam_recv[idx_recv[j] + k];
            idx_new[n_fam_new + 1] = idx_new[n_fam_new] + n1;
            n_fam_new += 1;
            j += 1;
          }
        }

        while (i < _n_fam) {
          int n0 = _family_idx[i+1] - _family_idx[i];
          for (k = 0; k < n0; k++)
            fam_new[idx_new[n_fam_new] + k] = _family[_family_idx[i] + k];
          idx_new[n_fam_new + 1] = idx_new[n_fam_new] + n0;
          n_fam_new += 1;
          i += 1;
        }
        while (j < sizes[0]) {
          int n1 = idx_recv[j+1] - idx_recv[j];
          for (k = 0; k < n1; k++)
            fam_new[idx_new[n_fam_new] + k] = fam_recv[idx_recv[j] + k];
          idx_new[n_fam_new + 1] = idx_new[n_fam_new] + n1;
          n_fam_new += 1;
          j += 1;
        }

        BFT_FREE(fam_recv);
        BFT_FREE(idx_recv);

        BFT_REALLOC(idx_new, n_fam_new + 1, int);
        BFT_REALLOC(fam_new, idx_new[n_fam_new], int);

        BFT_FREE(_family_idx);
        BFT_FREE(_family);

        _family_idx = idx_new;
        _family = fam_new;
        _n_fam = n_fam_new;

      } /* if (sizes[0] > 0) */

    }

    /* Odd ranks send once, then are finished for the merge step */

    else if (rank_id % 2 == 1) {

      int send_rank = (rank_id-1)*rank_mult;

      sizes[0] = _n_fam;
      sizes[1] = _family_idx[_n_fam];
      MPI_Send(sizes, 2, MPI_INT, send_rank, 1, cs_glob_mpi_comm);

      if (sizes[0] > 0) {
        MPI_Send(_family_idx, sizes[0] + 1, MPI_INT, send_rank, 2,
                 cs_glob_mpi_comm);
        MPI_Send(_family, sizes[1], MPI_INT, send_rank, 3, cs_glob_mpi_comm);
      }

      break;

    }

  }

  /* Now rank 0 broadcasts */

  sizes[0] = _n_fam;
  sizes[1] = _family_idx[_n_fam];
  MPI_Bcast(sizes, 2, MPI_INT, 0, cs_glob_mpi_comm);

  _n_fam = sizes[0];

  if (cs_glob_rank_id != 0) {
    BFT_REALLOC(_family_idx, sizes[0] + 1, int);
    BFT_REALLOC(_family, sizes[1], int);
  }

  MPI_Bcast(_family_idx, sizes[0] + 1, MPI_INT, 0, cs_glob_mpi_comm);
  MPI_Bcast(_family, sizes[1], MPI_INT, 0, cs_glob_mpi_comm);

  /* Finally generate renumbering array */

  BFT_MALLOC(renum, *n_fam, int);

  for (i = 0; i < *n_fam; i++) {

    int start_id, end_id, mid_id;
    int cmp_ret = 1;
    int n1 = (*family_idx)[i+1] - (*family_idx)[i];

    /* Use binary search to find entry */

    start_id = 0;
    end_id = _n_fam - 1;
    mid_id = start_id + ((end_id -start_id) / 2);

    while (start_id <= end_id) {
      int n0 = _family_idx[mid_id+1] - _family_idx[mid_id];
      cmp_ret = _sync_compare_families(n0,
                                       n1,
                                       _family + _family_idx[mid_id],
                                       (*family) + (*family_idx)[i]);
      if (cmp_ret < 0)
        start_id = mid_id + 1;
      else if (cmp_ret > 0)
        end_id = mid_id - 1;
      else
        break;
      mid_id = start_id + ((end_id -start_id) / 2);
    }

    assert(cmp_ret == 0);

    renum[i] = mid_id;
  }

  BFT_FREE(*family_idx);
  BFT_FREE(*family);

  *n_fam = _n_fam;
  *family_idx = _family_idx;
  *family = _family;

  return renum;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Descend binary tree for the ordering of a mesh's groups.
 *
 * \param[in, out]  mesh   pointer to mesh structure
 * \param[in]       level  level of the binary tree to descend
 * \param[in]       n      number of groups in the binary tree to descend
 * \param[in, out]  order  pre-allocated ordering table
 */
/*----------------------------------------------------------------------------*/

inline static void
_groups_descend_tree(const cs_mesh_t  *mesh,
                     size_t            level,
                     const size_t      n,
                     int               order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (n/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (strcmp(mesh->group + mesh->group_idx[i1],
                 mesh->group + mesh->group_idx[i2]) > 0)
        lv_cur++;
    }

    if (lv_cur >= n) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (strcmp(mesh->group + mesh->group_idx[i1],
               mesh->group + mesh->group_idx[i2]) >= 0)
      break;

    order[level] = order[lv_cur];
    level = lv_cur;
  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Order mesh groups.
 *
 * \param[in, out]  mesh   pointer to mesh structure
 * \param[out]      order  pre-allocated ordering table
 */
/*----------------------------------------------------------------------------*/

static void
_order_groups(const cs_mesh_t  *mesh,
              int               order[])
{
  int    o_save;
  size_t i;
  size_t n = mesh->n_groups;

  /* Initialize ordering array */

  for (i = 0; i < n; i++)
    order[i] = i;

  if (n < 2)
    return;

  /* Create binary tree */

  i = (n / 2);
  do {
    i--;
    _groups_descend_tree(mesh, i, n, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _groups_descend_tree(mesh, 0, i, order);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add mesh group definition if needed.
 *
 * \param[in, out]  mesh  pointer to mesh structure to modify
 * \param[in]       name  group name to add
 *
 * \return  id of group
 */
/*----------------------------------------------------------------------------*/

static int
_add_group(cs_mesh_t   *mesh,
           const char  *name)
{
  /* Check if already present */

  for (int i = 0; i < mesh->n_groups; i++) {
    const char *g_cur = mesh->group + mesh->group_idx[i];
    int d = strcmp(g_cur, name);
      if (d == 0)
        return i;
  }

  /* Add in case of empty groups */

  if (mesh->n_groups == 0) {
    mesh->n_groups = 1;
    BFT_MALLOC(mesh->group_idx, 2, int);
    mesh->group_idx[0] = 0;
    mesh->group_idx[1] = strlen(name) + 1;
    BFT_MALLOC(mesh->group, mesh->group_idx[1], char);
    strcpy(mesh->group, name);
    return 0;
  }

  /* Add to end then renumber */

  else {

    size_t l = strlen(name) + 1;
    int n_groups_o = mesh->n_groups;
    mesh->n_groups += 1;
    BFT_REALLOC(mesh->group_idx, mesh->n_groups + 1, int);
    BFT_REALLOC(mesh->group,
                mesh->group_idx[n_groups_o] + l,
                char);
    strcpy(mesh->group + mesh->group_idx[n_groups_o], name);
    mesh->group_idx[mesh->n_groups] = mesh->group_idx[n_groups_o] + l;
    cs_mesh_group_clean(mesh);

    for (int i = 0; i < mesh->n_groups; i++) {
      const char *g_cur = mesh->group + mesh->group_idx[i];
      int d = strcmp(g_cur, name);
      if (d == 0)
        return i;
    }

  }

  /* We should not arrive here */

  assert(0);
  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a group class for a group with a given name.
 *
 * If a group class with that group an no other is already present,
 * its id is returned.
 *
 * \param[in, out]  mesh   pointer to mesh structure to modify
 * \param[in]       name   group name to add
 *
 * \return          id of matching or added group class
 */
/*----------------------------------------------------------------------------*/

static int
_add_gc(cs_mesh_t   *mesh,
        const char  *name)
{
  int g_id = _add_group(mesh, name);
  int g_id_cmp = -g_id - 1;

  /* Check if a group class definition with only the selected group
     exists */

  int gc_id = -1;
  int gc_id_shift = 0;

  if (mesh->n_max_family_items > 1) {
    for (int i = 0; i < mesh->n_families; i++) {
      if (   mesh->family_item[i] == g_id_cmp
          && mesh->family_item[mesh->n_families + i] == 0) {
        gc_id = i;
        break;
      }
    }
  }

  else if (mesh->n_max_family_items == 1) {
    for (int i = 0; i < mesh->n_families; i++) {
      if (mesh->family_item[i] == g_id_cmp) {
        gc_id = i;
        break;
      }
    }
  }

  /* Add a family if needed */

  if (gc_id < 0) {

    int n_f_prv = mesh->n_families;
    int *f_prv = NULL;

    if (n_f_prv*mesh->n_max_family_items > 0) {
      BFT_MALLOC(f_prv, n_f_prv * mesh->n_max_family_items, int);
      memcpy(f_prv,
             mesh->family_item,
             (n_f_prv * mesh->n_max_family_items)*sizeof(int));
    }

    gc_id = mesh->n_families + gc_id_shift;

    mesh->n_families += 1;

    if (mesh->n_max_family_items == 0) {
      mesh->n_max_family_items = 1;
      BFT_REALLOC(mesh->family_item,
                  mesh->n_families*mesh->n_max_family_items,
                  int);
      for (int i = 0; i < mesh->n_families; i++)
        mesh->family_item[i] = 0;
      mesh->family_item[gc_id] = g_id_cmp;
    }
    else {
      BFT_REALLOC(mesh->family_item,
                  mesh->n_families*mesh->n_max_family_items,
                  int);
      for (int j = 0; j < mesh->n_max_family_items; j++) {
        for (int i = 0; i < n_f_prv; i++)
          mesh->family_item[mesh->n_families*j+i] = f_prv[n_f_prv*j + i];
      }
      mesh->family_item[gc_id - gc_id_shift] = g_id_cmp;
      for (int j = 1; j < mesh->n_max_family_items; j++)
        mesh->family_item[mesh->n_families*j + gc_id - gc_id_shift] = 0;
    }

    BFT_FREE(f_prv);
  }

  return gc_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to mesh entities, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh             pointer to mesh structure to modify
 * \param[in]       name             group name to assign to selected elements
 * \param[in]       n_selected_elts  number of selected elements
 * \param[in]       selected_elt_id  selected element ids (size: n_selected_elts)
 * \param[in, out]  gc_id            element group class ids (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

static void
_mesh_group_set(cs_mesh_t        *mesh,
                const char       *name,
                cs_lnum_t         n_selected_elts,
                const cs_lnum_t   selected_elt_id[],
                int               gc_id[])
{
  int _gc_id = _add_gc(mesh, name);

  for (cs_lnum_t i = 0; i < n_selected_elts; i++)
    gc_id[selected_elt_id[i]] = _gc_id + 1;

  if (mesh->class_defs != NULL)
    cs_mesh_update_selectors(mesh);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected mesh entities to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh             pointer to mesh structure to modify
 * \param[in]       name             group name to assign to selected elements
 * \param[in]       n_elts           number of elements in mesh entity
 * \param[in]       n_selected_elts  number of selected elements
 * \param[in]       selected_elt_id  selected element ids (size: n_selected_elts)
 * \param[in, out]  gc_id            element group class ids (size: n_elts)
 */
/*----------------------------------------------------------------------------*/

static void
_mesh_group_add(cs_mesh_t        *mesh,
                const char       *name,
                cs_lnum_t         n_elts,
                cs_lnum_t         n_selected_elts,
                const cs_lnum_t   selected_elt_id[],
                int               gc_id[])
{
  int _gc_id = _add_gc(mesh, name);

  int null_family = 0;
  if (mesh->n_families > 0) {
    if (mesh->family_item[0] == 0)
      null_family = 1;
  }

  /* Build index on entities (previous group class for elements
     not selected, previous + new for those selected) */

  cs_lnum_t *gc_tmp_idx = NULL, *gc_tmp = NULL;

  BFT_MALLOC(gc_tmp_idx, n_elts + 1, cs_lnum_t);
  gc_tmp_idx[0] = 0;

  for (cs_lnum_t i = 0; i < n_elts; i++)
    gc_tmp_idx[i+1] = 1;
  for (cs_lnum_t i = 0; i < n_selected_elts; i++) {
    cs_lnum_t j= selected_elt_id[i];
    if (gc_id[j] != null_family)
      gc_tmp_idx[j+1] += 1;
  }

  /* Transform count to index */

  for (cs_lnum_t i = 0; i < n_elts; i++)
    gc_tmp_idx[i+1] += gc_tmp_idx[i];

  /* Now assign multiple group classes */

  BFT_MALLOC(gc_tmp, gc_tmp_idx[n_elts], int);
  for (cs_lnum_t i = 0; i < n_elts; i++)
    gc_tmp[gc_tmp_idx[i]] = gc_id[i];
  for (cs_lnum_t i = 0; i < n_selected_elts; i++) {
    cs_lnum_t j = selected_elt_id[i];
    if (gc_id[j] != null_family)
      gc_tmp[gc_tmp_idx[j] + 1] = _gc_id + 1;
    else
      gc_tmp[gc_tmp_idx[j]] = _gc_id + 1;
  }

  /* Merge definitions */

  cs_mesh_group_combine_classes(mesh, n_elts, gc_tmp_idx, gc_tmp, gc_id);

  /* Cleanup */

  BFT_FREE(gc_tmp_idx);
  BFT_FREE(gc_tmp);

  if (mesh->class_defs != NULL)
    cs_mesh_update_selectors(mesh);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clean mesh group definitions.
 *
 * \param[in]  mesh  pointer to mesh structure to modify
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_clean(cs_mesh_t  *mesh)
{
  int i;
  size_t j;
  int n_groups = 0;
  size_t size_tot = 0;
  char *g_prev = NULL, *g_cur = NULL, *g_lst = NULL;
  int *order = NULL, *renum = NULL;

  if (mesh->n_groups < 1)
    return;

  /* Order group names */

  BFT_MALLOC(renum, mesh->n_groups, int);
  BFT_MALLOC(order, mesh->n_groups, int);

  _order_groups(mesh, order);

  /* Build compact group copy */

  BFT_MALLOC(g_lst, mesh->group_idx[mesh->n_groups], char);

  g_cur = mesh->group + mesh->group_idx[order[0]];
  g_prev = g_cur;
  n_groups += 1;
  strcpy(g_lst, g_cur);
  size_tot += strlen(g_cur) + 1;
  g_lst[size_tot - 1] = '\0';
  renum[order[0]] = 0;

  for (i = 1; i < mesh->n_groups; i++) {
    g_cur = mesh->group + mesh->group_idx[order[i]];
    if (strcmp(g_cur, g_prev) != 0) {
      g_prev = g_cur;
      strcpy(g_lst + size_tot, g_cur);
      n_groups += 1;
      size_tot += strlen(g_cur) + 1;
      g_lst[size_tot - 1] = '\0';
    }
    renum[order[i]] = n_groups - 1;
  }

  BFT_FREE(order);

  BFT_REALLOC(mesh->group_idx, n_groups + 1, cs_int_t);
  BFT_REALLOC(mesh->group, size_tot, char);

  mesh->n_groups = n_groups;
  memcpy(mesh->group, g_lst, size_tot);

  mesh->group_idx[0] = 0;
  for (i = 0; i < mesh->n_groups; i++) {
    j = strlen(mesh->group + mesh->group_idx[i]) + 1;
    mesh->group_idx[i + 1] = mesh->group_idx[i] + j;
  }

  BFT_FREE(g_lst);

  /* Now renumber groups in group class description */

  size_tot = mesh->n_families * mesh->n_max_family_items;

  for (j = 0; j < size_tot; j++) {
    int gc_id = mesh->family_item[j];
    if (gc_id < 0)
      mesh->family_item[j] = - renum[-gc_id - 1] - 1;
  }

  BFT_FREE(renum);

  /* Remove empty group if present (it should appear first now) */

  if (mesh->n_groups > 1) {

    if ((mesh->group_idx[1] - mesh->group_idx[0]) == 1) {

      size_t new_lst_size = (  mesh->group_idx[mesh->n_groups]
                             - mesh->group_idx[1]);
      for (i = 0; i < mesh->n_groups; i++)
        mesh->group_idx[i] = mesh->group_idx[i+1] - 1;
      mesh->n_groups -= 1;
      memmove(mesh->group, mesh->group+1, new_lst_size);

      BFT_REALLOC(mesh->group_idx, mesh->n_groups + 1, int);
      BFT_REALLOC(mesh->group, new_lst_size, char);

      for (j = 0; j < size_tot; j++) {
        int gc_id = mesh->family_item[j];
        if (gc_id < 0)
          mesh->family_item[j] += 1;
      }

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Combine mesh group classes.
 *
 * \param[in, out]  mesh          pointer to mesh structure to modify
 * \param[in]       n_elts        number of local elements
 * \param[in]       gc_id_idx     element group class index (size: n_elts +1)
 * \param[in]       gc_id         input element group classes
 *                                (size: gc_id_idx[n_elts])
 * \param[in]       gc_id_merged  output element group classes (size: n_elts)
 *
 * \return  array of new element group class ids
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_combine_classes(cs_mesh_t   *mesh,
                              cs_lnum_t    n_elts,
                              cs_lnum_t    gc_id_idx[],
                              int          gc_id[],
                              int          gc_id_merged[])
{
  int   n_fam = 0;
  int  *_gc_id_idx = NULL;
  int  *_gc_id = NULL;

  /* First, build families locally */

  if (n_elts > 0) {

    cs_lnum_t i, j, j_prev, n_prev;
    cs_lnum_t *order = NULL;
    cs_gnum_t *tmp_gc_id = NULL;
    const cs_lnum_t n_gc_id_values = gc_id_idx[n_elts];

    /* Build ordering of elements by associated families */

    BFT_MALLOC(tmp_gc_id, n_gc_id_values, cs_gnum_t);

    for (i = 0; i < n_gc_id_values; i++)
      tmp_gc_id[i] = gc_id[i] + 1;

    order = cs_order_gnum_i(NULL, tmp_gc_id, gc_id_idx, n_elts);

    BFT_FREE(tmp_gc_id);

    /* Build new family array, merging identical definitions */

    BFT_MALLOC(_gc_id_idx, n_elts + 1, int);
    BFT_MALLOC(_gc_id, gc_id_idx[n_elts], int);

    _gc_id_idx[0] = 0;

    j_prev = -1;
    n_prev = -1;

    for (i = 0; i < n_elts; i++) {
      cs_lnum_t k, l, n;
      bool is_same = true;
      j = order[i];
      n = gc_id_idx[j+1] - gc_id_idx[j];
      if (n != n_prev)
        is_same = false;
      else {
        for (k = gc_id_idx[j], l = gc_id_idx[j_prev];
             k < gc_id_idx[j + 1];
             k++, l++) {
          if (gc_id[k] != gc_id[l])
            is_same = false;
        }
      }
      if (is_same)
        gc_id_merged[j] = gc_id_merged[j_prev];
      else if (n == 0)
        gc_id_merged[j] = 0;
      else if (n == 1)
        gc_id_merged[j] = gc_id[gc_id_idx[j]];
      else {
        gc_id_merged[j] = mesh->n_families + 1 + n_fam;
        for (k = gc_id_idx[j], l = _gc_id_idx[n_fam];
             k < gc_id_idx[j + 1];
             k++, l++)
          _gc_id[l] = gc_id[k];
        _gc_id_idx[n_fam+1] = _gc_id_idx[n_fam] + n;
        n_fam += 1;
      }
      j_prev = j;
      n_prev = n;
    }

    BFT_FREE(order);

    BFT_REALLOC(_gc_id_idx, n_fam + 1, int);
    BFT_REALLOC(_gc_id, _gc_id_idx[n_fam], int);

  }
  else {
    BFT_MALLOC(_gc_id_idx, 1, int);
    _gc_id_idx[0] = 0;
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    cs_lnum_t i;
    int *renum = _sync_family_combinations(&n_fam,
                                           &_gc_id_idx,
                                           &_gc_id);

    for (i = 0; i < n_elts; i++) {
      int j = gc_id_merged[i] - mesh->n_families - 1;
      if (j >= 0)
        gc_id_merged[i] = mesh->n_families + renum[j] + 1;
    }

    BFT_FREE(renum);
  }

#endif

  /* Update mesh definition */

  {
    int i, j, k, l;
    int n_max_family_items = 0;

    for (i = 0; i < n_fam; i++) {
      int n_family_items = 0;
      for (j = _gc_id_idx[i]; j < _gc_id_idx[i+1]; j++) {
        int f_id = _gc_id[j] - 1;
        for (k = 0; k < mesh->n_max_family_items; k++) {
          if (mesh->family_item[mesh->n_families*k + f_id] != 0)
            n_family_items++;
        }
      }
      if (n_family_items > n_max_family_items)
        n_max_family_items = n_family_items;
    }

    /* Increase maximum number of definitions and pad it necessary */

    if (n_max_family_items > mesh->n_max_family_items) {
      BFT_REALLOC(mesh->family_item,
                  mesh->n_families*n_max_family_items,
                  cs_lnum_t);
      for (i = mesh->n_max_family_items;
           i < n_max_family_items;
           i++) {
        for (j = 0; j < mesh->n_families; j++)
          mesh->family_item[mesh->n_families*i + j] = 0;
      }
      mesh->n_max_family_items = n_max_family_items;
    }

    /* Increase number of families */

    mesh->n_families += n_fam;

    BFT_REALLOC(mesh->family_item,
                mesh->n_families * mesh->n_max_family_items,
                int);
    for (j = mesh->n_max_family_items - 1; j > 0; j--) {
      for (i = mesh->n_families - n_fam - 1; i > -1; i--)
        mesh->family_item[mesh->n_families*j + i]
          = mesh->family_item[(mesh->n_families - n_fam)*j + i];
    }
    for (i = mesh->n_families - n_fam, j = 0; i < mesh->n_families; i++, j++) {
      int n_family_items = 0;
      for (k = _gc_id_idx[j]; k < _gc_id_idx[j+1]; k++) {
        int f_id = _gc_id[k] - 1;
        for (l = 0; l < mesh->n_max_family_items; l++) {
          if (mesh->family_item[mesh->n_families*l + f_id] != 0) {
            mesh->family_item[mesh->n_families*n_family_items + i]
              = mesh->family_item[mesh->n_families*l + f_id];
            n_family_items++;
          }
        }
      }
      for (k = n_family_items; k < mesh->n_max_family_items; k++)
        mesh->family_item[mesh->n_families*k + i] = 0;
    }

  }

  BFT_FREE(_gc_id_idx);
  BFT_FREE(_gc_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to cells, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected cells
 * \param[in]       n_selected_cells  number of selected cells
 * \param[in]       selected_cell_id  selected cell ids (size: n_selected_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_cells_set(cs_mesh_t        *mesh,
                        const char       *name,
                        cs_lnum_t         n_selected_cells,
                        const cs_lnum_t   selected_cell_id[])
{
  _mesh_group_set(mesh,
                  name,
                  n_selected_cells,
                  selected_cell_id,
                  mesh->cell_family);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to interior faces, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_i_faces_set(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[])
{
  _mesh_group_set(mesh,
                  name,
                  n_selected_faces,
                  selected_face_id,
                  mesh->i_face_family);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a given group to boundary faces, removing those entities
 *        from previous groups if present.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_b_faces_set(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[])
{
  _mesh_group_set(mesh,
                  name,
                  n_selected_faces,
                  selected_face_id,
                  mesh->b_face_family);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected cells to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected cells
 * \param[in]       n_selected_cells  number of selected cells
 * \param[in]       selected_cell_id  selected cell ids (size: n_selected_cells)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_cells_add(cs_mesh_t        *mesh,
                        const char       *name,
                        cs_lnum_t         n_selected_cells,
                        const cs_lnum_t   selected_cell_id[])
{
  _mesh_group_add(mesh,
                  name,
                  mesh->n_cells,
                  n_selected_cells,
                  selected_cell_id,
                  mesh->cell_family);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected interior faces to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_i_faces_add(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[])
{
  _mesh_group_add(mesh,
                  name,
                  mesh->n_i_faces,
                  n_selected_faces,
                  selected_face_id,
                  mesh->i_face_family);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add selected boundary faces to a given group.
 *
 * The group is created if necessary.
 *
 * \param[in, out]  mesh              pointer to mesh structure to modify
 * \param[in]       name              group name to assign to selected faces
 * \param[in]       n_selected_faces  number of selected faces
 * \param[in]       selected_face_id  selected face ids (size: n_selected_faces)
 */
/*----------------------------------------------------------------------------*/

void
cs_mesh_group_b_faces_add(cs_mesh_t        *mesh,
                          const char       *name,
                          cs_lnum_t         n_selected_faces,
                          const cs_lnum_t   selected_face_id[])
{
  _mesh_group_add(mesh,
                  name,
                  mesh->n_b_faces,
                  n_selected_faces,
                  selected_face_id,
                  mesh->b_face_family);
}

/*---------------------------------------------------------------------------*/

END_C_DECLS
