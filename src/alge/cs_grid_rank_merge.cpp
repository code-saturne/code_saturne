/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(__STDC_VERSION__)      /* size_t */
#if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#include <stdlib.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_dispatch.h"
#include "base/cs_halo.h"
#include "base/cs_halo_perio.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "base/cs_reducers.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_matrix_tuning.h"
#include "alge/cs_matrix_util.h"
#include "base/cs_order.h"
#include "base/cs_parall.h"
#include "base/cs_reducers.h"
#include "alge/cs_sles.h"
#include "base/cs_sort.h"

#include "fvm/fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_grid.h"
#include "alge/cs_grid_priv.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file cs_grid.cpp
        Coarse grid construction for multigrid solver.

  The multigrid solver is based on an aggregation approach, with optional
  relaxation parameter.

  The idea behind the relaxation is that for an unstructured finite volume
  approach on arbitrary shaped cells, it is possible to improve the scaled
  Galerkin approach based on piecewise constant interpolation for the R0
  restriction and P0 prolongation operators. The off-diagonal entries of the
  P0 Galerkin coarse mesh matrix are rescaled by a parameter, which takes
  into account the mesh spacing ratio between the fine and coarse mesh in
  the vicinity of the coarse mesh cell boundaries.
  THis is detailed in \cite MeFoH09 .
*/

/*----------------------------------------------------------------------------*/

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
 * Return to pointer to root of a given grid.
 *
 * parameters:
 *   g  <-- Pointer to current grid
 *
 * return:
 *   pointer to grid's root in current hierarchy
 *----------------------------------------------------------------------------*/

static const cs_grid_t *
_root_grid(cs_grid_t  *g)
{
  const cs_grid_t  *r = g;

  while (r->parent != nullptr) {
    r = r->parent;
  }

  return r;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Empty a halo that is only useful for global synchronization
 *
 * parameters:
 *   h <-> pointer to halo structure being emptyied
 *----------------------------------------------------------------------------*/

static void
_empty_halo(cs_halo_t  *h)
{
  if (h == nullptr)
    return;

  h->n_c_domains = 0;
  CS_FREE(h->c_domain_rank);

  h->n_send_elts[0] = 0;
  h->n_send_elts[1] = 0;
  h->n_elts[0] = 0;
  h->n_elts[1] = 0;

  CS_FREE(h->send_list);
  CS_FREE(h->send_index);
  CS_FREE(h->send_perio_lst);
  CS_FREE(h->index);
  CS_FREE(h->perio_lst);
}

/*----------------------------------------------------------------------------
 * Append halo info for grid grouping and merging.
 *
 * parameters:
 *   g           <-> pointer to grid structure being merged
 *   new_cell_id <-> new cell ids for local cells
 *                   in: defined for local cells
 *                   out: updated also for halo cells
 *----------------------------------------------------------------------------*/

static void
_append_halos(cs_grid_t   *g,
              cs_lnum_t   *new_cell_id)
{
  int *recv_count = nullptr;
  cs_lnum_t *new_src_cell_id = nullptr;

  cs_halo_t *h = g->_halo;

  if (h == nullptr)
    return;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  const int n_transforms = h->n_transforms;
  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'h';

  if (g->merge_sub_size == 0) {
    _empty_halo(h);
    return;
  }

  int counts[3] = {h->n_c_domains, h->n_local_elts, h->n_elts[0]};

  /* local copy of halo sizes */

  int n_c_domains = h->n_c_domains;
  cs_lnum_t n_local_elts = h->n_local_elts;
  cs_lnum_t n_elts = h->n_elts[0];

  /* Exchange counters needed for concatenation */

  if (g->merge_sub_rank == 0) {
    CS_MALLOC(recv_count, g->merge_sub_size*3, int);
    for (int ii = 0; ii < 3; ii++)
      recv_count[ii] = counts[ii];
    for (int rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      MPI_Recv(recv_count + 3*rank_id, 3, MPI_INT, dist_rank,
               tag, comm, &status);
      for (int ii = 0; ii < 3; ii++)
        counts[ii] += recv_count[rank_id*3 + ii];
    }
  }
  else
    MPI_Send(counts, 3, MPI_INT, g->merge_sub_root, tag, comm);

  /* In case of periodic transforms, transpose perio list
     so as to have blocs of fixed n_transforms size, easier
     to work with for append + merge operations */

  int *perio_lst_tr = nullptr;
  if (n_transforms > 0) {
    CS_MALLOC(perio_lst_tr, counts[0]*n_transforms*2, cs_lnum_t);
    for (int rank_idx = 0; rank_idx < h->n_c_domains; rank_idx++) {
      for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
        for (int k = 0; k < 2; k++)
          perio_lst_tr[(n_transforms*rank_idx + tr_id)*2 + k]
            = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx + k];
      }
    }
    CS_FREE(h->perio_lst);
  }

  /* Adjust numbering */

  int *c_domain_rank;
  int n_domain_ranks = (g->merge_sub_rank == 0) ? counts[0] : h->n_c_domains;
  CS_MALLOC(c_domain_rank, n_domain_ranks, int);

  for (int rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    int init_rank_id = h->c_domain_rank[rank_id];
    c_domain_rank[rank_id]
      = init_rank_id - (init_rank_id % g->next_merge_stride);
  }

  if (g->merge_sub_rank == 0) {

    CS_MALLOC(new_src_cell_id, counts[2], cs_lnum_t);
    for (cs_lnum_t ii = g->n_rows, jj = 0;
         ii < g->n_cols_ext;
         ii++, jj++)
      new_src_cell_id[jj] = new_cell_id[ii];

    /* Reallocate index for receiving append data (local part is unmodified) */
    CS_REALLOC(h->index, counts[0]*2 + 1, cs_lnum_t);

    for (int rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      int n_c_domains_r = recv_count[rank_id*3 + 0];
      cs_lnum_t n_recv = recv_count[rank_id*3 + 2];

      cs_lnum_t index_shift = h->index[2*n_c_domains];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(c_domain_rank + n_c_domains, n_c_domains_r,
               MPI_INT, dist_rank, tag, comm, &status);
      MPI_Recv(new_src_cell_id + n_elts, n_recv,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);
      MPI_Recv(h->index + n_c_domains*2, n_c_domains_r*2+1,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);

      for (cs_lnum_t ii = 0, jj = n_c_domains*2;
           ii < n_c_domains_r*2+1;
           ii++, jj++)
        h->index[jj] += index_shift;

      if (n_transforms > 0) {
        MPI_Recv(perio_lst_tr + n_c_domains*n_transforms*2,
                 n_c_domains_r*n_transforms*2,
                 CS_MPI_LNUM, dist_rank, tag, comm, &status);
        for (cs_lnum_t ii = 0, jj = n_c_domains*n_transforms*2;
             ii < n_c_domains_r*n_transforms*2;
             ii+=2, jj+=2)
          perio_lst_tr[jj] += index_shift;
      }

      /* Update halo sizes */

      n_local_elts += recv_count[rank_id*3 + 1];
      n_c_domains += n_c_domains_r;
      n_elts += n_recv;
    }

  }
  else if (g->merge_sub_size > 1) {

    CS_MALLOC(new_src_cell_id, h->n_elts[0], cs_lnum_t);
    for (cs_lnum_t ii = g->n_rows, jj = 0; ii < g->n_cols_ext; ii++, jj++)
      new_src_cell_id[jj] = new_cell_id[ii];

    MPI_Send(c_domain_rank, h->n_c_domains, MPI_INT,
             g->merge_sub_root, tag, comm);
    MPI_Send(new_src_cell_id, h->n_elts[0], CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    MPI_Send(h->index, h->n_c_domains*2+1, CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);

    if (n_transforms > 0)
      MPI_Send(perio_lst_tr, h->n_c_domains*n_transforms*2,
               CS_MPI_LNUM, g->merge_sub_root, tag, comm);

  }

  /* Order list by rank, transform, and new element number */

  if (g->merge_sub_rank == 0) {

    const int  n_c_domains_ini = n_c_domains;
    const cs_lnum_t n_elts_ini = n_elts;

    int *elt_rank_idx;
    cs_lnum_t *elt_id;
    CS_MALLOC(elt_rank_idx, n_elts_ini, int);
    CS_MALLOC(elt_id, n_elts_ini, cs_lnum_t);

    int16_t *elt_tr_id = nullptr;
    if (n_transforms > 0)
      CS_MALLOC(elt_tr_id, n_elts_ini, int16_t);

    cs_lnum_t tr_mult = n_transforms+1;
    cs_lnum_t  *tmp_num = nullptr;
    CS_MALLOC(tmp_num, n_elts_ini*2, cs_lnum_t);

    for (int rank_idx = 0; rank_idx < n_c_domains_ini; rank_idx++) {
      int c_rank_id = c_domain_rank[rank_idx];
      /* use -2 instead of -1 for local rank so that division of
         c_rank_id*tr_mult + tr_id by tr_mult is < 0 */
      if (c_rank_id == cs_glob_rank_id) {
        c_rank_id = -2;
        c_domain_rank[rank_idx] = -1;
      }

      for (cs_lnum_t ii = h->index[rank_idx*2];
           ii < h->index[rank_idx*2+2];
           ii++)
        tmp_num[ii*2] = c_rank_id * tr_mult;
    }

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {
      tmp_num[ii*2 + 1] = new_src_cell_id[ii];
    }

    if (n_transforms > 0) {
      for (int rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
        for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
          cs_lnum_t s_id = perio_lst_tr[n_transforms*2*rank_id + 2*tr_id];
          cs_lnum_t e_id =   s_id
                           + perio_lst_tr[n_transforms*2*rank_id + 2*tr_id + 1];
          cs_lnum_t tr_num_shift = (tr_id + 1);
          for (cs_lnum_t ii = s_id; ii < e_id; ii++)
            tmp_num[ii*2] += tr_num_shift;
        }
      }
      CS_FREE(perio_lst_tr);
    }

    /* Now order elements based on rank, transform, and id */

    cs_lnum_t *order;
    CS_MALLOC(order, n_elts_ini*2, cs_lnum_t);

    cs_order_lnum_allocated_s(nullptr, tmp_num, 2, order, n_elts_ini);

    /* Now build ordered list and renumbering array */

    int prev_rank_tr_id = -1;
    cs_lnum_t prev_src_id = -1, n_elts_m = 0;

    cs_lnum_t *halo_o2n;
    CS_MALLOC(halo_o2n, n_elts_ini, int);

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {
      bool is_same = true;
      cs_lnum_t cur_id = order[ii];

      cs_lnum_t rank_tr_id = tmp_num[cur_id*2];
      int rank_id = rank_tr_id / tr_mult;
      int tr_id = (int16_t)(rank_tr_id % tr_mult) - 1;
      if (rank_tr_id < 0) {
        rank_id = -1;
        tr_id = (int16_t)((2*tr_mult+rank_tr_id) % tr_mult) - 1;
      }
      cs_lnum_t src_id = tmp_num[cur_id*2 + 1];

      /* If element is local, not in halo anymore */

      if (rank_id == -1 && tr_id == -1) {
        halo_o2n[cur_id] = src_id;
      }

      /* Otherwise, element in new halo */

      else {
        if (rank_tr_id != prev_rank_tr_id || src_id != prev_src_id)
          is_same = false;

        if (is_same == false) {
          elt_rank_idx[n_elts_m] = rank_id;
          elt_id[n_elts_m] = src_id;
          if (elt_tr_id != nullptr)
            elt_tr_id[n_elts_m] = tr_id;

          n_elts_m++;
        }

        halo_o2n[cur_id] = n_local_elts + n_elts_m - 1;

      }

      prev_rank_tr_id = rank_tr_id;
      prev_src_id = src_id;
    }

    /* Update halo section of cell renumbering array */

    cs_lnum_t n_send = recv_count[2];
    cs_lnum_t send_shift = n_send;

    for (cs_lnum_t ii = 0; ii < n_send; ii++)
      new_cell_id[g->n_rows + ii] = halo_o2n[ii];

    for (int rank_idx = 1; rank_idx < g->merge_sub_size; rank_idx++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_idx;
      int n_send_r = recv_count[rank_idx*3 + 2];
      MPI_Send(halo_o2n + send_shift, n_send_r, CS_MPI_LNUM,
               dist_rank, tag, comm);
      send_shift += n_send_r;
    }

    CS_FREE(halo_o2n);
    CS_FREE(order);
    CS_FREE(tmp_num);
    CS_FREE(recv_count);

    /* Now we do not need the unmerged halo data anymore, so can rebuild
       the new (merged) halo */

    const fvm_periodicity_t *periodicity = h->periodicity;
    cs_halo_destroy(&h);

    cs_rank_neighbors_t *rn
      = cs_rank_neighbors_create(n_c_domains, c_domain_rank);

    if (rn->size > 0) { // expected to be always true here.
      int r_idx = 0, r_id = rn->rank[0];
      for (cs_lnum_t i = 0; i < n_elts_m; i++) {
        int elt_r_id = elt_rank_idx[i];
        while (elt_r_id > r_id) {
          r_idx++;
          assert(r_idx < rn->size);
          r_id = rn->rank[r_idx];
        }
        elt_rank_idx[i] = r_idx;
      }

      if (rn->rank[0] < 0)
        rn->rank[0] = cs_glob_rank_id;
    }

    h = cs_halo_create_from_rank_neighbors(rn,
                                           n_local_elts,
                                           n_elts_m,
                                           elt_rank_idx,
                                           elt_id,
                                           elt_tr_id,
                                           periodicity);

    cs_rank_neighbors_destroy(&rn);

    CS_FREE(elt_tr_id);
    CS_FREE(elt_id);
    CS_FREE(elt_rank_idx);

  }

  else {
    CS_FREE(perio_lst_tr);

    if (g->merge_sub_size > 1)
      MPI_Recv(new_cell_id + g->n_rows, g->n_cols_ext - g->n_rows,
               CS_MPI_LNUM, g->merge_sub_root, tag, comm, &status);

    cs_halo_destroy(&h);
  }

  CS_FREE(new_src_cell_id);
  CS_FREE(c_domain_rank);

  g->halo = h;
  g->_halo = h;
}

/*----------------------------------------------------------------------------
 * Append cell arrays and counts for grid grouping and merging.
 *
 * parameters:
 *   g <-> pointer to grid structure being merged
 *----------------------------------------------------------------------------*/

static void
_append_cell_data(cs_grid_t  *g)
{
  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  const cs_lnum_t db_stride = g->db_size * g->db_size;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'c';

  /* Reallocate arrays for receiving rank and append data */

  if (g->merge_sub_rank == 0) {

    g->n_rows = g->merge_cell_idx[g->merge_sub_size];
    g->n_cols_ext = g->n_rows + g->halo->n_elts[0];

    if (g->relaxation > 0) {
      CS_REALLOC(g->_cell_cen, g->n_cols_ext, cs_real_3_t);
      CS_REALLOC(g->_cell_vol, g->n_cols_ext, cs_real_t);
    }

    CS_REALLOC(g->_da, g->n_cols_ext*db_stride, cs_real_t);

    for (int rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = (  g->merge_cell_idx[rank_id+1]
                          - g->merge_cell_idx[rank_id]);
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      if (g->relaxation > 0) {
        MPI_Recv(g->_cell_cen + g->merge_cell_idx[rank_id], n_recv*3,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->_cell_vol + g->merge_cell_idx[rank_id], n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);
      }

      MPI_Recv(g->_da + g->merge_cell_idx[rank_id]*db_stride,
               n_recv*db_stride,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

    }

  }
  else {
    if (g->relaxation > 0) {
      MPI_Send(g->_cell_cen, g->n_rows*3, CS_MPI_REAL, g->merge_sub_root,
               tag, comm);
      CS_FREE(g->_cell_cen);

      MPI_Send(g->_cell_vol, g->n_rows, CS_MPI_REAL, g->merge_sub_root,
               tag, comm);
      CS_FREE(g->_cell_vol);
    }

    MPI_Send(g->_da, g->n_rows*db_stride, CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    CS_FREE(g->_da);

    g->n_rows = 0;
    g->n_cols_ext = 0;
  }

  if (g->relaxation > 0) {
    g->cell_cen = g->_cell_cen;
    g->cell_vol = g->_cell_vol;
  }

  g->da = g->_da;
}

/*----------------------------------------------------------------------------
 * Synchronize cell array values for grid grouping and merging.
 *
 * parameters:
 *   g <-> pointer to grid structure being merged
 *----------------------------------------------------------------------------*/

static void
_sync_merged_cell_data(cs_grid_t  *g)
{
  if (g->halo != nullptr) {

    if (g->relaxation > 0) {
      cs_real_t *cell_cen = (cs_real_t *)(g->_cell_cen);
      cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD, cell_cen, 3);
      if (g->halo->n_transforms > 0)
        cs_halo_perio_sync_coords(g->halo, CS_HALO_STANDARD, cell_cen);
      cs_halo_sync(g->halo, CS_HALO_STANDARD, false, g->_cell_vol);
    }

    cs_lnum_t db_stride = g->db_size*g->db_size;
    cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD,
                             g->_da, db_stride);

  }
}

/*----------------------------------------------------------------------------
 * Append face arrays and counts for grid grouping and merging.
 *
 * parameters:
 *   g         <-> pointer to grid structure being merged
 *   n_faces   <-- number of faces to append
 *   face_list <-> list of faces to append
 *----------------------------------------------------------------------------*/

static void
_append_face_data(cs_grid_t   *g,
                  cs_lnum_t    n_faces,
                  cs_lnum_t   *face_list)
{
  cs_lnum_t *recv_count = nullptr;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'f';

  /* Exchange counters needed for concatenation */

  if (g->merge_sub_rank == 0) {
    CS_MALLOC(recv_count, g->merge_sub_size, cs_lnum_t);
    recv_count[0] = g->n_faces;
    for (int rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      MPI_Recv(recv_count + rank_id, 1, CS_MPI_LNUM, dist_rank,
               tag, comm, &status);
    }
  }
  else
    MPI_Send(&n_faces, 1, CS_MPI_LNUM, g->merge_sub_root, tag, comm);

  /* Reallocate arrays for receiving rank and append data */

  CS_FREE(g->coarse_face);

  if (g->merge_sub_rank == 0) {

    cs_lnum_t n_faces_tot = 0;

    for (int rank_id = 0; rank_id < g->merge_sub_size; rank_id++)
      n_faces_tot += recv_count[rank_id];

    CS_REALLOC(g->_face_cell, n_faces_tot, cs_lnum_2_t);

    if (g->symmetric == true)
      CS_REALLOC_HD(g->_xa, n_faces_tot, cs_real_t, g->alloc_mode);
    else
      CS_REALLOC_HD(g->_xa, n_faces_tot*2, cs_real_t, g->alloc_mode);

    if (g->relaxation > 0) {

      if (g->face_normal != nullptr)
        CS_REALLOC(g->face_normal, n_faces_tot, cs_real_3_t);

      CS_REALLOC(g->_xa0, n_faces_tot, cs_real_t);
      CS_REALLOC(g->xa0ij, n_faces_tot*3, cs_real_t);

    }

    g->n_faces = recv_count[0];

    for (int rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = recv_count[rank_id];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_face_cell + g->n_faces, n_recv*2,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);

      if (g->symmetric == true)
        MPI_Recv(g->_xa + g->n_faces, n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);
      else
        MPI_Recv(g->_xa + g->n_faces*2, n_recv*2,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

      if (g->relaxation > 0) {

        if (g->face_normal != nullptr)
          MPI_Recv(g->face_normal + g->n_faces, n_recv*3,
                   CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->_xa0 + g->n_faces, n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->xa0ij + g->n_faces*3, n_recv*3,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

      }

      g->n_faces += recv_count[rank_id];
    }

    CS_FREE(recv_count);

  }
  else {

    /* Filter face connectivity then send it */

    for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t p_face_id = face_list[face_id];
      assert(face_id <= p_face_id);
      g->_face_cell[face_id][0] = g->_face_cell[p_face_id][0];
      g->_face_cell[face_id][1] = g->_face_cell[p_face_id][1];
      if (g->symmetric == true)
        g->_xa[face_id] = g->_xa[p_face_id];
      else {
        g->_xa[face_id*2] = g->_xa[p_face_id*2];
        g->_xa[face_id*2+1] = g->_xa[p_face_id*2+1];
      }
      if (g->relaxation > 0) {
        if (g->face_normal != nullptr) {
          for (cs_lnum_t k = 0; k < 3; k++)
            g->face_normal[face_id][k] = g->face_normal[p_face_id][k];
        }
        g->_xa0[face_id] = g->_xa0[p_face_id];
        g->xa0ij[face_id*3] = g->xa0ij[p_face_id*3];
        g->xa0ij[face_id*3 + 1] = g->xa0ij[p_face_id*3 + 1];
        g->xa0ij[face_id*3 + 2] = g->xa0ij[p_face_id*3 + 2];
      }
    }

    MPI_Send(g->_face_cell, n_faces*2, CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    CS_FREE(g->_face_cell);

    if (g->symmetric == true)
      MPI_Send(g->_xa, n_faces, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    else
      MPI_Send(g->_xa, n_faces*2, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    CS_FREE(g->_xa);

    if (g->relaxation > 0) {
      if (g->face_normal != nullptr) {
        MPI_Send(g->face_normal, n_faces*3, CS_MPI_REAL,
                 g->merge_sub_root, tag, comm);
        CS_FREE(g->face_normal);
      }

      MPI_Send(g->_xa0, n_faces, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
      CS_FREE(g->_xa0);

      MPI_Send(g->xa0ij, n_faces*3, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
      CS_FREE(g->xa0ij);
    }

    g->n_faces = 0;
  }

  g->face_cell = g->_face_cell;
  g->xa = g->_xa;
  if (g->relaxation > 0) {
    g->xa0 = g->_xa0;
  }
}

/*----------------------------------------------------------------------------
 * Compute stride and ranks when merging grids on several ranks to a smaller
 * number of ranks.
 *
 * parameters:
 *   g            <-- Pointer to grid structure
 *   merge_stride <-- Associated merge stride
 *----------------------------------------------------------------------------*/

static void
_compute_merge_stride_and_ranks(cs_grid_t  *g,
                                int         merge_stride)
{
  int base_rank = cs_glob_rank_id;

  g->merge_sub_size = merge_stride;
  g->merge_stride = g->next_merge_stride;
  g->next_merge_stride *= merge_stride;

  if (base_rank % g->merge_stride != 0) {
    g->merge_sub_size = 0;
    g->merge_sub_root = -1;
    g->merge_sub_rank = -1;
  }
  else {
    int merge_rank = base_rank / g->merge_stride;
    int merge_size = (cs_glob_n_ranks/g->merge_stride);
    int merge_sub_root = (merge_rank / merge_stride) * merge_stride;
    if (cs_glob_n_ranks % g->merge_stride > 0)
      merge_size += 1;
    g->merge_sub_root = merge_sub_root * g->merge_stride;
    g->merge_sub_rank = merge_rank % merge_stride;
    if (merge_sub_root + g->merge_sub_size > merge_size)
      g->merge_sub_size = merge_size - merge_sub_root;
  }
}

#endif // defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Build coarse grid matrix in the MSR format given the 4 necessary arrays.
 *
 * parameters:
 *   c           <-> coarse grid structure
 *   symmetric   <-- symmetry flag
 *   c_row_index <-- MSR row index (0 to n-1)
 *   c_col_id    <-- MSR column id (0 to n-1)
 *   c_d_val     <-- diagonal values
 *   c_x_val     <-- extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_build_coarse_matrix_msr(cs_grid_t *c,
                         bool       symmetric,
                         cs_lnum_t *c_row_index,
                         cs_lnum_t *c_col_id,
                         cs_real_t *c_d_val,
                         cs_real_t *c_x_val)
{
  cs_matrix_structure_t  *ms
    = cs_matrix_structure_create_msr(CS_MATRIX_MSR,
                                     true, /* transfer */
                                     true, /* have_diag */
                                     true, /* ordered */
                                     c->n_rows,
                                     c->n_cols_ext,
                                     &c_row_index,
                                     &c_col_id,
                                     c->halo,
                                     nullptr);

  c->matrix_struct = ms;

  c->_matrix = cs_matrix_create(c->matrix_struct);
  c->matrix = c->_matrix;

  const cs_lnum_t *_c_row_index, *_c_col_id;

  cs_matrix_get_msr_arrays(c->matrix,
                           &_c_row_index, &_c_col_id,
                           nullptr, nullptr);

  cs_matrix_transfer_coefficients_msr(c->_matrix,
                                      symmetric,
                                      c->db_size,
                                      c->eb_size,
                                      _c_row_index,
                                      _c_col_id,
                                      &c_d_val,
                                      &c_x_val);
}

/*----------------------------------------------------------------------------
 * Create a distributed matrix structure from a native matrix structure.
 *
 * Note that the structure created maps global cell numbers to the given
 * existing face -> cell connectivity array, so it must be destroyed before
 * this array (usually the code's global cell numbering) is freed.
 *
 * parameters:
 *   n     <-- number of elements
 *   index <-> count (shifted by 1) in, index out (size: n+1)
 *----------------------------------------------------------------------------*/

static void
_count_to_index(cs_lnum_t  n,
                cs_lnum_t  index[])
{
  for (cs_lnum_t i = 0; i < n; i++)
    index[i+1] += index[i];
}

/*----------------------------------------------------------------------------
 * Prune MSR arrays to remove duplicate columns.
 *
 * This should occur only after discarding periodicity info, which only makes
 * sense for scalar matrices.
 *
 * parameters:
 *   matrix     <-- pointer to grid structure
 *   o2n        <-- optional old to new column ids, or null
 *   row_index  --> matrix row index
 *   col_id     --> matrix column ids
 *   d_val      --> matrix diagonal values
 *   x_val      --> matrix extradiagonal values
 *----------------------------------------------------------------------------*/

static void
_matrix_pruned_msr_arrays(cs_grid_t *          g,
                          const cs_lnum_t *    o2n,
                          cs_lnum_t *         &row_index,
                          cs_lnum_t *         &col_id,
                          cs_real_t *         &d_val,
                          cs_real_t *         &x_val)
{
  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  cs_matrix_t *matrix = g->_matrix;
  cs_alloc_mode_t alloc_mode = cs_matrix_get_alloc_mode(matrix);

  cs_assert(cs_matrix_get_diag_block_size(matrix) == 1);

  const cs_lnum_t n_rows = cs_matrix_get_n_rows(matrix);

  cs_matrix_release_msr_arrays(matrix,
                               &row_index,
                               &col_id,
                               &d_val,
                               &x_val);

  if (d_val == nullptr)
    d_val = g->_da;
  if (x_val == nullptr)
    x_val = g->_xa;

  if (row_index == nullptr && col_id == nullptr) {
    cs_matrix_structure_release_msr_arrays(g->matrix_struct,
                                           &row_index,
                                           &col_id);
  }

  g->da = nullptr;
  g->xa = nullptr;
  if (d_val != g->_da)
    CS_FREE(g->_da);
  else
    g->_da = nullptr;
  g->_da = nullptr;
  if (x_val != g->_xa)
    CS_FREE(g->_xa);
  else
    g->_xa = nullptr;

  cs_lnum_t *c_row_index;
  CS_MALLOC_HD(c_row_index, n_rows+1, cs_lnum_t, alloc_mode);
  c_row_index[0] = 0;

  int n_threads = cs_parall_n_threads(n_rows, CS_THR_MIN);
  bool need_compact = false;

  if (o2n != nullptr) {
    cs_lnum_t nnz = row_index[n_rows];

#   pragma omp parallel for  num_threads(n_threads)
    for (cs_lnum_t i = 0; i < nnz; i++) {
      col_id[i] = o2n[col_id[i]];
    }
  }

# pragma omp parallel for  num_threads(n_threads)
  for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
    cs_lnum_t *restrict row_col_id = col_id + row_index[row_id];
    cs_real_t *restrict row_val = x_val + row_index[row_id];
    cs_lnum_t n_row_elts = row_index[row_id+1] - row_index[row_id];
    bool unordered = false;
    for (cs_lnum_t i = 1; i < n_row_elts; i++) {
      if (row_col_id[i] <= row_col_id[i-1])
        unordered = true;
    }
    c_row_index[row_id + 1] = n_row_elts;

    if (unordered) {
      cs_sort_coupled_shell(0, n_row_elts, row_col_id, row_val);
      cs_lnum_t j = 0;
      for (cs_lnum_t i = 1; i < n_row_elts; i++) {
        if (row_col_id[i] > row_col_id[j]) {
          j += 1;
          row_col_id[j] = row_col_id[i];
          row_val[j] = row_val[i];
        }
        else
          row_val[j] += row_val[i];
      }
      c_row_index[row_id + 1] = j+1;
      if (j+1 < n_row_elts)
        need_compact = true;  /* no issue here for possible thread race, as
                                 this can only take one value */
    }
  }

  if (need_compact) {
    _count_to_index(n_rows, c_row_index);
    cs_lnum_t nnz = c_row_index[n_rows];

    cs_lnum_t *c_col_id;
    cs_real_t *c_x_val;
    CS_MALLOC_HD(c_col_id, nnz, cs_lnum_t, alloc_mode);
    CS_MALLOC_HD(c_x_val, nnz, cs_real_t, alloc_mode);

#   pragma omp parallel for  num_threads(n_threads)
    for (cs_lnum_t row_id = 0; row_id < n_rows; row_id++) {
      const cs_lnum_t src_shift = row_index[row_id];
      const cs_lnum_t dst_shift = c_row_index[row_id];
      const cs_lnum_t n_row_elts = c_row_index[row_id+1] - c_row_index[row_id];
      for (cs_lnum_t i = 0; i < n_row_elts; i++) {
        c_col_id[dst_shift + i] = col_id[src_shift + i];
        c_x_val[dst_shift + i] = x_val[src_shift + i];
      }
    }

    CS_FREE(row_index);
    CS_FREE(col_id);
    CS_FREE(x_val);

    row_index = c_row_index;
    col_id = c_col_id;
    x_val = c_x_val;
  }
  else
    CS_FREE(c_row_index);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      tm_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(tm_stop - tm_start);
    printf("%d: %s", cs_glob_rank_id, __func__);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Merge bottom grids on several ranks to a grid on a single rank
 *
 * parameters:
 *   g            <-- Pointer to grid structure
 *   merge_stride <-- Associated merge stride
 *   verbosity    <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_merge_bottom_grids_to_single(cs_grid_t  *g,
                              int         verbosity)
{
  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  if (g->db_size > 1 || g->level == 0)
    return;

  if (g->halo == nullptr)
    return;

  cs_alloc_mode_t alloc_mode = cs_matrix_get_alloc_mode(g->matrix);

  const cs_lnum_t n_rows = g->n_rows;
  cs_lnum_t base_shift = 0;

  int base_rank = 0, n_ranks = cs_glob_n_ranks;

#if defined(HAVE_MPI)
  if (g->comm != MPI_COMM_NULL) {
    MPI_Comm_rank(g->comm, &base_rank);
    MPI_Comm_size(g->comm, &n_ranks);
    MPI_Exscan(&n_rows, &base_shift, 1, CS_MPI_LNUM, MPI_SUM, g->comm);
  }
#endif

  /* Compute and exchange new row_ids */

  cs_lnum_t *new_row_id;
  CS_MALLOC(new_row_id, g->n_cols_ext, cs_lnum_t);
  for (cs_lnum_t j = 0; j < n_rows; j++)
    new_row_id[j] = base_shift + j;

  cs_halo_sync_untyped(g->halo,
                       CS_HALO_STANDARD,
                       sizeof(cs_lnum_t),
                       new_row_id);

  /* Access matrix MSR vectors */

  cs_lnum_t  nnz = 0;
  cs_lnum_t  *row_index = nullptr, *col_id = nullptr;
  cs_real_t  *d_val = nullptr, *x_val = nullptr;

  if (g->_matrix != nullptr) {
    _matrix_pruned_msr_arrays(g,
                              new_row_id,
                              row_index,
                              col_id,
                              d_val,
                              x_val);
    nnz = row_index[n_rows];
  }

  CS_FREE(new_row_id);

  cs_matrix_destroy(&(g->_matrix));
  g->matrix = nullptr;
  cs_matrix_structure_destroy(&(g->matrix_struct));
  if (g->_halo != nullptr)
    cs_halo_destroy(&(g->_halo));
  g->halo = nullptr;

  /* Discard geometric data */

  cs_grid_free_quantities(g);
  CS_FREE(g->_face_cell);
  CS_FREE(g->cell_face);
  CS_FREE(g->cell_face_sgn);

  /* Determine rank in merged group */

#if defined(HAVE_MPI)
  _compute_merge_stride_and_ranks(g, g->n_ranks);

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid:\n"
               "      merge_stride = %d; new ranks = %d\n",
               g->level, g->merge_stride, g->n_ranks);
    if (verbosity > 3)
      bft_printf("      sub_root = %d; sub_rank = %d; sub_size = %d\n",
                 g->merge_sub_root, g->merge_sub_rank, g->merge_sub_size);
  }

#else

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid periodic columns\n",
               g->level);
  }

#endif

  /* Compute number of coarse cells in new grid;
     local cell numbers will be shifted by cell_shift in the merged grid */

  cs_lnum_t *rank_index = nullptr;

  cs_lnum_t counts[2] = {g->n_rows, nnz};

#if defined(HAVE_MPI)

  if (base_rank == 0) {
    CS_MALLOC(g->merge_cell_idx, g->merge_sub_size + 1, cs_lnum_t);
    CS_MALLOC(rank_index, (n_ranks+1)*2, cs_lnum_t);
    rank_index[0] = 0; rank_index[1] = 0;
  }

  if (g->comm != MPI_COMM_NULL) {
    MPI_Gather(counts, 2, CS_MPI_LNUM, rank_index+2, 2, CS_MPI_LNUM,
               0, g->comm);
  }

  if (base_rank == 0) {
    rank_index[2] = g->n_rows; rank_index[3] = nnz;
    g->merge_cell_idx[0] = 0; g->merge_cell_idx[1] = g->n_rows;
    for (int i = 1; i < n_ranks; i++) {
      rank_index[(i+1)*2] += rank_index[i*2];
      rank_index[(i+1)*2+1] += rank_index[i*2 + 1];
      g->merge_cell_idx[i + 1] = rank_index[(i+1)*2];
    }
  }

  /* Now return computed info to grids that will be merged */

  if (g->comm != MPI_COMM_NULL) {
    MPI_Scatter(rank_index, 2, CS_MPI_LNUM, counts, 2, CS_MPI_LNUM,
                0, g->comm);
  }

#else

  {
    CS_MALLOC(rank_index, (n_ranks+1)*2, cs_lnum_t);
    rank_index[0] = 0; rank_index[1] = 0;
    rank_index[2] = g->n_rows; rank_index[3] = nnz;
  }

#endif

  /* Now update row index based on local shifts */

  if (base_rank > 0) {
    for (cs_lnum_t j = 0; j < n_rows+1; j++)
      row_index[j] += counts[1];
  }

  /* Now send matrix data to root rank */

  int *recvcounts = nullptr, *displs = nullptr;
  int n_send = 0;

  if (base_rank == 0) {
    g->n_rows = rank_index[n_ranks*2];
    g->n_cols_ext = rank_index[n_ranks*2];

    assert(g->n_g_rows == (cs_gnum_t)(g->n_rows));

    CS_REALLOC_HD(row_index, rank_index[n_ranks*2] + 1, cs_lnum_t, alloc_mode);
    CS_REALLOC_HD(col_id, rank_index[n_ranks*2 + 1], cs_lnum_t, alloc_mode);
    CS_REALLOC_HD(d_val, rank_index[n_ranks*2], cs_real_t, alloc_mode);
    CS_REALLOC_HD(x_val, rank_index[n_ranks*2 + 1], cs_real_t, alloc_mode);

    CS_MALLOC(recvcounts, n_ranks, int);
    CS_MALLOC(displs, n_ranks, int);

    recvcounts[0] = 0;
    displs[0] = 0;
    for (int i = 1; i < n_ranks; i++) {
      recvcounts[i] = rank_index[2*(i+1)] - rank_index[2*i];
      displs[i] = rank_index[2*i];
    }
  }
  else {
    n_send = g->n_rows;
  }

#if defined(HAVE_MPI)
  if (g->comm != MPI_COMM_NULL) {
    MPI_Gatherv(row_index+1, n_send, CS_MPI_LNUM,
                row_index+1, recvcounts, displs, CS_MPI_LNUM,
                0, g->comm);

    MPI_Gatherv(d_val, n_send, CS_MPI_REAL,
                d_val, recvcounts, displs, CS_MPI_REAL,
                0, g->comm);
  }
#endif

  if (base_rank == 0) {
    for (int i = 1; i < n_ranks; i++) {
      recvcounts[i] = rank_index[2*(i+1) + 1] - rank_index[2*i + 1];
      displs[i] = rank_index[2*i + 1];
    }
    CS_FREE(rank_index);
  }
  else {
    n_send = nnz;
  }

#if defined(HAVE_MPI)
  if (g->comm != MPI_COMM_NULL) {
    MPI_Gatherv(col_id, n_send, CS_MPI_LNUM,
                col_id, recvcounts, displs, CS_MPI_LNUM,
                0, g->comm);

    MPI_Gatherv(x_val, n_send, CS_MPI_REAL,
                x_val, recvcounts, displs, CS_MPI_REAL,
                0, g->comm);
  }
#endif

  CS_FREE(recvcounts);
  CS_FREE(displs);

  if (base_rank > 0) {
    g->n_rows = 0;
    g->n_cols_ext = 0;

    CS_FREE(row_index);
    CS_MALLOC_HD(row_index, 1, cs_lnum_t, g->alloc_mode);
    row_index[0] = 0;
    CS_FREE(col_id);
    CS_FREE(d_val);
    CS_FREE(x_val);
  }

  _build_coarse_matrix_msr(g,
                           g->symmetric,
                           row_index,
                           col_id,
                           d_val,
                           x_val);

  /* The matrix is now on a single rank */

#if defined(HAVE_MPI)

  g->n_ranks = 1;
  g->comm = MPI_COMM_NULL;

#endif

  if (verbosity > 3)
    bft_printf("      merged to %ld (from %ld) rows\n\n",
               (long)g->n_rows, (long)g->n_elts_r[0]);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      tm_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(tm_stop - tm_start);
    printf("%d: %s (level %d)", cs_glob_rank_id, __func__, g->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*============================================================================
 * Semi-private function definitions
 *
 * The following functions are intended to be used by the multigrid layer
 * (cs_multigrid.c), not directly by the user, so they are no more
 * documented than private static functions)
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Merge grids on several ranks to a grid on a single rank
 *
 * The coarse grid must previously have been initialized with _coarse_init()
 * and its coarsening indicator determined (at least for the local cells;
 * it is extended to halo cells here if necessary).
 *
 * parameters:
 *   g            <-- Pointer to grid structure
 *   merge_stride <-- Associated merge stride
 *   verbosity    <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_grid_merge_ranks(cs_grid_t  *g,
                    int         merge_stride,
                    int         verbosity)
{
  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  cs_lnum_t cell_shift = 0;
  cs_lnum_t n_faces = 0;
  cs_lnum_t *new_cell_id = nullptr, *face_list = nullptr;
  bool  *halo_cell_flag = nullptr;
  MPI_Comm comm = cs_glob_mpi_comm;
  MPI_Status status;

  static const int tag = 'm'+'e'+'r'+'g'+'e';

  if (merge_stride < 2)
    return;

  /* Determine rank in merged group */

  _compute_merge_stride_and_ranks(g, merge_stride);

  if (g->next_merge_stride > 1) {
    if (g->n_ranks % merge_stride)
      g->n_ranks = (g->n_ranks/merge_stride) + 1;
    else
      g->n_ranks = (g->n_ranks/merge_stride);
    g->comm = cs_base_get_rank_step_comm_recursive(g->comm, merge_stride);
  }

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid:\n"
               "      merge_stride = %d; new ranks = %d\n",
               g->level, g->merge_stride, g->n_ranks);
    if (verbosity > 3)
      bft_printf("      sub_root = %d; sub_rank = %d; sub_size = %d\n",
                 g->merge_sub_root, g->merge_sub_rank, g->merge_sub_size);
  }

  /* Compute number of coarse cells in new grid;
     local cell numbers will be shifted by cell_shift in the merged grid */

  if (g->merge_sub_size > 1) {

    if (g->merge_sub_rank == 0) {
      CS_MALLOC(g->merge_cell_idx, g->merge_sub_size + 1, cs_lnum_t);
      g->merge_cell_idx[0] = 0; g->merge_cell_idx[1] = g->n_rows;
      for (int i = 1; i < g->merge_sub_size; i++) {
        cs_lnum_t recv_val;
        int dist_rank = g->merge_sub_root + g->merge_stride*i;
        MPI_Recv(&recv_val, 1, CS_MPI_LNUM, dist_rank, tag, comm, &status);
        g->merge_cell_idx[i + 1] = recv_val + g->merge_cell_idx[i];
      }
    }
    else {
      cs_lnum_t send_val = g->n_rows;
      MPI_Send(&send_val, 1, CS_MPI_LNUM, g->merge_sub_root, tag, comm);
    }

    /* Now return computed info to grids that will be merged */

    if (g->merge_sub_rank == 0) {
      for (int i = 1; i < g->merge_sub_size; i++) {
        cs_lnum_t send_val;
        send_val = g->merge_cell_idx[i];
        MPI_Send(&send_val, 1, CS_MPI_LNUM,
                 g->merge_sub_root + g->merge_stride*i, tag, comm);
      }
    }
    else
      MPI_Recv(&cell_shift, 1, CS_MPI_LNUM,
               g->merge_sub_root, tag, comm, &status);
  }

  /* Compute and exchange new cell numbers */

  CS_MALLOC(new_cell_id, g->n_cols_ext, cs_lnum_t);
  for (cs_lnum_t j = 0; j < g->n_rows; j++)
    new_cell_id[j] = cell_shift + j;
  for (cs_lnum_t j = g->n_rows; j < g->n_cols_ext; j++)
    new_cell_id[j] = -1;

  cs_halo_sync_untyped(g->halo,
                       CS_HALO_STANDARD,
                       sizeof(cs_lnum_t),
                       new_cell_id);

  /* Now build face filter list (before halo is modified) */

  if (g->merge_sub_size > 1 && g->merge_sub_rank > 0 && g->n_faces > 0) {

    cs_lnum_t n_ghost_cells = g->n_cols_ext - g->n_rows;

    /* Mark which faces should be merged: to avoid duplicates, a face
       connected to a cell on a lower rank in the same merge set is
       discarded, as it has already been accounted for by that rank. */

    CS_MALLOC(face_list, g->n_faces, cs_lnum_t);
    CS_MALLOC(halo_cell_flag, n_ghost_cells, bool);
    for (cs_lnum_t j = 0; j < n_ghost_cells; j++)
      halo_cell_flag[j] = false;

    for (cs_lnum_t i = 0; i < g->halo->n_c_domains; i++) {
      int rank_id = g->halo->c_domain_rank[i];
      if (rank_id >= g->merge_sub_root && rank_id < cs_glob_rank_id) {
        for (int t_id = 0; t_id < g->halo->n_transforms; t_id++) {
          int t_shift = 4 * g->halo->n_c_domains * t_id;
          cs_lnum_t t_start = g->halo->perio_lst[t_shift + 4*i];
          cs_lnum_t t_end = t_start + g->halo->perio_lst[t_shift + 4*i + 1];
          for (cs_lnum_t j = t_start; j < t_end; j++)
            halo_cell_flag[j] = true;
        }
      }
      else {
        for (cs_lnum_t j = g->halo->index[2*i]; j < g->halo->index[2*i+2]; j++)
          halo_cell_flag[j] = true;
      }
    }

    for (cs_lnum_t face_id = 0; face_id < g->n_faces; face_id++) {
      bool use_face = true;
      cs_lnum_t ii = g->face_cell[face_id][0] - g->n_rows;
      cs_lnum_t jj = g->face_cell[face_id][1] - g->n_rows;
      if (ii > -1) {
        if (halo_cell_flag[ii] == false)
          use_face = false;
      }
      else if (jj > -1) {
        if (halo_cell_flag[jj] == false)
          use_face = false;
      }
      if (use_face == true)
        face_list[n_faces++] = face_id;
    }

    CS_FREE(halo_cell_flag);
  }

  /* Append and merge halos */

  _append_halos(g, new_cell_id);

  /* Update face ->cells connectivity */

  for (cs_lnum_t face_id = 0; face_id < g->n_faces; face_id++) {
    cs_lnum_t ii = g->face_cell[face_id][0];
    cs_lnum_t jj = g->face_cell[face_id][1];
    assert(ii != jj && new_cell_id[ii] != new_cell_id[jj]);
    g->_face_cell[face_id][0] = new_cell_id[ii];
    g->_face_cell[face_id][1] = new_cell_id[jj];
  }

  CS_FREE(new_cell_id);

  /* Merge cell and face data */

  if (g->merge_sub_size > 1) {
    _append_cell_data(g);
    _append_face_data(g, n_faces, face_list);
    CS_FREE(g->cell_face);
    CS_FREE(g->cell_face_sgn);
  }
  _sync_merged_cell_data(g);

  CS_FREE(face_list);

  if (verbosity > 3)
    bft_printf("      merged to %ld (from %ld) rows\n\n",
               (long)g->n_rows, (long)g->n_elts_r[0]);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      tm_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(tm_stop - tm_start);
    printf("%d: %s (level %d)", cs_glob_rank_id, __func__, g->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Merge bottom grid over adjacent ranks, to reduce communicator size.
 *
 * If this would lead to a mean number of rows exceeding a given threshold,
 * no merging is done.
 *
 * parameters:
 *   g                <-- Grid structure
 *   verbosity        <-- Verbosity level
 *   n_max_ranks      <-- Maximum number of MPI ranks for bottom grid
 *   max_row_factor   <-- maximum acceptable mean ratio of merged rows
 *                        (per MPI rank) to finest rows.
 *----------------------------------------------------------------------------*/

void
cs_grid_merge_bottom(cs_grid_t      *g,
                     int             verbosity,
                     int             n_max_ranks,
                     float           max_row_factor)
{
#if defined(HAVE_MPI)

  n_max_ranks = cs::min(g->n_ranks, n_max_ranks);
  float n_mean_rows_at_root;
  {
    const cs_grid_t * rg = _root_grid(g);
    n_mean_rows_at_root = (float)(rg->n_g_rows) / (float)(rg->n_ranks);
  }

  float n_mean_rows = (float)(g->n_g_rows) / (float)(n_max_ranks);
  if (n_mean_rows / n_mean_rows_at_root > max_row_factor)
    return;

  if (n_max_ranks == 1) {
    _merge_bottom_grids_to_single(g, verbosity);
  }

  else {
    bft_error(__FILE__, __LINE__, 0,
              "%s: not yet implemented for merge to multiple ranks.\n",
              __func__);
  }

#else

  /* Without MPI, merging can still be useful to replace periodic halo values
     with regular values, simplifying matrix). */

  if (g->halo != nullptr)
    _merge_bottom_grids_to_single(g, verbosity);

#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
