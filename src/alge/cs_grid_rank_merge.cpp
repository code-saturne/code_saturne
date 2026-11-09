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
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cmath>

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
#include "base/cs_profiling.h"
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
  \file cs_grid_rank_lerge.cpp
        Coarse grid merging over multiple ranks for multigrid solver.
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
 * Return a reduced communicator restricted to the first node of a given
 * parent communicator.
 *
 * The lowest ranks of the reduced communicator are used, so as to map
 * to a smaller number of nodes, assuming ranks on adjacent nodes have
 * neighboring ids (instead of actually checking the topology).
 *
 * parameters:
 *   parent_comm     <-- associated parent communicator (must be either
 *                       cs_glob_mpi_comm or a communicator returned by a
 *                       previous merge step)
 *   n_max_ranks     --> maximum number of ranks for communicator
 *   rank_step       --> associated multiple of ranks of parent communicator
 *   merge_rank_step --> rank to which grid on parent communicator is merged
 *----------------------------------------------------------------------------*/

static MPI_Comm
_get_rank_step_comm_compact(MPI_Comm   parent_comm,
                            int        n_max_ranks,
                            int       &rank_step,
                            int       &merge_rank_id)
{
  MPI_Comm sh_comm = MPI_COMM_NULL;

  rank_step = 1;
  merge_rank_id = -1;

#if (MPI_VERSION >= 3)
  /* Check which MPI ranks are on root node */

  int parent_n_ranks;
  MPI_Comm_size(parent_comm, &parent_n_ranks);

  int parent_rank_id;
  MPI_Comm_rank(parent_comm, &parent_rank_id);

  MPI_Comm_split_type(cs_glob_mpi_comm, MPI_COMM_TYPE_SHARED, 0,
                      MPI_INFO_NULL, &sh_comm);

  int node_rank_id = -1, node_n_ranks = 0;
  MPI_Comm_rank(sh_comm, &node_rank_id);
  MPI_Comm_size(sh_comm, &node_n_ranks);

  int parent_min_rank_id = parent_rank_id;
  MPI_Allreduce(MPI_IN_PLACE, &parent_min_rank_id, 1, MPI_INT, MPI_MIN,
                sh_comm);

  if (parent_min_rank_id > 0) {
    MPI_Comm_free(&sh_comm);
    sh_comm = MPI_COMM_NULL;
  }
  else {
    node_rank_id = -1;
  }

  MPI_Bcast(&node_n_ranks, 1, MPI_INT, 0, parent_comm);

  int *node_rank_ids;
  CS_MALLOC(node_rank_ids, node_n_ranks, int);
  if (sh_comm != MPI_COMM_NULL)
    MPI_Gather(&parent_rank_id, 1, MPI_INT,
               node_rank_ids, 1, MPI_INT, 0, sh_comm);

  MPI_Bcast(&node_rank_ids, node_n_ranks, MPI_INT, 0, parent_comm);

  int parent_rank_stride = 1;
  if (n_max_ranks < node_n_ranks)
    parent_rank_stride = node_n_ranks / n_max_ranks;

  int comm_n_ranks = cs::min(n_max_ranks, node_n_ranks);
  rank_step = parent_n_ranks / comm_n_ranks;
  if (parent_n_ranks % comm_n_ranks)
    rank_step += 1;

  merge_rank_id = node_rank_ids[(parent_rank_id/rank_step)*parent_rank_stride];

  CS_FREE(node_rank_ids);

#else
  cs_assert(0);  // Requires at least MPI-3
#endif

  return sh_comm;
}

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
 *   g                  <-> pointer to grid structure being merged
 *   new_cell_id        <-> new cell ids for local cells
 *                          in: defined for local cells
 *                          out: updated also for halo cells
 *----------------------------------------------------------------------------*/

static void
_append_halos(cs_grid_t  *g,
              cs_lnum_t  *new_cell_id)
{
  int *recv_count = nullptr;
  cs_lnum_t *new_src_cell_id = nullptr;

  cs_halo_t *h = g->_halo;

  if (h == nullptr)
    return;

  auto *grm = g->rank_merge;

  MPI_Comm  comm = cs_glob_mpi_comm;

  const int n_transforms = h->n_transforms;
  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'h';

  if (grm->sub_size[0] == 0) {
    _empty_halo(h);
    return;
  }

  /* local copy of halo sizes */

  int n_c_domains = 0;
  cs_lnum_t n_local_elts = 0;
  cs_lnum_t n_elts = 0;

  int n_requests = 1;
  bool is_receiver = (grm->rank_id == cs_glob_rank_id) ? true : false;
  if (is_receiver)
    n_requests += grm->sub_size[0];

  int  request_count = 0;

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  /* Allocate if necessary */

  if (n_requests*4 > 128) {
    CS_MALLOC(request, n_requests*4, MPI_Request);
    CS_MALLOC(status, n_requests*4, MPI_Status);
  }

  /* Exchange counters needed for concatenation */

  int recv_counts[3] = {0, 0, 0};
  int send_counts[3] = {h->n_c_domains, h->n_local_elts, h->n_elts[0]};

  if (is_receiver) {
    CS_MALLOC(recv_count, grm->sub_size[0]*3, int);
    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;
      MPI_Irecv(recv_count + 3*rank_idx, 3, MPI_INT, dist_rank,
                tag, comm, &request[request_count++]);
    }
  }
  MPI_Isend(send_counts, 3, MPI_INT, grm->rank_id,
            tag, comm, &request[request_count++]);

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  if (is_receiver) {
    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      for (int ii = 0; ii < 3; ii++)
        recv_counts[ii] += recv_count[rank_idx*3 + ii];
    }
  }

  /* In case of periodic transforms, transpose perio list
     so as to have blocs of fixed n_transforms size, easier
     to work with for append + merge operations */

  int *perio_lst_tr_s = nullptr, *perio_lst_tr_r = nullptr;
  if (n_transforms > 0) {
    CS_MALLOC(perio_lst_tr_s, send_counts[0]*n_transforms*2, cs_lnum_t);
    CS_MALLOC(perio_lst_tr_r, recv_counts[0]*n_transforms*2, cs_lnum_t);
    for (int rank_idx = 0; rank_idx < h->n_c_domains; rank_idx++) {
      for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
        for (int k = 0; k < 2; k++)
          perio_lst_tr_s[(n_transforms*rank_idx + tr_id)*2 + k]
            = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx + k];
      }
    }
    CS_FREE(h->perio_lst);
  }

  /* Adjust numbering */

  int *c_domain_rank_recv = nullptr, *c_domain_rank_send = nullptr;
  cs_lnum_t *index = nullptr;

  if (is_receiver) {

    CS_MALLOC(c_domain_rank_recv, recv_counts[0], int);
    CS_MALLOC(new_src_cell_id, recv_counts[2], cs_lnum_t);

    /* Allocate index for receiving append data */
    CS_MALLOC(index, recv_counts[0]*2 + 1, cs_lnum_t);
    index[0] = 0;

    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      int n_c_domains_r = recv_count[rank_idx*3 + 0];
      cs_lnum_t n_recv = recv_count[rank_idx*3 + 2];

      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;

      MPI_Irecv(c_domain_rank_recv + n_c_domains, n_c_domains_r,
                MPI_INT, dist_rank, tag, comm, &request[request_count++]);
      MPI_Irecv(new_src_cell_id + n_elts, n_recv,
                CS_MPI_LNUM, dist_rank, tag, comm, &request[request_count++]);
      MPI_Irecv(index + n_c_domains*2 + 1, n_c_domains_r*2,
                CS_MPI_LNUM, dist_rank, tag, comm, &request[request_count++]);

      if (n_transforms > 0) {
        MPI_Irecv(perio_lst_tr_r + n_c_domains*n_transforms*2,
                  n_c_domains_r*n_transforms*2,
                  CS_MPI_LNUM, dist_rank, tag, comm, &request[request_count++]);
      }

      /* Update halo sizes */

      n_local_elts += recv_count[rank_idx*3 + 1];
      n_c_domains += n_c_domains_r;
      n_elts += n_recv;
    }

  }

  {
    CS_MALLOC(c_domain_rank_send, h->n_c_domains, int);

    for (int rank_idx = 0; rank_idx < h->n_c_domains; rank_idx++) {
      int c_rank_idx = h->c_domain_rank[rank_idx] / g->rank_stride;
      c_domain_rank_send[rank_idx] = c_rank_idx  * g->rank_stride;
    }

    assert(h->index[0] == 0);
    MPI_Isend(c_domain_rank_send, h->n_c_domains, MPI_INT,
              grm->rank_id, tag, comm, &request[request_count++]);
    MPI_Isend(new_cell_id + g->n_rows, h->n_elts[0], CS_MPI_LNUM,
              grm->rank_id, tag, comm, &request[request_count++]);
    MPI_Isend(h->index+1, h->n_c_domains*2, CS_MPI_LNUM,
              grm->rank_id, tag, comm, &request[request_count++]);

    if (n_transforms > 0)
      MPI_Isend(perio_lst_tr_s, h->n_c_domains*n_transforms*2,
                CS_MPI_LNUM, grm->rank_id, tag, comm,
                &request[request_count++]);
  }

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  CS_FREE(c_domain_rank_send);

  cs_lnum_t *halo_o2n = nullptr;

  if (is_receiver) {

    n_c_domains = 0;
    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      int n_c_domains_r = recv_count[rank_idx*3 + 0];
      cs_lnum_t index_shift = index[2*n_c_domains];

      for (cs_lnum_t ii = 0, jj = n_c_domains*2+1;
           ii < n_c_domains_r*2;
           ii++, jj++)
        index[jj] += index_shift;

      if (n_transforms > 0) {
        for (cs_lnum_t ii = 0, jj = n_c_domains*n_transforms*2;
             ii < n_c_domains_r*n_transforms*2;
             ii+=2, jj+=2)
          perio_lst_tr_r[jj] += index_shift;
      }

      n_c_domains += n_c_domains_r;
    }

    /* Order list by rank, transform, and new element number */

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
      int c_rank_id = c_domain_rank_recv[rank_idx];
      /* use -2 instead of -1 for local rank so that division of
         c_rank_id*tr_mult + tr_id by tr_mult is < 0 */
      if (c_rank_id == cs_glob_rank_id) {
        c_rank_id = -2;
        c_domain_rank_recv[rank_idx] = -1;
      }

      for (cs_lnum_t ii = index[rank_idx*2];
           ii < index[rank_idx*2+2];
           ii++) {
        tmp_num[ii*2] = c_rank_id * tr_mult;
      }
    }

    CS_FREE(index);

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {
      tmp_num[ii*2 + 1] = new_src_cell_id[ii];
    }

    if (n_transforms > 0) {
      for (int rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
        for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
          cs_lnum_t s_id = perio_lst_tr_r[n_transforms*2*rank_id + 2*tr_id];
          cs_lnum_t e_id
            = s_id + perio_lst_tr_r[n_transforms*2*rank_id + 2*tr_id + 1];
          cs_lnum_t tr_num_shift = (tr_id + 1);
          for (cs_lnum_t ii = s_id; ii < e_id; ii++) {
            tmp_num[ii*2] += tr_num_shift;
          }
        }
      }
      CS_FREE(perio_lst_tr_r);
    }

    /* Now order elements based on rank, transform, and id */

    cs_lnum_t *order;
    CS_MALLOC(order, n_elts_ini*2, cs_lnum_t);

    cs_order_lnum_allocated_s(nullptr, tmp_num, 2, order, n_elts_ini);

    /* Now build ordered list and renumbering array */

    int prev_rank_tr_id = -1;
    cs_lnum_t prev_src_id = -1, n_elts_m = 0;

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

    cs_lnum_t send_shift = 0;
    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;
      int n_send_r = recv_count[rank_idx*3 + 2];
      MPI_Isend(halo_o2n + send_shift, n_send_r, CS_MPI_LNUM,
                dist_rank, tag, comm, &request[request_count++]);
      send_shift += n_send_r;
    }

    CS_FREE(order);
    CS_FREE(tmp_num);
    CS_FREE(recv_count);

    /* Now we do not need the unmerged halo data anymore, so can rebuild
       the new (merged) halo */

    const fvm_periodicity_t *periodicity = h->periodicity;
    cs_halo_destroy(&h);

    cs_rank_neighbors_t *rn
      = cs_rank_neighbors_create(n_c_domains, c_domain_rank_recv);

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
    CS_FREE(new_src_cell_id);
    CS_FREE(c_domain_rank_recv);
  }

  CS_FREE(perio_lst_tr_s);

  MPI_Irecv(new_cell_id + g->n_rows, g->n_cols_ext - g->n_rows,
            CS_MPI_LNUM, grm->rank_id, tag, comm,
            &request[request_count++]);

  MPI_Waitall(request_count, request, status);

  CS_FREE(halo_o2n);

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }

  if (is_receiver == false)
    cs_halo_destroy(&h);

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
  auto *grm = g->rank_merge;

  MPI_Comm  comm = cs_glob_mpi_comm;

  const cs_lnum_t db_stride = g->db_size * g->db_size;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'c';

  int n_requests = 1;
  bool is_receiver = (grm->rank_id == cs_glob_rank_id) ? true : false;
  if (is_receiver)
    n_requests += grm->sub_size[0];

  int  request_count = 0;

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  /* Allocate if necessary */

  if (n_requests*6 > 128) {
    CS_MALLOC(request, n_requests*6, MPI_Request);
    CS_MALLOC(status, n_requests*6, MPI_Status);
  }

  /* Reallocate arrays for receiving rank and append data */

  cs_real_t *da = nullptr;
  cs_real_3_t *cell_cen = nullptr;
  cs_real_t *cell_vol = nullptr;

  g->n_rows = 0;
  g->n_cols_ext = 0;

  if (is_receiver) {

    g->n_rows = grm->row_idx[grm->sub_size[0]];
    g->n_cols_ext = g->n_rows + g->halo->n_elts[0];

    if (g->relaxation > 0) {
      CS_MALLOC(cell_cen, g->n_cols_ext, cs_real_3_t);
      CS_MALLOC(cell_vol, g->n_cols_ext, cs_real_t);
    }

    CS_MALLOC(da, g->n_cols_ext*db_stride, cs_real_t);

    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {

      cs_lnum_t n_recv = (  grm->row_idx[rank_idx+1]
                          - grm->row_idx[rank_idx]);
      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;

      if (g->relaxation > 0) {
        MPI_Irecv(cell_cen + grm->row_idx[rank_idx], n_recv*3,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

        MPI_Irecv(cell_vol + grm->row_idx[rank_idx], n_recv,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);
      }

      MPI_Irecv(da + grm->row_idx[rank_idx]*db_stride,
                n_recv*db_stride,
                CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

    }

  }

  {
    cs_lnum_t n_send = g->n_elts_r[0];

    if (g->relaxation > 0) {
      MPI_Isend(g->_cell_cen, n_send*3, CS_MPI_REAL, grm->rank_id,
                tag, comm, &request[request_count++]);

      MPI_Isend(g->_cell_vol, n_send, CS_MPI_REAL, grm->rank_id,
                tag, comm, &request[request_count++]);
    }

    MPI_Isend(g->_da, n_send*db_stride, CS_MPI_REAL,
              grm->rank_id, tag, comm, &request[request_count++]);
  }

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }

  CS_FREE(g->_cell_cen);
  CS_FREE(g->_cell_vol);
  CS_FREE(g->_da);

  g->_cell_cen = cell_cen;
  g->_cell_vol = cell_vol;
  g->cell_cen = g->_cell_cen;
  g->cell_vol = g->_cell_vol;

  g->_da = da;
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
 *   g                  <-> pointer to grid structure being merged
 *   n_faces            <-- number of faces to append
 *   face_list          <-> list of faces to append
 *----------------------------------------------------------------------------*/

static void
_append_face_data(cs_grid_t   *g,
                  cs_lnum_t    n_faces,
                  cs_lnum_t   *face_list)
{
  auto *grm = g->rank_merge;

  cs_lnum_t *recv_count = nullptr;

  MPI_Comm  comm = cs_glob_mpi_comm;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'f';

  int n_requests = 1;
  bool is_receiver = (grm->rank_id == cs_glob_rank_id) ? true : false;
  if (is_receiver)
    n_requests += grm->sub_size[0];

  int  request_count = 0;

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  /* Allocate if necessary */

  if (n_requests*10 > 128) {
    CS_MALLOC(request, n_requests*10, MPI_Request);
    CS_MALLOC(status, n_requests*10, MPI_Status);
  }

  /* Exchange counters needed for concatenation */

  if (is_receiver) {
    CS_MALLOC(recv_count, grm->sub_size[0], cs_lnum_t);
    recv_count[0] = g->n_faces;
    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {
      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;
      MPI_Irecv(recv_count + rank_idx, 1, CS_MPI_LNUM, dist_rank,
                tag, comm, &request[request_count++]);
    }
  }
  MPI_Isend(&n_faces, 1, CS_MPI_LNUM, grm->rank_id, tag, comm,
            &request[request_count++]);

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  /* Reallocate arrays for receiving rank and append data */

  CS_FREE(g->coarse_face);

  cs_lnum_t faces_count = 0;
  cs_lnum_2_t *face_cell = nullptr;
  cs_real_t *xa = nullptr, *xa0 = nullptr, *xa0ij = nullptr;
  cs_real_3_t *face_normal = nullptr;

  if (is_receiver) {

    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++)
      faces_count += recv_count[rank_idx];

    CS_MALLOC(face_cell, faces_count, cs_lnum_2_t);

    if (g->symmetric == true)
      CS_REALLOC(xa, faces_count, cs_real_t);
    else
      CS_REALLOC(xa, faces_count*2, cs_real_t);

    if (g->relaxation > 0) {

      if (g->face_normal != nullptr)
        CS_MALLOC(face_normal, faces_count, cs_real_3_t);

      CS_MALLOC(xa0, faces_count, cs_real_t);
      CS_MALLOC(xa0ij, faces_count*3, cs_real_t);

    }

    faces_count = 0;

    for (int rank_idx = 0; rank_idx < grm->sub_size[0]; rank_idx++) {

      cs_lnum_t n_recv = recv_count[rank_idx];
      int dist_rank = grm->sub_root + grm->rank_stride*rank_idx;

      MPI_Irecv(face_cell + faces_count, n_recv*2,
                CS_MPI_LNUM, dist_rank, tag, comm, &request[request_count++]);

      if (g->symmetric == true)
        MPI_Irecv(xa + faces_count, n_recv,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);
      else
        MPI_Irecv(xa + faces_count*2, n_recv*2,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

      if (g->relaxation > 0) {

        if (g->face_normal != nullptr)
          MPI_Irecv(face_normal + faces_count, n_recv*3,
                    CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

        MPI_Irecv(xa0 + faces_count, n_recv,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

        MPI_Irecv(xa0ij + faces_count*3, n_recv*3,
                  CS_MPI_REAL, dist_rank, tag, comm, &request[request_count++]);

      }

      faces_count += recv_count[rank_idx];
    }

    CS_FREE(recv_count);

  }

  /* Filter face connectivity then send it */
  {
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

    MPI_Isend(g->_face_cell, n_faces*2, CS_MPI_LNUM,
              grm->rank_id, tag, comm, &request[request_count++]);

    if (g->symmetric == true)
      MPI_Isend(g->_xa, n_faces, CS_MPI_REAL,
                grm->rank_id, tag, comm, &request[request_count++]);
    else
      MPI_Isend(g->_xa, n_faces*2, CS_MPI_REAL,
                grm->rank_id, tag, comm, &request[request_count++]);
    if (g->relaxation > 0) {
      if (g->face_normal != nullptr) {
        MPI_Isend(g->face_normal, n_faces*3, CS_MPI_REAL,
                  grm->rank_id, tag, comm, &request[request_count++]);
      }

      MPI_Isend(g->_xa0, n_faces, CS_MPI_REAL,
                grm->rank_id, tag, comm, &request[request_count++]);

      MPI_Isend(g->xa0ij, n_faces*3, CS_MPI_REAL,
                grm->rank_id, tag, comm, &request[request_count++]);
    }
  }

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }

  CS_FREE(g->_face_cell);
  CS_FREE(g->face_normal);
  CS_FREE(g->_xa);
  CS_FREE(g->_xa0);
  CS_FREE(g->xa0ij);

  CS_FREE(recv_count);

  g->n_faces = faces_count;

  g->_face_cell = face_cell;
  g->face_cell = g->_face_cell;
  g->face_normal = face_normal;
  g->_xa = xa;
  g->xa = g->_xa;
  g->_xa0 = xa0;
  g->xa0 = g->_xa0;
  g->xa0ij = xa0ij;
}

/*----------------------------------------------------------------------------
 * Compute stride and ranks when merging grids on several ranks to a smaller
 * number of ranks, and associate rank merge object with grid.
 *
 * parameters:
 *   g           <-- Pointer to grid structure
 *   sub_size    <-- Rank subset size
 *   rank_stride <-- Associated destination rank stride (automatic if < 1)
 *----------------------------------------------------------------------------*/

static void
_compute_merge_stride_and_ranks(cs_grid_t  *g,
                                int         sub_size,
                                int         rank_stride = -1)
{
  assert(g->rank_merge == nullptr);

  /* Communicator creation is collective, so must be always called. */

  MPI_Comm parent_comm = g->comm;
  int parent_n_ranks = g->n_ranks;
  int parent_rank_stride = g->rank_stride;

  {
    /* Initial rank stride should match parent, but be safe. */
    if (g->parent != nullptr) {
      parent_rank_stride = g->parent->rank_stride;
      parent_n_ranks = g->parent->n_ranks;
      parent_comm = g->parent->comm;
    }

    int n_active_ranks = parent_n_ranks / sub_size;
    if (parent_n_ranks % sub_size)
      n_active_ranks += 1;

    if (rank_stride == -1)
      rank_stride = parent_rank_stride*sub_size;

    g->comm = cs_base_get_sub_comm(0,
                                   n_active_ranks*rank_stride,
                                   rank_stride);
  }

  if (parent_comm == MPI_COMM_NULL) {
    g->n_ranks = 1;
    return;
  }

  g->rank_merge = new cs::grid_rank_merge;

  auto *grm = g->rank_merge;

  /* Grid should be initialized with parent values;
     if parent is available, ensure this; otherwise assume it. */

  int n_ranks = parent_n_ranks;

  if (rank_stride < 1)
    g->rank_stride = parent_rank_stride * sub_size;
  else
    g->rank_stride = rank_stride;
  g->n_ranks = n_ranks / sub_size;
  if (n_ranks % sub_size > 0)
    g->n_ranks += 1;

  int base_rank = cs_glob_rank_id;
  int sub_id = (base_rank / parent_rank_stride) / sub_size;

  grm->rank_id = sub_id * g->rank_stride;
  grm->rank_stride = parent_rank_stride;
  grm->sub_root = sub_id * parent_rank_stride * sub_size;
  grm->sub_rank = (base_rank / parent_rank_stride) % sub_size;
  grm->sub_size[0] = sub_size;
  grm->sub_size[1] = sub_size;

  if (  grm->sub_root + grm->sub_size[0]*parent_rank_stride
      > parent_n_ranks*parent_rank_stride)
    grm->sub_size[0] = (parent_n_ranks*parent_rank_stride - grm->sub_root)
                       / parent_rank_stride;
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

  if (alloc_mode > CS_ALLOC_HOST) {
    alloc_mode = CS_ALLOC_HOST;
    cs_matrix_set_alloc_mode(const_cast<cs_matrix_t *>(g->_matrix),
                             alloc_mode);
    if (g->matrix_struct != nullptr)
      cs_matrix_structure_set_alloc_mode(g->matrix_struct, alloc_mode);
  }

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
  CS_MALLOC(c_row_index, n_rows+1, cs_lnum_t);
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
  CS_PROFILE_FUNC_RANGE();

  if (g->db_size > 1 || g->level == 0)
    return;

  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  const cs_lnum_t n_rows = g->n_rows;
  cs_lnum_t base_shift = 0;

  int base_rank = 0, n_ranks = 1;

#if defined(HAVE_MPI)

  /* Contrary to _compute_merge_stride_and_ranks, _merge_bottom_grids_to_single
     operates on an existing bottom grid g whose rows, halo, and
     matrix are already partitioned across g->comm, so directly use g->comm. */

  MPI_Comm comm = g->comm;

  if (comm != MPI_COMM_NULL) {
    MPI_Comm_rank(comm, &base_rank);
    MPI_Comm_size(comm, &n_ranks);
    MPI_Exscan(&n_rows, &base_shift, 1, CS_MPI_LNUM, MPI_SUM, comm);

    /* According to standard, base_shift on rank 0 is undefined
       after MPI_Exscan. */
    if (base_rank == 0)
      base_shift = 0;
  }
#endif

  /* Determine rank in merged group */

#if defined(HAVE_MPI)
  if (n_ranks > 1) {
    _compute_merge_stride_and_ranks(g, n_ranks);

    if (verbosity > 2) {
      bft_printf("\n"
                 "    merging level %2d grid:\n"
                 "      merge_stride = %d; new ranks = %d\n",
                 g->level, g->rank_stride, g->n_ranks);
      if (verbosity > 3)
        bft_printf("    rank_merge:\n"
                   "      rank_id = %d\n"
                   "      rank_stride = %d\n"
                   "      sub_root = %d\n"
                   "      sub_rank = %d\n"
                   "      sub_size = {%d, %d}\n",
                   g->rank_merge->rank_id, g->rank_stride,
                   g->rank_merge->sub_root, g->rank_merge->sub_rank,
                   g->rank_merge->sub_size[0], g->rank_merge->sub_size[1]);
    }
  }

#else

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid periodic columns\n",
               g->level);
  }

#endif

  if (g->halo == nullptr)
    return;

  cs_alloc_mode_t alloc_mode = CS_ALLOC_HOST;
  if (g->matrix != nullptr)
    alloc_mode = cs_matrix_get_alloc_mode(g->matrix);

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
  g->face_cell = nullptr;
  g->n_faces = 0;
  CS_FREE(g->cell_face);
  CS_FREE(g->cell_face_sgn);

  /* Compute number of coarse cells in new grid;
     local cell numbers will be shifted by cell_shift in the merged grid */

#if defined(HAVE_MPI)

  if (n_ranks > 1) {

    cs_lnum_t *rank_index = nullptr;

    cs_lnum_t counts[2] = {g->n_rows, nnz};

    auto *grm = g->rank_merge;

    if (base_rank == 0) {
      CS_MALLOC(grm->row_idx, grm->sub_size[0] + 1, cs_lnum_t);
      CS_MALLOC(rank_index, (n_ranks+1)*2, cs_lnum_t);
      rank_index[0] = 0; rank_index[1] = 0;
    }

    if (comm != MPI_COMM_NULL) {
      MPI_Gather(counts, 2, CS_MPI_LNUM, rank_index+2, 2, CS_MPI_LNUM,
                 0, comm);
    }

    if (base_rank == 0) {
      rank_index[2] = g->n_rows; rank_index[3] = nnz;
      grm->row_idx[0] = 0; grm->row_idx[1] = g->n_rows;
      for (int i = 1; i < n_ranks; i++) {
        rank_index[(i+1)*2] += rank_index[i*2];
        rank_index[(i+1)*2+1] += rank_index[i*2 + 1];
        grm->row_idx[i + 1] = rank_index[(i+1)*2];
      }
    }

    /* Now return computed info to grids that will be merged */

    if (comm != MPI_COMM_NULL) {
      MPI_Scatter(rank_index, 2, CS_MPI_LNUM, counts, 2, CS_MPI_LNUM,
                  0, comm);
    }

    /* Now update row index based on local shifts */

    if (base_rank > 0 && row_index != nullptr) {
      for (cs_lnum_t j = 0; j < n_rows+1; j++)
        row_index[j] += counts[1];
    }

    /* Now send matrix data to root rank */

    int *recvcounts = nullptr, *displs = nullptr;
    int n_send = 0;

    cs_lnum_t nnz_recv = 0;

    if (base_rank == 0) {
      g->n_rows = rank_index[n_ranks*2];
      g->n_cols_ext = rank_index[n_ranks*2];
      nnz_recv = rank_index[n_ranks*2 + 1];

      assert(g->n_g_rows == (cs_gnum_t)(g->n_rows));

      CS_REALLOC(row_index, g->n_rows + 1, cs_lnum_t);
      CS_REALLOC(col_id, nnz_recv, cs_lnum_t);
      CS_REALLOC(d_val, g->n_rows, cs_real_t);
      CS_REALLOC(x_val, nnz_recv, cs_real_t);

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

    if (comm != MPI_COMM_NULL) {
      cs_lnum_t *s_row_idx = (row_index != nullptr) ? row_index + 1 : nullptr;
      MPI_Gatherv(s_row_idx, n_send, CS_MPI_LNUM,
                  row_index+1, recvcounts, displs, CS_MPI_LNUM,
                  0, comm);

      MPI_Gatherv(d_val, n_send, CS_MPI_REAL,
                  d_val, recvcounts, displs, CS_MPI_REAL,
                  0, comm);
    }

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

    if (comm != MPI_COMM_NULL) {
      MPI_Gatherv(col_id, n_send, CS_MPI_LNUM,
                  col_id, recvcounts, displs, CS_MPI_LNUM,
                  0, comm);

      MPI_Gatherv(x_val, n_send, CS_MPI_REAL,
                  x_val, recvcounts, displs, CS_MPI_REAL,
                  0, comm);
    }

    CS_FREE(recvcounts);
    CS_FREE(displs);

    if (base_rank > 0) {
      g->n_rows = 0;
      g->n_cols_ext = 0;

      CS_FREE(row_index);
      CS_MALLOC_HD(row_index, 1, cs_lnum_t, alloc_mode);
      row_index[0] = 0;
      CS_FREE(col_id);
      CS_FREE(d_val);
      CS_FREE(x_val);
    }
    else if (alloc_mode > CS_ALLOC_HOST) {
      CS_REALLOC_HD(row_index, g->n_rows + 1, cs_lnum_t, alloc_mode);
      CS_REALLOC_HD(col_id, nnz_recv, cs_lnum_t, alloc_mode);
      CS_REALLOC_HD(d_val, g->n_rows, cs_real_t, alloc_mode);
      CS_REALLOC_HD(x_val, nnz_recv, cs_real_t, alloc_mode);
    }

  } // n_ranks > 1

#endif  // HAVE_MPI

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

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * cs::grid_rank_merge custom new and delete operators
 *----------------------------------------------------------------------------*/

void *
cs::grid_rank_merge::operator new
(
  std::size_t size   /*!<[in] Size to allocate */
)
{
  return cs_mem_malloc(1, size, "cs_grid_rank_merge::new",
                       __FILE__, __LINE__);
}

void
cs::grid_rank_merge::operator delete
(
  void *ptr   /*!<[in] Size to allocate */
)
{
  cs_mem_free(ptr, "cs_grid_rank_merge::delete",
              __FILE__, __LINE__);
}

/*----------------------------------------------------------------------------
 * Destructor
 *----------------------------------------------------------------------------*/

cs::grid_rank_merge::~grid_rank_merge()
{
  CS_FREE(this->buf);
  CS_FREE(this->row_idx);
}

/*----------------------------------------------------------------------------
 * Dump grid rank merge structure
 *----------------------------------------------------------------------------*/

void
cs::grid_rank_merge::dump()
{
  bft_printf("\n"
             "  rank_merge:           %p\n"
             "    rank_id:            %d\n"
             "    sub_root:           %d\n"
             "    sub_rank:           %d\n"
             "    sub_size:           {%d, %d}\n",
             (const void *)this, this->rank_id,
             this->sub_root, this->sub_rank,
             this->sub_size[0], this->sub_size[1]);

  if (this->row_idx != nullptr) {
    bft_printf("    row_idx\n");
    for (int i = 0; i < this->sub_size[0] + 1; i++)
      bft_printf("      %ld: %ld\n", (long)i, (long)this->row_idx[i]);
  }

  bft_printf("\n"
             "    buf:                %p\n"
             "    buf_size:           %lu\n",
             (const void *)this->buf,
             (unsigned long)this->buf_size);
}

#endif // HAVE_MPI

/*============================================================================
 * Semi-private function definitions
 *
 * The following functions are intended to be used by the multigrid layer
 * (cs_multigrid.cpp), not directly by the user, so they are no more
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
 *   g           <-- Pointer to grid structure
 *   sub_size    <-- Associated merge subset size
 *   verbosity   <-- verbosity level
 *   rank_stride <-- destination rank stride, or -1 for default
 *----------------------------------------------------------------------------*/

void
cs_grid_merge_ranks(cs_grid_t  *g,
                    int         sub_size,
                    int         verbosity,
                    int         rank_stride)
{
  if (sub_size < 2)
    return;

  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  /* Determine rank in merged group */

  _compute_merge_stride_and_ranks(g, sub_size, rank_stride);

  auto grm = g->rank_merge;

  if (grm == nullptr)
    return;

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid:\n"
               "      sub_size = %d; new ranks = %d\n",
               g->level, g->rank_stride, g->n_ranks);
    if (verbosity > 3)
      bft_printf("      sub_root = %d; sub_rank = %d; sub_size = {%d, %d}\n",
                 grm->sub_root, grm->sub_rank, grm->sub_size[0],
                 grm->sub_size[1]);
  }

  cs_lnum_t cell_shift = 0;
  cs_lnum_t *new_cell_id = nullptr, *face_list = nullptr;
  bool  *halo_cell_flag = nullptr;
  MPI_Comm comm = cs_glob_mpi_comm;

  static const int tag = 'm'+'e'+'r'+'g'+'e';

  int n_requests = 1;
  bool is_receiver = (grm->rank_id == cs_glob_rank_id) ? true : false;
  if (is_receiver)
    n_requests += grm->sub_size[0];

  int  request_count = 0;

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  /* Allocate if necessary */

  if (n_requests > 128) {
    CS_MALLOC(request, n_requests, MPI_Request);
    CS_MALLOC(status, n_requests, MPI_Status);
  }

  /* Compute number of coarse cells in new grid;
     local cell numbers will be shifted by cell_shift in the merged grid */

  if (is_receiver) {
    CS_MALLOC(grm->row_idx, grm->sub_size[0] + 1, cs_lnum_t);
    grm->row_idx[0] = 0;
    for (int i = 0; i < grm->sub_size[0]; i++) {
      int dist_rank = grm->sub_root + grm->rank_stride*i;
      MPI_Irecv(&(grm->row_idx[i + 1]), 1, CS_MPI_LNUM, dist_rank, tag,
                comm, &request[request_count++]);
    }
  }
  {
    cs_lnum_t send_val = g->n_rows;
    MPI_Isend(&send_val, 1, CS_MPI_LNUM, grm->sub_root, tag, comm,
              &request[request_count++]);
  }

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  if (is_receiver) {
    for (int i = 0; i < grm->sub_size[0]; i++)
      grm->row_idx[i + 1] += grm->row_idx[i];

    /* Now return computed info to grids that will be merged */

    if (grm->sub_rank == 0) {
      for (int i = 0; i < grm->sub_size[0]; i++) {
        MPI_Isend(&(grm->row_idx[i]), 1, CS_MPI_LNUM,
                  grm->sub_root + grm->rank_stride*i, tag,
                  comm, &request[request_count++]);
      }
    }

  }
  MPI_Irecv(&cell_shift, 1, CS_MPI_LNUM,
            grm->sub_root, tag, comm, &request[request_count++]);

  MPI_Waitall(request_count, request, status);
  request_count = 0;

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

  CS_MALLOC(face_list, g->n_faces, cs_lnum_t);

  cs_lnum_t n_faces = (grm->sub_rank == 0) ? g->n_faces : 0;

  if (grm->sub_size[0] > 1 && grm->sub_rank > 0 && g->n_faces > 0) {

    cs_lnum_t n_ghost_cells = g->n_cols_ext - g->n_rows;

    /* Mark which faces should be merged: to avoid duplicates, a face
       connected to a cell on a lower rank in the same merge set is
       discarded, as it has already been accounted for by that rank. */

    CS_MALLOC(halo_cell_flag, n_ghost_cells, bool);
    for (cs_lnum_t j = 0; j < n_ghost_cells; j++)
      halo_cell_flag[j] = false;

    for (cs_lnum_t i = 0; i < g->halo->n_c_domains; i++) {
      int rank_id = g->halo->c_domain_rank[i];
      /* Check if rank_id belongs to the same merge group
         and is lower than self. */
      if (   rank_id >= grm->sub_root
          && rank_id < cs_glob_rank_id
          && (rank_id - grm->sub_root) % grm->rank_stride == 0) {
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

  else {

    /* For first rank of merge set, all faces are used */
    cs_lnum_t _n_faces = g->n_faces;
    for (cs_lnum_t face_id = 0; face_id < _n_faces; face_id++)
      face_list[face_id] = face_id;

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

  if (grm->sub_size[0] > 1) {
    _append_cell_data(g);
    _append_face_data(g, n_faces, face_list);
    CS_FREE(g->cell_face);
    CS_FREE(g->cell_face_sgn);
  }
  _sync_merged_cell_data(g);

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }

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
  if (g->n_ranks == 1 && g->halo == nullptr)
    return;

  /* Cannot be called if this level already has merging info;
     in that case, should be called instead of previous merge
     (which may require re-building the cuurent grid from its parent). */

  if (g->rank_merge != nullptr)
    return;

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

  if (n_max_ranks == 1 && g->halo->n_rotations == 0) {
    _merge_bottom_grids_to_single(g, verbosity);
  }

  else {
    int rank_stride = std::max(cs_glob_node_n_ranks / n_max_ranks, 1);
    int sub_size = g->n_ranks / n_max_ranks;
    if (g->n_ranks % n_max_ranks)
      sub_size += 1;
    cs_grid_merge_ranks(g,
                        sub_size,
                        verbosity,
                        rank_stride);
  }

#else

  /* Without MPI, merging can still be useful to replace periodic halo values
     with regular values, simplifying matrix). */

  if (g->halo != nullptr)
    _merge_bottom_grids_to_single(g, verbosity);

#endif
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Gather coarse cell data in case of merged grids
 *
 * Though the grid is <const> here, its `merge_buf` member may be allocated
 * and used here.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_rank_merge_gather_row_data
(
  const cs_grid_t  *g,            /*!< grid structure */
  cs_datatype_t     data_type,    /*!< associated data type */
  size_t            data_stride,  /*!< associated data stride */
  void             *data          /*!< associated data buffer */
)
{
  assert(g != nullptr);

  auto *grm = g->rank_merge;
  if (grm == nullptr)
    return;

  MPI_Comm  comm = cs_glob_mpi_comm;
  int rank_id = cs_glob_rank_id;

  /* If grid merging has taken place, gather coarse data */

  MPI_Datatype mpi_type = cs_datatype_to_mpi[data_type];

  size_t elt_size = cs_datatype_size[data_type] * data_stride;
  cs_byte_t *send_buf = static_cast<cs_byte_t *>(data);
  cs_byte_t *data_p = static_cast<cs_byte_t *>(data);

  constexpr int tag = 'r'+'e'+'s'+'t'+'r'+'i'+'c'+'t';

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  int r_s_idx = -1, r_e_idx = -1;
  int n_requests = 0;
  if (rank_id == grm->rank_id) {
    r_s_idx = (grm->rank_id == grm->sub_root) ? 1 : 0;
    r_e_idx = grm->sub_size[0];
    n_requests = r_e_idx - r_s_idx;
  }

  if (r_s_idx != 1)
    n_requests++;

  /* Allocate if necessary */

  if (n_requests > 128) {
    CS_MALLOC(request, n_requests*2, MPI_Request);
    CS_MALLOC(status, n_requests*2, MPI_Status);
  }

  /* If ranks overlap, an extra buffer is needed */

  if (r_s_idx == 0) {
    cs::grid_rank_merge *_grm = const_cast<cs::grid_rank_merge *>(grm);
    size_t merge_buf_size = (size_t)(g->n_elts_r[0]) * elt_size;
    if (merge_buf_size > _grm->buf_size) {
      CS_FREE(_grm->buf);
      CS_MALLOC(_grm->buf, merge_buf_size, cs_byte_t);
    }
    memcpy(_grm->buf, data, merge_buf_size);
    send_buf = static_cast<cs_byte_t *>(_grm->buf);
  }

  /* Exchange data */

  int request_count = 0;

  for (int rank_idx = r_s_idx; rank_idx < r_e_idx; rank_idx++) {
    cs_lnum_t n_recv = (  grm->row_idx[rank_idx+1]
                        - grm->row_idx[rank_idx]) * data_stride;
    int dist_rank = grm->sub_root + rank_idx*grm->rank_stride;
    size_t s_shift = (size_t)(grm->row_idx[rank_idx]) * elt_size;
    MPI_Irecv(data_p + s_shift, n_recv, mpi_type,
              dist_rank, tag, comm, &request[request_count++]);
  }

  if (r_s_idx != 1) {
    MPI_Isend(send_buf, g->n_elts_r[0]*data_stride, mpi_type,
              grm->rank_id, tag, comm, &request[request_count++]);
  }

  MPI_Waitall(request_count, request, status);

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Scatter coarse cell data in case of merged grids
 *
 * Though the grid is <const> here, its `merge_buf` member may be allocated
 * and used here.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_rank_merge_scatter_row_data
(
  const cs_grid_t  *g,            /*!< grid structure */
  cs_datatype_t     data_type,    /*!< associated data type */
  size_t            data_stride,  /*!< associated data stride */
  void             *data          /*!< associated data buffer */
)
{
  assert(g != nullptr);

  auto *grm = g->rank_merge;
  if (grm == nullptr)
    return;

  MPI_Comm  comm = cs_glob_mpi_comm;
  int rank_id = cs_glob_rank_id;

  /* If grid merging has taken place, scatter coarse data */

  MPI_Datatype mpi_type = cs_datatype_to_mpi[data_type];

  size_t elt_size = cs_datatype_size[data_type] * data_stride;
  cs_byte_t *send_buf = static_cast<cs_byte_t *>(data);

  constexpr int tag = 'p'+'r'+'o'+'l'+'o'+'n'+'g';

  MPI_Request _request[128];
  MPI_Request *request = _request;
  MPI_Status _status[128];
  MPI_Status *status = _status;

  int r_s_idx = -1, r_e_idx = -1;
  int n_requests = 0;
  if (rank_id == grm->rank_id) {
    r_s_idx = (grm->rank_id == grm->sub_root) ? 1 : 0;
    r_e_idx = grm->sub_size[0];
    n_requests = r_e_idx - r_s_idx;
  }

  if (r_s_idx != 1)
    n_requests++;

  /* Allocate if necessary */

  if (n_requests > 128) {
    CS_MALLOC(request, n_requests*2, MPI_Request);
    CS_MALLOC(status, n_requests*2, MPI_Status);
  }

  /* If ranks overlap, an extra buffer is needed */

  if (r_s_idx == 0) {
    cs::grid_rank_merge *_grm = const_cast<cs::grid_rank_merge *>(grm);
    size_t merge_buf_size = (size_t)(g->n_rows) * elt_size;
    if (merge_buf_size > grm->buf_size) {
      CS_FREE(_grm->buf);
      CS_MALLOC(_grm->buf, merge_buf_size, cs_byte_t);
    }
    memcpy(_grm->buf, data, merge_buf_size);
    send_buf = static_cast<cs_byte_t *>(_grm->buf);
  }

  /* Exchange data */

  int request_count = 0;

  if (r_s_idx != 1) {
    MPI_Irecv(data, g->n_elts_r[0]*data_stride, mpi_type,
              grm->rank_id, tag, comm, &request[request_count++]);
  }

  for (int rank_idx = r_s_idx; rank_idx < r_e_idx; rank_idx++) {
    cs_lnum_t n_send = (  grm->row_idx[rank_idx+1]
                        - grm->row_idx[rank_idx]) * data_stride;
    int dist_rank = grm->sub_root + rank_idx*grm->rank_stride;
    size_t s_shift = (size_t)(grm->row_idx[rank_idx]) * elt_size;
    MPI_Isend(send_buf + s_shift, n_send, mpi_type,
              dist_rank, tag, comm, &request[request_count++]);
  }

  MPI_Waitall(request_count, request, status);

  if (request != _request) {
    CS_FREE(request);
    CS_FREE(status);
  }
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
