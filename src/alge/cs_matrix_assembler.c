/*============================================================================
 * Incremental or general construction of matrix.
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
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

#include "cs_defs.h"
#include "cs_rank_neighbors.h"
#include "cs_sort.h"
#include "cs_timer.h"

#include "cs_matrix.h"
#include "cs_matrix_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_assembler.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_matrix_assembler.c

  \brief Incremental or general construction of matrix structure.

  The matrix assembler is intended to assist the building of matrices
  in general parallel conditions, assuming each parallel rank is assigned a
  given number of rows, in increasing order. Column elements may refer to rows
  assigned to other ranks. This allows a good mapping to external libraries
  such as PETSc and Hypre.

  Moreover, in some cases, some global matrix elements may be computed on ranks
  to whom neither the row nor the column is assigned. This may happen for
  example when rows are assigned to mesh vertices, and a given cell
  provides a contribution to two vertices laying on parallel rank boundaries,
  and assigned to two different ranks, which are different from the cell's
  rank. Few such contributions are expected with a good partitionning, but
  they cannot always be avoided, or this would require complex partitioning
  constraints. libraries such as PETSc do not handle these cases, at least
  not efficiently (as the recommended preallocation required for local rows
  cannot be computed without the knowledge of those contributions), so
  the assembler should help manage these contributions in any case.

  The addition of a given non-zero to a matrix structure may be done
  multiple times (for example when looping on cells with an internal loop on
  vertices), but for better performance (and memory usage), it is recommended
  to avoid the same info multiple times, as duplicates may be removed only
  during the computation stage).
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Minimum size for OpenMP loops (needs benchmarking to adjust) */
#define THR_MIN 128

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure used to pre-build a matrix
 *----------------------------------------------------------------------------*/

struct _cs_matrix_assembler_t {

  bool        separate_diag;     /* is diagonal handled separately ? */

  bool        external_lib_ids;  /* if true, use global column ids instead
                                    of local ranks and ids (for use with
                                    external libraries such as PETSc, Hypre,
                                    or Lis) */

  cs_gnum_t   l_range[2];        /* local global row range */
  cs_gnum_t   n_g_rows;          /* global number of rows */
  cs_lnum_t   n_rows;            /* local number of rows */

  cs_lnum_t   size;              /* current insertion array size */
  cs_lnum_t   max_size;          /* maximum insertion array size */

  cs_lnum_t   *r_idx;            /* main row index (0 to n-1) */
  cs_lnum_t   *c_id;             /* main column ids (0 to n-1) */

  cs_gnum_t   *g_rc_id;          /* global row and column ids
                                    (local and distant) */

#if defined(HAVE_MPI)

  /* Distant columns associated with local rows */

  cs_lnum_t   *d_r_idx;         /* distant row index (0 to n-1) */
  cs_gnum_t   *d_g_c_id;        /* distant global column ids (0 to n-1) */

  /* Metadata for exchange of matrix coefficient values with other ranks */

  int           n_coeff_ranks;           /* number of MPI ranks with which
                                            coefficients are exchanged */
  int          *coeff_rank;              /* ranks with which coefficients are
                                            exchanged */

  cs_lnum_t     coeff_send_size;         /* number of coefficients to send */
  cs_lnum_t     coeff_recv_size;         /* number of coefficients to receive */

  cs_lnum_t     coeff_send_n_rows;       /* number of matching rows */
  cs_lnum_t    *coeff_send_index;        /* index of sent coefficient rows */
  cs_gnum_t    *coeff_send_row_g_id;     /* global ids matching rows (ordered) */
  cs_gnum_t    *coeff_send_col_g_id;     /* global ids matching columns
                                            (ordered) */

  cs_lnum_t    *coeff_rank_send_index;   /* index of data to send */
  cs_lnum_t    *coeff_rank_recv_index;   /* index of data to receive */

  cs_lnum_t    *coeff_recv_row_id;       /* local row ids associated with
                                            received data; */
  cs_lnum_t    *coeff_recv_col_id;       /* local column id associated with
                                            received data; the column numbering
                                            implicitely assumes local terms
                                            first, distant ones second */
  cs_gnum_t    *coeff_recv_col_g_id;     /* global column id couples
                                            associated with received data */

  /* Associated vector ghost element info */

  cs_halo_t    *halo;                    /* halo for associated vectors */

  /* Associated communicator */

  MPI_Comm     comm;             /* associated MPI communicator */

#endif /* HAVE_MPI */

  /* Statistics */

  int          n_ranks_init[2];  /* Number of ranks for initial exchange
                                    for distant rows then columns */

};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sort and compact local matrix elements.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_sort_and_compact_local(cs_matrix_assembler_t  *ma)
{
  /* Check if we are not already sorted */

  cs_lnum_t n_rows = ma->n_rows;
  bool ordered = true;

  for (cs_lnum_t i = 0; i < n_rows && ordered; i++) {
    cs_lnum_t *col_id = ma->c_id + ma->r_idx[i];
    cs_lnum_t n_cols = ma->r_idx[i+1] - ma->r_idx[i];
    for (cs_lnum_t j = 1; j < n_cols; j++) {
      if (col_id[j] <= col_id[j-1])
        ordered = false;
    }
  }

  if (ordered)
    return;

  /* Sort by row */

  bool direct_assembly = cs_sort_indexed(n_rows, ma->r_idx, ma->c_id);

  /* Compact elements if necessary */

  if (direct_assembly == false) {

    cs_lnum_t *tmpr_idx = NULL;

    BFT_MALLOC(tmpr_idx, n_rows+1, cs_lnum_t);
    memcpy(tmpr_idx, ma->r_idx, (n_rows+1)*sizeof(cs_lnum_t));

    cs_lnum_t  k = 0;

    for (cs_lnum_t i = 0; i < n_rows; i++) {
      cs_lnum_t *col_id = ma->c_id + ma->r_idx[i];
      cs_lnum_t n_cols = ma->r_idx[i+1] - ma->r_idx[i];
      cs_lnum_t col_id_prev = -1;
      ma->r_idx[i] = k;
      for (cs_lnum_t j = 0; j < n_cols; j++) {
        if (col_id_prev != col_id[j]) {
          ma->c_id[k++] = col_id[j];
          col_id_prev = col_id[j];
        }
      }
    }
    ma->r_idx[n_rows] = k;

    assert(ma->r_idx[n_rows] < tmpr_idx[n_rows]);

    BFT_FREE(tmpr_idx);
    BFT_REALLOC(ma->c_id, (ma->r_idx[n_rows]), cs_lnum_t);

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sort and compact distant matrix elements.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_sort_and_compact_distant(cs_matrix_assembler_t  *ma)
{
  /* Check if we are not already sorted */

  cs_lnum_t n_rows = ma->n_rows;
  bool ordered = true;

  for (cs_lnum_t i = 0; i < n_rows && ordered; i++) {
    cs_gnum_t *col_id = ma->d_g_c_id + ma->r_idx[i];
    cs_lnum_t n_cols = ma->r_idx[i+1] - ma->r_idx[i];
    for (cs_lnum_t j = 1; j < n_cols; j++) {
      if (col_id[j] <= col_id[j-1])
        ordered = false;
    }
  }

  if (ordered)
    return;

  /* Sort by row */

  bool direct_assembly
    = cs_sort_indexed_gnum(n_rows, ma->d_r_idx, ma->d_g_c_id);

  /* Compact elements if necessary */

  if (direct_assembly == false) {

    cs_lnum_t *tmpr_idx = NULL;

    BFT_MALLOC(tmpr_idx, n_rows+1, cs_lnum_t);
    memcpy(tmpr_idx, ma->d_r_idx, (n_rows+1)*sizeof(cs_lnum_t));

    cs_lnum_t  k = 0;

    for (cs_lnum_t i = 0; i < n_rows; i++) {
      cs_gnum_t *g_col_id = ma->d_g_c_id + ma->d_r_idx[i];
      cs_lnum_t n_cols = ma->d_r_idx[i+1] - ma->d_r_idx[i];
      ma->d_r_idx[i] = k;
      if (n_cols > 0)
        ma->d_g_c_id[k++] = g_col_id[0];
      for (cs_lnum_t j = 1; j < n_cols; j++) {
        if (g_col_id[j] != g_col_id[j-1])
          ma->d_g_c_id[k++] = g_col_id[j];
      }
    }
    ma->d_r_idx[n_rows] = k;

    assert(ma->d_r_idx[n_rows] < tmpr_idx[n_rows]);

    BFT_FREE(tmpr_idx);
    BFT_REALLOC(ma->d_g_c_id, (ma->d_r_idx[n_rows]), cs_gnum_t);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add terms to distant  matrix elements.
 *
 * \param[in, out]  ma    pointer to matrix assembler structure
 * \param[in]       n     number of added column and row couples
 * \param[in]       g_ij  added couples
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_complete_distant(cs_matrix_assembler_t  *ma,
                  cs_lnum_t               n,
                  cs_gnum_t               g_ij[])
{
  /* Count maximum terms to add for each row */

  cs_lnum_t n_rows = ma->n_rows;

  cs_lnum_t *d_c_count, *d_r_idx;
  BFT_MALLOC(d_c_count, n_rows, cs_lnum_t);
  BFT_MALLOC(d_r_idx, n_rows+1, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++)
    d_c_count[i] = 0;

  /* Only local rows should be encountered here; local column elements
     are ignored here (we will later check they do not require a local
     insertion, which is not expected, as we assume the local rank
     should already have provided a contribution to a local entry) */

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_gnum_t g_r_id = g_ij[i*2];
    cs_gnum_t g_c_id = g_ij[i*2 + 1];
    assert (g_r_id >=  ma->l_range[0] && g_r_id < ma->l_range[1]);
    if (g_c_id <  ma->l_range[0] || g_c_id >= ma->l_range[1])
      d_c_count[g_r_id - ma->l_range[0]] += 1;
  }

  d_r_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++)
    d_r_idx[i+1] = d_r_idx[i] + d_c_count[i];

  /* Expand matrix, starting copies from end to avoid overwrites
     (first line untouched) */

  BFT_REALLOC(ma->d_g_c_id, ma->d_r_idx[n_rows] + d_r_idx[n_rows], cs_gnum_t);

  for (cs_lnum_t i = n_rows-1; i > 0; i--) {
    cs_lnum_t n_cols = ma->d_r_idx[i+1] - ma->d_r_idx[i];
    cs_gnum_t *g_col_id_d = ma->d_g_c_id + ma->d_r_idx[i] + d_r_idx[i];
    const cs_gnum_t *g_col_id_s = ma->d_g_c_id + ma->d_r_idx[i];
    d_c_count[i] = n_cols;
    ma->d_r_idx[i+1] += d_r_idx[i+1];
    for (cs_lnum_t j = n_cols-1; j >= 0; j--)
      g_col_id_d[j] = g_col_id_s[j];
  }
  d_c_count[0] = ma->d_r_idx[1];
  ma->d_r_idx[1] += d_r_idx[1];

  BFT_FREE(d_r_idx);

  /* Now add terms */

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_gnum_t g_r_id = g_ij[i*2];
    cs_gnum_t g_c_id = g_ij[i*2 + 1];
    if (g_c_id <  ma->l_range[0] || g_c_id >= ma->l_range[1]) {
      cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
      ma->d_g_c_id[ma->d_r_idx[l_r_id] + d_c_count[l_r_id]] = g_c_id;
      d_c_count[l_r_id] += 1;
    }
  }
  BFT_FREE(d_c_count);

  /* Sort and remove duplicates */

  _sort_and_compact_distant(ma);
}

/*----------------------------------------------------------------------------
 * Determine assumed rank associated with a given global id.
 *
 * parameters:
 *   n_ranks <-- number of ranks
 *   n_g     <-- global number of ids
 *   g_id    <-- array of global ids for local rank (size: n)
 *
 * returns:
 *   assumed rank id for g_id
 *----------------------------------------------------------------------------*/

static inline int
_assumed_rank(int              n_ranks,
              cs_gnum_t        n_g,
              const cs_gnum_t  g_id)
{
  int rank_id;

  cs_gnum_t n_per_rank = n_g / n_ranks;
  cs_lnum_t rmdr = n_g - n_per_rank * (cs_gnum_t)n_ranks;

  assert(g_id < n_g);

  if (rmdr == 0)
    rank_id = g_id / n_per_rank;
  else {
    cs_gnum_t n_ranks_rmdr = n_ranks - rmdr;
    cs_gnum_t n_ranks_n_per_rank = n_ranks_rmdr * n_per_rank;
    if (g_id  <  n_ranks_n_per_rank)
      rank_id = g_id / n_per_rank;
    else
      rank_id = (g_id + n_ranks_rmdr) / (n_per_rank + 1);
  }

  return rank_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build Assumed rank neighbors info
 *
 * This operation is collective on communicator comm.
 *
 * \param[in]  n        number of global ids for local rank
 * \param[in]  n_g      global number of ids
 * \param[in]  l_range  global id range [min, max[ for local rank
 * \param[in]  g_id     array of global ids for local rank (size: n*stride)
 * \param[in]  comm     associated communicator
 *
 * \return  assumed rank neighbors structure
 */
/*----------------------------------------------------------------------------*/

static cs_rank_neighbors_t *
_assumed_rank_neighbors(cs_lnum_t        n,
                        cs_gnum_t        n_g,
                        cs_gnum_t        l_range[2],
                        const cs_gnum_t  g_id[],
                        MPI_Comm         comm)
{
  /* Prepare for determination of assumed rank;
     we will send both a partial local range description and global ids to the
     assumed partition, using a common array, so as to minimize
     all-to-all communication. */

  int n_ranks, l_rank;
  int *a_rank;

  MPI_Comm_size(comm, &n_ranks);
  MPI_Comm_rank(comm, &l_rank);

  /* Build real to assumed partition mapping */

  cs_lnum_t l_range_size = l_range[1] - l_range[0];

  BFT_MALLOC(a_rank, l_range_size + n, int);

  cs_lnum_t n_a_neighbors = 0, count = 0;

  /* Ranks to send real partition info to */

  int r_prev = -1;

  for (cs_lnum_t i = 0; i < l_range_size; i++) {
    cs_gnum_t _g_id = l_range[0] + (cs_gnum_t)i;
    int r_dest = _assumed_rank(n_ranks, n_g, _g_id);
    if (r_dest != r_prev && r_dest != l_rank) {
      a_rank[n_a_neighbors++] = r_dest;
      r_prev = r_dest;
    }
  }

  count = n_a_neighbors;

  /* Ranks to send real data to; assuming ids are mostly
     ordered, we limit the size of sent data. */

  for (cs_lnum_t i = 0; i < n; i++) {
    int r_dest = _assumed_rank(n_ranks, n_g, g_id[i]);
    if (r_dest != r_prev && r_dest != l_rank) {
      a_rank[count++] = r_dest;
      r_prev = r_dest;
    }
  }

  /* Build communicating rank neighbors structure;
     symmetrize the structure to tell ranks from which other
     ranks they may receive info */

  cs_rank_neighbors_t *arn = cs_rank_neighbors_create(count, a_rank);

  BFT_FREE(a_rank);

  cs_rank_neighbors_symmetrize(arn, comm);

  return arn;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build Assumed rank-helvetica-medium- neighbors info
 *
 * This operation is collective on communicator comm. The caller is
 * responsible for freeng the returned array.
 *
 * \param[in]  arn      info on ranks with which we will communicate
 * \param[in]  l_range  global id range [min, max[ for local rank
 * \param[in]  comm     associated communicator
 *
 * \return  true rank neighbors ranges array (size: arn->size*2)
 */
/*----------------------------------------------------------------------------*/

static cs_gnum_t *
_assumed_rank_exchange(cs_rank_neighbors_t  *arn,
                       cs_gnum_t             l_range[2],
                       MPI_Comm              comm)
{
  MPI_Request *request = NULL;
  MPI_Status *status = NULL;
  cs_gnum_t *d_ranges = NULL;

  BFT_MALLOC(d_ranges, arn->size*2, cs_gnum_t);
  BFT_MALLOC(request, arn->size*2, MPI_Request);
  BFT_MALLOC(status, arn->size*2, MPI_Status);

  /* Prepare for determination of assumed rank;
     we will send both a partial local range description and global ids to the
     assumed partition, using a common array, so as to minimize
     all-to-all communication. */

  /* Exchange local range with neighbor ranks */

  int request_count = 0;
  const int local_rank = cs_glob_rank_id;

  for (int i = 0; i < arn->size; i++) {
    MPI_Irecv(d_ranges + i*2,
              2,
              CS_MPI_GNUM,
              arn->rank[i],
              local_rank,
              comm,
              &(request[request_count++]));
  }

  for (int i = 0; i < arn->size; i++) {
    MPI_Isend(l_range,
              2,
              CS_MPI_GNUM,
              arn->rank[i],
              arn->rank[i],
              comm,
              &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);

  return d_ranges;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Use binary search to determine to which range a global id belongs
 *
 * \param[in]   n_ranges  number of distant ranges
 * \param[in]   d_ranges  global id ranges [min, max[ for neighboring ranks
 * \param[in ]  g_id      global id to search for
 *
 * \return  id of matching distant range, or -2 if not found
 */
/*----------------------------------------------------------------------------*/

static inline int
_g_id_rank(int              n_ranges,
           const cs_gnum_t  d_ranges[],
           cs_gnum_t        g_id)
{
  int start_id = 0;
  int end_id = n_ranges - 1;
  int mid_id = (end_id -start_id) / 2;

  while (start_id <= end_id) {
    if (d_ranges[mid_id*2+1] <= g_id)
      start_id = mid_id + 1;
    else if (d_ranges[mid_id*2] > g_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  if (g_id < d_ranges[mid_id*2] || g_id >= d_ranges[mid_id*2+1])
    mid_id = -2;

  return mid_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Process and exchange info relative to data that will
 *        be sent to neighboring ranks.
 *
 * At this stage, ma->rank will contain the index of the rank in the
 * rank neighborhood structure, not the true MPI rank.
 *
 * \param[in, out]  ma           pointer to matrix assembler structure
 * \param[in]       e_g_ij_size  size of coefficient data to send
 * \param[in]       e_g_ij       coefficient data (g_row_id, g_col_id couples)
 *                               to send
 */
/*----------------------------------------------------------------------------*/

static void
_process_assembly_data(cs_matrix_assembler_t  *ma,
                       cs_gnum_t               e_g_ij_size,
                       const cs_gnum_t        *e_g_ij)
{
  ma->coeff_send_size = e_g_ij_size;

  /* Count rows */

  ma->coeff_send_n_rows = 0;

  cs_gnum_t g_r_id_prev = ma->l_range[0]; /* impossible value here */

  for (cs_lnum_t i = 0; i < ma->coeff_send_size; i++) {
    cs_gnum_t g_r_id = e_g_ij[i*2];
    if (g_r_id != g_r_id_prev) {
      ma->coeff_send_n_rows += 1;
      g_r_id_prev = g_r_id;
    }
  }

  if (ma->coeff_send_size > 0) {

    cs_lnum_t row_count = 0;

    BFT_MALLOC(ma->coeff_send_index, ma->coeff_send_n_rows+1, cs_lnum_t);
    BFT_MALLOC(ma->coeff_send_row_g_id, ma->coeff_send_n_rows, cs_gnum_t);
    BFT_MALLOC(ma->coeff_send_col_g_id, ma->coeff_send_size, cs_gnum_t);

    g_r_id_prev = ma->l_range[0]; /* impossible value here */

    for (cs_lnum_t i = 0; i < ma->coeff_send_size; i++) {
      cs_gnum_t g_r_id = e_g_ij[i*2];
      if (g_r_id != g_r_id_prev) {
        ma->coeff_send_index[row_count] = i;
        ma->coeff_send_row_g_id[row_count] = g_r_id;
        row_count++;
        g_r_id_prev = g_r_id;
      }
      cs_gnum_t g_c_id = e_g_ij[i*2+1];
      ma->coeff_send_col_g_id[i] = g_c_id;
    }
    ma->coeff_send_index[row_count] = ma->coeff_send_size;

  }

  /* Determine ranks we may be communicating with for this stage */

  cs_rank_neighbors_t *arn = _assumed_rank_neighbors(ma->coeff_send_n_rows,
                                                     ma->n_g_rows,
                                                     ma->l_range,
                                                     ma->coeff_send_row_g_id,
                                                     ma->comm);

  ma->n_ranks_init[0] = arn->size;

  cs_gnum_t *d_ranges = _assumed_rank_exchange(arn,
                                               ma->l_range,
                                               ma->comm);

  cs_lnum_t *counts;
  BFT_MALLOC(counts, arn->size*2, cs_lnum_t);

  for (int i = 0; i < arn->size*2; i++)
    counts[i] = 0;

  for (cs_lnum_t i = 0; i < ma->coeff_send_n_rows; i++) {
    cs_gnum_t g_r_id = ma->coeff_send_row_g_id[i];
    int r_rank_id = _g_id_rank(arn->size, d_ranges, g_r_id);
    assert(r_rank_id >= 0 && r_rank_id < arn->size);
    assert(   g_r_id >= d_ranges[2*r_rank_id]
           && g_r_id < d_ranges[2*r_rank_id+1]);
    counts[r_rank_id] += ma->coeff_send_index[i+1] - ma->coeff_send_index[i];
  }

  BFT_FREE(d_ranges);

  /* Exchange counts */

  MPI_Request *request = NULL;
  MPI_Status *status = NULL;

  BFT_MALLOC(request, arn->size*2, MPI_Request);
  BFT_MALLOC(status, arn->size*2, MPI_Status);

  const int local_rank = cs_glob_rank_id;

  int request_count = 0;

  for (int i = 0; i < arn->size; i++) {
    MPI_Irecv(counts + arn->size + i,
              1,
              CS_MPI_LNUM,
              arn->rank[i],
              local_rank,
              ma->comm,
              &(request[request_count++]));
  }

  for (int i = 0; i < arn->size; i++) {
    MPI_Isend(counts + i,
              1,
              CS_MPI_LNUM,
              arn->rank[i],
              arn->rank[i],
              ma->comm,
              &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  /* Save exchange rank info */

  ma->n_coeff_ranks = 0;
  for (int i = 0; i < arn->size; i++) {
    if (counts[i] + counts[arn->size + i] > 0)
      ma->n_coeff_ranks += 1;
  }

  BFT_MALLOC(ma->coeff_rank, ma->n_coeff_ranks, int);

  if (ma->n_coeff_ranks > 0) {

    BFT_MALLOC(ma->coeff_rank_send_index, ma->n_coeff_ranks+1, cs_lnum_t);
    BFT_MALLOC(ma->coeff_rank_recv_index, ma->n_coeff_ranks+1, cs_lnum_t);

    ma->coeff_rank_send_index[0] = 0;
    ma->coeff_rank_recv_index[0] = 0;

    ma->n_coeff_ranks = 0;
    for (int i = 0; i < arn->size; i++) {
      if (counts[i] + counts[arn->size + i] > 0) {
        int j = ma->n_coeff_ranks;
        ma->coeff_rank[j] = arn->rank[i];
        ma->coeff_rank_send_index[j+1]
          = ma->coeff_rank_send_index[j] + counts[i];
        ma->coeff_rank_recv_index[j+1]
          = ma->coeff_rank_recv_index[j] + counts[arn->size + i];
        ma->n_coeff_ranks += 1;
      }
    }

  }

  BFT_FREE(counts);

  /* We don't need the ranks neighbors structure anymore
     (we'll use the saved and compacted info from now on) */

  cs_rank_neighbors_destroy(&arn);

  if (ma->n_coeff_ranks > 0) {

    /* Exchange data */

    ma->coeff_recv_size = ma->coeff_rank_recv_index[ma->n_coeff_ranks];

    cs_gnum_t  *recv_data;
    BFT_MALLOC(recv_data, ma->coeff_recv_size*2, cs_gnum_t);

    request_count = 0;

    for (int i = 0; i < ma->n_coeff_ranks; i++) {
      cs_lnum_t l_size =   ma->coeff_rank_recv_index[i+1]
                         - ma->coeff_rank_recv_index[i];
      if (l_size > 0) {
        cs_lnum_t recv_shift = ma->coeff_rank_recv_index[i]*2;
        MPI_Irecv(recv_data + recv_shift,
                  l_size*2,
                  CS_MPI_GNUM,
                  ma->coeff_rank[i],
                  local_rank,
                  ma->comm,
                  &(request[request_count++]));
      }
    }

    for (int i = 0; i < ma->n_coeff_ranks; i++) {
      cs_lnum_t l_size =   ma->coeff_rank_send_index[i+1]
                         - ma->coeff_rank_send_index[i];
      if (l_size > 0) {
        cs_lnum_t send_shift = ma->coeff_rank_send_index[i]*2;
        MPI_Isend(e_g_ij + send_shift,
                  l_size*2,
                  CS_MPI_GNUM,
                  ma->coeff_rank[i],
                  ma->coeff_rank[i],
                  ma->comm,
                  &(request[request_count++]));
      }
    }

    MPI_Waitall(request_count, request, status);

    /* Insert terms at local rows */

    if (ma->coeff_recv_size > 0)
      _complete_distant(ma, ma->coeff_recv_size, recv_data);

    /* Now we build the mapping from received terms to their positions
       in the array. For indexing, we implicitely consider that each row
       is built of 2 sub-rows, with the local columns first, and distant
       columns second. */

    BFT_MALLOC(ma->coeff_recv_row_id, ma->coeff_recv_size, cs_lnum_t);

    if (ma->external_lib_ids == false) {

      BFT_MALLOC(ma->coeff_recv_col_id, ma->coeff_recv_size, cs_lnum_t);

      for (cs_lnum_t i = 0; i < ma->coeff_recv_size; i++) {

        cs_lnum_t l_r_id = recv_data[i*2] - ma->l_range[0];
        cs_gnum_t g_c_id = recv_data[i*2+1];

        ma->coeff_recv_row_id[i] = l_r_id;

        /* Local part */

        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {
          const cs_lnum_t *col_id = ma->c_id + ma->r_idx[l_r_id];
          cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
          cs_lnum_t start_id = 0;
          cs_lnum_t end_id = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id] - 1;
          cs_lnum_t mid_id = (end_id -start_id) / 2;
          while (start_id <= end_id) {
            if (col_id[mid_id] < l_c_id)
              start_id = mid_id + 1;
            else if (col_id[mid_id] > l_c_id)
              end_id = mid_id - 1;
            else
              break;
            mid_id = start_id + ((end_id -start_id) / 2);
          }
          if (col_id[mid_id] == l_c_id)
            ma->coeff_recv_col_id[i] = mid_id;
          else {
            /* special case for diagonal */
            if (ma->separate_diag && l_c_id == l_r_id)
              ma->coeff_recv_col_id[i] = -1;
            /* unexpected case */
            else
              bft_error
                (__FILE__, __LINE__, 0,
                 "%s:\n"
                 "  Unexpected and unhandled use case:\n"
                 "    a distant rank has defined a matrix relation between\n"
                 "    elements %d and %d assigned to the current rank, with no\n"
                 "    such relation previously defined on the current rank.\n",
                 __func__, (int)l_r_id, (int)l_c_id);
          }
        }

        /* Distant part */

        else {

          const cs_gnum_t *g_col_id = ma->d_g_c_id + ma->d_r_idx[l_r_id];
          cs_lnum_t start_id = 0;
          cs_lnum_t end_id = ma->d_r_idx[l_r_id+1] - ma->d_r_idx[l_r_id] - 1;
          cs_lnum_t mid_id = (end_id -start_id) / 2;
          while (start_id <= end_id) {
            if (g_col_id[mid_id] < g_c_id)
              start_id = mid_id + 1;
            else if (g_col_id[mid_id] > g_c_id)
              end_id = mid_id - 1;
            else
              break;
            mid_id = start_id + ((end_id -start_id) / 2);
          }
          assert(g_col_id[mid_id] == g_c_id);
          /* column ids start and end of local row, so add r_idx[row_id+1] */
          ma->coeff_recv_col_id[i] = mid_id + ma->r_idx[l_r_id] + 1;
        }

      }

    }
    else { /* For most external libraries, use global column ids */

      BFT_MALLOC(ma->coeff_recv_col_g_id, ma->coeff_recv_size, cs_gnum_t);

      for (cs_lnum_t i = 0; i < ma->coeff_recv_size; i++) {
        ma->coeff_recv_row_id[i] = recv_data[i*2] - ma->l_range[0];
        ma->coeff_recv_col_g_id[i] = recv_data[i*2+1];
      }

    }

    BFT_FREE(recv_data);
  }

  BFT_FREE(request);
  BFT_FREE(status);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build vector halo structure associated with a matrix assembler.
 *
 * This assumes the ordered and compact list of external global ids
 * has been built.

 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_compute_halo(cs_matrix_assembler_t  *ma)
{
  /* Build neighbor rank info */

  cs_lnum_t   n_e_g_ids = 0;
  cs_gnum_t  *e_g_id = NULL;

  cs_lnum_t n = ma->d_r_idx[ma->n_rows];

  BFT_MALLOC(e_g_id, n, cs_gnum_t);

  /* Initial pass to reduce data a bit before sorting */

  cs_gnum_t g_id_prev = ma->l_range[0];

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_gnum_t g_id = ma->d_g_c_id[i];
    if (   (g_id < ma->l_range[0] || g_id >= ma->l_range[1])
        && g_id != g_id_prev) {
      e_g_id[n_e_g_ids++] = g_id;
      g_id_prev = g_id;
    }
  }

  n_e_g_ids = cs_sort_and_compact_gnum(n_e_g_ids, e_g_id);
  BFT_REALLOC(e_g_id, n_e_g_ids, cs_gnum_t);

  /* Now determine assumed rank neighbors based on compact info */

  cs_rank_neighbors_t *arn = _assumed_rank_neighbors(n_e_g_ids,
                                                     ma->n_g_rows,
                                                     ma->l_range,
                                                     e_g_id,
                                                     ma->comm);

  ma->n_ranks_init[1] = arn->size;

  /* Exchange info on distant ranks */

  cs_gnum_t *d_ranges = _assumed_rank_exchange(arn,
                                               ma->l_range,
                                               ma->comm);

  /* Identifiy rank and local id associated with each element */

  int        *rank_id;
  cs_lnum_t  *r_loc_id;
  BFT_MALLOC(rank_id, n_e_g_ids, int);
  BFT_MALLOC(r_loc_id, n_e_g_ids, cs_lnum_t);

# pragma omp parallel if(n_e_g_ids > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_e_g_ids; i++) {
    cs_gnum_t g_id = e_g_id[i];
    int r_id = _g_id_rank(arn->size, d_ranges, g_id);
    assert(r_id > -1);
    assert(g_id >= d_ranges[r_id*2]);
    r_loc_id[i] = g_id - d_ranges[r_id*2];
    rank_id[i] = r_id;
  }

  BFT_FREE(d_ranges);

  ma->halo = cs_halo_create_from_rank_neighbors(arn,
                                                 ma->n_rows,
                                                 n_e_g_ids,
                                                 rank_id,
                                                 r_loc_id);

  BFT_FREE(r_loc_id);
  BFT_FREE(rank_id);

  cs_rank_neighbors_destroy(&arn);
  BFT_FREE(e_g_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler
 *        in parallel mode.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_compute_mpi(cs_matrix_assembler_t  *ma)
{
  /* Compute number of rows */

  cs_lnum_t  n_rows = 0;
  if (ma->l_range[1] > ma->l_range[0])
    n_rows = ma->l_range[1] - ma->l_range[0];

  MPI_Allreduce(ma->l_range+1, &(ma->n_g_rows), 1, CS_MPI_GNUM,
                MPI_MAX, ma->comm);

  ma->n_rows = n_rows;

  /* Separate local fom distant rows: counting step */

  cs_lnum_t  e_size = 0;

  cs_lnum_t  *l_c_count, *d_c_count;
  BFT_MALLOC(ma->r_idx, n_rows+1, cs_lnum_t);
  BFT_MALLOC(ma->d_r_idx, n_rows+1, cs_lnum_t);
  BFT_MALLOC(l_c_count, n_rows, cs_lnum_t);
  BFT_MALLOC(d_c_count, n_rows, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    l_c_count[i] = 0;
    d_c_count[i] = 0;
  }

  for (cs_lnum_t i = 0; i < ma->size; i++) {
    cs_gnum_t g_r_id = ma->g_rc_id[i*2];
    cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
    if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
      cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
      if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1])
        l_c_count[l_r_id] += 1;
      else
        d_c_count[l_r_id] += 1;
    }
    else
      e_size += 1;
  }

  /* Build index and reset count */

  ma->r_idx[0] = 0;
  ma->d_r_idx[0] = 0;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    ma->r_idx[i+1] = ma->r_idx[i] + l_c_count[i];
    ma->d_r_idx[i+1] = ma->d_r_idx[i] + d_c_count[i];
    l_c_count[i] = 0;
    d_c_count[i] = 0;
  }

  /* Allocate structures holding data */

  cs_gnum_t *e_g_ij;
  BFT_MALLOC(ma->c_id, ma->r_idx[n_rows], cs_lnum_t);
  BFT_MALLOC(ma->d_g_c_id, ma->d_r_idx[n_rows], cs_gnum_t);
  BFT_MALLOC(e_g_ij, e_size*2, cs_gnum_t);

  /* Now fill data: local part is determined (will be cleaned),
     and list of distant rows defined; in cases where sizeof(cs_lnum_t),
     is smaller than sizeof(cs_gnum_t), this data, though
     unsorted, will already be smaller than the initial couples. */

  cs_lnum_t k = 0;

  for (cs_lnum_t i = 0; i < ma->size; i++) {
    cs_gnum_t g_r_id = ma->g_rc_id[i*2];
    cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
    if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
      cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
      if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {
        cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
        ma->c_id[ma->r_idx[l_r_id] + l_c_count[l_r_id]] = l_c_id;
        l_c_count[l_r_id] += 1;
      }
      else {
        ma->d_g_c_id[ma->d_r_idx[l_r_id] + d_c_count[l_r_id]] = g_c_id;
        d_c_count[l_r_id] += 1;
      }
    }
    else {
      e_g_ij[k++] = g_r_id;
      e_g_ij[k++] = g_c_id;
    }
  }

  BFT_FREE(d_c_count);
  BFT_FREE(l_c_count);

  /* We can free the initial insertion structure */

  BFT_FREE(ma->g_rc_id);

  assert(k == e_size*2);

  /* Rows and column ids in the local range can be finalized now;
     for most matrix structures, we will need to expand the
     rows with columns referencing distant rows (except when storing
     those terms separatly or passing them to an external library),
     but updating it now allows optimizing for space, and will make
     future updates simpler. */

  _sort_and_compact_local(ma);
  _sort_and_compact_distant(ma);

  /* We can also order and compact data that will be sent to external
     ranks now, as it will be needed at some point, and is simpler
     to do on global id couples than on local ids with higher stride.
     This may also reduce the size of the data we work on for the
     following steps. */

  e_size = cs_sort_and_compact_gnum_2(e_size, e_g_ij);

  _process_assembly_data(ma, e_size, e_g_ij);

  BFT_FREE(e_g_ij);

  _matrix_assembler_compute_halo(ma);
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler
 *        in serial mode.
 *
 * \param[in, out]  ma   pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_compute_local(cs_matrix_assembler_t  *ma)
{
  /* Separate local fom distant rows: counting step */

  cs_lnum_t n_rows = 0;
  if (ma->l_range[1] > ma->l_range[0])
    n_rows = ma->l_range[1] - ma->l_range[0];

  cs_lnum_t  *c_count;
  BFT_MALLOC(ma->r_idx, n_rows+1, cs_lnum_t);
  BFT_MALLOC(c_count, n_rows, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++)
    c_count[i] = 0;

  for (cs_lnum_t i = 0; i < ma->size; i++) {
    cs_gnum_t g_r_id = ma->g_rc_id[i*2];
    cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
    assert(g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]);
    cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
    assert(g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]);
    c_count[l_r_id] += 1;
  }

  /* Build index and reset count */

  ma->r_idx[0] = 0;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    ma->r_idx[i+1] = ma->r_idx[i] + c_count[i];
    c_count[i] = 0;
  }

  /* Allocate structures holding data */

  BFT_MALLOC(ma->c_id, ma->r_idx[n_rows], cs_lnum_t);

  /* Now fill data: local part is determined (will be cleaned),
     and list of distant rows defined */

  for (cs_lnum_t i = 0; i < ma->size; i++) {
    cs_lnum_t l_r_id = ma->g_rc_id[i*2] - ma->l_range[0];
    ma->c_id[ma->r_idx[l_r_id] + c_count[l_r_id]]
      = ma->g_rc_id[i*2+1] - ma->l_range[0];
    c_count[l_r_id] += 1;
  }

  /* Set global number of rows (will be updated in parallel) */

  ma->n_g_rows = n_rows;
  ma->n_rows = n_rows;

  /* Sort and compact elements for final structure */

  _sort_and_compact_local(ma);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix assembler structure.
 *
 * The associated matrix structure data will initially be empty, though
 * the range of rows associated with the current rank must be defined
 * immediately.
 *
 * \param[in]  l_range        global id range [min, max[ for local rank
 * \param[in]  separate_diag  if true, diagonal terms are handled separately
 *
 * \return  pointer to created matrix assembler structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create(cs_gnum_t  l_range[2],
                           bool       separate_diag)
{
  cs_matrix_assembler_t *ma = NULL;

  BFT_MALLOC(ma, 1, cs_matrix_assembler_t);

  ma->separate_diag = separate_diag;

  ma->external_lib_ids = false;

  ma->l_range[0] = l_range[0];
  ma->l_range[1] = l_range[1];

  ma->n_g_rows = 0;
  ma->n_rows = 0;

  ma->size = 0;
  ma->max_size = 0;

  ma->r_idx = NULL;
  ma->c_id = NULL;

  ma->g_rc_id = NULL;

#if defined(HAVE_MPI)

  ma->d_r_idx = NULL;
  ma->d_g_c_id = NULL;

  ma->n_coeff_ranks = 0;
  ma->coeff_rank = NULL;

  ma->coeff_send_size = 0;
  ma->coeff_recv_size = 0;

  ma->coeff_send_n_rows = 0;
  ma->coeff_send_index = NULL;
  ma->coeff_send_row_g_id = NULL;
  ma->coeff_send_col_g_id = NULL;

  ma->coeff_rank_send_index = NULL;
  ma->coeff_rank_recv_index = NULL;

  ma->coeff_recv_row_id = NULL;
  ma->coeff_recv_col_id = NULL;
  ma->coeff_recv_col_g_id = NULL;

  ma->halo = NULL;

  ma->comm = cs_glob_mpi_comm;

#endif /* HAVE_MPI */

  ma->n_ranks_init[0] = 0;
  ma->n_ranks_init[1] = 0;

  return ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix assembler structure.
 *
 * \param[in, out]  pointer to matrix assembler structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_destroy(cs_matrix_assembler_t  **ma)
{
  if (ma != NULL) {
    cs_matrix_assembler_t *_ma = *ma;

#if defined(HAVE_MPI)

    if (_ma->halo != NULL)
      cs_halo_destroy(&(_ma->halo));

    BFT_FREE(_ma->coeff_recv_col_g_id);
    BFT_FREE(_ma->coeff_recv_col_id);
    BFT_FREE(_ma->coeff_recv_row_id);

    BFT_FREE(_ma->coeff_rank_recv_index);
    BFT_FREE(_ma->coeff_rank_send_index);

    BFT_FREE(_ma->coeff_send_col_g_id);
    BFT_FREE(_ma->coeff_send_row_g_id);
    BFT_FREE(_ma->coeff_send_index);
    BFT_FREE(_ma->coeff_rank);

    BFT_FREE(_ma->d_g_c_id);
    BFT_FREE(_ma->d_r_idx);

#endif /* HAVE_MPI */

    BFT_FREE(_ma->g_rc_id);

    BFT_FREE(_ma->c_id);
    BFT_FREE(_ma->r_idx);

    BFT_FREE(*ma);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add entries to a matrix assembler structure.
 *
 * This function should be called by a single thread for a given assembler.
 *
 * \param[in, out]  ma        pointer to matrix assembler structure
 * \param[in]       n         number of entries
 * \param[in]       g_col_id  global column ids associated with entries
 * \param[in]       g_row_id  global row ids associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_add_ids(cs_matrix_assembler_t  *ma,
                            cs_lnum_t               n,
                            cs_gnum_t               g_row_id[],
                            cs_gnum_t               g_col_id[])
{
  /* Reallocate if needed;
     note that we could use optimized structures to avoid
     inserting duplicate values; A simple and relatively
     efficient variant may be to preallocate a local matrix
     with a fixed estimated number of columns, so as to
     handle the bulk of data, and add additional bins for
     the rest; sorting and removal of duplicates could occur
     when storage of preallocated sections is full.
     This may be done at a future stage, depending on timing
     and memory consumption measures for the current implementation. */

  if (ma->size + n >= ma->max_size) {
    if (ma->size == 0)
      ma->max_size = 4;
    while (ma->size + n >= ma->max_size)
      ma->max_size *= 2;
    BFT_REALLOC(ma->g_rc_id, ma->max_size*2, cs_gnum_t);
  }

  cs_gnum_t *_g_rc_id = ma->g_rc_id + ma->size*2;

  if (ma->separate_diag == false) {
    for (cs_lnum_t i = 0; i < n; i++) {
      _g_rc_id[i*2]   = g_row_id[i];
      _g_rc_id[i*2+1] = g_col_id[i];
    }
    ma->size += n;
  }
  else {
    for (cs_lnum_t i = 0; i < n; i++) {
      if (   g_row_id[i] != g_col_id[i]
          || g_row_id[i] <  ma->l_range[0]
          || g_row_id[i] >= ma->l_range[1]) {
        _g_rc_id[i*2]   = g_row_id[i];
        _g_rc_id[i*2+1] = g_col_id[i];
        ma->size += 1;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler.
 *
  * The associated vector halo is also computed at this stage.
 *
* \param[in, out]  ma    pointer to matrix assembler structure
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_compute_structure(cs_matrix_assembler_t  *ma)
{
#if defined(HAVE_MPI)
  if (ma->comm != MPI_COMM_NULL && ma->comm != MPI_COMM_SELF) {
    _matrix_assembler_compute_mpi(ma);
    return;
  }
#endif /* HAVE_MPI */

  _matrix_assembler_compute_local(ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign specific MPI communicator to matrix assembler.
 *
 * This must be called before \ref cs_matrix_assembler_compute.
 *
 * \param[in, out]  ma    pointer to matrix assembler structure
 * \param[in]       comm  assigned MPI communicator
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_matrix_assembler_set_comm(cs_matrix_assembler_t  *ma,
                             MPI_Comm                comm)
{
  assert(ma->n_ranks_init[1] == 0);

  ma->comm = comm;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to the halo structure associated with a
 *        matrix assembler.
 *
 * \param[in]  ma   pointer to matrix assembler structure
 *
 * \return  pointer to halo structure
 */
/*----------------------------------------------------------------------------*/

const cs_halo_t *
cs_matrix_assembler_get_halo(const cs_matrix_assembler_t  *ma)
{
#if defined(HAVE_MPI)

  return ma->halo;

#else

  return NULL;

#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
