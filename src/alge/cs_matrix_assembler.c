/*============================================================================
 * Incremental or general construction of matrix.
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

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_rank_neighbors.h"
#include "cs_sort.h"
#include "cs_timer.h"

#include "cs_matrix.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix_assembler_priv.h"
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

/* Fixed coefficient buffer size for accumulation
   (a reasonably small fixed size has the advantage of being easily usable
   on the stack and in a threading context, and that size should still
   be large enough to amortize calls to lower-level functions */

#undef  COEFF_GROUP_SIZE
#define COEFF_GROUP_SIZE 256

#undef  HISTOGRAM_SUBS
#define HISTOGRAM_SUBS  5

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print distribution of counter over ranks info to a given log type.
 *
 * \param[in]  log_id  log file type
 * \param[in]  count   counter value for current rank
 */
/*----------------------------------------------------------------------------*/

static void
_display_rank_histogram(cs_log_t  log_id,
                        int       count)
{
  int  i, j, k, count_max, count_min;
  double step;

  int h_count[HISTOGRAM_SUBS];
  int n_steps = HISTOGRAM_SUBS;

  const int n_counts = cs_glob_n_ranks;

  int *r_count;
  BFT_MALLOC(r_count, cs_glob_n_ranks, int);

  r_count[0] = count;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    MPI_Allgather(&count, 1, MPI_INT,
                  r_count, 1, MPI_INT, cs_glob_mpi_comm);
#endif

  /* Compute local min and max */

  count_max = r_count[0];
  count_min = r_count[0];
  for (i = 1; i < n_counts; i++) {
    count_min = CS_MIN(count_min, r_count[i]);
    count_max = CS_MAX(count_max, r_count[i]);
  }

  cs_log_printf(log_id, _("    minimum count =         %10d\n"), count_min);
  cs_log_printf(log_id, _("    maximum count =         %10d\n\n"), count_max);

  /* Define axis subdivisions */

  for (j = 0; j < n_steps; j++)
    h_count[j] = 0;

  if (count_max - count_min > 0) {

    if (count_max-count_min < n_steps)
      n_steps = CS_MAX(1, floor(count_max-count_min));

    step = (double)(count_max - count_min) / n_steps;

    /* Loop on counts */

    for (i = 0; i < n_counts; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (r_count[i] < count_min + k*step)
          break;
      }
      h_count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      cs_log_printf(log_id, "    %3d : [ %10d ; %10d [ = %10d\n",
                    i+1,
                    (int)(count_min + i*step),
                    (int)(count_min + j*step),
                    h_count[i]);

    cs_log_printf(log_id, "    %3d : [ %10d ; %10d ] = %10d\n",
                  n_steps,
                  (int)(count_min + (n_steps - 1)*step),
                  count_max,
                  h_count[n_steps - 1]);

  }

  else { /* if (count_max == count_min) */
    cs_log_printf(log_id, "    %3d : [ %10d ; %10d ] = %10d\n",
                  1, count_min, count_max, n_counts);
  }

  BFT_FREE(r_count);
}

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
    cs_lnum_t *col_id = ma->_c_id + ma->_r_idx[i];
    cs_lnum_t n_cols = ma->_r_idx[i+1] - ma->_r_idx[i];
    for (cs_lnum_t j = 1; j < n_cols; j++) {
      if (col_id[j] <= col_id[j-1])
        ordered = false;
    }
  }

  if (ordered)
    return;

  /* Sort by row */

  bool direct_assembly = cs_sort_indexed(n_rows, ma->_r_idx, ma->_c_id);

  /* Compact elements if necessary */

  if (direct_assembly == false) {

    cs_lnum_t *tmpr_idx = NULL;

    BFT_MALLOC(tmpr_idx, n_rows+1, cs_lnum_t);
    memcpy(tmpr_idx, ma->_r_idx, (n_rows+1)*sizeof(cs_lnum_t));

    cs_lnum_t  k = 0;

    for (cs_lnum_t i = 0; i < n_rows; i++) {
      cs_lnum_t *col_id = ma->_c_id + ma->_r_idx[i];
      cs_lnum_t n_cols = ma->_r_idx[i+1] - ma->_r_idx[i];
      cs_lnum_t col_id_prev = -1;
      ma->_r_idx[i] = k;
      for (cs_lnum_t j = 0; j < n_cols; j++) {
        if (col_id_prev != col_id[j]) {
          ma->_c_id[k++] = col_id[j];
          col_id_prev = col_id[j];
        }
      }
    }
    ma->_r_idx[n_rows] = k;

    assert(ma->_r_idx[n_rows] < tmpr_idx[n_rows]);

    BFT_FREE(tmpr_idx);
    BFT_REALLOC(ma->_c_id, (ma->_r_idx[n_rows]), cs_lnum_t);
    ma->c_id = ma->_c_id;

  }
}

#if defined(HAVE_MPI)

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
_g_id_rank_index(int              n_ranges,
                 const cs_gnum_t  d_ranges[],
                 cs_gnum_t        g_id)
{
  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = n_ranges - 1;
  cs_lnum_t mid_id = (end_id -start_id) / 2;

  while (start_id < end_id) {
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
    cs_gnum_t *col_id = ma->d_g_c_id + ma->d_r_idx[i];
    cs_lnum_t n_cols = ma->d_r_idx[i+1] - ma->d_r_idx[i];
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
      cs_gnum_t *col_g_id = ma->d_g_c_id + ma->d_r_idx[i];
      cs_lnum_t n_cols = ma->d_r_idx[i+1] - ma->d_r_idx[i];
      ma->d_r_idx[i] = k;
      if (n_cols > 0)
        ma->d_g_c_id[k++] = col_g_id[0];
      for (cs_lnum_t j = 1; j < n_cols; j++) {
        if (col_g_id[j] != col_g_id[j-1])
          ma->d_g_c_id[k++] = col_g_id[j];
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
 * \brief Add terms to local matrix elements defined by distant ranks
 *
 * \param[in, out]  ma    pointer to matrix assembler structure
 * \param[in]       n     number of added column and row couples
 * \param[in]       l_ij  added couples
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_complete_local(cs_matrix_assembler_t  *ma,
                cs_lnum_t               n,
                cs_lnum_2_t             l_ij[])
{
  /* Count maximum terms to add for each row */

  cs_lnum_t n_rows = ma->n_rows;

  cs_lnum_t *l_c_count, *l_r_idx;
  BFT_MALLOC(l_c_count, n_rows, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++)
    l_c_count[i] = 0;

  for (cs_lnum_t i = 0; i < n; i++)
    l_c_count[l_ij[i][0]] += 1;

  BFT_MALLOC(l_r_idx, n_rows+1, cs_lnum_t);

  l_r_idx[0] = 0;
  for (cs_lnum_t i = 0; i < n_rows; i++)
    l_r_idx[i+1] = l_r_idx[i] + l_c_count[i];

  /* Expand matrix, starting copies from end to avoid overwrites
     (first line untouched) */

  BFT_REALLOC(ma->_c_id, ma->r_idx[n_rows] + l_r_idx[n_rows], cs_lnum_t);
  ma->c_id = ma->_c_id;

  for (cs_lnum_t i = n_rows-1; i > 0; i--) {
    cs_lnum_t n_cols = ma->_r_idx[i+1] - ma->_r_idx[i];
    cs_lnum_t *col_id = ma->_c_id + ma->_r_idx[i] + l_r_idx[i];
    const cs_lnum_t *col_id_s = ma->_c_id + ma->_r_idx[i];
    l_c_count[i] = n_cols;
    ma->_r_idx[i+1] += l_r_idx[i+1];
    for (cs_lnum_t j = n_cols-1; j >= 0; j--)
      col_id[j] = col_id_s[j];
  }
  l_c_count[0] = ma->_r_idx[1];
  ma->_r_idx[1] += l_r_idx[1];

  BFT_FREE(l_r_idx);

  /* Now add terms */

  for (cs_lnum_t i = 0; i < n; i++) {
    cs_lnum_t r_id = l_ij[i][0];
    cs_gnum_t c_id = l_ij[i][1];
    ma->_c_id[ma->_r_idx[r_id] + l_c_count[r_id]] = c_id;
    l_c_count[r_id] += 1;
  }
  BFT_FREE(l_c_count);

  /* Sort and remove duplicates */

  _sort_and_compact_local(ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add terms to distant matrix elements.
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
    cs_gnum_t *col_g_id_d = ma->d_g_c_id + ma->d_r_idx[i] + d_r_idx[i];
    const cs_gnum_t *col_g_id_s = ma->d_g_c_id + ma->d_r_idx[i];
    d_c_count[i] = n_cols;
    ma->d_r_idx[i+1] += d_r_idx[i+1];
    for (cs_lnum_t j = n_cols-1; j >= 0; j--)
      col_g_id_d[j] = col_g_id_s[j];
  }
  if (n_rows > 0) {
    d_c_count[0] = ma->d_r_idx[1];
    ma->d_r_idx[1] += d_r_idx[1];
  }

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Append distant row ids to local row ids.
 *
 * As column ids remain ordered in each row, local columns appear first,
 * and ids matching distant columns second.
 *
 * The resulting structure is expected to be best suited for matrix-vector
 * products, or row extractions, as we do not need to handle separate
 * matrices, and multiple loops.
 *
 * For assembly, the extended part is also maintained separately. For a given
 * column appearing on row i at extended position j (i.e. ma->d_r_idx[i] + j),
 * its matching position in the extended matrix will be:
 * ma->r_idx[i+1] - (ma->d_r_idx[i+1] - ma->d_r_idx[i] + j.
 *
 * \param[in, out]  ma         pointer to matrix assembler structure
 * \param[in]       n_e_g_ids  number of unique external global ids
 * \param[in]       e_g_id     ordered unique external global ids
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_append_local_and_distant(cs_matrix_assembler_t  *ma,
                          cs_lnum_t               n_e_g_ids,
                          const cs_gnum_t         e_g_id[])
{
  /* Count terms to add for each row */

  cs_lnum_t n_rows = ma->n_rows;

  cs_lnum_t *c_count;
  BFT_MALLOC(c_count, n_rows, cs_lnum_t);

  /* Expand matrix, starting copies from end to avoid overwrites
     (first line untouched) */

  BFT_REALLOC(ma->_c_id, ma->_r_idx[n_rows] + ma->d_r_idx[n_rows], cs_lnum_t);
  ma->c_id = ma->_c_id;

  for (cs_lnum_t i = n_rows-1; i > 0; i--) {
    cs_lnum_t n_cols = ma->_r_idx[i+1] - ma->_r_idx[i];
    cs_lnum_t *col_id_d = ma->_c_id + ma->_r_idx[i] + ma->d_r_idx[i];
    const cs_lnum_t *col_id_s = ma->_c_id + ma->_r_idx[i];
    c_count[i] = n_cols;
    ma->_r_idx[i+1] += ma->d_r_idx[i+1];
    for (cs_lnum_t j = n_cols-1; j >= 0; j--)
      col_id_d[j] = col_id_s[j];
  }
  if (n_rows > 0) {
    c_count[0] = ma->_r_idx[1];
    ma->_r_idx[1] += ma->d_r_idx[1];
  }

  /* Now insert terms */

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    cs_lnum_t n_d_cols = ma->d_r_idx[i+1] - ma->d_r_idx[i];
    cs_lnum_t *col_id_d = ma->_c_id + ma->_r_idx[i+1] - n_d_cols;
    const cs_gnum_t *col_g_id_s = ma->d_g_c_id + ma->d_r_idx[i];
    for (cs_lnum_t j = 0; j < n_d_cols; j++) {
      cs_gnum_t g_c_id = col_g_id_s[j];
      cs_lnum_t k = _g_id_binary_find(n_e_g_ids, g_c_id, e_g_id);
      col_id_d[j] = n_rows + k;
    }
  }

  BFT_FREE(c_count);
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
 * \brief Exchange range info for building assumed rank neighbors info
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
_rank_ranges_exchange(cs_rank_neighbors_t  *arn,
                      cs_gnum_t             l_range[2],
                      MPI_Comm              comm)
{
  MPI_Request *request = NULL;
  MPI_Status *status = NULL;
  cs_gnum_t *d_ranges = NULL;

  BFT_MALLOC(d_ranges, arn->size*2, cs_gnum_t);
  BFT_MALLOC(request, arn->size*2, MPI_Request);
  BFT_MALLOC(status, arn->size*2, MPI_Status);

  /* Prepare for determination of assumed rank */

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
 * \brief Build true rank info based on assumed rank info
 *
 * This operation is collective on communicator comm.
 *
 * The caller is responsible for freeing the returned array
 *
 * \param[in]  arn       info on ranks with which we initially communicate
 * \param[in]  n         number of global ids for local rank
 * \param[in]  n_g       global number of ids
 * \param[in]  l_range   global id range [min, max[ for local rank
 * \param[in]  d_ranges  global id ranges [min, max[ for neighboring ranks
 * \param[in]  g_id      sorted array of global ids for local rank
 *                       (size: n)
 * \param[in]  comm      associated communicator
 *
 * \return  effective ranks matching given global ids (size n)
 */
/*----------------------------------------------------------------------------*/

static int *
_assumed_to_true_rank(cs_rank_neighbors_t  *arn,
                      cs_lnum_t             n,
                      cs_gnum_t             n_g,
                      const cs_gnum_t       l_range[2],
                      const cs_gnum_t       d_ranges[],
                      const cs_gnum_t       g_id[],
                      MPI_Comm              comm)
{
  int *d_rank;
  BFT_MALLOC(d_rank, n, int);

  MPI_Request *request = NULL;
  MPI_Status *status = NULL;

  const int local_rank = cs_glob_rank_id;

  /* Prepare for determination of assumed rank;
     we will send both a partial local range description and global ids to the
     assumed partition, using a common array, so as to minimize
     all-to-all communication. */

  int n_ranks, l_rank;
  int *a_rank;

  MPI_Comm_size(comm, &n_ranks);
  MPI_Comm_rank(comm, &l_rank);

  BFT_MALLOC(request, arn->size*2, MPI_Request);
  BFT_MALLOC(status, arn->size*2, MPI_Status);

  /* Determine ranks with which we will first exchange;
     we filter and immediately handle elements whose assumed
     rank matches the current rank */

  BFT_MALLOC(a_rank, n, int);

  cs_lnum_t n_range_match = 0;
  cs_lnum_t n_range_dist = 0;

  int a_rank_prev = -1;
  for (cs_lnum_t i = 0; i < n; i++) {
    int r = _assumed_rank(n_ranks, n_g, g_id[i]);
    if (r == local_rank) {
      int r_idx = _g_id_rank_index(arn->size, d_ranges, g_id[i]);
      assert(r_idx >= 0);
      d_rank[i] = arn->rank[r_idx];
      n_range_match += 1;
    }
    else {
      d_rank[i] = -1; /* flag for later */
      a_rank[n_range_dist] = r;
      n_range_dist += 1;
    }
    assert(r >= a_rank_prev);
    a_rank_prev = r;
  }

  cs_rank_neighbors_to_index(arn, n_range_dist, a_rank, a_rank);

  /* Count number of exchanges with each assumed rank
     and build index, first handled as counts;
     local values have index -1 */

  cs_lnum_t *send_index;
  BFT_MALLOC(send_index, arn->size+1, cs_lnum_t);

  send_index[0] = 0;
  cs_rank_neighbors_count(arn, n_range_dist, a_rank, send_index + 1);

  BFT_FREE(a_rank);

  /* Exchange counts */

  cs_lnum_t *recv_index;
  BFT_MALLOC(recv_index, arn->size + 1, cs_lnum_t);

  recv_index[0] = 0;

  int request_count = 0;

  for (int i = 0; i < arn->size; i++) {
    MPI_Irecv(recv_index + i + 1,
              1,
              CS_MPI_LNUM,
              arn->rank[i],
              local_rank,
              comm,
              &(request[request_count++]));
  }

  for (int i = 0; i < arn->size; i++) {
    MPI_Isend(send_index + i + 1,
              1,
              CS_MPI_LNUM,
              arn->rank[i],
              arn->rank[i],
              comm,
              &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  /* Now transform counts to indexes */

  for (int i = 0; i < arn->size; i++)
    send_index[i+1] += send_index[i];

  for (int i = 0; i < arn->size; i++)
    recv_index[i+1] += recv_index[i];

  /* Now exchange global ids */

  const cs_lnum_t recv_size = recv_index[arn->size];

  cs_gnum_t *recv_g_id;
  BFT_MALLOC(recv_g_id, recv_size, cs_gnum_t);

  request_count = 0;

  for (int i = 0; i < arn->size; i++) {
    if (recv_index[i+1] > recv_index[i])
      MPI_Irecv(recv_g_id + recv_index[i],
                recv_index[i+1] - recv_index[i],
                CS_MPI_GNUM,
                arn->rank[i],
                local_rank,
                comm,
                &(request[request_count++]));
  }

  for (int i = 0; i < arn->size; i++) {
    if (send_index[i+1] > send_index[i]) {
      cs_lnum_t shift = (arn->rank[i] > local_rank) ? n_range_match : 0;
      MPI_Isend(g_id + send_index[i] + shift,
                send_index[i+1] - send_index[i],
                CS_MPI_GNUM,
                arn->rank[i],
                arn->rank[i],
                comm,
                &(request[request_count++]));
    }
  }

  MPI_Waitall(request_count, request, status);

  /* Now for each global id, determine its rank */

  int *d_rank_ar;
  BFT_MALLOC(d_rank_ar, recv_size, int);

  for (cs_lnum_t i = 0; i < recv_size; i++) {
    cs_gnum_t _g_id = recv_g_id[i];
    if (l_range[0] <= _g_id && l_range[1] > _g_id)
      d_rank_ar[i] = local_rank;
    else {
      int r_idx = _g_id_rank_index(arn->size, d_ranges, recv_g_id[i]);
      assert(r_idx >= 0);
      d_rank_ar[i] = arn->rank[r_idx];
      assert(d_rank_ar[i] >= 0);
    }
  }

  BFT_FREE(recv_g_id);

  /* Send data back to original rank;
     remember part of the data is already known */

  request_count = 0;

  for (int i = 0; i < arn->size; i++) {
    if (send_index[i+1] > send_index[i]) {
      cs_lnum_t shift = (arn->rank[i] > local_rank) ? n_range_match : 0;
      MPI_Irecv(d_rank + send_index[i] + shift,
                send_index[i+1] - send_index[i],
                MPI_INT,
                arn->rank[i],
                arn->rank[i],
                comm,
                &(request[request_count++]));
    }
  }

  for (int i = 0; i < arn->size; i++) {
    if (recv_index[i+1] > recv_index[i])
      MPI_Isend(d_rank_ar + recv_index[i],
                recv_index[i+1] - recv_index[i],
                MPI_INT,
                arn->rank[i],
                local_rank,
                comm,
                &(request[request_count++]));
  }

  MPI_Waitall(request_count, request, status);

  /* Note: we could add an additional exchange here, for each
     connected rank, to indicate which ranks may need
     to communicate with it; this would avoid the need for a
     symmetrization of the rank neighbors structure built from
     the d_rank[] array, as ranks would also know who needs to
     send some info. An all-to-all communication with just
     1 value per rank could thus be replaced by a point-to-point
     communication, but with slightly larger messages.
     We are not sure the added complexity would be worth it,
     so at this stage, we just add this placeholder. */

  BFT_FREE(d_rank_ar);

  BFT_FREE(recv_index);
  BFT_FREE(send_index);

  BFT_FREE(request);
  BFT_FREE(status);

  return d_rank;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine rank neighbors based on global id info
 *
 * This assumes the ordered and compact list of external global ids
 * has been built.
 *
 * \param[in]   n_e_g_ids     number of unique external global ids
 * \param[in]   n_g           global number of ids
 * \param[in]   l_range       global id range [min, max[ for local rank
 * \param[in]   e_g_id        ordered unique external global ids
 * \param[out]  n_init_ranks  number of ranks in initial assumed
 *                            neighborhood (for logging), or NULL
 * \param[in]   comm          associated communicator
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static cs_rank_neighbors_t *
_rank_neighbors(cs_lnum_t          n_e_g_ids,
                cs_gnum_t          n_g,
                cs_gnum_t          l_range[2],
                const cs_gnum_t    e_g_id[],
                int               *n_init_ranks,
                MPI_Comm           comm)
{
  /* Determine assumed rank neighbors based on compact info */

  cs_rank_neighbors_t *arn = _assumed_rank_neighbors(n_e_g_ids,
                                                     n_g,
                                                     l_range,
                                                     e_g_id,
                                                     comm);

  if (n_init_ranks != NULL)
    *n_init_ranks = arn->size;

  cs_gnum_t *d_ranges = _rank_ranges_exchange(arn,
                                              l_range,
                                              comm);

  int *e_rank_id = _assumed_to_true_rank(arn,
                                         n_e_g_ids,
                                         n_g,
                                         l_range,
                                         d_ranges,
                                         e_g_id,
                                         comm);

  BFT_FREE(d_ranges);

  /* Now rebuild rank neighbors */

  cs_rank_neighbors_destroy(&arn);

  arn = cs_rank_neighbors_create(n_e_g_ids, e_rank_id);

  BFT_FREE(e_rank_id);

  cs_rank_neighbors_symmetrize(arn, comm);

  return arn;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Process and exchange info relative to data that will be sent to
 *        neighboring ranks.
 *
 * At this stage, ma->rank will contain the index of the rank in the rank
 * neighborhood structure, not the true MPI rank.
 *
 * \param[in, out]  ma           pointer to matrix assembler structure
 * \param[in]       e_g_ij_size  size of coefficient data to send
 * \param[in]       e_g_ij       coefficient data (row_g_id, col_g_id couples)
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

  cs_lnum_t  n_e_g_ids = 0;
  cs_gnum_t *e_g_id = NULL;

  for (cs_lnum_t i = 0; i < ma->coeff_send_size; i++) {
    cs_gnum_t g_r_id = e_g_ij[i*2];
    if (g_r_id != g_r_id_prev) {
      n_e_g_ids += 1;
      ma->coeff_send_n_rows += 1;
      g_r_id_prev = g_r_id;
    }
    cs_gnum_t g_c_id = e_g_ij[i*2+1];
    if (g_c_id < ma->l_range[0] || g_c_id >= ma->l_range[1])
      n_e_g_ids += 1;
  }

  if (ma->coeff_send_size > 0) {

    BFT_MALLOC(e_g_id, n_e_g_ids, cs_gnum_t);

    n_e_g_ids = 0;

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
        e_g_id[n_e_g_ids++] = g_r_id;
        row_count++;
        g_r_id_prev = g_r_id;
      }
      cs_gnum_t g_c_id = e_g_ij[i*2+1];
      ma->coeff_send_col_g_id[i] = g_c_id;
      if (g_c_id < ma->l_range[0] || g_c_id >= ma->l_range[1]) {
        e_g_id[n_e_g_ids++] = g_c_id;
      }
    }
    ma->coeff_send_index[row_count] = ma->coeff_send_size;

    n_e_g_ids = cs_sort_and_compact_gnum(n_e_g_ids, e_g_id);

  }

  /* Determine ranks we may be communicating with for this stage */

  cs_rank_neighbors_t *arn = _rank_neighbors(n_e_g_ids,
                                             ma->n_g_rows,
                                             ma->l_range,
                                             e_g_id,
                                             &(ma->n_ranks_init[0]),
                                             ma->comm);

  BFT_FREE(e_g_id);

  cs_gnum_t *d_ranges = _rank_ranges_exchange(arn, ma->l_range, ma->comm);

  /* Prepare counts */

  cs_lnum_t *counts;
  BFT_MALLOC(counts, arn->size*2, cs_lnum_t);

  for (int i = 0; i < arn->size*2; i++)
    counts[i] = 0;

  for (cs_lnum_t i = 0; i < ma->coeff_send_n_rows; i++) {
    cs_gnum_t g_r_id = ma->coeff_send_row_g_id[i];
    int r_rank_id = _g_id_rank_index(arn->size, d_ranges, g_r_id);
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

    /* Now we build the mapping from received terms to their positions in the
       array. For indexing, we implicitly consider that each row is built of 2
       sub-rows, with the local columns first, and distant columns second. */

    BFT_MALLOC(ma->coeff_recv_row_id, ma->coeff_recv_size, cs_lnum_t);

    if (ma->flags & CS_MATRIX_DISTANT_ROW_USE_COL_IDX) {

      cs_lnum_t    n_l_insert = 0, n_l_insert_max = 0;
      cs_lnum_2_t *l_insert = NULL;

      BFT_MALLOC(ma->coeff_recv_col_idx, ma->coeff_recv_size, cs_lnum_t);

      /* First pass: determine local row ids, and check if insertion of
         previously unknown entries is required; column ids are set for local
         entries */

      for (cs_lnum_t i = 0; i < ma->coeff_recv_size; i++) {

        cs_lnum_t l_r_id = recv_data[i*2] - ma->l_range[0];
        cs_gnum_t g_c_id = recv_data[i*2+1];

        ma->coeff_recv_row_id[i] = l_r_id;

        /* Local part */

        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {

          cs_lnum_t n_cols = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];
          cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
          const cs_lnum_t *col_id = ma->c_id + ma->r_idx[l_r_id];

          ma->coeff_recv_col_idx[i] = _l_id_binary_search(n_cols,
                                                          l_c_id,
                                                          col_id);

          /* Special case for separate diagonal, other cases require insertion */

          if (   ma->coeff_recv_col_idx[i] < 0
              && (!ma->separate_diag || l_c_id != l_r_id)) {

            if (n_l_insert >= n_l_insert_max) {
              n_l_insert_max = CS_MAX(n_l_insert_max*2, 16);
              BFT_REALLOC(l_insert, n_l_insert_max, cs_lnum_2_t);
            }

            l_insert[n_l_insert][0] = l_r_id;
            l_insert[n_l_insert][1] = l_c_id;
            n_l_insert++;

          }

        }

      }

      /* Insert additional local terms if required */

      if (n_l_insert > 0) {

        _complete_local(ma, n_l_insert, l_insert);

        BFT_FREE(l_insert);
        n_l_insert_max = 0;

      }

      /* Second pass: determine local column ids, and check if insertion of
         previously unknown entries is required */

      for (cs_lnum_t i = 0; i < ma->coeff_recv_size; i++) {

        cs_lnum_t l_r_id = recv_data[i*2] - ma->l_range[0];
        cs_gnum_t g_c_id = recv_data[i*2+1];

        /* Local part */

        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {

          if (n_l_insert == 0)  /* Already up-to-date if no insertion */
            continue;

          cs_lnum_t n_cols = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];
          cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
          const cs_lnum_t *col_id = ma->c_id + ma->r_idx[l_r_id];

          ma->coeff_recv_col_idx[i] = _l_id_binary_search(n_cols,
                                                          l_c_id,
                                                          col_id);

          assert(   ma->coeff_recv_col_idx[i] >= 0
                 || (ma->separate_diag && l_c_id == l_r_id));

        }

        /* Distant part */

        else {

          cs_lnum_t n_d_cols = ma->d_r_idx[l_r_id+1] - ma->d_r_idx[l_r_id];
          cs_lnum_t n_cols = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];
          const cs_gnum_t *col_g_id = ma->d_g_c_id + ma->d_r_idx[l_r_id];

          cs_lnum_t d_c_idx = _g_id_binary_find(n_d_cols,
                                                g_c_id,
                                                col_g_id);

          /* column ids start and end of local row, so add n_cols
             (local part only, matrices are not appened at this stage) */
          ma->coeff_recv_col_idx[i] = d_c_idx + n_cols;

        }

      }

    }

    if (ma->flags & CS_MATRIX_DISTANT_ROW_USE_COL_G_ID) {

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

 * \param[in, out]  ma         pointer to matrix assembler structure
 * \param[in]       n_e_g_ids  number of unique external global ids
 * \param[in]       e_g_id     ordered unique external global ids
 *
 * This function should be called by a single thread for a given assembler.
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_compute_halo(cs_matrix_assembler_t  *ma,
                               cs_lnum_t               n_e_g_ids,
                               const cs_gnum_t         e_g_id[])
{
  /* Determine assumed rank neighbors based on compact info */

  cs_rank_neighbors_t *arn = _rank_neighbors(n_e_g_ids,
                                             ma->n_g_rows,
                                             ma->l_range,
                                             e_g_id,
                                             &(ma->n_ranks_init[1]),
                                             ma->comm);

  cs_gnum_t *d_ranges = _rank_ranges_exchange(arn, ma->l_range, ma->comm);

  /* Identifiy rank and local id associated with each element */

  int        *rank_id;
  cs_lnum_t  *r_loc_id;
  BFT_MALLOC(rank_id, n_e_g_ids, int);
  BFT_MALLOC(r_loc_id, n_e_g_ids, cs_lnum_t);

# pragma omp parallel if(n_e_g_ids > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_e_g_ids; i++) {
    cs_gnum_t g_id = e_g_id[i];
    int r_id = _g_id_rank_index(arn->size, d_ranges, g_id);
    assert(r_id > -1);
    assert(g_id >= d_ranges[r_id*2]);
    r_loc_id[i] = g_id - d_ranges[r_id*2];
    rank_id[i] = r_id;
  }

  BFT_FREE(d_ranges);

  ma->_halo = cs_halo_create_from_rank_neighbors(arn,
                                                 ma->n_rows,
                                                 n_e_g_ids,
                                                 rank_id,
                                                 r_loc_id);
  ma->halo = ma->_halo;

  BFT_FREE(r_loc_id);
  BFT_FREE(rank_id);

  cs_rank_neighbors_destroy(&arn);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute internal structures required by a matrix assembler in
 *        parallel mode.
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

  BFT_MALLOC(ma->_r_idx, n_rows+1, cs_lnum_t);
  BFT_MALLOC(ma->d_r_idx, n_rows+1, cs_lnum_t);
  ma->r_idx = ma->_r_idx;

  BFT_MALLOC(l_c_count, n_rows, cs_lnum_t);
  BFT_MALLOC(d_c_count, n_rows, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    l_c_count[i] = 0;
    d_c_count[i] = 0;
  }

  if (ma->separate_diag) {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
        cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {
          if (g_c_id != g_r_id)
            l_c_count[l_r_id] += 1;
        }
        else
          d_c_count[l_r_id] += 1;
      }
      else
        e_size += 1;
    }

  }
  else {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
        cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {
          l_c_count[l_r_id] += 1;
        }
        else
          d_c_count[l_r_id] += 1;
      }
      else
        e_size += 1;
    }

  }

  /* Build index and reset count */

  ma->_r_idx[0] = 0;
  ma->d_r_idx[0] = 0;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    ma->_r_idx[i+1] = ma->_r_idx[i] + l_c_count[i];
    ma->d_r_idx[i+1] = ma->d_r_idx[i] + d_c_count[i];
    l_c_count[i] = 0;
    d_c_count[i] = 0;
  }

  /* Allocate structures holding data */

  BFT_MALLOC(ma->_c_id, ma->_r_idx[n_rows], cs_lnum_t);
  BFT_MALLOC(ma->d_g_c_id, ma->d_r_idx[n_rows], cs_gnum_t);
  ma->c_id = ma->_c_id;

  cs_gnum_t *e_g_ij;
  BFT_MALLOC(e_g_ij, e_size*2, cs_gnum_t);

  /* Now fill data: local part is determined (will be cleaned),
     and list of distant rows defined; in cases where sizeof(cs_lnum_t),
     is smaller than sizeof(cs_gnum_t), this data, though
     unsorted, will already be smaller than the initial couples. */

  cs_lnum_t k = 0;

  if (ma->separate_diag) {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
        cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
        if (   g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]
            && g_c_id != g_r_id) {
          cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
          ma->_c_id[ma->_r_idx[l_r_id] + l_c_count[l_r_id]] = l_c_id;
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

  }
  else {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      if (g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]) {
        cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
        if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {
          cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
          ma->_c_id[ma->_r_idx[l_r_id] + l_c_count[l_r_id]] = l_c_id;
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

  /* Build compact external global id array */

  cs_lnum_t   n_e_g_ids = 0;
  cs_gnum_t  *e_g_id = NULL;

  cs_lnum_t d_size = ma->d_r_idx[ma->n_rows];

  BFT_MALLOC(e_g_id, d_size, cs_gnum_t);

  /* Initial pass to reduce data a bit before sorting */

  cs_gnum_t g_id_prev = ma->l_range[0];

  for (cs_lnum_t i = 0; i < d_size; i++) {
    cs_gnum_t g_id = ma->d_g_c_id[i];
    if (   (g_id < ma->l_range[0] || g_id >= ma->l_range[1])
        && g_id != g_id_prev) {
      e_g_id[n_e_g_ids++] = g_id;
      g_id_prev = g_id;
    }
  }

  n_e_g_ids = cs_sort_and_compact_gnum(n_e_g_ids, e_g_id);
  BFT_REALLOC(e_g_id, n_e_g_ids, cs_gnum_t);

  /* Finally, append distant row to local rows, and build matching halo */

  _append_local_and_distant(ma, n_e_g_ids, e_g_id);

  if (! (ma->flags & CS_MATRIX_EXTERNAL_HALO)) {
    _matrix_assembler_compute_halo(ma, n_e_g_ids, e_g_id);
    assert(ma->halo->n_elts[0] == n_e_g_ids);
  }

  /* Assign array to structure */

  ma->n_e_g_ids = n_e_g_ids;
  ma->e_g_id = e_g_id; e_g_id = NULL;
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

#if defined(HAVE_CUDA)
  cs_alloc_mode_t alloc_mode = CS_ALLOC_HOST_DEVICE_SHARED;
#else
  cs_alloc_mode_t alloc_mode = CS_ALLOC_HOST;
#endif

  CS_MALLOC_HD(ma->_r_idx, n_rows+1, cs_lnum_t, alloc_mode);

  ma->r_idx = ma->_r_idx;

  BFT_MALLOC(c_count, n_rows, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n_rows; i++)
    c_count[i] = 0;

  if (ma->separate_diag) {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      assert(g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]);
      cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
      assert(g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]);
      if (g_c_id != g_r_id)
        c_count[l_r_id] += 1;
    }

  }
  else {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_gnum_t g_r_id = ma->g_rc_id[i*2];
      cs_gnum_t g_c_id = ma->g_rc_id[i*2 + 1];
      assert(g_r_id >= ma->l_range[0] && g_r_id < ma->l_range[1]);
      cs_lnum_t l_r_id = g_r_id - ma->l_range[0];
      assert(g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]);
      c_count[l_r_id] += 1;
    }

  }

  /* Build index and reset count */

  ma->_r_idx[0] = 0;

  for (cs_lnum_t i = 0; i < n_rows; i++) {
    ma->_r_idx[i+1] = ma->_r_idx[i] + c_count[i];
    c_count[i] = 0;
  }

  /* Allocate structures holding data */

  CS_MALLOC_HD(ma->_c_id, ma->_r_idx[n_rows], cs_lnum_t, alloc_mode);
  ma->c_id = ma->_c_id;

  /* Now fill data: local part is determined (will be cleaned),
     and list of distant rows defined */

  if (ma->separate_diag) {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_lnum_t l_r_id = ma->g_rc_id[i*2] - ma->l_range[0];
      if (ma->g_rc_id[i*2+1] != ma->g_rc_id[i*2]) {
        ma->_c_id[ma->_r_idx[l_r_id] + c_count[l_r_id]]
          = ma->g_rc_id[i*2+1] - ma->l_range[0];
        c_count[l_r_id] += 1;
      }
    }

  }
  else {

    for (cs_lnum_t i = 0; i < ma->size; i++) {
      cs_lnum_t l_r_id = ma->g_rc_id[i*2] - ma->l_range[0];
      ma->_c_id[ma->_r_idx[l_r_id] + c_count[l_r_id]]
        = ma->g_rc_id[i*2+1] - ma->l_range[0];
      c_count[l_r_id] += 1;
    }

  }

  BFT_FREE(c_count);

  /* Set global number of rows (will be updated in parallel) */

  ma->n_g_rows = n_rows;
  ma->n_rows = n_rows;

  /* Sort and compact elements for final structure */

  _sort_and_compact_local(ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build diagonal index array when conversion between included
 *        and separate diagonal is required between a values assembler
 *        and its matching matrix assembler.
 *
 * \param[in, out]  pointer to matrix values assembler structure
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_diag_idx(cs_matrix_assembler_values_t  *mav)
{
  if (mav->diag_idx != NULL)
    return;

  const cs_matrix_assembler_t *ma = mav->ma;

  if (ma->separate_diag == mav->separate_diag)
    return;

  BFT_MALLOC(mav->diag_idx, ma->n_rows, cs_lnum_t);

  /*
    Two cases:
    - assembler has separate diagonal but values assembler does not:
      determine after which index rows should be shifted
    - assembler has included diagonal but values assembler does not:
      determine after which index rows should be shifted
  */

  if (ma->separate_diag) {
    for (cs_lnum_t i = 0; i < ma->n_rows; i++) {
      cs_lnum_t j = ma->r_idx[i], k = ma->r_idx[i+1];
      while (j < k) {
        if (ma->c_id[j] > i)
          k = j;
        j++;
      }
      mav->diag_idx[i] = k - ma->r_idx[i];
    }
  }
  else if (mav->separate_diag) {
    for (cs_lnum_t i = 0; i < ma->n_rows; i++) {
      cs_lnum_t j = ma->r_idx[i], k = ma->r_idx[i+1];
      while (j < k) {
        if (ma->c_id[j] == i)
          k = j;
        j++;
      }
      mav->diag_idx[i] = k - ma->r_idx[i];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add coefficient values based on row indexes when the diagonal
 *        must be shifted between the assembler and the values assembler.
 *
 * This function may be called by different threads if the caller
 * ensures different threads will not write to the same final parts of
 * the matrix coefficients.
 *
 * \param[in, out]  mav      pointer to matrix assembler values structure
 * \param[in]       n        number of values to add
 * \param[in]       stride   matrix block stride
 * \param[in]       row_id   matrix row ids
 * \param[in]       col_idx  matrix assembler column indexes
 * \param[in]       val      values to add
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_add_cnv_idx(cs_matrix_assembler_values_t  *mav,
                                     cs_lnum_t                      n,
                                     cs_lnum_t                      stride,
                                     const cs_lnum_t                row_id[],
                                     const cs_lnum_t                col_idx[],
                                     const cs_real_t                val[])
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  /* Case where the assembler has a separate diagonal, and values have
     an included diagonal */

  if (ma->separate_diag) {

    assert(mav->separate_diag == false);

    cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

    for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

      cs_lnum_t b_size = COEFF_GROUP_SIZE;
      if (i + COEFF_GROUP_SIZE > n)
        b_size = n - i;

      for (cs_lnum_t j = 0; j < b_size; j++) {
        cs_lnum_t r_id = row_id[i+j];
        cs_lnum_t c_idx = col_idx[i+j];
        if (r_id < 0)
          continue;
        if (c_idx == -1)
          s_col_idx[j] = mav->diag_idx[r_id];
        else if (c_idx < mav->diag_idx[r_id])
          s_col_idx[j] = c_idx;
        else
          s_col_idx[j] = c_idx + 1;
      }

      mav->add_values(mav->matrix,
                      b_size,
                      stride,
                      row_id + i,
                      s_col_idx,
                      val + (i*stride));

    }

  }

  /* Case where the assembler has an included diagonal, and values have
     a separate diagonal */

  else {

    assert(ma->separate_diag == false);
    assert(mav->separate_diag == true);

    cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

    for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

      cs_lnum_t b_size = COEFF_GROUP_SIZE;
      if (i + COEFF_GROUP_SIZE > n)
        b_size = n - i;

      for (cs_lnum_t j = 0; j < b_size; j++) {
        cs_lnum_t r_id = row_id[i+j];
        cs_lnum_t c_idx = col_idx[i+j];
        if (r_id < 0)
          continue;
        if (c_idx < mav->diag_idx[r_id])
          s_col_idx[j] = c_idx;
        else if (c_idx == mav->diag_idx[r_id])
          s_col_idx[j] = -1;
        else
          s_col_idx[j] = c_idx - 1;
      }

      mav->add_values(mav->matrix,
                      b_size,
                      stride,
                      row_id + i,
                      s_col_idx,
                      val + (i*stride));

    }

  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add coefficient values based on local row ids and global column ids,
 *        using a local addition function.
 *
 * This function may be called by different threads if the caller
 * ensures different threads will not write to the same final parts of
 * the matrix coefficients.
 *
 * \param[in, out]  mav       pointer to matrix assembler values structure
 * \param[in]       n         number of values to add
 * \param[in]       stride    matrix block stride
 * \param[in]       row_id    matrix row ids
 * \param[in]       col_g_id  global column ids
 * \param[in]       val       values to add
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_add_lg(cs_matrix_assembler_values_t  *mav,
                                cs_lnum_t                      n,
                                cs_lnum_t                      stride,
                                const cs_lnum_t                row_id[],
                                const cs_gnum_t                col_g_id[],
                                const cs_real_t                val[])
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  assert(mav->add_values != NULL);

  cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

  for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

    cs_lnum_t b_size = COEFF_GROUP_SIZE;
    if (i + COEFF_GROUP_SIZE > n)
      b_size = n - i;

    for (cs_lnum_t j = 0; j < b_size; j++) {

      cs_lnum_t k = i+j;

      cs_lnum_t l_r_id = row_id[k];
      cs_gnum_t g_c_id = col_g_id[k];

      /* Local part */

      if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {

        cs_lnum_t l_c_id = g_c_id - ma->l_range[0];
        cs_lnum_t n_cols =  ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];

        /* remark: we could slightly reduce the binary search by subtracting
           the last (ma->d_r_idx[l_r_id+1] - ma->d_r_idx[l_r_id]) from
           the column count, but this would require a separate function
           for the serial case where that array may not exist */

        s_col_idx[j] = _l_id_binary_search(n_cols,
                                           l_c_id,
                                           ma->c_id + ma->r_idx[l_r_id]);

        assert(s_col_idx[j] > -1 || (ma->separate_diag && l_c_id == l_r_id));

      }

      /* Distant part */

#if defined(HAVE_MPI)

      else {

        cs_lnum_t n_cols = ma->d_r_idx[l_r_id+1] - ma->d_r_idx[l_r_id];

        cs_lnum_t d_c_idx = _g_id_binary_find(n_cols,
                                              g_c_id,
                                              ma->d_g_c_id + ma->d_r_idx[l_r_id]);

        /* column ids start and end of local row, so add r_idx[row_id+1] */
        s_col_idx[j] = d_c_idx + ma->r_idx[l_r_id] + 1;

      }

#endif /* HAVE_MPI */

    }

    /* Case where assembler and values assembler have the same
       "separate diagonal" behavior */

    if (ma->separate_diag == mav->separate_diag)
      mav->add_values(mav->matrix,
                      b_size,
                      stride,
                      row_id + i,
                      s_col_idx,
                      val + (i*stride));

    else
      _matrix_assembler_values_add_cnv_idx(mav,
                                           b_size,
                                           stride,
                                           row_id + i,
                                           s_col_idx,
                                           val + (i*stride));

  }
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add coefficient values based on local row ids and column ids,
 *        using a global addition function.
 *
 * This function may be called by different threads if the caller
 * ensures different threads will not write to the same final parts of
 * the matrix coefficients.
 *
 * \param[in, out]  mav     pointer to matrix assembler values structure
 * \param[in]       n       number of values to add
 * \param[in]       stride  matrix block stride
 * \param[in]       row_id  matrix row ids
 * \param[in]       col_id  matrix assembler column ids
 * \param[in]       val     values to add
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_add_ll_g(cs_matrix_assembler_values_t  *mav,
                                  cs_lnum_t                      n,
                                  cs_lnum_t                      stride,
                                  const cs_lnum_t                row_id[],
                                  const cs_lnum_t                col_id[],
                                  const cs_real_t                val[])

{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_gnum_t row_g_id[COEFF_GROUP_SIZE];
  cs_gnum_t col_g_id[COEFF_GROUP_SIZE];

  for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

    cs_lnum_t b_size = COEFF_GROUP_SIZE;
    if (i + COEFF_GROUP_SIZE > n)
      b_size = n - i;

    for (cs_lnum_t j = 0; j < b_size; j++) {

      cs_lnum_t k = i+j;

      row_g_id[j] = row_id[k] + ma->l_range[0];

      cs_lnum_t l_c_id = col_id[k];
      if (l_c_id < ma->n_rows)
        col_g_id[j] = row_id[i+j] + ma->l_range[0];
      else
        col_g_id[j] = ma->e_g_id[l_c_id - ma->n_rows];

    }

    mav->add_values_g(mav->matrix,
                      b_size,
                      stride,
                      row_g_id,
                      col_g_id,
                      val + (i*stride));

  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add coefficient values based on local row ids and column indexes,
 *        using a global addition function.
 *
 * This function may be called by different threads if the caller ensures
 * different threads will not write to the same final parts of the matrix
 * coefficients.
 *
 * \param[in, out]  mav      pointer to matrix assembler values structure
 * \param[in]       n        number of values to add
 * \param[in]       stride   matrix block stride
 * \param[in]       row_id   matrix row ids
 * \param[in]       col_idx  matrix assembler column indexes
 * \param[in]       val      values to add
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_add_llx_g(cs_matrix_assembler_values_t  *mav,
                                   cs_lnum_t                      n,
                                   cs_lnum_t                      stride,
                                   const cs_lnum_t                row_id[],
                                   const cs_lnum_t                col_idx[],
                                   const cs_real_t                val[])

{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_gnum_t row_g_id[COEFF_GROUP_SIZE];
  cs_gnum_t col_g_id[COEFF_GROUP_SIZE];

  for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

    cs_lnum_t b_size = COEFF_GROUP_SIZE;
    if (i + COEFF_GROUP_SIZE > n)
      b_size = n - i;

    for (cs_lnum_t j = 0; j < b_size; j++) {

      cs_lnum_t k = i+j;
      cs_lnum_t l_r_id = row_id[k];
      cs_lnum_t l_c_id;

      if (col_idx[k] > -1) {
        l_c_id = ma->c_id[ma->r_idx[l_r_id] + col_idx[k]];
        assert(l_c_id > -1);
      }
      else
        l_c_id = l_r_id;

      row_g_id[j] = l_r_id + ma->l_range[0];
      if (l_c_id < ma->n_rows)
        col_g_id[j] = l_c_id + ma->l_range[0];
      else
        col_g_id[j] = ma->e_g_id[l_c_id - ma->n_rows];

    }

    mav->add_values_g(mav->matrix,
                      b_size,
                      stride,
                      row_g_id,
                      col_g_id,
                      val + (i*stride));

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add coefficient values based on local row ids and global column ids,
 *        using a global addition function.
 *
 * This function may be called by different threads if the caller
 * ensures different threads will not write to the same final parts of
 * the matrix coefficients.
 *
 * \param[in, out]  mav       pointer to matrix assembler values structure
 * \param[in]       n         number of values to add
 * \param[in]       stride    matrix block stride
 * \param[in]       row_id    matrix row ids
 * \param[in]       col_g_id  matrix global column ids
 * \param[in]       val       values to add
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_assembler_values_add_lg_g(cs_matrix_assembler_values_t  *mav,
                                  cs_lnum_t                      n,
                                  cs_lnum_t                      stride,
                                  const cs_lnum_t                row_id[],
                                  const cs_gnum_t                col_g_id[],
                                  const cs_real_t                val[])
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_gnum_t row_g_id[COEFF_GROUP_SIZE];

  for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

    cs_lnum_t b_size = COEFF_GROUP_SIZE;
    if (i + COEFF_GROUP_SIZE > n)
      b_size = n - i;

    for (cs_lnum_t j = 0; j < b_size; j++)
      row_g_id[j] = row_id[i+j] + ma->l_range[0];

    mav->add_values_g(mav->matrix,
                      b_size,
                      stride,
                      row_g_id,
                      col_g_id + i,
                      val + (i*stride));

  }
}

#endif /* HAVE_MPI */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
 * \return  pointer to created matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create(const cs_gnum_t  l_range[2],
                           bool             separate_diag)
{
  cs_matrix_assembler_t *ma = NULL;

  BFT_MALLOC(ma, 1, cs_matrix_assembler_t);

  ma->separate_diag = separate_diag;

  ma->flags = CS_MATRIX_DISTANT_ROW_USE_COL_IDX;

  ma->l_range[0] = l_range[0];
  ma->l_range[1] = l_range[1];

  ma->n_g_rows = 0;
  ma->n_rows = 0;

  ma->size = 0;
  ma->max_size = 0;

  ma->r_idx = NULL;
  ma->c_id = NULL;
  ma->_r_idx = NULL;
  ma->_c_id = NULL;

  ma->d_r_idx = NULL;
  ma->d_g_c_id = NULL;

  ma->g_rc_id = NULL;

#if defined(HAVE_MPI)

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
  ma->coeff_recv_col_idx = NULL;
  ma->coeff_recv_col_g_id = NULL;

  ma->comm = cs_glob_mpi_comm;

  ma->n_ranks_init[0] = 0;
  ma->n_ranks_init[1] = 0;

#endif /* HAVE_MPI */

  ma->halo = NULL;
  ma->_halo = NULL;

  ma->n_e_g_ids = 0;
  ma->e_g_id = NULL;

  return ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a matrix assembler structure based on a given connectivity
 *        and associated halo structure.
 *
 * The assembler may not be later modified when built this way.
 *
 * The connectivity arrays and halo will be shared by the caller, so their
 * life cycle must be at least as long as the matrix assembler structure.
 *
 * \param[in]  n_rows         number fo local rows
 * \param[in]  separate_diag  if true, diagonal terms are handled separately
 * \param[in]  row_idx        matrix row index
 * \param[in]  col_id         matrix column indexes
 * \param[in]  halo           associated halo structure
 *
 * \return  pointer to created matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_t *
cs_matrix_assembler_create_from_shared(cs_lnum_t         n_rows,
                                       bool              separate_diag,
                                       const cs_lnum_t   row_idx[],
                                       const cs_lnum_t   col_id[],
                                       const cs_halo_t  *halo)
{
  cs_gnum_t l_range[2] = {0, n_rows};
  cs_gnum_t n_g_rows = n_rows;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t g_shift;
    cs_gnum_t l_shift = n_rows;
    MPI_Scan(&l_shift, &g_shift, 1, CS_MPI_GNUM, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&l_shift, &n_g_rows, 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
    l_range[0] = g_shift - l_shift;
    l_range[1] = g_shift;
  }
#endif

  /* Base creation */

  cs_matrix_assembler_t *ma = cs_matrix_assembler_create(l_range,
                                                         separate_diag);

  ma->n_g_rows = n_g_rows;
  ma->n_rows = n_rows;

  ma->r_idx = row_idx;
  ma->c_id = col_id;

  ma->halo = halo;

  /* Additional required parallel data */

  if (ma->halo != NULL) {

    cs_gnum_t *t_g_id;
    BFT_MALLOC(ma->e_g_id, ma->halo->n_elts[0], cs_gnum_t);
    BFT_MALLOC(t_g_id, ma->n_rows + ma->halo->n_elts[0], cs_gnum_t);

    for (cs_lnum_t i = 0; i < ma->n_rows; i++)
      t_g_id[i] = (cs_gnum_t)i + ma->l_range[0];

    cs_halo_sync_untyped(ma->halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_gnum_t),
                         t_g_id);

    ma->n_e_g_ids = ma->halo->n_elts[0];
    for (cs_lnum_t i = 0; i < ma->n_e_g_ids; i++)
      ma->e_g_id[i] = t_g_id[ma->n_rows + i];

    BFT_FREE(t_g_id);

    /* Build distant columns index and global ids */

    BFT_MALLOC(ma->d_r_idx, ma->n_rows+1, cs_lnum_t);

    ma->d_r_idx[0] = 0;
    for (cs_lnum_t i = 0; i < ma->n_rows; i++) {
      cs_lnum_t n_e_cols = 0;
      cs_lnum_t n_cols = ma->r_idx[i+1] - ma->r_idx[i];
      const cs_lnum_t *c_id = ma->c_id + ma->r_idx[i];
      for (cs_lnum_t j = 0; j < n_cols; j++) {
        if (c_id[j] >= ma->n_rows)
          n_e_cols++;
      }
      ma->d_r_idx[i+1] = n_e_cols;
    }

    for (cs_lnum_t i = 0; i < ma->n_rows; i++)
      ma->d_r_idx[i+1] += ma->d_r_idx[i];

    BFT_MALLOC(ma->d_g_c_id, ma->d_r_idx[ma->n_rows], cs_gnum_t);
    for (cs_lnum_t i = 0; i < ma->n_rows; i++) {
      cs_lnum_t offset = ma->d_r_idx[i];
      cs_lnum_t n_cols = ma->r_idx[i+1] - ma->r_idx[i];
      const cs_lnum_t *c_id = ma->c_id + ma->r_idx[i];
      for (cs_lnum_t j = 0; j < n_cols; j++) {
        if (c_id[j] >= ma->n_rows) {
          assert(c_id[j] - ma->n_rows < ma->n_e_g_ids);
          ma->d_g_c_id[offset++] = ma->e_g_id[c_id[j] - ma->n_rows];
        }
      }
    }

  }

  return ma;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a matrix assembler structure.
 *
 * \param[in, out]  ma  pointer to matrix assembler structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_destroy(cs_matrix_assembler_t  **ma)
{
  if (ma != NULL && *ma != NULL) {
    cs_matrix_assembler_t *_ma = *ma;

    BFT_FREE(_ma->e_g_id);

    if (_ma->_halo != NULL)
      cs_halo_destroy(&(_ma->_halo));

#if defined(HAVE_MPI)

    BFT_FREE(_ma->coeff_recv_col_g_id);
    BFT_FREE(_ma->coeff_recv_col_idx);
    BFT_FREE(_ma->coeff_recv_row_id);

    BFT_FREE(_ma->coeff_rank_recv_index);
    BFT_FREE(_ma->coeff_rank_send_index);

    BFT_FREE(_ma->coeff_send_col_g_id);
    BFT_FREE(_ma->coeff_send_row_g_id);
    BFT_FREE(_ma->coeff_send_index);
    BFT_FREE(_ma->coeff_rank);

#endif /* HAVE_MPI */

    BFT_FREE(_ma->g_rc_id);

    BFT_FREE(_ma->d_g_c_id);
    BFT_FREE(_ma->d_r_idx);

    CS_FREE_HD(_ma->_c_id);
    CS_FREE_HD(_ma->_r_idx);

    BFT_FREE(*ma);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set option flags for a given matrix assembler structure.
 *
 * When used, this function should be called before defining entries
 * (in case of future option additions) and in any case before computing
 * the final structure.
 *
 * Flags are defined as a sum (bitwise or) of constants, described below:
 * - \ref CS_MATRIX_DISTANT_ROW_USE_COL_IDX indicates that column indexes
 *   matching distant row data are computed and maintained, so the
 *   associated coefficient contributions can be added more efficiently.
 * - \ref CS_MATRIX_DISTANT_ROW_USE_COL_G_ID  indicates that columns
 *   global ids matching distant row data are maintained directly,
 *   so no lookup to a global id description of matrix columns is needed;
 *   this option is useful only when using an external library allowing
 *   incremental setting or accumulation of coefficients using global
 *   column ids, such as with PETSc's MatSetValues. When building matrices
 *   locally first (which includes most cases, whether internal matrices,
 *   or using library conversion or import functions such as PETSc's
 *   MatCreateMPIAIJWithArrays) this is slightly less efficient, as
 *   column ids will need to be matched to each row's columns a second
 *   time for those (exchanged) coefficients.
 * - \ref CS_MATRIX_EXTERNAL_HALO indicates that we do not need to
 *   build an associated halo structure, as it will be built externally
 *   (i.e. by an external library such as PETSc, HYPRE, ...) using its
 *   own equivalent structures.
 *
 * \param[in, out]  ma     pointer to matrix assembler structure
 * \param[in]       flags  sum of matrix assembler flag constants
 *                         (bitwise or)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_set_options(cs_matrix_assembler_t  *ma,
                                int                     flags)
{
  ma->flags = flags;

  /* Ensure at least one distant row exchange option is set */

  if (   !(ma->flags & CS_MATRIX_DISTANT_ROW_USE_COL_IDX)
      && !(ma->flags & CS_MATRIX_DISTANT_ROW_USE_COL_G_ID))
    ma->flags = ma->flags | CS_MATRIX_DISTANT_ROW_USE_COL_IDX;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the option flags for a given matrix assembler structure.
 *
 * Flags are defined as a sum (bitwise or) of constants, described in
 * \ref cs_matrix_assembler_set_options.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  option flags (sum of integer constants) used y this structure
 */
/*----------------------------------------------------------------------------*/

int
cs_matrix_assembler_get_options(const cs_matrix_assembler_t  *ma)
{
  return ma->flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add entries to a matrix assembler structure.
 *
 * This function should be called by a single thread for a given assembler.
 *
 * \param[in, out]  ma        pointer to matrix assembler structure
 * \param[in]       n         number of entries
 * \param[in]       col_g_id  global column ids associated with entries
 * \param[in]       row_g_id  global row ids associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_add_g_ids(cs_matrix_assembler_t  *ma,
                              cs_lnum_t               n,
                              const cs_gnum_t         row_g_id[],
                              const cs_gnum_t         col_g_id[])
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
      _g_rc_id[i*2]   = row_g_id[i];
      _g_rc_id[i*2+1] = col_g_id[i];
    }
    ma->size += n;
  }
  else {
    cs_lnum_t j = 0;
    for (cs_lnum_t i = 0; i < n; i++) {
      if (   row_g_id[i] != col_g_id[i]
          || row_g_id[i] <  ma->l_range[0]
          || row_g_id[i] >= ma->l_range[1]) {
        _g_rc_id[j*2]   = row_g_id[i];
        _g_rc_id[j*2+1] = col_g_id[i];
        j += 1;
      }
    }
    ma->size += j;
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
cs_matrix_assembler_compute(cs_matrix_assembler_t  *ma)
{
#if defined(HAVE_MPI)
  if (ma->comm != MPI_COMM_NULL && ma->comm != MPI_COMM_SELF) {
    _matrix_assembler_compute_mpi(ma);
    return;
  }
#endif /* HAVE_MPI */

  _matrix_assembler_compute_local(ma);

  /* Non-null array to simplify binary search (rare diagonal-only case) */

  if (ma->c_id == NULL) {
    BFT_MALLOC(ma->_c_id, 1, cs_lnum_t);
    ma->c_id = ma->_c_id;
    ma->_c_id[0] = -1;
  }
}

#if defined(HAVE_MPI)

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
 * \brief Return a pointer to local global row range associated with a
 *        matrix assembler.
 *
 * \param[in]  ma   pointer to matrix assembler structure
 *
 * \return  pointer to local range
 */
/*----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_assembler_get_l_range(const cs_matrix_assembler_t  *ma)
{
  return ma->l_range;
}

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
  return ma->halo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if the matrix assembler is based on a separate diagonal.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  true if the associated structure has a separate diagonal,
 *          false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_matrix_assembler_get_separate_diag(const cs_matrix_assembler_t  *ma)
{
  return ma->separate_diag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of rows associated with a matrix assembler.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  number of rows associated with matrix assembler
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_assembler_get_n_rows(const cs_matrix_assembler_t  *ma)
{
  return ma->n_rows;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the global number of rows associated with a matrix assembler.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  number of rows associated with matrix assembler
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_matrix_assembler_get_n_g_rows(const cs_matrix_assembler_t  *ma)
{
  return ma->n_g_rows;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of columns associated with a matrix assembler.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  number of columns associated with matrix assembler
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_matrix_assembler_get_n_columns(const cs_matrix_assembler_t  *ma)
{
#if defined(HAVE_MPI)
  return (ma->n_rows + ma->n_e_g_ids);
#else
  return ma->n_rows;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a row index associated with a matrix assembler.
 *
 * This index is a CSR type index relative to all columns, including or
 * excluding the diagonal depending on the \c separate_diag option used
 * when creating the matrix assembler. The matching column ids can be
 * obtained using \ref cs_matrix_assembler_get_col_ids.
 *
 * \warning the returned index is valid only as long as the matrix assembly
 * structure exists.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  pointer to matrix structure row index
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_assembler_get_row_index(const cs_matrix_assembler_t  *ma)
{
  return ma->r_idx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a column ids associated with a matrix assembler.
 *
 * These ids relative to all columns of successive rows, with columns
 * ordered by local id, so local columns appear first, distant columns
 * after. Depending on the \c separate_diag option used when creating the
 * matrix assembler, the columns may include the diagonal or not.
 * The matching index can be obtained using
 * \ref cs_matrix_assembler_get_row_index.
 *
 * \warning the returned index is valid only as long as the matrix assembly
 * structure exists.
 *
 * \param[in]  ma  pointer to matrix assembler structure
 *
 * \return  pointer to matrix structure row index
 */
/*----------------------------------------------------------------------------*/

const cs_lnum_t *
cs_matrix_assembler_get_col_ids(const cs_matrix_assembler_t  *ma)
{
  return ma->c_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return info on the number of neighbor ranks a matrix assembler
 *        may communicate with.
 *
 * \param[in]   ma   pointer to matrix assembler structure
 * \param[out]  rc   rank counts; the 4 values are:
 *                   - 0 number of communicating ranks for vector update (halo)
 *                   - 1 number of communicating ranks for matrix assembly
 *                   - 2 maximum number of communicating ranks for halo
 *                       construction (assumed ranks algorithm)
 *                   - 3 maximum number of communicating ranks for matrix
 *                       assembly determination (assumed ranks algorithm)
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_get_rank_counts(const cs_matrix_assembler_t  *ma,
                                    int                           rc[4])
{
  rc[0] = 0;
  if (ma->halo != NULL)
    rc[0] = ma->halo->n_c_domains;
#if defined(HAVE_MPI)
  rc[1] = ma->n_coeff_ranks;
  rc[2] = ma->n_ranks_init[1];
  rc[3] = ma->n_ranks_init[0];
#else
  rc[1] = 0;
  rc[2] = 0;
  rc[3] = 0;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log rank counts for a given matrix assembler.
 *
 * \param[in]  ma      pointer to matrix assembler structure
 * \param[in]  log_id  log type
 * \param[in]  name    name of this assembler
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_log_rank_counts(const cs_matrix_assembler_t  *ma,
                                    cs_log_t                      log_id,
                                    const char                   *name)
{
  cs_log_printf(log_id,
                _("\nNeighbor rank counts for matrix assembler: %s\n"
                  "-----------------------------------------\n"),
                name);

  const char *count_name[]
    = {N_("Neighbor ranks for vector update (halo)"),
       N_("Neighbor ranks for matrix assembly"),
       N_("Maximum neighbor ranks during halo construction"),
       N_("Maximum neighbor ranks during assembly determination")};

  int counts[4];

  cs_matrix_assembler_get_rank_counts(ma, counts);

  for (int i = 0; i < 4; i++) {

    char ul[120];
    size_t j = 0;
    size_t u_count = cs_log_strlen(_(count_name[i]));
    for (j = 0; j < u_count && j < 119; j++)
      ul[j] = '-';
    ul[j] = '\0';

    cs_log_printf(log_id, "\n  %s:\n  %s\n\n", _(count_name[i]), ul);

    _display_rank_histogram(log_id, counts[i]);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize a matrix assembler values structure.
 *
 * The associated values will initially be set to zero.
 *
 * This is a low-level function, which should be called by a simpler
 * function (\ref cs_matrix_assembler_values_init) which provides
 * the necessary function pointers.
 *
 * \warning  The matrix pointer must point to valid data when the selection
 *           function is called, so the life cycle of the data pointed to
 *           should be at least as long as that of the assembler values
 *           structure. In a similar manner, the life cycle of the associated
 *           matrix assembler must also be at least as long as that
 *           of the assembler values structure.
 *
 * \param[in]       ma        associated matrix assembly structure
 * \param[in]       sep_diag  true if diagonal terms are stored separately
 * \param[in]       db_size   diagonal block sizes
 * \param[in]       eb_size   extra-diagonal block sizes
 * \param[in, out]  matrix    untyped pointer to matrix description structure
 * \param[in]       init      pointer to matrix coefficients
 *                            initialization function
 * \param[in]       add       pointer to matrix coefficients addition from
 *                            local ids function
 * \param[in]       add_g     pointer to matrix coefficients addition from
 *                            global ids function
 * \param[in]       begin     pointer to matrix coefficients assembly
 *                            start function (optional)
 * \param[in]       end       pointer to matrix coefficients assembly
 *                            end function (optional)
 *
 * \return  pointer to initialized matrix assembler structure;
 */
/*----------------------------------------------------------------------------*/

cs_matrix_assembler_values_t *
cs_matrix_assembler_values_create(const cs_matrix_assembler_t          *ma,
                                  bool                                  sep_diag,
                                  cs_lnum_t                             db_size,
                                  cs_lnum_t                             eb_size,
                                  void                                 *matrix,
                                  cs_matrix_assembler_values_init_t    *init,
                                  cs_matrix_assembler_values_add_t     *add,
                                  cs_matrix_assembler_values_add_g_t   *add_g,
                                  cs_matrix_assembler_values_begin_t   *begin,
                                  cs_matrix_assembler_values_end_t     *end)
{
  cs_matrix_assembler_values_t *mav;

  BFT_MALLOC(mav, 1, cs_matrix_assembler_values_t);

  mav->ma = ma;

  mav->separate_diag = sep_diag;
  mav->final_assembly = false;

  mav->db_size = db_size;
  mav->eb_size = eb_size;

  mav->diag_idx = NULL;

  mav->matrix = matrix;

  mav->init = init;
  mav->add_values = add;
  mav->add_values_g = add_g;
  mav->assembly_begin = begin;
  mav->assembly_end = end;

  const cs_lnum_t eb_size_2 = eb_size*eb_size;

#if defined(HAVE_MPI)

  cs_lnum_t  alloc_size = ma->coeff_send_size * eb_size_2;

  BFT_MALLOC(mav->coeff_send, alloc_size, cs_real_t);

  for (cs_lnum_t i = 0; i < alloc_size; i++)
    mav->coeff_send[i] = 0;

#endif /* HAVE_MPI */

  /* Build diagonal index if it will be needed */

  if (mav->separate_diag != ma->separate_diag)
    _matrix_assembler_values_diag_idx(mav);

  if (mav->init != NULL)
    mav->init(mav->matrix, mav->db_size, mav->eb_size);

  return mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize matrix assembler values structure.
 *
 * When this function returns, the assembler values structure has been
 * destroyed, and the associated matrix is fully assembled (i.e. ready to
 * use).
 *
 * \param[in, out]  mav  pointer to matrix assembler values structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_finalize(cs_matrix_assembler_values_t  **mav)
{
  if (mav == NULL)
    return;

  cs_matrix_assembler_values_t  *_mav = *mav;
  if (_mav == NULL)
    return;

  if (_mav->final_assembly == false)
    cs_matrix_assembler_values_done(_mav);

  if (_mav->assembly_end != NULL)
    _mav->assembly_end(_mav->matrix);

  BFT_FREE(*mav);
  *mav = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using local
 *        row and column ids.
 *
 * If the matching matrix coefficients use a block structure, the
 * values passed to this function should use the same block size.
 *
 * Note also that if the matrix has different diagonal and extradiagonal
 * block sizes, separate calls to this function should be used for both
 * types of coefficients; in addition, in this case, diagonal coefficients
 * should only be provided by the owning rank (this also impacts how
 * the associated matrix assembler structure is defined).
 *
 * This function may be called by different threads, as long those threads
 * do not add contributions to the same rows (assuming caution is taken
 * in the case of external libraries so that their builder functions
 * have the same property).
 *
 * \param[in, out]  mav     pointer to matrix assembler values structure
 * \param[in]       n       number of entries
 * \param[in]       col_id  local column ids associated with entries
 * \param[in]       row_id  local row ids associated with entries
 * \param[in]       val     values associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_add(cs_matrix_assembler_values_t  *mav,
                               cs_lnum_t                      n,
                               const cs_lnum_t                row_id[],
                               const cs_lnum_t                col_id[],
                               const cs_real_t                val[])
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_lnum_t stride = 0;

  if (n < 1)
    return;

  /* Base stride on first type of value encountered */

  if (row_id[0] == col_id[0])
    stride = mav->db_size * mav->db_size;
  else
    stride = mav->eb_size * mav->eb_size;

  /* Case where we compute a column index first */

  if (mav->add_values != NULL) {

    cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

    for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

      cs_lnum_t b_size = COEFF_GROUP_SIZE;
      if (i + COEFF_GROUP_SIZE > n)
        b_size = n - i;

      for (cs_lnum_t j = 0; j < b_size; j++) {

        cs_lnum_t k = i+j;

        cs_lnum_t l_r_id = row_id[k];
        cs_lnum_t l_c_id = col_id[k];

        cs_lnum_t n_cols = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];

        s_col_idx[j] = _l_id_binary_search(n_cols,
                                           col_id[k],
                                           ma->c_id + ma->r_idx[l_r_id]);

        assert(s_col_idx[j] > -1 || (ma->separate_diag && l_c_id == l_r_id));

      }

      if (ma->separate_diag == mav->separate_diag)
        mav->add_values(mav->matrix,
                        b_size,
                        stride,
                        row_id + i,
                        s_col_idx,
                        val + (i*stride));

      else
        _matrix_assembler_values_add_cnv_idx(mav,
                                             b_size,
                                             stride,
                                             row_id + i,
                                             s_col_idx,
                                             val + (i*stride));

    }

  }

  else { /* use a global id-based function */

    assert(mav->add_values_g != NULL);

    _matrix_assembler_values_add_ll_g(mav,
                                      n,
                                      stride,
                                      row_id,
                                      col_id,
                                      val);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add values to a matrix assembler values structure using global
 *        row and column ids.
 *
 * If the matching matrix coefficients use a block structure, the
 * values passed to this function should use the same block size.
 *
 * Note also that if the matrix has different diagonal and extradiagonal
 * block sizes, separate calls to this function should be used for both
 * types of coefficients; in addition, in this case, diagonal coefficients
 * should only be provided by the owning rank (this also impacts how
 * the associated matrix assembler structure is defined).
 *
 * This function may be called by different threads, as long those threads
 * do not add contributions to the same rows (assuming caution is taken
 * in the case of external libraries so that their builder functions
 * have the same property).
 *
 * \param[in, out]  mav       pointer to matrix assembler values structure
 * \param[in]       n         number of entries
 * \param[in]       col_g_id  global column ids associated with entries
 * \param[in]       row_g_id  global row ids associated with entries
 * \param[in]       val       values associated with entries
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_add_g(cs_matrix_assembler_values_t  *mav,
                                 cs_lnum_t                      n,
                                 const cs_gnum_t                row_g_id[],
                                 const cs_gnum_t                col_g_id[],
                                 const cs_real_t                val[])
{
  const cs_matrix_assembler_t  *ma = mav->ma;

  cs_lnum_t stride = 0;

  if (n < 1)
    return;

  /* Base stride on first type of value encountered */

  if (row_g_id[0] == col_g_id[0])
    stride = mav->db_size * mav->db_size;
  else
    stride = mav->eb_size * mav->eb_size;

  cs_gnum_t s_g_row_id[COEFF_GROUP_SIZE];
  cs_gnum_t s_g_col_id[COEFF_GROUP_SIZE];

  for (cs_lnum_t i = 0; i < n; i+= COEFF_GROUP_SIZE) {

    cs_lnum_t b_size = COEFF_GROUP_SIZE;
    if (i + COEFF_GROUP_SIZE > n)
      b_size = n - i;

    for (cs_lnum_t j = 0; j < b_size; j++) {

      cs_lnum_t k = i+j;

      cs_gnum_t g_r_id = row_g_id[k];

#if defined(HAVE_MPI)

      /* Case where coefficient is handled by other rank */

      if (g_r_id < ma->l_range[0] || g_r_id >= ma->l_range[1]) {

        s_g_row_id[j] = ma->l_range[1];
        s_g_col_id[j] = 0;

        cs_lnum_t e_r_id = _g_id_binary_find(ma->coeff_send_n_rows,
                                             g_r_id,
                                             ma->coeff_send_row_g_id);

        cs_lnum_t r_start = ma->coeff_send_index[e_r_id];
        cs_lnum_t n_e_rows = ma->coeff_send_index[e_r_id+1] - r_start;

        cs_lnum_t e_id =   r_start
                         +_g_id_binary_find(n_e_rows,
                                            col_g_id[k],
                                            ma->coeff_send_col_g_id + r_start);

        /* Now add values to send coefficients */

        for (cs_lnum_t l = 0; l < stride; l++)
#         pragma omp atomic
          mav->coeff_send[e_id*stride + l] += val[k*stride + l];

      }

      /* Standard case */

      else {

        s_g_row_id[j] = g_r_id;
        s_g_col_id[j] = col_g_id[k];

      }

#else

      s_g_row_id[j] = g_r_id;
      s_g_col_id[j] = col_g_id[k];

#endif /* HAVE_MPI */

    }

    if (mav->add_values_g != NULL) /* global id-based assembler function */
      mav->add_values_g(mav->matrix,
                        b_size,
                        stride,
                        s_g_row_id,
                        s_g_col_id,
                        val + (i*stride));

    else if (ma->d_r_idx != NULL) { /* local-id-based function, need to adapt */

      cs_lnum_t s_row_id[COEFF_GROUP_SIZE];
      cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

      for (cs_lnum_t j = 0; j < b_size; j++) {

        if (s_g_row_id[j] == ma->l_range[1]) { /* filter other-rank values */
          s_row_id[j] = -1;
          s_col_idx[j] = -1;
        }

        else { /* Handle values assigned to local rank */

          cs_lnum_t l_r_id = s_g_row_id[j] - ma->l_range[0];
          cs_gnum_t g_c_id = s_g_col_id[j];

          s_row_id[j] = l_r_id;

          cs_lnum_t n_l_cols = (  ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id]
                                - ma->d_r_idx[l_r_id+1] + ma->d_r_idx[l_r_id]);

          /* Local part */

          if (g_c_id >= ma->l_range[0] && g_c_id < ma->l_range[1]) {

            cs_lnum_t l_c_id = g_c_id - ma->l_range[0];

            s_col_idx[j] = _l_id_binary_search(n_l_cols,
                                               l_c_id,
                                               ma->c_id + ma->r_idx[l_r_id]);

            assert(s_col_idx[j] > -1 || (ma->separate_diag && l_c_id == l_r_id));

          }

          /* Distant part */

          else {

            cs_lnum_t n_cols = ma->d_r_idx[l_r_id+1] - ma->d_r_idx[l_r_id];

            cs_lnum_t d_c_idx = _g_id_binary_find(n_cols,
                                                  g_c_id,
                                                  ma->d_g_c_id + ma->d_r_idx[l_r_id]);

            /* column ids start and end of local row, so add n_l_cols */
            s_col_idx[j] = d_c_idx + n_l_cols;

          }

        }

      }

      if (ma->separate_diag == mav->separate_diag) {
        mav->add_values(mav->matrix,
                        b_size,
                        stride,
                        s_row_id,
                        s_col_idx,
                        val + (i*stride));
      }

      else
        _matrix_assembler_values_add_cnv_idx(mav,
                                             b_size,
                                             stride,
                                             s_row_id,
                                             s_col_idx,
                                             val + (i*stride));

    }

    else { /* local-only case or matrix part */

      assert(ma->d_r_idx == NULL);

      cs_lnum_t s_row_id[COEFF_GROUP_SIZE];
      cs_lnum_t s_col_idx[COEFF_GROUP_SIZE];

      for (cs_lnum_t j = 0; j < b_size; j++) {

        assert(s_g_row_id[j] < ma->l_range[1]);
        assert(s_g_col_id[j] < ma->l_range[1]);

        cs_lnum_t l_r_id = s_g_row_id[j] - ma->l_range[0];
        cs_lnum_t l_c_id = s_g_col_id[j] - ma->l_range[0];

        assert(   s_g_col_id[j] >= ma->l_range[0]
               && s_g_col_id[j] < ma->l_range[1]);

        s_row_id[j] = l_r_id;

        cs_lnum_t n_cols = ma->r_idx[l_r_id+1] - ma->r_idx[l_r_id];

        s_col_idx[j] = _l_id_binary_search(n_cols,
                                           l_c_id,
                                           ma->c_id + ma->r_idx[l_r_id]);

        assert(s_col_idx[j] > -1 || (ma->separate_diag && l_c_id == l_r_id));

      }

      if (ma->separate_diag == mav->separate_diag)
        mav->add_values(mav->matrix,
                        b_size,
                        stride,
                        s_row_id,
                        s_col_idx,
                        val + (i*stride));
      else
        _matrix_assembler_values_add_cnv_idx(mav,
                                             b_size,
                                             stride,
                                             s_row_id,
                                             s_col_idx,
                                             val + (i*stride));

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Start assembly of matrix values structure.
 *
 * The call to this function is always optional, and indicates that assembly
 * may begin. Calling it before \ref cs_matrix_assembler_values_finalize
 * may be useful when the underlying libraries can overlap communication
 * (or other operations such as matrix trasformation) and computation.
 *
 * \remark When values have been added to rows belonging to another rank,
 *         communication will occur here. Splitting this function could
 *         add another opportunity for overlap of computation and
 *         communication, but is not done as of 2016, as it would add
 *         complexity, with probably limited returns, as the effective
 *         asynchonous progress of MPI libraries is usually limited.
 *
 * \param[in, out]  mav  pointer to matrix assembler data structure
 */
/*----------------------------------------------------------------------------*/

void
cs_matrix_assembler_values_done(cs_matrix_assembler_values_t  *mav)
{
  if (mav == NULL)
    return;

  /* Exchange row data with other ranks if required */

#if defined(HAVE_MPI)

  const cs_matrix_assembler_t  *ma = mav->ma;

  if (ma->n_coeff_ranks > 0) {

    cs_real_t  *recv_coeffs = NULL;

    cs_lnum_t stride = mav->eb_size * mav->eb_size;

    MPI_Request *request = NULL;
    MPI_Status *status = NULL;

    BFT_MALLOC(recv_coeffs, ma->coeff_recv_size*stride, cs_real_t);

    BFT_MALLOC(request, ma->n_coeff_ranks*2, MPI_Request);
    BFT_MALLOC(status, ma->n_coeff_ranks*2, MPI_Status);

    int request_count = 0;
    int local_rank = cs_glob_rank_id;

    for (int i = 0; i < ma->n_coeff_ranks; i++) {
      cs_lnum_t l_size = (  ma->coeff_rank_recv_index[i+1]
                          - ma->coeff_rank_recv_index[i]) * stride;
      if (l_size > 0) {
        cs_lnum_t recv_shift = ma->coeff_rank_recv_index[i]*stride;
        MPI_Irecv(recv_coeffs + recv_shift,
                  l_size,
                  CS_MPI_REAL,
                  ma->coeff_rank[i],
                  local_rank,
                  ma->comm,
                  &(request[request_count++]));
      }
    }

    for (int i = 0; i < ma->n_coeff_ranks; i++) {
      cs_lnum_t l_size = (  ma->coeff_rank_send_index[i+1]
                          - ma->coeff_rank_send_index[i]) * stride;
      if (l_size > 0) {
        cs_lnum_t send_shift = ma->coeff_rank_send_index[i]*stride;
        MPI_Isend(mav->coeff_send + send_shift,
                  l_size,
                  CS_MPI_REAL,
                  ma->coeff_rank[i],
                  ma->coeff_rank[i],
                  ma->comm,
                  &(request[request_count++]));
      }
    }

    MPI_Waitall(request_count, request, status);

    BFT_FREE(request);
    BFT_FREE(status);

    /* Now add coefficients to local rows */

    if (ma->coeff_recv_size > 0) {

      /* Using local id-based function */

      if (mav->add_values != NULL) {

        assert(ma->coeff_recv_row_id != NULL);
        assert(ma->coeff_recv_col_idx != NULL);

        /* Case where columns are already pre-indexed */

        if (ma->coeff_recv_col_idx != NULL) {
          if (ma->separate_diag == mav->separate_diag)
            mav->add_values(mav->matrix,
                            ma->coeff_recv_size,
                            stride,
                            ma->coeff_recv_row_id,
                            ma->coeff_recv_col_idx,
                            recv_coeffs);
          else
            _matrix_assembler_values_add_cnv_idx(mav,
                                                 ma->coeff_recv_size,
                                                 stride,
                                                 ma->coeff_recv_row_id,
                                                 ma->coeff_recv_col_idx,
                                                 recv_coeffs);
        }
        else {
          assert(ma->coeff_recv_col_g_id != NULL);
          _matrix_assembler_values_add_lg(mav,
                                          ma->coeff_recv_size,
                                          stride,
                                          ma->coeff_recv_row_id,
                                          ma->coeff_recv_col_g_id,
                                          recv_coeffs);
        }

      }

      /* Using global id-based function */

      else {
        assert(mav->add_values_g != NULL);

        if (ma->coeff_recv_col_g_id != NULL)
          _matrix_assembler_values_add_lg_g(mav,
                                            ma->coeff_recv_size,
                                            stride,
                                            ma->coeff_recv_row_id,
                                            ma->coeff_recv_col_g_id,
                                            recv_coeffs);

        else {
          assert(ma->coeff_recv_col_idx);
          _matrix_assembler_values_add_llx_g(mav,
                                             ma->coeff_recv_size,
                                             stride,
                                             ma->coeff_recv_row_id,
                                             ma->coeff_recv_col_idx,
                                             recv_coeffs);
        }
      }

    }

    BFT_FREE(recv_coeffs);

  }

  BFT_FREE(mav->coeff_send);

#endif /* HAVE_MPI */

  BFT_FREE(mav->diag_idx);

  mav->final_assembly = true;

  if (mav->assembly_begin != NULL)
    mav->assembly_begin(mav->matrix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
