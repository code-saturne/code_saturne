/*============================================================================
 * Determine and update geometrical neighborhood information.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include <limits.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_box.h"
#include "fvm_box_tree.h"

#include "cs_order.h"
#include "cs_part_to_block.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_neighborhood.h"

#include "cs_all_to_all.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro and Type definitions
 *============================================================================*/

/* Neighborhood statistics (using box tree) */
/*------------------------------------------*/

typedef struct {

  int         dim;                     /* Layout dimension */

  /* The following fields have 3 global values:
     mean on ranks, minimum on ranks, and maximum on ranks */

  int         depth[3];                /* Tree depth */
  cs_lnum_t   n_leaves[3];             /* Number of leaves */
  cs_lnum_t   n_boxes[3];              /* Number of associated boxes */
  cs_lnum_t   n_threshold_leaves[3];   /* Number of leaves over threshold */
  cs_lnum_t   n_leaf_boxes[3];         /* Number of boxes per leaf */
  size_t      mem_used[3];             /* Memory used */
  size_t      mem_required[3];         /* Memory temporarily required */

} _box_tree_stats_t;

/* Main neighborhood structure */
/*-----------------------------*/

struct _fvm_neighborhood_t {

  cs_lnum_t         n_elts;          /* Number of elements */

  cs_gnum_t        *elt_num;         /* Global numbers associated with
                                        elements in local block
                                        (size: n_elts) */
  cs_lnum_t        *neighbor_index;  /* Start index of neighbors
                                        (size: n_elts + 1) */
  cs_gnum_t        *neighbor_num;    /* Global element neighbor numbers
                                        (size: neighbor_index[n_elts]) */

#if defined(HAVE_MPI)
  MPI_Comm  comm;                    /* Associated MPI communicator */
#endif

  /* Algorithm-related options */

  int  max_tree_depth;               /* Maximum search tree depth */
  int  leaf_threshold;               /* Maximum number of boxes which can
                                        be related to a leaf of the tree if
                                        level < max_tree_depth */
  float  max_box_ratio;              /* Stop adding levels to tree when
                                        (  n_linked_boxes
                                         > max_box_ratio*n_initial_boxes) */
  float  max_box_ratio_distrib;      /* In parallel, max_box_ratio for
                                        initial coarse tree used for
                                        distribution */

  _box_tree_stats_t  bt_stats;       /* Statistics associated with the
                                        box-trees used for search */

  /* Timings */

  double  cpu_time[2];   /* CPU time for tree construction and query */
  double  wtime[2];      /* Wall clock time for tree construction and query */

};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize box_tree statistics
 *
 * parameters:
 *   bts --> pointer to box tree statistics structure
 *---------------------------------------------------------------------------*/

static void
_init_bt_statistics(_box_tree_stats_t  *bts)
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}
/*----------------------------------------------------------------------------
 * Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * parameters:
 *   bts   <-> pointer to box tree statistics structure
 *   boxes <-- pointer to box tree structure
 *---------------------------------------------------------------------------*/

static void
_update_bt_statistics(_box_tree_stats_t     *bts,
                      const fvm_box_tree_t  *bt)
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = fvm_box_tree_get_stats(bt,
                               bts->depth,
                               bts->n_leaves,
                               bts->n_boxes,
                               bts->n_threshold_leaves,
                               bts->n_leaf_boxes,
                               bts->mem_used,
                               mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = CS_MAX(bts->mem_required[i], mem_required[i]);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Distribute bounding boxes over the ranks according to a Morton encoding
 * index. Try to get a well-balanced distribution and spatially coherent.
 *
 * parameters:
 *   n     <-> pointer to neighborhood management structure
 *   boxes <-> box set to redistribute
 *---------------------------------------------------------------------------*/

static void
_redistribute_boxes(fvm_neighborhood_t  *n,
                    fvm_box_set_t       *boxes)
{
  fvm_box_tree_t  *coarse_tree = NULL;
  fvm_box_distrib_t  *distrib = NULL;

  /* Sanity checks */

  assert(boxes != NULL);

  coarse_tree = fvm_box_tree_create(n->max_tree_depth,
                                    n->leaf_threshold,
                                    n->max_box_ratio_distrib);

  /* Build a tree and associate boxes */

  fvm_box_tree_set_boxes(coarse_tree,
                         boxes,
                         FVM_BOX_TREE_SYNC_LEVEL);

  _update_bt_statistics(&(n->bt_stats), coarse_tree);

  /*
    Compute an index based on Morton encoding to ensure a good distribution
    of bounding boxes among the ranks.
  */

  distrib = fvm_box_tree_get_distrib(coarse_tree, boxes);

  fvm_box_tree_destroy(&coarse_tree);

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  fvm_box_distrib_dump_statistics(distrib, n->comm);
#endif

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  fvm_box_set_redistribute(distrib, boxes);

  /* Delete intermediate structures */

  fvm_box_distrib_destroy(&distrib);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * using a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

static inline void
_gnum_shellsort(cs_lnum_t   l,
                cs_lnum_t   r,
                cs_gnum_t   a[])
{
  cs_lnum_t i, j, h;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1);

  /* Sort array */
  for (; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      cs_gnum_t   v = a[i];

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
 * Remove multiple element neighbor occurences
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

static void
_clean_neighbor_nums(fvm_neighborhood_t  *n)
{
  cs_lnum_t   i, j, start_id, end_id, saved_id, n_elts, n_neighbors;

  cs_lnum_t   n_count = 0;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;
  n_neighbors = n->neighbor_index[n_elts];

  /* Remove redundant elements */

  saved_id = n->neighbor_index[0];

  for (i = 0; i < n_elts; i++) {

    start_id = saved_id;
    end_id = n->neighbor_index[i+1];

    if (end_id - start_id > 0) {

      _gnum_shellsort(start_id, end_id, n->neighbor_num);

      n->neighbor_num[n_count++] = n->neighbor_num[start_id];

      for (j = start_id + 1; j < end_id; j++) {
        if (n->neighbor_num[j] != n->neighbor_num[j-1])
          n->neighbor_num[n_count++] = n->neighbor_num[j];
      }

    }

    saved_id = end_id;
    n->neighbor_index[i+1] = n_count;

  } /* End of loop on elements */

  if (n_count < n_neighbors)
    BFT_REALLOC(n->neighbor_num, n_count, cs_gnum_t);
}

/*----------------------------------------------------------------------------
 * Order a neighborhood based on the global numbering of elements.
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

static void
_order_neighborhood(fvm_neighborhood_t  *n)
{
  cs_lnum_t   i, j, k, order_id, shift, e_count;
  cs_lnum_t   n_elts, n_neighbors, n_elt_neighbors;
  cs_gnum_t   prev_num, cur_num;

  cs_lnum_t   *order = NULL, *old_index = NULL;
  cs_gnum_t   *old_e_num = NULL, *old_n_num = NULL;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;
  n_neighbors = n->neighbor_index[n_elts];

  BFT_MALLOC(order, n_elts, cs_lnum_t);
  BFT_MALLOC(old_e_num, n_elts, cs_gnum_t);
  BFT_MALLOC(old_index, n_elts + 1, cs_lnum_t);
  BFT_MALLOC(old_n_num, n_neighbors, cs_gnum_t);

  memcpy(old_e_num, n->elt_num, n_elts*sizeof(cs_gnum_t));
  memcpy(old_index, n->neighbor_index, (n_elts + 1)*sizeof(cs_lnum_t));
  memcpy(old_n_num, n->neighbor_num, n_neighbors*sizeof(cs_gnum_t));

  /* Order elt_num */

  cs_order_gnum_allocated(NULL, old_e_num, order, n_elts);

  /* Reshape according to the new ordering */

  /* Add first element */

  order_id = order[0];
  shift = 0;

  n->elt_num[0] = old_e_num[order_id];
  prev_num = n->elt_num[0];

  n->neighbor_index[0] = 0;
  n->neighbor_index[1] = old_index[order_id+1] - old_index[order_id];

  /* Loop on second-to last elements, merging data if an element has
     already appeared */

  for (i = 1, e_count = 1; i < n_elts; i++) {

    order_id = order[i];

    n_elt_neighbors = old_index[order_id+1] - old_index[order_id];

    shift = n->neighbor_index[i];

    cur_num = old_e_num[order_id];

    if (cur_num != prev_num) {
      n->elt_num[e_count] = cur_num;
      n->neighbor_index[e_count+1] = (  n->neighbor_index[e_count]
                                      + n_elt_neighbors);
      e_count += 1;
      prev_num = cur_num;
    }
    else
      n->neighbor_index[e_count] += n_elt_neighbors;

    for (j = old_index[order_id], k = 0; k < n_elt_neighbors; j++, k++)
      n->neighbor_num[shift + k] = old_n_num[j];

  } /* End of loop on elements */

  assert(n->neighbor_index[e_count] == n_neighbors);

  /* Free temporary memory */

  BFT_FREE(order);
  BFT_FREE(old_e_num);
  BFT_FREE(old_index);
  BFT_FREE(old_n_num);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Synchronize a neighborhood management structure and distribute the
 * resulting data over ranks by block.
 *
 * parameters:
 *   n        <-- pointer to neighborhood management structure
 *   n_g_elts <-- global number of elements
 *----------------------------------------------------------------------------*/

static void
_sync_by_block(fvm_neighborhood_t  *n,
               cs_gnum_t            n_g_elts)
{
  assert(n != NULL);

  if (n_g_elts == 0 || n->comm == MPI_COMM_NULL)
    return;

  /* Initialization */

  int  rank_id, n_ranks;

  MPI_Comm_rank(n->comm, &rank_id);
  MPI_Comm_size(n->comm, &n_ranks);

  cs_block_dist_info_t bi = cs_block_dist_compute_sizes(rank_id,
                                                        n_ranks,
                                                        0,
                                                        0,
                                                        n_g_elts);

  cs_all_to_all_t *d
    = cs_all_to_all_create_from_block(n->n_elts,
                                      0, /* flags */
                                      n->elt_num,
                                      bi,
                                      n->comm);

  cs_gnum_t *send_meta;
  BFT_MALLOC(send_meta, n->n_elts*2, cs_gnum_t);

  for (cs_lnum_t i = 0; i < n->n_elts; i++) {
    send_meta[i*2] = n->elt_num[i];
    send_meta[i*2+1] = n->neighbor_index[i+1] - n->neighbor_index[i];
  }

  cs_gnum_t *recv_meta
    = cs_all_to_all_copy_array(d,
                               CS_GNUM_TYPE,
                               2,
                               false, /* reverse */
                               send_meta,
                               NULL);

  BFT_FREE(send_meta);

  cs_lnum_t n_recv = cs_all_to_all_n_elts_dest(d);

  cs_gnum_t *recv_gnum;
  BFT_MALLOC(recv_gnum, n_recv, cs_gnum_t);
  cs_lnum_t *recv_index;
  BFT_MALLOC(recv_index, n_recv+1, cs_lnum_t);

  recv_index[0] = 0;
  for (cs_lnum_t i = 0; i < n_recv; i++) {
    recv_gnum[i] = recv_meta[i*2];
    recv_index[i+1] = recv_index[i] + recv_meta[i*2+1];
  }

  BFT_FREE(recv_meta);

  cs_gnum_t *recv_neighbor
    = cs_all_to_all_copy_indexed(d,
                                 CS_GNUM_TYPE,
                                 false, /* reverse */
                                 n->neighbor_index,
                                 n->neighbor_num,
                                 recv_index,
                                 NULL);

  /* Build arrays corresponding to block distribution of neighborhood */

  n->n_elts = bi.gnum_range[1] - bi.gnum_range[0];

  BFT_FREE(n->elt_num);
  BFT_FREE(n->neighbor_index);
  BFT_FREE(n->neighbor_num);

  BFT_MALLOC(n->elt_num, n->n_elts, cs_gnum_t);
  BFT_MALLOC(n->neighbor_index, n->n_elts + 1, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n->n_elts; i++) {
    n->elt_num[i] = bi.gnum_range[0] + i;
    n->neighbor_index[i] = 0;
  }
  n->neighbor_index[n->n_elts] = 0;

  /* Count element neighbors in block distribution */

  for (cs_lnum_t i = 0; i < n_recv; i++) {
    cs_lnum_t elt_id = recv_gnum[i] - bi.gnum_range[0];
    n->neighbor_index[elt_id + 1] += recv_index[i+1] - recv_index[i];
  }

  /* Transform element neighbor count to index */

  n->neighbor_index[0] = 0;
  for (cs_lnum_t i = 0; i < n->n_elts; i++)
    n->neighbor_index[i+1] += n->neighbor_index[i];

  BFT_MALLOC(n->neighbor_num, n->neighbor_index[n->n_elts], cs_gnum_t);

  /* Fill element neighbors in block distribution */

  cs_lnum_t *counter;
  BFT_MALLOC(counter, n->n_elts, cs_lnum_t);

  for (cs_lnum_t i = 0; i < n->n_elts; i++)
    counter[i] = 0;

  for (cs_lnum_t i = 0; i < n_recv; i++) {

    cs_lnum_t elt_id = recv_gnum[i] - bi.gnum_range[0];

    cs_lnum_t s_id = recv_index[i];
    cs_lnum_t e_id = recv_index[i+1];

    cs_lnum_t shift = n->neighbor_index[elt_id] + counter[elt_id];

    for (cs_lnum_t j = s_id; j < e_id; j++)
      n->neighbor_num[shift++] = recv_neighbor[j];

    counter[elt_id] += e_id - s_id;

  } /* End of loop on ranks */

  BFT_FREE(recv_gnum);
  BFT_FREE(counter);
  BFT_FREE(recv_neighbor);
  BFT_FREE(recv_index);

  cs_all_to_all_destroy(&d);

  /* Remove redundant data */

  _clean_neighbor_nums(n);
}

#endif /* defined(HAVE_MPI) */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a neighborhood_t structure and initialize it.
 *
 * parameters:
 *   comm  <-- associated MPI communicator
 *
 * returns:
 *   pointer to an empty fvm_box_tree_t structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
fvm_neighborhood_t *
fvm_neighborhood_create(MPI_Comm  comm)
#else
fvm_neighborhood_t *
fvm_neighborhood_create(void)
#endif
{
  double  w_start, w_end, cpu_start, cpu_end;

  fvm_neighborhood_t *n = NULL;

  /* Timer start */

  w_start = cs_timer_wtime();
  cpu_start = cs_timer_cpu_time();

  /* Allocate and initialize */

  BFT_MALLOC(n, 1, fvm_neighborhood_t);

  n->n_elts = 0;
  n->elt_num = NULL;
  n->neighbor_index = NULL;
  n->neighbor_num = NULL;

#if defined(HAVE_MPI)
  n->comm = comm;
#endif

  /* Algorithm options */

  n->max_tree_depth = 30; /* defaults */
  n->leaf_threshold = 30;
  n->max_box_ratio = 10.0;
  n->max_box_ratio_distrib = 6.0;

  _init_bt_statistics(&(n->bt_stats));

  /* Timer end */

  w_end = cs_timer_wtime();
  cpu_end = cs_timer_cpu_time();

  n->cpu_time[0] = cpu_end - cpu_start;  /* build time */
  n->wtime[0] = w_end - w_start;
  n->cpu_time[1] = 0.0;                  /* query */
  n->wtime[1] = 0.0;

  /* Return structure */

  return n;
}

/*----------------------------------------------------------------------------
 * Destroy a neighborhood_t structure.
 *
 * parameters:
 *   n <-> pointer to pointer to fvm_neighborhood_t structure to destroy.
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_destroy(fvm_neighborhood_t  **n)
{
  if (n != NULL) {
    fvm_neighborhood_t *_n = *n;
    if (_n != NULL) {
      if (_n->elt_num != NULL)
        BFT_FREE(_n->elt_num);
      if (_n->neighbor_index != NULL)
        BFT_FREE(_n->neighbor_index);
      if (_n->neighbor_num != NULL)
        BFT_FREE(_n->neighbor_num);
    }
    BFT_FREE(*n);
  }
}

/*----------------------------------------------------------------------------
 * Set non-default algorithm parameters for neighborhood management structure.
 *
 * parameters:
 *   n                     <-> pointer to neighborhood management structure
 *   max_tree_depth        <-- maximum search tree depth
 *   leaf_threshold        <-- maximum number of boxes which can be related to
 *                             a leaf of the tree if level < max_tree_depth
 *   max_box_ratio         <-- stop adding levels to tree when
 *                             (n_linked_boxes > max_box_ratio*n_init_boxes)
 *   max_box_ratio_distrib <-- maximum box ratio when computing coarse
 *                             tree prior to parallel distribution
 *---------------------------------------------------------------------------*/

void
fvm_neighborhood_set_options(fvm_neighborhood_t  *n,
                             int                  max_tree_depth,
                             int                  leaf_threshold,
                             float                max_box_ratio,
                             float                max_box_ratio_distrib)
{
  if (n == NULL)
    return;

  n->max_tree_depth = max_tree_depth;
  n->leaf_threshold = leaf_threshold;
  n->max_box_ratio = max_box_ratio;
  n->max_box_ratio_distrib = max_box_ratio_distrib;
}

/*----------------------------------------------------------------------------
 * Retrieve pointers to of arrays from a neighborhood_t structure.
 *
 * Arrays remain the property of the neighborhood_t structure, and must not
 * be modified by the caller.
 *
 * parameters:
 *   n              <-> pointer to fvm_neighborhood_t structure.
 *   n_elts         --> number of elements with neighbors in block
 *                      associated with local rank
 *   elt_num        --> global element numbers in local block (size: n_elts)
 *   neighbor_index --> start index of neighbors (size: n_elts + 1)
 *   neighbor_num   --> global element neighbor numbers
 *                      (size: neighbor_index[n_elts])
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_get_data(const fvm_neighborhood_t         *n,
                          cs_lnum_t                        *n_elts,
                          cs_gnum_t                 **const elt_num,
                          cs_lnum_t                 **const neighbor_index,
                          cs_gnum_t                 **const neighbor_num)
{
  if (n != NULL) {

    if (n_elts != NULL)
      *n_elts = n->n_elts;

    if (elt_num != NULL)
      *elt_num = n->elt_num;

    if (neighbor_index != NULL)
      *neighbor_index = n->neighbor_index;

    if (neighbor_num != NULL)
      *neighbor_num = n->neighbor_num;
  }
}

/*----------------------------------------------------------------------------
 * Transfer ownership of arrays from a neighborhood_t structure to the
 * calling program.
 *
 * Arrays that are transferred are removed from the structure, so its
 * use after calling this function is limited to querying timing information
 * until a new neighborhood is computed.
 *
 * parameters:
 *   n              <-> pointer to fvm_neighborhood_t structure.
 *   n_elts         --> number of elements with neighbors in block
 *                      associated with local rank
 *   elt_num        --> global element numbers in local block (size: n_elts)
 *   neighbor_index --> start index of neighbors (size: n_elts + 1)
 *   neighbor_num   --> global element neighbor numbers
 *                      (size: neighbor_index[n_elts])
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_transfer_data(fvm_neighborhood_t   *n,
                               cs_lnum_t            *n_elts,
                               cs_gnum_t           **elt_num,
                               cs_lnum_t           **neighbor_index,
                               cs_gnum_t           **neighbor_num)
{
  if (n != NULL) {

    if (n_elts != NULL)
      *n_elts = n->n_elts;

    if (elt_num != NULL) {
      *elt_num = n->elt_num;
      n->elt_num = NULL;
    }
    if (neighbor_index != NULL) {
      *neighbor_index = n->neighbor_index;
      n->neighbor_index = NULL;
    }
    if (neighbor_num != NULL) {
      *neighbor_num = n->neighbor_num;
      n->neighbor_num = NULL;
    }
  }
}

/*----------------------------------------------------------------------------
 * Determine intersecting boxes.
 *
 * Box global numbers and extents may be either copied for the structure's
 * internal use from the caller, or transferred to the neighborhood management
 * structure: both the box_gnum and extents arguments have an "assigned"
 * variant, in which case a pointer to a pointer is provided, and the
 * argument's property is transferred to the neighborhod management structure.
 * The unused variant of an argument should be set to NULL.
 *
 * Boxes may be distributed among processors, so their intersections are
 * determined using a block distribution, and defined using their
 * global numbers.
 *
 * All global numbers appearing in box_gnum[] will have a matching entry in
 * the neighborhod structure. To remove global numbers entries with no
 * neighbors from the structure, the fvm_neighborhood_prune() function
 * may be used.
 *
 * parameters:
 *   n                 <-> pointer to neighborhood management structure
 *   dim               <-- spatial dimension
 *   n_boxes           <-- local number of boxes
 *   box_gnum          <-- global numbering of boxes
 *   extents           <-- coordinate extents (size: n_boxes*dim*2, as
 *                         xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   box_gnum_assigned <-> as box_gnum, ownership transferred (NULL on return)
 *   extents_assigned  <-> as extents, ownership transferred (NULL on return)
 *---------------------------------------------------------------------------*/

void
fvm_neighborhood_by_boxes(fvm_neighborhood_t  *n,
                          int                  dim,
                          cs_lnum_t            n_boxes,
                          const cs_gnum_t     *box_gnum,
                          const cs_coord_t    *extents,
                          cs_gnum_t          **box_gnum_assigned,
                          cs_coord_t         **extents_assigned)
{
  double  clock_start, clock_end, cpu_start, cpu_end;

  fvm_box_tree_t  *bt = NULL;
  fvm_box_set_t  *boxes = NULL;

  const cs_gnum_t   *_box_gnum = box_gnum;
  const cs_coord_t  *_extents = extents;

  int  n_ranks = 1;

  clock_start = cs_timer_wtime();
  cpu_start = cs_timer_cpu_time();

  /* Transfer data if necessary */

  if (box_gnum_assigned != NULL)
    _box_gnum = *box_gnum_assigned;
  if (extents_assigned != NULL)
    _extents = *extents_assigned;

  /* Reset structure if necessary */

  n->n_elts = 0;
  if (n->elt_num != NULL)
    BFT_FREE(n->elt_num);
  if (n->neighbor_index != NULL)
    BFT_FREE(n->neighbor_index);
  if (n->neighbor_num != NULL)
    BFT_FREE(n->neighbor_num);

  /* Allocate fvm_box_set_t structure and initialize it */

#if defined(HAVE_MPI)

  if (n->comm != MPI_COMM_NULL)
    MPI_Comm_size(n->comm, &n_ranks);

  boxes = fvm_box_set_create(dim,
                             1,  /* normalize */
                             1,  /* allow_projection */
                             n_boxes,
                             _box_gnum,
                             _extents,
                             n->comm);

  if (n_ranks > 1)
    _redistribute_boxes(n, boxes);

#else

  boxes = fvm_box_set_create(dim,
                             1,  /* normalize */
                             1,  /* allow_projection */
                             n_boxes,
                             _box_gnum,
                             _extents);

#endif

  /* Free transferred data if applicable */

  if (box_gnum_assigned != NULL) {
    _box_gnum = NULL;
    BFT_FREE(*box_gnum_assigned);
  }

  if (extents_assigned != NULL) {
    _extents = NULL;
    BFT_FREE(*extents_assigned);
  }

  /* Build a tree structure and use it to order bounding boxes */

  /* Create and initialize a box tree structure */

  bt = fvm_box_tree_create(n->max_tree_depth,
                           n->leaf_threshold,
                           n->max_box_ratio);

  /* Build a tree and put bounding boxes */

  fvm_box_tree_set_boxes(bt,
                         boxes,
                         FVM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics((&n->bt_stats), bt);

  /* Update construction times. */

  clock_end = cs_timer_wtime();
  cpu_end = cs_timer_cpu_time();

  n->cpu_time[0] = cpu_end - cpu_start;
  n->wtime[0] = clock_end - clock_start;

  clock_start = clock_end;
  cpu_start = cpu_end;

  /* Allocate structure to store intersections between boxes */

  n->n_elts = fvm_box_set_get_size(boxes);

  BFT_MALLOC(n->elt_num, n->n_elts, cs_gnum_t);
  if (n->n_elts > 0)
    memcpy(n->elt_num,
           fvm_box_set_get_g_num(boxes),
           n->n_elts*sizeof(cs_gnum_t));

  fvm_box_tree_get_intersects(bt,
                              boxes,
                              &(n->neighbor_index),
                              &(n->neighbor_num));

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  fvm_box_tree_dump(bt);
  fvm_box_set_dump(boxes, 1);
#endif

  /* Destroy the associated box tree */

  fvm_box_tree_destroy(&bt);

  /* Compact intersections list, delete redundancies and order intersections */

  _order_neighborhood(n);

#if defined(HAVE_MPI)

  /* Synchronize list of intersections for each element of the list
     and distribute it by block over the ranks */

  if (n_ranks > 1)
    _sync_by_block(n, fvm_box_set_get_global_size(boxes));

#endif /* HAVE_MPI */

  /* Destroy the box set structures */

  fvm_box_set_destroy(&boxes);

  _clean_neighbor_nums(n);

  /* Update query times. */

  clock_end = cs_timer_wtime();
  cpu_end = cs_timer_cpu_time();

  n->cpu_time[1] = cpu_end - cpu_start;
  n->wtime[1] = clock_end - clock_start;
}

/*----------------------------------------------------------------------------
 * Prune a neighborhood (remove entries with no neighbors).
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_prune(fvm_neighborhood_t  *n)
{
  cs_lnum_t   i, start_id, end_id, saved_id, n_elts;

  cs_lnum_t   e_count = 0;

  assert(n != NULL);

  if (n->n_elts == 0)
    return;

  n_elts = n->n_elts;

  /* Remove elements with no neighbors */

  saved_id = n->neighbor_index[0];

  for (i = 0; i < n_elts; i++) {

    start_id = saved_id;
    end_id = n->neighbor_index[i+1];

    if (end_id - start_id > 0) {

      n->elt_num[e_count] = n->elt_num[i];

      saved_id = end_id;
      n->neighbor_index[e_count+1] = end_id;

      e_count += 1;

    }

  }

  if (e_count < n_elts) {
    n->n_elts = e_count;
    BFT_REALLOC(n->elt_num, e_count, cs_gnum_t);
    BFT_REALLOC(n->neighbor_index, e_count + 1, cs_lnum_t);
  }
}

/*----------------------------------------------------------------------------
 * Get global statistics relative to the search structures used
 * by fvm_neighborhood_by_boxes().
 *
 * All fields returned are optional: if their argument is set to NULL,
 * the corresponding information will not be returned.
 *
 * For each field not set to NULL, 3 values are always returned:
 * the mean on all ranks (rounded to the closest integer), the minimum,
 * and the maximum value respectively.
 *
 * In serial mode, the mean, minimum, and maximum will be identical for most
 * fields, but all 3 values are returned nonetheless.
 *
 * Note that the final memory use is only relative to the final search
 * structure, and the comparison of the total (or average) with the minima
 * and maxima may give an indication on load balancing.

 * The mem_required field only takes into account the theoretical maximum
 * memory footprint of the main search structure during its construction phase,
 * and of that of temporary structures before load balancing in parallel mode,
 * but does not include minor work arrays or buffers used during the algorithm.
 *
 * Neither of the 2 memory fields include the footprint of the arrays
 * containing the query results.
 *
 * parameters:
 *   n                  <-- pointer to neighborhood management structure
 *   dim                --> layout dimension (3, 2, or 1)
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_final          --> theoretical memory for final search structure
 *   mem_required       --> theoretical maximum memory for main structures
 *                          used during the algorithm
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
fvm_neighborhood_get_box_stats(const fvm_neighborhood_t  *n,
                               int                        depth[3],
                               cs_lnum_t                  n_leaves[3],
                               cs_lnum_t                  n_boxes[3],
                               cs_lnum_t                  n_threshold_leaves[3],
                               cs_lnum_t                  n_leaf_boxes[3],
                               size_t                     mem_final[3],
                               size_t                     mem_required[3])
{
  size_t i;

  if (n == NULL)
    return 0;

  for (i = 0; i < 3; i++) {

    if (depth != NULL)
      depth[i] = n->bt_stats.depth[i];

    if (n_leaves != NULL)
      n_leaves[i] = n->bt_stats.n_leaves[i];

    if (n_boxes != NULL)
      n_boxes[i] = n->bt_stats.n_boxes[i];

    if (n_threshold_leaves != NULL)
      n_threshold_leaves[i] = n->bt_stats.n_threshold_leaves[i];

    if (n_leaf_boxes != NULL)
      n_leaf_boxes[i] = n->bt_stats.n_leaf_boxes[i];

    if (mem_final != NULL)
      mem_final[i] = n->bt_stats.mem_used[i];

    if (mem_required != NULL)
      mem_required[i] = n->bt_stats.mem_required[i];
  }
  return n->bt_stats.dim;
}

/*----------------------------------------------------------------------------
 * Return timing information.
 *
 * parameters:
 *   n              <-- pointer to neighborhood management structure
 *   build_wtime    --> initialization Wall-clock time (or NULL)
 *   build_cpu_time --> initialization CPU time (or NULL)
 *   query_wtime    --> query Wall-clock time (or NULL)
 *   query_cpu_time --> query CPU time (or NULL)
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_get_times(const fvm_neighborhood_t  *n,
                           double                    *build_wtime,
                           double                    *build_cpu_time,
                           double                    *query_wtime,
                           double                    *query_cpu_time)
{
  if (n == NULL)
    return;

  if (build_wtime != NULL)
    *build_wtime = n->wtime[0];
  if (build_cpu_time != NULL)
    *build_cpu_time = n->cpu_time[0];

  if (query_wtime != NULL)
    *query_wtime = n->wtime[1];
  if (query_cpu_time != NULL)
    *query_cpu_time = n->cpu_time[1];
}

/*----------------------------------------------------------------------------
 * Dump a neighborhood management structure.
 *
 * parameters:
 *   n <-- pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_dump(const fvm_neighborhood_t  *n)
{
  cs_lnum_t   i, j;

  bft_printf("\n"
             "Neighborhood information: %p\n\n", (const void *)n);

  if (n == NULL)
    return;

  bft_printf("number of elements: %10d\n"
             "list size:          %10d\n\n",
             (int)(n->n_elts), (int)(n->neighbor_index[n->n_elts]));

  bft_printf("max tree depth:     %d\n"
             "leaf threshold:     %d\n"
             "max box ratio       %f\n\n",
             n->max_tree_depth, n->leaf_threshold, n->max_box_ratio);

#if defined(HAVE_MPI)
  if (n->comm != MPI_COMM_NULL)
    bft_printf("\n"
               "Associated MPI communicator: %ld\n",
               (long)(n->comm));
#endif

  bft_printf("CPU time:           %f construction, %f query\n"
             "Wall-clock time:    %f construction, %f query\n\n",
             n->cpu_time[0], n->cpu_time[1],
             n->wtime[0], n->wtime[1]);

  for (i = 0; i < n->n_elts; i++) {

    int n_neighbors = (n->neighbor_index[i+1] - n->neighbor_index[i]);

    bft_printf("global num.: %10llu | n_neighbors : %3d |",
               (unsigned long long)(n->elt_num[i]), n_neighbors);

    for (j = n->neighbor_index[i]; j < n->neighbor_index[i+1]; j++)
      bft_printf("  %10llu ", (unsigned long long)(n->neighbor_num[j]));
    bft_printf("\n");

  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
