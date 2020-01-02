#ifndef __FVM_NEIGHBORHOOD_H__
#define __FVM_NEIGHBORHOOD_H__

/*============================================================================
 * Search octrees and quadtrees of boxes.
 *============================================================================*/

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvm_neighborhood_t fvm_neighborhood_t;

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
fvm_neighborhood_create(MPI_Comm  comm);

#else

fvm_neighborhood_t *
fvm_neighborhood_create(void);

#endif

/*----------------------------------------------------------------------------
 * Destroy a neighborhood_t structure.
 *
 * parameters:
 *   n <-> pointer to pointer to fvm_neighborhood_t structure to destroy.
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_destroy(fvm_neighborhood_t  **n);

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
                             float                max_box_ratio_distrib);

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
                          cs_gnum_t                 **const neighbor_num);

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
                               cs_gnum_t           **neighbor_num);

/*----------------------------------------------------------------------------
 * Determine intersecting boxes.
 *
 * Box global numbers and extents may be either copied for the structure's
 * internal use from the caller, or transferred to the neighborhood management
 * structure: both the box_gnum and extents arguments have an "assigned"
 * variant, in which cas a pointer to a pointer is provided, and the
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
 *   box_gnum_transfer <-> as box_gnum, ownership transferred (NULL on return)
 *   extents_transfer  <-> as extents, ownership transferred (NULL on return)
 *   comm       <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
fvm_neighborhood_by_boxes(fvm_neighborhood_t  *n,
                          int                  dim,
                          cs_lnum_t            n_boxes,
                          const cs_gnum_t     *box_gnum,
                          const cs_coord_t    *extents,
                          cs_gnum_t          **box_gnum_assigned,
                          cs_coord_t         **extents_assigned);

/*----------------------------------------------------------------------------
 * Prune a neighborhood (remove entries with no neighbors).
 *
 * parameters:
 *   n <-> pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_prune(fvm_neighborhood_t  *n);

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
                               size_t                     mem_required[3]);

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
                           double                    *query_cpu_time);

/*----------------------------------------------------------------------------
 * Dump a neighborhood management structure.
 *
 * parameters:
 *   n <-- pointer to neighborhood management structure
 *----------------------------------------------------------------------------*/

void
fvm_neighborhood_dump(const fvm_neighborhood_t  *n);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_NEIGHBORHOOD_H__ */
