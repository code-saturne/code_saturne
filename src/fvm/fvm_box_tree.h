#ifndef __FVM_BOX_TREE_H__
#define __FVM_BOX_TREE_H__

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_box.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _fvm_box_tree_t fvm_box_tree_t;

typedef enum {

  FVM_BOX_TREE_ASYNC_LEVEL,  /* Boxes are placed according to tree parameters,
                                and potentially at different levels */
  FVM_BOX_TREE_SYNC_LEVEL    /* All boxes are placed for all ranks at the
                                same level */

} fvm_box_tree_sync_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a fvm_box_tree_t structure and initialize it.
 *
 * parameters:
 *  max_level     <-- max possible level
 *  threshold     <-- max number of  boxes linked to an octant if
 *                    max_level is not reached
 *  max_box_ratio <-- max n_linked_boxes / n_boxes ratio
 *
 * returns:
 *   pointer to an empty fvm_box_tree_t structure.
 *----------------------------------------------------------------------------*/

fvm_box_tree_t *
fvm_box_tree_create(int    max_level,
                    int    threshold,
                    float  max_box_ratio);

/*----------------------------------------------------------------------------
 * Destroy a fvm_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to fvm_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

void
fvm_box_tree_destroy(fvm_box_tree_t  **bt);

/*----------------------------------------------------------------------------
 * Get the deepest level allowed by the tree structure.
 *
 * parameters:
 *   bt <-- pointer to fvm_box_tree_t structure.
 *
 * returns:
 *   deepest allowed level of the tree
 *----------------------------------------------------------------------------*/

int
fvm_box_tree_get_max_level(const fvm_box_tree_t  *bt);

/*----------------------------------------------------------------------------
 * Assign a set of boxes to an empty fvm_box_tree_t structure.
 *
 * The box tree structure must have been created using to fvm_tree_create().
 *
 * The depth of the tree is adjusted so that a maximum of max_n_elts boxes
 * will be assigned to each leaf, unless this would require going beyond
 * the tree's maximum level.
 *
 * If max_level = -1, the highest level reachable is FVM_TREE_MAX_LEVEL but
 * there is no defined target level.
 *
 * parameters:
 *   bt         <-> pointer to fvm_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   build_type <-- layout variant for building the tree structure
 *----------------------------------------------------------------------------*/

void
fvm_box_tree_set_boxes(fvm_box_tree_t       *bt,
                       const fvm_box_set_t  *boxes,
                       fvm_box_tree_sync_t   build_type);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute an index based on Morton encoding to ensure a good distribution
 * of boxes among the participating ranks.
 *
 * parameters:
 *   bt         <-> pointer to fvm_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *
 * returns:
 *   pointer to newly created fvm_box_distrib_t structure.
 *----------------------------------------------------------------------------*/

fvm_box_distrib_t *
fvm_box_tree_get_distrib(fvm_box_tree_t        *bt,
                         const fvm_box_set_t   *boxes);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Build an indexed list on boxes to list intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   boxes     <-- pointer to a associated box set
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxes
 *----------------------------------------------------------------------------*/

void
fvm_box_tree_get_intersects(fvm_box_tree_t       *bt,
                            const fvm_box_set_t  *boxes,
                            cs_lnum_t            *box_index[],
                            cs_gnum_t            *box_g_num[]);

/*----------------------------------------------------------------------------
 * Get global box tree statistics.
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
 * Note that the theoretical memory use includes that of the associated
 * box set.
 *
 * parameters:
 *   bt                 <-- pointer to box tree structure
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_used           --> theoretical used memory
 *   mem_allocated      --> theoretical allocated memory
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
fvm_box_tree_get_stats(const fvm_box_tree_t  *bt,
                       int                    depth[3],
                       cs_lnum_t              n_leaves[3],
                       cs_lnum_t              n_boxes[3],
                       cs_lnum_t              n_threshold_leaves[3],
                       cs_lnum_t              n_leaf_boxes[3],
                       size_t                 mem_used[3],
                       size_t                 mem_allocated[3]);

/*----------------------------------------------------------------------------
 * Display local statistics about a fvm_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
fvm_box_tree_dump_statistics(const fvm_box_tree_t  *bt);

/*----------------------------------------------------------------------------
 * Dump an fvm_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
fvm_box_tree_dump(fvm_box_tree_t  *bt);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_BOX_TREE_H__ */
