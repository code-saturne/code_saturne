#ifndef __FVM_BOX_H__
#define __FVM_BOX_H__

/*============================================================================
 * Handle boxes aligned with Cartesian axes.
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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"
#include "fvm_morton.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

/* Collection of boxes */

typedef struct _fvm_box_set_t fvm_box_set_t;

/* Distribution on octree or quadtree */

typedef struct _fvm_box_distrib_t fvm_box_distrib_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a set of boxes and initialize it.
 *
 * parameters:
 *   dim              <-- spatial dimension
 *   normalize        <-- 1 if boxes are to be normalized, 0 otherwize
 *   allow_projection <-- if 1, project to lower dimension if all boxes
 *                        are cut by the median plane of the set.
 *   n_boxes          <-- number of elements to create
 *   box_gnum         <-- global numbering of boxes
 *   extents          <-- coordinate extents (size: n_boxes*dim*2, as
 *                        xmin1, ymin1, .. xmax1, ymax1, ..., xmin2, ...)
 *   comm             <-- associated MPI communicator
 *
 * returns:
 *   a new allocated pointer to a fvm_box_set_t structure.
 *---------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
fvm_box_set_t *
fvm_box_set_create(int                dim,
                   int                normalize,
                   int                allow_projection,
                   cs_lnum_t          n_boxes,
                   const cs_gnum_t   *box_gnum,
                   const cs_coord_t  *box_extents,
                   MPI_Comm           comm);
#else
fvm_box_set_t *
fvm_box_set_create(int                dim,
                   int                normalize,
                   int                allow_projection,
                   cs_lnum_t          n_boxes,
                   const cs_gnum_t   *box_gnum,
                   const cs_coord_t  *box_extents);

#endif

/*----------------------------------------------------------------------------
 * Destroy a fvm_box_set_t structure.
 *
 * parameters:
 *   boxes <-> pointer to pointer to the fvm_box_set_t structure to delete
 *----------------------------------------------------------------------------*/

void
fvm_box_set_destroy(fvm_box_set_t  **boxes);

/*----------------------------------------------------------------------------
 * Return the dimension associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   associated spatial dimension
 *---------------------------------------------------------------------------*/

int
fvm_box_set_get_dim(const fvm_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return the local number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

cs_lnum_t
fvm_box_set_get_size(const fvm_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return the global number of boxes in a set.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   local number of boxes
 *---------------------------------------------------------------------------*/

cs_gnum_t
fvm_box_set_get_global_size(const fvm_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return extents associated with a set of boxes.
 *
 * The extents array is organized in the following fashion:
 * {x_min_0, y_min_0, ..., x_max_0, y_max_0, ...
 *  x_min_n, y_min_n, ..., x_max_n, y_max_n, ...}
 *
 * Its size is thus: n_boxes * dim * 2.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to extents array
 *---------------------------------------------------------------------------*/

const cs_coord_t *
fvm_box_set_get_extents(fvm_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Return global numbers associated with a set of boxes.
 *
 * parameters:
 *   boxes <-- pointer to set of boxes
 *
 * returns:
 *   pointer to global box numbers array
 *---------------------------------------------------------------------------*/

const cs_gnum_t *
fvm_box_set_get_g_num(fvm_box_set_t  *boxes);

/*----------------------------------------------------------------------------
 * Build a Morton_index to get a well-balanced distribution of the boxes.
 *
 * parameters:
 *  boxes      <-- pointer to associated fvm_box_set_t structure
 *  distrib    <-> pointer to a fvm_box_distrib_t structure
 *  n_leaves   <-- number of leaves with weight > 0
 *  leaf_codes <-- Morton code for each leaf
 *  weight     <-- number of boxes related to each leaf
 *---------------------------------------------------------------------------*/

void
fvm_box_set_build_morton_index(const fvm_box_set_t  *boxes,
                               fvm_box_distrib_t    *distrib,
                               cs_lnum_t             n_leaves,
                               fvm_morton_code_t    *leaf_codes,
                               cs_lnum_t            *weight);

/*----------------------------------------------------------------------------
 * Redistribute boxes over the ranks according to the Morton index to
 * assume a better balanced distribution of the boxes.
 *
 * parameters:
 *  box_distrib <--  data structure on box distribution
 *  box_set     <->  pointer to the structure to redistribute
 *---------------------------------------------------------------------------*/

void
fvm_box_set_redistribute(const fvm_box_distrib_t  *box_distrib,
                         fvm_box_set_t            *boxes);

/*----------------------------------------------------------------------------
 * Dump a fvm_box_set_t structure.
 *
 * parameters:
 *   box_set   <-- pointer to the fvm_box_t structure
 *   verbosity <-- verbosity level (0 or 1)
 *----------------------------------------------------------------------------*/

void
fvm_box_set_dump(const fvm_box_set_t  *boxes,
                 int                   verbosity);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a fvm_box_distrib_t structure.
 *
 * parameters:
 *   n_boxes   <-- number of boxes
 *   n_g_boxes <-- global number of boxes
 *   max_level <-- max level reached locally in the related tree
 *   comm      <-- MPI comm. on which distribution takes place
 *
 * returns:
 *   a pointer to a new allocated fvm_box_distrib_t structure.
 *---------------------------------------------------------------------------*/

fvm_box_distrib_t *
fvm_box_distrib_create(cs_lnum_t  n_boxes,
                       cs_gnum_t  n_g_boxes,
                       int        max_level,
                       MPI_Comm   comm);

/*----------------------------------------------------------------------------
 * Destroy a fvm_box_distrib_t structure.
 *
 * parameters:
 *   distrib <-> pointer to pointer to the structure to destroy
 *---------------------------------------------------------------------------*/

void
fvm_box_distrib_destroy(fvm_box_distrib_t  **distrib);

/*----------------------------------------------------------------------------
 * Delete redundancies in box distribution
 *
 * parameters:
 *   distrib <-> pointer to the fvm_box_distrib_t structure
 *---------------------------------------------------------------------------*/

void
fvm_box_distrib_clean(fvm_box_distrib_t  *distrib);

/*----------------------------------------------------------------------------
 * Display a histogramm on leaves associated to the boxes and
 * several other pieces of information (min, max, ...)
 *
 * parameters:
 *   distrib <-- pointer to the fvm_box_distrib_t structure
 *   comm    <-- associated MPI communicator
 *---------------------------------------------------------------------------*/

void
fvm_box_distrib_dump_statistics(const fvm_box_distrib_t  *distrib,
                                MPI_Comm                  comm);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __FVM_BOX_H__ */
