#ifndef __CS_RENUMBER_H__
#define __CS_RENUMBER_H__

/*============================================================================
 * Optional mesh renumbering
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Renumbering algorithms */

typedef enum {

  CS_RENUMBER_CELLS_SCOTCH_PART,     /* SCOTCH partitioning */
  CS_RENUMBER_CELLS_SCOTCH_ORDER,    /* SCOTCH ordering */
  CS_RENUMBER_CELLS_METIS_PART,      /* METIS partitioning */
  CS_RENUMBER_CELLS_METIS_ORDER,     /* METIS ordering */
  CS_RENUMBER_CELLS_MORTON,          /* Morton space filling curve */
  CS_RENUMBER_CELLS_HILBERT,         /* Hilbert space filling curve */
  CS_RENUMBER_CELLS_RCM,             /* Reverse Cuthill-McKee */
  CS_RENUMBER_CELLS_NONE             /* No cells renumbering */

} cs_renumber_cells_type_t;

typedef enum {

  CS_RENUMBER_I_FACES_BLOCK,         /* No shared cell in block */
  CS_RENUMBER_I_FACES_MULTIPASS,     /* Use multipass face numbering */
  CS_RENUMBER_I_FACES_SIMD,          /* Renumber for vector (SIMD) operations */
  CS_RENUMBER_I_FACES_NONE           /* No interior face numbering */

} cs_renumber_i_faces_type_t;

typedef enum {

  CS_RENUMBER_B_FACES_THREAD,        /* No cell shared between threads */
  CS_RENUMBER_B_FACES_SIMD,          /* Renumber for vector (SIMD) operations */
  CS_RENUMBER_B_FACES_NONE           /* No interior face numbering */

} cs_renumber_b_faces_type_t;

typedef enum {

  CS_RENUMBER_VERTICES_NONE          /* No vertex numbering */

} cs_renumber_vertices_type_t;

/* Ordering options for adjacency arrays */

typedef enum {

  CS_RENUMBER_ADJACENT_LOW,         /* Lowest adjacent id first */
  CS_RENUMBER_ADJACENT_HIGH         /* Highest adjacent id first */

} cs_renumber_ordering_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set the target number of threads for mesh renumbering.
 *
 * By default, the target number of threads is set to cs_glob_n_threads,
 * but the value may be forced using this function. This is mainly useful
 * for testing purposes.
 *
 * parameters:
 *   n_threads <-- target number of threads for mesh numbering
 *----------------------------------------------------------------------------*/

void
cs_renumber_set_n_threads(int  n_threads);

/*----------------------------------------------------------------------------
 * Return the target number of threads for mesh renumbering.
 *
 * returns:
 *   the target number of threads for mesh numbering
 *----------------------------------------------------------------------------*/

int
cs_renumber_get_n_threads(void);

/*----------------------------------------------------------------------------
 * Set the minimum sunset sizes when renumbering for threads.
 *
 * parameters:
 *   min_i_subset_size <-- minimum number of interior faces per
 *                         thread per group
 *   min_b_subset_size <-- minimum number of boundary faces per
 *                         thread per group
 *----------------------------------------------------------------------------*/

void
cs_renumber_set_min_subset_size(cs_lnum_t  min_i_subset_size,
                                cs_lnum_t  min_b_subset_size);

/*----------------------------------------------------------------------------
 * Get the minimum sunset sizes when renumbering for threads.
 *
 *   min_i_subset_size --> minimum number of interior faces per
 *                         thread per group, or NULL
 *   min_b_subset_size --> minimum number of boundary faces per
 *                         thread per group, or NULL
 *----------------------------------------------------------------------------*/

void
cs_renumber_get_min_subset_size(cs_lnum_t  *min_i_subset_size,
                                cs_lnum_t  *min_b_subset_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select the algorithm for mesh renumbering.
 *
 * \param[in]  halo_adjacent_cells_last  if true, cells adjacent to ghost cells
 *                                       will be placed last
 *                                       (after pre-numbering)
 * \param[in]  halo_adjacent_faces_last  if true, interior faces adjacent to
 *                                       ghost cells will be placed last
 *                                       (after pre-numbering)
 * \param[in]  i_faces_base_ordering     pre-ordering of interior faces by
 *                                       lowest or highest adjacent cell id
 * \param[in]  cells_pre_numbering       algorithm for cells pre-numbering
 * \param[in]  cells_numbering           algorithm for cells numbering
 * \param[in]  i_faces_numbering         algorithm for interior faces numbering
 * \param[in]  b_faces_numbering         algorithm for boundary faces numbering
 * \param[in]  vertices_numbering        algorithm for vertices numbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_set_algorithm(bool                         halo_adjacent_cells_last,
                          bool                         halo_adjacent_faces_last,
                          cs_renumber_ordering_t       i_faces_base_ordering,
                          cs_renumber_cells_type_t     cells_pre_numbering,
                          cs_renumber_cells_type_t     cells_numbering,
                          cs_renumber_i_faces_type_t   i_faces_numbering,
                          cs_renumber_b_faces_type_t   b_faces_numbering,
                          cs_renumber_vertices_type_t  vertices_numbering);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the algorithms for mesh renumbering.
 *
 * Any argument may be passed NULL if this option is not queried.
 *
 * \param[out]  halo_adjacent_cells_last  if true, cells adjacent to ghost cells
 *                                        will be placed last
 *                                        (after pre-numbering)
 * \param[out]  halo_adjacent_faces_last  if true, interior faces adjacent to
 *                                        ghost cells will be placed last
 *                                        (after pre-numbering)
 * \param[out]  i_faces_base_ordering     pre-ordering of interior faces by
 *                                        lowest or highest adjacent cell id
 * \param[out]  cells_pre_numbering       algorithm for cells pre-numbering
 * \param[out]  cells_numbering           algorithm for cells numbering
 * \param[out]  i_faces_numbering         algorithm for interior faces numbering
 * \param[out]  b_faces_numbering         algorithm for boundary faces numbering
 * \param[out]  vertices_numbering        algorithm for vertices numbering
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_get_algorithm(bool                        *halo_adjacent_cells_last,
                          bool                        *halo_adjacent_faces_last,
                          cs_renumber_ordering_t      *i_faces_base_ordering,
                          cs_renumber_cells_type_t    *cells_pre_numbering,
                          cs_renumber_cells_type_t    *cells_numbering,
                          cs_renumber_i_faces_type_t  *i_faces_numbering,
                          cs_renumber_b_faces_type_t  *b_faces_numbering,
                          cs_renumber_vertices_type_t *vertices_numbering);

/*----------------------------------------------------------------------------
 * Renumber mesh elements for vectorization or threading depending on code
 * options and target machine.
 *
 * Renumbering cells may also allow improving locality (and favor faces
 * renumbering).
 * It is also possible to place cells connected to ghost cells last,
 * which may be useful to enable computation/communication overlap.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_mesh(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber cells depending on code options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_cells(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Renumber interior faces for vectorization or threading depending on code
 * options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_i_faces(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Renumber interior faces by global number.
 *
 * This effectively resets the interior faces to their initial numbering.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_i_faces_by_gnum(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Renumber boundary faces for vectorization or threading depending on code
 * options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_b_faces(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------
 * Renumber boundary faces by global number.
 *
 * This effectively resets the boundary faces to their initial numbering.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *----------------------------------------------------------------------------*/

void
cs_renumber_b_faces_by_gnum(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Renumber vertices depending on code options and target machine.
 *
 * parameters:
 *   mesh  <->  pointer to global mesh structure
 *
 * \param[in, out]  mesh  pointer to global mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_renumber_vertices(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RENUMBER_H__ */
