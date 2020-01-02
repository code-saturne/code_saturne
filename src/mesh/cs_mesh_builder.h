#ifndef __CS_MESH_BUILDER_H__
#define __CS_MESH_BUILDER_H__

/*============================================================================
 * Auxiliary structure used to read, write, and partition mesh data.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_interface.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Auxiliary and temporary structure used to build or distribute mesh */
/* ------------------------------------------------------------------ */

typedef struct {

  /* Global dimensions */

  cs_gnum_t     n_g_faces;             /* Number of faces */
  cs_gnum_t     n_g_face_connect_size; /* Size of face connectivity */

  int           n_perio;               /* Number of periodicities */

  bool          have_cell_rank;        /* True if cell_rank array is defined */
  bool          have_face_r_gen;       /* True if face level is defined */

  /* Temporary mesh data */

  cs_gnum_t    *face_cells;
  cs_lnum_t    *face_vertices_idx;
  cs_gnum_t    *face_vertices;
  int          *cell_gc_id;
  int          *face_gc_id;
  cs_real_t    *vertex_coords;

  /* Refinement features */

  char         *face_r_gen;

  /* Periodic features */

  int          *periodicity_num;         /* Periodicity numbers */
  cs_lnum_t    *n_per_face_couples;      /* Nb. face couples per periodicity */
  cs_gnum_t    *n_g_per_face_couples;    /* Global nb. couples per periodicity */

  cs_gnum_t   **per_face_couples;        /* Periodic face couples list. */

  /* Optional partitioning info */

  int          *cell_rank;               /* Partition id for each cell */

  /* Block ranges for parallel distribution */

  int                    min_rank_step;  /* Minimum block rank step */
  cs_block_dist_info_t   cell_bi;        /* Block info for cell data */
  cs_block_dist_info_t   face_bi;        /* Block info for face data */
  cs_block_dist_info_t   vertex_bi;      /* Block info for vertex data */
  cs_block_dist_info_t  *per_face_bi;    /* Block info for parallel face
                                            couples */

} cs_mesh_builder_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

extern cs_mesh_builder_t  *cs_glob_mesh_builder; /* Pointer to builder mesh
                                                    structure */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create an empty mesh builder structure.
 *
 * returns:
 *   A pointer to a mesh builder structure
 *----------------------------------------------------------------------------*/

cs_mesh_builder_t *
cs_mesh_builder_create(void);

/*----------------------------------------------------------------------------*
 * Destroy a cs_mesh_builder_t structure.
 *
 * parameters:
 *   mb <-> pointer to pointer of structure to destroy
 *----------------------------------------------------------------------------*/

void
cs_mesh_builder_destroy(cs_mesh_builder_t  **mb);

/*----------------------------------------------------------------------------
 * Define block distribution sizes for mesh builder.
 *
 * parameters:
 *   mb             <-> mesh builder
 *   rank_id        <-- id of local rank
 *   n_ranks        <-- number of associated ranks
 *   min_rank_step  <-- minimum rank step between blocks
 *   min_block_size <-- minimum number of entities per block
 *   n_g_cells      <-- global number of cells
 *   n_g_faces      <-- global number of faces
 *   n_g_vertices   <-- global number of vertices
 *----------------------------------------------------------------------------*/

void
cs_mesh_builder_define_block_dist(cs_mesh_builder_t  *mb,
                                  int                 rank_id,
                                  int                 n_ranks,
                                  int                 min_rank_step,
                                  int                 min_block_size,
                                  cs_gnum_t           n_g_cells,
                                  cs_gnum_t           n_g_faces,
                                  cs_gnum_t           n_g_vertices);

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Extract periodic face connectivity information from faces interface
 * for mesh builder whn running in parallel mode.
 *
 * parameters:
 *   n_init_perio <-- number of initial periodicities
 *   mesh         <-- pointer to mesh structure
 *   mb           <-> pointer to mesh builder structure
 *   periodicity  <--  periodicity information
 *   face_gnum    <--  global face numbers, or NULL
 *   face_ifs     <-- parallel and periodic faces interfaces set
 *----------------------------------------------------------------------------*/

void
cs_mesh_builder_extract_periodic_faces_g(int                        n_init_perio,
                                         cs_mesh_builder_t         *mb,
                                         fvm_periodicity_t         *periodicity,
                                         const cs_gnum_t           *face_gnum,
                                         const cs_interface_set_t  *face_ifs);

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_BUILDER_H__ */
