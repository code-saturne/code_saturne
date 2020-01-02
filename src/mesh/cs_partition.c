/*============================================================================
 * Partitioner
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

#include <assert.h>
#include <errno.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * METIS library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HAVE_PARMETIS)
#include <parmetis.h>
#endif

#if defined(HAVE_METIS)
#include <metis.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

/*----------------------------------------------------------------------------
 * SCOTCH library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_SCOTCH)
#include <stdint.h>
#include <scotch.h>
#elif defined(HAVE_PTSCOTCH)
#include <stdint.h>
#include <ptscotch.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"

#include "cs_all_to_all.h"
#include "cs_base.h"
#include "cs_block_dist.h"
#include "cs_block_to_part.h"
#include "cs_file.h"
#include "cs_io.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_builder.h"
#include "cs_part_to_block.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_partition.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Public Type documentation
 *============================================================================*/

/*
 * \enum cs_partition_stage_t cs_partition.h
 *
 * Partitioning is always done just after reading the mesh, unless a
 * partitioning input file is available, in which case the partitioning
 * read replaces this stage.
 *
 * When a mesh modification implying a change of cell connectivity graph
 * is expected, the mesh may be re-partitioned after the pre-processing
 * stage, prior to calculation. By default, re-partitioning is only done
 * if the partitioning algorithm chosen for that stage is expected to
 * produce different results due to the connectivity change. This is
 * the case for graph-based algorithms such as those of METIS or SCOTCH,
 * when mesh joining is defined, or additional periodic matching is defined
 * (and the algorithm is not configured to ignore periodicity information).
 *
 * There are thus two possible partitioning stages:
 *
 * - CS_PARTITION_FOR_PREPROCESS, which is optional, and occurs
 *   just  after reading the mesh.
 * - CS_PARTITION_MAIN, which occurs just after reading the mesh if
 *   it is the only stage,, or after mesh preprocessing (and before
 *   computation), if the partitioning for preprocessing stage is
 *   activated.
 *
 * The number of partitioning stages is determined automatically based on
 * information provided through \ref cs_partition_set_preprocess_hints,
 * but re-partitioning may also be forced or inhibited using the
 * \ref cs_partition_set_preprocess function.
 */

/*
 * \enum cs_partition_algorithm_t cs_partition.h
 *
 * Partitioning algorithm type.

 * If the default algorithm is selected, the choice will be based on the
 * following priority, depending on available libraries:
 * -  PT-SCOTCH (or SCOTCH if partitioning on one rank);
 * -  ParMETIS (or METIS if partitioning on one rank);
 * -  Morton space-filling curve (in bounding box)
 *
 * If both partitioning stages are active, the default for the preprocessing
 * stage will be based on the Morton space-filling curve (in bounding box),
 * as this should be cheaper, and the initial cell connectivity graph
 * is usually expected to be modified during preprocessing.
 *
 * \var CS_PARTITION_DEFAULT           Default partitioning (based on stage)
 * \var CS_PARTITION_SFC_MORTON_BOX    Morton (Z) curve in bounding box
 * \var CS_PARTITION_SFC_MORTON_CUBE   Morton (Z) curve in bounding cube
 * \var CS_PARTITION_SFC_HILBERT_BOX   Peano-Hilbert curve in bounding box
 * \var CS_PARTITION_SFC_HILBERT_CUBE  Peano-Hilbert curve in bounding cube
 * \var CS_PARTITION_SCOTCH            PT-SCOTCH or SCOTCH
 * \var CS_PARTITION_METIS             ParMETIS or METIS
 * \var CS_PARTITION_BLOCK             Unoptimized (naive) block partitioning
 * \var CS_PARTITION_NONE              No repartitioning (for computation
 *                                     stage after preprocessing)
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Local Type definitions
 *============================================================================*/

typedef double  _vtx_coords_t[3];

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_partition_algorithm_t   _part_algorithm[2] = {CS_PARTITION_DEFAULT,
                                                        CS_PARTITION_DEFAULT};
static int                        _part_rank_step[2] = {1, 1};
static bool                       _part_ignore_perio[2] = {false, false};

static int                        _part_compute_join_hint = false;
static int                        _part_compute_perio_hint = false;
static int                        _part_preprocess_active = 1; /* 0: inactive;
                                                                  1: default;
                                                                  2: active */
static int                        _part_write_output = 1;
static int                        _part_n_extra_partitions = 0;
static int                       *_part_extra_partitions_list = NULL;

static bool                       _part_uniform_sfc_block_size = false;

#if defined(WIN32) || defined(_WIN32)
static const char _dir_separator = '\\';
#else
static const char _dir_separator = '/';
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a reduced communicator
 *
 * parameters:
 *   part_step <-- step between active ranks
 *
 * returns:
 *   handle to reduced communicator
 *----------------------------------------------------------------------------*/

static MPI_Comm
_init_reduced_communicator(int  part_step)
{
  int n_ranks;
  int ranges[1][3];
  MPI_Group old_group, new_group;
  MPI_Comm part_comm = MPI_COMM_NULL;

  n_ranks = cs_glob_n_ranks;

  MPI_Barrier(cs_glob_mpi_comm); /* For debugging */

  MPI_Comm_size(cs_glob_mpi_comm, &n_ranks);
  MPI_Comm_group(cs_glob_mpi_comm, &old_group);

  ranges[0][0] = 0;
  ranges[0][1] = n_ranks - 1;
  ranges[0][2] = part_step;

  MPI_Group_range_incl(old_group, 1, ranges, &new_group);
  MPI_Comm_create(cs_glob_mpi_comm, new_group, &part_comm);
  MPI_Group_free(&new_group);

  MPI_Group_free(&old_group);

  MPI_Barrier(cs_glob_mpi_comm); /* For debugging */

  if (cs_glob_rank_id > -1 && cs_glob_rank_id % part_step) {
    /* MPI_Comm_free(&part_comm) does not work outside of the group */
    part_comm = MPI_COMM_NULL;
  }

  return part_comm;
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Display the distribution of cells per partition in serial mode
 *
 * parameters:
 *   cell_range <-- first and past-the-last cell numbers for this rank
 *   n_parts    <-- number of partitions
 *   part       <-- cell partition number
 *----------------------------------------------------------------------------*/

static void
_cell_part_histogram(cs_gnum_t   cell_range[2],
                     int         n_parts,
                     const   int part[])
{
  int i, k;
  size_t j;
  double step;

  cs_lnum_t *n_part_cells;
  cs_lnum_t n_min, n_max;
  size_t count[10];
  int n_steps = 10;
  size_t n_cells = 0;

  if (cell_range[1] > cell_range[0])
    n_cells = cell_range[1] - cell_range[0];

  if (n_parts <= 1) /* Should never happen */
    return;

  bft_printf(_("  Number of cells per domain (histogramm):\n"));

  BFT_MALLOC(n_part_cells, n_parts, cs_lnum_t);

  for (i = 0; i < n_parts; i++)
    n_part_cells[i] = 0;

  for (j = 0; j < n_cells; j++)
    n_part_cells[part[j]] += 1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    cs_lnum_t *n_part_cells_sum;
    BFT_MALLOC(n_part_cells_sum, n_parts, cs_lnum_t);
    MPI_Allreduce(n_part_cells, n_part_cells_sum, n_parts,
                  CS_MPI_LNUM, MPI_SUM, cs_glob_mpi_comm);
    BFT_FREE(n_part_cells);
    n_part_cells = n_part_cells_sum;
    n_part_cells_sum = NULL;
  }

#endif /* defined(HAVE_MPI) */

  /* Compute min and max */

  n_min = n_part_cells[0];
  n_max = n_part_cells[0];

  for (i = 1; i < n_parts; i++) {
    if (n_part_cells[i] > n_max)
      n_max = n_part_cells[i];
    else if (n_part_cells[i] < n_min)
      n_min = n_part_cells[i];
  }

  /* Define axis subdivisions */

  for (i = 0; i < n_steps; i++)
    count[i] = 0;

  if (n_max - n_min > 0) {

    if (n_max-n_min < n_steps)
      n_steps = n_max-n_min > 0 ? n_max-n_min : 1;

    step = (double)(n_max - n_min) / n_steps;

    /* Loop on partitions */

    for (i = 0; i < n_parts; i++) {

      /* Associated subdivision */

      for (j = 0, k = 1; k < n_steps; j++, k++) {
        if (n_part_cells[i] < n_min + k*step)
          break;
      }
      count[j] += 1;

    }

    for (i = 0, j = 1; i < n_steps - 1; i++, j++)
      bft_printf("    [ %10d ; %10d [ = %10d\n",
                 (int)(n_min + i*step),
                 (int)(n_min + j*step),
                 (int)(count[i]));

    bft_printf("    [ %10d ; %10d ] = %10d\n",
               (int)(n_min + (n_steps - 1)*step),
               (int)n_max,
               (int)(count[n_steps - 1]));

  }

  else { /* if (n_max == n_min) */
    bft_printf("    [ %10d ; %10d ] = %10d\n",
               (int)(n_min), (int)n_max, (int)n_parts);
  }

  BFT_FREE(n_part_cells);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Compute cell centers using minimal local data in parallel mode
 *
 * parameters:
 *   n_cells      <-- number of cells
 *   n_faces      <-- number of faces
 *   face_cells   <-- face -> cells connectivity
 *   face_vtx_idx <-- face -> vertices connectivity index
 *   face_vtx     <-- face -> vertices connectivity
 *   vtx_coord    <-- vertex coordinates
 *   cell_center  --> cell centers
 *----------------------------------------------------------------------------*/

static void
_cell_center_g(cs_lnum_t         n_cells,
               cs_lnum_t         n_faces,
               const cs_lnum_t   face_cells[],
               const cs_lnum_t   face_vtx_idx[],
               const cs_lnum_t   face_vtx[],
               const cs_real_t   vtx_coord[],
               cs_coord_t        cell_center[])
{
  cs_lnum_t i, j;
  cs_lnum_t vtx_id, face_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3], face_center[3];

  cs_lnum_t n_max_face_vertices = 0;

  _vtx_coords_t *face_vtx_coord = NULL;
  cs_coord_t *weight = NULL;

  const double surf_epsilon = 1e-24;

  assert(face_vtx_idx[0] == 0);

  BFT_MALLOC(weight, n_cells, cs_coord_t);

  for (i = 0; i < n_cells; i++) {
    weight[i] = 0.0;
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] = 0.0;
  }

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (face_id = 0; face_id < n_faces; face_id++) {
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices, _vtx_coords_t);

  /* Loop on each face */

  for (face_id = 0; face_id < n_faces; face_id++) {

    /* Initialization */

    cs_lnum_t tri_id;

    cs_lnum_t cell_id_0 = face_cells[face_id*2] -1;
    cs_lnum_t cell_id_1 = face_cells[face_id*2 + 1] -1;
    cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
    cs_coord_t face_surface = 0.0;

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id];
    end_id = face_vtx_idx[face_id + 1];

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      cs_lnum_t shift = 3 * (face_vtx[vtx_id] - 1);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    /* Compute the barycenter of the face vertices */

    for (i = 0; i < 3; i++) {
      vtx_cog[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        vtx_cog[i] += face_vtx_coord[vtx_id][i];
      vtx_cog[i] /= n_face_vertices;
    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycenter) */

    for (i = 0; i < 3; i++) {
      ref_normal[i] = 0.;
      face_center[i] = 0.0;
    }

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      cs_coord_t tri_surface;
      cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

      cs_lnum_t id0 = tri_id;
      cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

      /* Normal for each triangle */

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
        vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
      }

      tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
      tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
      tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

      if (tri_id == 0) {
        for (i = 0; i < 3; i++)
          ref_normal[i] = tri_normal[i];
      }

      /* Center of gravity for a triangle */

      for (i = 0; i < 3; i++) {
        tri_center[i] = (  vtx_cog[i]
                         + face_vtx_coord[id0][i]
                         + face_vtx_coord[id1][i]) / 3.0;
      }

      tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                         + tri_normal[1]*tri_normal[1]
                         + tri_normal[2]*tri_normal[2]) * 0.5;

      if ((  tri_normal[0]*ref_normal[0]
           + tri_normal[1]*ref_normal[1]
           + tri_normal[2]*ref_normal[2]) < 0.0)
        tri_surface *= -1.0;

      /* Now compute contribution to face center and surface */

      face_surface += tri_surface;

      for (i = 0; i < 3; i++) {
        face_center[i] += tri_surface * tri_center[i];
        unweighted_center[i] = tri_center[i];
      }

    } /* End of loop  on triangles of the face */

    if (face_surface > surf_epsilon) {
      for (i = 0; i < 3; i++)
        face_center[i] /= face_surface;
    }
    else {
      face_surface = surf_epsilon;
      for (i = 0; i < 3; i++)
        face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
    }

    /* Now contribute to cell centers */

    assert(cell_id_0 > -2 && cell_id_1 > -2);

    if (cell_id_0 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_0*3 + i] += face_center[i]*face_surface;
      weight[cell_id_0] += face_surface;
    }

    if (cell_id_1 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_1*3 + i] += face_center[i]*face_surface;
      weight[cell_id_1] += face_surface;
    }

  } /* End of loop on faces */

  BFT_FREE(face_vtx_coord);

  for (i = 0; i < n_cells; i++) {
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] /= weight[i];
  }

  BFT_FREE(weight);
}

/*----------------------------------------------------------------------------
 * Compute cell centers using block data in parallel mode.
 *
 * parameters:
 *   mb          <-- pointer to mesh builder helper structure
 *   cell_center --> cell centers array
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_precompute_cell_center_g(const cs_mesh_builder_t  *mb,
                          cs_coord_t                cell_center[],
                          MPI_Comm                  comm)
{
  cs_lnum_t i;
  int n_ranks = 0;

  cs_datatype_t gnum_type = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
  cs_datatype_t real_type = (sizeof(cs_real_t) == 8) ? CS_DOUBLE : CS_FLOAT;

  cs_lnum_t _n_cells = 0;
  cs_lnum_t _n_faces = 0;
  cs_lnum_t _n_vertices = 0;

  cs_gnum_t *_cell_num = NULL;
  cs_gnum_t *_face_num = NULL;
  cs_gnum_t *_vtx_num = NULL;
  cs_gnum_t *_face_gcells = NULL;
  cs_gnum_t *_face_gvertices = NULL;

  cs_lnum_t *_face_cells = NULL;
  cs_lnum_t *_face_vertices_idx = NULL;
  cs_lnum_t *_face_vertices = NULL;

  cs_real_t *_vtx_coord = NULL;

  cs_block_to_part_t *d = NULL;

  /* Initialization */

  MPI_Comm_size(comm, &n_ranks);

  assert((sizeof(cs_lnum_t) == 4) || (sizeof(cs_lnum_t) == 8));

  _n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

  BFT_MALLOC(_cell_num, _n_cells, cs_gnum_t);

  for (i = 0; i < _n_cells; i++)
    _cell_num[i] = mb->cell_bi.gnum_range[0] + i;

  /* Distribute faces */
  /*------------------*/

  d = cs_block_to_part_create_by_adj_s(comm,
                                       mb->face_bi,
                                       mb->cell_bi,
                                       2,
                                       mb->face_cells,
                                       NULL,
                                       NULL);

  _n_faces = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_face_gcells, _n_faces*2, cs_gnum_t);

  /* Face -> cell connectivity */

  cs_block_to_part_copy_array(d,
                              gnum_type,
                              2,
                              mb->face_cells,
                              _face_gcells);

  /* Now convert face -> cell connectivity to local cell numbers */

  BFT_MALLOC(_face_cells, _n_faces*2, cs_lnum_t);

  cs_block_to_part_global_to_local(_n_faces*2,
                                   1,
                                   _n_cells,
                                   true,
                                   _cell_num,
                                   _face_gcells,
                                   _face_cells);

  BFT_FREE(_cell_num);
  BFT_FREE(_face_gcells);

  /* Face connectivity */

  BFT_MALLOC(_face_vertices_idx, _n_faces + 1, cs_lnum_t);

  cs_block_to_part_copy_index(d,
                              mb->face_vertices_idx,
                              _face_vertices_idx);

  BFT_MALLOC(_face_gvertices, _face_vertices_idx[_n_faces], cs_gnum_t);

  cs_block_to_part_copy_indexed(d,
                                gnum_type,
                                mb->face_vertices_idx,
                                mb->face_vertices,
                                _face_vertices_idx,
                                _face_gvertices);

  _face_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Vertices */

  d = cs_block_to_part_create_adj(comm,
                                  mb->vertex_bi,
                                  _face_vertices_idx[_n_faces],
                                  _face_gvertices);

  _n_vertices = cs_block_to_part_get_n_part_ents(d);

  BFT_MALLOC(_vtx_coord, _n_vertices*3, cs_real_t);

  cs_block_to_part_copy_array(d,
                              real_type,
                              3,
                              mb->vertex_coords,
                              _vtx_coord);

  _vtx_num = cs_block_to_part_transfer_gnum(d);

  cs_block_to_part_destroy(&d);

  /* Now convert face -> vertex connectivity to local vertex numbers */

  BFT_MALLOC(_face_vertices, _face_vertices_idx[_n_faces], cs_lnum_t);

  cs_block_to_part_global_to_local(_face_vertices_idx[_n_faces],
                                   1,
                                   _n_vertices,
                                   true,
                                   _vtx_num,
                                   _face_gvertices,
                                   _face_vertices);

  BFT_FREE(_face_gvertices);

  _cell_center_g(_n_cells,
                 _n_faces,
                 _face_cells,
                 _face_vertices_idx,
                 _face_vertices,
                 _vtx_coord,
                 cell_center);

  BFT_FREE(_vtx_coord);
  BFT_FREE(_vtx_num);

  BFT_FREE(_face_cells);

  BFT_FREE(_face_vertices_idx);
  BFT_FREE(_face_vertices);

  BFT_FREE(_face_num);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Compute cell centers using block data in serial mode.
 *
 * parameters:
 *   mb          <-- pointer to mesh builder helper structure
 *   cell_center --> cell centers array
 *----------------------------------------------------------------------------*/

static void
_precompute_cell_center_l(const cs_mesh_builder_t  *mb,
                          cs_coord_t                cell_center[])
{
  cs_lnum_t i, j;
  cs_lnum_t vtx_id, face_id, start_id, end_id;
  cs_lnum_t n_face_vertices;
  cs_coord_t ref_normal[3], vtx_cog[3], face_center[3];

  cs_lnum_t n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];
  cs_lnum_t n_faces = mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0];

  const cs_gnum_t  *face_cells = mb->face_cells;
  const cs_lnum_t  *face_vtx_idx = mb->face_vertices_idx;
  const cs_gnum_t  *face_vtx = mb->face_vertices;
  const cs_real_t  *vtx_coord = mb->vertex_coords;

  cs_lnum_t n_max_face_vertices = 0;

  _vtx_coords_t *face_vtx_coord = NULL;
  cs_coord_t *weight = NULL;

  const double surf_epsilon = 1e-24;

  assert(face_vtx_idx[0] == 0);

  BFT_MALLOC(weight, n_cells, cs_coord_t);

  for (i = 0; i < n_cells; i++) {
    weight[i] = 0.0;
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] = 0.0;
  }

  /* Counting and allocation */

  n_max_face_vertices = 0;

  for (face_id = 0; face_id < n_faces; face_id++) {
    n_face_vertices = face_vtx_idx[face_id + 1] - face_vtx_idx[face_id];
    if (n_max_face_vertices <= n_face_vertices)
      n_max_face_vertices = n_face_vertices;
  }

  BFT_MALLOC(face_vtx_coord, n_max_face_vertices, _vtx_coords_t);

  /* Loop on each face */

  for (face_id = 0; face_id < n_faces; face_id++) {

    /* Initialization */

    cs_lnum_t tri_id;

    cs_lnum_t cell_id_0 = face_cells[face_id*2] -1;
    cs_lnum_t cell_id_1 = face_cells[face_id*2 + 1] -1;
    cs_coord_t unweighted_center[3] = {0.0, 0.0, 0.0};
    cs_coord_t face_surface = 0.0;

    n_face_vertices = 0;

    start_id = face_vtx_idx[face_id];
    end_id = face_vtx_idx[face_id + 1];

    /* Define the polygon (P) according to the vertices (Pi) of the face */

    for (vtx_id = start_id; vtx_id < end_id; vtx_id++) {

      cs_lnum_t shift = 3 * (face_vtx[vtx_id] - 1);
      for (i = 0; i < 3; i++)
        face_vtx_coord[n_face_vertices][i] = vtx_coord[shift + i];
      n_face_vertices++;

    }

    /* Compute the barycenter of the face vertices */

    for (i = 0; i < 3; i++) {
      vtx_cog[i] = 0.0;
      for (vtx_id = 0; vtx_id < n_face_vertices; vtx_id++)
        vtx_cog[i] += face_vtx_coord[vtx_id][i];
      vtx_cog[i] /= n_face_vertices;
    }

    /* Loop on the triangles of the face (defined by an edge of the face
       and its barycenter) */

    for (i = 0; i < 3; i++) {
      ref_normal[i] = 0.;
      face_center[i] = 0.0;
    }

    for (tri_id = 0 ; tri_id < n_face_vertices ; tri_id++) {

      cs_coord_t tri_surface;
      cs_coord_t vect1[3], vect2[3], tri_normal[3], tri_center[3];

      cs_lnum_t id0 = tri_id;
      cs_lnum_t id1 = (tri_id + 1)%n_face_vertices;

      /* Normal for each triangle */

      for (i = 0; i < 3; i++) {
        vect1[i] = face_vtx_coord[id0][i] - vtx_cog[i];
        vect2[i] = face_vtx_coord[id1][i] - vtx_cog[i];
      }

      tri_normal[0] = vect1[1] * vect2[2] - vect2[1] * vect1[2];
      tri_normal[1] = vect2[0] * vect1[2] - vect1[0] * vect2[2];
      tri_normal[2] = vect1[0] * vect2[1] - vect2[0] * vect1[1];

      if (tri_id == 0) {
        for (i = 0; i < 3; i++)
          ref_normal[i] = tri_normal[i];
      }

      /* Center of gravity for a triangle */

      for (i = 0; i < 3; i++) {
        tri_center[i] = (  vtx_cog[i]
                         + face_vtx_coord[id0][i]
                         + face_vtx_coord[id1][i]) / 3.0;
      }

      tri_surface = sqrt(  tri_normal[0]*tri_normal[0]
                         + tri_normal[1]*tri_normal[1]
                         + tri_normal[2]*tri_normal[2]) * 0.5;

      if ((  tri_normal[0]*ref_normal[0]
           + tri_normal[1]*ref_normal[1]
           + tri_normal[2]*ref_normal[2]) < 0.0)
        tri_surface *= -1.0;

      /* Now compute contribution to face center and surface */

      face_surface += tri_surface;

      for (i = 0; i < 3; i++) {
        face_center[i] += tri_surface * tri_center[i];
        unweighted_center[i] = tri_center[i];
      }

    } /* End of loop  on triangles of the face */

    if (face_surface > surf_epsilon) {
      for (i = 0; i < 3; i++)
        face_center[i] /= face_surface;
    }
    else {
      face_surface = surf_epsilon;
      for (i = 0; i < 3; i++)
        face_center[i] = unweighted_center[i] * face_surface / n_face_vertices;
    }

    /* Now contribute to cell centers */

    assert(cell_id_0 > -2 && cell_id_1 > -2);

    if (cell_id_0 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_0*3 + i] += face_center[i]*face_surface;
      weight[cell_id_0] += face_surface;
    }

    if (cell_id_1 > -1) {
      for (i = 0; i < 3; i++)
        cell_center[cell_id_1*3 + i] += face_center[i]*face_surface;
      weight[cell_id_1] += face_surface;
    }

  } /* End of loop on faces */

  BFT_FREE(face_vtx_coord);

  for (i = 0; i < n_cells; i++) {
    for (j = 0; j < 3; j++)
      cell_center[i*3 + j] /= weight[i];
  }

  BFT_FREE(weight);
}

/*----------------------------------------------------------------------------
 * Define cell ranks using a space-filling curve.
 *
 * parameters:
 *   n_g_cells   <-- global number of cells
 *   n_ranks     <-- number of ranks in partition
 *   mb          <-- pointer to mesh builder helper structure
 *   sfc_type    <-- type of space-filling curve
 *   cell_rank   --> cell rank (1 to n numbering)
 *   comm        <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

static void
_cell_rank_by_sfc(cs_gnum_t                 n_g_cells,
                  int                       n_ranks,
                  const cs_mesh_builder_t  *mb,
                  fvm_io_num_sfc_t          sfc_type,
                  int                       cell_rank[],
                  MPI_Comm                  comm)

#else

static void
_cell_rank_by_sfc(cs_gnum_t                 n_g_cells,
                  int                       n_ranks,
                  const cs_mesh_builder_t  *mb,
                  fvm_io_num_sfc_t          sfc_type,
                  int                       cell_rank[])

#endif
{
  cs_lnum_t i;
  cs_timer_t  start_time, end_time;
  cs_timer_counter_t dt;

  cs_lnum_t n_cells = 0, block_size = 0;

  cs_coord_t *cell_center = NULL;
  fvm_io_num_t *cell_io_num = NULL;
  const cs_gnum_t *cell_num = NULL;

  bft_printf(_("\n Partitioning by space-filling curve: %s.\n"),
             _(fvm_io_num_sfc_type_name[sfc_type]));

  start_time = cs_timer_time();

  n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];
  block_size = mb->cell_bi.block_size;

  BFT_MALLOC(cell_center, n_cells*3, cs_coord_t);

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    _precompute_cell_center_g(mb, cell_center, comm);
#endif
  if (n_ranks == 1)
    _precompute_cell_center_l(mb, cell_center);

  end_time = cs_timer_time();
  dt = cs_timer_diff(&start_time, &end_time);
  start_time = end_time;

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  precompute cell centers:    %.3g s\n"),
                (double)(dt.wall_nsec)/1.e9);

  cell_io_num = fvm_io_num_create_from_sfc(cell_center,
                                           3,
                                           n_cells,
                                           sfc_type);

  BFT_FREE(cell_center);

  cell_num = fvm_io_num_get_global_num(cell_io_num);

  block_size = n_g_cells / n_ranks;
  if (n_g_cells % n_ranks)
    block_size += 1;

  /* Determine rank based on global numbering with SFC ordering; */

  if (_part_uniform_sfc_block_size == false) {

    cs_gnum_t cells_per_rank = n_g_cells / n_ranks;
    cs_lnum_t rmdr = n_g_cells - cells_per_rank * (cs_gnum_t)n_ranks;

    if (rmdr == 0) {
      for (i = 0; i < n_cells; i++)
        cell_rank[i] = (cell_num[i] - 1) / cells_per_rank;
    }
    else {
      cs_gnum_t n_ranks_rmdr = n_ranks - rmdr;
      cs_gnum_t n_ranks_cells_per_rank = n_ranks_rmdr * cells_per_rank;
      for (i = 0; i < n_cells; i++) {
        if ((cell_num[i] - 1)  <  n_ranks_cells_per_rank)
          cell_rank[i] = (cell_num[i] - 1) / cells_per_rank;
        else
          cell_rank[i] = (cell_num[i] + n_ranks_rmdr - 1) / (cells_per_rank + 1);
      }
    }

  }

  else {

    /* Plan for case where we would need a fixed block size,
       for example, using an external linear solver assuming this.
       This may not work at high process counts, where the last
       ranks will have no data (a solution to this would be
       to build a slightly smaller MPI communicator). */

    for (i = 0; i < n_cells; i++) {
      cell_rank[i] = ((cell_num[i] - 1) / block_size);
      assert(cell_rank[i] > -1 && cell_rank[i] < n_ranks);
    }

  }

  cell_io_num = fvm_io_num_destroy(cell_io_num);

  end_time = cs_timer_time();
  dt = cs_timer_diff(&start_time, &end_time);

  if (sfc_type < FVM_IO_NUM_SFC_HILBERT_BOX)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Morton (Z) curve:           %.3g s\n"),
                  (double)(dt.wall_nsec)/1.e9);
  else
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  Peano-Hilbert curve:        %.3g s\n"),
                  (double)(dt.wall_nsec)/1.e9);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Update global face -> cells adjacency for periodicity in parallel mode.
 *
 * If face Fi is periodic with face Fj, and face Fi is adjacent to cell Ci,
 * while face Fj is adjacent to face Cj, we add Cj to Fi's adjacent cells,
 * and Ci to Fj's adjacent cells, just as if Fi and Fj were regular interior
 * faces (this ignores the geometric transformation, but is suffient to
 * build the cells -> cells graph for domain partitioning).
 *
 * This function should be called when faces are distributed by blocks,
 * that is prior to redistribution based on face -> cells adjacency.
 *
 * arguments:
 *   bi                 <-- block size and range info for faces
 *   g_face_cells       <-> global face->cells adjacency to update
 *   n_periodic_couples <-- number of periodic couples
 *   periodic_couples   <-- array indicating periodic couples (interlaced,
 *                          using global numberings)
 *   comm               <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_add_perio_to_face_cells_g(cs_block_dist_info_t  bi,
                           cs_gnum_t             g_face_cells[],
                           cs_lnum_t             n_periodic_couples,
                           const cs_gnum_t       periodic_couples[],
                           MPI_Comm              comm)
{
  /* Initialization */

  int flags = 0;

  const cs_lnum_t n_send = n_periodic_couples*2;

  cs_all_to_all_t *d = cs_all_to_all_create_from_block(n_send,
                                                       flags,
                                                       periodic_couples,
                                                       bi,
                                                       comm);

  /* Distribute to blocks */

  cs_gnum_t *b_data = cs_all_to_all_copy_array(d,
                                               CS_GNUM_TYPE,
                                               1,
                                               false, /* reverse */
                                               periodic_couples,
                                               NULL);

  cs_lnum_t  n_b = cs_all_to_all_n_elts_dest(d);

  /* Receive buffer contains global cell face whose cell adjacency
     is defined on the local rank, and we replace the received value
     by the adjacent cell number, for return exchange */

  for (cs_lnum_t i = 0; i < n_b; i++) {
    cs_gnum_t g_face_id = b_data[i] - bi. gnum_range[0];
    cs_gnum_t c_num_0 = g_face_cells[g_face_id*2];
    const cs_gnum_t c_num_1 = g_face_cells[g_face_id*2 + 1];
    assert(c_num_0 == 0 || c_num_1 == 0);
    /* assign c_num_0 or c_num_1 depending on which is nonzero
       (or 0 if both are 0, which should not happen) */
    b_data[i] = c_num_0 + c_num_1;
  }

  cs_gnum_t *send_adj;
  BFT_MALLOC(send_adj, n_send*2, cs_gnum_t);

  cs_gnum_t *r_data = cs_all_to_all_copy_array(d,
                                               CS_GNUM_TYPE,
                                               1,
                                               true, /* reverse */
                                               b_data,
                                               NULL);

  BFT_FREE(b_data);

  /* Now r_data contains the global cell number matching a given face;
     Send global face number and cell number adjacent with its
     periodic face to rank defining the adjacency for this face */

  for (cs_lnum_t i = 0; i < n_send; i++) {

    cs_lnum_t k = (i + 1) % 2;  /* 1 for first value, 0 for
                                   second (permutation) */

    send_adj[i*2]     = periodic_couples[i];
    send_adj[i*2 + 1] = r_data[i + k];

  }

  BFT_FREE(r_data);

  b_data = cs_all_to_all_copy_array(d,
                                    CS_GNUM_TYPE,
                                    2,
                                    false, /* reverse */
                                    send_adj,
                                    NULL);

  /* Update face -> cell connectivity */

  for (cs_lnum_t i = 0; i < n_b; i++) {
    const cs_gnum_t g_face_id = b_data[2*i] - bi. gnum_range[0];
    const cs_gnum_t g_cell_num = b_data[2*i + 1];
    if (g_face_cells[g_face_id*2] == 0)
      g_face_cells[g_face_id*2] = g_cell_num;
    else {
      assert(g_face_cells[g_face_id*2 + 1] == 0);
      g_face_cells[g_face_id*2 + 1] = g_cell_num;
    }
  }

  BFT_FREE(b_data);
  BFT_FREE(send_adj);

  cs_all_to_all_destroy(&d);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Update global face -> cells adjacency for periodicity in local mode.
 *
 * If face Fi is periodic with face Fj, and face Fi is adjacent to cell Ci,
 * while face Fj is adjacent to face Cj, we add Cj to Fi's adjacent cells,
 * and Ci to Fj's adjacent cells, just as if Fi and Fj were regular interior
 * faces (this ignores the geometric transformation, but is suffient to
 * build the cells -> cells graph for domain partitioning).
 *
 * This function should be called when faces are distributed by blocks,
 * that is prior to redistribution based on face -> cells adjacency.
 *
 * arguments:
 *   bi                 <-- block size and range info for faces
 *   g_face_cells       <-> global face->cells adjacency to update
 *   n_periodic_couples <-- number of periodic couples
 *   periodic_couples   <-- array indicating periodic couples (interlaced,
 *                          using global numberings)
 *
 * returns:
 *   initialized partition to block distributor
 *----------------------------------------------------------------------------*/

static void
_add_perio_to_face_cells_l(cs_gnum_t         g_face_cells[],
                           cs_lnum_t         n_periodic_couples,
                           const cs_gnum_t   periodic_couples[])
{
  cs_lnum_t j;

  /* Update face -> cell connectivity */

  for (j = 0; j < n_periodic_couples; j++) {

    cs_gnum_t f_id_0 = periodic_couples[j*2] - 1;
    cs_gnum_t f_id_1 = periodic_couples[j*2 + 1] - 1;
    cs_gnum_t c_num_00 = g_face_cells[f_id_0*2];
    cs_gnum_t c_num_01 = g_face_cells[f_id_0*2 + 1];
    cs_gnum_t c_num_10 = g_face_cells[f_id_1*2];
    cs_gnum_t c_num_11 = g_face_cells[f_id_1*2 + 1];

    assert(c_num_00 == 0 || c_num_01 == 0);
    assert(c_num_10 == 0 || c_num_11 == 0);

    /* assign c_num_i0 or c_num_i1 depending on which is nonzero
       (or 0 if both are 0, which should not happen) */

    if (g_face_cells[f_id_0*2] == 0)
      g_face_cells[f_id_0*2] = c_num_10 + c_num_11;
    else
      g_face_cells[f_id_0*2 + 1] = c_num_00 + c_num_01;

    if (g_face_cells[f_id_1*2] == 0)
      g_face_cells[f_id_1*2] = c_num_00 + c_num_01;
    else
      g_face_cells[f_id_1*2 + 1] = c_num_10 + c_num_11;
  }
}

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

/*----------------------------------------------------------------------------
 * Build cell -> cell connectivity
 *
 * parameters:
 *   n_cells        <-- number of cells in mesh
 *   n_faces        <-- number of cells in mesh
 *   start_cell     <-- number of first cell for the curent rank
 *   face_cells     <-- face->cells connectivity
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *----------------------------------------------------------------------------*/

static void
_metis_cell_cells(size_t       n_cells,
                  size_t       n_faces,
                  cs_gnum_t    start_cell,
                  cs_gnum_t   *face_cells,
                  idx_t      **cell_idx,
                  idx_t      **cell_neighbors)
{
  size_t i, id_0, id_1;

  cs_gnum_t  c_num[2];

  idx_t  *n_neighbors;
  idx_t  *_cell_idx;
  idx_t  *_cell_neighbors;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, idx_t);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    assert(   c_num[0] - start_cell < n_cells
           || c_num[1] - start_cell < n_cells);

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells)
      n_neighbors[c_num[0] - start_cell] += 1;
    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells)
      n_neighbors[c_num[1] - start_cell] += 1;
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, idx_t);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++)
    _cell_idx[i + 1] = _cell_idx[i] + n_neighbors[i];

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_cells], idx_t);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells) {
      id_0 = c_num[0] - start_cell;
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = c_num[1] - 1;
      n_neighbors[id_0] += 1;
    }

    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells) {
      id_1 = c_num[1] - start_cell;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = c_num[0] - 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_FREE(n_neighbors);

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
}

/*----------------------------------------------------------------------------
 * Compute partition using METIS
 *
 * parameters:
 *   n_cells       <-- number of cells in mesh
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *----------------------------------------------------------------------------*/

static void
_part_metis(size_t   n_cells,
            int      n_parts,
            idx_t   *cell_idx,
            idx_t   *cell_neighbors,
            int     *cell_part)
{
  size_t i;
  double  start_time, end_time;

  idx_t   _n_constraints = 1;

  idx_t    edgecut    = 0; /* <-- Number of faces on partition */

  idx_t   _n_cells = n_cells;
  idx_t   _n_parts = n_parts;
  idx_t  *_cell_part = NULL;

  start_time = cs_timer_wtime();

  if (sizeof(idx_t) == sizeof(int))
    _cell_part = (idx_t *)cell_part;

  else
    BFT_MALLOC(_cell_part, n_cells, idx_t);

  if (n_parts < 8) {

    bft_printf(_("\n"
                 " Partitioning %llu cells to %d domains"
                 "   (METIS_PartGraphRecursive).\n"),
               (unsigned long long)n_cells, n_parts);

    METIS_PartGraphRecursive(&_n_cells,
                             &_n_constraints,
                             cell_idx,
                             cell_neighbors,
                             NULL,       /* vwgt:   cell weights */
                             NULL,       /* vsize:  size of the vertices */
                             NULL,       /* adjwgt: face weights */
                             &_n_parts,
                             NULL,       /* tpwgts */
                             NULL,       /* ubvec: load imbalance tolerance */
                             NULL,       /* options */
                             &edgecut,
                             _cell_part);

  }

  else {

    bft_printf(_("\n"
                 " Partitioning %llu cells to %d domains\n"
                 "  (METIS_PartGraphKway).\n"),
               (unsigned long long)n_cells, n_parts);

    METIS_PartGraphKway(&_n_cells,
                        &_n_constraints,
                        cell_idx,
                        cell_neighbors,
                        NULL,       /* vwgt:   cell weights */
                        NULL,       /* vsize:  size of the vertices */
                        NULL,       /* adjwgt: face weights */
                        &_n_parts,
                        NULL,       /* tpwgts */
                        NULL,       /* ubvec: load imbalance tolerance */
                        NULL,       /* options */
                        &edgecut,
                        _cell_part);

  }

  end_time = cs_timer_wtime();

  bft_printf(_("\n"
               "  Total number of faces on parallel boundaries: %llu\n"
               "  wall-clock time: %f s\n\n"),
             (unsigned long long)edgecut,
             (double)(end_time - start_time));

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  METIS_PartGraphKway:        %.3g s\n",
                (double)(end_time - start_time));

  if (sizeof(idx_t) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
}

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_PARMETIS)

/*----------------------------------------------------------------------------
 * Compute partition using ParMETIS
 *
 * parameters:
 *   n_g_cells     <-- global number of cells
 *   cell_range    <-- first and past-the-last cell numbers for this rank
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_part_parmetis(cs_gnum_t   n_g_cells,
               cs_gnum_t   cell_range[2],
               int         n_parts,
               idx_t      *cell_idx,
               idx_t      *cell_neighbors,
               int        *cell_part,
               MPI_Comm    comm)
{
  size_t i;
  double  start_time, end_time;

  unsigned long long  edgecut    = 0; /* <-- Number of faces on partition */

  int       n_ranks;
  size_t    n_cells = cell_range[1] - cell_range[0];
  idx_t     vtxstart = cell_range[0] - 1;
  idx_t     vtxend = cell_range[1] - 1;
  idx_t    *vtxdist = NULL;
  idx_t    *_cell_part = NULL;
  MPI_Datatype mpi_idx_t = MPI_DATATYPE_NULL;

  start_time = cs_timer_wtime();

  MPI_Comm_size(comm, &n_ranks);

  /* Adjust mpi_idx_t if necessary */

  if (sizeof(idx_t) == sizeof(short))
    mpi_idx_t = MPI_SHORT;
  else if (sizeof(idx_t) == sizeof(int))
    mpi_idx_t = MPI_INT;
  else if (sizeof(idx_t) == sizeof(long))
    mpi_idx_t = MPI_LONG; /* standard ParMETIS 3.1.1 only short or int */
  else {
    assert(0); /* porting error, possible with future or modified ParMETIS */
  }

  if (sizeof(idx_t) == sizeof(int))
    _cell_part = cell_part;

  else
    BFT_MALLOC(_cell_part, n_cells, idx_t);

  bft_printf(_("\n"
               " Partitioning %llu cells to %d domains on %d ranks\n"
               "   (ParMETIS_V3_PartKway).\n"),
             (unsigned long long)n_g_cells, n_parts, n_ranks);

  /* Build vtxdist */

  BFT_MALLOC(vtxdist, n_ranks + 1, idx_t);

  MPI_Allgather(&vtxstart, 1, mpi_idx_t,
                vtxdist, 1, mpi_idx_t, comm);
  MPI_Allreduce(&vtxend, vtxdist + n_ranks, 1, mpi_idx_t, MPI_MAX, comm);

  /* Call ParMETIS partitioning */

  {
    int      j;
    idx_t  _edgecut = 0;
    idx_t  _n_parts = n_parts;
    idx_t  ncon     = 1; /* number of weights for each vertex */
    idx_t  options[3] = {0, 1, 15}; /* By default if options[0] = 0 */
    idx_t  numflag  = 0; /* 0 to n-1 numbering (C type) */
    idx_t  wgtflag  = 0; /* No weighting for faces or cells */

    real_t wgt = 1.0/n_parts;
    real_t ubvec[]  = {1.5};
    real_t *tpwgts = NULL;

    BFT_MALLOC(tpwgts, n_parts, real_t);

    for (j = 0; j < n_parts; j++)
      tpwgts[j] = wgt;

    int retval = ParMETIS_V3_PartKway
                   (vtxdist,
                    cell_idx,
                    cell_neighbors,
                    NULL,       /* vwgt:   cell weights */
                    NULL,       /* adjwgt: face weights */
                    &wgtflag,
                    &numflag,
                    &ncon,
                    &_n_parts,
                    tpwgts,
                    ubvec,
                    options,
                    &_edgecut,
                    _cell_part,
                    &comm);

    BFT_FREE(tpwgts);

    edgecut = _edgecut;

    if (retval != METIS_OK)
      bft_error(__FILE__, __LINE__, 0,
                _("Error in ParMETIS partitioning.\n"));

  }

  end_time = cs_timer_wtime();

  BFT_FREE(vtxdist);

  if (edgecut > 0)
    bft_printf(_("\n"
                 "  Total number of faces on parallel boundaries: %llu\n"
                 "  wall-clock time: %f s\n\n"),
               (unsigned long long)edgecut,
               (double)(end_time - start_time));
  else
    bft_printf(_("\n"
                 "  wall-clock time: %f s\n\n"),
               (double)(end_time - start_time));

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  ParMETIS_V3_PartKway:       %.3g s\n",
                (double)(end_time - start_time));

  if (sizeof(idx_t) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }
}

#endif /* defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Sort an array "a" between its left bound "l" and its right bound "r"
 * thanks to a shell sort (Knuth algorithm).
 *
 * parameters:
 *   l <-- left bound
 *   r <-- right bound
 *   a <-> array to sort
 *---------------------------------------------------------------------------*/

static void
_scotch_sort_shell(SCOTCH_Num  l,
                   SCOTCH_Num  r,
                   SCOTCH_Num  a[])
{
  int i, j, h;

  /* Compute stride */
  for (h = 1; h <= (r-l)/9; h = 3*h+1) ;

  /* Sort array */
  for (; h > 0; h /= 3) {

    for (i = l+h; i < r; i++) {

      SCOTCH_Num v = a[i];

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
 * Print an error message and exit with an EXIT_FAILURE code.
 *
 * An implementation of this function is required by libScotch.
 *
 * parameters:
 *   errstr <-- format string, as printf() and family.
 *   ...    <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

void
SCOTCH_errorPrint(const char  *errstr,
                  ...)
{
  if (cs_glob_rank_id < 1) {

    va_list  errlist;

    fflush(stdout);

    fprintf(stderr, "\n");

    fprintf(stderr, _("\nFatal error encountered by libScotch.\n\n"));

    va_start(errlist, errstr);
    vfprintf(stderr, errstr, errlist);
    va_end(errlist);
    fprintf(stderr, "\n\n");
    fflush(stderr);
  }

  assert(0);

#if defined(HAVE_MPI)
  {
    int mpi_flag;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag != 0)
      MPI_Abort(cs_glob_mpi_comm, EXIT_FAILURE);
  }
#endif /* HAVE_MPI */

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Print a warning message.
 *
 * An implementation of this function is required by libScotch.
 *
 * parameters:
 *   errstr <-- format string, as printf() and family.
 *   ...    <-- variable arguments based on format string.
 *----------------------------------------------------------------------------*/

void
SCOTCH_errorPrintW (const char *errstr,
                    ...)
{
  if (cs_glob_rank_id < 1) {

    va_list  errlist;

    fflush(stdout);

    fprintf(stdout, "\n");

    fprintf(stdout, _("\nWarning (libScotch):\n\n"));

    va_start(errlist, errstr);
    vfprintf(stdout, errstr, errlist);
    va_end(errlist);
    fprintf(stdout, "\n\n");
    fflush(stdout);
  }
}

/*----------------------------------------------------------------------------
 * Build cell -> cell connectivity
 *
 * parameters:
 *   n_cells        <-- number of cells in mesh
 *   n_faces        <-- number of cells in mesh
 *   start_cell     <-- number of first cell for the curent rank
 *   face_cells     <-- face->cells connectivity
 *   cell_idx       --> cell->cells index
 *   cell_neighbors --> cell->cells connectivity
 *----------------------------------------------------------------------------*/

static void
_scotch_cell_cells(size_t        n_cells,
                   size_t        n_faces,
                   cs_gnum_t     start_cell,
                   cs_gnum_t    *face_cells,
                   SCOTCH_Num  **cell_idx,
                   SCOTCH_Num  **cell_neighbors)
{
  size_t i;
  cs_gnum_t  id_0, id_1, c_num[2];
  SCOTCH_Num  start_id, end_id, c_id;

  SCOTCH_Num  *n_neighbors;
  SCOTCH_Num  *_cell_idx;
  SCOTCH_Num  *_cell_neighbors;

  /* Count and allocate arrays */

  BFT_MALLOC(n_neighbors, n_cells, SCOTCH_Num);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    assert(   c_num[0] - start_cell < n_cells
           || c_num[1] - start_cell < n_cells);

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells)
      n_neighbors[c_num[0] - start_cell] += 1;
    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells)
      n_neighbors[c_num[1] - start_cell] += 1;
  }

  BFT_MALLOC(_cell_idx, n_cells + 1, SCOTCH_Num);

  _cell_idx[0] = 0;

  for (i = 0; i < n_cells; i++)
    _cell_idx[i + 1] = _cell_idx[i] + n_neighbors[i];

  BFT_MALLOC(_cell_neighbors, _cell_idx[n_cells], SCOTCH_Num);

  for (i = 0; i < n_cells; i++)
    n_neighbors[i] = 0;

  for (i = 0; i < n_faces; i++) {

    c_num[0] = face_cells[i*2];
    c_num[1] = face_cells[i*2 + 1];

    if (c_num[0] == 0 || c_num[1] == 0 || c_num[0] == c_num[1])
      continue;

    if (c_num[0] >= start_cell && c_num[0] - start_cell < n_cells) {
      id_0 = c_num[0] - start_cell;
      _cell_neighbors[_cell_idx[id_0] + n_neighbors[id_0]] = c_num[1] - 1;
      n_neighbors[id_0] += 1;
    }

    if (c_num[1] >= start_cell && c_num[1] - start_cell < n_cells) {
      id_1 = c_num[1] - start_cell;
      _cell_neighbors[_cell_idx[id_1] + n_neighbors[id_1]] = c_num[0] - 1;
      n_neighbors[id_1] += 1;
    }
  }

  BFT_FREE(n_neighbors);

  /* Clean graph */

  c_id = 0;
  start_id = _cell_idx[0]; /* also = 0 */
  end_id = 0;

  for (i = 0; i < n_cells; i++) {

    SCOTCH_Num j, n_prev;

    end_id = _cell_idx[i+1];

    _scotch_sort_shell(start_id, end_id, _cell_neighbors);

    n_prev = _cell_neighbors[start_id];
    _cell_neighbors[c_id] = n_prev;
    c_id += 1;

    for (j = start_id + 1; j < end_id; j++) {
      if (_cell_neighbors[j] != n_prev) {
        n_prev = _cell_neighbors[j];
        _cell_neighbors[c_id] = n_prev;
        c_id += 1;
      }
    }

    start_id = end_id;
    _cell_idx[i+1] = c_id;

  }

  if (c_id < end_id)
    BFT_REALLOC(_cell_neighbors, c_id, SCOTCH_Num);

  /* Set return values */

  *cell_idx = _cell_idx;
  *cell_neighbors = _cell_neighbors;
}

/*----------------------------------------------------------------------------
 * Compute partition using SCOTCH
 *
 * parameters:
 *   n_cells       <-- number of cells in mesh
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *----------------------------------------------------------------------------*/

static void
_part_scotch(SCOTCH_Num   n_cells,
             int          n_parts,
             SCOTCH_Num  *cell_idx,
             SCOTCH_Num  *cell_neighbors,
             int         *cell_part)
{
  SCOTCH_Num  i;
  SCOTCH_Graph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;

  double  start_time, end_time;

  int     retval = 0;

  SCOTCH_Num    edgecut = 0; /* <-- Number of faces on partition */
  SCOTCH_Num  *_cell_part = NULL;

  /* Initialization */

  start_time = cs_timer_wtime();

  if (sizeof(SCOTCH_Num) == sizeof(int))
    _cell_part = (SCOTCH_Num *)cell_part;
  else
    BFT_MALLOC(_cell_part, n_cells, SCOTCH_Num);

  bft_printf(_("\n"
               " Partitioning %llu cells to %d domains\n"
               "   (SCOTCH_graphPart).\n"),
             (unsigned long long)n_cells, n_parts);

  /* Partition using libScotch */

  SCOTCH_graphInit(&grafdat);

  retval
    = SCOTCH_graphBuild(&grafdat,
                        0,                  /* baseval; 0 to n -1 numbering */
                        n_cells,            /* vertnbr */
                        cell_idx,           /* verttab */
                        NULL,               /* vendtab: verttab + 1 or NULL */
                        NULL,               /* velotab: vertex weights */
                        NULL,               /* vlbltab; vertex labels */
                        cell_idx[n_cells],  /* edgenbr */
                        cell_neighbors,     /* edgetab */
                        NULL);              /* edlotab */

  if (retval == 0) {

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_graphCheck(&grafdat) == 0)
      retval = SCOTCH_graphPart(&grafdat, n_parts, &stradat, _cell_part);

    SCOTCH_stratExit(&stradat);
  }

  SCOTCH_graphExit(&grafdat);

  /* Shift cell_part values to 1 to n numbering and free possible temporary */

  if (sizeof(SCOTCH_Num) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }

  /* Compute edge cut */

  if (retval == 0) {

    SCOTCH_Num cell_id, edgenum, commcut;

    commcut = 0;

    for (cell_id = 0; cell_id < n_cells; cell_id++) {
      SCOTCH_Num  edgennd,  partval;
      partval = cell_part[cell_id];
      for (edgenum = cell_idx[cell_id], edgennd = cell_idx[cell_id + 1];
           edgenum < edgennd;
           edgenum++) {
        if (cell_part[cell_neighbors[edgenum]] != partval)
          commcut++;
      }
    }

    edgecut = commcut / 2;
  }

  /* Finalization */

  end_time = cs_timer_wtime();

  bft_printf(_("\n"
               "  Total number of faces on parallel boundaries: %llu\n"
               "  wall-clock time: %f s\n\n"),
             (unsigned long long)edgecut,
             (double)(end_time - start_time));

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  SCOTCH_graphPart:           %.3g s\n",
                (double)(end_time - start_time));
}

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

#if defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Compute partition using PT-SCOTCH
 *
 * parameters:
 *   n_g_cells     <-- global number of cells
 *   cell_range    <-- first and past-the-last cell numbers for this rank
 *   n_parts       <-- number of partitions
 *   cell_cell_idx <-- cell->cells index
 *   cell_cell     <-- cell->cells connectivity
 *   cell_part     --> cell partition
 *   comm          <-- associated MPI communicator
 *----------------------------------------------------------------------------*/

static void
_part_ptscotch(cs_gnum_t    n_g_cells,
               cs_gnum_t    cell_range[2],
               int          n_parts,
               SCOTCH_Num  *cell_idx,
               SCOTCH_Num  *cell_neighbors,
               int         *cell_part,
               MPI_Comm     comm)
{
  int  n_ranks;
  SCOTCH_Num  i;
  SCOTCH_Dgraph  grafdat;  /* Scotch graph object to interface with libScotch */
  SCOTCH_Strat  stradat;

  double  start_time, end_time;

  int     retval = 0;

  SCOTCH_Num    n_cells = cell_range[1] - cell_range[0];
  SCOTCH_Num  *_cell_part = NULL;

  /* Initialization */

  start_time = cs_timer_wtime();

  MPI_Comm_size(comm, &n_ranks);

  if (sizeof(SCOTCH_Num) == sizeof(int))
    _cell_part = (SCOTCH_Num *)cell_part;
  else
    BFT_MALLOC(_cell_part, n_cells, SCOTCH_Num);

  bft_printf(_("\n"
               " Partitioning %llu cells to %d domains on %d ranks\n"
               "   (SCOTCH_dgraphPart).\n"),
             (unsigned long long)n_g_cells, n_parts, n_ranks);

  /* Partition using libScotch */

  retval = SCOTCH_dgraphInit(&grafdat, comm);

  if (retval == 0) {
    retval = SCOTCH_dgraphBuild
               (&grafdat,
                0,                  /* baseval; 0 to n -1 numbering */
                n_cells,            /* vertlocnbr */
                n_cells,            /* vertlocmax (= vertlocnbr) */
                cell_idx,           /* vertloctab */
                NULL,               /* vendloctab: vertloctab + 1 or NULL */
                NULL,               /* veloloctab: vertex weights */
                NULL,               /* vlblloctab; vertex labels */
                cell_idx[n_cells],  /* edgelocnbr */
                cell_idx[n_cells],  /* edgelocsiz */
                cell_neighbors,     /* edgeloctab */
                NULL,               /* edgegstab */
                NULL);              /* edloloctab */
  }

  if (retval == 0) {

    SCOTCH_stratInit(&stradat);

    if (SCOTCH_dgraphCheck(&grafdat) == 0)
      retval = SCOTCH_dgraphPart(&grafdat, n_parts, &stradat, _cell_part);

    SCOTCH_stratExit(&stradat);
  }

  SCOTCH_dgraphExit(&grafdat);

  /* Shift cell_part values to 1 to n numbering and free possible temporary */

  if (sizeof(SCOTCH_Num) != sizeof(int)) {
    for (i = 0; i < n_cells; i++)
      cell_part[i] = _cell_part[i];
    BFT_FREE(_cell_part);
  }

  /* Finalization */

  end_time = cs_timer_wtime();

  bft_printf(_("\n"
               "  wall-clock time: %f s\n\n"),
             (double)(end_time - start_time));

  cs_log_printf(CS_LOG_PERFORMANCE,
                "  SCOTCH_dgraphPart:          %.3g s\n",
                (double)(end_time - start_time));
}

#endif /* defined(HAVE_PTSCOTCH) */

/*----------------------------------------------------------------------------
 * Prepare input from mesh builder for use by partitioner.
 *
 * parameters:
 *   mesh           <-- pointer to mesh structure
 *   mb             <-- pointer to mesh builder structure
 *   rank_step      <-- Step between active partitioning ranks
 *                      (1 in basic case, > 1 if we seek to partition on a
 *                      reduced number of ranks)
 *   ignore_perio   <-- ignore periodicity information if true
 *   cell_range     <-- first and past-the-last cell numbers for this rank
 *   n_faces        <-- number of local faces for current rank
 *   g_face_cells   <-> global face -> cells connectivity
 *----------------------------------------------------------------------------*/

static void
_prepare_input(const cs_mesh_t           *mesh,
               const cs_mesh_builder_t   *mb,
               int                        rank_step,
               bool                       ignore_perio,
               cs_gnum_t                  cell_range[2],
               cs_lnum_t                 *n_faces,
               cs_gnum_t                **g_face_cells)
{
  int rank_id = cs_glob_rank_id;
  int n_ranks = cs_glob_n_ranks;

  cs_gnum_t *_g_face_cells = NULL;

#if defined(HAVE_MPI)
  cs_block_to_part_t *d = NULL;
#endif

  cs_block_dist_info_t cell_bi =  cs_block_dist_compute_sizes(rank_id,
                                                              n_ranks,
                                                              rank_step,
                                                              0,
                                                              mesh->n_g_cells);

  /* By default, the face -> cells connectivity is that of the mesh builder */

  *g_face_cells = mb->face_cells;

  /* In case of periodicity, update global face -> cells connectivity:
     if face Fi is periodic with face Fj, and face Fi is connected
     to cell Ci, while face Fj is connected to face Cj, we add Cj
     to the cells connected to Fi, and Ci to the cells connected to Fj,
     just as if Fi and Fj were regular interior faces (this ignores
     the geometric transformation, but is sufficient to build the
     cells -> cells graph for domain partitioning). */

  if (ignore_perio == false && mb->n_perio > 0) {

    int perio_id;
    cs_lnum_t n_per_face_couples = 0;
    cs_gnum_t n_g_per_face_couples = 0;
    cs_gnum_t _n_b_faces = 0;
    cs_gnum_t *_per_face_couples = NULL;
    const cs_gnum_t *per_face_couples = NULL;

    bft_printf(_("\n"
                 " Ignoring periodicity for graph-based partitioning.\n"));

    if (mb->face_bi.gnum_range[1] > mb->face_bi.gnum_range[0])
      _n_b_faces = mb->face_bi.gnum_range[1] - mb->face_bi.gnum_range[0];

    BFT_MALLOC(_g_face_cells, _n_b_faces*2, cs_gnum_t);
    *g_face_cells = _g_face_cells;

    memcpy(_g_face_cells, mb->face_cells, sizeof(cs_gnum_t)*_n_b_faces*2);

    if (mb->n_perio > 1) {

      /* Assemble faces from multiple periodicities into one */

      for (perio_id = 0; perio_id < mb->n_perio; perio_id++) {
        n_per_face_couples += mb->n_per_face_couples[perio_id];
        n_g_per_face_couples += mb->n_g_per_face_couples[perio_id];
      }

      BFT_MALLOC(_per_face_couples, n_per_face_couples*2, cs_gnum_t);
      per_face_couples = _per_face_couples;

      n_per_face_couples = 0;

      for (perio_id = 0; perio_id < mb->n_perio; perio_id++) {
        cs_gnum_t *per_face_couples_p = _per_face_couples + 2*n_per_face_couples;
        memcpy(per_face_couples_p,
               mb->per_face_couples[perio_id],
               sizeof(cs_gnum_t)*mb->n_per_face_couples[perio_id]*2);
        n_per_face_couples += mb->n_per_face_couples[perio_id];
      }

    }
    else { /* if mb->n_perio == 1 */
      n_per_face_couples = mb->n_per_face_couples[0];
      n_g_per_face_couples = mb->n_g_per_face_couples[0];
      per_face_couples = mb->per_face_couples[0];
    }

    if (n_g_per_face_couples > 0) {

#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1)
        _add_perio_to_face_cells_g(mb->face_bi,
                                   _g_face_cells,
                                   n_per_face_couples,
                                   per_face_couples,
                                   cs_glob_mpi_comm);
#endif

      if (n_ranks == 1)
        _add_perio_to_face_cells_l(_g_face_cells,
                                   n_per_face_couples,
                                   per_face_couples);

    }

    if (_per_face_couples != NULL)
      BFT_FREE(_per_face_couples);
  }

  /* Distribute faces if necessary */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_datatype_t gnum_type
      = (sizeof(cs_gnum_t) == 8) ? CS_UINT64 : CS_UINT32;
    cs_gnum_t *g_face_cells_tmp = NULL;

    d = cs_block_to_part_create_by_adj_s(cs_glob_mpi_comm,
                                         mb->face_bi,
                                         cell_bi,
                                         2,
                                         *g_face_cells,
                                         NULL,
                                         NULL);

    *n_faces = cs_block_to_part_get_n_part_ents(d);

    BFT_MALLOC(g_face_cells_tmp, (*n_faces)*2, cs_gnum_t);

    /* Face -> cell connectivity */

    cs_block_to_part_copy_array(d,
                                gnum_type,
                                2,
                                *g_face_cells,
                                g_face_cells_tmp);

    if (_g_face_cells != NULL) /* in case of periodicity */
      BFT_FREE(_g_face_cells);

    *g_face_cells = g_face_cells_tmp;

    cs_block_to_part_destroy(&d);
  }

#endif /* defined(HAVE_MPI) */

  if (cs_glob_n_ranks == 1)
    *n_faces = mb->n_g_faces;

  /* Prepare auxiliary return values */

  cell_range[0] = cell_bi.gnum_range[0];
  cell_range[1] = cell_bi.gnum_range[1];
}

#if   defined(HAVE_METIS) || defined(HAVE_PARMETIS) \
   || defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

/*----------------------------------------------------------------------------
 * Distribute partitioning info so as to match mesh builder block info.
 *
 * parameters:
 *   mb           <-- pointer to mesh builder structure
 *   rank_step    <-- Step between active partitioning ranks
 *                    (1 in basic case, > 1 if we seek to partition on a
 *                    reduced number of ranks)
 *   cell_range   <-- first and past-the-last cell numbers for this rank
 *   cell_rank    <-> pointer to pointer to cell rank
 *----------------------------------------------------------------------------*/

static void
_distribute_output(const cs_mesh_builder_t   *mb,
                   int                        rank_step,
                   const cs_gnum_t            cell_range[2],
                   int                      **cell_rank)
{
  /* Distribute partitioning info if necessary */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1 && (mb->cell_bi.rank_step != rank_step)) {

    cs_gnum_t i;
    cs_gnum_t n_b_cells = 0, n_p_cells = 0;
    cs_datatype_t int_type = (sizeof(int) == 8) ? CS_INT64 : CS_INT32;

    cs_part_to_block_t *d = NULL;
    cs_gnum_t *global_cell_num = NULL;

    int *_cell_rank = NULL;

    if (mb->cell_bi.gnum_range[1] > mb->cell_bi.gnum_range[0])
      n_b_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

    if (cell_range[1] > cell_range[0])
      n_p_cells = cell_range[1] - cell_range[0];

    BFT_MALLOC(_cell_rank, n_b_cells, int);
    BFT_MALLOC(global_cell_num, n_p_cells, cs_gnum_t);

    for (i = 0; i < n_p_cells; i++)
      global_cell_num[i] = cell_range[0] + i;

    d = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                        mb->cell_bi,
                                        n_p_cells,
                                        global_cell_num);
    cs_part_to_block_transfer_gnum(d, global_cell_num);
    global_cell_num = NULL;

    /* Cell rank id */

    cs_part_to_block_copy_array(d,
                                int_type,
                                1,
                                *cell_rank,
                                _cell_rank);

    cs_part_to_block_destroy(&d);

    BFT_FREE(*cell_rank);
    *cell_rank = _cell_rank;
  }

#endif /* defined(HAVE_MPI) */
}

#endif /*    defined(HAVE_METIS) || defined(HAVE_PARMETIS) \
          || defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

/*----------------------------------------------------------------------------
 * Write output file.
 *
 * parameters:
 *   n_g_cells   <-- global number of cells
 *   cell_range  <-- first and past-the-last cell numbers for this rank
 *   n_ranks     <-- number of ranks corresonding to output file
 *   domain_rank <-- domain rank array (size: n_cells)
 *----------------------------------------------------------------------------*/

static void
_write_output(cs_gnum_t  n_g_cells,
              cs_gnum_t  cell_range[2],
              int        n_ranks,
              const int  domain_rank[])
{
  size_t i;
  int n_ranks_size;
  cs_file_access_t method;
  char *filename = NULL;
  cs_io_t *fh = NULL;
  int *domain_num = NULL;
  cs_datatype_t datatype_gnum = CS_DATATYPE_NULL;
  cs_datatype_t datatype_int = CS_DATATYPE_NULL;

  const char dir[] = "partition_output";
  const char magic_string[] = "Domain partitioning, R0";

  if (sizeof(int) == 4)
    datatype_int = CS_INT32;
  else if (sizeof(int) == 8)
    datatype_int = CS_INT64;
  else {
    assert(0);
  }

  if (sizeof(cs_gnum_t) == 4)
    datatype_gnum = CS_UINT32;
  else if (sizeof(cs_gnum_t) == 8)
    datatype_gnum = CS_UINT64;
  else {
    assert(0);
  }

  if (cell_range[1] > cell_range[0]) {
    cs_gnum_t _n_cells = cell_range[1] - cell_range[0];
    BFT_MALLOC(domain_num, _n_cells, int);
    for (i = 0; i < _n_cells; i++)
      domain_num[i] = domain_rank[i] + 1;
  }

  /* Create directory if required */

  if (cs_glob_rank_id < 1) {
    if (cs_file_isdir(dir) != 1) {
      if (cs_file_mkdir_default(dir) != 0)
        bft_error(__FILE__, __LINE__, errno,
                  _("The partitioning directory cannot be created"));
    }
  }

  /* Open file */

  for (i = n_ranks, n_ranks_size = 1;
       i >= 10;
       i /= 10, n_ranks_size += 1);

  BFT_MALLOC(filename,
             strlen(dir) + strlen("domain_number_") + n_ranks_size + 2,
             char);

  sprintf(filename,
          "%s%cdomain_number_%d",
          dir, _dir_separator, n_ranks);

#if defined(HAVE_MPI)
  {
    MPI_Info  hints;
    MPI_Comm  block_comm, comm;
    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method, &hints);
    cs_file_get_default_comm(NULL, NULL, &block_comm, &comm);
    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);
    fh = cs_io_initialize(filename,
                          magic_string,
                          CS_IO_MODE_WRITE,
                          method,
                          CS_IO_ECHO_OPEN_CLOSE,
                          hints,
                          block_comm,
                          comm);
  }
#else
  {
    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method);
    fh = cs_io_initialize(filename,
                          magic_string,
                          CS_IO_MODE_WRITE,
                          method,
                          CS_IO_ECHO_OPEN_CLOSE);
  }
#endif

  BFT_FREE(filename);

  /* Write headers */

  cs_io_write_global("n_cells",
                     1,
                     1,
                     0,
                     1,
                     datatype_gnum,
                     &n_g_cells,
                     fh);

  cs_io_write_global("n_ranks",
                     1,
                     1,
                     0,
                     1,
                     datatype_int,
                     &n_ranks,
                     fh);

  cs_io_write_block_buffer("cell:domain number",
                           n_g_cells,
                           cell_range[0],
                           cell_range[1],
                           1,
                           0,
                           1,
                           datatype_int,
                           domain_num,
                           fh);

  cs_io_finalize(&fh);

  BFT_FREE(domain_num);
}

/*----------------------------------------------------------------------------
 * Read cell rank if available
 *
 * parameters:
 *   mesh <-- pointer to mesh structure
 *   mb   <-> pointer to mesh builder helper structure
 *   echo <-- echo (verbosity) level
 *----------------------------------------------------------------------------*/

static void
_read_cell_rank(cs_mesh_t          *mesh,
                cs_mesh_builder_t  *mb,
                long                echo)
{
  char file_name[64]; /* more than enough for
                         "partition_input/domain_number_<n_ranks>" */
  size_t  i;
  cs_file_access_t  method;
  cs_io_sec_header_t  header;

  cs_io_t  *rank_pp_in = NULL;
  cs_lnum_t   n_ranks = 0;
  cs_gnum_t   n_elts = 0;
  cs_gnum_t   n_g_cells = 0;

  const char magic_string[] = "Domain partitioning, R0";
  const char  *unexpected_msg = N_("Section of type <%s> on <%s>\n"
                                   "unexpected or of incorrect size");

  if (n_ranks == 1)
    return;

#if (__STDC_VERSION__ < 199901L)
  sprintf(file_name,
          "partition_input%cdomain_number_%d",
          _dir_separator, cs_glob_n_ranks);
#else
  snprintf(file_name, 64,
           "partition_input%cdomain_number_%d",
           _dir_separator, cs_glob_n_ranks);
#endif
  file_name[63] = '\0'; /* Just in case; processor counts would need to be
                           in the exa-range for this to be necessary. */

  /* Test if file exists */

  if (! cs_file_isreg(file_name)) {
    bft_printf(_(" No \"%s\" file available;\n"), file_name);
    return;
  }

  /* Open file */

#if defined(HAVE_MPI)
  {
    MPI_Info           hints;
    MPI_Comm           block_comm, comm;
    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method, &hints);
    cs_file_get_default_comm(NULL, NULL, &block_comm, &comm);
    assert(comm == cs_glob_mpi_comm || comm == MPI_COMM_NULL);
    rank_pp_in = cs_io_initialize(file_name,
                                  magic_string,
                                  CS_IO_MODE_READ,
                                  method,
                                  echo,
                                  hints,
                                  block_comm,
                                  comm);
  }
#else
  {
    cs_file_get_default_access(CS_FILE_MODE_WRITE, &method);
    rank_pp_in = cs_io_initialize(file_name,
                                  magic_string,
                                  CS_IO_MODE_READ,
                                  method,
                                  echo);
  }
#endif

  if (echo > 0)
    bft_printf("\n");

  /* Loop on read sections */

  while (rank_pp_in != NULL) {

    /* Receive headers */

    cs_io_read_header(rank_pp_in, &header);

    /* Treatment according to the header name */

    if (strncmp(header.sec_name, "n_cells",
                CS_IO_NAME_LEN) == 0) {

      if (header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_cs_gnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_g_cells, rank_pp_in);
        if (n_g_cells != mesh->n_g_cells)
          bft_error(__FILE__, __LINE__, 0,
                    _("The number of cells reported by file\n"
                      "\"%s\" (%llu)\n"
                      "does not correspond to those of the mesh (%llu)."),
                    cs_io_get_name(rank_pp_in),
                    (unsigned long long)(n_g_cells),
                    (unsigned long long)(mesh->n_g_cells));
      }

    }
    else if (strncmp(header.sec_name, "n_ranks",
                     CS_IO_NAME_LEN) == 0) {

      if (header.n_vals != 1)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        cs_io_set_cs_lnum(&header, rank_pp_in);
        cs_io_read_global(&header, &n_ranks, rank_pp_in);
        if (n_ranks != cs_glob_n_ranks)
          bft_error(__FILE__, __LINE__, 0,
                    _("The number of ranks reported by file\n"
                      "\"%s\" (%d) does not\n"
                      "correspond to the current number of ranks (%d)."),
                    cs_io_get_name(rank_pp_in), (int)n_ranks,
                    (int)cs_glob_n_ranks);
      }

    }
    else if (strncmp(header.sec_name, "cell:domain number",
                     CS_IO_NAME_LEN) == 0) {

      n_elts = mesh->n_g_cells;
      if (header.n_vals != (cs_file_off_t)n_elts)
        bft_error(__FILE__, __LINE__, 0,
                  _(unexpected_msg), header.sec_name,
                  cs_io_get_name(rank_pp_in));
      else {
        mb->have_cell_rank = true;
        cs_io_set_cs_lnum(&header, rank_pp_in);
        if (mb->cell_bi.gnum_range[0] > 0)
          n_elts = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];
        BFT_MALLOC(mb->cell_rank, n_elts, cs_lnum_t);
        cs_io_read_block(&header,
                         mb->cell_bi.gnum_range[0],
                         mb->cell_bi.gnum_range[1],
                         mb->cell_rank, rank_pp_in);
        for (i = 0; i < n_elts; i++) /* Convert 1 to n to 0 to n-1 */
          mb->cell_rank[i] -= 1;
      }
      cs_io_finalize(&rank_pp_in);
      rank_pp_in = NULL;

    }

    else
      bft_error(__FILE__, __LINE__, 0,
                _("Section of type <%s> on <%s> is unexpected."),
                header.sec_name, cs_io_get_name(rank_pp_in));
  }

  if (rank_pp_in != NULL)
    cs_io_finalize(&rank_pp_in);
}

/*----------------------------------------------------------------------------*
 * Define a naive partitioning by blocks.
 *
 * parameters:
 *   mesh      <-- pointer to mesh structure
 *   mb        <-- pointer to mesh builder structure
 *   cell_part --> assigned cell partition
 *----------------------------------------------------------------------------*/

static void
_block_partititioning(const cs_mesh_t          *mesh,
                      const cs_mesh_builder_t  *mb,
                      int                      *cell_part)
{
  cs_lnum_t i;

  int  n_ranks = cs_glob_n_ranks;
  cs_lnum_t block_size = mesh->n_g_cells / n_ranks;
  cs_lnum_t n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

  if (mesh->n_g_cells % n_ranks)
    block_size += 1;

  /* Use variable block size if necessary so that all ranks have some
     cells assigned if possible */

  if ((cs_lnum_t)((mesh->n_g_cells - 1) % block_size) < n_ranks - 1) {

    cs_gnum_t cells_per_rank = mesh->n_g_cells / n_ranks;
    cs_lnum_t rmdr = mesh->n_g_cells - cells_per_rank * (cs_gnum_t)n_ranks;

    if (rmdr == 0) {
      for (i = 0; i < n_cells; i++) {
        cs_gnum_t cell_num = mb->cell_bi.gnum_range[0] + i;
        cell_part[i] = (cell_num - 1) / cells_per_rank;
      }
    }
    else {
      cs_gnum_t n_ranks_rmdr = n_ranks - rmdr;
      cs_gnum_t n_ranks_cells_per_rank = n_ranks_rmdr * cells_per_rank;
      for (i = 0; i < n_cells; i++) {
        cs_gnum_t cell_num = mb->cell_bi.gnum_range[0] + i;
        if ((cell_num - 1)  <  n_ranks_cells_per_rank)
          cell_part[i] = (cell_num - 1) / cells_per_rank;
        else
          cell_part[i] = (cell_num + n_ranks_rmdr - 1) / (cells_per_rank + 1);
      }
    }

  }

  else {
    for (i = 0; i < n_cells; i++) {
      cs_gnum_t cell_num = mb->cell_bi.gnum_range[0] + i;
      cell_part[i] = ((cell_num - 1) / block_size);
    }
  }

}

/*----------------------------------------------------------------------------*
 * Indicate if a usable graph partitioning algorithm is available
 * for a given partitioning stage.
 *
 * If a non-default partitioning algorithm for this stage has ben defined,
 * it is considered usable here (checking may be done later).
 *
 * parameters:
 *   stage <-- associated partitioning stage
 *
 * returns:
 *   true if graph partitioning is available, false otherwise.
 *----------------------------------------------------------------------------*/

static bool
_have_usable_graph_partitioning(cs_partition_stage_t  stage)
{
  bool retval = false;

  cs_partition_algorithm_t a = _part_algorithm[stage];

  int n_part_ranks = cs_glob_n_ranks / _part_rank_step[stage];
  if (n_part_ranks < 1)
    n_part_ranks = 1;

  if (a >= CS_PARTITION_SCOTCH && a <= CS_PARTITION_METIS)
    retval = true;

#if defined(HAVE_PTSCOTCH)
  if (a == CS_PARTITION_DEFAULT)
    retval = true;
#elif defined(HAVE_SCOTCH)
  if (n_part_ranks == 1 && a == CS_PARTITION_DEFAULT)
    retval = true;
#elif defined(HAVE_PARMETIS)
  if (a == CS_PARTITION_DEFAULT)
    retval = true;
#elif defined(HAVE_METIS)
  if (n_part_ranks == 1 && a == CS_PARTITION_DEFAULT)
    retval = true;
#endif

  return retval;
}

/*----------------------------------------------------------------------------*
 * Return default algorithm for domain partitioning.
 *
 * parameters:
 *   stage <-- associated partitioning stage
 *
 * returns:
 *   default partitioning algorithm for this stage
 *----------------------------------------------------------------------------*/

static cs_partition_algorithm_t
_select_algorithm(cs_partition_stage_t  stage)
{
  cs_partition_algorithm_t retval = _part_algorithm[stage];

  if (retval == CS_PARTITION_DEFAULT) {

    int n_part_ranks;

    n_part_ranks = cs_glob_n_ranks / _part_rank_step[stage];
    if (n_part_ranks < 1)
      n_part_ranks = 1;

    /* last stage of 1 or 2 */

#if defined(HAVE_PTSCOTCH)
    if (retval == CS_PARTITION_DEFAULT)
      retval = CS_PARTITION_SCOTCH;
#elif defined(HAVE_SCOTCH)
    if (n_part_ranks == 1 && retval == CS_PARTITION_DEFAULT)
      retval = CS_PARTITION_SCOTCH;
#elif defined(HAVE_PARMETIS)
    retval = CS_PARTITION_METIS;
#elif defined(HAVE_METIS)
    if (n_part_ranks == 1 && retval == CS_PARTITION_DEFAULT)
      retval = CS_PARTITION_METIS;
#endif
    if (retval == CS_PARTITION_DEFAULT)
      retval = CS_PARTITION_SFC_MORTON_BOX;

    /* 1st stage of 2:
       If 2nd stage uses a space-filling curve, use same curve by default;
       otherwise, use default space-filling curve. */

    if (stage == CS_PARTITION_FOR_PREPROCESS) {

      if (   _part_algorithm[CS_PARTITION_MAIN] >= CS_PARTITION_SFC_MORTON_BOX
          && _part_algorithm[CS_PARTITION_MAIN] <= CS_PARTITION_SFC_HILBERT_CUBE)
        retval = _part_algorithm[CS_PARTITION_MAIN];
      else
        retval = CS_PARTITION_SFC_MORTON_BOX;

    }

  }

  return retval;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print information on external libraries
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_external_library_info(void)
{
  int  n_ext_libs = 0;

  if (cs_glob_rank_id >= 1)
    return;

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  n_ext_libs++;
#endif
#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)
  n_ext_libs++;
#endif

  if (n_ext_libs == 1)
    bft_printf(_("\n  External libraries for partitioning:\n"));
  else if (n_ext_libs > 1)
    bft_printf(_("\n  External libraries for partitioning:\n"));

#if defined(HAVE_METIS)
#if defined(HAVE_METIS_H) && defined(METIS_VER_MAJOR)
  bft_printf("    METIS %d.%d.%d\n",
             METIS_VER_MAJOR, METIS_VER_MINOR, METIS_VER_SUBMINOR);
#else
  bft_printf("    METIS\n");
#endif
#elif defined(HAVE_PARMETIS)
#if defined(PARMETIS_MAJOR_VERSION) && defined(PARMETIS_SUBMINOR_VERSION)
  bft_printf("    ParMETIS %d.%d.%d\n",
             PARMETIS_MAJOR_VERSION, PARMETIS_MINOR_VERSION,
             PARMETIS_SUBMINOR_VERSION);
#elif defined(PARMETIS_MAJOR_VERSION)
  bft_printf("    ParMETIS %d.%d\n",
             PARMETIS_MAJOR_VERSION, PARMETIS_MINOR_VERSION);
#else
  bft_printf("    ParMETIS\n");
#endif
#endif

#if defined(HAVE_SCOTCH)
#if defined(SCOTCH_VERSION) && defined(SCOTCH_RELEASE)
  bft_printf("    SCOTCH %d.%d.%d\n",
             SCOTCH_VERSION, SCOTCH_RELEASE, SCOTCH_PATCHLEVEL);
#else
  bft_printf("    SCOTCH\n");
#endif
#elif defined(HAVE_PTSCOTCH)
#if defined(SCOTCH_VERSION) && defined(SCOTCH_RELEASE)
  bft_printf("    PT-SCOTCH %d.%d.%d\n",
             SCOTCH_VERSION, SCOTCH_RELEASE, SCOTCH_PATCHLEVEL);
#else
  bft_printf("    PT-SCOTCH\n");
#endif
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set algorithm for domain partitioning.
 *
 * \param[in]  stage         associated partitioning stage
 * \param[in]  algorithm     partitioning algorithm choice
 * \param[in]  rank_step     if > 1, partitioning done on at most
 *                           n_ranks / rank_step processes
 *                           (for graph-based partitioning only)
 * \param[in]  ignore_perio  if true, ignore periodicity information when
 *                           present (for graph-based partitioning only)
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_set_algorithm(cs_partition_stage_t      stage,
                           cs_partition_algorithm_t  algorithm,
                           int                       rank_step,
                           bool                      ignore_perio)
{
  int n_part_ranks = cs_glob_n_ranks / rank_step;

  if (n_part_ranks < 1) {
    rank_step = cs_glob_n_ranks;
    n_part_ranks = 1;
  }

  /* Check consistency of choice */

#if !defined(HAVE_PTSCOTCH)
  if (algorithm == CS_PARTITION_SCOTCH) {
#if defined(HAVE_SCOTCH)
    if (n_part_ranks > 1) {
      rank_step = cs_glob_n_ranks;
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(
        _("Partitioning with %s requested, but %s is not available,\n"
           "so serial %s will be used."), "LibSCOTCH", "PT-SCOTCH", "SCOTCH");
    }
#else
    bft_error(__FILE__, __LINE__, 0,
              _("Partitioning with %s required but neither\n"
                "%s nor %s is available."), "LibSCOTCH", "PT-SCOTCH", "SCOTCH");
#endif
  }
#endif /* defined(HAVE_PTSCOTCH) */

#if !defined(HAVE_PARMETIS)
  if (algorithm == CS_PARTITION_METIS) {
#if defined(HAVE_METIS)
    if (n_part_ranks > 1) {
      rank_step = cs_glob_n_ranks;
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(
         _("Partitioning with %s requested, but %s is not available,\n"
           "so serial %s will be used."), "METIS", "ParMETIS", "METIS");
    }
#else
    bft_error(__FILE__, __LINE__, 0,
              _("Partitioning with %s required but neither\n"
                "%s nor %s is available."), "METIS", "ParMETIS", "METIS");
#endif
  }
#endif /* defined(HAVE_PARMETIS) */

  /* Set options */

  _part_algorithm[stage] = algorithm;
  _part_rank_step[stage] = rank_step;
  _part_ignore_perio[stage] = ignore_perio;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set partitioning write to file option.
 *
 * Partitioning information for subsequent calculations is written to file
 * after the last partitioning stage depending on the output level.
 *
 * Note that partitioning information for additional partitionings is
 * always written to file, regardless of this option.
 *
 * \param[in]  write_flag  option to save partitioning information:
 *                         0: never
 *                         1: for graph-based partitioning only
 *                         2: always
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_set_write_level(int  write_flag)
{
  _part_write_output = write_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define hints indicating if initial partitioning fo a preprocessing
 *        stage is required.
 *
 * \param[in]  join           true if a mesh joining operation is planned.
 * \param[in]  join_periodic  true if a mesh periodic matching operation
 *                            is planned
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_set_preprocess_hints(bool  join,
                                  bool  join_periodic)
{
  _part_compute_join_hint = join;
  _part_compute_perio_hint = join_periodic;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate or deactivate initial partitioning for preprocessing.
 *
 * \param[in]  active  true to activate pre-partitiong for the preprocessing
 *                     stage, false to de-activate it
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_set_preprocess(bool  active)
{
  if (active)
    _part_preprocess_active = 2;
  else
    _part_preprocess_active = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Indicate if re-partitiong for the computation stage is required.
 *
 * \return  true if initial partitioning for preprocessing is active,
 *          false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_partition_get_preprocess(void)
{
  bool retval = false;

  if (_part_preprocess_active < 1)
    retval = false;

  else if (_part_preprocess_active == 1) {

    if (_have_usable_graph_partitioning(CS_PARTITION_MAIN)) {
      if (_part_compute_join_hint)
        retval = true;
      if (   _part_compute_perio_hint
          && _part_ignore_perio[CS_PARTITION_MAIN] == false)
        retval = true;
    }
  }

  else
    retval = true;

  if (cs_glob_n_ranks < 2)
    retval = false;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define list of extra partitionings to build.
 *
 * Partitionings in this list will be output to file, and may be used for
 * subsequent calculations.
 *
 * When partitioning for both preprocessing and calculation stages, output to
 * file of partioning data or generation of additional partitionings
 * (see \ref cs_partition_add_partitions) will only be done for the
 * second stage.
 *
 * \param[in]  n_extra_partitions     number of extra partitionings to compute
 * \param[in]  extra_partitions_list  list of extra partitions to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_partition_add_partitions(int  n_extra_partitions,
                            int  extra_partitions_list[])
{
  _part_n_extra_partitions = n_extra_partitions;

  BFT_REALLOC(_part_extra_partitions_list, n_extra_partitions, int);

  if (n_extra_partitions > 0)
    memcpy(_part_extra_partitions_list,
           extra_partitions_list,
           sizeof(int)*n_extra_partitions);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Partition mesh based on current options.
 *
 * \param[in]       mesh   pointer to mesh structure
 * \param[in, out]  mb     pointer to mesh builder structure
 * \param[in]       stage  associated partitioning stage
 */
/*----------------------------------------------------------------------------*/

void
cs_partition(cs_mesh_t             *mesh,
             cs_mesh_builder_t     *mb,
             cs_partition_stage_t   stage)
{
  cs_timer_t t0, t1;
  cs_timer_counter_t dt;

  int n_part_ranks = cs_glob_n_ranks / _part_rank_step[stage];
  cs_partition_algorithm_t _algorithm = _select_algorithm(stage);

  bool write_output = false;
  int  n_extra_partitions = 0;

  int  *cell_part = NULL;

  cs_gnum_t  cell_range[2] = {0, 0};
  cs_lnum_t  n_cells = 0;
  cs_lnum_t  n_faces = 0;
  cs_gnum_t  *face_cells = NULL;

  /* Initialize local options */

  if (stage == CS_PARTITION_MAIN) {
    if (   (   _algorithm == CS_PARTITION_METIS
            || _algorithm == CS_PARTITION_SCOTCH)
        && _part_write_output > 0)
      write_output = true;
    else if (_part_write_output > 1)
      write_output = true;
    n_extra_partitions = _part_n_extra_partitions;
  }

  if (n_part_ranks < 1)
    n_part_ranks = 1;

  /* Free previous cell rank info if present */

  if (mb->cell_rank != NULL)
    BFT_FREE(mb->cell_rank);

  /* Read cell rank data if available */

  if (cs_glob_n_ranks > 1) {
    if (   stage != CS_PARTITION_MAIN
        || cs_partition_get_preprocess() == false) {
      _read_cell_rank(mesh, mb, CS_IO_ECHO_OPEN_CLOSE);
      if (mb->have_cell_rank)
        return;
    }
  }
  else { /* if (cs_glob_n_ranks == 1) */
    if (stage != CS_PARTITION_MAIN || n_extra_partitions < 1)
      return;
  }

  (void)cs_timer_wtime();

  /* Print header */

  bft_printf("\n ----------------------------------------------------------\n");

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Partitioning:\n\n"));

  /* Start timer */

  t0 = cs_timer_time();

  /* Adapt builder data for partitioning */

  if (_algorithm == CS_PARTITION_METIS || _algorithm == CS_PARTITION_SCOTCH) {

    _prepare_input(mesh,
                   mb,
                   _part_rank_step[stage],
                   _part_ignore_perio[stage],
                   cell_range,
                   &n_faces,
                   &face_cells);

    n_cells = cell_range[1] - cell_range[0];

  }
  else {

    n_part_ranks = cs_glob_n_ranks;
    n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

  }

  /* Build and partition graph */

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)

  if (_algorithm == CS_PARTITION_METIS) {

    int  i;
    cs_timer_t  t2;
    idx_t  *cell_idx = NULL, *cell_neighbors = NULL;

    _metis_cell_cells(n_cells,
                      n_faces,
                      cell_range[0],
                      face_cells,
                      &cell_idx,
                      &cell_neighbors);

    if (face_cells != mb->face_cells)
      BFT_FREE(face_cells);

    t2 = cs_timer_time();
    dt = cs_timer_diff(&t0, &t2);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  preparing graph:            %.3g s\n"),
                  (double)(dt.wall_nsec)/1.e9);

#if defined(HAVE_PARMETIS)

    if (n_part_ranks > 1) {

      MPI_Comm part_comm = cs_glob_mpi_comm;

      if (_part_rank_step[stage] > 1)
        part_comm = _init_reduced_communicator(_part_rank_step[stage]);

      for (i = 0; i < n_extra_partitions + 1; i++) {

        int  n_ranks = cs_glob_n_ranks;

        if (i < n_extra_partitions) {
          n_ranks = _part_extra_partitions_list[i];
          if (n_ranks == cs_glob_n_ranks) {
            write_output = true;
            continue;
          }
        }

        if (n_ranks < 2)
          continue;

        BFT_REALLOC(cell_part, n_cells, int);

        if (cs_glob_rank_id % _part_rank_step[stage] == 0)
          _part_parmetis(mesh->n_g_cells,
                         cell_range,
                         n_ranks,
                         cell_idx,
                         cell_neighbors,
                         cell_part,
                         part_comm);

        _distribute_output(mb,
                           _part_rank_step[stage],
                           cell_range,
                           &cell_part);

        _cell_part_histogram(mb->cell_bi.gnum_range, n_ranks, cell_part);

        if (write_output || i < n_extra_partitions)
          _write_output(mesh->n_g_cells,
                        mb->cell_bi.gnum_range,
                        n_ranks,
                        cell_part);
      }

      if (part_comm != cs_glob_mpi_comm && part_comm != MPI_COMM_NULL)
        MPI_Comm_free(&part_comm);
    }

#endif

    if (n_part_ranks == 1) {

      for (i = 0; i < n_extra_partitions + 1; i++) {

        int  n_ranks = cs_glob_n_ranks;

        if (i < n_extra_partitions) {
          n_ranks = _part_extra_partitions_list[i];
          if (n_ranks == cs_glob_n_ranks) {
            write_output = true;
            continue;
          }
        }

        if (n_ranks < 2)
          continue;

        BFT_REALLOC(cell_part, n_cells, int);

        if (cs_glob_rank_id < 0 || (cs_glob_rank_id % _part_rank_step[stage] == 0))
          _part_metis(n_cells,
                      n_ranks,
                      cell_idx,
                      cell_neighbors,
                      cell_part);

        _distribute_output(mb,
                           _part_rank_step[stage],
                           cell_range,
                           &cell_part);

        _cell_part_histogram(mb->cell_bi.gnum_range, n_ranks, cell_part);

        if (write_output || i < n_extra_partitions)
          _write_output(mesh->n_g_cells,
                        mb->cell_bi.gnum_range,
                        n_ranks,
                        cell_part);
      }
    }

    BFT_FREE(cell_idx);
    BFT_FREE(cell_neighbors);
  }

#endif /* defined(HAVE_METIS) || defined(HAVE_PARMETIS) */

#if defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH)

  if (_algorithm == CS_PARTITION_SCOTCH) {

    int  i;
    cs_timer_t  t2;
    SCOTCH_Num  *cell_idx = NULL, *cell_neighbors = NULL;

    _scotch_cell_cells(n_cells,
                       n_faces,
                       cell_range[0],
                       face_cells,
                       &cell_idx,
                       &cell_neighbors);

    if (face_cells != mb->face_cells)
      BFT_FREE(face_cells);

    t2 = cs_timer_time();
    dt = cs_timer_diff(&t0, &t2);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  preparing graph:            %.3g s\n"),
                  (double)(dt.wall_nsec)/1.e9);

#if defined(HAVE_PTSCOTCH)

    if (n_part_ranks > 1) {

      MPI_Comm part_comm = cs_glob_mpi_comm;

      if (_part_rank_step[stage] > 1)
        part_comm = _init_reduced_communicator(_part_rank_step[stage]);

      for (i = 0; i < n_extra_partitions + 1; i++) {

        int  n_ranks = cs_glob_n_ranks;

        if (i < n_extra_partitions) {
          n_ranks = _part_extra_partitions_list[i];
          if (n_ranks == cs_glob_n_ranks) {
            write_output = true;
            continue;
          }
        }

        if (n_ranks < 2)
          continue;

        BFT_REALLOC(cell_part, n_cells, int);

        if (cs_glob_rank_id % _part_rank_step[stage] == 0)
          _part_ptscotch(mesh->n_g_cells,
                         cell_range,
                         n_ranks,
                         cell_idx,
                         cell_neighbors,
                         cell_part,
                         part_comm);

        _distribute_output(mb,
                           _part_rank_step[stage],
                           cell_range,
                           &cell_part);

        _cell_part_histogram(mb->cell_bi.gnum_range, n_ranks, cell_part);

        if (write_output || i < n_extra_partitions)
          _write_output(mesh->n_g_cells,
                        mb->cell_bi.gnum_range,
                        n_ranks,
                        cell_part);
      }

      if (part_comm != cs_glob_mpi_comm && part_comm != MPI_COMM_NULL)
        MPI_Comm_free(&part_comm);
    }

#endif

    if (n_part_ranks == 1) {

      for (i = 0; i < n_extra_partitions + 1; i++) {

        int  n_ranks = cs_glob_n_ranks;

        if (i < n_extra_partitions) {
          n_ranks = _part_extra_partitions_list[i];
          if (n_ranks == cs_glob_n_ranks) {
            write_output = true;
            continue;
          }
        }

        if (n_ranks < 2)
          continue;

        BFT_REALLOC(cell_part, n_cells, int);

        if (cs_glob_rank_id < 0 || (cs_glob_rank_id % _part_rank_step[stage] == 0))
          _part_scotch(n_cells,
                       n_ranks,
                       cell_idx,
                       cell_neighbors,
                       cell_part);

        _distribute_output(mb,
                           _part_rank_step[stage],
                           cell_range,
                           &cell_part);

        _cell_part_histogram(mb->cell_bi.gnum_range, n_ranks, cell_part);

        if (write_output || i < n_extra_partitions)
          _write_output(mesh->n_g_cells,
                        mb->cell_bi.gnum_range,
                        n_ranks,
                        cell_part);
      }
    }

    BFT_FREE(cell_idx);
    BFT_FREE(cell_neighbors);
  }

#endif /* defined(HAVE_SCOTCH) || defined(HAVE_PTSCOTCH) */

  if (   _algorithm >= CS_PARTITION_SFC_MORTON_BOX
      && _algorithm <= CS_PARTITION_SFC_HILBERT_CUBE) {

    int i;
    fvm_io_num_sfc_t sfc_type = _algorithm - CS_PARTITION_SFC_MORTON_BOX;

    BFT_MALLOC(cell_part, n_cells, int);

    for (i = 0; i < n_extra_partitions + 1; i++) {

      int  n_ranks = cs_glob_n_ranks;

      if (i < n_extra_partitions) {
        n_ranks = _part_extra_partitions_list[i];
        if (n_ranks == cs_glob_n_ranks) {
          write_output = true;
          continue;
        }
      }

      if (n_ranks < 2)
        continue;

#if defined(HAVE_MPI)
      _cell_rank_by_sfc(mesh->n_g_cells,
                        n_ranks,
                        mb,
                        sfc_type,
                        cell_part,
                        cs_glob_mpi_comm);
#else
      _cell_rank_by_sfc(mesh->n_g_cells, n_ranks, mb, sfc_type, cell_part);
#endif

      _cell_part_histogram(mb->cell_bi.gnum_range, n_ranks, cell_part);

      if (write_output || i < n_extra_partitions)
        _write_output(mesh->n_g_cells,
                      mb->cell_bi.gnum_range,
                      n_ranks,
                      cell_part);
    }

  }

  /* Naive partitioner */

  else if (_algorithm == CS_PARTITION_BLOCK) {

    BFT_MALLOC(cell_part, n_cells, int);

    _block_partititioning(mesh, mb, cell_part);

  }

  /* Reset extra partitions list if used */

  if (n_extra_partitions > 0) {
    BFT_FREE(_part_extra_partitions_list);
    _part_n_extra_partitions = 0;
  }

  /* Copy to mesh builder */

  mb->have_cell_rank = true;
  mb->cell_rank = cell_part;

  n_cells = mb->cell_bi.gnum_range[1] - mb->cell_bi.gnum_range[0];

  /* End of this section for log file */

  t1 = cs_timer_time();

  dt = cs_timer_diff(&t0, &t1);

  bft_printf(_("\n"
               " Partitioning finished (%.3g s)\n"),
             (double)(dt.wall_nsec)/1.e9);
  bft_printf_flush();

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("  wall clock time:            %.3g s\n\n"),
                (double)(dt.wall_nsec)/1.e9);

  cs_log_separator(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
