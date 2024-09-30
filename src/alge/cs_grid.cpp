/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <chrono>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_tuning.h"
#include "cs_matrix_util.h"
#include "cs_order.h"
#include "cs_prototypes.h"
#include "cs_sles.h"
#include "cs_sort.h"

#include "fvm_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_grid.h"

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define EPZERO  1.E-12

#if !defined(HUGE_VAL)
#define HUGE_VAL  1.E+12
#endif

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/* Structure associated with grid */
/*--------------------------------*/

struct _cs_grid_t {

  int                 level;        /* Level in multigrid hierarchy */

  cs_alloc_mode_t     alloc_mode;   /* Preferred allocation mode */

  bool                conv_diff;    /* true if convection/diffusion case,
                                       false otherwise */
  bool                symmetric;    /* Symmetric matrix coefficients
                                       indicator */
  bool                use_faces;    /* True if face information is present */

  cs_lnum_t           db_size;      /* Block sizes for diagonal */
  cs_lnum_t           eb_size;      /* Block sizes for extra diagonal */

  cs_gnum_t           n_g_rows;     /* Global number of rows */

  cs_lnum_t           n_rows;       /* Local number of rows */
  cs_lnum_t           n_cols_ext;   /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */

  cs_lnum_t           n_elts_r[2];  /* Size of array used for restriction
                                       operations ({n_rows, n_cols_ext} when
                                       no grid merging has taken place) */

  /* Grid hierarchy information */

  const struct _cs_grid_t  *parent; /* Pointer to parent (finer) grid */

  /* Connectivity information */

  cs_lnum_t           n_faces;      /* Local number of faces */
  const cs_lnum_2_t  *face_cell;    /* Face -> cells connectivity (1 to n) */
  cs_lnum_2_t        *_face_cell;   /* Face -> cells connectivity
                                       (private array) */

  /* Restriction from parent to current level */

  cs_lnum_t          *coarse_row;   /* Fine -> coarse row connectivity;
                                       size: parent n_cols_ext
                                       -2: undetermined (initial value)
                                       -1: penalized
                                        {0, n-1}: standard mapping */
  cs_lnum_t          *coarse_face;  /* Fine -> coarse face connectivity
                                       (1 to n, signed:
                                       = 0 fine face inside coarse cell
                                       > 0 orientation same as parent
                                       < 0 orientation opposite as parent);
                                       size: parent n_faces */

  /* Face to cell date (if owner); uses same index as MSR matrix */

  cs_lnum_t    *cell_face;           /* Cell to faces adjacency */
  short int    *cell_face_sgn;       /* Cell to faces orientation */

  /* Geometric data */

  cs_real_t         relaxation;     /* P0/P1 relaxation parameter */
  const cs_real_t  *cell_cen;       /* Cell center (shared) */
  cs_real_t        *_cell_cen;      /* Cell center (private) */

  const cs_real_t  *cell_vol;       /* Cell volume (shared) */
  cs_real_t        *_cell_vol;      /* Cell volume (private) */

  const cs_real_t  *face_normal;    /* Surface normal of internal faces.
                                       (shared; L2 norm = face area) */
  cs_real_t        *_face_normal;   /* Surface normal of internal faces.
                                       (private; L2 norm = face area) */

  /* Parallel / periodic halo */

  const cs_halo_t  *halo;           /* Halo for this connectivity (shared) */
  cs_halo_t        *_halo;          /* Halo for this connectivity (private) */

  /* Matrix-related data */

  const cs_real_t  *da;             /* Diagonal (shared) */
  cs_real_t        *_da;            /* Diagonal (private) */

  const cs_real_t  *xa;             /* Extra-diagonal (shared) */
  cs_real_t        *_xa;            /* Extra-diagonal (private) */
  cs_real_t        *xa_conv;        /* Extra-diagonal (except level 0) */
  cs_real_t        *xa_diff;        /* Extra-diagonal (except level 0) */

  const cs_real_t  *xa0;            /* Symmetrized extra-diagonal (shared) */
  cs_real_t        *_xa0;           /* Symmetrized extra-diagonal (private) */
  cs_real_t        *xa0_diff;       /* Symmetrized extra-diagonal */

  cs_real_t        *xa0ij;

  cs_matrix_structure_t   *matrix_struct;  /* Associated matrix structure */
  const cs_matrix_t       *matrix;         /* Associated matrix (shared) */
  cs_matrix_t             *_matrix;        /* Associated matrix (private) */

#if defined(HAVE_MPI)

  /* Additional fields to allow merging grids */

  int               merge_sub_root;    /* sub-root when merging */
  int               merge_sub_rank;    /* sub-rank when merging
                                          (0 <= sub-rank < merge_sub_size) */
  int               merge_sub_size;    /* current number of merged ranks
                                          for this subset */
  int               merge_stride;      /* total number of ranks over which
                                          merging occurred at previous levels */
  int               next_merge_stride; /* total number of ranks over which
                                          merging occurred at current level */

  cs_lnum_t        *merge_cell_idx;    /* start cell_id for each sub-rank
                                          when merge_sub_rank = 0
                                          (size: merge_size + 1) */

  int               n_ranks;           /* Number of active ranks */
  MPI_Comm          comm;              /* Associated communicator */

#endif
};

/* Structure associated with traversal and update of graph */
/*---------------------------------------------------------*/

typedef struct _cs_graph_m_ptr_t {

  cs_lnum_t           m_min;        /* Minimum local number of columns */
  cs_lnum_t           m_max;        /* Maximum local number of columns */

  cs_lnum_t          *m_head;       /* pointer to head of rows list for each
                                       value of m */

  cs_lnum_t          *next;         /* pointer to next row */
  cs_lnum_t          *prev;         /* pointer to previous row */

} cs_graph_m_ptr_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

cs_real_t _penalization_threshold = 1e4;

 /* Threshold under which diagonal dominant rows are ignored in the
    aggregation scheme (assuming that dominance allows strong enough
    convergence for those rows by itself). Y Notay recommends using
    a value of 5. */

cs_real_t _dd_threshold_pw = 5;

/* Names for coarsening options */

const char *cs_grid_coarsening_type_name[]
  = {N_("default"),
     N_("SPD, diag/extra-diag ratio based"),
     N_("SPD, max extra-diag ratio based"),
     N_("SPD, (multiple) pairwise aggregation"),
     N_("convection + diffusion")};

/* Select tuning options */

static int _grid_tune_max_level = 0;
static int *_grid_tune_max_fill_level = NULL;
static cs_matrix_variant_t **_grid_tune_variant = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Binary search for a given local id in a given array of
 *        ordered local ids, when the id might not be present
 *
 * We assume the id is present in the array.
 *
 * \param[in]  l_id_array size  array_size
 * \param[in]  l_id             local id to search for
 * \param[in]  l_id_array       ordered unique local ids array
 *
 * \return  index of l_id in l_id_array, or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline cs_lnum_t
_l_id_binary_search(cs_lnum_t        l_id_array_size,
                    cs_lnum_t        l_id,
                    const cs_lnum_t  l_id_array[])
{
  if (l_id_array_size < 1)
    return -1;

  cs_lnum_t start_id = 0;
  cs_lnum_t end_id = l_id_array_size - 1;
  cs_lnum_t mid_id = (end_id -start_id) / 2;
  while (start_id < end_id) {
    if (l_id_array[mid_id] < l_id)
      start_id = mid_id + 1;
    else if (l_id_array[mid_id] > l_id)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }
  if (l_id_array[mid_id] != l_id)
    mid_id = -1;

  return mid_id;
}

/*----------------------------------------------------------------------------
 * Reduce block values to a single equivalent value for computation
 * of aggregation criteria.
 *
 * parameters:
 *   n_blocs   <-- number of blocks
 *   b_size    <-- block dimensions
 *   b_vals    <-- block values
 *   b_vals    --> scalar values
 *----------------------------------------------------------------------------*/

static void
_reduce_block(cs_lnum_t         n_blocks,
              const cs_lnum_t   b_size,
              const cs_real_t   b_vals[],
              cs_real_t         s_vals[])
{
  const cs_real_t b_div = 1.0 / b_size;
  const cs_lnum_t b_stride = b_size*b_size;

# pragma omp parallel for if(n_blocks > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_blocks; i++) {
    cs_real_t s = 0;
    for (cs_lnum_t j = 0; j < b_size; j++)
      s += b_vals[i*b_stride + b_size*j + j];
    s_vals[i] = s * b_div;
  }
}

/*----------------------------------------------------------------------------
 * Allocate empty grid structure
 *
 * - Periodic faces are not handled yet
 *
 * returns:
 *   empty grid structure
 *----------------------------------------------------------------------------*/

static cs_grid_t *
_create_grid(void)
{
  cs_grid_t *g;

  BFT_MALLOC(g, 1, cs_grid_t);

  g->conv_diff = false;
  g->symmetric = false;

  g->db_size = 1;
  g->eb_size = 1;

  g->level = 0;

  g->n_g_rows = 0;
  g->n_rows = 0;
  g->n_cols_ext = 0;

  g->n_faces = 0;

  g->n_elts_r[0] = 0;
  g->n_elts_r[1] = 0;

  g->parent = NULL;
  g->conv_diff = false;

  g->relaxation = 0;

  g->face_cell = NULL;
  g->_face_cell = NULL;

  g->coarse_row = NULL;
  g->coarse_face = NULL;

  g->cell_face = nullptr;
  g->cell_face_sgn = nullptr;

  g->cell_cen = NULL;
  g->_cell_cen = NULL;
  g->cell_vol = NULL;
  g->_cell_vol = NULL;
  g->face_normal = NULL;
  g->_face_normal = NULL;

  g->halo = NULL;
  g->_halo = NULL;

  g->da = NULL;
  g->_da = NULL;
  g->xa = NULL;
  g->_xa = NULL;
  g->xa_conv = NULL;
  g->xa_diff = NULL;
  g->xa0 = NULL;
  g->_xa0 = NULL;
  g->xa0_diff = NULL;

  g->xa0ij = NULL;

  g->matrix_struct = NULL;
  g->matrix = NULL;
  g->_matrix = NULL;

#if defined(HAVE_MPI)

  g->merge_sub_root = 0;
  g->merge_sub_rank = 0;
  g->merge_sub_size = 1;
  g->merge_stride = 0;
  g->next_merge_stride = 1;

  g->merge_cell_idx = NULL;

  g->n_ranks = cs_glob_n_ranks;
  g->comm = cs_glob_mpi_comm;

#endif
  return g;
}

/*----------------------------------------------------------------------------
 * Initialize coarse grid from fine grid
 *
 * This creates au quasi-empty grid structure, with only symmetry and
 * level information initialized, and coarsening array allocated and
 * zeroed.
 *
 * After this function is called, the coarsening array must be determined,
 * then _coarsen must be called to complete the coarse grid.
 *
 * parameters:
 *   f          <-- Pointer to fine (parent) grid structure
 *   alloc_mode <-- allocation mode
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

static cs_grid_t *
_coarse_init(const cs_grid_t  *f,
             cs_alloc_mode_t   alloc_mode)
{
  cs_grid_t *c = _create_grid();

  c->parent = f;
  c->alloc_mode = alloc_mode;

  c->level = f->level + 1;
  c->symmetric = f->symmetric;
  c->conv_diff = f->conv_diff;
  c->use_faces = f->use_faces;
  c->db_size = f->db_size;
  c->eb_size = f->eb_size;

  cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(f->matrix);
  CS_MALLOC_HD(c->coarse_row, f->n_cols_ext, cs_lnum_t, amode);

# pragma omp parallel for if(f->n_cols_ext > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < f->n_cols_ext; ii++)
    c->coarse_row[ii] = -2;

#if defined(HAVE_MPI)
  c->merge_stride = f->merge_stride;
  c->next_merge_stride = f->next_merge_stride;
  c->n_ranks = f->n_ranks;
  c->comm = f->comm;
#endif

  return c;
}

/*----------------------------------------------------------------------------
 * Log aggregation information using the face to cells adjacency.
 *
 * parameters:
 *   f         <-- Fine grid structure
 *   c         <-- Coarse grid structure
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_aggregation_stats_log(const cs_grid_t  *f,
                       const cs_grid_t  *c,
                       int               verbosity)
{
  const cs_lnum_t f_n_rows = f->n_rows;
  const cs_lnum_t c_n_rows = c->n_rows;

  cs_gnum_t c_n_rows_g = c->n_rows;

  const cs_lnum_t *c_coarse_row = (const cs_lnum_t *)(c->coarse_row);

  cs_lnum_t *c_aggr_count;

  BFT_MALLOC(c_aggr_count, c_n_rows, cs_lnum_t);
  for (cs_lnum_t i = 0; i < c_n_rows; i++)
    c_aggr_count[i] = 0;

  for (cs_lnum_t i = 0; i < f_n_rows; i++) {
    cs_lnum_t j = c_coarse_row[i];
    if (j > -1)
      c_aggr_count[j] += 1;
  }

  cs_lnum_t aggr_min = f->n_rows, aggr_max = 0;
  cs_gnum_t aggr_tot = 0;

  for (cs_lnum_t i = 0; i < c_n_rows; i++) {
    aggr_min = CS_MIN(aggr_min, c_aggr_count[i]);
    aggr_max = CS_MAX(aggr_max, c_aggr_count[i]);
    aggr_tot += c_aggr_count[i];
  }

#if defined(HAVE_MPI)
  if (f->comm != MPI_COMM_NULL) {
    cs_lnum_t _aggr_min = aggr_min, _aggr_max = aggr_max;
    cs_gnum_t _tot[2] = {aggr_tot, c_n_rows_g}, tot[2] = {0, 0};
    MPI_Allreduce(&_aggr_min, &aggr_min, 1, CS_MPI_LNUM, MPI_MIN, f->comm);
    MPI_Allreduce(&_aggr_max, &aggr_max, 1, CS_MPI_LNUM, MPI_MAX, f->comm);
    MPI_Allreduce(&_tot, &tot, 2, CS_MPI_GNUM, MPI_SUM, f->comm);
    aggr_tot = tot[0];
    c_n_rows_g = tot[1];
  }
#endif

  bft_printf("       aggregation min = %ld; max = %ld; mean = %8.2f\n",
             (long)aggr_min, (long)aggr_max,
             (double)aggr_tot/(double)c_n_rows_g);

  cs_lnum_t aggr_count = aggr_max - aggr_min + 1;

  if (verbosity > 4 && c_n_rows_g > 0 && aggr_count > 0) {

    /* Histogram */

    bft_printf("       histogram\n");

    cs_lnum_t *histogram;
    BFT_MALLOC(histogram, aggr_count, cs_lnum_t);
    for (cs_lnum_t i = 0; i < aggr_count; i++)
      histogram[i] = 0;
    for (cs_lnum_t ic = 0; ic < c_n_rows; ic++) {
      for (cs_lnum_t i = 0; i < aggr_count; i++) {
        if (c_aggr_count[ic] == aggr_min + i)
          histogram[i] += 1;
      }
    }
#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
    if (f->comm != MPI_COMM_NULL)
      MPI_Allreduce(MPI_IN_PLACE, histogram, aggr_count, CS_MPI_LNUM, MPI_SUM,
                    f->comm);
#endif
    for (cs_lnum_t i = 0; i < aggr_count; i++) {
      double epsp = 100. * histogram[i] / c_n_rows_g;
      bft_printf("         aggregation %ld = %8.2f %%\n",
                 (long)(aggr_min + i), epsp);
    }
    BFT_FREE(histogram);
  }

  BFT_FREE(c_aggr_count);
}

/*----------------------------------------------------------------------------
 * Build fine -> coarse face connectivity and coarse grid face->cell
 * connectivity (also determining the number of coarse faces) from
 * fine grid and cell coarsening array.
 *
 * Also, it is the caller's responsibility to free arrays coarse_face[]
 * and coarse_face_cell_id[] when they are no longer used.
 *
 * parameters:
 *   fine             <-- Pointer to fine (parent) grid structure
 *   coarse_row       <-- Fine -> coarse row id
 *                        size: fine->n_cols_ext
 *   n_coarse_rows    <-- Number of coarse cells
 *   n_coarse_faces   --> Number of coarse faces
 *   coarse_face      --> Fine -> coarse face connectivity (1 to n, signed:
 *                        = 0 fine face inside coarse cell
 *                        > 0 orientation same as parent
 *                        < 0 orientation opposite as parent);
 *                        size: fine->n_faces
 *   coarse_face_cell --> coarse face -> cell connectivity
 *                        size: n_coarse_faces
 *----------------------------------------------------------------------------*/

static void
_coarsen_faces(const cs_grid_t    *fine,
               const cs_lnum_t    *restrict coarse_row,
               cs_lnum_t           n_coarse_rows,
               cs_lnum_t          *n_coarse_faces,
               cs_lnum_t         **coarse_face,
               cs_lnum_2_t       **coarse_face_cell)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_lnum_t  ii, jj, face_id, connect_size;

  cs_lnum_t  *restrict c_cell_cell_cnt = NULL;
  cs_lnum_t  *restrict c_cell_cell_idx = NULL;
  cs_lnum_t  *restrict c_cell_cell_id = NULL;
  cs_lnum_t  *restrict c_cell_cell_face = NULL;

  cs_lnum_t    *restrict _coarse_face = NULL;
  cs_lnum_2_t  *restrict _c_face_cell = NULL;

  cs_lnum_t   c_n_faces = 0;

  const cs_lnum_t c_n_cells = n_coarse_rows;
  const cs_lnum_t f_n_faces = fine->n_faces;
  const cs_lnum_2_t *restrict f_face_cell = fine->face_cell;

  /* Pre-allocate return values
     (coarse face->cell connectivity is over-allocated) */

  BFT_MALLOC(_coarse_face, f_n_faces, cs_lnum_t);
  BFT_MALLOC(_c_face_cell, f_n_faces, cs_lnum_2_t);

# pragma omp parallel for if(f_n_faces > CS_THR_MIN)
  for (face_id = 0; face_id < f_n_faces; face_id++)
    _coarse_face[face_id] = 0;

  /* Allocate index */

  BFT_MALLOC(c_cell_cell_idx, c_n_cells + 1, cs_lnum_t);

# pragma omp parallel for if(c_n_cells > CS_THR_MIN)
  for (ii = 0; ii <= c_n_cells; ii++)
    c_cell_cell_idx[ii] = 0;

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    ii = coarse_row[f_face_cell[face_id][0]];
    jj = coarse_row[f_face_cell[face_id][1]];

    if (ii < jj && ii > -1)
      c_cell_cell_idx[ii+1] += 1;
    else if (ii > jj && jj > -1)
      c_cell_cell_idx[jj+1] += 1;

  }

  /* Transform to index */

  for (ii = 0; ii < c_n_cells; ii++)
    c_cell_cell_idx[ii+1] += c_cell_cell_idx[ii];

  BFT_MALLOC(c_cell_cell_id,
             c_cell_cell_idx[c_n_cells],
             cs_lnum_t);

  BFT_MALLOC(c_cell_cell_face,
             c_cell_cell_idx[c_n_cells],
             cs_lnum_t);

  connect_size = c_cell_cell_idx[c_n_cells];

# pragma omp parallel for if(connect_size > CS_THR_MIN)
  for (ii = 0; ii < connect_size; ii++)
    c_cell_cell_face[ii] = 0;

  /* Use a counter for array population, as array will usually
     not be fully populated */

  BFT_MALLOC(c_cell_cell_cnt, c_n_cells, cs_lnum_t);

# pragma omp parallel for if(c_n_cells > CS_THR_MIN)
  for (ii = 0; ii < c_n_cells; ii++)
    c_cell_cell_cnt[ii] = 0;

  /* Build connectivity */

  c_n_faces = 0;

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    cs_lnum_t kk, start_id, end_id;
    cs_lnum_t sign = 1;

    ii = coarse_row[f_face_cell[face_id][0]];
    jj = coarse_row[f_face_cell[face_id][1]];

    if (ii != jj && ii > -1 && jj > -1) {

      if (ii > jj) {
        sign = -1;
        kk = ii;
        ii = jj;
        jj = kk;
      }

      start_id = c_cell_cell_idx[ii];
      end_id   = c_cell_cell_idx[ii] + c_cell_cell_cnt[ii];

      for (kk = start_id; kk < end_id; kk++) {
        if (c_cell_cell_id[kk] == jj)
          break;
      }
      if (kk == end_id) {
        c_cell_cell_id[kk] = jj;
        c_cell_cell_face[kk] = c_n_faces + 1;
        _c_face_cell[c_n_faces][0] = ii;
        _c_face_cell[c_n_faces][1] = jj;
        c_n_faces++;
        c_cell_cell_cnt[ii] += 1;
      }
      _coarse_face[face_id] = sign*c_cell_cell_face[kk];

    }

  }

  /* Free temporary arrays */

  BFT_FREE(c_cell_cell_cnt);
  BFT_FREE(c_cell_cell_face);
  BFT_FREE(c_cell_cell_id);
  BFT_FREE(c_cell_cell_idx);

  /* Set return values */

  BFT_REALLOC(_c_face_cell, c_n_faces, cs_lnum_2_t);

  *n_coarse_faces = c_n_faces;
  *coarse_face = _coarse_face;
  *coarse_face_cell = _c_face_cell;

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           fine->level, fine->level - 1);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Send prepared coarsening ids (based on halo send lists) and receive
 * similar values to a cell coarsening array halo.
 *
 * parameters:
 *   halo          <->  pointer to halo structure
 *   coarse_send   <->  values defined on halo send lists
 *   coarse_row    <->  pointer to local row coarsening array
 *                      whose halo values are to be received
 *----------------------------------------------------------------------------*/

static void
_exchange_halo_coarsening(const cs_halo_t  *halo,
                          cs_lnum_t         coarse_send[],
                          cs_lnum_t         coarse_row[])
{
  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int  request_count = 0;
    const int  local_rank = cs_glob_rank_id;

    MPI_Request _request[128];
    MPI_Request *request = _request;
    MPI_Status _status[128];
    MPI_Status *status = _status;

    /* Allocate if necessary */

    if (halo->n_c_domains*2 > 128) {
      BFT_MALLOC(request, halo->n_c_domains*2, MPI_Request);
      BFT_MALLOC(status, halo->n_c_domains*2, MPI_Status);
    }

    /* Receive data from distant ranks */

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      cs_lnum_t start = halo->index[2*rank_id];
      cs_lnum_t length = halo->index[2*rank_id + 2] - halo->index[2*rank_id];

      if (halo->c_domain_rank[rank_id] != local_rank) {

        MPI_Irecv(coarse_row + halo->n_local_elts + start,
                  length,
                  CS_MPI_LNUM,
                  halo->c_domain_rank[rank_id],
                  halo->c_domain_rank[rank_id],
                  cs_glob_mpi_comm,
                  &(request[request_count++]));

      }
      else
        local_rank_id = rank_id;

    }

    /* We wait for posting all receives (often recommended) */

    // MPI_Barrier(comm);

    /* Send data to distant ranks */

    for (int rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank_id] != local_rank) {

        cs_lnum_t start = halo->send_index[2*rank_id];
        cs_lnum_t length =   halo->send_index[2*rank_id + 2]
                           - halo->send_index[2*rank_id];

        MPI_Isend(coarse_send + start,
                  length,
                  CS_MPI_LNUM,
                  halo->c_domain_rank[rank_id],
                  local_rank,
                  cs_glob_mpi_comm,
                  &(request[request_count++]));

      }

    }

    /* Wait for all exchanges */

    MPI_Waitall(request_count, request, status);

    if (request != _request) {
      BFT_FREE(request);
      BFT_FREE(status);
    }
  }

#endif /* defined(HAVE_MPI) */

  /* Copy local values in case of periodicity */

  if (halo->n_transforms > 0) {

    if (local_rank_id > -1) {

      cs_lnum_t *_coarse_row
        = coarse_row + halo->n_local_elts + halo->index[2*local_rank_id];

      cs_lnum_t start = halo->send_index[2*local_rank_id];
      cs_lnum_t length =   halo->send_index[2*local_rank_id + 2]
                         - halo->send_index[2*local_rank_id];

#     pragma omp parallel for if(length > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < length; i++)
        _coarse_row[i] = coarse_send[start + i];

    }

  }
}

/*----------------------------------------------------------------------------
 * Build a halo for a coarse grid from a fine grid and its halo.
 *
 * The coarse grid must previously have been initialized with _coarse_init()
 * and its coarsening indicator determined for the local cells.
 *
 * parameters:
 *   f <-- Pointer to fine (parent) grid structure
 *   c <-> Pointer to coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_coarsen_halo(const cs_grid_t   *f,
              cs_grid_t         *c)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  int domain_id, tr_id, section_id;
  cs_lnum_t ii, jj;
  cs_lnum_t start_id, end_id, sub_count;

  cs_lnum_t *start_end_id = NULL;
  cs_lnum_t *sub_num = NULL;
  cs_lnum_t  *coarse_send = NULL;

  cs_lnum_t *restrict coarse_row = c->coarse_row;

  cs_halo_t *c_halo = NULL;
  const cs_halo_t *f_halo = f->halo;

  const cs_lnum_t c_n_rows = c->n_rows;

  const int stride = f_halo->n_c_domains*4;
  const int n_sections = f_halo->n_transforms + 1;
  const int n_f_send = f_halo->n_send_elts[1]; /* Size of full list */

  c_halo = cs_halo_create_from_ref(f_halo);

  c->halo = c_halo;
  c->_halo = c_halo;

  /* Initialize coarse halo counters */

  c_halo->n_local_elts = c_n_rows;
  c_halo->n_send_elts[0] = 0;
  c_halo->n_send_elts[1] = 0;

# pragma omp parallel for if(f->n_cols_ext > CS_THR_MIN)
  for (ii = f->n_rows; ii < f->n_cols_ext; ii++)
    coarse_row[ii] = -1;

  /* Allocate and initialize counters */

  BFT_MALLOC(start_end_id, n_sections*2, cs_lnum_t);
  BFT_MALLOC(sub_num, c_n_rows + 1, cs_lnum_t);
  BFT_MALLOC(coarse_send, n_f_send, cs_lnum_t);

  /* Sub_num values shifted by 1, so 0 can be used for handling
     of removed (penalized) rows) */

  sub_num[0] = -2;
# pragma omp parallel for if(c_n_rows > CS_THR_MIN)
  for (ii = 1; ii <= c_n_rows; ii++)
    sub_num[ii] = -1;

# pragma omp parallel for if(n_f_send > CS_THR_MIN)
  for (ii = 0; ii < n_f_send; ii++)
    coarse_send[ii] = -1;

  /* Counting and marking pass */
  /*---------------------------*/

  /* Proceed halo section by halo section */

  for (domain_id = 0; domain_id < f_halo->n_c_domains; domain_id++) {

    /* Handle purely parallel cells */

    start_end_id[0] = f_halo->send_index[domain_id*2];

    if (f_halo->n_transforms == 0)
      start_end_id[1] = f_halo->send_index[(domain_id+1)*2];
    else
      start_end_id[1] = f_halo->send_perio_lst[4*domain_id];

    /* Now handle periodic cells transformation by transformation */

    for (tr_id = 0; tr_id < f_halo->n_transforms; tr_id++) {
      start_end_id[tr_id*2 + 2]
        = f_halo->send_perio_lst[stride*tr_id + 4*domain_id];
      start_end_id[tr_id*2 + 3]
        =   start_end_id[tr_id*2 + 2]
          + f_halo->send_perio_lst[stride*tr_id + 4*domain_id + 1];
    }

    /* We can now loop on sections independently of the transform */

    for (section_id = 0; section_id < n_sections; section_id++) {

      start_id = start_end_id[section_id*2];
      end_id = start_end_id[section_id*2 + 1];

      sub_count = 0;

      if (section_id > 0)
        c_halo->send_perio_lst[stride*(section_id-1) + 4*domain_id]
          = c_halo->n_send_elts[0];

      /* Build halo renumbering */

      for (ii = start_id; ii < end_id; ii++) {
        jj = coarse_row[f_halo->send_list[ii]] + 1;
        if (sub_num[jj] == -1) {
          sub_num[jj] = sub_count;
          sub_count += 1;
        }
        coarse_send[ii] = sub_num[jj];
      }

      /* Update send_list indexes and counts */

      c_halo->n_send_elts[0] += sub_count;

      if (section_id > 0) {
        c_halo->send_perio_lst[stride*(section_id-1) + 4*domain_id + 1]
          = sub_count;
        c_halo->send_perio_lst[stride*(section_id-1) + 4*domain_id + 2]
          = c_halo->n_send_elts[0];
        c_halo->send_perio_lst[stride*(section_id-1) + 4*domain_id + 3]
          = 0;
      }

      /* Reset sub_num for next section or domain */

      for (ii = start_id; ii < end_id; ii++)
        sub_num[coarse_row[f_halo->send_list[ii]] + 1] = -1;
      sub_num[0] = -2;

    }

    /* Update send index (initialized at halo creation) */

    c_halo->send_index[domain_id*2 + 1] = c_halo->n_send_elts[0];
    c_halo->send_index[domain_id*2 + 2] = c_halo->n_send_elts[0];

  }

  /* Exchange and update coarsening list halo */
  /*------------------------------------------*/

  _exchange_halo_coarsening(f_halo, coarse_send, coarse_row);

  BFT_FREE(coarse_send);

  /* Proceed halo section by halo section */

  c_halo->n_elts[0] = 0;
  c_halo->n_elts[1] = 0;
  c_halo->index[0] = 0;

  for (domain_id = 0; domain_id < f_halo->n_c_domains; domain_id++) {

    /* Handle purely parallel cells */

    start_end_id[0] = f_halo->index[domain_id*2];

    if (f_halo->n_transforms == 0)
      start_end_id[1] = f_halo->index[(domain_id+1)*2];
    else
      start_end_id[1] = f_halo->perio_lst[4*domain_id];

    /* Now handle periodic cells transformation by transformation */

    for (tr_id = 0; tr_id < f_halo->n_transforms; tr_id++) {
      start_end_id[tr_id*2 + 2]
        = f_halo->perio_lst[stride*tr_id + 4*domain_id];
      start_end_id[tr_id*2 + 3]
        =   start_end_id[tr_id*2 + 2]
          + f_halo->perio_lst[stride*tr_id + 4*domain_id + 1];
    }

    /* We can now loop on sections independently of the transform */

    for (section_id = 0; section_id < n_sections; section_id++) {

      start_id = f_halo->n_local_elts + start_end_id[section_id*2];
      end_id = f_halo->n_local_elts + start_end_id[section_id*2 + 1];

      sub_count = 0;

      /* Build halo renumbering */

      for (ii = start_id, jj = -1; ii < end_id; ii++) {
        if (coarse_row[ii] > -1) {
          if (coarse_row[ii] > jj)
            jj += 1;
          coarse_row[ii] += c_n_rows + c_halo->n_elts[0];
        }
      }

      if (section_id > 0) {
        c_halo->perio_lst[stride*(section_id-1) + 4*domain_id]
          = c_halo->n_elts[0];
        c_halo->perio_lst[stride*(section_id-1) + 4*domain_id + 1]
          = jj + 1;
      }

      c_halo->n_elts[0] += jj + 1;

    }

    for (tr_id = 0; tr_id < f_halo->n_transforms; tr_id++) {
      c_halo->perio_lst[stride*tr_id + 4*domain_id + 2]
        = c_halo->n_elts[0];
      c_halo->perio_lst[stride*tr_id + 4*domain_id + 3]
        = 0;
    }

    /* Update index */

    c_halo->n_elts[1] = c_halo->n_elts[0];

    c_halo->index[domain_id*2 + 1] = c_halo->n_elts[0];
    c_halo->index[domain_id*2 + 2] = c_halo->n_elts[0];
  }

  /* Pass to build send list */
  /*-------------------------*/

  CS_MALLOC_HD(c_halo->send_list, c_halo->n_send_elts[0], cs_lnum_t,
               cs_check_device_ptr(c_halo->send_index));

  c_halo->n_send_elts[0] = 0;

  /* Proceed halo section by halo section */

  for (domain_id = 0; domain_id < f_halo->n_c_domains; domain_id++) {

    /* Handle purely parallel cells */

    start_end_id[0] = f_halo->send_index[domain_id*2];

    if (f_halo->n_transforms == 0)
      start_end_id[1] = f_halo->send_index[(domain_id+1)*2];
    else
      start_end_id[1] = f_halo->send_perio_lst[4*domain_id];

    /* Now handle periodic cells transformation by transformation */

    for (tr_id = 0; tr_id < f_halo->n_transforms; tr_id++) {
      start_end_id[tr_id*2 + 2]
        = f_halo->send_perio_lst[stride*tr_id + 4*domain_id];
      start_end_id[tr_id*2 + 3]
        =   start_end_id[tr_id*2 + 2]
          + f_halo->send_perio_lst[stride*tr_id + 4*domain_id + 1];
    }

    /* We can now loop on sections independently of the transform */

    for (section_id = 0; section_id < n_sections; section_id++) {

      start_id = start_end_id[section_id*2];
      end_id = start_end_id[section_id*2 + 1];

      sub_count = 0;

      if (section_id > 0)
        c_halo->send_perio_lst[stride*(section_id-1) + 4*domain_id]
          = c_halo->n_send_elts[0];

      /* Build halo renumbering */

      for (ii = start_id; ii < end_id; ii++) {
        jj = coarse_row[f_halo->send_list[ii]] + 1;
        if (sub_num[jj] == -1) {
          sub_num[jj] = sub_count;
          c_halo->send_list[c_halo->n_send_elts[0] + sub_count] = jj - 1;
          sub_count += 1;
        }
      }

      c_halo->n_send_elts[0] += sub_count;

      /* Reset sub_num for next section or domain */

      for (ii = start_id; ii < end_id; ii++)
        sub_num[coarse_row[f_halo->send_list[ii]] + 1] = -1;
      sub_num[0] = -2;

    }

  }

  c_halo->n_send_elts[1] = c_halo->n_send_elts[0];

  /* Prune halo */
  /*------------*/

  /* Proceed halo section by halo section */

  int domain_count = 0;

  for (domain_id = 0; domain_id < c_halo->n_c_domains; domain_id++) {

    cs_lnum_t recv_size =   c_halo->index[2*domain_id + 2]
                          - c_halo->index[2*domain_id];
    cs_lnum_t send_size =   c_halo->send_index[2*domain_id + 2]
                          - c_halo->send_index[2*domain_id];

    if (recv_size + send_size == 0) {
      c_halo->c_domain_rank[domain_id] = -1;
    }

    else {
      for (int i = 0; i < 2; i++) {
        c_halo->index[2*domain_count + i]
          = c_halo->index[2*domain_id + i];
        c_halo->send_index[2*domain_count + i]
          = c_halo->send_index[2*domain_id + i];
      }

      domain_count++;
    }

  }

  c_halo->index[2*domain_count] = c_halo->index[2*c_halo->n_c_domains];
  c_halo->send_index[2*domain_count] = c_halo->send_index[2*c_halo->n_c_domains];

  CS_REALLOC_HD(c_halo->index, domain_count*2+1, cs_lnum_t,
                cs_check_device_ptr(c_halo->index));
  CS_REALLOC_HD(c_halo->send_index, domain_count*2+1, cs_lnum_t,
                cs_check_device_ptr(c_halo->send_index));

  if (domain_count < c_halo->n_c_domains && n_sections > 0) {

    const cs_lnum_t  stride_o = 4*c_halo->n_c_domains;
    const cs_lnum_t  stride_n = 4*domain_count;

    cs_lnum_t *send_perio_lst, *perio_lst;
    BFT_MALLOC(send_perio_lst, c_halo->n_transforms*stride_n, cs_lnum_t);
    CS_MALLOC_HD(perio_lst,
                 c_halo->n_transforms*stride_n,
                 cs_lnum_t,
                 c->alloc_mode);

    for (int i = 0; i < c_halo->n_transforms; i++) {

      int j = 0;

      for (int k = 0; k < c_halo->n_c_domains; k++) {

        if (c_halo->c_domain_rank[k] > -1) {

          for (int l = 0; l < 4; l++)
            send_perio_lst[i*stride_n + 4*j + l]
              = c_halo->send_perio_lst[i*stride_o + 4*k + l];

          for (int l = 0; l < 4; l++)
            perio_lst[i*stride_n + 4*j + l]
              = c_halo->perio_lst[i*stride_o + 4*k + l];

          j += 1;

        }

      }

      assert(j == domain_count);

    }

    BFT_FREE(c_halo->send_perio_lst);
    CS_FREE(c_halo->perio_lst);
    c_halo->send_perio_lst = send_perio_lst;
    c_halo->perio_lst = perio_lst;

  }

  domain_count = 0;
  for (int i = 0; i < c_halo->n_c_domains; i++) {
    if (c_halo->c_domain_rank[i] > -1) {
      c_halo->c_domain_rank[domain_count] = c_halo->c_domain_rank[i];
      domain_count++;
    }
  }
  c_halo->n_c_domains = domain_count;
  BFT_REALLOC(c_halo->c_domain_rank, c_halo->n_c_domains, int);

  if (cs_check_device_ptr(c_halo->send_index) > CS_ALLOC_HOST) {
    cs_sync_h2d(c_halo->send_index);
    cs_sync_h2d(c_halo->send_list);
  }

  /* Free memory */

  BFT_FREE(coarse_send);
  BFT_FREE(sub_num);
  BFT_FREE(start_end_id);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, c->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Build coarse grid from fine grid
 *
 * The coarse grid must previously have been initialized with _coarse_init()
 * and its coarsening indicator determined (at least for the local cells;
 * it is extended to halo cells here if necessary).
 *
 * - Periodic faces are not handled yet
 *
 * parameters:
 *   f <-- Pointer to fine (parent) grid structure
 *   c <-> Pointer to coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_coarsen(const cs_grid_t   *f,
         cs_grid_t         *c)
{
  cs_lnum_t  c_n_rows = c->n_rows;

  const cs_lnum_t f_n_faces = f->n_faces;
  const cs_lnum_t f_n_rows = cs_matrix_get_n_rows(f->matrix);
  const cs_lnum_2_t *restrict f_face_cell = f->face_cell;

  /* Sanity check */

  if (f_face_cell != NULL) {
#   pragma omp parallel for if(f_n_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < f_n_faces; face_id++) {
      cs_lnum_t ii = f_face_cell[face_id][0];
      cs_lnum_t jj = f_face_cell[face_id][1];
      if (ii == jj)
        bft_error(__FILE__, __LINE__, 0,
                  _("Connectivity error:\n"
                    "Face %d has same cell %d on both sides."),
                  (int)(face_id+1), (int)(ii+1));
    }
  }

  /* Compute number of coarse rows */

  if (c_n_rows < 0) {
    c_n_rows = 0;
    for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
      if (c->coarse_row[ii] >= c_n_rows)
        c_n_rows = c->coarse_row[ii] + 1;
    }
    c->n_rows = c_n_rows;
  }

  c->n_g_rows = c_n_rows;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    if (f->comm != MPI_COMM_NULL) {
      cs_gnum_t _c_n_rows = c_n_rows;
      MPI_Allreduce(&_c_n_rows, &(c->n_g_rows), 1, CS_MPI_GNUM, MPI_SUM,
                    f->comm);
    }
  }
#endif

  /* Prolong mesh coarsening indicator to halo rows and build
     coarse mesh halos if necessary */

  if (f->halo != NULL) {
    _coarsen_halo(f, c);
    c->n_cols_ext = c->n_rows + c->halo->n_elts[0];
  }
  else
    c->n_cols_ext = c_n_rows;

  c->n_elts_r[0] = c->n_rows;
  c->n_elts_r[1] = c->n_cols_ext;

  /* Build face coarsening and coarse grid face -> cells connectivity */

  if (  f->face_cell != NULL
      && (   c->relaxation > 0
          || cs_matrix_get_type(f->matrix) == CS_MATRIX_NATIVE)) {
    _coarsen_faces(f,
                   c->coarse_row,
                   c->n_rows,
                   &(c->n_faces),
                   &(c->coarse_face),
                   &(c->_face_cell));

    c->face_cell = (const cs_lnum_2_t  *)(c->_face_cell);
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Merge halo info after appending arrays.
 *
 * parameters:
 *   h                 <-> Pointer to halo structure
 *   alloc_mode        <-- Memory allocation mode
 *   new_src_cell_num  <-> new list of cells the senders should provide
 *----------------------------------------------------------------------------*/

static void
_rebuild_halo_send_lists(cs_halo_t        *h,
                         cs_alloc_mode_t   alloc_mode,
                         cs_lnum_t         new_src_cell_id[])
{
  cs_lnum_t start, length;

  int rank_id, tr_id;
  int n_sections = 1 + h->n_transforms;
  int request_count = 0;
  cs_lnum_t *send_buf = NULL, *recv_buf = NULL;
  MPI_Status *status = NULL;
  MPI_Request *request = NULL;

  BFT_MALLOC(status, h->n_c_domains*2, MPI_Status);
  BFT_MALLOC(request, h->n_c_domains*2, MPI_Request);
  BFT_MALLOC(send_buf, h->n_c_domains*n_sections, cs_lnum_t);
  BFT_MALLOC(recv_buf, h->n_c_domains*n_sections, cs_lnum_t);

  /* Exchange sizes */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++)
    MPI_Irecv(recv_buf + rank_id*n_sections,
              n_sections,
              CS_MPI_LNUM,
              h->c_domain_rank[rank_id],
              h->c_domain_rank[rank_id],
              cs_glob_mpi_comm,
              &(request[request_count++]));

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    send_buf[rank_id*n_sections]
      = h->index[2*(rank_id+1)] - h->index[2*rank_id];
    for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
      send_buf[rank_id*n_sections + tr_id + 1]
        = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 1];
    }
    MPI_Isend(send_buf + rank_id*n_sections,
              n_sections,
              CS_MPI_LNUM,
              h->c_domain_rank[rank_id],
              cs_glob_rank_id,
              cs_glob_mpi_comm,
              &(request[request_count++]));
  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  /* Update sizes */

  CS_MALLOC_HD(h->send_index, h->n_c_domains*2 + 1, cs_lnum_t, CS_ALLOC_HOST);
  h->send_index[0] = 0;
  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    h->send_index[rank_id*2 + 1]
      = h->send_index[rank_id*2] + recv_buf[rank_id*n_sections];
    h->send_index[rank_id*2 + 2] = h->send_index[rank_id*2 + 1];
  }

  /* Update send_perio_lst in case of transforms */

  if (h->n_transforms > 0) {
    BFT_MALLOC(h->send_perio_lst, h->n_c_domains*h->n_transforms*4, cs_lnum_t);
    for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
      cs_lnum_t n_cur_vals = recv_buf[rank_id*n_sections];
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++)
        n_cur_vals -= recv_buf[rank_id*n_sections + 1 + tr_id];
      n_cur_vals += h->send_index[rank_id*2];
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
        cs_lnum_t n_tr_vals = recv_buf[rank_id*n_sections + 1 + tr_id];
        h->send_perio_lst[h->n_c_domains*4*tr_id + 4*rank_id] = n_cur_vals;
        h->send_perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 1] = n_tr_vals;
        n_cur_vals += n_tr_vals;
        h->send_perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 2] = n_cur_vals;
        h->send_perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 3] = 0;
      }
    }
  }

  BFT_FREE(send_buf);
  BFT_FREE(recv_buf);

  h->n_send_elts[0] = h->send_index[h->n_c_domains*2];
  h->n_send_elts[1] = h->n_send_elts[0];

  CS_MALLOC_HD(h->send_list, h->n_send_elts[0], cs_lnum_t,
               alloc_mode);

  /* Receive data from distant ranks */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    start = h->send_index[2*rank_id];
    length = h->send_index[2*rank_id + 1] - h->send_index[2*rank_id];
    MPI_Irecv(h->send_list + start,
              length,
              CS_MPI_LNUM,
              h->c_domain_rank[rank_id],
              h->c_domain_rank[rank_id],
              cs_glob_mpi_comm,
              &(request[request_count++]));
  }

  /* Send data to distant ranks */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    start = h->index[2*rank_id];
    length = h->index[2*rank_id + 1] - h->index[2*rank_id];
    MPI_Isend(new_src_cell_id + start,
              length,
              CS_MPI_LNUM,
              h->c_domain_rank[rank_id],
              cs_glob_rank_id,
              cs_glob_mpi_comm,
              &(request[request_count++]));
  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);

  BFT_FREE(request);
  BFT_FREE(status);
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
  if (h == NULL)
    return;

  h->n_c_domains = 0;
  BFT_FREE(h->c_domain_rank);

  h->n_send_elts[0] = 0;
  h->n_send_elts[1] = 0;
  h->n_elts[0] = 0;
  h->n_elts[1] = 0;

  CS_FREE_HD(h->send_list);
  CS_FREE_HD(h->send_index);
  BFT_FREE(h->send_perio_lst);
  BFT_FREE(h->index);
  BFT_FREE(h->perio_lst);
}

/*----------------------------------------------------------------------------
 * Merge halo info after appending arrays.
 *
 * parameters:
 *   h                 <-> pointer to halo structure
 *   loc_rank_id       <-- local rank id
 *   n_new_cells       <-- number of new local cells
 *   new_src_cell_id   <-> in: new halo sender cell id matching halo cell;
 *                         out: new list of cells the senders should provide
 *                              (same numbers, merged and possibly reordered)
 *   new_halo_cell_id --> new cell id for each halo cell
 *                         (< n_new_cells for cells that have become local,
 *                         >= n_new_cells for cells still in halo)
 *----------------------------------------------------------------------------*/

static void
_merge_halo_data(cs_halo_t   *h,
                 int          loc_rank_id,
                 cs_lnum_t    n_new_cells,
                 cs_lnum_t    new_src_cell_id[],
                 cs_lnum_t    new_halo_cell_id[])
{
  const int  stride = (h->n_transforms > 0) ? 3 : 2;
  const int  n_c_domains_ini = h->n_c_domains;
  const int  n_sections = h->n_transforms + 1;

  cs_lnum_t   n_elts_ini = h->n_elts[0];

  if (h->n_elts[0] + h->n_send_elts[0] < 1) {
    _empty_halo(h);
    return;
  }

  /* Order list by rank, transform, and new element number */

  cs_gnum_t  *tmp_num = NULL;
  BFT_MALLOC(tmp_num, n_elts_ini*stride, cs_gnum_t);

  for (int rank_idx = 0; rank_idx < n_c_domains_ini; rank_idx++) {
    int c_rank_id = h->c_domain_rank[rank_idx];
    if (c_rank_id == loc_rank_id) {
      h->c_domain_rank[rank_idx] = -1; /* to appear first */
      c_rank_id = -1;
    }

    for (cs_lnum_t ii = h->index[rank_idx*2];
         ii < h->index[(rank_idx+1)*2];
         ii++)
      tmp_num[ii*stride] = c_rank_id + 1;
  }

  if (stride == 2) {
    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++)
      tmp_num[ii*2 + 1] = new_src_cell_id[ii];
  }
  else { /* if (stride == 3) */

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {
      tmp_num[ii*3 + 1] = 0;
      tmp_num[ii*3 + 2] = new_src_cell_id[ii];
    }

    for (int rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
      for (int tr_id = 0; tr_id < h->n_transforms; tr_id++) {
        cs_lnum_t ii_0 = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id];
        cs_lnum_t ii_1 = ii_0 + h->perio_lst[  h->n_c_domains*4*tr_id
                                             + 4*rank_id + 1];
        for (cs_lnum_t ii = ii_0; ii < ii_1; ii++)
          tmp_num[ii*3 + 1] = tr_id + 1;
      }
    }
  }

  /* Rebuild lists and build renumbering.

     As the send and receive indexes might not contain exactly the same
     list of distant ranks (at least in case of periodicity), we must be careful
     not to neighbor information here, so we only remove duplicates in the
     current list. */

  {
    int j = 0;
    cs_sort_int_shell(0, h->n_c_domains, h->c_domain_rank);

    for (int i = 1; i < h->n_c_domains; i++) {
      if (h->c_domain_rank[i] != h->c_domain_rank[j]) {
        j++;
        h->c_domain_rank[j] = h->c_domain_rank[i];
      }
    }

    h->n_c_domains = j+1;
    BFT_REALLOC(h->c_domain_rank, h->n_c_domains, int);
  }

  cs_lnum_t *section_idx = NULL;

  cs_lnum_t *order = cs_order_gnum_s(NULL, tmp_num, stride, n_elts_ini);

  for (int i = 0; i < h->n_c_domains*2 + 1; i++)
    h->index[i] = 0;
  h->n_elts[0] = 0;
  h->n_elts[1] = 0;

  cs_lnum_t prev_src_id = 0;

  if (stride == 2) {

    int prev_rank_id = -1, rank_idx = -1;

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {

      bool is_same = true;
      cs_lnum_t cur_id = order[ii];

      int rank_id = tmp_num[cur_id*2] - 1;
      cs_lnum_t src_id = tmp_num[cur_id*2 + 1];

      if (rank_id != -1) { /* local rank */

        if (rank_id != prev_rank_id) {
          if (rank_idx > -1)
            h->index[rank_idx*2+1] = h->n_elts[0];
          rank_idx++;
          while (   h->c_domain_rank[rank_idx] != rank_id
                 && rank_idx < h->n_c_domains)
            rank_idx++;
          assert(rank_idx < h->n_c_domains);
          is_same = false;
        }
        if (src_id != prev_src_id)
          is_same = false;

        if (is_same == false) {
          new_src_cell_id[h->n_elts[0]] = src_id;
          h->n_elts[0] += 1;
        }

        new_halo_cell_id[cur_id] = n_new_cells + h->n_elts[0] - 1;

        prev_rank_id = rank_id;
        prev_src_id = src_id;

      }
      else { /* local rank */
        new_halo_cell_id[cur_id] = tmp_num[cur_id*2 + 1];
        assert(new_halo_cell_id[cur_id] < n_new_cells);
      }
    }

    if (rank_idx > -1)
      h->index[rank_idx*2+1] = h->n_elts[0];

  }
  else { /* if (stride == 3) */

    int prev_rank_id = -2, rank_idx = -1;
    cs_lnum_t section_idx_size = 0;

    int prev_section_id = -1;

    /* Index will be initialized as count */

    for (cs_lnum_t ii = 0; ii < n_elts_ini; ii++) {

      bool is_same = true;
      cs_lnum_t cur_id = order[ii];

      int rank_id = tmp_num[cur_id*3] - 1;
      int section_id = tmp_num[cur_id*3 + 1];
      cs_lnum_t src_id = tmp_num[cur_id*3 + 2];

      if (rank_id != -1 || tmp_num[cur_id*3 + 1] != 0) {

        if (rank_id != prev_rank_id) {
          if (rank_idx > -1)
            h->index[rank_idx*2+1] = h->n_elts[0];
          rank_idx++;
          while (   h->c_domain_rank[rank_idx] != rank_id
                 && rank_idx < h->n_c_domains)
            rank_idx++;
          assert(rank_idx < h->n_c_domains);
          is_same = false;
        }
        if (section_id != prev_section_id)
          is_same = false;
        if (src_id != prev_src_id)
          is_same = false;

        if (is_same == false) {
          new_src_cell_id[h->n_elts[0]] = src_id;
          h->n_elts[0] += 1;
          while (rank_idx*n_sections + section_id + 1 >= section_idx_size) {
            cs_lnum_t section_idx_size_prv = section_idx_size;
            section_idx_size
              = (section_idx_size_prv > 0) ? section_idx_size_prv*2 : 16;
            BFT_REALLOC(section_idx, section_idx_size, cs_lnum_t);
            for (cs_lnum_t sid = section_idx_size_prv;
                 sid < section_idx_size;
                 sid++)
              section_idx[sid] = 0;
          };
          section_idx[rank_idx*n_sections + section_id + 1] += 1;
        }

        new_halo_cell_id[cur_id] = n_new_cells + h->n_elts[0] - 1;

        prev_rank_id = rank_id;
        prev_section_id = section_id;
        prev_src_id = src_id;

      }
      else { /* if (rank_id == -1 && tmp_num[cur_id*3 + 1] == 0) */
        new_halo_cell_id[cur_id] = tmp_num[cur_id*3 + 2];
        assert(new_halo_cell_id[cur_id] < n_new_cells);
      }

    }

    if (rank_idx > -1)
      h->index[rank_idx*2+1] = h->n_elts[0];

    /* Transform count to index */

    section_idx_size = n_sections * h->n_c_domains + 1;
    for (cs_lnum_t sid = 1; sid < section_idx_size; sid++)
      section_idx[sid] += section_idx[sid - 1];
  }

  BFT_FREE(order);
  BFT_FREE(tmp_num);

  /* Restore local rank in connected domains list */

  for (int i = 0; i < h->n_c_domains; i++) {
    if (h->c_domain_rank[i] == -1)
      (h->c_domain_rank[i] = loc_rank_id);
  }

  /* Update transforms list */

  if (h->n_transforms > 0) {

    for (int tr_id = 0; tr_id < h->n_transforms; tr_id++) {
      for (int rank_idx = 0; rank_idx < h->n_c_domains; rank_idx++) {
        int section_id = rank_idx*n_sections + tr_id + 1;
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx]
          = section_idx[section_id];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx + 1]
          = section_idx[section_id + 1] - section_idx[section_id];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx + 2]
          = section_idx[section_id + 1];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_idx + 3]
          = 0;
      }
    }

    BFT_FREE(section_idx);

  }

  /* Free and resize memory */

  BFT_REALLOC(h->index, h->n_c_domains*2+1, cs_lnum_t);
  if (h->n_transforms > 0)
    BFT_REALLOC(h->perio_lst,
                h->n_c_domains * h->n_transforms * 4,
                cs_lnum_t);

  h->n_elts[1] = h->n_elts[0];
  h->index[h->n_c_domains*2] = h->n_elts[0];
  for (int i = 1; i < h->n_c_domains; i++) {
    if (h->index[i*2+1] == 0)
      h->index[i*2+1] = h->index[i*2-1];
    h->index[i*2] = h->index[i*2-1];
  }
}

/*----------------------------------------------------------------------------
 * Append halo info for grid grouping and merging.
 *
 * parameters:
 *   g           <-> pointer to grid structure being merged
 *   new_cell_id <-> new cell ids for local cells
 *                   in: defined for local cells
 *                   out: updated also for halo cells
 *----------------------------------------------------------------------------*/

static void
_append_halos(cs_grid_t   *g,
              cs_lnum_t   *new_cell_id)
{
  cs_lnum_t ii, jj;
  int rank_id;
  int counts[3];

  int *recv_count = NULL;
  cs_lnum_t *new_src_cell_id = NULL, *new_halo_cell_id = NULL;

  cs_halo_t *h = g->_halo;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  const int n_transforms = h->n_transforms;
  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'h';

  /* Remove send indexes, which will be rebuilt */

  h->n_send_elts[0] = 0;
  h->n_send_elts[1] = 0;
  CS_FREE_HD(h->send_list);
  CS_FREE_HD(h->send_index);
  BFT_FREE(h->send_perio_lst);

  if (g->merge_sub_size == 0) {
    _empty_halo(h);
    return;
  }

  /* Adjust numbering */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    int init_rank_id = h->c_domain_rank[rank_id];
    h->c_domain_rank[rank_id]
      = init_rank_id - (init_rank_id % g->next_merge_stride);
  }

  counts[0] = h->n_c_domains;
  counts[1] = h->n_local_elts;
  counts[2] = h->n_elts[0];

  /* Exchange counters needed for concatenation */

  if (g->merge_sub_rank == 0) {
    BFT_MALLOC(recv_count, g->merge_sub_size*3, int);
    for (ii = 0; ii < 3; ii++)
      recv_count[ii] = counts[ii];
    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      MPI_Recv(recv_count + 3*rank_id, 3, MPI_INT, dist_rank,
               tag, comm, &status);
      for (ii = 0; ii < 3; ii++)
        counts[ii] += recv_count[rank_id*3 + ii];
    }
  }
  else
    MPI_Send(counts, 3, MPI_INT, g->merge_sub_root, tag, comm);

  /* In case of periodic transforms, transpose perio list
     so as to have blocs of fixed n_transforms size, easier
     to work with for append + merge operations */

  if (n_transforms > 0) {
    cs_lnum_t *perio_list_tr;
    BFT_MALLOC(perio_list_tr, h->n_c_domains*n_transforms*4, cs_lnum_t);
    for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
      for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
        for (int k = 0; k < 4; k++)
          perio_list_tr[n_transforms*4*rank_id + 4*tr_id + k]
            = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + k];
      }
    }
    memcpy(h->perio_lst,
           perio_list_tr,
           h->n_c_domains*n_transforms*4*sizeof(cs_lnum_t));
    BFT_FREE(perio_list_tr);
  }

  /* Reallocate arrays for receiving rank and append data */

  if (g->merge_sub_rank == 0) {

    BFT_MALLOC(new_src_cell_id, counts[2], cs_lnum_t);
    for (ii = g->n_rows, jj = 0; ii < g->n_cols_ext; ii++, jj++)
      new_src_cell_id[jj] = new_cell_id[ii];

    BFT_REALLOC(h->c_domain_rank, counts[0], int);
    BFT_REALLOC(h->index, counts[0]*2 + 1, cs_lnum_t);
    BFT_REALLOC(h->perio_lst, counts[0]*n_transforms*4, cs_lnum_t);

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      int n_c_domains_r = recv_count[rank_id*3 + 0];
      cs_lnum_t n_recv = recv_count[rank_id*3 + 2];

      cs_lnum_t index_shift = h->index[2*h->n_c_domains];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(h->c_domain_rank + h->n_c_domains, n_c_domains_r,
               MPI_INT, dist_rank, tag, comm, &status);
      MPI_Recv(new_src_cell_id + h->n_elts[0], n_recv,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);
      MPI_Recv(h->index + h->n_c_domains*2, n_c_domains_r*2+1,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);

      for (ii = 0, jj = h->n_c_domains*2;
           ii < n_c_domains_r*2+1;
           ii++, jj++)
        h->index[jj] += index_shift;

      if (n_transforms > 0) {
        MPI_Recv(h->perio_lst + h->n_c_domains*n_transforms*4,
                 n_c_domains_r*n_transforms*4,
                 CS_MPI_LNUM, dist_rank, tag, comm, &status);
        for (ii = 0, jj = h->n_c_domains*n_transforms*4;
             ii < n_c_domains_r*n_transforms*4;
             ii+=2, jj+=2)
          h->perio_lst[jj] += index_shift;
      }

      /* Update halo sizes */

      h->n_local_elts += recv_count[rank_id*3 + 1];
      h->n_c_domains += n_c_domains_r;
      h->n_elts[0] += n_recv;
      h->n_elts[1] = h->n_elts[0];
    }

  }
  else if (g->merge_sub_size > 1) {

    BFT_MALLOC(new_src_cell_id, h->n_elts[0], cs_lnum_t);
    for (ii = g->n_rows, jj = 0; ii < g->n_cols_ext; ii++, jj++)
      new_src_cell_id[jj] = new_cell_id[ii];

    MPI_Send(h->c_domain_rank, h->n_c_domains, MPI_INT,
             g->merge_sub_root, tag, comm);
    MPI_Send(new_src_cell_id, h->n_elts[0], CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    MPI_Send(h->index, h->n_c_domains*2+1, CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);

    if (n_transforms > 0)
      MPI_Send(h->perio_lst, h->n_c_domains*n_transforms*4,
               CS_MPI_LNUM, g->merge_sub_root, tag, comm);

    _empty_halo(h);
  }

  /* Cleanup halo and set pointer for coarsening (sub_root) ranks*/

  if (h != NULL) {

    /* In case of periodic transforms, transpose perio list back to its
       standard order */

    if (n_transforms > 0) {
      cs_lnum_t *perio_list_tr;
      BFT_MALLOC(perio_list_tr, h->n_c_domains*n_transforms*4, cs_lnum_t);
      for (int tr_id = 0; tr_id < n_transforms; tr_id++) {
        for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
          for (int k = 0; k < 4; k++)
            perio_list_tr[h->n_c_domains*4*tr_id + 4*rank_id + k]
              = h->perio_lst[n_transforms*4*rank_id + 4*tr_id + k];
        }
      }
      memcpy(h->perio_lst,
             perio_list_tr,
             h->n_c_domains*n_transforms*4*sizeof(cs_lnum_t));
      BFT_FREE(perio_list_tr);
    }

    /* Now merge data */

    BFT_MALLOC(new_halo_cell_id, h->n_elts[0], cs_lnum_t);

    _merge_halo_data(h,
                     cs_glob_rank_id,
                     counts[1],
                     new_src_cell_id,
                     new_halo_cell_id);

    _rebuild_halo_send_lists(h, g->alloc_mode, new_src_cell_id);

  }

  if (new_src_cell_id != NULL)
    BFT_FREE(new_src_cell_id);

  g->halo = h;
  g->_halo = h;

  /* Finally, update halo section of cell renumbering array */

  if (g->merge_sub_rank == 0) {

    cs_lnum_t n_send = recv_count[2];
    cs_lnum_t send_shift = n_send;

    for (ii = 0; ii < n_send; ii++)
      new_cell_id[g->n_rows + ii] = new_halo_cell_id[ii];

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      n_send = recv_count[rank_id*3 + 2];
      MPI_Send(new_halo_cell_id + send_shift, n_send, CS_MPI_LNUM,
               dist_rank, tag, comm);
      send_shift += n_send;
    }

    BFT_FREE(recv_count);
  }
  else if (g->merge_sub_size > 1)
    MPI_Recv(new_cell_id + g->n_rows, g->n_cols_ext - g->n_rows,
             CS_MPI_LNUM, g->merge_sub_root, tag, comm, &status);

  BFT_FREE(new_halo_cell_id);
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
  int rank_id;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  const cs_lnum_t db_stride = g->db_size * g->db_size;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'c';

  /* Reallocate arrays for receiving rank and append data */

  if (g->merge_sub_rank == 0) {

    g->n_rows = g->merge_cell_idx[g->merge_sub_size];
    g->n_cols_ext = g->n_rows + g->halo->n_elts[0];

    if (g->relaxation > 0) {
      BFT_REALLOC(g->_cell_cen, g->n_cols_ext * 3, cs_real_t);
      BFT_REALLOC(g->_cell_vol, g->n_cols_ext, cs_real_t);
    }

    BFT_REALLOC(g->_da, g->n_cols_ext*db_stride, cs_real_t);

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = (  g->merge_cell_idx[rank_id+1]
                          - g->merge_cell_idx[rank_id]);
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      if (g->relaxation > 0) {
        MPI_Recv(g->_cell_cen + g->merge_cell_idx[rank_id]*3, n_recv*3,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->_cell_vol + g->merge_cell_idx[rank_id], n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);
      }

      MPI_Recv(g->_da + g->merge_cell_idx[rank_id]*db_stride,
               n_recv*db_stride,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

    }

  }
  else {
    if (g->relaxation > 0) {
      MPI_Send(g->_cell_cen, g->n_rows*3, CS_MPI_REAL, g->merge_sub_root,
               tag, comm);
      BFT_FREE(g->_cell_cen);

      MPI_Send(g->_cell_vol, g->n_rows, CS_MPI_REAL, g->merge_sub_root,
               tag, comm);
      BFT_FREE(g->_cell_vol);
    }

    MPI_Send(g->_da, g->n_rows*db_stride, CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_da);

    g->n_rows = 0;
    g->n_cols_ext = 0;
  }

  if (g->relaxation > 0) {
    g->cell_cen = g->_cell_cen;
    g->cell_vol = g->_cell_vol;
  }

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
  if (g->halo != NULL) {

    if (g->relaxation > 0) {
      cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD, g->_cell_cen, 3);
      if (g->halo->n_transforms > 0)
        cs_halo_perio_sync_coords(g->halo, CS_HALO_STANDARD, g->_cell_cen);

      cs_halo_sync_var(g->halo, CS_HALO_STANDARD, g->_cell_vol);
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
 *   g           <-> pointer to grid structure being merged
 *   n_faces     <-- number of faces to append
 *   new_cel_num <-> new cell numbering for local cells
 *                   in: defined for local cells
 *                   out: updated also for halo cells
 *----------------------------------------------------------------------------*/

static void
_append_face_data(cs_grid_t   *g,
                  cs_lnum_t    n_faces,
                  cs_lnum_t   *face_list)
{
  int rank_id;

  cs_lnum_t *recv_count = NULL;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'f';

  /* Exchange counters needed for concatenation */

  if (g->merge_sub_rank == 0) {
    BFT_MALLOC(recv_count, g->merge_sub_size, cs_lnum_t);
    recv_count[0] = g->n_faces;
    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      MPI_Recv(recv_count + rank_id, 1, CS_MPI_LNUM, dist_rank,
               tag, comm, &status);
    }
  }
  else
    MPI_Send(&n_faces, 1, CS_MPI_LNUM, g->merge_sub_root, tag, comm);

  /* Reallocate arrays for receiving rank and append data */

  BFT_FREE(g->coarse_face);

  if (g->merge_sub_rank == 0) {

    cs_lnum_t n_faces_tot = 0;

    for (rank_id = 0; rank_id < g->merge_sub_size; rank_id++)
      n_faces_tot += recv_count[rank_id];

    BFT_REALLOC(g->_face_cell, n_faces_tot, cs_lnum_2_t);

    if (g->symmetric == true)
      BFT_REALLOC(g->_xa, n_faces_tot, cs_real_t);
    else
      BFT_REALLOC(g->_xa, n_faces_tot*2, cs_real_t);

    if (g->relaxation > 0) {

      BFT_REALLOC(g->_face_normal, n_faces_tot*3, cs_real_t);

      BFT_REALLOC(g->_xa0, n_faces_tot, cs_real_t);
      BFT_REALLOC(g->xa0ij, n_faces_tot*3, cs_real_t);

    }

    g->n_faces = recv_count[0];

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = recv_count[rank_id];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_face_cell + g->n_faces, n_recv*2,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);

      if (g->symmetric == true)
        MPI_Recv(g->_xa + g->n_faces, n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);
      else
        MPI_Recv(g->_xa + g->n_faces*2, n_recv*2,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

      if (g->relaxation > 0) {

        MPI_Recv(g->_face_normal + g->n_faces*3, n_recv*3,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->_xa0 + g->n_faces, n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

        MPI_Recv(g->xa0ij + g->n_faces*3, n_recv*3,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

      }

      g->n_faces += recv_count[rank_id];
    }

    BFT_FREE(recv_count);

  }
  else {

    cs_lnum_t face_id = 0;

    /* Filter face connectivity then send it */

    for (face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t p_face_id = face_list[face_id];
      g->_face_cell[face_id][0] = g->_face_cell[p_face_id][0];
      g->_face_cell[face_id][1] = g->_face_cell[p_face_id][1];
      if (g->symmetric == true)
        g->_xa[face_id] = g->_xa[p_face_id];
      else {
        g->_xa[face_id*2] = g->_xa[p_face_id*2];
        g->_xa[face_id*2+1] = g->_xa[p_face_id*2+1];
      }
      if (g->relaxation > 0) {
        g->_face_normal[face_id*3] = g->_face_normal[p_face_id*3];
        g->_face_normal[face_id*3 + 1] = g->_face_normal[p_face_id*3 + 1];
        g->_face_normal[face_id*3 + 2] = g->_face_normal[p_face_id*3 + 2];
        g->_xa0[face_id] = g->_xa0[p_face_id];
        g->xa0ij[face_id*3] = g->xa0ij[p_face_id*3];
        g->xa0ij[face_id*3 + 1] = g->xa0ij[p_face_id*3 + 1];
        g->xa0ij[face_id*3 + 2] = g->xa0ij[p_face_id*3 + 2];
      }
    }

    MPI_Send(g->_face_cell, n_faces*2, CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_face_cell);

    if (g->symmetric == true)
      MPI_Send(g->_xa, n_faces, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    else
      MPI_Send(g->_xa, n_faces*2, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    BFT_FREE(g->_xa);

    if (g->relaxation > 0) {
      MPI_Send(g->_face_normal, n_faces*3, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
      BFT_FREE(g->_face_normal);

      MPI_Send(g->_xa0, n_faces, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
      BFT_FREE(g->_xa0);

      MPI_Send(g->xa0ij, n_faces*3, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
      BFT_FREE(g->xa0ij);
    }

    g->n_faces = 0;
  }

  g->face_cell = (const cs_lnum_2_t  *)(g->_face_cell);
  g->xa = g->_xa;
  if (g->relaxation > 0) {
    g->face_normal = g->_face_normal;
    g->xa0 = g->_xa0;
  }
}

/*----------------------------------------------------------------------------
 * Merge grids on several ranks to a grid on a single rank
 *
 * The coarse grid must previously have been initialized with _coarse_init()
 * and its coarsening indicator determined (at least for the local cells;
 * it is extended to halo cells here if necessary).
 *
 * - Periodic faces are not handled yet
 *
 * parameters:
 *   g            <-- Pointer to grid structure
 *   merge_stride <-- Associated merge stride
 *   verbosity    <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_merge_grids(cs_grid_t  *g,
             int         merge_stride,
             int         verbosity)
{
  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  int i, rank_id, t_id;
  cs_lnum_t j, face_id;
  int base_rank = cs_glob_rank_id;
  cs_lnum_t cell_shift = 0;
  cs_lnum_t n_faces = 0;
  cs_lnum_t *new_cell_id = NULL, *face_list = NULL;
  bool  *halo_cell_flag = NULL;
  MPI_Comm comm = cs_glob_mpi_comm;
  MPI_Status status;

  static const int tag = 'm'+'e'+'r'+'g'+'e';

  if (merge_stride < 2)
    return;

  /* Determine rank in merged group */

  g->merge_sub_size = merge_stride;
  g->merge_stride = g->next_merge_stride;
  g->next_merge_stride *= merge_stride;

  if (base_rank % g->merge_stride != 0) {
    g->merge_sub_size = 0;
    g->merge_sub_root = -1;
    g->merge_sub_rank = -1;
  }
  else {
    int merge_rank = base_rank / g->merge_stride;
    int merge_size = (cs_glob_n_ranks/g->merge_stride);
    int merge_sub_root = (merge_rank / merge_stride) * merge_stride;
    if (cs_glob_n_ranks % g->merge_stride > 0)
      merge_size += 1;
    g->merge_sub_root = merge_sub_root * g->merge_stride;
    g->merge_sub_rank = merge_rank % merge_stride;
    if (merge_sub_root + g->merge_sub_size > merge_size)
      g->merge_sub_size = merge_size - merge_sub_root;
  }

  if (g->next_merge_stride > 1) {
    if (g->n_ranks % merge_stride)
      g->n_ranks = (g->n_ranks/merge_stride) + 1;
    else
      g->n_ranks = (g->n_ranks/merge_stride);
    g->comm = cs_base_get_rank_step_comm_recursive(g->comm, merge_stride);
  }

  if (verbosity > 2) {
    bft_printf("\n"
               "    merging level %2d grid:\n"
               "      merge_stride = %d; new ranks = %d\n",
               g->level, g->merge_stride, g->n_ranks);
    if (verbosity > 3)
      bft_printf("      sub_root = %d; sub_rank = %d; sub_size = %d\n",
                 g->merge_sub_root, g->merge_sub_rank, g->merge_sub_size);
  }

  /* Compute number of coarse cells in new grid;
     local cell numbers will be shifted by cell_shift in the merged grid */

  if (g->merge_sub_size > 1) {

    if (g->merge_sub_rank == 0) {
      BFT_MALLOC(g->merge_cell_idx, g->merge_sub_size + 1, cs_lnum_t);
      g->merge_cell_idx[0] = 0; g->merge_cell_idx[1] = g->n_rows;
      for (i = 1; i < g->merge_sub_size; i++) {
        cs_lnum_t recv_val;
        int dist_rank = g->merge_sub_root + g->merge_stride*i;
        MPI_Recv(&recv_val, 1, CS_MPI_LNUM, dist_rank, tag, comm, &status);
        g->merge_cell_idx[i + 1] = recv_val + g->merge_cell_idx[i];
      }
    }
    else {
      cs_lnum_t send_val = g->n_rows;
      MPI_Send(&send_val, 1, CS_MPI_LNUM, g->merge_sub_root, tag, comm);
    }

    /* Now return computed info to grids that will be merged */

    if (g->merge_sub_rank == 0) {
      for (i = 1; i < g->merge_sub_size; i++) {
        cs_lnum_t send_val;
        send_val = g->merge_cell_idx[i];
        MPI_Send(&send_val, 1, CS_MPI_LNUM,
                 g->merge_sub_root + g->merge_stride*i, tag, comm);
      }
    }
    else
      MPI_Recv(&cell_shift, 1, CS_MPI_LNUM,
               g->merge_sub_root, tag, comm, &status);
  }

  /* Compute and exchange new cell numbers */

  BFT_MALLOC(new_cell_id, g->n_cols_ext, cs_lnum_t);
  for (j = 0; j < g->n_rows; j++)
    new_cell_id[j] = cell_shift + j;
  for (j = g->n_rows; j < g->n_cols_ext; j++)
    new_cell_id[j] = -1;

  cs_halo_sync_untyped(g->halo,
                       CS_HALO_STANDARD,
                       sizeof(cs_lnum_t),
                       new_cell_id);

  /* Now build face filter list (before halo is modified) */

  if (g->merge_sub_size > 1 && g->merge_sub_rank > 0 && g->n_faces > 0) {

    cs_lnum_t n_ghost_cells = g->n_cols_ext - g->n_rows;

    /* Mark which faces should be merged: to avoid duplicates, a face
       connected to a cell on a lower rank in the same merge set is
       discarded, as it has already been accounted for by that rank. */

    BFT_MALLOC(face_list, g->n_faces, cs_lnum_t);
    BFT_MALLOC(halo_cell_flag, n_ghost_cells, bool);
    for (j = 0; j < n_ghost_cells; j++)
      halo_cell_flag[j] = false;

    for (i = 0; i < g->halo->n_c_domains; i++) {
      rank_id = g->halo->c_domain_rank[i];
      if (rank_id >= g->merge_sub_root && rank_id < cs_glob_rank_id) {
        for (t_id = 0; t_id < g->halo->n_transforms; t_id++) {
          int t_shift = 4 * g->halo->n_c_domains * t_id;
          cs_lnum_t t_start = g->halo->perio_lst[t_shift + 4*i];
          cs_lnum_t t_end = t_start + g->halo->perio_lst[t_shift + 4*i + 1];
          for (j = t_start; j < t_end; j++)
            halo_cell_flag[j] = true;
        }
      }
      else {
        for (j = g->halo->index[2*i]; j < g->halo->index[2*i+2]; j++)
          halo_cell_flag[j] = true;
      }
    }

    for (face_id = 0; face_id < g->n_faces; face_id++) {
      bool use_face = true;
      cs_lnum_t ii = g->face_cell[face_id][0] - g->n_rows + 1;
      cs_lnum_t jj = g->face_cell[face_id][1] - g->n_rows + 1;
      if (ii > 0) {
        if (halo_cell_flag[ii - 1] == false)
          use_face = false;
      }
      else if (jj > 0) {
        if (halo_cell_flag[jj - 1] == false)
          use_face = false;
      }
      if (use_face == true)
        face_list[n_faces++] = face_id;
    }

    BFT_FREE(halo_cell_flag);
  }

  /* Append and merge halos */

  _append_halos(g, new_cell_id);

  /* Update face ->cells connectivity */

  for (face_id = 0; face_id < g->n_faces; face_id++) {
    cs_lnum_t ii = g->face_cell[face_id][0];
    cs_lnum_t jj = g->face_cell[face_id][1];
    assert(ii != jj && new_cell_id[ii] != new_cell_id[jj]);
    g->_face_cell[face_id][0] = new_cell_id[ii];
    g->_face_cell[face_id][1] = new_cell_id[jj];
  }

  BFT_FREE(new_cell_id);

  /* Merge cell and face data */

  if (g->merge_sub_size > 1) {
    _append_cell_data(g);
    _append_face_data(g, n_faces, face_list);
  }
  _sync_merged_cell_data(g);

  BFT_FREE(face_list);

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

/*----------------------------------------------------------------------------
 * Scatter coarse cell integer data in case of merged grids
 *
 * parameters:
 *   g    <-- Grid structure
 *   num  <-> Variable defined on merged grid cells in,
 *            defined on scattered to grid cells out
 *----------------------------------------------------------------------------*/

static void
_scatter_row_int(const cs_grid_t  *g,
                 int              *num)
{
  assert(g != NULL);

  /* If grid merging has taken place, scatter coarse data */

  if (g->merge_sub_size > 1) {

    MPI_Comm  comm = cs_glob_mpi_comm;
    static const int tag = 'p'+'r'+'o'+'l'+'o'+'n'+'g';

    /* Append data */

    if (g->merge_sub_rank == 0) {
      int rank_id;
      assert(cs_glob_rank_id == g->merge_sub_root);
      for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
        cs_lnum_t n_send = (  g->merge_cell_idx[rank_id+1]
                            - g->merge_cell_idx[rank_id]);
        int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
        MPI_Send(num + g->merge_cell_idx[rank_id], n_send, CS_MPI_LNUM,
                 dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(num, g->n_elts_r[0], CS_MPI_LNUM,
               g->merge_sub_root, tag, comm, &status);
    }
  }
}

/*----------------------------------------------------------------------------
 * Scatter coarse cell integer data in case of merged grids
 *
 * parameters:
 *   g    <-- Grid structure
 *   num  <-> Variable defined on merged grid cells in,
 *            defined on scattered to grid cells out
 *----------------------------------------------------------------------------*/

static void
_scatter_row_num(const cs_grid_t  *g,
                 cs_lnum_t        *num)
{
  assert(g != NULL);

  /* If grid merging has taken place, scatter coarse data */

  if (g->merge_sub_size > 1) {

    MPI_Comm  comm = cs_glob_mpi_comm;
    static const int tag = 'p'+'r'+'o'+'l'+'o'+'n'+'g';

    /* Append data */

    if (g->merge_sub_rank == 0) {
      int rank_id;
      assert(cs_glob_rank_id == g->merge_sub_root);
      for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
        cs_lnum_t n_send = (  g->merge_cell_idx[rank_id+1]
                            - g->merge_cell_idx[rank_id]);
        int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
        MPI_Send(num + g->merge_cell_idx[rank_id], n_send, CS_MPI_LNUM,
                 dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(num, g->n_elts_r[0], CS_MPI_LNUM,
               g->merge_sub_root, tag, comm, &status);
    }
  }
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Remove an element from a list (queue) of elements listed by adjacency size.
 *
 * parameters:
 *   s      <-> structure associated with graph traversal and update
 *   a_m_e  <-- number of adjacencies considered for this element
 *   elt_id <-- id of element to remove
 *----------------------------------------------------------------------------*/

static inline void
_graph_m_ptr_remove_m(cs_graph_m_ptr_t  *s,
                      cs_lnum_t          a_m_e,
                      cs_lnum_t          elt_id)
{
  cs_lnum_t _m = a_m_e;

  cs_lnum_t s_p = s->prev[elt_id];
  cs_lnum_t s_n = s->next[elt_id];
  if (s_p > -1) { /* not head of list */
    s->next[s_p] = s_n;
    if (s_n > -1)
      s->prev[s_n] = s_p;
  }
  else { /* head of list */
    assert(s->m_head[_m] == elt_id);
    s->m_head[_m] = s_n;
    if (s_n > -1)
      s->prev[s_n] = -1;
    if (s->m_head[_m] < 0) { /* Update min and max m */
      if (_m == s->m_min) {
        while (s->m_min < s->m_max && s->m_head[s->m_min] < 0)
          s->m_min++;
      }
      else if (_m == s->m_max) {
        s->m_max -= 1;
        while (s->m_max > s->m_min) {
          if (s->m_head[s->m_max] > -1)
            break;
          else
            s->m_max--;
        }
      }
    }
  }

  s->prev[elt_id] = -1;
  s->next[elt_id] = -1;
}

/*----------------------------------------------------------------------------
 * Insert an element at the head of a list (queue) of elements listed
 * by adjacency size.
 *
 * Elements inserted through this function must have been removed from
 * a list with higher adjacency, so a_m_e < s->m_max.
 *
 * parameters:
 *   s      <-> structure associated with graph traversal and update
 *   a_m_e  <-> number of adjacencies considered for this element
 *   elt_id <-- id of element to remove
 *----------------------------------------------------------------------------*/

static inline void
_graph_m_ptr_insert_m(cs_graph_m_ptr_t  *s,
                      cs_lnum_t          a_m_e,
                      cs_lnum_t          elt_id)
{
  cs_lnum_t _m = a_m_e;

  cs_lnum_t s_p = s->m_head[_m];
  if (s_p > -1) { /* list not empty */
    assert(s->prev[s_p] == -1);
    s->prev[s_p] = elt_id;
    s->next[elt_id] = s_p;
  }
  else {  /* Update min and max m */
    if (_m < s->m_min)
      s->m_min = _m;
    else if (_m > s->m_max) {
      s->m_max = _m;
      if (s->m_head[s->m_min] < 0)
        s->m_min = _m;
    }
    assert(_m <= s->m_max);
  }
  s->m_head[_m] = elt_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reorder fine to coarse row mapping so that coarse rows are
 *        ordered by increasing first associated fine row.
 *
 * If called before _scan_f_c_row, where f_c_row used thread-local ordering
 * this function can be called locally for each thread.
 *
 * \param[in]       f_n_rows  number of fine rows for this thread
 * \param[in]       c_n_rows  number of coarse rows for this thread
 * \param           c_o2n     work array (size c_n_rows) for this thread
 * \param[in, out]  f_c_row   fine to coarse rows mapping for this thread.
 */
/*----------------------------------------------------------------------------*/

static void
_reorder_f_c_row(cs_lnum_t    f_n_rows,
                 cs_lnum_t    c_n_rows,
                 cs_lnum_t   *restrict c_o2n,
                 cs_lnum_t   *restrict f_c_row)
{
  for (cs_lnum_t i = 0; i < c_n_rows; i++) {
    c_o2n[i] = -1;
  }

  cs_lnum_t c_count = 0;
  for (cs_lnum_t i = 0; i < f_n_rows; i++) {
    cs_lnum_t c_row = f_c_row[i];
    if (c_row > -1) {
      if (c_o2n[c_row] < 0) {
        c_o2n[c_row] = c_count;
        f_c_row[i] = c_count;
        c_count++;
      }
      else
        f_c_row[i] = c_o2n[c_row];
    }
  }

  assert(c_count == c_n_rows);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update fine-to coarse row mapping using a scan type operation
 *        when that mapping has been computed locally on separate threads.
 *
 * \warning This operator calls  \ref cs_parall_thread_range with f_n_rows,
 * so the caller or upstream algorithm should have used to same thread
 * ranges.
 *
 * \param[in]       n_threads  number of threads
 * \param[in]       f_n_rows   number of rows in fine grid
 * \param[in, out]  t_c_scan   coarse row count for each thread
 * \param[in, out]  f_c_row    fine to coarse rows mapping (per-thread in,
 *                             shared out).
 *
 * \return  local number of resulting coarse rows
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_OPENMP)

static cs_lnum_t
_scan_f_c_row(int          n_threads,
              cs_lnum_t    f_n_rows,
              cs_lnum_t    t_c_scan[],
              cs_lnum_t   *restrict f_c_row)
{
  cs_lnum_t c_n_rows = 0;

  for (int i = 0; i < n_threads; i++) {
    cs_lnum_t n_t_rows = t_c_scan[i];
    t_c_scan[i] = c_n_rows;
    c_n_rows += n_t_rows;
  }

  #pragma omp parallel  num_threads(n_threads)
  {
    const int t_id = omp_get_thread_num();
    cs_lnum_t r_shift = t_c_scan[t_id];

    if (r_shift != 0) { /* No shift needed for thread 0 */
      cs_lnum_t t_s_id, t_e_id;
      cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        if (f_c_row[ii] > -1)
          f_c_row[ii] += r_shift;
      }
    }
  }

  return c_n_rows;
}

#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply one step of the pairwise aggregation algorithm for a
 *        matrix expected to be an M-matrix.
 *
 * This assumes positive diagonal entries, and negative off-diagonal entries.
 * The matrix is provided in scalar MSR (modified CSR with separate diagonal)
 * form. To use block matrices with this function, their blocks must be
 * condensed to equivalent scalars.
 *
 * \param[in]   f_n_rows      number of rows in fine grid
 * \param[in]   beta          aggregation criterion
 * \param[in]   dd_threshold  diagonal dominance threshold; if > 0, ignore rows
 *                            whose diagonal dominance is above this threshold
 * \param[in]   row_index     matrix row index (separate diagonal)
 * \param[in]   col_id        matrix column ids
 * \param[in]   d_val         matrix diagonal values (scalar)
 * \param[in]   x_val         matrix extradiagonal values (scalar)
 * \param[out]  f_c_row       fine to coarse rows mapping
 *
 * \return  local number of resulting coarse rows
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_pairwise_msr(cs_lnum_t                  f_n_rows,
              const cs_real_t            beta,
              const cs_real_t            dd_threshold,
              const cs_lnum_t  *restrict row_index,
              const cs_lnum_t  *restrict col_id,
              const cs_real_t  *restrict d_val,
              const cs_real_t  *restrict x_val,
              cs_lnum_t        *restrict f_c_row)
{
  cs_lnum_t c_n_rows = 0;

  /* Coarse row count scan array for shared-memory parallelism */

  cs_lnum_t *t_c_scan = nullptr;
  int n_loc_threads = cs_parall_n_threads(f_n_rows, CS_THR_MIN);
  if (n_loc_threads > 1) {
    BFT_MALLOC(t_c_scan, n_loc_threads, cs_lnum_t);
    for (int i = 0; i < n_loc_threads; i++)
      t_c_scan[i] = 0;
  }

  /* Allocate working arrays */

  short int  *a_m;   /* active m for row */
  cs_real_t  *a_max; /* max per line */
  cs_lnum_t *s_next;

  BFT_MALLOC(a_m, f_n_rows, short int);
  BFT_MALLOC(a_max, f_n_rows, cs_real_t);
  BFT_MALLOC(s_next, f_n_rows*2, cs_lnum_t);

  #pragma omp parallel shared(c_n_rows)  num_threads(n_loc_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    cs_lnum_t t_f_n_rows = t_e_id - t_s_id; // Number of rows for this thread.
    cs_lnum_t t_c_n_rows = 0;

    /* Mark all elements of fine to coarse rows as uninitialized */
    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++)
      f_c_row[ii] = -2;

    /* Computation of the maximum over line ii and test if the line ii is
     * ignored. Be careful that the sum has to be the sum of the
     * absolute value of every extra-diagonal coefficient, but the maximum
     * is only on the negative coefficient. */

    cs_lnum_t m_max = -1;
    cs_lnum_t n_remain = t_f_n_rows;

    if (dd_threshold > 0) {

      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {

        cs_real_t sum = 0.0;

        cs_lnum_t s_id = row_index[ii];
        cs_lnum_t e_id = row_index[ii+1];
        a_max[ii] = 0.0;

        for (cs_lnum_t jj = s_id; jj < e_id; jj++) {
          cs_real_t xv = x_val[jj];
          sum += CS_ABS(xv);
          if (xv < 0)
            a_max[ii] = cs_math_fmax(a_max[ii], -xv);
        }

        /* Check if the line seems ignored or not */

        if (d_val[ii] > dd_threshold * sum) {
          a_m[ii] = -1;
          n_remain -= 1;
          f_c_row[ii] = -1;
        }
        else {
          a_m[ii] = 0;
          for (cs_lnum_t jj = e_id-1; jj >= s_id; jj--) {
            if (col_id[jj] < f_n_rows && x_val[jj] < beta*sum)
              a_m[ii] += 1;
          }
        }

        if (m_max < a_m[ii])
          m_max = a_m[ii];

      }

    }
    else { /* variant with no diagonal dominance check */

      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {

        cs_lnum_t s_id = row_index[ii];
        cs_lnum_t e_id = row_index[ii+1];
        a_max[ii] = 0.0;

        for (cs_lnum_t jj = s_id; jj < e_id; jj++) {
          cs_real_t xv = x_val[jj];
          if (xv < 0)
            a_max[ii] = cs_math_fmax(a_max[ii], -xv);
        }

        a_m[ii] = 0;
        for (cs_lnum_t jj = e_id-1; jj >= s_id; jj--) {
          if (col_id[jj] < f_n_rows)
            a_m[ii] += 1;
        }

        if (m_max < a_m[ii])
          m_max = a_m[ii];

      }

    }

    if (m_max >= 0) {

      /* Build pointers to lists of rows by a_m
         (to allow access to row with lowest m) */

      cs_graph_m_ptr_t s;

      s.m_min = 0;
      s.m_max = m_max;

      BFT_MALLOC(s.m_head, s.m_max+1, cs_lnum_t);
      s.next = s_next;
      s.prev = s.next + f_n_rows;

      for (cs_lnum_t ii = 0; ii < s.m_max+1; ii++)
        s.m_head[ii] = -1;
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        s.next[ii] = -1;
        s.prev[ii] = -1;
      }

      for (cs_lnum_t ii = t_e_id-1; ii >= t_s_id; ii--) {
        short int _m = a_m[ii];
        if (_m >= 0) {
          cs_lnum_t prev_head = s.m_head[_m];
          s.m_head[_m] = ii;
          s.next[ii] = prev_head;
          if (prev_head > -1)
            s.prev[prev_head] = ii;
        }
      }

      /* Now build pairs */

      while (s.m_min < s.m_max) {
        if (s.m_head[s.m_min] < 0)
          s.m_min++;
        else
          break;
      }

      while (s.m_min > -1) {

        /* Select remaining ii with minimal a_m */

        cs_lnum_t ii = s.m_head[s.m_min];
        assert(ii >= t_s_id && ii < t_e_id);

        cs_lnum_t gg[2] = {ii, -1};

        /* Select remaining jj such that aij = min_over_k_aik */

        cs_lnum_t s_id = row_index[ii];
        cs_lnum_t e_id = row_index[ii+1];

        f_c_row[ii] = t_c_n_rows; /* Add i to "pair" in all cases */

        if (e_id > s_id) {

          cs_lnum_t jj = -1;
          cs_real_t _a_min = HUGE_VAL;
          for (cs_lnum_t kk_idx = s_id; kk_idx < e_id; kk_idx++) {
            cs_lnum_t kk = col_id[kk_idx];
            /* not aggregated yet */
            if (kk >= t_s_id && kk < t_e_id && f_c_row[kk] == -2) {
              cs_real_t xv = x_val[kk_idx];
              if (xv < _a_min) {
                _a_min = xv;
                jj = kk;
              }
            }
          }

          /* Keep jj only if within threshold */

          if (_a_min >= -beta*a_max[ii])
            jj = -1;
          else {
            f_c_row[jj] = t_c_n_rows; /* Add jj to "pair" */
            gg[1] = jj;
          }

        }

        t_c_n_rows++;

        /* Now update search set */

        for (int ip = 0; ip < 2; ip++) {
          cs_lnum_t i = gg[ip];
          if (i >= t_s_id && i < t_e_id) {
            cs_lnum_t _m = a_m[i];
            if (_m < 0)
              continue;

            _graph_m_ptr_remove_m(&s, _m, i);
            a_m[i] = -1;
            cs_lnum_t _s_id = row_index[i];
            cs_lnum_t _e_id = row_index[i+1];
            for (cs_lnum_t k = _e_id-1; k >= _s_id; k--) {
              cs_lnum_t j = col_id[k];
              if (j < t_s_id || j >= t_e_id)
                continue;
              _m = a_m[j];
              if (_m >= 0) {
                _graph_m_ptr_remove_m(&s, _m, j);
                if (_m > 0)
                  _graph_m_ptr_insert_m(&s, _m-1, j);
                a_m[j] = _m - 1;
              }
            }
            n_remain--;
          } /* i > -1 */
        }
        if (s.m_min >= s.m_max) { /* check if list has become empty */
          if (s.m_head[s.m_min] < 0) {
            s.m_min = -1;
            s.m_max = -1;
          }
        }

      }

      /* We might have remaining cells */

      for (int ii = t_s_id; ii < t_e_id; ii++) {
        if (f_c_row[ii] < -1)
          f_c_row[ii] = t_c_n_rows++;
      }

      BFT_FREE(s.m_head);

      /* Prepare combining thread results */

#if defined(HAVE_OPENMP)
      if (n_loc_threads > 1)
        t_c_scan[omp_get_thread_num()] = t_c_n_rows;
      else
        c_n_rows = t_c_n_rows;
#else
      c_n_rows = t_c_n_rows;
#endif

    } //  if (m_max >= 0)

  } // end of OpenMP section

  /* Combine thread results */

#if defined(HAVE_OPENMP)
  if (n_loc_threads > 1) {
    c_n_rows = _scan_f_c_row(n_loc_threads, f_n_rows, t_c_scan, f_c_row);
  }
#endif

  /* Free working arrays */

  BFT_FREE(s_next);
  BFT_FREE(a_max);
  BFT_FREE(a_m);

  BFT_FREE(t_c_scan);

  return c_n_rows;
}

 /*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using an
 * automatic criterion and pairwise aggregation variant 1,
 * with a matrix in MSR format.
 *
 * parameters:
 *   f                   <-- Fine grid structure
 *   verbosity           <-- Verbosity level
 *   f_c_row             --> Fine row -> coarse row connectivity
 *
 * return:
 *   number of coarse rows
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_automatic_aggregation_pw_msr(const cs_grid_t  *f,
                              int               verbosity,
                              cs_lnum_t        *f_c_row)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const cs_real_t beta = 0.25;
  const cs_real_t dd_threshold = (f->level == 0) ? _dd_threshold_pw : -1;
  const cs_lnum_t f_n_rows = f->n_rows;

  /* Access matrix MSR vectors */

  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;
  cs_real_t *_d_val = NULL, *_x_val = NULL;

  cs_matrix_get_msr_arrays(f->matrix,
                           &row_index,
                           &col_id,
                           &d_val,
                           &x_val);

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t eb_size = f->eb_size;

  if (db_size > 1) {
    BFT_MALLOC(_d_val, f_n_rows, cs_real_t);
    _reduce_block(f_n_rows, db_size, d_val, _d_val);
    d_val = _d_val;
  }

  if (eb_size > 1) {
    cs_lnum_t f_n_enz = row_index[f_n_rows];
    BFT_MALLOC(_x_val, f_n_enz, cs_real_t);
    _reduce_block(f_n_enz, eb_size, x_val, _x_val);
    x_val = _x_val;
  }

  if (verbosity > 3)
    bft_printf("\n     %s: beta %5.3e; diag_dominance_threshold: %5.3e\n",
               __func__, beta, dd_threshold);

  cs_lnum_t c_n_rows = _pairwise_msr(f_n_rows,
                                     beta,
                                     dd_threshold,
                                     row_index,
                                     col_id,
                                     d_val,
                                     x_val,
                                     f_c_row);

  /* Free working arrays */

  BFT_FREE(_d_val);
  BFT_FREE(_x_val);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, f->level-1);
    printf(", total = %ld\n", elapsed.count());
  }

  return c_n_rows;
}

/*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using an
 * automatic criterion and pairwise aggregation variant 1,
 * with a matrix in native format.
 *
 * parameters:
 *   f                   <-- Fine grid structure
 *   max_aggregation     <-- Max fine rows per coarse row
 *   verbosity           <-- Verbosity level
 *   f_c_row             --> Fine row -> coarse row connectivity
 *
 * return:
 *   number of coarse rows
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_automatic_aggregation_mx_native(const cs_grid_t  *f,
                                 cs_lnum_t         max_aggregation,
                                 int               verbosity,
                                 cs_lnum_t        *f_c_row)
{
  const cs_lnum_t f_n_rows = f->n_rows;

  /* Algorithm parameters */
  const cs_real_t beta = 0.25; /* 0.5 for HHO */
  const int ncoarse = 8;
  int npass_max = 10;

  int _max_aggregation = 1, npass = 0;
  cs_lnum_t aggr_count = f_n_rows;
  cs_lnum_t c_n_rows = 0;

  cs_lnum_t *c_aggr_count = NULL;
  cs_real_t *maxi = NULL;

  /* Access matrix native arrays */

  cs_lnum_t n_edges = 0;
  const cs_lnum_2_t  *edges;
  const cs_real_t  *d_val, *x_val;
  cs_real_t *_d_val = NULL, *_x_val = NULL;

  bool symmetric = true;
  cs_lnum_t isym = 2;

  cs_matrix_get_native_arrays(f->matrix,
                              &symmetric,
                              &n_edges,
                              &edges,
                              &d_val,
                              &x_val);

  if (f->symmetric)
    isym = 1;

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t eb_size = f->eb_size;

  if (db_size > 1) {
    BFT_MALLOC(_d_val, f_n_rows, cs_real_t);
    _reduce_block(f_n_rows, db_size, d_val, _d_val);
    d_val = _d_val;
  }

  if (eb_size > 1) {
    BFT_MALLOC(_x_val, n_edges*isym, cs_real_t);
    _reduce_block(n_edges*isym, eb_size, x_val, _x_val);
    x_val = _x_val;
  }

  /* Allocate working arrays */

  BFT_MALLOC(c_aggr_count, f_n_rows, cs_lnum_t);
  BFT_MALLOC(maxi, f_n_rows, cs_real_t);

  for (cs_lnum_t ii = 0; ii < f_n_rows; ii++){
    c_aggr_count[ii] = 1;
  }

  /* Computation of the maximum over line ii and test if the line ii is
   * penalized. Be careful that the sum has to be the sum of the
   * absolute value of every extra-diagonal coefficient, but the maximum
   * is only on the negative coefficient. */

  cs_real_t *sum;
  BFT_MALLOC(sum, f_n_rows, cs_real_t);

  for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
    sum[ii] = 0;
    maxi[ii] = 0.0;
  }

  for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {
    cs_lnum_t ii = edges[e_id][0];
    cs_lnum_t jj = edges[e_id][1];
    cs_real_t xv0 = x_val[e_id*isym];
    cs_real_t xv1 = x_val[(e_id+1)*isym-1];
    if (ii < f_n_rows) {
      sum[ii] += CS_ABS(xv0);
      if (xv0 < 0.)
        maxi[ii] = cs_math_fmax(maxi[ii], -xv0);
    }
    if (jj < f_n_rows) {
      sum[jj] += CS_ABS(xv1);
      if (xv1 < 0.)
        maxi[jj] = cs_math_fmax(maxi[jj], -xv1);
    }
  }

  /* Check if the line seems penalized or not. */
  if (f->level == 0) {
    for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
      if (d_val[ii] > _penalization_threshold * sum[ii])
        f_c_row[ii] = -1;
    }
  }
  BFT_FREE(sum);

  /* Passes */

  if (verbosity > 3)
    bft_printf("\n     %s:\n", __func__);

  do {

    npass++;
    _max_aggregation++;
    _max_aggregation = CS_MIN(_max_aggregation, max_aggregation);

    /* Pairwise aggregation */

    for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

      cs_lnum_t ii = edges[e_id][0];
      cs_lnum_t jj = edges[e_id][1];

      /* Exclude rows on parallel or periodic boundary, so as not to */
      /* coarsen the grid across those boundaries (which would change */
      /* the communication pattern and require a more complex algorithm). */

      /* ii or jj is candidate to aggregation only if it is not penalized */

      if (ii >= f_n_rows || jj >= f_n_rows)
        continue;

      cs_lnum_t ii_c = f_c_row[ii];
      cs_lnum_t jj_c = f_c_row[jj];

      if (ii_c == -1 || jj_c == -1)  /* penalization */
        continue;

      cs_real_t xv = x_val[e_id*isym];

      if (isym == 2)
        xv = cs_math_fmax(xv, x_val[e_id*2 + 1]);

      /* Test if ii and jj are strongly negatively coupled and at */
      /* least one of them is not already in an aggregate. */

      if (xv < -beta*maxi[ii] && (ii_c == -2 || jj_c == -2)) {

        if (ii_c > -1 && jj_c < -1) {
          if (c_aggr_count[ii_c] < _max_aggregation +1) {
            f_c_row[jj] = ii_c;
            c_aggr_count[ii_c] += 1;
          }
        }
        else if (ii_c < -1 && jj_c > -1) {
          if (c_aggr_count[jj_c] < _max_aggregation +1) {
            f_c_row[ii] = jj_c;
            c_aggr_count[jj_c] += 1;
          }
        }
        else if (ii_c < -1 && jj_c < -1) {
          f_c_row[ii] = c_n_rows;
          f_c_row[jj] = c_n_rows;
          c_aggr_count[c_n_rows] += 1;
          c_n_rows++;
        }
      }

    }

    /* Check the number of coarse rows created */
    aggr_count = 0;
    for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
      if (f_c_row[ii] < -1)
        aggr_count++;
    }

    /* Additional passes if aggregation is insufficient */
    if (aggr_count == 0 || (c_n_rows + aggr_count)*ncoarse < f_n_rows)
      npass_max = npass;

  } while (npass < npass_max); /* Loop on passes */

  /* Finish assembly: rows that are not diagonally dominant and not in an
   * aggregate form their own aggregate */
  for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
    if (f_c_row[ii] < -1) {
      f_c_row[ii] = c_n_rows;
      c_n_rows++;
    }
  }

  /* Free working arrays */

  BFT_FREE(_d_val);
  BFT_FREE(_x_val);
  BFT_FREE(c_aggr_count);
  BFT_FREE(maxi);

  return c_n_rows;
}

/*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using an
 * automatic criterion and pairwise aggregation variant 1,
 * with a matrix in MSR format.
 *
 * parameters:
 *   f                   <-- Fine grid structure
 *   max_aggregation     <-- Max fine rows per coarse row
 *   verbosity           <-- Verbosity level
 *   f_c_row             --> Fine row -> coarse row connectivity
 *
 * return:
 *   number of coarse rows
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_automatic_aggregation_mx_msr(const cs_grid_t  *f,
                              cs_lnum_t         max_aggregation,
                              int               verbosity,
                              cs_lnum_t        *f_c_row)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const cs_lnum_t f_n_rows = f->n_rows;
  cs_lnum_t c_n_rows = -1;

  const int npass_max = 10;

  cs_lnum_t *c_aggr_count = NULL;
  cs_real_t *maxi = NULL;

  /* Algorithm parameters */
  const cs_real_t beta = 0.25; /* 0.5 for HHO */
  const int ncoarse = 8;

  if (verbosity > 3)
    bft_printf("\n     %s: npass_max: %d; n_coarse: %d;"
               " beta %5.3e; pena_thd: %5.3e\n",
               __func__, npass_max, ncoarse, beta, _penalization_threshold);

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t eb_size = f->eb_size;

  /* Access matrix MSR vectors */

  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;
  cs_real_t *_d_val = NULL, *_x_val = NULL;

  cs_matrix_get_msr_arrays(f->matrix,
                           &row_index,
                           &col_id,
                           &d_val,
                           &x_val);

  const cs_lnum_t f_nnz = row_index[f_n_rows];

  if (db_size > 1) {
    BFT_MALLOC(_d_val, f_n_rows, cs_real_t);
    _reduce_block(f_n_rows, db_size, d_val, _d_val);
    d_val = _d_val;
  }

  if (eb_size > 1) {
    cs_lnum_t f_n_enz = row_index[f_n_rows];
    BFT_MALLOC(_x_val, f_n_enz, cs_real_t);
    _reduce_block(f_n_enz, eb_size, x_val, _x_val);
    x_val = _x_val;
  }

  /* Coarse row count scan array for shared-memory parallelism */

  cs_lnum_t *t_c_scan = nullptr;
  int n_loc_threads = cs_parall_n_threads(f_n_rows, CS_THR_MIN);
  if (n_loc_threads > 1) {
    BFT_MALLOC(t_c_scan, n_loc_threads, cs_lnum_t);
    for (int i = 0; i < n_loc_threads; i++)
      t_c_scan[i] = 0;
  }

  /* Allocate working arrays */

  BFT_MALLOC(c_aggr_count, f_n_rows, cs_lnum_t);
  BFT_MALLOC(maxi, f_n_rows, cs_real_t);

  /* aggregation queue: local column id, index in matrix, and
     index of symmetric element if needed */
  cs_lnum_t *ag_queue;
  BFT_MALLOC(ag_queue, f_nnz*2, cs_lnum_t);

  /* Handle a block of rows per thread; aggregation will not be done
     across thread blocks, in a similar manner that it is not done
     across MPI ranks. */

  #pragma omp parallel shared(c_n_rows) num_threads(n_loc_threads)
  {

#if defined(HAVE_OPENMP)
    const int t_id = omp_get_thread_num();
#else
    const int t_id = 0;
#endif
    /* Use a simple thread range for now. Perhaps a more complex range
       computation could lead to better work balance between threads,
       as rows with low ids will tend to get more work than rows
       with high ids, given that we consider the upper triangular part */

    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    int _max_aggregation = 1;

    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      c_aggr_count[ii] = 1;
    }

    /* Computation of the maximum over line ii and test if the line ii is
     * penalized. Be careful that the sum has to be the sum of the absolute
     * value of every extra-diagonal coefficient, but the maximum is only on the
     * negative coefficient. */

    if (f->level == 0) {
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        maxi[ii] = 0.0;
        cs_real_t  sum = 0.0;

        for (cs_lnum_t jj = row_index[ii]; jj < row_index[ii+1]; jj++) {
          const cs_real_t  xv = x_val[jj];
          if (xv < 0.) {
            sum -= xv;
            maxi[ii] = cs_math_fmax(maxi[ii], -xv);
          }
          else {
            sum += xv;
          }
        }

        /* Check if the line seems penalized or not. */
        if (d_val[ii] > _penalization_threshold * sum)
          f_c_row[ii] = -1;
      }
    }
    else {
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        maxi[ii] = 0.0;
        for (cs_lnum_t jj = row_index[ii]; jj < row_index[ii+1]; jj++) {
          const cs_real_t  xv = x_val[jj];
          if (xv < 0.)
            maxi[ii] = cs_math_fmax(maxi[ii], -xv);
        }
      }
    }

    /* Indexes and counters */

    const cs_lnum_t t_r_s_id = row_index[t_s_id] * 2;
    cs_lnum_t t_r_e_id = t_r_s_id;

    cs_lnum_t t_c_n_rows = 0;
    cs_lnum_t r_n_faces = 0;

    /* First pass */

    {
      if (_max_aggregation < max_aggregation)
        _max_aggregation++;

      /* Pairwise aggregation */

      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {

        /* ii is candidate to aggregation only if it is not penalized */

        cs_lnum_t ii_c = f_c_row[ii];

        if (ii_c == -1) /* penalized */
          continue;

        const cs_real_t  row_criterion = -beta*maxi[ii];

        for (cs_lnum_t jidx = row_index[ii]; jidx < row_index[ii+1]; jidx++) {

          /* Exclude rows on parallel or periodic boundary, so as not to */
          /* coarsen the grid across those boundaries (which would change */
          /* the communication pattern and require a more complex algorithm). */

          cs_lnum_t jj = col_id[jidx];
          if (jj < t_s_id || jj >= t_e_id)
            continue;

          cs_lnum_t jj_c = f_c_row[jj];

          if (   (ii_c > -1 && jj_c > -1)   /* both already aggregated */
              || jj_c == -1)                /* one penalized */
            continue;

          /* Test if ii and jj are strongly negatively coupled */

          int count = 0;

          if (x_val[jidx] < row_criterion) {

            if (ii_c > -1 && jj_c < -1) {
              if (c_aggr_count[t_s_id + ii_c] < _max_aggregation +1) {
                f_c_row[jj] = ii_c;
                c_aggr_count[t_s_id + ii_c] += 1;
                count++;
              }
            }
            else if (ii_c < -1 && jj_c > -1) {
              if (c_aggr_count[t_s_id + jj_c] < _max_aggregation +1) {
                f_c_row[ii] = jj_c;
                c_aggr_count[t_s_id + jj_c] += 1;
                count++;
              }
            }
            else if (ii_c < -1 && jj_c < -1) {
              f_c_row[ii] = t_c_n_rows;
              f_c_row[jj] = t_c_n_rows;
              c_aggr_count[t_s_id + t_c_n_rows] += 1;
              t_c_n_rows++;
              count++;
            }
          }

          if (count == 0) {
            ag_queue[t_r_e_id] = ii;
            ag_queue[t_r_e_id+1] = jidx;
            t_r_e_id += 2;
            r_n_faces++;
          }

        } /* Loop on columns */

      } /* Loop on rows */
    }

    /* Other passes */

    for (int npass = 2; npass <= npass_max; npass++) {

      /* Exit loop on passes if aggregation is sufficient
         (note that each thread may loop independently). */

      cs_lnum_t aggr_count = 0; /* Count of cells not yet aggregated */
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        if (f_c_row[ii] < -1)
          aggr_count++;
      }

      if (   aggr_count == 0
          || (t_c_n_rows + aggr_count)*ncoarse < (t_e_id - t_s_id)
          || r_n_faces == 0) {
        break;
      }

      if (_max_aggregation < max_aggregation)
        _max_aggregation++;

      /* Re-initialize non-eliminated faces */
      r_n_faces = 0;

      /* Loop on non-eliminated faces (i.e. "aggregation queue").

         Faces adjacent to penalized rows, or row on parallel boundaries,
         were already excluded from the aggregation queue in the first pass,
         so we do not need to test for those anymore. */

      cs_lnum_t r_idx_n = t_r_s_id;

      for (cs_lnum_t r_idx = t_r_s_id; r_idx < t_r_e_id; r_idx += 2) {

        cs_lnum_t ii = ag_queue[r_idx];
        cs_lnum_t jidx = ag_queue[r_idx+1];

        cs_lnum_t jj = col_id[jidx];

        /* If both cells are already aggregated, pass */

        cs_lnum_t ii_c = f_c_row[ii];
        cs_lnum_t jj_c = f_c_row[jj];

        if (ii_c >= 0 && jj_c >= 0)
          continue;

        /* Test if ii and jj are strongly negatively coupled  */

        const cs_real_t  row_criterion = -beta*maxi[ii];

        int count = 0;

        if (x_val[jidx] < row_criterion) {

          if (ii_c > -1 && jj_c < -1) {
            if (c_aggr_count[t_s_id + ii_c] < _max_aggregation +1) {
              f_c_row[jj] = ii_c;
              c_aggr_count[t_s_id + ii_c] += 1;
              count++;
            }
          }
          else if (ii_c < -1 && jj_c > -1) {
            if (c_aggr_count[t_s_id + jj_c] < _max_aggregation +1) {
              f_c_row[ii] = jj_c;
              c_aggr_count[t_s_id + jj_c] += 1;
              count++;
            }
          }
          else if (ii_c < -1 && jj_c < -1) {
            f_c_row[ii] = t_c_n_rows;
            f_c_row[jj] = t_c_n_rows;
            c_aggr_count[t_s_id + t_c_n_rows] += 1;
            t_c_n_rows++;
            count++;
          }
        }

        /* Place faces where no merging occurs back in queue
           (we reuse the same location, knowing that
           t_r_s_id + n_r_faces*ag_queue_stride < r_id) */
        if (count == 0) {
          assert(r_idx_n < t_r_e_id);
          ag_queue[r_idx_n] = ii;
          ag_queue[r_idx_n+1] = jidx;
          r_idx_n += 0;
          r_n_faces++;
        }

      } /* Loop on aggregation queue */

      t_r_e_id = r_idx_n;

    } /* Loop on passes */

    /* Finish assembly: rows that are not diagonally dominant and not in an
     * aggregate form their own aggregate */
    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      if (f_c_row[ii] < -1) {
        f_c_row[ii] = t_c_n_rows;
        t_c_n_rows++;
      }
    }

    /* Prepare combining thread results */

    if (n_loc_threads > 1)
      t_c_scan[t_id] = t_c_n_rows;
    else
      c_n_rows = t_c_n_rows;

  } // End of OpenMP section

  /* Combine thread results */

#if defined(HAVE_OPENMP)
  if (n_loc_threads > 1) {
    c_n_rows = _scan_f_c_row(n_loc_threads, f_n_rows, t_c_scan, f_c_row);
  }
#endif

  /* Free working arrays */

  BFT_FREE(ag_queue);
  BFT_FREE(c_aggr_count);
  BFT_FREE(_d_val);
  BFT_FREE(_x_val);
  BFT_FREE(maxi);
  BFT_FREE(t_c_scan);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, f->level-1);
    printf(", total = %ld\n", elapsed.count());
  }

  return c_n_rows;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute coarse to fine grid adjacency.
 *
 * If n_f_threads > 1, it must be consistent with the value used for
 * coarsening, because in that case we assume that each coarse row is
 * assigned to a single thread.
 *
 * \param[in]   f_n_rows       number of rows in fine grid
 * \param[in]   c_n_rows       number of rows in fine grid
 * \param[in]   alloc_mode     allocation mode for cf_row_index and cs_row_ids
 * \param[in]   n_f_threads    number of threads for fine rows
 * \param[in]   f_c_row        fine to coarse row map
 * \param[in]   f_row_index    MSR row index for fine grid (to count columns)
 * \param[out]  c_f_row_index  coarse to fine rows index
 * \param[out]  c_f_row_ids    coarse to fine row ids
 * \param[out]  c_row_index_0  index based on maximum column per coarse row,
 *                             or nullptr
 */
/*----------------------------------------------------------------------------*/

static void
_coarse_to_fine_adjacency_msr(cs_lnum_t         f_n_rows,
                              cs_lnum_t         c_n_rows,
                              cs_alloc_mode_t   alloc_mode,
                              int               n_f_threads,
                              const cs_lnum_t  *restrict f_c_row,
                              const cs_lnum_t  *restrict f_row_index,
                              cs_lnum_t        **c_f_row_index,
                              cs_lnum_t        **c_f_row_ids,
                              cs_lnum_t        **c_row_index_0)
{
  int n_c_threads = cs_parall_n_threads(c_n_rows, CS_THR_MIN);

  cs_lnum_t *cf_r_idx, *c_r_idx_0, *cf_r_shift;
  CS_MALLOC_HD(cf_r_idx, c_n_rows+1, cs_lnum_t, alloc_mode);
  BFT_MALLOC(c_r_idx_0, c_n_rows+1, cs_lnum_t);
  BFT_MALLOC(cf_r_shift, c_n_rows, cs_lnum_t);

# pragma omp parallel for  num_threads(n_c_threads)
  for (cs_lnum_t i = 0; i < c_n_rows; i++) {
    cf_r_idx[i] = 0;
    cf_r_shift[i] = 0;
    c_r_idx_0[i] = 0;
  }
  cf_r_idx[c_n_rows] = 0;
  c_r_idx_0[c_n_rows] = 0;

  /* First loop for counting */

  cs_lnum_t  c_nnz_0 = 0;

  #pragma omp parallel shared(c_nnz_0) num_threads(n_f_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    cs_lnum_t *restrict cf_count = cf_r_idx+1;
    cs_lnum_t *restrict c_r_count = c_r_idx_0+1;

    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      cs_lnum_t jj = f_c_row[ii];
      if (jj > -1) {
        cf_count[jj] += 1;
        c_r_count[jj] += f_row_index[ii+1] - f_row_index[ii];
      }
    }
  }

  /* Pre-scan for coarse rows */

  cs_lnum_t *t_c_scan = nullptr;
  if (n_c_threads > 1) {
    BFT_MALLOC(t_c_scan, n_f_threads*2, cs_lnum_t);
    for (int i = 0; i < n_f_threads; i++) {
      t_c_scan[i*2] = 0;
      t_c_scan[i*2+1] = 0;
    }
  }

  #pragma omp parallel shared(c_nnz_0) num_threads(n_c_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(c_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);
    t_s_id += 1;  /* Shift by 1 as we are dealing with an index */
    t_e_id += 1;

    cs_lnum_t t_c_n_rows = 0;
    cs_lnum_t t_c_nnz_0 = 0;

    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      t_c_n_rows += cf_r_idx[ii];
      cf_r_idx[ii] = t_c_n_rows;

      t_c_nnz_0 += c_r_idx_0[ii];
      c_r_idx_0[ii] = t_c_nnz_0;
    }

#if defined(HAVE_OPENMP)
    if (t_c_scan != nullptr) {
      int t_id = omp_get_thread_num();
      t_c_scan[t_id*2] = t_c_n_rows;
      t_c_scan[t_id*2+1] = t_c_nnz_0;
    }
#endif
  }

#if defined(HAVE_OPENMP)

  /* Scan operation for threads */

  if (t_c_scan != nullptr) {

    c_nnz_0 = 0;
    cs_lnum_t cf_count = 0;
    for (int i = 0; i < n_c_threads; i++) {
      cf_count += t_c_scan[i*2];
      c_nnz_0 += t_c_scan[i*2+1];
      t_c_scan[i*2] = cf_count;
      t_c_scan[i*2+1] = c_nnz_0;
    }

    /* Now call post-scan */

    #pragma omp parallel shared(c_nnz_0) num_threads(n_c_threads)
    {
      int t_id = omp_get_thread_num();

      if (t_id > 0) {
        cs_lnum_t t_s_id, t_e_id;
        cs_parall_thread_range(c_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);
        t_s_id += 1;  /* Shift by 1 as we are dealing with an index */
        t_e_id += 1;

        cs_lnum_t n_rows_shift = t_c_scan[(t_id-1)*2];
        cs_lnum_t nnz_0_shift = t_c_scan[(t_id-1)*2 + 1];

        for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
          cf_r_idx[ii] += n_rows_shift;
          c_r_idx_0[ii] += nnz_0_shift;
        }
      }
    }

  }

#endif /* defined(HAVE_OPENMP) */

  BFT_FREE(t_c_scan);

  /* Now assign rows */

  cs_lnum_t *cf_r_ids;
  CS_MALLOC_HD(cf_r_ids, cf_r_idx[c_n_rows], cs_lnum_t, alloc_mode);

  #pragma omp parallel num_threads(n_f_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      cs_lnum_t jj = f_c_row[ii];
      if (jj > -1) {
        cs_lnum_t kk = cf_r_idx[jj] + cf_r_shift[jj];
        cf_r_shift[jj] += 1;
        cf_r_ids[kk] = ii;
      }
    }
  }

  BFT_FREE(cf_r_shift);
  BFT_FREE(t_c_scan);

  *c_f_row_index = cf_r_idx;
  *c_f_row_ids = cf_r_ids;

  if (c_row_index_0 != nullptr)
    *c_row_index_0 = c_r_idx_0;
  else
    BFT_FREE(c_r_idx_0);
}

/*----------------------------------------------------------------------------
 * Build a coarse grid msr structure from the previous level
 * with a matrix in MSR format.
 *
 * parameters:
 *   f_n_rows          <-- number of fine rows
 *   c_n_rows          <-- number of coarse rows
 *   alloc_mode        <-- allocation mode
 *   f_row_index       <-- fine matrix row index
 *   f_col_id          <-- fine matrix column ids
 *   f_c_row           <-- fine to coarse row mapping
 *   c_f_row_index     <-- coarse to fine row index
 *   c_f_row_ids       <-- coarse to fine row ids
 *   c_row_index_0     <-- oversized (pre-merging) row index for coarse matrix
 *   c_row_index       --> coarse row index
 *   c_col_ids         --> coarse column ids
 *----------------------------------------------------------------------------*/

static void
_coarse_msr_struct(cs_lnum_t          f_n_rows,
                   cs_lnum_t          c_n_rows,
                   cs_alloc_mode_t    alloc_mode,
                   const cs_lnum_t   *restrict f_row_index,
                   const cs_lnum_t   *restrict f_col_id,
                   const cs_lnum_t   *restrict f_c_row,
                   const cs_lnum_t   *restrict c_f_row_index,
                   const cs_lnum_t   *restrict c_f_row_ids,
                   const cs_lnum_t   *restrict c_row_index_0,
                   cs_lnum_t        **c_row_index,
                   cs_lnum_t        **c_col_ids)
{
  cs_lnum_t *restrict c_row_idx;
  BFT_MALLOC(c_row_idx, c_n_rows+1, cs_lnum_t);
  c_row_idx[0] = 0;

  cs_lnum_t *t_c_scan = nullptr;
  int n_loc_threads = cs_parall_n_threads(f_n_rows, CS_THR_MIN);
  if (n_loc_threads > 1) {
    BFT_MALLOC(t_c_scan, n_loc_threads*2, cs_lnum_t);
    for (int i = 0; i < n_loc_threads; i++) {
      t_c_scan[i*2] = 0;
      t_c_scan[i*2 + 1] = 0;
    }
  }

  /* Allocate working arrays

     c_row_idx_0 and c_col_id_0 are oversized, as they assume
     each coarse row can have as many columns as the sum of the
     column count of the aggregated fine rows.
     This ensures the working
     array for each thread is large enough and data from each thread
     is independent. Data will then need to be combined */

  cs_lnum_t c_nnz = 0;
  const cs_lnum_t c_nnz_0 = c_row_index_0[c_n_rows];

  cs_lnum_t *c_col_id_0;
  BFT_MALLOC(c_col_id_0, c_nnz_0, cs_lnum_t);

  /* Loop on coarse rows */

  #pragma omp parallel shared(c_nnz) num_threads(n_loc_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(c_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    cs_lnum_t *t_col_id_p = c_col_id_0 + c_row_index_0[t_s_id];
    cs_lnum_t t_nnz = 0;

    /* Loop on block's coarse rows */
    for (cs_lnum_t c_idx = t_s_id; c_idx < t_e_id; c_idx++) {
      cs_lnum_t *restrict l_buf = t_col_id_p + t_nnz;
      cs_lnum_t b_e = 0;
      cs_lnum_t r_s_id = c_f_row_index[c_idx];
      cs_lnum_t r_e_id = c_f_row_index[c_idx+1];
      for (cs_lnum_t r_idx = r_s_id; r_idx < r_e_id; r_idx++) {
        cs_lnum_t ii = c_f_row_ids[r_idx];
        cs_lnum_t s_id = f_row_index[ii];
        cs_lnum_t e_id = f_row_index[ii+1];
        for (cs_lnum_t jj = s_id; jj < e_id; jj++) {
          cs_lnum_t kk = f_c_row[f_col_id[jj]];
          if (kk > -1 && kk != c_idx) {
            l_buf[b_e++] = kk;
          }
        }
      }
      cs_sort_lnum(l_buf, b_e);
      cs_lnum_t r_count = 0;
      if (b_e > 0) {
        r_count = 1;
        for (cs_lnum_t jj = 1; jj < b_e; jj++) {
          if (l_buf[jj] > l_buf[jj-1]) {
            l_buf[r_count] = l_buf[jj];
            r_count++;
          }
        }
      }
      t_nnz += r_count;
      c_row_idx[c_idx + 1] = t_nnz; // final value set here
    }

#if defined(HAVE_OPENMP)
    if (t_c_scan != nullptr) {
      int t_id = omp_get_thread_num();
      t_c_scan[t_id*2] = c_row_index_0[t_s_id];
      t_c_scan[t_id*2+1] = t_nnz;
    }
    else
      c_nnz = t_nnz;
#else
    c_nnz = t_nnz;
#endif

  } // End of OpenMP section

  cs_lnum_t *c_col_id = nullptr;

#if defined(HAVE_OPENMP)

  /* Assemble thread submatrices */

  if (n_loc_threads > 1) {

    c_nnz = 0;
    for (int i = 0; i < n_loc_threads; i++) {
      c_nnz += t_c_scan[i*2+1];
      t_c_scan[i*2+1] = c_nnz;
    }

    CS_MALLOC_HD(c_col_id, c_nnz, cs_lnum_t, alloc_mode);

    #pragma omp parallel num_threads(n_loc_threads)
    {
      int t_id = omp_get_thread_num();
      cs_lnum_t t_s_id, t_e_id;
      cs_parall_thread_range(c_n_rows, sizeof(cs_real_t), t_id, n_loc_threads,
                             &t_s_id, &t_e_id);

      cs_lnum_t c_s_id = 0, c_e_id = t_c_scan[t_id*2], t_nnz = t_c_scan[1];
      if (t_id > 0) {
        c_s_id = t_c_scan[(t_id-1)*2+1];
        c_e_id = t_c_scan[(t_id)*2+1];
        t_nnz = c_e_id - c_s_id;
      }
      memcpy(c_col_id + c_s_id,
             c_col_id_0 + c_row_index_0[t_s_id],
             sizeof(cs_lnum_t)*t_nnz);
      if (t_id > 0) {
        cs_lnum_t shift =  t_c_scan[(t_id-1)*2+1];
        for (cs_lnum_t c_idx = t_s_id; c_idx < t_e_id; c_idx++) {
          c_row_idx[c_idx + 1] += shift;
        }
      }
    }

    BFT_FREE(c_col_id_0);

  }

#endif

  if (n_loc_threads == 1) {
    CS_REALLOC_HD(c_col_id_0, c_nnz, cs_lnum_t, alloc_mode);
    c_col_id = c_col_id_0;
  }

  /* Free working arrays */

  BFT_FREE(t_c_scan);

  *c_row_index = c_row_idx;
  *c_col_ids = c_col_id;
}

/*----------------------------------------------------------------------------
 * Compute MSR face adjacency for a give grid.
 *
 * We choose to orient faces in order of increasing local row id.
 *
 * Note that in the case of distributed (MPI) computation, this means
 * that faces on parallel boundaries will not have the same orientation
 * for each rank, which is not an issue here as we should never need this
 * hypothesis in the multigrid infrastructure (contrary to other portions
 * of the code where summing on faces may be done).
 *
 * parameters:
 *   g             <-- pointer to grid structure.
 *   row_index     <-- matrix row index
 *   col_id        <-- matrix column ids
 *----------------------------------------------------------------------------*/

static void
_msr_face_adjacency(cs_grid_t        *g,
                    const cs_lnum_t  *restrict row_index,
                    const cs_lnum_t  *restrict col_id)
{
  std::chrono::high_resolution_clock::time_point tm_start;
  if (cs_glob_timer_kernels_flag > 0)
    tm_start = std::chrono::high_resolution_clock::now();

  const cs_lnum_t n_rows = g->n_rows;
  const cs_lnum_t c_nnz = row_index[g->n_rows];

  cs_lnum_t *cell_to_face;
  short int *cell_to_face_sgn;
  CS_MALLOC_HD(cell_to_face, c_nnz, cs_lnum_t, g->alloc_mode);
  CS_MALLOC_HD(cell_to_face_sgn, c_nnz, short int, g->alloc_mode);

  cs_lnum_t *t_f_scan = nullptr;
  int n_loc_threads = cs_parall_n_threads(c_nnz, CS_THR_MIN);

  cs_lnum_t n_faces = 0;

  if (n_loc_threads > 1) {
    BFT_MALLOC(t_f_scan, n_loc_threads*2, cs_lnum_t);
    for (int i = 0; i < n_loc_threads; i++) {
      t_f_scan[i] = 0;
    }
  }

  /* First stage: determine coarse face ids in direct orientation */

  #pragma omp parallel shared(n_faces) num_threads(n_loc_threads)
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    cs_lnum_t t_n_faces = 0;

    /* Loop on block's coarse rows */
    for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++) {
      const cs_lnum_t s_id = row_index[c_id];
      const cs_lnum_t e_id = row_index[c_id+1];

      for (cs_lnum_t idx = s_id; idx < e_id; idx++) {
        cs_lnum_t c_id_a = col_id[idx];
        if (c_id_a < c_id) {
          cell_to_face[idx] = -1;
          cell_to_face_sgn[idx] = -1;
        }
        else {
          cell_to_face[idx] = t_n_faces;
          cell_to_face_sgn[idx] = 1;
          t_n_faces++;
        }
      }
    }

#if defined(HAVE_OPENMP)
    if (t_f_scan != nullptr) {
      int t_id = omp_get_thread_num();
      t_f_scan[t_id] = t_n_faces;
    }
    else
      n_faces = t_n_faces;
#else
    n_faces = t_n_faces;
#endif
  } // End of OpenMP section

  /* Second stage: scan operation in case of multiple threads */

#if defined(HAVE_OPENMP)

  if (n_loc_threads > 1) {

    n_faces = 0;
    for (int i = 0; i < n_loc_threads; i++) {
      n_faces += t_f_scan[i];
      t_f_scan[i] = n_faces;
    }

    #pragma omp parallel num_threads(n_loc_threads)
    {
      int t_id = omp_get_thread_num();

      if (t_id > 0) {
        cs_lnum_t t_s_id, t_e_id;
        cs_parall_thread_range(n_rows, sizeof(cs_real_t),
                               t_id, n_loc_threads, &t_s_id, &t_e_id);

        const cs_lnum_t face_shift = t_f_scan[t_id-1];

        const cs_lnum_t s_idx = row_index[t_s_id];
        const cs_lnum_t e_idx = row_index[t_e_id];

        for (cs_lnum_t idx = s_idx; idx < e_idx; idx++) {
          if (cell_to_face[idx] > -1)
            cell_to_face[idx] += face_shift;
        }
      }

    }

  }

#endif /* defined(HAVE_OPENMP) */

  BFT_FREE(t_f_scan);

  /* Third stage: determine coarse face ids in opposite orientation */

  #pragma omp parallel for num_threads(n_loc_threads)
  for (cs_lnum_t c_id = 0; c_id < n_rows; c_id++) {
    const cs_lnum_t s_id = row_index[c_id];
    const cs_lnum_t e_id = row_index[c_id+1];
    for (cs_lnum_t idx = s_id; idx < e_id; idx++) {
      if (cell_to_face[idx] == -1) {
        cs_lnum_t c_id_a = col_id[idx];
        const cs_lnum_t s_id_a = row_index[c_id_a];
        const cs_lnum_t e_id_a = row_index[c_id_a];
        const cs_lnum_t idx_a = _l_id_binary_search(e_id_a - s_id_a,
                                                    c_id,
                                                    col_id + s_id_a);
        cell_to_face[idx] = cell_to_face[s_id_a + idx_a];
      }
    }
  }

  g->n_faces = n_faces;
  g->cell_face = cell_to_face;
  g->cell_face_sgn = cell_to_face_sgn;

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

/*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using
 * an automatic criterion, using the face to cells adjacency.
 *
 * parameters:
 *   f                    <-- Fine grid structure
 *   coarsening_type      <-- Coarsening type
 *   max_aggregation      <-- Max fine cells per coarse cell
 *   relaxation_parameter <-- P0/P1 relaxation factor
 *   verbosity            <-- Verbosity level
 *   f_c_cell             --> Fine cell -> coarse cell connectivity
 *
 * return:
 *   number of coarse rows
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_automatic_aggregation_fc(const cs_grid_t       *f,
                          cs_grid_coarsening_t   coarsening_type,
                          cs_lnum_t              max_aggregation,
                          double                 relaxation_parameter,
                          int                    verbosity,
                          cs_lnum_t             *f_c_cell)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_lnum_t n_faces;

  cs_lnum_t isym = 2;
  int ncoarse = 8, npass_max = 10, inc_nei = 0;
  int _max_aggregation = 1, npass = 0;

  cs_lnum_t f_n_cells = f->n_rows;
  cs_lnum_t f_n_cells_ext = f->n_cols_ext;
  cs_lnum_t f_n_faces = f->n_faces;

  cs_lnum_t aggr_count = f_n_cells;

  cs_lnum_t r_n_faces = f_n_faces;
  cs_lnum_t c_n_cells = 0;

  cs_real_t epsilon = 1.e-6;

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t eb_size = f->eb_size;
  const cs_lnum_2_t *f_face_cells = f->face_cell;
  const cs_real_t *f_da = f->da;
  const cs_real_t *f_xa = f->xa;

  /* Reduce block to scalar equivalents */

  const cs_real_t *_f_da = f_da;
  const cs_real_t *_f_xa = f_xa;

  cs_real_t *s_da = NULL, *s_xa = NULL;

  if (db_size > 1) {
    BFT_MALLOC(s_da, f_n_cells, cs_real_t);
    _reduce_block(f_n_cells, db_size, f_da, s_da);
    _f_da = s_da;
  }

  if (eb_size > 1) {
    BFT_MALLOC(s_xa, f_n_faces*isym, cs_real_t);
    _reduce_block(f_n_faces*isym, eb_size, f_xa, s_xa);
    _f_xa = s_xa;
  }

  /* Allocate working arrays */

  cs_lnum_t *i_work_array = NULL;
  BFT_MALLOC(i_work_array, f_n_cells_ext*2 + f_n_faces*3, cs_lnum_t);

  cs_lnum_t *c_cardinality = i_work_array;
  cs_lnum_t *c_aggr_count = i_work_array + f_n_cells_ext;
  cs_lnum_t *f_c_face = i_work_array + 2*f_n_cells_ext; /* fine -> coarse face */
  cs_lnum_t *merge_flag = i_work_array + 2*f_n_cells_ext + f_n_faces;

  /* Initialization */

  if (f->symmetric == true)
    isym = 1;

# pragma omp parallel for if(f_n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < f_n_cells_ext; ii++) {
    c_cardinality[ii] = -1;
    c_aggr_count[ii] = 1;
  }

  cs_real_t *xv_sum = nullptr;
  if (f->level == 0) {
    BFT_MALLOC(xv_sum, f_n_cells_ext, cs_real_t);
#   pragma omp parallel for if(f_n_cells_ext > CS_THR_MIN)
    for (cs_lnum_t ii = 0; ii < f_n_cells_ext; ii++) {
      xv_sum[ii] = 0;
    }
  }

# pragma omp parallel for if(f_n_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < f_n_faces; face_id++) {
    merge_flag[face_id] = face_id +1;
    f_c_face[face_id] = 0;
  }

  /* Compute cardinality (number of neighbors for each cell -1) */

  if (f->level == 0) {
    for (cs_lnum_t face_id = 0; face_id < f_n_faces; face_id++) {
      cs_lnum_t ii = f_face_cells[face_id][0];
      cs_lnum_t jj = f_face_cells[face_id][1];

      c_cardinality[ii] += 1;
      c_cardinality[jj] += 1;

      xv_sum[ii] += fabs(_f_xa[face_id*isym]);
      xv_sum[jj] += fabs(_f_xa[face_id*isym + isym-1]);
    }

    /* Also check for penalization */

    for (cs_lnum_t ii = 0; ii < f_n_cells; ii++) {
      if (_f_da[ii] > _penalization_threshold * xv_sum[ii])
        f_c_cell[ii] = -1;
    }

    BFT_FREE(xv_sum);
  }

  else { /* if f->level > 0) */
    for (cs_lnum_t face_id = 0; face_id < f_n_faces; face_id++) {
      cs_lnum_t ii = f_face_cells[face_id][0];
      cs_lnum_t jj = f_face_cells[face_id][1];

      c_cardinality[ii] += 1;
      c_cardinality[jj] += 1;
    }
  }

  /* Passes */

  if (verbosity > 3)
    bft_printf("\n     %s:\n", __func__);

  do {

    inc_nei++;
    npass++;
    n_faces = r_n_faces;
    _max_aggregation++;
    _max_aggregation = CS_MIN(_max_aggregation, max_aggregation);

#   pragma omp parallel for if(n_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
      f_c_face[face_id] = merge_flag[face_id];
      merge_flag[face_id] = 0;
    }

    if (verbosity > 3)
      bft_printf("       pass %3d; r_n_faces = %10ld; aggr_count = %10ld\n",
                 npass, (long)r_n_faces, (long)aggr_count);

    /* Initialize non-eliminated faces */
    r_n_faces = 0;

    /* Loop on non-eliminated faces */

    cs_real_t ag_threshold = 1. - epsilon;
    if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX)
      ag_threshold = (1. - epsilon) * pow(relaxation_parameter, npass);

    for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {

      cs_lnum_t c_face = f_c_face[face_id] -1;

      cs_lnum_t ii = f_face_cells[c_face][0];
      cs_lnum_t jj = f_face_cells[c_face][1];

      /* Exclude faces on parallel or periodic boundary, so as not to */
      /* coarsen the grid across those boundaries (which would change */
      /* the communication pattern and require a more complex algorithm). */

      if (   (ii < f_n_cells && jj < f_n_cells)
          && (f_c_cell[ii] < -1 || f_c_cell[jj] < -1)) {

        cs_lnum_t count = 0;

        cs_lnum_t ix0 = c_face*isym, ix1 = (c_face +1)*isym -1;

        cs_real_t f_da0_da1_inv
          =   ((c_cardinality[ii]+inc_nei)*(c_cardinality[jj]+inc_nei))
            / (_f_da[ii] * _f_da[jj]);

        cs_real_t aggr_crit;

        if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX) {
          cs_real_t f_xa0 = cs_math_fmax(-_f_xa[ix0], 0.);
          cs_real_t f_xa1 = cs_math_fmax(-_f_xa[ix1], 0.);
          aggr_crit = cs_math_fmax(f_xa0, f_xa1) * sqrt(f_da0_da1_inv);
        }
        else {
          cs_real_t f_xa0_xa1 =  _f_xa[ix0] * _f_xa[ix1];
          aggr_crit = f_da0_da1_inv * f_xa0_xa1;
        }

        if (aggr_crit > ag_threshold) {

          if (f_c_cell[ii] > -1 && f_c_cell[jj] < -1 ) {
            if (c_aggr_count[f_c_cell[ii]] < _max_aggregation +1) {
              f_c_cell[jj] = f_c_cell[ii];
              c_aggr_count[f_c_cell[ii]] += 1;
              count++;
            }
          }
          else if (f_c_cell[ii] < -1 && f_c_cell[jj] > -1) {
            if (c_aggr_count[f_c_cell[jj]] < _max_aggregation +1) {
              f_c_cell[ii] = f_c_cell[jj];
              c_aggr_count[f_c_cell[jj]] += 1;
              count++;
            }
          }
          else if (f_c_cell[ii] < -1 && f_c_cell[jj] < -1) {
            f_c_cell[ii] = c_n_cells;
            f_c_cell[jj] = c_n_cells;
            c_aggr_count[c_n_cells] += 1;
            c_n_cells++;
            count++;
          }
        }

        if (count == 0 && (f_c_cell[ii] < -1 || f_c_cell[jj] < -1)) {
          merge_flag[r_n_faces] = c_face +1;
          r_n_faces++;
        }

      }

    }

    /* Check the number of coarse cells created */
    aggr_count = 0;
#   pragma omp parallel for reduction(+:aggr_count) if(f_n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < f_n_cells; i++) {
      if (f_c_cell[i] < -1)
        aggr_count++;
    }

    /* Additional passes if aggregation is insufficient */

    if (   aggr_count == 0 || (c_n_cells + aggr_count)*ncoarse < f_n_cells
        || r_n_faces == 0)
      npass_max = npass;

  } while (npass < npass_max); /* Loop on passes */

  BFT_FREE(s_da);
  BFT_FREE(s_xa);

  /* Finish assembly */
  for (cs_lnum_t i = 0; i < f_n_cells; i++) {
    if (f_c_cell[i] < -1) {
      f_c_cell[i] = c_n_cells;
      c_n_cells++;
    }
  }

  /* Free working arrays */

  BFT_FREE(i_work_array);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, f->level-1);
    printf(", total = %ld\n", elapsed.count());
  }

  return c_n_cells;
}

/*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using
 * an automatic criterion, using a MSR matrix adjacency.
 *
 * parameters:
 *   f                    <-- Fine grid structure
 *   coarsening_type      <-- Coarsening type
 *   max_aggregation      <-- Max fine cells per coarse cell
 *   relaxation_parameter <-- P0/P1 relaxation factor
 *   verbosity            <-- Verbosity level
 *   f_c_row              --> Fine row -> coarse row mapping
 *
 * return:
 *   number of coarse rows
 *----------------------------------------------------------------------------*/

static cs_lnum_t
_automatic_aggregation_dx_msr(const cs_grid_t       *f,
                              cs_grid_coarsening_t   coarsening_type,
                              cs_lnum_t              max_aggregation,
                              double                 relaxation_parameter,
                              int                    verbosity,
                              cs_lnum_t             *f_c_row)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  const int npass_max = 10;

  cs_lnum_t isym = 2;
  if (f->symmetric == true)
    isym = 1;

  const cs_lnum_t f_n_rows = f->n_rows;
  const cs_lnum_t f_n_faces = f->n_faces;

  cs_lnum_t c_n_rows = -1;

  cs_real_t epsilon = 1.e-6;

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t eb_size = f->eb_size;

  /* Access matrix MSR vectors */

  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;
  cs_real_t *_d_val = NULL, *_x_val = NULL;

  cs_matrix_get_msr_arrays(f->matrix,
                           &row_index,
                           &col_id,
                           &d_val,
                           &x_val);

  const cs_lnum_t f_nnz = row_index[f_n_rows];

  if (db_size > 1) {
    BFT_MALLOC(_d_val, f_n_rows, cs_real_t);
    _reduce_block(f_n_rows, db_size, d_val, _d_val);
    d_val = _d_val;
  }

  if (eb_size > 1) {
    BFT_MALLOC(_x_val, f_nnz, cs_real_t);
    _reduce_block(f_nnz, eb_size, x_val, _x_val);
    x_val = _x_val;
  }

  /* Coarse row count scan array for shared-memory parallelism */

  cs_lnum_t *t_c_scan = nullptr;
  int n_loc_threads = cs_parall_n_threads(f_n_rows, CS_THR_MIN);
  if (n_loc_threads > 1) {
    BFT_MALLOC(t_c_scan, n_loc_threads, cs_lnum_t);
    for (int i = 0; i < n_loc_threads; i++)
      t_c_scan[i] = 0;
  }

  /* Allocate working arrays */

  short *c_cardinality, *c_aggr_count;
  BFT_MALLOC(c_cardinality, f_n_rows, short);
  BFT_MALLOC(c_aggr_count, f_n_rows, short);

  /* aggregation queue: local column id, index in matrix, and
     index of symmetric element if needed */
  const cs_lnum_t ag_queue_stride = (1+isym);
  cs_lnum_t ag_work_size = CS_MAX(f_nnz*ag_queue_stride, f_n_rows);
  cs_lnum_t *ag_work;
  BFT_MALLOC(ag_work, ag_work_size, cs_lnum_t);
  cs_lnum_t *ag_queue = ag_work;

  cs_lnum_t *log_counts = nullptr;
  if (verbosity > 3) {
    BFT_MALLOC(log_counts, n_loc_threads*npass_max*2, cs_lnum_t);
    memset(log_counts, 0, n_loc_threads*npass_max*2*sizeof(cs_lnum_t));
  }

  /* Handle a block of rows per thread; aggregation will not be done
     across thread blocks, in a similar manner that it is not done
     across MPI ranks. */

  #pragma omp parallel shared(c_n_rows) num_threads(n_loc_threads)
  {

#if defined(HAVE_OPENMP)
    const int t_id = omp_get_thread_num();
#else
    const int t_id = 0;
#endif

    /* Use a simple thread range for now. Perhaps a more complex range
       computation could lead to better work balance between threads,
       as rows with low ids will tend to get more work than rows
       with high ids, given that we consider the upper triangular part */

    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(f_n_rows, sizeof(cs_real_t), &t_s_id, &t_e_id);

    int ncoarse = 8, inc_nei = 0;
    int _max_aggregation = 1;

    /* Computation of the cardinality and test if the line is
     * penalized, in which case its cardinality is set to 0. */

    if (f->level == 0) {
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        cs_lnum_t t_r_s_id = row_index[ii];
        cs_lnum_t t_r_e_id = row_index[ii+1];

        c_cardinality[ii] = t_r_e_id - t_r_s_id + 1;
        c_aggr_count[ii] = 1;

        cs_real_t  sum = 0.0;
        for (cs_lnum_t jj = t_r_s_id; jj < t_r_e_id; jj++) {
          const cs_real_t  xv = x_val[jj];
          if (xv < 0.)
            sum -= xv;
          else
            sum += xv;
        }

        /* Check if the line seems penalized or not. */
        if (d_val[ii] > _penalization_threshold * sum) {
          f_c_row[ii] = -1;
          c_aggr_count[ii] = 0;
        }
      }
    }
    else {
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        c_cardinality[ii] = row_index[ii+1] - row_index[ii] + 1;
        c_aggr_count[ii] = 1;
      }
    }

    /* Indexes and counterss */

    const cs_lnum_t t_r_s_id = row_index[t_s_id] * ag_queue_stride;
    cs_lnum_t t_r_e_id = t_r_s_id;

    cs_lnum_t t_c_n_rows = 0;
    cs_lnum_t r_n_faces = 0;

    /* Passes */

    if (verbosity > 3 && t_s_id == 0)
      bft_printf("\n     %s:\n", __func__);

    /* First pass */

    {
      int npass = 1;

      if (_max_aggregation < max_aggregation)
        _max_aggregation++;

      if (verbosity > 3 && t_id == 0) {
        bft_printf("       pass 1; r_n_faces = %ld;"
                   " aggr_count = %ld\n",
                   (long)f_n_faces, (long)f_n_rows);
      }

      /* Loop on non-eliminated faces */

      cs_real_t ag_threshold = 1. - epsilon;
      if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX)
        ag_threshold = (1. - epsilon) * pow(relaxation_parameter, npass);

      /* For each thread, build indexes for subsequent
         passes; these are based on the matrix adjacency (oversized
         since only positively oriented, and we only need to hande
         faces not merged on the first pass, but avoids reallocation
         or extra counting passes) */

      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        cs_lnum_t row_s_id = row_index[ii];
        cs_lnum_t row_e_id = row_index[ii+1];
        for (cs_lnum_t ix0 = row_s_id; ix0 < row_e_id; ix0++) {

          cs_lnum_t jj = col_id[ix0];

          /* Only consider faces in upper triangular part, as
             each face (graph edge) should appear once each for
             triangular part in the matrix.

             Also exclude aggregation on parallel or periodic boundary, so as
             not to coarsen the grid across those boundaries (which would
             change the communication pattern and require a more complex
             algorithm), and even across "artificial" local thread boundaries,
             for simpler multithreading.

             So we ignore this entry if jj < ii, or if jj is outside
             the [t_s_id, t_e_id[ bounds. Since ii >= t_s_id, this
             simplifies to the test below. */

          if (jj < ii || jj >= t_e_id)
            continue;

          cs_lnum_t ii_c = f_c_row[ii];
          cs_lnum_t jj_c = f_c_row[jj];

          /* If both rows are already aggregated or one is penalized, pass */

          if (ii_c == -1 || jj_c == -1 || (ii_c >= 0 && jj_c >= 0))
            continue;

          cs_lnum_t ix1 = ix0;
          if (isym == 2) {
            ix1 = row_index[jj];
            cs_lnum_t ix1_e = row_index[jj+1];
            while (col_id[ix1] != ii && ix1 < ix1_e)
              ix1++;
            if (ix1 >= ix1_e) // non-symmetric structure; pass (should not occur)
              continue;
          }

          cs_real_t f_da0_da1_inv =   (c_cardinality[ii]*c_cardinality[jj])
                                    / (d_val[ii] * d_val[jj]);

          cs_real_t aggr_crit;

          if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX) {
            cs_real_t f_xa0 = cs_math_fmax(-x_val[ix0], 0.);
            cs_real_t f_xa1 = cs_math_fmax(-x_val[ix1], 0.);
            aggr_crit = cs_math_fmax(f_xa0, f_xa1) * sqrt(f_da0_da1_inv);
          }
          else {
            cs_real_t f_xa0_xa1 =  x_val[ix0] * x_val[ix1];
            aggr_crit = f_da0_da1_inv * f_xa0_xa1;
          }

          int count = 0;

          if (aggr_crit > ag_threshold) {
            if (ii_c > -1 && jj_c < -1 ) {
              if (c_aggr_count[t_s_id + ii_c] < _max_aggregation +1) {
                f_c_row[jj] = ii_c;
                c_aggr_count[t_s_id + ii_c] += 1;
                count++;
              }
            }
            else if (ii_c < -1 && jj_c > -1) {
              if (c_aggr_count[t_s_id + jj_c] < _max_aggregation +1) {
                f_c_row[ii] = jj_c;
                c_aggr_count[t_s_id + jj_c] += 1;
                count++;
              }
            }
            else if (ii_c < -1 && jj_c < -1) {
              f_c_row[ii] = t_c_n_rows;
              f_c_row[jj] = t_c_n_rows;
              c_aggr_count[t_s_id + t_c_n_rows] += 1;
              t_c_n_rows++;
              count++;
            }
          }

          if (count == 0) {
            ag_queue[t_r_e_id] = ii;
            ag_queue[t_r_e_id+1] = ix0;
            if (ag_queue_stride == 3)
              ag_queue[t_r_e_id+1] = ix1;
            t_r_e_id += ag_queue_stride;
            r_n_faces++;
          }

        }  /* Loop on columns */

      }  /* Loop on rows */
    }

    /* Other passes */

    for (int npass = 2; npass <= npass_max; npass++) {

      /* Exit loop on passes if aggregation is sufficient
         (note that each thread may loop independently). */

      cs_lnum_t aggr_count = 0; /* Count of cells not yet aggregated */
      for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
        if (f_c_row[ii] < -1)
          aggr_count++;
      }

      if (   aggr_count == 0
          || (t_c_n_rows + aggr_count)*ncoarse < (t_e_id - t_s_id)
          || r_n_faces == 0) {
        break;
      }

      if (_max_aggregation < max_aggregation)
        _max_aggregation++;

      if (verbosity > 3) {
        log_counts[t_id*(npass-1)*2] = r_n_faces;
        log_counts[t_id*(npass-1)*2 + 1] = aggr_count;
      }

      /* Re-initialize non-eliminated faces */
      r_n_faces = 0;

      inc_nei += 1;

      /* Loop on non-eliminated faces (i.e. "aggregation queue").

         Faces adjacent to penalized rows, or row on parallel boundaries,
         were already excluded from the aggregation queue in the first pass,
         so we do not need to test for those anymore. */

      cs_lnum_t r_idx_n = t_r_s_id;

      for (cs_lnum_t r_idx = t_r_s_id;
           r_idx < t_r_e_id;
           r_idx += ag_queue_stride) {

        cs_real_t ag_threshold = 1. - epsilon;
        if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX)
          ag_threshold = (1. - epsilon) * pow(relaxation_parameter, npass);

        cs_lnum_t ii = ag_queue[r_idx];
        cs_lnum_t ix0 = ag_queue[r_idx+1];
        cs_lnum_t ix1 = (ag_queue_stride == 3) ? ag_queue[r_idx+2] : ix0;

        cs_lnum_t jj = col_id[ix0];

        /* If both cells are already aggregated, pass */

        cs_lnum_t ii_c = f_c_row[ii];
        cs_lnum_t jj_c = f_c_row[jj];

        if (ii_c >= 0 && jj_c >= 0)
          continue;

        cs_lnum_t count = 0;

        cs_real_t f_da0_da1_inv = (   (c_cardinality[ii]+inc_nei)
                                    * (c_cardinality[jj]+inc_nei))
                                  / (d_val[ii] * d_val[jj]);

        cs_real_t aggr_crit;

        if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX) {
          cs_real_t f_xa0 = cs_math_fmax(-x_val[ix0], 0.);
          cs_real_t f_xa1 = cs_math_fmax(-x_val[ix1], 0.);
          aggr_crit = cs_math_fmax(f_xa0, f_xa1) * sqrt(f_da0_da1_inv);
        }
        else {
          cs_real_t f_xa0_xa1 =  x_val[ix0] * x_val[ix1];
          aggr_crit = f_da0_da1_inv * f_xa0_xa1;
        }

        if (aggr_crit > ag_threshold) {
          if (ii_c > -1 && jj_c < -1 ) {
            if (c_aggr_count[t_s_id + ii_c] < _max_aggregation +1) {
              f_c_row[jj] = ii_c;
              c_aggr_count[t_s_id + ii_c] += 1;
              count++;
            }
          }
          else if (ii_c < -1 && jj_c > -1) {
            if (c_aggr_count[t_s_id + jj_c] < _max_aggregation +1) {
              f_c_row[ii] = jj_c;
              c_aggr_count[t_s_id + jj_c] += 1;
              count++;
            }
          }
          else if (ii_c < -1 && jj_c < -1) {
            f_c_row[ii] = t_c_n_rows;
            f_c_row[jj] = t_c_n_rows;
            c_aggr_count[t_s_id + t_c_n_rows] += 1;
            t_c_n_rows++;
            count++;
          }
        }

        /* Place faces where no merging occurs back in queue
           (we reuse the same location, knowing that
           t_r_s_id + n_r_faces*ag_queue_stride < r_id) */
        if (count == 0) {
          assert(r_idx_n < t_r_e_id);
          ag_queue[r_idx_n] = ii;
          ag_queue[r_idx_n+1] = ix0;
          if (ag_queue_stride == 3)
            ag_queue[r_idx_n+2] = ix1;
          r_idx_n += ag_queue_stride;
          r_n_faces++;
        }

      } /* Loop on aggregation queue */

      t_r_e_id = r_idx_n;

    } /* Loop on passes */

    /* Finish assembly */
    for (cs_lnum_t ii = t_s_id; ii < t_e_id; ii++) {
      if (f_c_row[ii] < -1) {
        f_c_row[ii] = t_c_n_rows;
        t_c_n_rows++;
      }
    }

    /* We need to wait for all threads if we want to reuse
       the work array, as the portion we use is not the same
       (we could use the same portion if we were absolutely sure that
        (row_index[t_s_id] - row_index[t_e_id]) / 2 > (t_e_id - t_s_id),
       which is always the case if there is at least one non-diagonal
       term per row, but could be false in some strage configurations). */

    #pragma omp barrier

    _reorder_f_c_row(t_e_id - t_s_id,
                     t_c_n_rows,
                     ag_work + t_s_id,
                     f_c_row + t_s_id);

    /* Prepare combining thread results */

    if (n_loc_threads > 1)
      t_c_scan[t_id] = t_c_n_rows;
    else
      c_n_rows = t_c_n_rows;

  } // End of OpenMP section

  if (verbosity > 3) {
    for (int npass = 2; npass <= npass_max; npass++) {
      for (int t_id = 0; t_id < n_loc_threads; t_id++) {
        cs_lnum_t r_n_faces = log_counts[t_id*(npass-1)*2];
        cs_lnum_t aggr_count = log_counts[t_id*(npass-1)*2 + 1];
        if (r_n_faces > 0 && aggr_count > 0)
          bft_printf("       pass %d, thread %d; r_n_faces = %ld;"
                     " aggr_count = %ld\n",
                     npass, t_id, (long)r_n_faces, (long)aggr_count);
      }
    }
    BFT_FREE(log_counts);
  }

  /* Combine thread results */

#if defined(HAVE_OPENMP)
  if (n_loc_threads > 1) {
    c_n_rows = _scan_f_c_row(n_loc_threads, f_n_rows, t_c_scan, f_c_row);
  }
#endif

  /* Free working arrays */

  BFT_FREE(ag_work);
  BFT_FREE(c_cardinality);
  BFT_FREE(c_aggr_count);

  BFT_FREE(_d_val);
  BFT_FREE(_x_val);

  BFT_FREE(t_c_scan);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, f->level-1);
    printf(", total = %ld\n", elapsed.count());
  }

  return c_n_rows;
}

/*----------------------------------------------------------------------------
 * Compute volume and center of coarse cells.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-> Coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_cell_quantities(const cs_grid_t  *fine_grid,
                                cs_grid_t        *coarse_grid)
{
  cs_lnum_t f_n_cells = fine_grid->n_rows;

  cs_lnum_t c_n_cells = coarse_grid->n_rows;
  cs_lnum_t c_n_cells_ext = coarse_grid->n_cols_ext;

  cs_lnum_t *c_coarse_row = coarse_grid->coarse_row;

  cs_real_t *c_cell_vol = coarse_grid->_cell_vol;
  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;

  const cs_real_t *f_cell_vol = fine_grid->cell_vol;
  const cs_real_t *f_cell_cen = fine_grid->cell_cen;

  /* Compute volume and center of coarse cells */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t ic = 0; ic < c_n_cells_ext; ic++) {
    c_cell_vol[ic] = 0.;
    c_cell_cen[3*ic]    = 0.;
    c_cell_cen[3*ic +1] = 0.;
    c_cell_cen[3*ic +2] = 0.;
  }

  for (cs_lnum_t ii = 0; ii < f_n_cells; ii++) {
    cs_lnum_t ic = c_coarse_row[ii];
    if (ic > -1) {
      c_cell_vol[ic] += f_cell_vol[ii];
      c_cell_cen[3*ic]    += f_cell_vol[ii]*f_cell_cen[3*ii];
      c_cell_cen[3*ic +1] += f_cell_vol[ii]*f_cell_cen[3*ii +1];
      c_cell_cen[3*ic +2] += f_cell_vol[ii]*f_cell_cen[3*ii +2];
    }
  }

# pragma omp parallel for if(c_n_cells > CS_THR_MIN)
  for (cs_lnum_t ic = 0; ic < c_n_cells; ic++) {
    c_cell_cen[3*ic]    /= c_cell_vol[ic];
    c_cell_cen[3*ic +1] /= c_cell_vol[ic];
    c_cell_cen[3*ic +2] /= c_cell_vol[ic];
  }
}

/*----------------------------------------------------------------------------
 * Verification for matrix quantities.
 *
 * parameters:
 *   g <-- pointer to grid structure
 *----------------------------------------------------------------------------*/

static void
_verify_matrix(const cs_grid_t  *g)
{
  const cs_lnum_t n_cols = cs_matrix_get_n_columns(g->matrix);
  const cs_lnum_t n_rows = cs_matrix_get_n_rows(g->matrix);
  const cs_lnum_t db_size = g->db_size;

  cs_real_t *val;
  BFT_MALLOC(val, n_cols*db_size, cs_real_t);

  /* Evaluate fine and coarse grids diagonal dominance */

  cs_matrix_diag_dominance(g->matrix, val);

  cs_real_t vmin = HUGE_VAL, vmax = -HUGE_VAL;

  const cs_lnum_t n_ents = n_rows*db_size;

  for (cs_lnum_t i = 0; i < n_ents; i++) {
    if (val[i] < vmin)
      vmin = val[i];
    else if (val[i] > vmax)
      vmax = val[i];
  }

  BFT_FREE(val);

#if defined(HAVE_MPI)
  if (cs_glob_mpi_comm != MPI_COMM_NULL) {
    cs_real_t _vmin = vmin, _vmax = vmax;
    MPI_Allreduce(&_vmin, &vmin, 1, CS_MPI_REAL, MPI_MIN, g->comm);
    MPI_Allreduce(&_vmax, &vmax, 1, CS_MPI_REAL, MPI_MAX, g->comm);
  }
#endif

  bft_printf(_("       grid level %2d diag. dominance: min = %12.5e\n"
               "                                      max = %12.5e\n\n"),
             g->level, vmin, vmax);
}

/*----------------------------------------------------------------------------
 * Verification for coarse quantities computed from fine quantities.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-- Coarse grid structure
 *   n_clips_min <-- number of clippings to minimum value
 *   n_clips_max <-- number of clippings to maximum value
 *   interp      <-- 0 for no intepolation, > 0 for interpolation
 *----------------------------------------------------------------------------*/

static void
_verify_coarse_quantities(const cs_grid_t  *fine_grid,
                          const cs_grid_t  *coarse_grid,
                          cs_gnum_t         n_clips_min,
                          cs_gnum_t         n_clips_max,
                          int               interp)
{
  cs_lnum_t ic, jc, ii, jj, c_face, face_id;

  int isym = 2;

  cs_lnum_t f_n_cells = fine_grid->n_rows;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cols_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells = coarse_grid->n_rows;
  cs_lnum_t c_n_cells_ext = coarse_grid->n_cols_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  const cs_real_t *c_xa0 = coarse_grid->_xa0;
  const cs_real_t *c_xa = coarse_grid->_xa;

  cs_real_t *w1 = NULL;

  const cs_lnum_t db_size = fine_grid->db_size;
  const cs_lnum_t db_stride = db_size*db_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_xa = fine_grid->xa;

  BFT_MALLOC(w1, f_n_cells_ext*db_stride, cs_real_t);

  if (fine_grid->symmetric == true)
    isym = 1;

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
  MPI_Comm comm = fine_grid->comm;
  if (comm != MPI_COMM_NULL) {
    cs_gnum_t n_clips[2] = {n_clips_min, n_clips_max};
    MPI_Allreduce(MPI_IN_PLACE, n_clips, 2, CS_MPI_GNUM, MPI_SUM, comm);
    n_clips_min = n_clips[0];
    n_clips_max = n_clips[0];
  }
#endif

  if (n_clips_min+n_clips_max > 0)
    bft_printf("\n     %s:\n"
               "       coarse_matrix < xag0 on %10llu faces\n"
               "                     > 0    on %10llu faces\n",
               __func__,
               (unsigned long long)n_clips_min,
               (unsigned long long)n_clips_max);

  cs_real_t *w2, *w3, *w4;
  double anmin[2] = {HUGE_VAL, HUGE_VAL};
  double anmax[2] = {-HUGE_VAL, -HUGE_VAL};

  BFT_MALLOC(w2, f_n_cells_ext*db_stride, cs_real_t);
  BFT_MALLOC(w3, c_n_cells_ext*db_stride, cs_real_t);
  BFT_MALLOC(w4, c_n_cells_ext*db_stride, cs_real_t);

  /* Evaluate anisotropy of fine and coarse grids */

  for (ii = 0; ii < f_n_cells_ext; ii++) {
    w1[ii] = -HUGE_VAL;
    w2[ii] = HUGE_VAL;
  }

  for (ic = 0; ic < c_n_cells_ext; ic++) {
    w3[ic] = -HUGE_VAL;
    w4[ic] = HUGE_VAL;
  }

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cell[face_id][0];
    jj = f_face_cell[face_id][1];
    w1[ii] = cs_math_fmax(fabs(f_xa[face_id*isym]), w1[ii]);
    w2[ii] = cs_math_fmin(fabs(f_xa[face_id*isym]), w2[ii]);
    w1[jj] = cs_math_fmax(fabs(f_xa[(face_id +1)*isym -1]), w1[jj]);
    w2[jj] = cs_math_fmin(fabs(f_xa[(face_id +1)*isym -1]), w2[jj]);
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    w3[ic] = cs_math_fmax(fabs(c_xa[c_face*isym]), w3[ic]);
    w4[ic] = cs_math_fmin(fabs(c_xa[c_face*isym]), w4[ic]);
    w3[jc] = cs_math_fmax(fabs(c_xa[(c_face +1)*isym -1]), w3[jc]);
    w4[jc] = cs_math_fmin(fabs(c_xa[(c_face +1)*isym -1]), w4[jc]);
  }

  for (ii = 0; ii < f_n_cells; ii++)
    w1[ii] = w2[ii] / w1[ii];

  for (ic = 0; ic < c_n_cells; ic++)
    w3[ic] = w4[ic] / w3[ic];

  anmin[0] = HUGE_VAL; anmin[1] = HUGE_VAL;
  anmax[0] = -HUGE_VAL; anmax[1] = -HUGE_VAL;

  for (ii = 0; ii < f_n_cells; ii++) {
    if (w1[ii] < anmin[0])
      anmin[0] = w1[ii];
    else if (w1[ii] > anmax[0])
      anmax[0] = w1[ii];
  }

  for (ic = 0; ic < c_n_cells; ic++) {
    if (w3[ic] < anmin[1])
      anmin[1] = w3[ic];
    else if (w3[ic] > anmax[1])
      anmax[1] = w3[ic];
  }

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
  if (comm != MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE, anmin, 2, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, anmax, 2, MPI_DOUBLE, MPI_MAX, comm);
  }
#endif

  bft_printf(_("       fine   grid anisotropy: min      = %12.5e\n"
               "                               max      = %12.5e\n"
               "       coarse grid anisotropy: min      = %12.5e\n"
               "                               max      = %12.5e\n"),
             anmin[0], anmax[0], anmin[1], anmax[1]);

  BFT_FREE(w2);
  BFT_FREE(w4);

  if (interp == 1) {
    double rmin = HUGE_VAL, rmax = -HUGE_VAL;
    for (c_face = 0; c_face < c_n_faces; c_face++) {
      rmin = cs_math_fmin(rmin, c_xa[c_face*isym] / c_xa0[c_face]);
      rmax = cs_math_fmax(rmax, c_xa[c_face*isym] / c_xa0[c_face]);
    }
#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
    if (comm != MPI_COMM_NULL) {
      MPI_Allreduce(MPI_IN_PLACE, &rmin, 1, MPI_DOUBLE, MPI_MIN, comm);
      MPI_Allreduce(MPI_IN_PLACE, &rmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    }
#endif
    bft_printf(_("       minimum xag_p1 / xag_p0          = %12.5e\n"
                 "       maximum xag_p1 / xag_p0          = %12.5e\n"),
               rmin, rmax);
  }

  BFT_FREE(w3);
  BFT_FREE(w1);
}

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
                                     c->n_rows,
                                     c->n_cols_ext,
                                     &c_row_index,
                                     &c_col_id,
                                     c->halo,
                                     NULL);

  c->matrix_struct = ms;

  c->_matrix = cs_matrix_create(c->matrix_struct);
  c->matrix = c->_matrix;

  const cs_lnum_t *_c_row_index, *_c_col_id;

  cs_matrix_get_msr_arrays(c->matrix,
                           &_c_row_index, &_c_col_id,
                           NULL, NULL);

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
 * Build a coarse level from a finer level using face->cells connectivity
 *
 * parameters:
 *   fine_grid   <-- fine grid structure
 *   coarse_grid <-> coarse grid structure
 *   verbosity   <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities_native(const cs_grid_t  *fine_grid,
                                  cs_grid_t        *coarse_grid,
                                  int               verbosity)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_lnum_t ii, jj, kk, face_id;

  cs_real_t dsigjg, dsxaij, agij;

  int isym = 2;

  cs_lnum_t f_n_cells = fine_grid->n_rows;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cols_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells_ext = coarse_grid->n_cols_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  cs_lnum_t *c_coarse_row = coarse_grid->coarse_row;
  cs_lnum_t *c_coarse_face = coarse_grid->coarse_face;

  cs_gnum_t n_clips_min = 0, n_clips_max = 0;

  cs_real_t *f_xa0ij = fine_grid->xa0ij;

  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;
  cs_real_t *c_face_normal = coarse_grid->_face_normal;

  cs_real_t *c_xa0 = coarse_grid->_xa0;
  cs_real_t *c_xa0ij = coarse_grid->xa0ij;
  cs_real_t *c_da = coarse_grid->_da;
  cs_real_t *c_xa = coarse_grid->_xa;

  const cs_lnum_t db_size = fine_grid->db_size;
  const cs_lnum_t db_stride = db_size*db_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_face_normal = fine_grid->face_normal;
  const cs_real_t *f_xa0 = fine_grid->xa0;
  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_xa = fine_grid->xa;

  const cs_real_t relax_param = coarse_grid->relaxation;

  if (fine_grid->symmetric == true)
    isym = 1;

  /*  Finalize computation of matrix in c_da, c_xa
   *
   *  relax_param:
   *    <= 0 : P0 restriction / P0 prolongation => c_xa = c_xa0
   *    > 0  : P0 restriction / P1 prolongation => c_xa = c_xa0ij/icjc
   */

  /* Extradiagonal terms */

  cs_lnum_t c_n_vals = c_n_faces*isym;

# pragma omp parallel for if(c_n_vals*6 > CS_THR_MIN)
  for (cs_lnum_t c_val = 0; c_val < c_n_vals; c_val++) {
    c_xa[c_val] = 0.;
  }

  if (relax_param <= 0) {

    if (isym == 1) {

      for (face_id = 0; face_id < f_n_faces; face_id++) {
        cs_lnum_t c_face = c_coarse_face[face_id];
        if (c_face > 0) {
          c_face = c_face - 1;
          c_xa[c_face] += f_xa[face_id];
        }
        else if (c_face < 0) {
          c_face = -c_face - 1;
          c_xa[c_face] += f_xa[face_id];
        }
      }

    }

    else if (isym == 2) {

      for (face_id = 0; face_id < f_n_faces; face_id++) {
        cs_lnum_t c_face = c_coarse_face[face_id];
        if (c_face > 0) {
          c_face = c_face - 1;
          c_xa[c_face*2] += f_xa[face_id*2];
          c_xa[c_face*2 + 1] += f_xa[face_id*2 + 1];
        }
        else if (c_face < 0) {
          c_face = -c_face - 1;
          c_xa[c_face*2] += f_xa[face_id*2];
          c_xa[c_face*2 + 1] += f_xa[face_id*2 + 1];
        }
      }

    }

  }

  else if (relax_param > 0) {

    /* P0 restriction of matrices, "interior" surface: */
    /* xag0(nfacg), surfag(3,nfacgl), xagxg0(2,nfacg) */

#   pragma omp parallel for if(c_n_faces*6 > CS_THR_MIN)
    for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++) {
      c_xa0[c_face] = 0.;
      c_face_normal[3*c_face]    = 0.;
      c_face_normal[3*c_face +1] = 0.;
      c_face_normal[3*c_face +2] = 0.;
      c_xa0ij[3*c_face]    = 0.;
      c_xa0ij[3*c_face +1] = 0.;
      c_xa0ij[3*c_face +2] = 0.;
    }

    if (f_face_normal != NULL) {

      for (face_id = 0; face_id < f_n_faces; face_id++) {

        if (c_coarse_face[face_id] > 0 ) {
          cs_lnum_t c_face = c_coarse_face[face_id] -1;

          c_xa0[c_face] += f_xa0[face_id];
          c_face_normal[3*c_face]    += f_face_normal[3*face_id];
          c_face_normal[3*c_face +1] += f_face_normal[3*face_id +1];
          c_face_normal[3*c_face +2] += f_face_normal[3*face_id +2];
          c_xa0ij[3*c_face]    += f_xa0ij[3*face_id];
          c_xa0ij[3*c_face +1] += f_xa0ij[3*face_id +1];
          c_xa0ij[3*c_face +2] += f_xa0ij[3*face_id +2];
        }
        else if (c_coarse_face[face_id] < 0) {
          cs_lnum_t c_face = -c_coarse_face[face_id] -1;

          c_xa0[c_face] += f_xa0[face_id];
          c_face_normal[3*c_face]    -= f_face_normal[3*face_id];
          c_face_normal[3*c_face +1] -= f_face_normal[3*face_id +1];
          c_face_normal[3*c_face +2] -= f_face_normal[3*face_id +2];
          c_xa0ij[3*c_face]    -= f_xa0ij[3*face_id];
          c_xa0ij[3*c_face +1] -= f_xa0ij[3*face_id +1];
          c_xa0ij[3*c_face +2] -= f_xa0ij[3*face_id +2];
        }

      }

    }
    else { /* f_face_normal = NULL */

      for (face_id = 0; face_id < f_n_faces; face_id++) {

        if (c_coarse_face[face_id] > 0 ) {
          cs_lnum_t c_face = c_coarse_face[face_id] -1;
          cs_lnum_t ic = c_face_cell[c_face][0];
          cs_lnum_t jc = c_face_cell[c_face][1];

          c_xa0[c_face] += f_xa0[face_id];
          c_face_normal[3*c_face]    += c_cell_cen[3*jc]   - c_cell_cen[3*ic];
          c_face_normal[3*c_face +1] += c_cell_cen[3*jc+1] - c_cell_cen[3*ic+1];
          c_face_normal[3*c_face +2] += c_cell_cen[3*jc+2] - c_cell_cen[3*ic+2];
          c_xa0ij[3*c_face]    += f_xa0ij[3*face_id];
          c_xa0ij[3*c_face +1] += f_xa0ij[3*face_id +1];
          c_xa0ij[3*c_face +2] += f_xa0ij[3*face_id +2];
        }
        else if (c_coarse_face[face_id] < 0) {
          cs_lnum_t c_face = -c_coarse_face[face_id] -1;
          cs_lnum_t ic = c_face_cell[c_face][0];
          cs_lnum_t jc = c_face_cell[c_face][1];

          c_xa0[c_face] += f_xa0[face_id];
          c_face_normal[3*c_face]    -= c_cell_cen[3*jc]   - c_cell_cen[3*ic];
          c_face_normal[3*c_face +1] -= c_cell_cen[3*jc+1] - c_cell_cen[3*ic+1];
          c_face_normal[3*c_face +2] -= c_cell_cen[3*jc+2] - c_cell_cen[3*ic+2];
          c_xa0ij[3*c_face]    -= f_xa0ij[3*face_id];
          c_xa0ij[3*c_face +1] -= f_xa0ij[3*face_id +1];
          c_xa0ij[3*c_face +2] -= f_xa0ij[3*face_id +2];
        }

      }

    }

    /* Matrix initialized to c_xa0 */

    if (isym == 1) {
#     pragma omp parallel for if(c_n_faces > CS_THR_MIN)
      for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++)
        c_xa[c_face] = c_xa0[c_face];
    }
    else {
#     pragma omp parallel for if(c_n_faces*2 > CS_THR_MIN)
      for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++) {
        c_xa[2*c_face]    = c_xa0[c_face];
        c_xa[2*c_face +1] = c_xa0[c_face];
      }
    }

    for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++) {

      cs_lnum_t ic = c_face_cell[c_face][0];
      cs_lnum_t jc = c_face_cell[c_face][1];

      dsigjg =   (  c_cell_cen[3*jc]
                  - c_cell_cen[3*ic])    * c_face_normal[3*c_face]
               + (  c_cell_cen[3*jc +1]
                  - c_cell_cen[3*ic +1]) * c_face_normal[3*c_face +1]
               + (  c_cell_cen[3*jc +2]
                  - c_cell_cen[3*ic +2]) * c_face_normal[3*c_face +2];

      dsxaij =   c_xa0ij[3*c_face]    * c_face_normal[3*c_face]
               + c_xa0ij[3*c_face +1] * c_face_normal[3*c_face +1]
               + c_xa0ij[3*c_face +2] * c_face_normal[3*c_face +2];

      if (fabs(dsigjg) > EPZERO) {

        agij = dsxaij/dsigjg;

        /* Standard */
        c_xa[c_face*isym] = agij;
        c_xa[(c_face +1)*isym -1] = agij;

        /* Clipped matrix */
        if (agij < c_xa0[c_face] || agij > 0.) {
          c_xa[c_face*isym] = c_xa0[c_face];
          c_xa[(c_face +1)*isym -1] = c_xa0[c_face];
          if (agij < c_xa0[c_face]) n_clips_min++;
          if (agij > 0.) n_clips_max++;
        }

      }
    }

    /* Possible P1 matrix / P0 matrix relaxation defined
       by the user in cs_user_parameters.c
       (using cs_multigrid_set_coarsening_options) */

    if (isym == 1) {
      for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++)
        c_xa[c_face] = relax_param*c_xa[c_face]
                     + (1. - relax_param)*c_xa0[c_face];
    }
    else {
      for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++) {
        c_xa[2*c_face]    =         relax_param  * c_xa[2*c_face]
                            + (1. - relax_param) * c_xa0[c_face];
        c_xa[2*c_face +1] =         relax_param  * c_xa[2*c_face +1]
                            + (1. - relax_param) * c_xa0[c_face];
      }
    }

  } /* relax_param > 0 */

  /* Initialize non differential fine grid term saved in w1 */

  cs_real_t *w1 = NULL;
  BFT_MALLOC(w1, f_n_cells_ext*db_stride, cs_real_t);

  if (db_size == 1) {
#   pragma omp parallel for if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++)
      w1[ii] = f_da[ii];
  }
  else {
#   pragma omp parallel for private(jj, kk) if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      for (jj = 0; jj < db_size; jj++) {
        for (kk = 0; kk < db_size; kk++)
          w1[ii*db_stride + db_size*jj + kk]
            = f_da[ii*db_stride + db_size*jj + kk];
      }
    }
  }
# pragma omp parallel for if(f_n_cells_ext - f_n_cells > CS_THR_MIN)
  for (ii = f_n_cells*db_stride; ii < f_n_cells_ext*db_stride; ii++)
    w1[ii] = 0.;

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cell[face_id][0];
    jj = f_face_cell[face_id][1];
    for (kk = 0; kk < db_size; kk++) {
      w1[ii*db_stride + db_size*kk + kk] += f_xa[face_id*isym];
      w1[jj*db_stride + db_size*kk + kk] += f_xa[(face_id +1)*isym -1];
    }
  }

  /* Diagonal term */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t ic = 0; ic < c_n_cells_ext*db_stride; ic++)
    c_da[ic] = 0.;

  if (db_size == 1) {
    for (ii = 0; ii < f_n_cells; ii++) {
      cs_lnum_t ic = c_coarse_row[ii];
      if (ic > -1)
        c_da[ic] += w1[ii];
    }
  }
  else {
    for (ii = 0; ii < f_n_cells; ii++) {
      cs_lnum_t ic = c_coarse_row[ii];
      if (ic > -1) {
        for (jj = 0; jj < db_size; jj++) {
          for (kk = 0; kk < db_size; kk++)
            c_da[ic*db_stride + db_size*jj + kk]
              += w1[ii*db_stride + db_size*jj + kk];
        }
      }
    }
  }

  for (cs_lnum_t c_face = 0; c_face < c_n_faces; c_face++) {
    cs_lnum_t ic = c_face_cell[c_face][0];
    cs_lnum_t jc = c_face_cell[c_face][1];
    for (kk = 0; kk < db_size; kk++) {
      c_da[ic*db_stride + db_size*kk + kk] -= c_xa[c_face*isym];
      c_da[jc*db_stride + db_size*kk + kk] -= c_xa[(c_face +1)*isym -1];
    }
  }

  BFT_FREE(w1);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           fine_grid->level, coarse_grid->level);
    printf(", total = %ld\n", elapsed.count());
  }

  /* Optional verification */

  if (verbosity > 3) {
    int interp = (relax_param > 0) ? 1 : 0;
    _verify_coarse_quantities(fine_grid,
                              coarse_grid,
                              n_clips_min,
                              n_clips_max,
                              interp);
  }
}

/*----------------------------------------------------------------------------
 * Build a coarse level from a finer level for convection/diffusion case.
 *
 * parameters:
 *   fine_grid   <-- fine grid structure
 *   coarse_grid <-> coarse grid structure
 *   verbosity   <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities_conv_diff(const cs_grid_t  *fine_grid,
                                     cs_grid_t        *coarse_grid,
                                     int               verbosity)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  cs_lnum_t ic, jc, ii, jj, c_face, face_id;

  cs_real_t dsigjg, dsxaij, agij;

  cs_lnum_t f_n_cells = fine_grid->n_rows;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cols_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells_ext = coarse_grid->n_cols_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  cs_lnum_t *c_coarse_row = coarse_grid->coarse_row;
  cs_lnum_t *c_coarse_face = coarse_grid->coarse_face;

  cs_gnum_t n_clips_min = 0, n_clips_max = 0;

  cs_real_t *f_xa0ij = fine_grid->xa0ij;

  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;
  cs_real_t *c_face_normal = coarse_grid->_face_normal;

  cs_real_t *c_xa0 = coarse_grid->_xa0;
  cs_real_t *c_xa0_diff = coarse_grid->xa0_diff;
  cs_real_t *c_xa0ij = coarse_grid->xa0ij;
  cs_real_t *c_da = coarse_grid->_da;
  cs_real_t *c_xa = coarse_grid->_xa;
  cs_real_t *c_xa_conv = coarse_grid->xa_conv;
  cs_real_t *c_xa_diff = coarse_grid->xa_diff;

  cs_real_t *w1 = NULL;

  const cs_lnum_t db_size = fine_grid->db_size;
  const cs_lnum_t db_stride = db_size*db_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_face_normal = fine_grid->face_normal;
  const cs_real_t *f_xa0 = fine_grid->xa0;
  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_xa = fine_grid->xa;
  const cs_real_t *f_xa_conv = fine_grid->xa_conv;
  const cs_real_t *f_xa0_diff = fine_grid->xa0_diff;
  const cs_real_t *f_xa_diff = fine_grid->xa_diff;

  BFT_MALLOC(w1, 2*f_n_cells_ext*db_stride, cs_real_t);

  assert(fine_grid->symmetric == false);

  /* P0 restriction of matrices, "interior" surface: */
  /* xag0(nfacg), surfag(3,nfacgl), xagxg0(2,nfacg) */

# pragma omp parallel for if(c_n_faces*6 > CS_THR_MIN)
  for (c_face = 0; c_face < c_n_faces; c_face++) {
    c_xa0[2*c_face]         = 0.;
    c_xa0[2*c_face +1]      = 0.;
    c_xa0_diff[c_face]      = 0.;
    c_face_normal[3*c_face]    = 0.;
    c_face_normal[3*c_face +1] = 0.;
    c_face_normal[3*c_face +2] = 0.;
    c_xa0ij[3*c_face]    = 0.;
    c_xa0ij[3*c_face +1] = 0.;
    c_xa0ij[3*c_face +2] = 0.;
    c_xa_conv[2*c_face]    = 0.;
    c_xa_conv[2*c_face +1] = 0.;
  }

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    if (c_coarse_face[face_id] > 0 ) {
      c_face = c_coarse_face[face_id] -1;

      c_xa0[2*c_face]         += f_xa0[2*face_id];
      c_xa0[2*c_face +1]      += f_xa0[2*face_id +1];
      c_xa0_diff[c_face]      += f_xa0_diff[face_id];
      c_face_normal[3*c_face]    += f_face_normal[3*face_id];
      c_face_normal[3*c_face +1] += f_face_normal[3*face_id +1];
      c_face_normal[3*c_face +2] += f_face_normal[3*face_id +2];
      c_xa0ij[3*c_face]    += f_xa0ij[3*face_id];
      c_xa0ij[3*c_face +1] += f_xa0ij[3*face_id +1];
      c_xa0ij[3*c_face +2] += f_xa0ij[3*face_id +2];
    }
    else if (c_coarse_face[face_id] < 0) {
      c_face = -c_coarse_face[face_id] -1;

      c_xa0[2*c_face]         += f_xa0[2*face_id +1];
      c_xa0[2*c_face +1]      += f_xa0[2*face_id];
      c_xa0_diff[c_face]      += f_xa0_diff[face_id];
      c_face_normal[3*c_face]    -= f_face_normal[3*face_id];
      c_face_normal[3*c_face +1] -= f_face_normal[3*face_id +1];
      c_face_normal[3*c_face +2] -= f_face_normal[3*face_id +2];
      c_xa0ij[3*c_face]    -= f_xa0ij[3*face_id];
      c_xa0ij[3*c_face +1] -= f_xa0ij[3*face_id +1];
      c_xa0ij[3*c_face +2] -= f_xa0ij[3*face_id +2];
    }

  }

  /*  Finalize computation of matrix in c_da, c_xa */
  /*  interp = 0 : P0 restriction / P0 prolongation => c_xa = c_xa0 */
  /*  interp = 1 : P0 restriction / P1 prolongation => c_xa = c_xa0ij/icjc */

  /*  Initialization */

  const cs_real_t relax_param = coarse_grid->relaxation;

  int interp = 1;

  /* Initialize non differential fine grid term saved in w1 */

  if (db_size == 1) {
#   pragma omp parallel for if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      w1[ii] = f_da[ii];
    }
  }
  else {
#   pragma omp parallel for private(jj) if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      for (jj = 0; jj < db_size; jj++) {
        for (cs_lnum_t kk = 0; kk < db_size; kk++)
          w1[ii*db_stride + db_size*jj + kk]
            = f_da[ii*db_stride + db_size*jj + kk];
      }
    }
  }

# pragma omp parallel for if(f_n_cells_ext - f_n_cells > CS_THR_MIN)
  for (ii = f_n_cells*db_stride; ii < f_n_cells_ext*db_stride; ii++)
    w1[ii] = 0.;

  /* Finest grid: we need to separate convective and diffusive coefficients.
     By construction, extra-diagonal coefficients are expected to be
     negative, with symmetric diffusion and pure upwind convection
     components. This allows separating both parts. */

  if (fine_grid->level == 0) {
    for (face_id = 0; face_id < f_n_faces; face_id++) {
      ii = f_face_cell[face_id][0];
      jj = f_face_cell[face_id][1];
      cs_real_t f_xa_conv_0, f_xa_conv_1, f_xa_diff_s;
      if (f_xa[2*face_id] < f_xa[2*face_id+1]) {
        f_xa_diff_s = f_xa[2*face_id+1];
        f_xa_conv_0 = f_xa[2*face_id] - f_xa_diff_s;
        f_xa_conv_1 = -0.;
      }
      else {
        f_xa_diff_s = f_xa[2*face_id];
        f_xa_conv_0 = -0;
        f_xa_conv_1 = f_xa[2*face_id+1] - f_xa_diff_s;
      }
      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        w1[ii*db_stride + db_size*kk + kk] += f_xa_conv_0 + f_xa_diff_s;
        w1[jj*db_stride + db_size*kk + kk] += f_xa_conv_1 + f_xa_diff_s;
      }
      if (c_coarse_face[face_id] > 0 ) {
        c_face = c_coarse_face[face_id] -1;
        c_xa_conv[2*c_face]    += f_xa_conv_0;
        c_xa_conv[2*c_face +1] += f_xa_conv_1;
      }
      else if (c_coarse_face[face_id] < 0) {
        c_face = -c_coarse_face[face_id] -1;
        c_xa_conv[2*c_face]    += f_xa_conv_1;
        c_xa_conv[2*c_face +1] += f_xa_conv_0;
      }
    }
  }

  /* Coarser grids */

  else {
    for (face_id = 0; face_id < f_n_faces; face_id++) {
      ii = f_face_cell[face_id][0];
      jj = f_face_cell[face_id][1];
      for (cs_lnum_t kk = 0; kk < db_size; kk++) {
        w1[ii*db_stride + db_size*kk + kk] +=   f_xa_conv[2*face_id]
                                              + f_xa_diff[face_id];
        w1[jj*db_stride + db_size*kk + kk] +=   f_xa_conv[2*face_id +1]
                                              + f_xa_diff[face_id];
      }
      if (c_coarse_face[face_id] > 0 ) {
        c_face = c_coarse_face[face_id] -1;
        c_xa_conv[2*c_face]    += f_xa_conv[2*face_id];
        c_xa_conv[2*c_face +1] += f_xa_conv[2*face_id +1];
      }
      else if (c_coarse_face[face_id] < 0) {
        c_face = -c_coarse_face[face_id] -1;
        c_xa_conv[2*c_face]    += f_xa_conv[2*face_id +1];
        c_xa_conv[2*c_face +1] += f_xa_conv[2*face_id];
      }
    }
  }

  /* Initialize coarse matrix storage on (c_da, c_xa) */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (ic = 0; ic < c_n_cells_ext*db_stride; ic++) {
    c_da[ic] = 0.;
  }

  /* Extradiagonal terms */
  /* (symmetric matrices for now, even with non symmetric storage isym = 2) */

  /* Matrix initialized to c_xa0 (interp=0) */

# pragma omp parallel for if(c_n_faces*2 > CS_THR_MIN)
  for (c_face = 0; c_face < c_n_faces; c_face++) {
    c_xa[2*c_face]         = c_xa0[2*c_face];
    c_xa[2*c_face +1]      = c_xa0[2*c_face +1];
    c_xa_diff[c_face]      = c_xa0_diff[c_face];
  }

  if (interp == 1) {
   for (c_face = 0; c_face < c_n_faces; c_face++) {

      ic = c_face_cell[c_face][0];
      jc = c_face_cell[c_face][1];

      dsigjg =   (  c_cell_cen[3*jc]
                  - c_cell_cen[3*ic])    * c_face_normal[3*c_face]
               + (  c_cell_cen[3*jc +1]
                  - c_cell_cen[3*ic +1]) * c_face_normal[3*c_face +1]
               + (  c_cell_cen[3*jc +2]
                  - c_cell_cen[3*ic +2]) * c_face_normal[3*c_face +2];

      dsxaij =   c_xa0ij[3*c_face]    * c_face_normal[3*c_face]
               + c_xa0ij[3*c_face +1] * c_face_normal[3*c_face +1]
               + c_xa0ij[3*c_face +2] * c_face_normal[3*c_face +2];

      if (fabs(dsigjg) > EPZERO) {

        agij = dsxaij/dsigjg;

        /* Standard */
        c_xa_diff[c_face] = agij;

        /* Clipped matrix */
        if (agij < c_xa0_diff[c_face] || agij > 0.) {
          c_xa_diff[c_face] = c_xa0_diff[c_face];
          if (agij < c_xa0_diff[c_face]) n_clips_min++;
          if (agij > 0.) n_clips_max++;
        }

      }
    }

    /* Possible P1 matrix / P0 matrix relaxation defined
       by the user in cs_user_parameters.c
       (using cs_multigrid_set_coarsening_options) */
    for (c_face = 0; c_face < c_n_faces; c_face++) {
      c_xa_diff[c_face] =         relax_param  * c_xa_diff[c_face]
                          + (1. - relax_param) * c_xa0_diff[c_face];
    }

  } /* endif interp == 1 */

  if (interp != 0 && interp != 1)
    bft_error(__FILE__, __LINE__, 0, "interp incorrectly defined.");

  /* Diagonal term */

  if (db_size == 1) {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_row[ii];
      if (ic > -1)
        c_da[ic] += w1[ii];
    }
  }
  else {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_row[ii];
      if (ic > -1) {
        for (jj = 0; jj < db_size; jj++) {
          for (cs_lnum_t kk = 0; kk < db_size; kk++)
            c_da[ic*db_stride + db_size*jj + kk]
              += w1[ii*db_stride + db_size*jj + kk];
        }
      }
    }
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    for (cs_lnum_t kk = 0; kk < db_size; kk++) {
      c_da[ic*db_stride + db_size*kk + kk] -=   c_xa_conv[2*c_face]
                                              + c_xa_diff[c_face];
      c_da[jc*db_stride + db_size*kk + kk] -=   c_xa_conv[2*c_face +1]
                                              + c_xa_diff[c_face];
    }
  }

  /* Convection/diffusion matrix */

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    c_xa[2*c_face]    = c_xa_conv[2*c_face]    + c_xa_diff[c_face];
    c_xa[2*c_face +1] = c_xa_conv[2*c_face +1] + c_xa_diff[c_face];
  }

  BFT_FREE(w1);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           fine_grid->level, coarse_grid->level);
    printf(", total = %ld\n", elapsed.count());
  }

  /* Optional verification */

  if (verbosity > 3)
    _verify_coarse_quantities(fine_grid,
                              coarse_grid,
                              n_clips_min,
                              n_clips_max,
                              interp);
}

/*----------------------------------------------------------------------------
 * Build a coarse level from a finer level with an MSR matrix.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-> Coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities_msr(const cs_grid_t  *fine_grid,
                               cs_grid_t        *coarse_grid)

{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  int n_f_threads = cs_parall_n_threads(fine_grid->n_rows, CS_THR_MIN);

  const cs_lnum_t db_size = fine_grid->db_size;
  const cs_lnum_t db_stride = db_size*db_size;

  const cs_lnum_t eb_size = fine_grid->eb_size;
  const cs_lnum_t eb_stride = eb_size*eb_size;

  const cs_lnum_t f_n_rows = fine_grid->n_rows;

  const cs_lnum_t c_n_rows = coarse_grid->n_rows;
  const cs_lnum_t *f_c_row = coarse_grid->coarse_row;

  /* Fine matrix in the MSR format */

  const cs_lnum_t  *f_row_index, *f_col_id;
  const cs_real_t  *f_d_val, *f_x_val;

  cs_matrix_get_msr_arrays(fine_grid->matrix,
                           &f_row_index,
                           &f_col_id,
                           &f_d_val,
                           &f_x_val);

  /* Determine reverse coarse to fine adjacency */

  cs_lnum_t *c_f_row_index = nullptr, *c_f_row_ids = nullptr;
  cs_lnum_t *c_row_index_0 = nullptr;

  _coarse_to_fine_adjacency_msr(f_n_rows,
                                c_n_rows,
                                coarse_grid->alloc_mode,
                                n_f_threads,
                                f_c_row,
                                f_row_index,
                                &c_f_row_index,
                                &c_f_row_ids,
                                &c_row_index_0);

  /* Coarse matrix elements in the MSR format */

  /* Build structure
     --------------- */

  cs_lnum_t *c_row_index,  *c_col_id;

  _coarse_msr_struct(f_n_rows,
                     c_n_rows,
                     coarse_grid->alloc_mode,
                     f_row_index,
                     f_col_id,
                     f_c_row,
                     c_f_row_index,
                     c_f_row_ids,
                     c_row_index_0,
                     &c_row_index,
                     &c_col_id);

  /* Assign values
     ------------- */

  cs_lnum_t c_nnz = c_row_index[c_n_rows];

  cs_real_t *restrict c_d_val, *restrict c_x_val;
  CS_MALLOC_HD(c_d_val, c_n_rows*db_stride, cs_real_t, coarse_grid->alloc_mode);
  CS_MALLOC_HD(c_x_val, c_nnz*eb_stride, cs_real_t, coarse_grid->alloc_mode);

  /* Scalar case */

  if (db_size == 1) {

    #pragma omp parallel for num_threads(n_f_threads)
    for (cs_lnum_t ic = 0; ic < c_n_rows; ic++) {

      const cs_lnum_t s_id = c_row_index[ic];
      const cs_lnum_t n_cols = c_row_index[ic+1] - s_id;
      const cs_lnum_t n_vals = n_cols;
      cs_real_t *restrict row_x_vals = c_x_val + s_id;

      c_d_val[ic] = 0;
      for (cs_lnum_t ll = 0; ll < n_vals; ll++)
        row_x_vals[ll] = 0;

      cs_lnum_t r_s_id = c_f_row_index[ic];
      cs_lnum_t r_e_id = c_f_row_index[ic+1];
      for (cs_lnum_t r_idx = r_s_id; r_idx < r_e_id; r_idx++) {
        cs_lnum_t ii = c_f_row_ids[r_idx];

        /* Diagonal term */

        c_d_val[ic] += f_d_val[ii];

        /* Extra-diagonal terms */

        cs_lnum_t f_s_id = f_row_index[ii];
        cs_lnum_t f_e_id = f_row_index[ii+1];
        for (cs_lnum_t jj = f_s_id; jj < f_e_id; jj++) {
          cs_lnum_t jc = f_c_row[f_col_id[jj]];
          if (jc > -1) {
            if (ic != jc) {
              /* ids are sorted, so binary search possible */
              cs_lnum_t k = _l_id_binary_search(n_cols, jc, c_col_id + s_id);
              row_x_vals[k] += f_x_val[jj];
            }
            else { /* ic == jc */
              c_d_val[ic] += f_x_val[jj];
            }
          } /* If fine row has matching coarse row */
        } /* Loop in fine columns */
      } /* Loop on fine rows */

    } /* OpenMP loop on coarse rows */

  }

  /* General case */

  else {

    #pragma omp parallel for num_threads(n_f_threads)
    for (cs_lnum_t ic = 0; ic < c_n_rows; ic++) {

      const cs_lnum_t s_id = c_row_index[ic];
      const cs_lnum_t n_cols = c_row_index[ic+1] - s_id;
      const cs_lnum_t n_vals = n_cols*eb_stride;
      cs_real_t *restrict row_x_vals = c_x_val + s_id*eb_stride;

      for (cs_lnum_t ll = 0; ll < db_stride; ll++) {
        c_d_val[ic*db_stride + ll] = 0;
      }

      for (cs_lnum_t ll = 0; ll < n_vals; ll++)
        row_x_vals[ll] = 0;

      cs_lnum_t r_s_id = c_f_row_index[ic];
      cs_lnum_t r_e_id = c_f_row_index[ic+1];
      for (cs_lnum_t r_idx = r_s_id; r_idx < r_e_id; r_idx++) {
        cs_lnum_t ii = c_f_row_ids[r_idx];

        /* Diagonal terms */

        for (cs_lnum_t jj = 0; jj < db_stride; jj++) {
          c_d_val[ic*db_stride + jj] += f_d_val[ii*db_stride + jj];
        }

        /* Extra-diagonal terms */

        cs_lnum_t f_s_id = f_row_index[ii];
        cs_lnum_t f_e_id = f_row_index[ii+1];
        for (cs_lnum_t jj = f_s_id; jj < f_e_id; jj++) {
          cs_lnum_t jc = f_c_row[f_col_id[jj]];
          if (jc > -1) {
            if (ic != jc) {
              /* ids are sorted, so binary search possible */
              cs_lnum_t k = _l_id_binary_search(n_cols, jc, c_col_id + s_id);
              assert(k > -1);
              for (cs_lnum_t l = 0; l < eb_stride; l++)
                row_x_vals[k*eb_stride + l] += f_x_val[jj*eb_stride + l];
            }
            else { /* ic == jc */
              if (eb_size == 1) {
                for (cs_lnum_t kk = 0; kk < db_size; kk++) {
                  /* contribution to diagonal by isotropic extra-diagonal */
                  c_d_val[ic*db_stride + db_size*kk + kk] += f_x_val[jj];
                }
              }
              else {
                for (cs_lnum_t kk = 0; kk < eb_stride; kk++) {
                  /* contribution to diagonal by complete extra-diagonal */
                  c_d_val[ic*db_stride + kk] += f_x_val[jj*eb_stride + kk];
                }
              }
            }
          } /* If fine row has matching coarse row */
        } /* Loop in fine columns */
      } /* Loop on fine rows */

    } /* OPenMP loop on coarse rows */

  }

  /* Now build matrix */

  _build_coarse_matrix_msr(coarse_grid, fine_grid->symmetric,
                           c_row_index, c_col_id,
                           c_d_val, c_x_val);

  /* Free working arrays */

  CS_FREE(c_row_index_0);
  CS_FREE(c_f_row_ids);
  CS_FREE(c_f_row_index);

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           fine_grid->level, coarse_grid->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief compute coarse mesh quantities from fine quantities with MSR matrix.
 *
 * Also build direct a fine to coarse MSR matrix index.
 *
 * Though thie can be computed on the fly, in cases whare multiple passes
 * are required, this will avoid extra searches.
 *
 * Whether running searches multiple times on the coarse structure,
 * or requiring more bandwith but with a simpler access (choice here)
 * provides best performance would need exploring.
 *
 * The caller is responsible for freeing the returned array.
 *
 * \param[in]       f              fine grid structure
 * \param[in, out]  c              coarse grid structure
 * \param[in]       f_row_index    fine matrix row index
 * \param[in]       f_col_id       fine matrix column ids
 * \param[in]       c_f_row_index  coarse to fine rows index
 * \param[in]       c_f_row_ids    coarse to fine row ids
 * \param[in]       c_row_index    coarse row index
 * \param[in]       c_col_id       coarse column ids
 *
 * \return  fine to coarse index mapping, with -1 for unmapped faces,
 *          and -1 - coarse_row_id when an index maps to the diagonal
 *          (i.e. if the index associates 2 merged rows).
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t *
_coarse_mesh_quantities_msr(const cs_grid_t  *f,
                            cs_grid_t        *c,
                            const cs_lnum_t  *f_row_index,
                            const cs_lnum_t  *f_col_id,
                            const cs_lnum_t  *c_f_row_index,
                            const cs_lnum_t  *c_f_row_ids,
                            const cs_lnum_t  *c_row_index,
                            const cs_lnum_t  *c_col_id)
{
  std::chrono::high_resolution_clock::time_point t_0;
  std::chrono::high_resolution_clock::time_point t_1;
  if (cs_glob_timer_kernels_flag > 0)
    t_0 = std::chrono::high_resolution_clock::now();

  const cs_lnum_t f_n_rows = f->n_rows;

  const cs_lnum_t c_n_rows = c->n_rows;
  const cs_lnum_t c_n_cols = c->n_cols_ext;
  const cs_lnum_t c_n_faces = c->n_faces;

  const cs_lnum_t *f_c_row = c->coarse_row;

  const cs_alloc_mode_t alloc_mode = f->alloc_mode;

  const cs_lnum_t   *f_cell_idx;
  const cs_lnum_t   *f_cell_face;
  const short int   *f_cell_face_sgn;
  const cs_real_3_t *f_cell_cen;
  const cs_real_3_t *f_face_normal;
  const cs_real_t   *f_cell_vol;

  cs_matrix_get_mesh_association(f->matrix,
                                 &f_cell_idx,
                                 &f_cell_face,
                                 &f_cell_face_sgn,
                                 &f_cell_cen,
                                 &f_cell_vol,
                                 &f_face_normal);

  const cs_lnum_t  *c_cell_face = c->cell_face;

  int n_loc_threads = cs_parall_n_threads(f_n_rows, CS_THR_MIN);

  cs_real_3_t  *c_cell_cen;
  cs_real_t    *c_cell_vol;
  cs_real_3_t  *c_face_normal;

  CS_MALLOC_HD(c_cell_cen, c_n_cols, cs_real_3_t, alloc_mode);
  CS_MALLOC_HD(c_cell_vol, c_n_cols, cs_real_t, alloc_mode);
  CS_MALLOC_HD(c_face_normal, c_n_faces, cs_real_3_t, alloc_mode);

  #pragma omp parallel num_threads(n_loc_threads)
  for (cs_lnum_t i = 0; i < c_n_faces; i++) {
    c_face_normal[i][0] = 0.;
    c_face_normal[i][1] = 0.;
    c_face_normal[i][2] = 0.;
  }

  /* Allocate mapping */

  cs_lnum_t f_nnz = f_row_index[f_n_rows];
  cs_lnum_t *fc_m_idx;
  CS_MALLOC_HD(fc_m_idx, f_nnz, cs_lnum_t, f->alloc_mode);

  /* Loop on coarse rows */

  #pragma omp parallel for  num_threads(n_loc_threads)
  for (cs_lnum_t ii_c = 0; ii_c < c_n_rows; ii_c++) {

    const cs_lnum_t c_s_idx = c_row_index[ii_c];
    const cs_lnum_t n_c_cols = c_row_index[ii_c+1] - c_s_idx;

    cs_real_t *restrict _c_cell_cen = c_cell_cen[ii_c];

    for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++)
      _c_cell_cen[coo_id] = 0.;
    c_cell_vol[ii_c] = 0;

    /* Loop on matching fine rows */
    cs_lnum_t r_s_id = c_f_row_index[ii_c];
    cs_lnum_t r_e_id = c_f_row_index[ii_c+1];
    for (cs_lnum_t r_idx = r_s_id; r_idx < r_e_id; r_idx++) {
      cs_lnum_t ii = c_f_row_ids[r_idx];

      /* Cell quantities */
      for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++)
        _c_cell_cen[coo_id] += f_cell_cen[ii][coo_id];
      c_cell_vol[ii_c] += f_cell_vol[ii];

      cs_lnum_t f_s_idx = f_row_index[ii];
      cs_lnum_t f_e_idx = f_row_index[ii+1];
      for (cs_lnum_t f_idx = f_s_idx; f_idx < f_e_idx; f_idx++) {
        cs_lnum_t jj = f_col_id[f_idx];
        cs_lnum_t jj_c = f_c_row[jj];
        if (jj_c < 0)
          fc_m_idx[f_idx] = -1;
        else if (jj_c == ii_c)
          fc_m_idx[f_idx] = -1 - jj_c;
        else {
          cs_lnum_t k = _l_id_binary_search(n_c_cols,
                                            jj_c,
                                            c_col_id + c_s_idx);
          assert(k > -1);
          cs_lnum_t c_idx = c_s_idx + k;
          fc_m_idx[f_idx] = c_idx;

          if (jj_c > ii_c) {
            int sgn = f_cell_face_sgn[f_idx];
            if (ii > jj)
              sgn *= -1;

            /* Note: coarse face id is handled by a single
               thread, since it is associated with the coarse
               upper diagonal matrix. */

            cs_lnum_t f_face_id = f_cell_face[f_idx];
            cs_lnum_t c_face_id = c_cell_face[c_idx];

            for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++)
              c_face_normal[c_face_id][coo_id]
                += sgn * f_face_normal[f_face_id][coo_id];
          }
        }
      } /* Loop on fine columns */

    } /* Loop on fine rows associated to current coarse row */

  } /* Loop on coarse rows */

  /* Add to coarse grid */

  c->_cell_cen = (cs_real_t *)c_cell_cen;
  c->cell_cen = c->cell_cen;

  c->_cell_vol = (cs_real_t *)c_cell_vol;
  c->cell_vol = c->_cell_vol;

  c->_face_normal = (cs_real_t *)c_face_normal;
  c->face_normal = c->_face_normal;

  if (cs_glob_timer_kernels_flag > 0) {
    t_1 = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_1 - t_0);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           f->level, c->level);
    printf(", total = %ld\n", elapsed.count());
  }

  /* Synchronize grid's geometric quantities */

  if (c->halo != NULL) {
    cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_cell_cen, 3);
    if (c->halo->n_transforms > 0)
      cs_halo_perio_sync_coords(c->halo, CS_HALO_STANDARD, c->_cell_cen);

    cs_halo_sync_var(c->halo, CS_HALO_STANDARD, c->_cell_vol);
  }

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_2 = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_2 - t_1);
    printf("%d: %s halo exchange (level %d)", cs_glob_rank_id, __func__,
           c->level);
    printf(", total = %ld\n", elapsed.count());
  }

  return fc_m_idx;
}

/*----------------------------------------------------------------------------
 * Build a coarse level from a finer level with an MSR matrix.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-> Coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities_msr_with_faces(const cs_grid_t  *fine_grid,
                                          cs_grid_t        *coarse_grid)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  int n_f_threads = cs_parall_n_threads(fine_grid->n_rows, CS_THR_MIN);

  const cs_lnum_t f_n_rows = fine_grid->n_rows;

  const cs_lnum_t c_n_rows = coarse_grid->n_rows;
  const cs_lnum_t *f_c_row = coarse_grid->coarse_row;

  /* Fine matrix in the MSR format */

  const cs_lnum_t  *f_row_index, *f_col_id;
  const cs_real_t  *f_d_val, *f_x_val;

  cs_matrix_get_msr_arrays(fine_grid->matrix,
                           &f_row_index,
                           &f_col_id,
                           &f_d_val,
                           &f_x_val);

  /* Determine reverse coarse to fine adjacency */

  cs_lnum_t *c_f_row_index = nullptr, *c_f_row_ids = nullptr;
  cs_lnum_t *c_row_index_0 = nullptr;

  _coarse_to_fine_adjacency_msr(f_n_rows,
                                c_n_rows,
                                coarse_grid->alloc_mode,
                                n_f_threads,
                                f_c_row,
                                f_row_index,
                                &c_f_row_index,
                                &c_f_row_ids,
                                &c_row_index_0);

  /* Coarse matrix elements in the MSR format */

  /* Build structure
     --------------- */

  cs_lnum_t *c_row_index,  *c_col_id;

  _coarse_msr_struct(f_n_rows,
                     c_n_rows,
                     coarse_grid->alloc_mode,
                     f_row_index,
                     f_col_id,
                     f_c_row,
                     c_f_row_index,
                     c_f_row_ids,
                     c_row_index_0,
                     &c_row_index,
                     &c_col_id);

  CS_FREE(c_row_index_0);

  /* Compute face adjacency
     ---------------------- */

  _msr_face_adjacency(coarse_grid,
                      c_row_index,
                      c_col_id);

  /* Compute fine face->coarse adjacency */

  cs_lnum_t *fine_to_coarse_idx
    = _coarse_mesh_quantities_msr(fine_grid,
                                  coarse_grid,
                                  f_row_index,
                                  f_col_id,
                                  c_f_row_index,
                                  c_f_row_ids,
                                  c_row_index,
                                  c_col_id);

  CS_FREE(fine_to_coarse_idx);

#if 0
  {
  const cs_lnum_t f_n_faces = fine_grid->n_faces;

  int n_loc_threads = cs_parall_n_threads(f_n_faces, CS_THR_MIN);

  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;
  cs_real_t *c_face_normal = coarse_grid->_face_normal;

  cs_real_t *c_xa0 = coarse_grid->_xa0;
  cs_real_t *c_xa0ij = coarse_grid->xa0ij;

  /* Build fine-> coarse face mapping */

  cs_lnum_t *fc_face;
  CS_MALLOC_HD(fc_face, f_n_faces, cs_lnum_t, fine_grid->alloc_mode);

  #pragma omp parallel num_threads(n_loc_threads)
  for (cs_lnum_t ii = 0; ii < f_n_faces; ii++) {
    fc_face[ii] = 0;
  }

  #pragma omp parallel for num_threads(n_loc_threads)
  for (cs_lnum_t ic = 0; ic < c_n_rows; ic++) {

    const cs_lnum_t s_id = c_row_index[ic];
    const cs_lnum_t n_cols = c_row_index[ic+1] - s_id;
    const cs_lnum_t n_vals = n_cols;
    cs_real_t *restrict row_x_vals = c_x_val + s_id;

    c_d_val[ic] = 0;
    for (cs_lnum_t ll = 0; ll < n_vals; ll++)
      row_x_vals[ll] = 0;

    cs_lnum_t r_s_id = c_f_row_index[ic];
    cs_lnum_t r_e_id = c_f_row_index[ic+1];
    for (cs_lnum_t r_idx = r_s_id; r_idx < r_e_id; r_idx++) {
      cs_lnum_t ii = c_f_row_ids[r_idx];

      cs_lnum_t f_s_id = f_row_index[ii];
      cs_lnum_t f_e_id = f_row_index[ii+1];
      for (cs_lnum_t jj = f_s_id; jj < f_e_id; jj++) {
        cs_lnum_t jc = f_c_row[f_col_id[jj]];
        if (jc > ic) { /* Has matching coarse row, + upper diagonal part
                          (matching coarse face numbering) */
          /* ids are sorted, so binary search possible */
          cs_lnum_t k = _l_id_binary_search(n_cols, jc, c_col_id + s_id);
          cs_lnum_t f_face_id = f_c2f[jj];
          cs_lnum_t f_face_sgn = f_c2f_sgn[jj];
          f_c_face[s_id + k] = (f_face_id+1) * f_face_sgn;
        }
      } /* Loop in fine columns */
    } /* Loop on fine rows */

  } /* OpenMP loop on coarse rows */
#endif

  /* Assign values
     ------------- */

  /* TODO */

  CS_FREE(c_col_id);
  CS_FREE(c_row_index);

  /* Free working arrays */

  CS_FREE(c_f_row_ids);
  CS_FREE(c_f_row_index);

  /* Now build matrix */

#if 0
  _build_coarse_matrix_msr(coarse_grid, fine_grid->symmetric,
                           c_row_index, c_col_id,
                           c_d_val, c_x_val);
#endif

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d -> %d)", cs_glob_rank_id, __func__,
           fine_grid->level, coarse_grid->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Build edge-based matrix values from MSR matrix.
 *
 * parameters:
 *   grid <-> grid structure
 *----------------------------------------------------------------------------*/

static void
_native_from_msr(cs_grid_t  *g)
{
  const cs_lnum_t db_stride = g->db_size*g->db_size;
  const cs_lnum_t eb_stride = g->eb_size*g->eb_size;

  const cs_lnum_t n_rows = g->n_rows;
  const cs_lnum_t n_cols_ext = g->n_cols_ext;

  /* Matrix in the MSR format */

  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(g->matrix,
                           &row_index,
                           &col_id,
                           &d_val,
                           &x_val);

  {
    CS_FREE(g->_da);
    BFT_MALLOC(g->_da, db_stride*n_cols_ext, cs_real_t);
    g->da = g->_da;

    for (cs_lnum_t i = 0; i < n_rows; i++) {
      for (cs_lnum_t l = 0; l < eb_stride; l++)
        g->_da[i*db_stride + l] = d_val[i*db_stride + l];
    }

    for (cs_lnum_t i = n_rows; i < n_cols_ext; i++) {
      for (cs_lnum_t l = 0; l < eb_stride; l++)
        g->_da[i*db_stride + l] = 0;
    }
  }

  if (g->symmetric) {
    BFT_REALLOC(g->_face_cell, row_index[n_rows], cs_lnum_2_t);
    BFT_REALLOC(g->_xa, eb_stride*row_index[n_rows], cs_real_t);

    cs_lnum_t n_edges = 0;

    for (cs_lnum_t i = 0; i < n_rows; i++) {
      cs_lnum_t s_id = row_index[i];
      cs_lnum_t e_id = row_index[i+1];
      for (cs_lnum_t j = s_id; j < e_id; j++) {
        cs_lnum_t k = col_id[j];
        if (k <= i)
          continue;
        g->_face_cell[n_edges][0] = i;
        g->_face_cell[n_edges][1] = k;
        for (cs_lnum_t l = 0; l < eb_stride; l++)
          g->_xa[n_edges*eb_stride + l] = x_val[j*eb_stride + l];
        n_edges += 1;
      }
    }

    g->n_faces = n_edges;
    BFT_REALLOC(g->_face_cell, n_edges, cs_lnum_2_t);
    BFT_REALLOC(g->_xa, eb_stride*n_edges, cs_real_t);
    g->face_cell = (const cs_lnum_2_t *)(g->_face_cell);
    g->xa = g->_xa;
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: currently only implemented for symmetric cases.",
              __func__);

  cs_matrix_destroy(&(g->_matrix));
  g->matrix = NULL;
  cs_matrix_structure_destroy(&(g->matrix_struct));
}

/*----------------------------------------------------------------------------
 * Build matrix from edge-based matrix values.
 *
 * parameters:
 *   grid <-> grid structure
 *----------------------------------------------------------------------------*/

static void
_matrix_from_native(cs_matrix_type_t   cm_type,
                    cs_grid_t         *g)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  g->matrix_struct = cs_matrix_structure_create(cm_type,
                                                g->n_rows,
                                                g->n_cols_ext,
                                                g->n_faces,
                                                g->face_cell,
                                                g->halo,
                                                NULL);

  g->_matrix = cs_matrix_create(g->matrix_struct);
  cs_matrix_set_alloc_mode(g->_matrix, g->alloc_mode);

  cs_matrix_set_coefficients(g->_matrix,
                             g->symmetric,
                             g->db_size,
                             g->eb_size,
                             g->n_faces,
                             g->face_cell,
                             g->da,
                             g->xa);

  g->matrix = g->_matrix;

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d)", cs_glob_rank_id, __func__, g->level);
    printf(", total = %ld\n", elapsed.count());
  }
}

/*----------------------------------------------------------------------------
 * Project coarse grid row numbers to parent grid.
 *
 * parameters:
 *   c  <-- coarse grid structure
 *----------------------------------------------------------------------------*/

static void
_project_coarse_row_to_parent(cs_grid_t  *c)
{
  assert(c != NULL);

  const cs_grid_t  *f = c->parent;

  assert(f != NULL);

  if (f->parent != NULL) { /* level > 0*/

    cs_lnum_t *c_coarse_row = c->coarse_row;
    cs_lnum_t *f_coarse_row = f->coarse_row;

    cs_lnum_t f_n_cols = f->parent->n_cols_ext;

    for (cs_lnum_t ii = 0; ii < f_n_cols; ii++) {
      cs_lnum_t ic = f_coarse_row[ii];
      if (ic >= 0)
        f_coarse_row[ii] = c_coarse_row[ic];
    }

  }
}

/*----------------------------------------------------------------------------
 * Build locally empty matrix.
 *
 * parameters:
 *   c           <-> coarse grid structure
 *   matrix_type <-- matrix type
 *----------------------------------------------------------------------------*/

static void
_build_coarse_matrix_null(cs_grid_t         *c,
                          cs_matrix_type_t   matrix_type)
{
  cs_matrix_structure_t  *ms
    = cs_matrix_structure_create(matrix_type,
                                 0,
                                 0,
                                 0,
                                 NULL,
                                 NULL,
                                 NULL);

  c->matrix_struct = ms;

  c->_matrix = cs_matrix_create(c->matrix_struct);
  c->matrix = c->_matrix;

  cs_matrix_set_alloc_mode(c->_matrix, c->alloc_mode);

  cs_matrix_set_coefficients(c->_matrix,
                             c->symmetric,
                             c->db_size,
                             c->eb_size,
                             0,
                             NULL,
                             NULL,
                             NULL);
}

/*----------------------------------------------------------------------------
 * Compute fine row integer values from coarse row values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_num   --> Variable (rank) defined on coarse grid rows
 *   f_num   <-- Variable (rank) defined on fine grid rows
 *----------------------------------------------------------------------------*/

static void
_prolong_row_int(const cs_grid_t  *c,
                 const cs_grid_t  *f,
                 int              *c_num,
                 int              *f_num)
{
  const cs_lnum_t *coarse_row;
  const int *_c_num = c_num;

  cs_lnum_t f_n_rows = f->n_rows;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_row != NULL || f_n_rows == 0);
  assert(f_num != NULL);
  assert(c_num != NULL);

#if defined(HAVE_MPI)
  _scatter_row_int(c, c_num);
#endif

  /* Set fine values (possible penalization at first level) */

  coarse_row = c->coarse_row;

# pragma omp parallel for if(f_n_rows > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < f_n_rows; ii++) {
    cs_lnum_t i = coarse_row[ii];
    if (i >= 0)
      f_num[ii] = _c_num[i];
    else
      f_num[ii] = -1;
  }
}

/*============================================================================
 * Semi-private function definitions
 *
 * The following functions are intended to be used by the multigrid layer
 * (cs_multigrid.c), not directly by the user, so they are no more
 * documented than private static functions)
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Create base grid by mapping from shared mesh values.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   n_faces        <-- Local number of faces
 *   db_size        <-- Block sizes for diagonal
 *   eb_size        <-- Block sizes for diagonal
 *   face_cell      <-- Face -> cells connectivity
 *   a              <-- Associated matrix
 *   conv_diff      <-- Convection-diffusion mode
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_shared(cs_lnum_t              n_faces,
                           cs_lnum_t              db_size,
                           cs_lnum_t              eb_size,
                           const cs_lnum_2_t     *face_cell,
                           const cs_matrix_t     *a,
                           bool                   conv_diff)
{
  cs_grid_t *g = NULL;

  /* Create empty structure and map base data */

  g = _create_grid();

  g->alloc_mode = cs_matrix_get_alloc_mode(a);
  g->level = 0;
  g->conv_diff = conv_diff;
  g->symmetric = cs_matrix_is_symmetric(a);
  g->use_faces = false;

  g->db_size = db_size;
  g->eb_size = eb_size;

  g->n_rows = cs_matrix_get_n_rows(a);
  g->n_cols_ext = cs_matrix_get_n_columns(a);
  g->n_faces = n_faces;
  g->n_g_rows = g->n_rows;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t _n_rows = g->n_rows;
    MPI_Allreduce(&_n_rows, &(g->n_g_rows), 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  if (cs_matrix_is_mapped_from_native(a)) {
    g->face_cell = face_cell;
    g->use_faces = true;
  }

  g->relaxation = 0;

  const cs_real_t *cell_vol;
  const cs_real_3_t *cell_cen, *face_normal;

  cs_matrix_get_mesh_association(a,
                                 NULL,
                                 NULL,
                                 NULL,
                                 &cell_cen,
                                 &cell_vol,
                                 &face_normal);

  g->cell_cen = (const cs_real_t *)cell_cen;
  g->cell_vol = (const cs_real_t *)cell_vol;
  g->face_normal = (const cs_real_t *)face_normal;

  g->halo = cs_matrix_get_halo(a);

  /* Set shared matrix coefficients */

  if (cs_matrix_is_mapped_from_native(a)) {
    g->da = cs_matrix_get_diagonal(a);
    g->xa= cs_matrix_get_extra_diagonal(a);
  }

  if (g->face_cell != NULL) {

    /* Build symmetrized extra-diagonal terms if necessary,
       or point to existing terms if already symmetric */

    if (g->symmetric == true) {
      g->xa0 = g->xa;
      g->_xa0 = NULL;
    }
    else if (g->conv_diff) {
      g->xa0  = g->xa;
      g->_xa0 = NULL;
    }
    else {
      BFT_MALLOC(g->_xa0, n_faces, cs_real_t);
      for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++)
        g->_xa0[face_id] = 0.5 * (g->xa[face_id*2] + g->xa[face_id*2+1]);
      g->xa0 = g->_xa0;
    }

    /* Compute multigrid-specific terms */

    BFT_MALLOC(g->xa0ij, n_faces*3, cs_real_t);

    const cs_real_t *restrict g_xa0 = g->xa0;

    if (g->conv_diff) {
      BFT_MALLOC(g->xa0_diff, n_faces, cs_real_t);

      /* Finest grid: we need to separate convective and diffusive coefficients.
         By construction, extra-diagonal coefficients are expected to be
         negative, with symmetric diffusion and pure upwind convection
         components. This allows separating both parts. */

#     pragma omp parallel for  if(n_faces > CS_THR_MIN)
      for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
        if (g->xa[2*face_id] < g->xa[2*face_id+1]) {
          g->xa0_diff[face_id] = g->xa[2*face_id+1];
        }
        else {
          g->xa0_diff[face_id] = g->xa[2*face_id];
        }
      }

      g_xa0 = g->xa0_diff;
    }

#   pragma omp parallel for  if(n_faces > CS_THR_MIN)
    for (cs_lnum_t face_id = 0; face_id < n_faces; face_id++) {
      cs_lnum_t i0 = face_cell[face_id][0];
      cs_lnum_t i1 = face_cell[face_id][1];
      for (cs_lnum_t kk = 0; kk < 3; kk++) {
        g->xa0ij[face_id*3 + kk] =   g_xa0[face_id]
                                   * (  cell_cen[i1][kk]
                                      - cell_cen[i0][kk]);
      }
    }

  }

  g->matrix_struct = NULL;
  g->matrix = a;
  g->_matrix = NULL;

  return g;
}

/*----------------------------------------------------------------------------
 * Create base grid by mapping from parent (possibly shared) matrix.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   a       <-- associated matrix
 *   n_ranks <-- number of active ranks (<= 1 to restrict to local values)
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_parent(const cs_matrix_t  *a,
                           int                 n_ranks)
{
  cs_grid_t *g = NULL;

  /* Create empty structure and map base data */

  g = _create_grid();

  bool local = true;
  const cs_halo_t *h = cs_matrix_get_halo(a);
  if (h != NULL) {
    local = false;
    if (h->n_c_domains == 1) {
      if (h->c_domain_rank[0] == cs_glob_rank_id)
        local = true;
    }
  }

  if (n_ranks > 1 || local) {
    g->matrix = a;
#if defined(HAVE_MPI)
    g->comm = cs_base_get_rank_step_comm(g->merge_stride);
    MPI_Comm_size(g->comm, &(g->n_ranks));
#endif
  }
  else {
    g->_matrix = cs_matrix_create_by_local_restrict(a);
    g->matrix = g->_matrix;
#if defined(HAVE_MPI)
    g->comm = cs_base_get_rank_step_comm(g->merge_stride);
    MPI_Comm_size(g->comm, &(g->n_ranks));
#endif
  }

  g->level = 0;
  g->symmetric = cs_matrix_is_symmetric(g->matrix);
  g->use_faces = false;

  const cs_lnum_t db_size = cs_matrix_get_diag_block_size(g->matrix);
  const cs_lnum_t eb_size = cs_matrix_get_extra_diag_block_size(g->matrix);

  g->db_size = db_size;
  g->eb_size = eb_size;

  g->n_rows = cs_matrix_get_n_rows(g->matrix);
  g->n_cols_ext = cs_matrix_get_n_columns(g->matrix);
  g->halo = cs_matrix_get_halo(g->matrix);

  g->n_g_rows = g->n_rows;

#if defined(HAVE_MPI)
  if (g->halo != NULL && g->comm != MPI_COMM_NULL) {
    cs_gnum_t _g_n_rows = g->n_rows;
    MPI_Allreduce(&_g_n_rows, &(g->n_g_rows), 1, CS_MPI_GNUM, MPI_SUM, g->comm);
  }
#endif

  return g;
}

/*----------------------------------------------------------------------------
 * Destroy a grid structure.
 *
 * parameters:
 *   grid <-> Pointer to grid structure pointer
 *----------------------------------------------------------------------------*/

void
cs_grid_destroy(cs_grid_t **grid)
{
  if (grid != NULL && *grid != NULL) {

    cs_grid_t *g = *grid;
    cs_grid_free_quantities(g);

    BFT_FREE(g->_face_cell);

    CS_FREE_HD(g->coarse_row);

    if (g->_halo != NULL)
      cs_halo_destroy(&(g->_halo));

    BFT_FREE(g->_da);
    BFT_FREE(g->_xa);

    cs_matrix_destroy(&(g->_matrix));
    cs_matrix_structure_destroy(&(g->matrix_struct));

    CS_FREE(g->cell_face);
    CS_FREE(g->cell_face_sgn);

#if defined(HAVE_MPI)
    BFT_FREE(g->merge_cell_idx);
#endif

    BFT_FREE(*grid);
  }
}

/*----------------------------------------------------------------------------
 * Free a grid structure's associated quantities.
 *
 * The quantities required to compute a coarser grid with relaxation from a
 * given grid are not needed after that stage, so may be freed.
 *
 * parameters:
 *   g <-> Pointer to grid structure
 *----------------------------------------------------------------------------*/

void
cs_grid_free_quantities(cs_grid_t  *g)
{
  assert(g != NULL);

  if (cs_matrix_get_type(g->matrix) != CS_MATRIX_NATIVE) {
    BFT_FREE(g->_face_cell);
    g->face_cell = NULL;
    BFT_FREE(g->_xa);
    g->xa = NULL;
    if (cs_matrix_get_type(g->matrix) == CS_MATRIX_CSR) {
      BFT_FREE(g->_da);
      g->da = NULL;
    }
  }

  BFT_FREE(g->coarse_face);

  BFT_FREE(g->_cell_cen);
  BFT_FREE(g->_cell_vol);
  BFT_FREE(g->_face_normal);

  BFT_FREE(g->xa_conv);
  BFT_FREE(g->xa_diff);
  BFT_FREE(g->_xa0);
  BFT_FREE(g->xa0_diff);

  BFT_FREE(g->xa0ij);
}

/*----------------------------------------------------------------------------
 * Get grid information.
 *
 * parameters:
 *   g          <-- Grid structure
 *   level      --> Level in multigrid hierarchy (or NULL)
 *   symmetric  --> Symmetric matrix coefficients indicator (or NULL)
 *   db_size    --> Size of the diagonal block (or NULL)
 *   eb_size    --> Size of the extra diagonal block (or NULL)
 *   n_ranks    --> number of ranks with data (or NULL)
 *   n_rows     --> Number of local rows (or NULL)
 *   n_cols_ext --> Number of columns including ghosts (or NULL)
 *   n_entries  --> Number of entries (or NULL)
 *   n_g_rows   --> Number of global rows (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 bool             *symmetric,
                 cs_lnum_t        *db_size,
                 cs_lnum_t        *eb_size,
                 int              *n_ranks,
                 cs_lnum_t        *n_rows,
                 cs_lnum_t        *n_cols_ext,
                 cs_lnum_t        *n_entries,
                 cs_gnum_t        *n_g_rows)
{
  assert(g != NULL);

  if (level != NULL)
    *level = g->level;

  if (symmetric != NULL)
    *symmetric = g->symmetric;

  if (db_size != NULL)
    *db_size = g->db_size;

  if (eb_size != NULL)
    *eb_size = g->eb_size;

  if (n_ranks != NULL) {
#if defined(HAVE_MPI)
    *n_ranks = g->n_ranks;
#else
    *n_ranks = 1;
#endif
  }

  if (n_rows != NULL)
    *n_rows = g->n_rows;
  if (n_cols_ext != NULL)
    *n_cols_ext = g->n_cols_ext;
  assert(g->n_rows <= g->n_cols_ext);
  if (n_entries != NULL) {
    if (g->matrix != NULL)
      *n_entries = cs_matrix_get_n_entries(g->matrix);
    else
      *n_entries = 0;
  }

  if (n_g_rows != NULL)
    *n_g_rows = g->n_g_rows;
}

/*----------------------------------------------------------------------------
 * Get memory allocation mode corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   memory allocation typt
 *----------------------------------------------------------------------------*/

cs_alloc_mode_t
cs_grid_get_alloc_mode(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->alloc_mode;
}

/*----------------------------------------------------------------------------
 * Get number of rows corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of rows of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_rows(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_rows;
}

/*----------------------------------------------------------------------------
 * Get number of extended (local + ghost) columns corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of extended rows of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cols_ext(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_cols_ext;
}

/*----------------------------------------------------------------------------
 * Get maximum number of extended (local + ghost) columns corresponding to
 * a grid, both with and without merging between ranks
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   maximum number of extended cells of grid structure, with or without
 *   merging
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cols_max(const cs_grid_t  *g)
{
  cs_lnum_t retval = 0;

  if (g != NULL)
    retval = CS_MAX(g->n_cols_ext, g->n_elts_r[1]);

  return retval;
}

/*----------------------------------------------------------------------------
 * Get global number of rows corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   global number of rows of grid structure
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_grid_get_n_g_rows(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_g_rows;
}

/*----------------------------------------------------------------------------
 * Get grid's associated matrix information.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   pointer to matrix structure
 *----------------------------------------------------------------------------*/

const cs_matrix_t *
cs_grid_get_matrix(const cs_grid_t  *g)
{
  const cs_matrix_t *m = NULL;

  assert(g != NULL);

  m = g->matrix;

  return m;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Get the MPI subcommunicator for a given grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
cs_grid_get_comm(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->comm;
}

/*----------------------------------------------------------------------------
 * Get the MPI subcommunicator for a given merge stride.
 *
 * parameters:
 *   parent       <-- parent MPI communicator
 *   merge_stride <-- associated merge stride
 *
 * returns:
 *   MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
cs_grid_get_comm_merge(MPI_Comm  parent,
                       int       merge_stride)
{
  MPI_Comm comm = MPI_COMM_NULL;

  if (parent != MPI_COMM_NULL) {
    int size;
    MPI_Comm_size(parent, &size);
    int rank_step = cs_glob_n_ranks / size;
    if (cs_glob_n_ranks %size > 0)
      rank_step += 1;
    rank_step *= merge_stride;
    comm = cs_base_get_rank_step_comm(rank_step);
  }

  return comm;
}

#endif

/*----------------------------------------------------------------------------
 * Create coarse grid from fine grid.
 *
 * parameters:
 *   f                          <-- Fine grid structure
 *   alloc_mode                 <-- Memory allocation mode
 *   coarsening_type            <-- Coarsening criteria type
 *   aggregation_limit          <-- Maximum allowed fine rows per coarse rows
 *   verbosity                  <-- Verbosity level
 *   merge_stride               <-- Associated merge stride
 *   merge_rows_mean_threshold  <-- mean number of rows under which
 *                                  merging should be applied
 *   merge_rows_glob_threshold  <-- global number of rows under which
 *                                  merging should be applied
 *   relaxation_parameter       <-- P0/P1 relaxation factor
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen(const cs_grid_t      *f,
                cs_alloc_mode_t       alloc_mode,
                cs_grid_coarsening_t  coarsening_type,
                int                   aggregation_limit,
                int                   verbosity,
                int                   merge_stride,
                int                   merge_rows_mean_threshold,
                cs_gnum_t             merge_rows_glob_threshold,
                double                relaxation_parameter)
{
  std::chrono::high_resolution_clock::time_point t_start;
  if (cs_glob_timer_kernels_flag > 0)
    t_start = std::chrono::high_resolution_clock::now();

  int recurse = 0;
  cs_lnum_t isym = 2;
  bool conv_diff = f->conv_diff;

  /* By default, always use MSR structure, as it usually provides the
     best performance, and is required for the hybrid Gauss-Seidel-Jacobi
     smoothers. In multithreaded case, we also prefer to use a matrix
     structure allowing threading without a specific renumbering, as
     structures are rebuilt often (so only CSR and MSR can be considered) */

  cs_matrix_type_t fine_matrix_type = cs_matrix_get_type(f->matrix);
  cs_matrix_type_t coarse_matrix_type = CS_MATRIX_MSR;

  cs_grid_t *c = NULL;

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t db_stride = db_size * db_size;

  assert(f != NULL);

  /* Initialization */

  c = _coarse_init(f, alloc_mode);
  c->n_rows = -1;

  if (f->symmetric == true)
    isym = 1;

  c->relaxation = relaxation_parameter;
  if (f->use_faces == false && c->relaxation > 0)
    c->relaxation = 0;

  /* Ensure default is available */

  if (coarsening_type == CS_GRID_COARSENING_DEFAULT) {
    if (f->use_faces) {
      if (f->conv_diff == false)
        coarsening_type = CS_GRID_COARSENING_SPD_DX;
      else
        coarsening_type = CS_GRID_COARSENING_CONV_DIFF_DX;
    }
    else
      coarsening_type = CS_GRID_COARSENING_SPD_MX;
  }
  else if (coarsening_type == CS_GRID_COARSENING_SPD_DX) {
    /* closest altenative */
    if (f->face_cell == NULL)
      coarsening_type = CS_GRID_COARSENING_SPD_MX;
  }

  else if (coarsening_type == CS_GRID_COARSENING_SPD_PW) {
    /* closest altenative */
    if (fine_matrix_type == CS_MATRIX_NATIVE)
      coarsening_type = CS_GRID_COARSENING_SPD_MX;
  }

  else if (coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX) {
    /* closest altenative */
    if (f->face_cell == NULL)
      coarsening_type = CS_GRID_COARSENING_SPD_MX;
  }

  /* Determine fine->coarse cell connectivity (aggregation) */

  if (   coarsening_type == CS_GRID_COARSENING_SPD_DX
      || coarsening_type == CS_GRID_COARSENING_CONV_DIFF_DX) {
    /* Performance still seems slightly better for the legacy version
       of the "DX" aggregation algorithm in single-threaded mode
       compared to the MSR-based version. */
    if (   f->use_faces
        && (fine_matrix_type != CS_MATRIX_MSR || cs_glob_n_threads == 1)) {
      c->n_rows = _automatic_aggregation_fc(f,
                                            coarsening_type,
                                            aggregation_limit,
                                            c->relaxation,
                                            verbosity,
                                            c->coarse_row);
    }
    else if (fine_matrix_type == CS_MATRIX_MSR)
      c->n_rows = _automatic_aggregation_dx_msr(f,
                                                coarsening_type,
                                                aggregation_limit,
                                                c->relaxation,
                                                verbosity,
                                                c->coarse_row);
  }
  else if (coarsening_type == CS_GRID_COARSENING_SPD_MX) {
    switch (fine_matrix_type) {
    case CS_MATRIX_NATIVE:
      c->n_rows = _automatic_aggregation_mx_native(f,
                                                   aggregation_limit,
                                                   verbosity,
                                                   c->coarse_row);
      break;
    case CS_MATRIX_MSR:
      c->n_rows = _automatic_aggregation_mx_msr(f,
                                                aggregation_limit,
                                                verbosity,
                                                c->coarse_row);
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _("%s: selected coarsening type (%s)\n"
                  "not available for %s structured matrices."),
                __func__, cs_grid_coarsening_type_name[coarsening_type],
                _(cs_matrix_get_type_name(f->matrix)));
    }
  }
  else if (coarsening_type == CS_GRID_COARSENING_SPD_PW) {
    switch (fine_matrix_type) {
    case CS_MATRIX_MSR:
      c->n_rows = _automatic_aggregation_pw_msr(f, verbosity, c->coarse_row);
      if (aggregation_limit > 2)
        recurse = 2;
      break;
    default:
      bft_error(__FILE__, __LINE__, 0,
                _("%s: selected coarsening type (%s)\n"
                  "not available for %s structured matrices."),
                __func__, cs_grid_coarsening_type_name[coarsening_type],
                _(cs_matrix_get_type_name(f->matrix)));
    }
  }

  _coarsen(f, c);

  if (verbosity > 3)
    _aggregation_stats_log(f, c, verbosity);

  BFT_MALLOC(c->_da, c->n_cols_ext * db_stride, cs_real_t);
  c->da = c->_da;

  BFT_MALLOC(c->_xa, c->n_faces*isym, cs_real_t);
  c->xa = c->_xa;

  if (  (fine_matrix_type == CS_MATRIX_NATIVE || f->use_faces)
      && c->relaxation > 0) {

    /* Allocate permanent arrays in coarse grid */

    BFT_MALLOC(c->_cell_cen, c->n_cols_ext*3, cs_real_t);
    c->cell_cen = c->_cell_cen;

    BFT_MALLOC(c->_cell_vol, c->n_cols_ext, cs_real_t);
    c->cell_vol = c->_cell_vol;

    BFT_MALLOC(c->_face_normal, c->n_faces*3, cs_real_t);
    c->face_normal = c->_face_normal;

    if (conv_diff) {
      BFT_MALLOC(c->xa_conv, c->n_faces*2, cs_real_t);
      BFT_MALLOC(c->xa_diff, c->n_faces, cs_real_t);
    }

    /* We could have xa0 point to xa if symmetric, but this would require
       caution in CRSTGR to avoid overwriting. */

    BFT_MALLOC(c->_xa0, c->n_faces*isym, cs_real_t);
    c->xa0 = c->_xa0;

    if (conv_diff) {
      BFT_MALLOC(c->xa0_diff, c->n_faces, cs_real_t);
    }

    BFT_MALLOC(c->xa0ij, c->n_faces*3, cs_real_t);

    /* Matrix-related data */

    _compute_coarse_cell_quantities(f, c);

    /* Synchronize grid's geometric quantities */

    if (c->halo != NULL) {

      cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_cell_cen, 3);
      if (c->halo->n_transforms > 0)
        cs_halo_perio_sync_coords(c->halo, CS_HALO_STANDARD, c->_cell_cen);

      cs_halo_sync_var(c->halo, CS_HALO_STANDARD, c->_cell_vol);

    }

    /* New threaded variant, in progress */
#if 0
    if (   fine_matrix_type == CS_MATRIX_MSR
        && coarse_matrix_type == CS_MATRIX_MSR)
      _compute_coarse_quantities_msr_with_faces(f, c);
#endif

  }

  if (fine_matrix_type == CS_MATRIX_MSR && c->relaxation <= 0) {

   _compute_coarse_quantities_msr(f, c);

  }

  else if (f->use_faces) {

    if (conv_diff)
      _compute_coarse_quantities_conv_diff(f, c, verbosity);
    else
      _compute_coarse_quantities_native(f, c, verbosity);

    /* Synchronize matrix's geometric quantities */

    if (c->halo != NULL)
      cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_da, db_stride);

    _matrix_from_native(coarse_matrix_type, c);

  }

  /* Apply tuning if needed */

  if (_grid_tune_max_level > 0 && c->matrix != nullptr && recurse < 2) {

    cs_matrix_fill_type_t mft
      = cs_matrix_get_fill_type(f->symmetric,
                                f->db_size,
                                f->eb_size);

    if (_grid_tune_max_level > f->level) {
      int k = CS_MATRIX_N_FILL_TYPES*(f->level) + mft;
      cs_matrix_variant_t *coarse_mv = _grid_tune_variant[k];

      /* Create tuned variant upon first pass for this level and
         fill type */

      if  (   coarse_mv == NULL
           && _grid_tune_max_fill_level[mft] > f->level) {

        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("\n"
                        "Tuning for coarse matrices of level %d and type: %s\n"
                        "==========================\n"),
                      f->level + 1, cs_matrix_fill_type_name[mft]);

        int n_min_products = cs_matrix_get_tuning_runs();

        coarse_mv = cs_matrix_variant_tuned(c->matrix,
                                            1,
                                            n_min_products);

        _grid_tune_variant[k] = coarse_mv;

        if  (_grid_tune_max_fill_level[mft] == f->level + 1) {
          cs_log_printf(CS_LOG_PERFORMANCE, "\n");
          cs_log_separator(CS_LOG_PERFORMANCE);
        }

      }

      if (coarse_mv != NULL)
        cs_matrix_variant_apply_tuned(c->_matrix, coarse_mv);
    }

  }

  if (c->matrix == NULL) {
    assert(c->n_rows == 0);
    _build_coarse_matrix_null(c, coarse_matrix_type);
  }

  /* Recurse if necessary */

  if (recurse > 1) {

    /* Build coarser grid from coarse grid */
    cs_grid_t *cc = cs_grid_coarsen(c,
                                    c->alloc_mode,
                                    coarsening_type,
                                    aggregation_limit / recurse,
                                    verbosity,
                                    0, // merge_stride,
                                    0, // merge_rows_mean_threshold,
                                    0, //merge_rows_glob_threshold,
                                    relaxation_parameter);

    /* Project coarsening */

    _project_coarse_row_to_parent(cc);
    CS_FREE_HD(cc->coarse_row);
    cc->coarse_row = c->coarse_row;
    c->coarse_row = NULL;

    if (c->use_faces) {
      BFT_FREE(cc->coarse_face);
      BFT_FREE(cc->_face_cell);
      _coarsen_faces(f,
                     cc->coarse_row,
                     cc->n_rows,
                     &(cc->n_faces),
                     &(cc->coarse_face),
                     &(cc->_face_cell));
      cc->face_cell = (const cs_lnum_2_t  *)(cc->_face_cell);
    }

    /* Keep coarsest grid only */

    cc->level -=1;
    cc->parent = f;

    assert(cc->parent == c->parent);

    cs_grid_destroy(&c);
    c = cc;
  }

  /* Merge grids if we are below the threshold */

#if defined(HAVE_MPI)

  if (merge_stride > 1 && c->n_ranks > 1) {
    cs_gnum_t  _n_ranks = c->n_ranks;
    cs_gnum_t  _n_mean_g_rows = c->n_g_rows / _n_ranks;
    if (   _n_mean_g_rows < (cs_gnum_t)merge_rows_mean_threshold
        || c->n_g_rows < merge_rows_glob_threshold) {
      if (c->_xa == NULL && c->n_faces > 0)
        _native_from_msr(c);
      _merge_grids(c, merge_stride, verbosity);
      if (c->_matrix != nullptr) {
        cs_matrix_type_t cm_type = cs_matrix_get_type(c->matrix);
        cs_matrix_destroy(&(c->_matrix));
        cs_matrix_structure_destroy(&(c->matrix_struct));
        _matrix_from_native(cm_type, c);
      }
    }
  }

#endif

  if (f->use_faces)
    cs_matrix_set_mesh_association(c->_matrix,
                                   NULL,
                                   NULL,
                                   NULL,
                                   (const cs_real_3_t *)c->cell_cen,
                                   (const cs_real_t *)c->cell_vol,
                                   (const cs_real_3_t *)c->_face_normal);

  /* Optional verification */

  if (verbosity > 3) {
    if (f->level == 0)
      _verify_matrix(f);
    _verify_matrix(c);
  }

  if (cs_glob_timer_kernels_flag > 0) {
    std::chrono::high_resolution_clock::time_point
      t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds elapsed
      = std::chrono::duration_cast
          <std::chrono::microseconds>(t_stop - t_start);
    printf("%d: %s (level %d)", cs_glob_rank_id, __func__, f->level);
    printf(", total = %ld\n", elapsed.count());
  }

  /* Return new (coarse) grid */

  return c;
}

/*----------------------------------------------------------------------------
 * Create coarse grid with only one row per rank from fine grid.
 *
 * parameters:
 *   f            <-- Fine grid structure
 *   alloc_mode   <-- Memory allocation mode
 *   merge_stride <-- Associated merge stride
 *   verbosity    <-- Verbosity level
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen_to_single(const cs_grid_t  *f,
                          cs_alloc_mode_t   alloc_mode,
                          int               merge_stride,
                          int               verbosity)
{
  cs_lnum_t isym = 2;

  /* By default, always use MSR structure, as it often seems to provide the
     best performance, and is required for the hybrid Gauss-Seidel-Jacobi
     smoothers. In multithreaded case, we also prefer to use a matrix
     structure allowing threading without a specific renumbering, as
     structures are rebuilt often (so only CSR and MSR can be considered) */

  cs_matrix_type_t fine_matrix_type = cs_matrix_get_type(f->matrix);

  cs_grid_t *c = NULL;

  const cs_lnum_t db_size = f->db_size;
  const cs_lnum_t db_stride = db_size * db_size;

  assert(f != NULL);

  /* Initialization */

  c = _coarse_init(f, alloc_mode);

  if (f->symmetric == true)
    isym = 1;

  c->relaxation = 0;

  /* All to a single row */

  for (cs_lnum_t i = 0; i < f->n_rows; i++)
    c->coarse_row[i] = 0;

  _coarsen(f, c);

  if (verbosity > 3)
    _aggregation_stats_log(f, c, verbosity);

  if (fine_matrix_type == CS_MATRIX_MSR) {
    _compute_coarse_quantities_msr(f, c);

#if defined(HAVE_MPI)
    if (c->n_ranks > 1 && merge_stride > 1) {
      _native_from_msr(c);
      _merge_grids(c, merge_stride, verbosity);
      _matrix_from_native(CS_MATRIX_MSR, c);
    }
#endif
  }

  else if (f->use_faces) {

    BFT_MALLOC(c->_da, c->n_cols_ext * db_stride, cs_real_t);
    c->da = c->_da;

    BFT_MALLOC(c->_xa, c->n_faces*isym, cs_real_t);
    c->xa = c->_xa;

    _compute_coarse_quantities_native(f, c, verbosity);

    /* Synchronize matrix's geometric quantities */

    if (c->halo != NULL)
      cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_da, db_stride);

    /* Merge grids if we are below the threshold */

#if defined(HAVE_MPI)
    if (c->n_ranks > 1 && merge_stride > 1)
      _merge_grids(c, merge_stride, verbosity);
#endif

    _matrix_from_native(CS_MATRIX_MSR, c);
  }

  c->matrix = c->_matrix;

  /* Optional verification */

  if (verbosity > 3) {
    if (f->level == 0)
      _verify_matrix(f);
    _verify_matrix(c);
  }

  /* Return new (coarse) grid */

  return c;
}

/*----------------------------------------------------------------------------
 * Project coarse grid row numbers to base grid.
 *
 * If a global coarse grid row number is larger than max_num, its
 * value modulo max_num is used.
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   max_num     <-- Values of c_row_num = global_num % max_num
 *   c_row_num   --> Global coarse row number (modulo max_num)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_row_num(const cs_grid_t  *g,
                        cs_lnum_t         n_base_rows,
                        int               max_num,
                        int               c_row_num[])
{
  cs_lnum_t ii = 0;
  cs_gnum_t base_shift = 1;
  cs_gnum_t _max_num = max_num;
  cs_lnum_t n_max_rows = 0;
  cs_lnum_t *tmp_num_1 = NULL, *tmp_num_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(c_row_num != NULL);

  /* Initialize array */

  n_max_rows = g->n_rows;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_rows > n_max_rows)
      n_max_rows = _g->n_rows;
  }

  BFT_MALLOC(tmp_num_1, n_max_rows, cs_lnum_t);

  /* Compute local base starting row number in parallel mode */

#if defined(HAVE_MPI)
  if (g->comm != MPI_COMM_NULL) {
    cs_gnum_t local_shift = g->n_rows;
    cs_gnum_t global_shift = 0;
    MPI_Scan(&local_shift, &global_shift, 1, CS_MPI_GNUM, MPI_SUM,
             g->comm);
    base_shift = 1 + global_shift - g->n_rows;
  }
#endif

  for (ii = 0; ii < g->n_rows; ii++)
    tmp_num_1[ii] = (cs_gnum_t)(ii + base_shift) % _max_num;

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_num_2, n_max_rows, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_max_rows; i++)
      tmp_num_2[i] = -1;        /* singleton */

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_rows = _g->parent->n_rows;

#if defined(HAVE_MPI)
      _scatter_row_num(_g, tmp_num_1);
#endif /* defined(HAVE_MPI) */

      for (ii = 0; ii < n_parent_rows; ii++) {
        cs_lnum_t ic = _g->coarse_row[ii];
        if (ic >= 0)
          tmp_num_2[ii] = tmp_num_1[ic];
      }

      for (ii = 0; ii < n_parent_rows; ii++)
        tmp_num_1[ii] = tmp_num_2[ii];

    }

    assert(_g->level == 0);
    assert(_g->n_rows == n_base_rows);

    /* Free temporary arrays */

    BFT_FREE(tmp_num_2);
  }

  memcpy(c_row_num, tmp_num_1, n_base_rows*sizeof(int));

  BFT_FREE(tmp_num_1);
}

/*----------------------------------------------------------------------------
 * Project coarse grid row rank to base grid.
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   f_row_rank  --> Global coarse row rank projected to fine rows
 *----------------------------------------------------------------------------*/

void
cs_grid_project_row_rank(const cs_grid_t  *g,
                         cs_lnum_t         n_base_rows,
                         int               f_row_rank[])
{
  cs_lnum_t ii;
  cs_lnum_t n_max_rows = 0;
  int *tmp_rank_1 = NULL, *tmp_rank_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(f_row_rank != NULL || g->n_rows == 0);

  n_max_rows = g->n_rows;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_rows > n_max_rows)
      n_max_rows = _g->n_rows;
  }

  BFT_MALLOC(tmp_rank_1, n_max_rows, int);

  for (ii = 0; ii < g->n_rows; ii++)
    tmp_rank_1[ii] = cs_glob_rank_id;

  /* Project to finer levels */

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_rank_2, n_max_rows, int);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_rows = _g->parent->n_rows;

      _prolong_row_int(_g,
                       _g->parent,
                       tmp_rank_1,
                       tmp_rank_2);

      for (ii = 0; ii < n_parent_rows; ii++)
        tmp_rank_1[ii] = tmp_rank_2[ii];

    }

    assert(_g->level == 0);
    assert(_g->n_rows == n_base_rows);

    /* Free temporary arrays */

    BFT_FREE(tmp_rank_2);
  }

  memcpy(f_row_rank, tmp_rank_1, n_base_rows*sizeof(int));

  BFT_FREE(tmp_rank_1);
}

/*----------------------------------------------------------------------------
 * Project variable from coarse grid to base grid
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   c_var       <-- Variable on coarse grid
 *   f_var       --> Variable projected to fine grid
 *----------------------------------------------------------------------------*/

void
cs_grid_project_var(const cs_grid_t  *g,
                    cs_lnum_t         n_base_rows,
                    const cs_real_t   c_var[],
                    cs_real_t         f_var[])
{
  cs_lnum_t ii;
  int i;
  cs_lnum_t n_max_rows = 0;
  cs_real_t *tmp_var_1 = NULL, *tmp_var_2 = NULL;
  const cs_grid_t *_g = g;

  const cs_lnum_t db_size = g->db_size;

  assert(g != NULL);
  assert(c_var != NULL || g->n_rows == 0);
  assert(f_var != NULL);

  n_max_rows = g->n_rows;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_rows > n_max_rows)
      n_max_rows = _g->n_rows;
  }

  cs_alloc_mode_t amode = cs_check_device_ptr(f_var);
  CS_MALLOC_HD(tmp_var_1, n_max_rows*db_size, cs_real_t, amode);
  memcpy(tmp_var_1, c_var, g->n_rows*db_size*sizeof(cs_real_t));

  /* Project to finer levels */

  if (g->level > 0) {

    cs_dispatch_context ctx;
    if (amode == CS_ALLOC_HOST)
      ctx.set_use_gpu(false);

    /* Allocate temporary arrays */

    CS_MALLOC_HD(tmp_var_2, n_max_rows*db_size, cs_real_t, amode);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_rows = _g->parent->n_rows;

      cs_grid_prolong_row_var(ctx,
                              _g,
                              _g->parent,
                              false,
                              tmp_var_1,
                              tmp_var_2);

      for (ii = 0; ii < n_parent_rows; ii++)
        for (i = 0; i < db_size; i++)
          tmp_var_1[ii*db_size+i] = tmp_var_2[ii*db_size+i];

    }

    assert(_g->level == 0);
    assert(_g->n_rows == n_base_rows);

    /* Free temporary arrays */

    BFT_FREE(tmp_var_2);
  }

  memcpy(f_var, tmp_var_1, n_base_rows*db_size*sizeof(cs_real_t));

  BFT_FREE(tmp_var_1);
}

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric and project it to base grid
 *
 * parameters:
 *   g           <-- Grid structure
 *   n_base_rows <-- Number of rows in base grid
 *   diag_dom    --> Diagonal dominance metric (on fine grid)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_diag_dom(const cs_grid_t  *g,
                         cs_lnum_t         n_base_rows,
                         cs_real_t         diag_dom[])
{
  cs_real_t *dd = NULL;
  const cs_lnum_t db_size = g->db_size;
  const cs_lnum_t db_stride = db_size * db_size;

  assert(g != NULL);
  assert(diag_dom != NULL);

  if (g->level == 0)
    dd = diag_dom;
  else
    BFT_MALLOC(dd, g->n_cols_ext*db_stride, cs_real_t);

  /* Compute coarse diagonal dominance */

  cs_matrix_diag_dominance(g->matrix, dd);

  /* Now project to finer levels */

  if (dd != diag_dom) {
    cs_grid_project_var(g, n_base_rows, dd, diag_dom);
    BFT_FREE(dd);
  }
}

/*----------------------------------------------------------------------------
 * Finalize global info related to multigrid solvers
 *----------------------------------------------------------------------------*/

void
cs_grid_finalize(void)
{
  if (_grid_tune_max_level > 0) {

    for (int i = 0; i < _grid_tune_max_level; i++) {
      for (int j = 0; j < CS_MATRIX_N_FILL_TYPES; j++) {
        int k = CS_MATRIX_N_FILL_TYPES*i + j;
        if (_grid_tune_variant[k] != NULL)
          cs_matrix_variant_destroy(&(_grid_tune_variant[k]));
      }
    }

    BFT_FREE(_grid_tune_variant);
    BFT_FREE(_grid_tune_max_fill_level);

    _grid_tune_max_level = 0;
  }
}

/*----------------------------------------------------------------------------
 * Dump grid structure
 *
 * parameters:
 *   g <-- grid structure that should be dumped
 *----------------------------------------------------------------------------*/

void
cs_grid_dump(const cs_grid_t  *g)
{
  cs_lnum_t  i;

  if (g == NULL) {
    bft_printf("\n\n  grid: null\n");
    return;
  }
  bft_printf("\n"
             "  grid:          %p\n"
             "  level:         %d (parent: %p)\n"
             "  n_rows:        %d\n"
             "  n_cols_ext:    %d\n"
             "  n_faces:       %d\n"
             "  n_g_cells:     %d\n"
             "  n_elts_r:      [%d, %d]\n",
             (const void *)g, g->level, (const void *)(g->parent),
             (int)(g->n_rows), (int)(g->n_cols_ext),
             (int)(g->n_faces), (int)(g->n_g_rows),
             (int)(g->n_elts_r[0]), (int)(g->n_elts_r[1]));

#if defined(HAVE_MPI)

  bft_printf("\n"
             "  merge_sub_root:     %d\n"
             "  merge_sub_rank:     %d\n"
             "  merge_sub_size:     %d\n"
             "  merge_stride:       %d\n"
             "  next_merge_stride:  %d\n"
             "  n_ranks:            %d\n",
             g->merge_sub_root, g->merge_sub_rank, g->merge_sub_size,
             g->merge_stride, g->next_merge_stride, g->n_ranks);

  if (g->merge_cell_idx != NULL) {
    bft_printf("  merge_cell_idx\n");
    for (i = 0; i < g->merge_sub_size + 1; i++)
      bft_printf("    %ld: %ld\n", (long)i, (long)g->merge_cell_idx[i]);
  }

#endif
  bft_printf("\n"
             "  face_cell:      %p\n"
             "  _face_cell:     %p\n"
             "  coarse_row:     %p\n"
             "  coarse_face:    %p\n"
             "  halo:           %p\n",
             (const void *)g->face_cell, (const void *)g->_face_cell,
             (const void *)g->coarse_row, (const void *)g->coarse_face,
             (const void *)g->halo);

  if (g->face_cell != NULL) {
    bft_printf("\n"
               "  face -> cell connectivity;\n");
    for (i = 0; i < g->n_faces; i++)
      bft_printf("    %ld : %ld, %ld\n", (long)(i+1),
                 (long)(g->face_cell[i][0]), (long)(g->face_cell[i][1]));
  }

  if (g->coarse_row != NULL && g->parent != NULL) {
    bft_printf("\n"
               "  coarse_row;\n");
    for (i = 0; i < g->parent->n_rows; i++)
      bft_printf("    %ld : %ld\n",
                 (long)(i+1), (long)(g->coarse_row[i]));
  }

  if (g->coarse_face != NULL && g->parent != NULL) {
    bft_printf("\n"
               "  coarse_face;\n");
    for (i = 0; i < g->parent->n_faces; i++)
      bft_printf("    %ld : %ld\n",
                 (long)(i+1), (long)(g->coarse_face[i]));
  }

  cs_halo_dump(g->halo, 1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set matrix tuning behavior for multigrid coarse meshes.
 *
 * The finest mesh (level 0) is handled by the default tuning options,
 * so only coarser meshes are considered here.
 *
 * \param[in]  fill_type  associated matrix fill type
 * \param[in]  max_level  maximum level for which tuning is active
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_set_matrix_tuning(cs_matrix_fill_type_t  fill_type,
                          int                    max_level)
{
  if (_grid_tune_max_level < max_level) {

    if (_grid_tune_max_level == 0) {
      BFT_MALLOC(_grid_tune_max_fill_level, CS_MATRIX_N_FILL_TYPES, int);
      for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++)
        _grid_tune_max_fill_level[i] = 0;
    }

    BFT_REALLOC(_grid_tune_variant,
                CS_MATRIX_N_FILL_TYPES*max_level, cs_matrix_variant_t *);

    for (int i = _grid_tune_max_level; i < max_level; i++) {
      for (int j = 0; j < CS_MATRIX_N_FILL_TYPES; j++) {
        _grid_tune_variant[CS_MATRIX_N_FILL_TYPES*i + j] = NULL;
      }
    }

    _grid_tune_max_level = max_level;
  }

  _grid_tune_max_fill_level[fill_type] = max_level;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

/*----------------------------------------------------------------------------
 * Compute coarse row variable values from fine row values
 *
 * parameters:
 *   ctx     <-> Reference to dispatch context
 *   f       <-- Fine grid structure
 *   c       <-- Fine grid structure
 *   f_var   <-- Variable defined on fine grid rows
 *   c_var   --> Variable defined on coarse grid rows
 *----------------------------------------------------------------------------*/

void
cs_grid_restrict_row_var(cs_dispatch_context  &ctx,
                         const cs_grid_t      *f,
                         const cs_grid_t      *c,
                         const cs_real_t      *f_var,
                         cs_real_t            *c_var)
{
  cs_lnum_t f_n_rows = f->n_rows;
  cs_lnum_t c_n_cols_ext = c->n_elts_r[1];

  const cs_lnum_t *coarse_row;
  const cs_lnum_t db_size = f->db_size;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_row != NULL || f_n_rows == 0);
  assert(f_var != NULL || f_n_rows == 0);
  assert(c_var != NULL || c_n_cols_ext == 0);

  /* Set coarse values */

  coarse_row = c->coarse_row;

  cs_lnum_t _c_n_cols_ext = c_n_cols_ext*db_size;

  cs_dispatch_sum_type_t sum_type = CS_DISPATCH_SUM_SIMPLE;
  if (ctx.use_gpu())
    sum_type = CS_DISPATCH_SUM_ATOMIC;
#if defined(HAVE_OPENMP)
  else if (ctx.n_min_for_cpu_threads() <= f_n_rows)
    sum_type = CS_DISPATCH_SUM_ATOMIC;
#endif

  ctx.parallel_for(_c_n_cols_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
    c_var[ii] = 0.;
  });

  if (db_size == 1) {
    ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
      cs_lnum_t i = coarse_row[ii];
      if (i >= 0)
        cs_dispatch_sum(&(c_var[i]), f_var[ii], sum_type);
    });
  }
  else {
    switch(db_size) {
    case 3:
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t i = coarse_row[ii];
        if (i >= 0) {
          cs_real_t tmp[3];
          for (cs_lnum_t j = 0; j < 3; j++)
            tmp[j] = f_var[ii*db_size+j];
          cs_dispatch_sum<3>(&(c_var[i*db_size]), tmp, sum_type);
        }
      });
      break;
    default:
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t i = coarse_row[ii];
        if (i >= 0) {
          for (cs_lnum_t j = 0; j < db_size; j++)
            cs_dispatch_sum(&(c_var[i*db_size+j]),
                            f_var[ii*db_size+j],
                            sum_type);
        }
      });
    }
  }
  ctx.wait();

#if defined(HAVE_MPI)

  /* If grid merging has taken place, gather coarse data */

  if (c->merge_sub_size > 1) {

    MPI_Comm  comm = cs_glob_mpi_comm;
    static const int tag = 'r'+'e'+'s'+'t'+'r'+'i'+'c'+'t';

    /* Append data */

    if (c->merge_sub_rank == 0) {
      int rank_id;
      MPI_Status status;
      assert(cs_glob_rank_id == c->merge_sub_root);
      for (rank_id = 1; rank_id < c->merge_sub_size; rank_id++) {
        cs_lnum_t n_recv = (  c->merge_cell_idx[rank_id+1]
                            - c->merge_cell_idx[rank_id]);
        int dist_rank = c->merge_sub_root + c->merge_stride*rank_id;
        MPI_Recv(c_var + c->merge_cell_idx[rank_id]*db_size,
                 n_recv*db_size, CS_MPI_REAL, dist_rank, tag, comm, &status);
      }
    }
    else
      MPI_Send(c_var, c->n_elts_r[0]*db_size, CS_MPI_REAL,
               c->merge_sub_root, tag, comm);
  }

#endif /* defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------
 * Compute fine row variable values from coarse row values
 *
 * parameters:
 *   ctx       <-> Reference to dispatch context
 *   c         <-- Fine grid structure
 *   f         <-- Fine grid structure
 *   increment <-- if true, add value to f_var; otherwise, overwrite it
 *   c_var     <-- Variable defined on coarse grid rows
 *   f_var     <-> Variable defined on fine grid rows
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_row_var(cs_dispatch_context  &ctx,
                        const cs_grid_t      *c,
                        const cs_grid_t      *f,
                        bool                  increment,
                        cs_real_t            *c_var,
                        cs_real_t            *f_var)
{
  const cs_lnum_t *coarse_row;
  const cs_real_t *_c_var = c_var;

  const cs_lnum_t db_size = f->db_size;

  cs_lnum_t f_n_rows = f->n_rows;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_row != NULL || f_n_rows == 0);
  assert(f_var != NULL || f_n_rows == 0);
  assert(c_var != NULL || c->n_cols_ext == 0);

#if defined(HAVE_MPI)

  /* If grid merging has taken place, scatter coarse data */

  if (c->merge_sub_size > 1) {

    MPI_Comm  comm = cs_glob_mpi_comm;
    static const int tag = 'p'+'r'+'o'+'l'+'o'+'n'+'g';

    /* Append data */

    if (c->merge_sub_rank == 0) {
      int rank_id;
      assert(cs_glob_rank_id == c->merge_sub_root);
      for (rank_id = 1; rank_id < c->merge_sub_size; rank_id++) {
        cs_lnum_t n_send = (  c->merge_cell_idx[rank_id+1]
                            - c->merge_cell_idx[rank_id]);
        int dist_rank = c->merge_sub_root + c->merge_stride*rank_id;
        MPI_Send(c_var + c->merge_cell_idx[rank_id]*db_size,
                 n_send*db_size, CS_MPI_REAL, dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(c_var, c->n_elts_r[0]*db_size, CS_MPI_REAL,
               c->merge_sub_root, tag, comm, &status);
    }
  }

#endif /* defined(HAVE_MPI) */

  /* Set fine values (possible penalization at first level) */

  coarse_row = c->coarse_row;

  /* Increment mode */

  if (increment) {

    if (db_size == 1) {
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t ic = coarse_row[ii];
        if (ic >= 0) {
          f_var[ii] += _c_var[ic];
        }
      });
    }
    else {
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t ic = coarse_row[ii];
        if (ic >= 0) {
          for (cs_lnum_t i = 0; i < db_size; i++)
            f_var[ii*db_size+i] += _c_var[ic*db_size+i];
        }
      });
    }

  }

  /* Overwrite mode */

  else {
    if (db_size == 1) {
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t ic = coarse_row[ii];
        if (ic >= 0) {
          f_var[ii] = _c_var[ic];
        }
        else
          f_var[ii] = 0;
      });
    }
    else {
      ctx.parallel_for(f_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t ii) {
        cs_lnum_t ic = coarse_row[ii];
        if (ic >= 0) {
          for (cs_lnum_t i = 0; i < db_size; i++)
            f_var[ii*db_size+i] = _c_var[ic*db_size+i];
        }
        else {
          for (cs_lnum_t i = 0; i < db_size; i++)
            f_var[ii*db_size+i] = 0;
        }
      });
    }

  }

  ctx.wait();
}

/*----------------------------------------------------------------------------*/

