/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

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
#include "cs_order.h"
#include "cs_prototypes.h"
#include "cs_sles.h"

#include "fvm_defs.h"
#include "fvm_hilbert.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_grid.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
  bool                conv_diff;    /* true if convection/diffusion case,
                                       false otherwhise */
  bool                symmetric;    /* Symmetric matrix coefficients
                                       indicator */
  int                 diag_block_size[4]; /* Block sizes for diagonal, or NULL */
  int                 extra_diag_block_size[4]; /* Block sizes for
                                                   extra diagonal, or NULL */

  cs_lnum_t           n_cells;      /* Local number of cells */
  cs_lnum_t           n_cells_ext;  /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */
  cs_lnum_t           n_faces;      /* Local number of faces */
  cs_gnum_t           n_g_cells;    /* Global number of cells */

  cs_lnum_t           n_cells_r[2]; /* Size of array used for restriction
                                       operations ({ncells, n_cells_ext} when
                                       no grid merging has taken place) */

  /* Grid hierarchy information */

  const struct _cs_grid_t  *parent; /* Pointer to parent (finer) grid */

  /* Connectivity information */

  const cs_lnum_2_t  *face_cell;    /* Face -> cells connectivity (1 to n) */
  cs_lnum_2_t        *_face_cell;   /* Face -> cells connectivity
                                       (private array) */

  /* Restriction from parent to current level */

  cs_lnum_t          *coarse_cell;  /* Fine -> coarse cell connectivity
                                       (1 to n); size: parent n_cells_ext */
  cs_lnum_t          *coarse_face;  /* Fine -> coarse face connectivity
                                       (1 to n, signed:
                                       = 0 fine face inside coarse cell
                                       > 0 orientation same as parent
                                       < 0 orientation opposite as parent);
                                       size: parent n_faces */

  /* Geometric data */

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
  const cs_real_t  *da_conv;        /* Diagonal (shared) */
  cs_real_t        *_da_conv;       /* Diagonal (private) */
  const cs_real_t  *da_diff;        /* Diagonal (shared) */
  cs_real_t        *_da_diff;       /* Diagonal (private) */

  const cs_real_t  *xa;             /* Extra-diagonal (shared) */
  cs_real_t        *_xa;            /* Extra-diagonal (private) */
  const cs_real_t  *xa_conv;        /* Extra-diagonal (shared) */
  cs_real_t        *_xa_conv;       /* Extra-diagonal (private) */
  const cs_real_t  *xa_diff;        /* Extra-diagonal (shared) */
  cs_real_t        *_xa_diff;       /* Extra-diagonal (private) */

  const cs_real_t  *xa0;            /* Symmetrized extra-diagonal (shared) */
  cs_real_t        *_xa0;           /* Symmetrized extra-diagonal (private) */
  const cs_real_t  *xa0_diff;       /* Symmetrized extra-diagonal (shared) */
  cs_real_t        *_xa0_diff;      /* Symmetrized extra-diagonal (private) */

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
  int               comm_id;           /* Associated communicator
                                          (owner when merge_idx != NULL) */

  /* Distribution information */

#endif
};

/*============================================================================
 *  Global variables
 *============================================================================*/

#if defined(HAVE_MPI)

static cs_gnum_t   _grid_merge_mean_threshold = 300;
static cs_gnum_t   _grid_merge_glob_threshold = 500;
static int         _grid_merge_min_ranks = 1;
static int         _grid_merge_stride = 1;

static int        _n_grid_comms = 0;
static int       *_grid_ranks = NULL;
static MPI_Comm  *_grid_comm = NULL;

#endif /* defined(HAVE_MPI) */

/* Names for coarsening options */

const char *cs_grid_coarsening_type_name[]
  = {N_("algebraic, natural face traversal"),
     N_("algebraic, face traveral by criteria"),
     N_("algebraic, face traversal by Hilbert SFC")};

/* Select tuning options */

static int _grid_tune_max_level = 0;
static int *_grid_tune_max_fill_level = NULL;
static cs_matrix_variant_t **_grid_tune_variant = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

  g->diag_block_size[0] = 1;
  g->diag_block_size[1] = 1;
  g->diag_block_size[2] = 1;
  g->diag_block_size[3] = 1;

  g->extra_diag_block_size[0] = 1;
  g->extra_diag_block_size[1] = 1;
  g->extra_diag_block_size[2] = 1;
  g->extra_diag_block_size[3] = 1;

  g->level = 0;

  g->n_cells = 0;
  g->n_cells_ext = 0;
  g->n_faces = 0;
  g->n_g_cells = 0;

  g->n_cells_r[0] = 0;
  g->n_cells_r[1] = 0;

  g->parent = NULL;
  g->conv_diff = false;

  g->face_cell = NULL;
  g->_face_cell = NULL;

  g->coarse_cell = NULL;
  g->coarse_face = NULL;

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
  g->da_conv = NULL;
  g->_da_conv = NULL;
  g->da_diff = NULL;
  g->_da_diff = NULL;
  g->xa = NULL;
  g->_xa = NULL;
  g->xa_conv = NULL;
  g->_xa_conv = NULL;
  g->xa_diff = NULL;
  g->_xa_diff = NULL;
  g->xa0 = NULL;
  g->_xa0 = NULL;
  g->xa0_diff = NULL;
  g->_xa0_diff = NULL;

  g->xa0ij = NULL;

  g->matrix = NULL;

#if defined(HAVE_MPI)

  g->merge_sub_root = 0;
  g->merge_sub_rank = 0;
  g->merge_sub_size = 1;
  g->merge_stride = 0;
  g->next_merge_stride = 1;

  g->merge_cell_idx = NULL;

  g->n_ranks = cs_glob_n_ranks;
  g->comm_id = 0;

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
 *   f <-- Pointer to fine (parent) grid structure
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

static cs_grid_t *
_coarse_init(const cs_grid_t *f)
{
  cs_lnum_t ii;
  int i;
  cs_grid_t *c = NULL;

  c = _create_grid();

  c->parent = f;

  c->level = f->level + 1;
  c->symmetric = f->symmetric;
  c->conv_diff = f->conv_diff;
  for (i = 0; i < 4; i++)
    c->diag_block_size[i] = f->diag_block_size[i];

  for (i = 0; i < 4; i++)
    c->extra_diag_block_size[i] = f->extra_diag_block_size[i];

  BFT_MALLOC(c->coarse_cell, f->n_cells_ext, cs_lnum_t);

# pragma omp parallel for if(f->n_cells_ext > CS_THR_MIN)
  for (ii = 0; ii < f->n_cells_ext; ii++)
    c->coarse_cell[ii] = 0;

#if defined(HAVE_MPI)
  c->merge_stride = f->merge_stride;
  c->next_merge_stride = f->next_merge_stride;
  c->n_ranks = f->n_ranks;
  c->comm_id = f->comm_id;
#endif

  return c;
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
 *   coarse_cell      <-- Fine -> coarse cell id (1 to n numbering)
 *                        size: fine->n_cells_ext
 *   n_coarse_cells   <-- Number of coarse cells
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
               const cs_lnum_t    *restrict coarse_cell,
               cs_lnum_t           n_coarse_cells,
               cs_lnum_t          *n_coarse_faces,
               cs_lnum_t         **coarse_face,
               cs_lnum_2_t       **coarse_face_cell)
{
  cs_lnum_t  ii, jj, face_id, connect_size;

  cs_lnum_t  *restrict c_cell_cell_cnt = NULL;
  cs_lnum_t  *restrict c_cell_cell_idx = NULL;
  cs_lnum_t  *restrict c_cell_cell_id = NULL;
  cs_lnum_t  *restrict c_cell_cell_face = NULL;

  cs_lnum_t    *restrict _coarse_face = NULL;
  cs_lnum_2_t  *restrict _c_face_cell = NULL;

  cs_lnum_t   c_n_faces = 0;

  const cs_lnum_t c_n_cells = n_coarse_cells;
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

    ii = coarse_cell[f_face_cell[face_id][0]] - 1;
    jj = coarse_cell[f_face_cell[face_id][1]] - 1;

    if (ii < jj)
      c_cell_cell_idx[ii+1] += 1;
    else if (ii > jj)
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

    ii = coarse_cell[f_face_cell[face_id][0]] - 1;
    jj = coarse_cell[f_face_cell[face_id][1]] - 1;

    if (ii != jj) {

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
}

/*----------------------------------------------------------------------------
 * Send prepared coarsening ids (based on halo send lists) and receive
 * similar values to a cell coarsening array halo.
 *
 * parameters:
 *   halo          <->  pointer to halo structure
 *   coarse_send   <->  values defined on halo send lists
 *   coarse_cell   <->  pointer to local cell coarsening array
 *                      whose halo values are to be received
 *----------------------------------------------------------------------------*/

static void
_exchange_halo_coarsening(const cs_halo_t  *halo,
                          cs_lnum_t         coarse_send[],
                          cs_lnum_t         coarse_cell[])
{
  cs_lnum_t  i, start, length;

  int local_rank_id = (cs_glob_n_ranks == 1) ? 0 : -1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int rank_id;
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

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      start = halo->index[2*rank_id];
      length = halo->index[2*rank_id + 2] - halo->index[2*rank_id];

      if (halo->c_domain_rank[rank_id] != local_rank) {

        MPI_Irecv(coarse_cell + halo->n_local_elts + start,
                  length,
                  CS_MPI_INT,
                  halo->c_domain_rank[rank_id],
                  halo->c_domain_rank[rank_id],
                  cs_glob_mpi_comm,
                  &(request[request_count++]));

      }
      else
        local_rank_id = rank_id;

    }

    /* We wait for posting all receives (often recommended) */

    MPI_Barrier(cs_glob_mpi_comm);

    /* Send data to distant ranks */

    for (rank_id = 0; rank_id < halo->n_c_domains; rank_id++) {

      /* If this is not the local rank */

      if (halo->c_domain_rank[rank_id] != local_rank) {

        start = halo->send_index[2*rank_id];
        length =  halo->send_index[2*rank_id + 2] - halo->send_index[2*rank_id];

        MPI_Isend(coarse_send + start,
                  length,
                  CS_MPI_INT,
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

      cs_lnum_t *_coarse_cell
        = coarse_cell + halo->n_local_elts + halo->index[2*local_rank_id];

      start = halo->send_index[2*local_rank_id];
      length =   halo->send_index[2*local_rank_id + 2]
               - halo->send_index[2*local_rank_id];

#     pragma omp parallel for if(length > CS_THR_MIN)
      for (i = 0; i < length; i++)
        _coarse_cell[i] = coarse_send[start + i];

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
  int domain_id, tr_id, section_id;
  cs_lnum_t ii, jj;
  cs_lnum_t start_id, end_id, sub_count;

  cs_lnum_t *start_end_id = NULL;
  cs_lnum_t *sub_num = NULL;
  cs_lnum_t  *coarse_send = NULL;

  cs_lnum_t *restrict coarse_cell = c->coarse_cell;

  cs_halo_t *c_halo = NULL;
  const cs_halo_t *f_halo = f->halo;

  const cs_lnum_t c_n_cells = c->n_cells;

  const int stride = f_halo->n_c_domains*4;
  const int n_sections = f_halo->n_transforms + 1;
  const int n_f_send = f_halo->n_send_elts[1]; /* Size of full list */

  c_halo = cs_halo_create_from_ref(f_halo);

  c->halo = c_halo;
  c->_halo = c_halo;

  /* Initialize coarse halo counters */

  c_halo->n_local_elts = c_n_cells;
  c_halo->n_send_elts[0] = 0;
  c_halo->n_send_elts[1] = 0;

# pragma omp parallel for if(f->n_cells_ext > CS_THR_MIN)
  for (ii = f->n_cells; ii < f->n_cells_ext; ii++)
    coarse_cell[ii] = 0;

  /* Allocate and initialize counters */

  BFT_MALLOC(start_end_id, n_sections*2, cs_lnum_t);
  BFT_MALLOC(sub_num, c_n_cells, cs_lnum_t);
  BFT_MALLOC(coarse_send, n_f_send, cs_lnum_t);

# pragma omp parallel for if(c_n_cells > CS_THR_MIN)
  for (ii = 0; ii < c_n_cells; ii++)
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
        jj = coarse_cell[f_halo->send_list[ii]] - 1;
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
        sub_num[coarse_cell[f_halo->send_list[ii]] - 1] = -1;

    }

    /* Update send index (initialized at halo creation) */

    c_halo->send_index[domain_id*2 + 1] = c_halo->n_send_elts[0];
    c_halo->send_index[domain_id*2 + 2] = c_halo->n_send_elts[0];

  }

  /* Exchange and update coarsening list halo */
  /*------------------------------------------*/

  _exchange_halo_coarsening(f_halo, coarse_send, coarse_cell);

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
        if (coarse_cell[ii] > jj)
          jj += 1;
        coarse_cell[ii] += c_n_cells + c_halo->n_elts[0] + 1;
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

  BFT_MALLOC(c_halo->send_list, c_halo->n_send_elts[0], cs_lnum_t);

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
        jj = coarse_cell[f_halo->send_list[ii]] - 1;
        if (sub_num[jj] == -1) {
          sub_num[jj] = sub_count;
          c_halo->send_list[c_halo->n_send_elts[0] + sub_count] = jj;
          sub_count += 1;
        }
      }

      c_halo->n_send_elts[0] += sub_count;

      /* Reset sub_num for next section or domain */

      for (ii = start_id; ii < end_id; ii++)
        sub_num[coarse_cell[f_halo->send_list[ii]] - 1] = -1;

    }

  }

  c_halo->n_send_elts[1] = c_halo->n_send_elts[0];

  /* Free memory */

  BFT_FREE(coarse_send);
  BFT_FREE(sub_num);
  BFT_FREE(start_end_id);
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
  cs_lnum_t  c_n_cells = 0;

  const cs_lnum_t f_n_faces = f->n_faces;
  const cs_lnum_2_t *restrict f_face_cell = f->face_cell;

  /* Sanity check */

# pragma omp parallel for if(f_n_faces > CS_THR_MIN)
  for (cs_lnum_t face_id = 0; face_id < f_n_faces; face_id++) {
    cs_lnum_t ii = f_face_cell[face_id][0];
    cs_lnum_t jj = f_face_cell[face_id][1];
    if (ii == jj)
      bft_error(__FILE__, __LINE__, 0,
                _("Connectivity error:\n"
                  "Face %d has same cell %d on both sides."),
                (int)(face_id+1), (int)(ii+1));
  }

  /* Compute number of coarse cells */

  for (cs_lnum_t ii = 0; ii < f->n_cells; ii++) {
    if (c->coarse_cell[ii] > c_n_cells)
      c_n_cells = c->coarse_cell[ii];
  }

  c->n_cells = c_n_cells;
  c->n_g_cells = c_n_cells;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t _c_n_cells = c_n_cells;
    MPI_Allreduce(&_c_n_cells, &(c->n_g_cells), 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  /* Prolong mesh coarsening indicator to halo cells and build
     coarse mesh halos if necessary */

  if (f->halo != NULL) {
    _coarsen_halo(f, c);
    c->n_cells_ext = c->n_cells + c->halo->n_elts[0];
  }
  else
    c->n_cells_ext = c_n_cells;

  c->n_cells_r[0] = c->n_cells;
  c->n_cells_r[1] = c->n_cells_ext;

  /* Build face coarsening and coarse grid face -> cells connectivity */

  _coarsen_faces(f,
                 c->coarse_cell,
                 c->n_cells,
                 &(c->n_faces),
                 &(c->coarse_face),
                 &(c->_face_cell));

  c->face_cell = (const cs_lnum_2_t  *)(c->_face_cell);
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of an array of real values.
 *
 * parameters:
 *   val       <-- pointer to entities that should be ordered.
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_real_descend_tree(const cs_real_t   val[],
                         size_t            level,
                         const size_t      nb_ent,
                         cs_lnum_t         order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (val[i1] > val[i2]) lv_cur++;
    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (val[i1] >= val[i2]) break;

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of doubles.
 *
 * parameters:
 *   val      <-- pointer to entities that should be ordered.
 *   order    <-- pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_real_local(const cs_real_t   val[],
                  cs_lnum_t         order[],
                  const size_t      nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0; i < nb_ent; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2);
  do {
    i--;
    _order_real_descend_tree(val, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_real_descend_tree(val, 0, i, order);
  }
}

/*----------------------------------------------------------------------------
 * Descend binary tree for the ordering of Hilbert code couples.
 *
 * parameters:
 *   val       <-- pointer to entities that should be ordered.
 *   level     <-- level of the binary tree to descend
 *   nb_ent    <-- number of entities in the binary tree to descend
 *   order     <-> ordering array
 *----------------------------------------------------------------------------*/

inline static void
_order_hilbert_descend_tree(const fvm_hilbert_code_t  val[],
                            size_t                    level,
                            const size_t              nb_ent,
                            cs_lnum_t                 order[])
{
  size_t i_save, i1, i2, lv_cur;

  i_save = (size_t)(order[level]);

  while (level <= (nb_ent/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < nb_ent - 1) {

      i1 = (size_t)(order[lv_cur+1]);
      i2 = (size_t)(order[lv_cur]);

      if (val[i1*2] >= val[i2*2]) {
        if (val[i1*2] > val[i2*2] || val[i1*2+1] > val[i2*2+1])
          lv_cur++;
      }
    }

    if (lv_cur >= nb_ent) break;

    i1 = i_save;
    i2 = (size_t)(order[lv_cur]);

    if (val[i1*2] >= val[i2*2]) {
      if (val[i1*2] > val[i2*2] || val[i1*2+1] >= val[i2*2+1])
        break;
    }

    order[level] = order[lv_cur];
    level = lv_cur;

  }

  order[level] = i_save;
}

/*----------------------------------------------------------------------------
 * Order an array of Hilbert codes couples.
 *
 * parameters:
 *   val      <-- pointer to entities that should be ordered.
 *   order    <-- pre-allocated ordering table
 *   nb_ent   <-- number of entities considered
 *----------------------------------------------------------------------------*/

static void
_order_hilbert2_local(const fvm_hilbert_code_t  val[],
                      cs_lnum_t                 order[],
                      const size_t              nb_ent)
{
  size_t i;
  cs_lnum_t o_save;

  /* Initialize ordering array */

  for (i = 0; i < nb_ent; i++)
    order[i] = i;

  if (nb_ent < 2)
    return;

  /* Create binary tree */

  i = (nb_ent / 2);
  do {
    i--;
    _order_hilbert_descend_tree(val, i, nb_ent, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = nb_ent - 1; i > 0; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_hilbert_descend_tree(val, 0, i, order);
  }
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Create a set of reduced communicators
 *
 * parameters:
 *   grid_merge_stride <-- size multiple between communicators
 *----------------------------------------------------------------------------*/

static void
_init_reduced_communicators(int  grid_merge_stride)
{
  int comm_id;
  int n_ranks;
  int ranges[1][3];
  MPI_Group old_group, new_group;

  _grid_merge_stride = grid_merge_stride;

  /* Determine number of communicators */

  _n_grid_comms = 1;
  n_ranks = cs_glob_n_ranks;
  while (n_ranks > 1) {
    _n_grid_comms += 1;
    if (n_ranks % _grid_merge_stride == 0)
      n_ranks = n_ranks / _grid_merge_stride;
    else
      n_ranks = (n_ranks / _grid_merge_stride) + 1;
  }

  BFT_MALLOC(_grid_comm, _n_grid_comms, MPI_Comm);
  BFT_MALLOC(_grid_ranks, _n_grid_comms, cs_lnum_t);

  n_ranks = cs_glob_n_ranks;
  _grid_ranks[0] = cs_glob_n_ranks;
  _grid_comm[0] = cs_glob_mpi_comm;

  MPI_Barrier(cs_glob_mpi_comm); /* For debugging */

  MPI_Comm_size(cs_glob_mpi_comm, &n_ranks);
  MPI_Comm_group(cs_glob_mpi_comm, &old_group);

  ranges[0][0] = 0;
  ranges[0][1] = n_ranks - 1;
  ranges[0][2] = 1;

  for (comm_id = 1; comm_id < _n_grid_comms; comm_id++) {

    if (n_ranks % _grid_merge_stride == 0)
      n_ranks = n_ranks / _grid_merge_stride;
    else
      n_ranks = (n_ranks / _grid_merge_stride) + 1;
    _grid_ranks[comm_id] = n_ranks;

    ranges[0][2] *= _grid_merge_stride;

    MPI_Group_range_incl(old_group, 1, ranges, &new_group);
    MPI_Comm_create(cs_glob_mpi_comm, new_group, &(_grid_comm[comm_id]));
    MPI_Group_free(&new_group);

  }

  MPI_Group_free(&old_group);

  MPI_Barrier(cs_glob_mpi_comm); /* For debugging */
}

/*----------------------------------------------------------------------------
 * Destroy a set of reduced communicators
 *----------------------------------------------------------------------------*/

static void
_finalize_reduced_communicators(void)
{
  int comm_id;

  for (comm_id = 1; comm_id < _n_grid_comms; comm_id++) {
    if (_grid_comm[comm_id] != MPI_COMM_NULL)
      MPI_Comm_free(&(_grid_comm[comm_id]));
  }

  BFT_FREE(_grid_comm);
  BFT_FREE(_grid_ranks);

  _n_grid_comms = 0;
}

/*----------------------------------------------------------------------------
 * Merge halo info after appending arrays.
 *
 * parameters:
 *   h                 <-> Pointer to halo structure
 *   new_src_cell_num  <-> new list of cells the senders should provide
 *----------------------------------------------------------------------------*/

static void
_rebuild_halo_send_lists(cs_halo_t  *h,
                         cs_lnum_t   new_src_cell_id[])
{
  /* As halos are based on interior faces, every domain is both
     sender and receiver */

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

  BFT_MALLOC(h->send_index, h->n_c_domains*2 + 1, cs_lnum_t);
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

  BFT_MALLOC(h->send_list, h->n_send_elts[0], cs_lnum_t);

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

  BFT_FREE(h->send_list);
  BFT_FREE(h->send_index);
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
  int  rank_id, prev_rank_id, tr_id, prev_section_id;
  cs_lnum_t  ii, ii_0, ii_1, cur_id, section_id, src_id, prev_src_id;

  int   stride = (h->n_transforms > 0) ? 3 : 2;
  int   n_c_domains_ini = h->n_c_domains;

  cs_lnum_t   n_elts_ini = h->n_elts[0];
  cs_lnum_t  *order = NULL, *section_idx = NULL;
  cs_gnum_t  *tmp_num = NULL;

  const int  n_sections = h->n_transforms + 1;

  if (h->n_elts[0] < 1) {
    _empty_halo(h);
    return;
  }

  /* Order list by rank, transform, and new element number */

  BFT_MALLOC(tmp_num, n_elts_ini*stride, cs_gnum_t);

  for (rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
    for (ii = h->index[rank_id*2];
         ii < h->index[(rank_id+1)*2];
         ii++)
      tmp_num[ii*stride] = h->c_domain_rank[rank_id];
  }

  if (stride == 2) {
    for (ii = 0; ii < n_elts_ini; ii++)
      tmp_num[ii*2 + 1] = new_src_cell_id[ii];
  }
  else { /* if (stride == 3) */

    for (ii = 0; ii < n_elts_ini; ii++) {
      tmp_num[ii*3 + 1] = 0;
      tmp_num[ii*3 + 2] = new_src_cell_id[ii];
    }

    for (rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
        ii_0 = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id];
        ii_1 = ii_0 + h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 1];
        for (ii = ii_0; ii < ii_1; ii++)
          tmp_num[ii*3 + 1] = tr_id + 1;
      }
    }
  }

  order = cs_order_gnum_s(NULL, tmp_num, stride, n_elts_ini);

  /* Rebuild lists and build renumbering */

  cur_id = order[0];
  h->n_c_domains = 0;
  h->index[0] = 0;
  h->n_elts[0] = 0;
  h->n_elts[1] = 0;

  prev_rank_id = -1;
  prev_src_id = 0;

  if (stride == 2) {

    for (ii = 0; ii < n_elts_ini; ii++) {

      bool is_same = true;
      cur_id = order[ii];

      rank_id = tmp_num[cur_id*2];
      src_id = tmp_num[cur_id*2 + 1];

      if (rank_id != loc_rank_id) {

        if (rank_id != prev_rank_id) {
          h->c_domain_rank[h->n_c_domains] = rank_id;
          h->index[h->n_c_domains*2] = h->n_elts[0];
          h->n_c_domains += 1;
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
      else { /* if (rank_id == loc_rank_id) */
        new_halo_cell_id[cur_id] = tmp_num[cur_id*2 + 1];
        assert(new_halo_cell_id[cur_id] < n_new_cells);
      }
    }

  }
  else { /* if (stride == 3) */

    int   rank_idx = -1;
    cs_lnum_t section_idx_size = 0;

    prev_section_id = -1;

    /* Index will be initialized as count */

    for (ii = 0; ii < n_elts_ini; ii++) {

      bool is_same = true;
      cur_id = order[ii];

      rank_id = tmp_num[cur_id*3];
      section_id = tmp_num[cur_id*3 + 1];
      src_id = tmp_num[cur_id*3 + 2];

      if (rank_id != loc_rank_id || tmp_num[cur_id*3 + 1] != 0) {

        if (rank_id != prev_rank_id) {
          h->c_domain_rank[h->n_c_domains] = rank_id;
          h->index[h->n_c_domains*2] = h->n_elts[0];
          h->n_c_domains += 1;
          is_same = false;
          rank_idx += 1;
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
      else { /* if (rank_id == loc_rank_id && tmp_num[cur_id*3 + 1] == 0) */
        new_halo_cell_id[cur_id] = tmp_num[cur_id*3 + 2];
        assert(new_halo_cell_id[cur_id] < n_new_cells);
      }

    }

    /* Transform count to index */

    section_idx_size = n_sections * h->n_c_domains + 1;
    for (cs_lnum_t sid = 1; sid < section_idx_size; sid++)
      section_idx[sid] += section_idx[sid - 1];
  }

  if (h->n_transforms > 0) {

    for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
      for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
        section_id = rank_id*n_sections + tr_id + 1;
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id]
          = section_idx[section_id];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 1]
          = section_idx[section_id + 1] - section_idx[section_id];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 2]
          = section_idx[section_id + 1];
        h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 3]
          = 0;
      }
    }

    BFT_FREE(section_idx);

  }

  /* Free and resize memory */

  BFT_FREE(order);
  BFT_FREE(tmp_num);
  BFT_REALLOC(h->c_domain_rank, h->n_c_domains, int);
  BFT_REALLOC(h->index, h->n_c_domains*2+1, cs_lnum_t);
  if (h->n_transforms > 0)
    BFT_REALLOC(h->perio_lst,
                h->n_c_domains * h->n_transforms * 4,
                cs_lnum_t);

  h->n_elts[1] = h->n_elts[0];
  h->index[h->n_c_domains*2] = h->n_elts[0];
  for (ii = 0; ii < h->n_c_domains; ii++)
    h->index[ii*2+1] = h->index[ii*2+2];
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
  BFT_FREE(h->send_list);
  BFT_FREE(h->send_index);
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
    for (ii = g->n_cells, jj = 0; ii < g->n_cells_ext; ii++, jj++)
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
    for (ii = g->n_cells, jj = 0; ii < g->n_cells_ext; ii++, jj++)
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

  /* Cleanup halo and set pointer for coarsening (sub_rooot) ranks*/

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

    _rebuild_halo_send_lists(h, new_src_cell_id);

  }

  if (new_src_cell_id != NULL)
    BFT_FREE(new_src_cell_id);

  cs_halo_update_buffers(h);

  g->halo = h;
  g->_halo = h;

  /* Finally, update halo section of cell renumbering array */

  if (g->merge_sub_rank == 0) {

    cs_lnum_t n_send = recv_count[2];
    cs_lnum_t send_shift = n_send;

    for (ii = 0; ii < n_send; ii++)
      new_cell_id[g->n_cells + ii] = new_halo_cell_id[ii];

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
    MPI_Recv(new_cell_id + g->n_cells, g->n_cells_ext - g->n_cells,
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

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'c';

  /* Reallocate arrays for receiving rank and append data */

  if (g->merge_sub_rank == 0) {

    g->n_cells = g->merge_cell_idx[g->merge_sub_size];
    g->n_cells_ext = g->n_cells + g->halo->n_elts[0];

    BFT_REALLOC(g->_cell_cen, g->n_cells_ext * 3, cs_real_t);
    BFT_REALLOC(g->_cell_vol, g->n_cells_ext, cs_real_t);
    BFT_REALLOC(g->_da, g->n_cells_ext*g->diag_block_size[3], cs_real_t);

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = (  g->merge_cell_idx[rank_id+1]
                           - g->merge_cell_idx[rank_id]);
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_cell_cen + g->merge_cell_idx[rank_id]*3, n_recv*3,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->_cell_vol + g->merge_cell_idx[rank_id], n_recv,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->_da + g->merge_cell_idx[rank_id]*g->diag_block_size[3],
               n_recv*g->diag_block_size[3],
               CS_MPI_REAL, dist_rank, tag, comm, &status);

    }

  }
  else {
    MPI_Send(g->_cell_cen, g->n_cells*3, CS_MPI_REAL, g->merge_sub_root,
             tag, comm);
    BFT_FREE(g->_cell_cen);

    MPI_Send(g->_cell_vol, g->n_cells, CS_MPI_REAL, g->merge_sub_root,
             tag, comm);
    BFT_FREE(g->_cell_vol);

    MPI_Send(g->_da, g->n_cells*g->diag_block_size[3], CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_da);

    g->n_cells = 0;
    g->n_cells_ext = 0;
  }

  g->cell_cen = g->_cell_cen;
  g->cell_vol = g->_cell_vol;
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

    cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD, g->_cell_cen, 3);
    if (g->halo->n_transforms > 0)
      cs_halo_perio_sync_coords(g->halo, CS_HALO_STANDARD, g->_cell_cen);

    cs_halo_sync_var(g->halo, CS_HALO_STANDARD, g->_cell_vol);

    cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD,
                             g->_da, g->diag_block_size[3]);

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

    BFT_REALLOC(g->_face_normal, n_faces_tot*3, cs_real_t);

    if (g->symmetric == true)
      BFT_REALLOC(g->_xa, n_faces_tot, cs_real_t);
    else
      BFT_REALLOC(g->_xa, n_faces_tot*2, cs_real_t);

    BFT_REALLOC(g->_xa0, n_faces_tot, cs_real_t);
    BFT_REALLOC(g->xa0ij, n_faces_tot*3, cs_real_t);

    g->n_faces = recv_count[0];

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      cs_lnum_t n_recv = recv_count[rank_id];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_face_cell + g->n_faces, n_recv*2,
               CS_MPI_LNUM, dist_rank, tag, comm, &status);

      MPI_Recv(g->_face_normal + g->n_faces*3, n_recv*3,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      if (g->symmetric == true)
        MPI_Recv(g->_xa + g->n_faces, n_recv,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);
      else
        MPI_Recv(g->_xa + g->n_faces*2, n_recv*2,
                 CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->_xa0 + g->n_faces, n_recv,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->xa0ij + g->n_faces*3, n_recv*3,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

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
      g->_face_normal[face_id*3] = g->_face_normal[p_face_id*3];
      g->_face_normal[face_id*3 + 1] = g->_face_normal[p_face_id*3 + 1];
      g->_face_normal[face_id*3 + 2] = g->_face_normal[p_face_id*3 + 2];
      if (g->symmetric == true)
        g->_xa[face_id] = g->_xa[p_face_id];
      else {
        g->_xa[face_id*2] = g->_xa[p_face_id*2];
        g->_xa[face_id*2+1] = g->_xa[p_face_id*2+1];
      }
      g->_xa0[face_id] = g->_xa0[p_face_id];
      g->xa0ij[face_id*3] = g->xa0ij[p_face_id*3];
      g->xa0ij[face_id*3 + 1] = g->xa0ij[p_face_id*3 + 1];
      g->xa0ij[face_id*3 + 2] = g->xa0ij[p_face_id*3 + 2];
    }

    MPI_Send(g->_face_cell, n_faces*2, CS_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_face_cell);

    MPI_Send(g->_face_normal, n_faces*3, CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_face_normal);

    if (g->symmetric == true)
      MPI_Send(g->_xa, n_faces, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    else
      MPI_Send(g->_xa, n_faces*2, CS_MPI_REAL,
               g->merge_sub_root, tag, comm);
    BFT_FREE(g->_xa);

    MPI_Send(g->_xa0, n_faces, CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->_xa0);

    MPI_Send(g->xa0ij, n_faces*3, CS_MPI_REAL,
             g->merge_sub_root, tag, comm);
    BFT_FREE(g->xa0ij);

    g->n_faces = 0;
  }

  g->face_cell = (const cs_lnum_2_t  *)(g->_face_cell);
  g->face_normal = g->_face_normal;
  g->xa = g->_xa;
  g->xa0 = g->_xa0;
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
 *   g         <-- Pointer to grid structure
 *   verbosity <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_merge_grids(cs_grid_t  *g,
             int         verbosity)
{
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

  if (_grid_merge_stride > 1) {
     if (_n_grid_comms == 0)
       _init_reduced_communicators(_grid_merge_stride);
  }

  /* Determine rank in merged group */

  g->merge_sub_size = _grid_merge_stride;
  g->merge_stride = g->next_merge_stride;
  g->next_merge_stride *= _grid_merge_stride;

  if (base_rank % g->merge_stride != 0) {
    g->merge_sub_size = 0;
    g->merge_sub_root = -1;
    g->merge_sub_rank = -1;
  }
  else {
    int merge_rank = base_rank / g->merge_stride;
    int merge_size = (cs_glob_n_ranks/g->merge_stride);
    int merge_sub_root = (merge_rank / _grid_merge_stride) * _grid_merge_stride;
    if (cs_glob_n_ranks % g->merge_stride > 0)
      merge_size += 1;
    g->merge_sub_root = merge_sub_root * g->merge_stride;
    g->merge_sub_rank = merge_rank % _grid_merge_stride;
    if (merge_sub_root + g->merge_sub_size > merge_size)
      g->merge_sub_size = merge_size - merge_sub_root;
  }

  if (g->next_merge_stride > 1) {
    if (g->n_ranks % _grid_merge_stride)
      g->n_ranks = (g->n_ranks/_grid_merge_stride) + 1;
    else
      g->n_ranks = (g->n_ranks/_grid_merge_stride);
    g->comm_id = 0;
    while (   _grid_ranks[g->comm_id] != g->n_ranks
           && g->comm_id < _n_grid_comms)
      g->comm_id++;
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
      g->merge_cell_idx[0] = 0; g->merge_cell_idx[1] = g->n_cells;
      for (i = 1; i < g->merge_sub_size; i++) {
        cs_lnum_t recv_val;
        int dist_rank = g->merge_sub_root + g->merge_stride*i;
        MPI_Recv(&recv_val, 1, CS_MPI_LNUM, dist_rank, tag, comm, &status);
        g->merge_cell_idx[i + 1] = recv_val + g->merge_cell_idx[i];
      }
    }
    else {
      cs_lnum_t send_val = g->n_cells;
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

  BFT_MALLOC(new_cell_id, g->n_cells_ext, cs_lnum_t);
  for (j = 0; j < g->n_cells; j++)
    new_cell_id[j] = cell_shift + j;
  for (j = g->n_cells; j < g->n_cells_ext; j++)
    new_cell_id[j] = -1;

  cs_halo_sync_untyped(g->halo,
                       CS_HALO_STANDARD,
                       sizeof(cs_lnum_t),
                       new_cell_id);

  /* Now build face filter list (before halo is modified) */

  if (g->merge_sub_size > 1 && g->merge_sub_rank > 0 && g->n_faces > 0) {

    cs_lnum_t n_ghost_cells = g->n_cells_ext - g->n_cells;

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
      cs_lnum_t ii = g->face_cell[face_id][0] - g->n_cells + 1;
      cs_lnum_t jj = g->face_cell[face_id][1] - g->n_cells + 1;
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
    bft_printf("      merged to %d (from %d) cells\n\n",
               g->n_cells, g->n_cells_r[0]);
}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Determine face traversal order.
 *
 * parameters:
 *   g               <-- Grid structure
 *   coarsening_type <-- Coarsening type
 *   order           --> Face ordering
 *----------------------------------------------------------------------------*/

static void
_order_face_traversal(const cs_grid_t   *g,
                      int                coarsening_type,
                      cs_lnum_t          order[])
{
  cs_lnum_t ii;
  cs_lnum_t n_faces = g->n_faces;

  if (coarsening_type == 2) {

    cs_coord_t extents[6];
    fvm_hilbert_code_t  *c_h_code = NULL;
    fvm_hilbert_code_t  *f_h_code = NULL;

    BFT_MALLOC(f_h_code, g->n_faces*2, fvm_hilbert_code_t);
    BFT_MALLOC(c_h_code, g->n_cells_ext, fvm_hilbert_code_t);

#if defined(HAVE_MPI)
    fvm_hilbert_get_coord_extents(3, g->n_cells_ext, g->cell_cen, extents,
                                  MPI_COMM_NULL);
#else
    fvm_hilbert_get_coord_extents(3, g->n_cells_ext, g->cell_cen, extents);
#endif

    fvm_hilbert_encode_coords(3, extents, g->n_cells_ext, g->cell_cen,
                              c_h_code);

    for (ii = 0; ii < n_faces; ii++) {
      double c0 = c_h_code[g->face_cell[ii][0]];
      double c1 = c_h_code[g->face_cell[ii][1]];
      if (c0 < c1) {
        f_h_code[ii*2] = c0;
        f_h_code[ii*2+1] = c1;
      }
      else {
        f_h_code[ii*2] = c1;
        f_h_code[ii*2+1] = c0;
      }
    }

    BFT_FREE(c_h_code);

    _order_hilbert2_local(f_h_code, order, g->n_faces);

    BFT_FREE(f_h_code);
  }

  else {
    for (ii = 0; ii < n_faces; ii++)
      order[ii] = ii;
  }

}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Scatter coarse cell integer data in case of merged grids
 *
 * parameters:
 *   g    <-- Grid structure
 *   num  <-> Variable defined on merged grid cells in,
 *            defined on scattered to grid cells out
 *----------------------------------------------------------------------------*/

static void
_scatter_cell_num(const cs_grid_t  *g,
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
        MPI_Send(num + g->merge_cell_idx[rank_id], n_send, MPI_INT,
                 dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(num, g->n_cells_r[0], MPI_INT,
               g->merge_sub_root, tag, comm, &status);
    }
  }

}

#endif /* defined(HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Build a coarse grid level from the previous level using
 * an automatic criterion.
 *
 * parameters:
 *   f                    <-- Fine grid structure
 *   coarsening_type      <-- Coarsening type
 *   max_aggregation      <-- Max fine cells per coarse cell
 *   relaxation_parameter <-- P0/P1 relaxation factor
 *   verbosity            <-- Verbosity level
 *   f_c_cell             --> Fine cell -> coarse cell connectivity
 *----------------------------------------------------------------------------*/

static void
_automatic_aggregation(const cs_grid_t  *fine_grid,
                       int               coarsening_type,
                       cs_lnum_t         max_aggregation,
                       double            relaxation_parameter,
                       int               verbosity,
                       cs_lnum_t        *f_c_cell)
{
  cs_lnum_t ii, jj, kk, face_id, c_face, face_f;
  cs_lnum_t n_faces;
  cs_real_t f_xa1, f_xa2;
  cs_real_t trace1, trace2;

  int isym = 2;
  int ncoarse = 8, npass_max = 10, inc_nei = 1;
  int _max_aggregation = 1, npass = 0;

  cs_lnum_t f_n_cells = fine_grid->n_cells;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cells_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t aggr_count = f_n_cells, count = 0;

  cs_lnum_t r_n_faces = f_n_faces;
  cs_lnum_t c_n_cells = 0;

  cs_real_t epsilon = 1.e-6;

  cs_lnum_t *face_order = NULL;
  cs_lnum_t *c_cardinality = NULL, *f_c_face = NULL;
  cs_lnum_t *merge_flag = NULL, *c_aggr_count = NULL;
  cs_lnum_t *i_work_array = NULL;
  cs_real_t *r_work_array = NULL;

  cs_real_t *aggr_crit = NULL;

  const int *db_size = fine_grid->diag_block_size;
  const cs_lnum_2_t *f_face_cells = fine_grid->face_cell;
  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_xa = fine_grid->xa;

  /* Allocate working arrays */

  BFT_MALLOC(i_work_array, f_n_cells_ext*2 + f_n_faces*3, cs_lnum_t);
  BFT_MALLOC(r_work_array, f_n_faces, cs_real_t);

  c_cardinality = i_work_array;
  c_aggr_count = i_work_array + f_n_cells_ext;
  f_c_face = i_work_array + 2*f_n_cells_ext; /* fine face -> coarse face */
  merge_flag = i_work_array + 2*f_n_cells_ext + f_n_faces;
  face_order = merge_flag + f_n_faces;

  aggr_crit = r_work_array;

  /* Initialization */

  _order_face_traversal(fine_grid, coarsening_type, face_order);

  if (fine_grid->symmetric == true)
    isym = 1;

# pragma omp parallel for if(f_n_cells_ext > CS_THR_MIN)
  for (ii = 0; ii < f_n_cells_ext; ii++) {
    c_cardinality[ii] = -1;
    f_c_cell[ii] = 0;
    c_aggr_count[ii] = 1;
  }

# pragma omp parallel for if(f_n_faces > CS_THR_MIN)
  for (face_id = 0; face_id < f_n_faces; face_id++) {
    merge_flag[face_id] = face_id +1;
    f_c_face[face_id] = 0;
    aggr_crit[face_id] = 2.;
  }

  /* Compute cardinality (number of neighbors for each cell -1) */

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cells[face_id][0];
    jj = f_face_cells[face_id][1];

    c_cardinality[ii] += 1;
    c_cardinality[jj] += 1;
  }

  /* Passes */

  if (verbosity > 3)
    bft_printf("\n     %s:\n", __func__);

  do {

    npass++;
    n_faces = r_n_faces;
    _max_aggregation++;
    _max_aggregation = CS_MIN(_max_aggregation, max_aggregation);

#   pragma omp parallel for if(n_faces > CS_THR_MIN)
    for (face_id = 0; face_id < n_faces; face_id++) {
      f_c_face[face_id] = merge_flag[face_id];
      merge_flag[face_id] = 0;
    }

    if (n_faces < f_n_faces) {
#     pragma omp parallel for if(f_n_faces > CS_THR_MIN)
      for (face_id = n_faces; face_id < f_n_faces; face_id++) {
        merge_flag[face_id] = 0;
        f_c_face[face_id] = 0;
      }
    }

    if (verbosity > 3)
      bft_printf("       pass %3d; r_n_faces = %10ld; aggr_count = %10ld\n",
                 npass, (long)r_n_faces, (long)aggr_count);

    /* Increment number of neighbors */

#   pragma omp parallel for if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++)
      c_cardinality[ii] += inc_nei;

    /* Initialize non-eliminated faces */
    r_n_faces = 0;

    if (fine_grid->conv_diff) {

      for (face_id = 0; face_id < n_faces; face_id++) {

        c_face = f_c_face[face_id] -1;

        ii = f_face_cells[c_face][0];
        jj = f_face_cells[c_face][1];

        /* Exclude faces on parallel or periodic boundary, so as not to */
        /* coarsen the grid across those boundaries (which would change */
        /* the communication pattern and require a more complex algorithm). */

        if (ii < f_n_cells && jj < f_n_cells) {
          f_xa1 = f_xa[2*c_face];
          f_xa2 = f_xa[2*c_face +1];
          /* TODO: remove these tests, or adimensionalize them */
          f_xa1 = CS_MAX(-f_xa1, 1.e-15);
          f_xa2 = CS_MAX(-f_xa2, 1.e-15);

          aggr_crit[face_id] =   CS_MAX(f_xa1, f_xa2)
                               / CS_MAX(sqrt((f_da[ii] / c_cardinality[ii])
                                             * (f_da[jj] / c_cardinality[jj])),
                                        1.e-15);
        }
      }

    }
    else {
      for (face_id = 0; face_id < n_faces; face_id++) {

        c_face = f_c_face[face_id] -1;

        ii = f_face_cells[c_face][0];
        jj = f_face_cells[c_face][1];

        /* Exclude faces on parallel or periodic boundary, so as not to */
        /* coarsen the grid across those boundaries (which would change */
        /* the communication pattern and require a more complex algorithm). */

        if (ii < f_n_cells && jj < f_n_cells) {
          f_xa1 = f_xa[c_face*isym];//TODO if f_xa1 is a block
          f_xa2 = f_xa[(c_face +1)*isym -1];
          /* TODO: remove these tests, or adimensionalize them */
          f_xa1 = CS_MAX(-f_xa1, 1.e-15);
          f_xa2 = CS_MAX(-f_xa2, 1.e-15);

          trace1 = 0.;
          trace2 = 0.;
          for (kk = 0; kk < db_size[0]; kk++) {
            trace1 += f_da[ii*db_size[3] + db_size[2]*kk + kk];
            trace2 += f_da[jj*db_size[3] + db_size[2]*kk + kk];
          }
          aggr_crit[face_id] =   (trace1/c_cardinality[ii])
                               * (trace2/c_cardinality[jj])
                               / (f_xa1*db_size[0]*f_xa2*db_size[0]);
        }

      }
    }

    /* Order faces by criteria (0 to n-1 numbering) */

    if (coarsening_type == 1)
      _order_real_local(aggr_crit, face_order, n_faces);

    /* Loop on non-eliminated faces */

    cs_real_t ag_mult = 1.;
    cs_real_t ag_threshold = 1. - epsilon;
    if (fine_grid->conv_diff) {
      // ag_mult = -1.;
      // ag_threshold = - (1. - epsilon) * pow(relaxation_parameter, npass);
      ag_threshold = (1. - epsilon) * pow(relaxation_parameter, npass);
    }

    for (face_id = 0; face_id < n_faces; face_id++) {

      face_f = face_order[face_id];
      c_face = f_c_face[face_f] -1;
      ii = f_face_cells[c_face][0];
      jj = f_face_cells[c_face][1];

      /* Exclude faces on parallel or periodic boundary, so as not to */
      /* coarsen the grid across those boundaries (which would change */
      /* the communication pattern and require a more complex algorithm). */

      if (ii < f_n_cells && jj < f_n_cells) {
        count = 0;

        if (   ag_mult*aggr_crit[face_f] < ag_threshold
            && (f_c_cell[ii] < 1 || f_c_cell[jj] < 1)) {

          if (f_c_cell[ii] > 0 && f_c_cell[jj] < 1 ) {
            if (c_aggr_count[f_c_cell[ii] -1] < _max_aggregation +1) {
              f_c_cell[jj] = f_c_cell[ii];
              c_aggr_count[f_c_cell[ii] -1] += 1;
              count++;
            }
          }
          else if (f_c_cell[ii] < 1 && f_c_cell[jj] > 0) {
            if (c_aggr_count[f_c_cell[jj] -1] < _max_aggregation +1) {
              f_c_cell[ii] = f_c_cell[jj];
              c_aggr_count[f_c_cell[jj] -1] += 1;
              count++;
            }
          }
          else if (f_c_cell[ii] < 1 && f_c_cell[jj] < 1) {
            f_c_cell[ii] = c_n_cells +1;
            f_c_cell[jj] = c_n_cells +1;
            c_aggr_count[c_n_cells] += 1;
            c_n_cells++;
            count++;
          }
        }
        if (count == 0 && (f_c_cell[ii] < 1 || f_c_cell[jj] < 1)) {
          merge_flag[r_n_faces] = c_face +1;
          face_order[r_n_faces] = r_n_faces;
          r_n_faces++;
        }

      }

    }

    /* Check the number of coarse cells created */
    aggr_count = 0;
#   pragma omp parallel for reduction(+:aggr_count) if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      if (f_c_cell[ii] < 1)
        aggr_count++;
    }

    /* Additional passes if aggregation is insufficient */

    if (   aggr_count == 0 || (c_n_cells + aggr_count)*ncoarse < f_n_cells
        || r_n_faces == 0)
      npass_max = npass;

  } while (npass < npass_max); /* Loop on passes */

  /* Finish assembly */
  for (ii = 0; ii < f_n_cells; ii++) {
    if (f_c_cell[ii] < 1) {
      f_c_cell[ii] = c_n_cells + 1;
      c_n_cells++;
    }
  }

  /* Various checks */

  if (verbosity > 3) {

    cs_lnum_t ic;
    cs_lnum_t aggr_min = 2*f_n_cells, aggr_max = 0;
    cs_gnum_t c_n_cells_g = c_n_cells;

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
    MPI_Comm comm = cs_glob_mpi_comm;
    if (fine_grid->n_ranks > 1) {
      if (_grid_ranks != NULL)
        comm = _grid_comm[fine_grid->comm_id];
    }
#endif

    for (ii = 0; ii < c_n_cells; ii++) {
      aggr_min = CS_MIN(aggr_min, c_aggr_count[ii]);
      aggr_max = CS_MAX(aggr_max, c_aggr_count[ii]);
    }

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
    if (comm != MPI_COMM_NULL) {
      MPI_Allreduce(MPI_IN_PLACE, &aggr_min, 1, CS_MPI_LNUM, MPI_MIN, comm);
      MPI_Allreduce(MPI_IN_PLACE, &aggr_max, 1, CS_MPI_LNUM, MPI_MAX, comm);
      MPI_Allreduce(MPI_IN_PLACE, &c_n_cells_g, 1, CS_MPI_GNUM, MPI_SUM, comm);
    }
#endif

    bft_printf("       aggregation min = %10ld; max = %10ld; target = %10ld\n",
               (long)aggr_min, (long)aggr_max, (long)max_aggregation);

    bft_printf("       histogram\n");

    aggr_count = aggr_max - aggr_min + 1;
    if (aggr_count > 0) {
      cs_lnum_t i;
      cs_lnum_t *histogram;
      BFT_MALLOC(histogram, aggr_count, cs_lnum_t);
      for (i = 0; i < aggr_count; i++)
        histogram[i] = 0;
      for (ic = 0; ic < c_n_cells; ic++) {
        for (i = 0; i < aggr_count; i++) {
          if (c_aggr_count[ic] == aggr_min + i)
            histogram[i] += 1;
        }
      }
#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
      if (comm != MPI_COMM_NULL)
        MPI_Allreduce(MPI_IN_PLACE, histogram, aggr_count, CS_MPI_LNUM, MPI_SUM,
                      comm);
#endif
      if (c_n_cells_g > 0) {
        for (i = 0; i < aggr_count; i++) {
          double epsp = 100. * histogram[i] / c_n_cells_g;
          bft_printf("         aggregation %10ld = %12.5f %%\n",
                     (long)(aggr_min + i), epsp);
        }
      }
      BFT_FREE(histogram);
    }

    /* Additional checks */

    for (ii = 0; ii < f_n_cells; ii++)
      c_cardinality[ii] = 0;
    for (ii = 0; ii < f_n_cells; ii++)
      c_cardinality[f_c_cell[ii] - 1] += 1;

    ii = 0;
    jj = 2*f_n_faces;
    aggr_count = 0;
    for (ic = 0; ic < c_n_cells; ic++) {
      ii = CS_MAX(ii, c_cardinality[ic]);
      jj = CS_MIN(jj, c_cardinality[ic]);
      aggr_count += c_cardinality[ic];
    }

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
    if (comm != MPI_COMM_NULL) {
      MPI_Allreduce(MPI_IN_PLACE, &aggr_min, 1, CS_MPI_LNUM, MPI_MIN, comm);
      MPI_Allreduce(MPI_IN_PLACE, &aggr_max, 1, CS_MPI_LNUM, MPI_MAX, comm);
    }
#endif

    bft_printf("       cell aggregation min = %10ld; max = %10ld\n",
               (long)aggr_min, (long)aggr_max);

  /* Consistency check */

    if (aggr_count != f_n_cells)
      bft_error(__FILE__, __LINE__, 0,
                "Inconsistency in %s aggregation.", __func__);

  }

  /* Free working arrays */

  BFT_FREE(i_work_array);
  BFT_FREE(r_work_array);
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
  cs_lnum_t ic, ii;

  cs_lnum_t f_n_cells = fine_grid->n_cells;

  cs_lnum_t c_n_cells = coarse_grid->n_cells;
  cs_lnum_t c_n_cells_ext = coarse_grid->n_cells_ext;

  cs_lnum_t *c_coarse_cell = coarse_grid->coarse_cell;

  cs_real_t *c_cell_vol = coarse_grid->_cell_vol;
  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;

  const cs_real_t *f_cell_vol = fine_grid->cell_vol;
  const cs_real_t *f_cell_cen = fine_grid->cell_cen;

  /* Compute volume and center of coarse cells */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (ic = 0; ic < c_n_cells_ext; ic++) {
    c_cell_vol[ic] = 0.;
    c_cell_cen[3*ic]    = 0.;
    c_cell_cen[3*ic +1] = 0.;
    c_cell_cen[3*ic +2] = 0.;
  }

  for (ii = 0; ii < f_n_cells; ii++) {
    ic = c_coarse_cell[ii] -1;
    c_cell_vol[ic] += f_cell_vol[ii];
    c_cell_cen[3*ic]    += f_cell_vol[ii]*f_cell_cen[3*ii];
    c_cell_cen[3*ic +1] += f_cell_vol[ii]*f_cell_cen[3*ii +1];
    c_cell_cen[3*ic +2] += f_cell_vol[ii]*f_cell_cen[3*ii +2];
  }

# pragma omp parallel for if(c_n_cells > CS_THR_MIN)
  for (ic = 0; ic < c_n_cells; ic++) {
    c_cell_cen[3*ic]    /= c_cell_vol[ic];
    c_cell_cen[3*ic +1] /= c_cell_vol[ic];
    c_cell_cen[3*ic +2] /= c_cell_vol[ic];
  }
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
  cs_lnum_t ic, jc, ii, jj, kk, c_face, face_id;

  int isym = 2;

  cs_lnum_t f_n_cells = fine_grid->n_cells;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cells_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells = coarse_grid->n_cells;
  cs_lnum_t c_n_cells_ext = coarse_grid->n_cells_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  const cs_real_t *c_xa0 = coarse_grid->_xa0;
  const cs_real_t *c_da = coarse_grid->_da;
  const cs_real_t *c_xa = coarse_grid->_xa;

  cs_real_t *w1 = NULL;

  const int *db_size = fine_grid->diag_block_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_xa = fine_grid->xa;

  BFT_MALLOC(w1, f_n_cells_ext*db_size[3], cs_real_t);

  if (fine_grid->symmetric == true)
    isym = 1;

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
  MPI_Comm comm = cs_glob_mpi_comm;
  if (fine_grid->n_ranks > 1) {
    if (_grid_ranks != NULL)
      comm = _grid_comm[fine_grid->comm_id];
    if (comm != MPI_COMM_NULL) {
      cs_gnum_t n_clips[2] = {n_clips_min, n_clips_max};
      MPI_Allreduce(MPI_IN_PLACE, n_clips, 2, CS_MPI_GNUM, MPI_SUM, comm);
      n_clips_min = n_clips[0];
      n_clips_max = n_clips[0];
    }
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

  BFT_MALLOC(w2, f_n_cells_ext*db_size[3], cs_real_t);
  BFT_MALLOC(w3, c_n_cells_ext*db_size[3], cs_real_t);
  BFT_MALLOC(w4, c_n_cells_ext*db_size[3], cs_real_t);

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
    w1[ii] = CS_MAX(fabs(f_xa[face_id*isym]), w1[ii]);
    w2[ii] = CS_MIN(fabs(f_xa[face_id*isym]), w2[ii]);
    w1[jj] = CS_MAX(fabs(f_xa[(face_id +1)*isym -1]), w1[jj]);
    w2[jj] = CS_MIN(fabs(f_xa[(face_id +1)*isym -1]), w2[jj]);
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    w3[ic] = CS_MAX(fabs(c_xa[c_face*isym]), w3[ic]);
    w4[ic] = CS_MIN(fabs(c_xa[c_face*isym]), w4[ic]);
    w3[jc] = CS_MAX(fabs(c_xa[(c_face +1)*isym -1]), w3[jc]);
    w4[jc] = CS_MIN(fabs(c_xa[(c_face +1)*isym -1]), w4[jc]);
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
      rmin = CS_MIN(rmin, c_xa[c_face*isym] / c_xa0[c_face]);
      rmax = CS_MAX(rmax, c_xa[c_face*isym] / c_xa0[c_face]);
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

  /* Evaluate fine and coarse grids diagonal dominance */

  for (ii = 0; ii < f_n_cells; ii++) {
    for (jj = 0; jj < db_size[0]; jj++) {
      for (kk = 0; kk < db_size[0]; kk++)
        w1[ii*db_size[3] + db_size[2]*jj + kk]
          = fabs(f_da[ii*db_size[3] + db_size[2]*jj + kk]);
    }
  }
  for (ii = f_n_cells; ii < f_n_cells_ext; ii++) {
    for (jj = 0; jj < db_size[0]; jj++) {
      for (kk = 0; kk < db_size[0]; kk++)
        w1[ii*db_size[3] + db_size[2]*jj + kk] = 0.;
    }
  }

  for (ic = 0; ic < c_n_cells; ic++) {
    for (jj = 0; jj < db_size[0]; jj++) {
      for (kk = 0; kk < db_size[0]; kk++)
        w3[ic*db_size[3] + db_size[2]*jj + kk]
          = fabs(c_da[ic*db_size[3] + db_size[2]*jj + kk]);
    }
  }
  for (ic = c_n_cells; ic < c_n_cells_ext; ic++) {
    for (jj = 0; jj < db_size[0]; jj++) {
      for (kk = 0; kk < db_size[0]; kk++)
        w3[ic*db_size[3] + db_size[2]*jj + kk] = 0.;
    }
  }

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cell[face_id][0];
    jj = f_face_cell[face_id][1];
    for (kk = 0; kk < db_size[0]; kk++) {
      w1[ii*db_size[3] + db_size[2]*kk + kk]
        -= fabs(f_xa[face_id*isym]);
      w1[jj*db_size[3] + db_size[2]*kk + kk]
        -= fabs(f_xa[(face_id +1)*isym -1]);
    }
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    for (kk = 0; kk < db_size[0]; kk++) {
      w3[ic*db_size[3] + db_size[2]*kk + kk]
        -= fabs(c_xa[c_face*isym]);
      w3[jc*db_size[3] + db_size[2]*kk + kk]
        -= fabs(c_xa[(c_face +1)*isym -1]);
    }
  }

  for (ii = 0; ii < f_n_cells; ii++)
    w1[ii] /= fabs(f_da[ii]);

  for (ic = 0; ic < c_n_cells; ic++)
    w3[ic] /= fabs(c_da[ic]);

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

  BFT_FREE(w3);
  BFT_FREE(w1);

#if defined(HAVE_MPI) && defined(HAVE_MPI_IN_PLACE)
  if (comm != MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE, anmin, 2, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, anmax, 2, MPI_DOUBLE, MPI_MAX, comm);
  }
#endif

  bft_printf(_("       fine   grid diag. dominance: min = %12.5e\n"
               "                                    max = %12.5e\n"
               "       coarse grid diag. dominance: min = %12.5e\n"
               "                                    max = %12.5e\n\n"),
             anmin[0], anmax[0], anmin[1], anmax[1]);
}

/*----------------------------------------------------------------------------
 * Build a coarse level from a finer level.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-> Coarse grid structure
 *   relax_param <-- P0/P1 relaxation parameter
 *   verbosity   <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities(const cs_grid_t  *fine_grid,
                           cs_grid_t        *coarse_grid,
                           cs_real_t         relax_param,
                           int               verbosity)
{
  int interp;

  cs_lnum_t ic, jc, ii, jj, kk, c_face, face_id;

  cs_real_t dsigjg, dsxaij, agij;

  int isym = 2;

  cs_lnum_t f_n_cells = fine_grid->n_cells;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cells_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells_ext = coarse_grid->n_cells_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  cs_lnum_t *c_coarse_cell = coarse_grid->coarse_cell;
  cs_lnum_t *c_coarse_face = coarse_grid->coarse_face;

  cs_gnum_t n_clips_min = 0, n_clips_max = 0;

  cs_real_t *f_xa0ij = fine_grid->xa0ij;

  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;
  cs_real_t *c_face_normal = coarse_grid->_face_normal;

  cs_real_t *c_xa0 = coarse_grid->_xa0;
  cs_real_t *c_xa0ij = coarse_grid->xa0ij;
  cs_real_t *c_da = coarse_grid->_da;
  cs_real_t *c_xa = coarse_grid->_xa;

  cs_real_t *w1 = NULL;

  const int *db_size = fine_grid->diag_block_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_face_normal = fine_grid->face_normal;
  const cs_real_t *f_xa0 = fine_grid->xa0;
  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_xa = fine_grid->xa;

  BFT_MALLOC(w1, f_n_cells_ext*db_size[3], cs_real_t);

  if (fine_grid->symmetric == true)
    isym = 1;

  /* P0 restriction of matrixes, "interior" surface: */
  /* xag0(nfacg), surfag(3,nfacgl), xagxg0(2,nfacg) */

# pragma omp parallel for if(c_n_faces*6 > CS_THR_MIN)
  for (c_face = 0; c_face < c_n_faces; c_face++) {
    c_xa0[c_face] = 0.;
    c_face_normal[3*c_face]    = 0.;
    c_face_normal[3*c_face +1] = 0.;
    c_face_normal[3*c_face +2] = 0.;
    c_xa0ij[3*c_face]    = 0.;
    c_xa0ij[3*c_face +1] = 0.;
    c_xa0ij[3*c_face +2] = 0.;
  }

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    if (c_coarse_face[face_id] > 0 ) {
      c_face = c_coarse_face[face_id] -1;

      c_xa0[c_face] += f_xa0[face_id];
      c_face_normal[3*c_face]    += f_face_normal[3*face_id];
      c_face_normal[3*c_face +1] += f_face_normal[3*face_id +1];
      c_face_normal[3*c_face +2] += f_face_normal[3*face_id +2];
      c_xa0ij[3*c_face]    += f_xa0ij[3*face_id];
      c_xa0ij[3*c_face +1] += f_xa0ij[3*face_id +1];
      c_xa0ij[3*c_face +2] += f_xa0ij[3*face_id +2];
    }
    else if (c_coarse_face[face_id] < 0) {
      c_face = -c_coarse_face[face_id] -1;

      c_xa0[c_face] += f_xa0[face_id];
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

  interp = 1;

  /* Initialize non differential fine grid term saved in w1 */

  if (db_size[0] == 1) {
#   pragma omp parallel for if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++)
      w1[ii] = f_da[ii];
  }
  else {
#   pragma omp parallel for private(jj, kk) if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      for (jj = 0; jj < db_size[0]; jj++) {
        for (kk = 0; kk < db_size[0]; kk++)
          w1[ii*db_size[3] + db_size[2]*jj + kk]
            = f_da[ii*db_size[3] + db_size[2]*jj + kk];
      }
    }
  }
# pragma omp parallel for if(f_n_cells_ext - f_n_cells > CS_THR_MIN)
  for (ii = f_n_cells*db_size[3]; ii < f_n_cells_ext*db_size[3]; ii++)
    w1[ii] = 0.;

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cell[face_id][0];
    jj = f_face_cell[face_id][1];
    for (kk = 0; kk < db_size[0]; kk++) {
      w1[ii*db_size[3] + db_size[2]*kk + kk] += f_xa[face_id*isym];
      w1[jj*db_size[3] + db_size[2]*kk + kk] += f_xa[(face_id +1)*isym -1];
    }
  }

  /* Initialize coarse matrix storage on (c_da, c_xa) */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (ic = 0; ic < c_n_cells_ext*db_size[3]; ic++)
    c_da[ic] = 0.;

  /* Extradiagonal terms */
  /* (symmetric matrixes for now, even with non symmetric storage isym = 2) */

  /* Matrix initialized to c_xa0 (interp=0) */

  if (isym == 1) {
#   pragma omp parallel for if(c_n_faces > CS_THR_MIN)
    for (c_face = 0; c_face < c_n_faces; c_face++)
      c_xa[c_face] = c_xa0[c_face];
  }
  else {
#   pragma omp parallel for if(c_n_faces*2 > CS_THR_MIN)
    for (c_face = 0; c_face < c_n_faces; c_face++) {
      c_xa[2*c_face]    = c_xa0[c_face];
      c_xa[2*c_face +1] = c_xa0[c_face];
    }
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
      for (c_face = 0; c_face < c_n_faces; c_face++)
        c_xa[c_face] = relax_param*c_xa[c_face]
                     + (1. - relax_param)*c_xa0[c_face];
    }
    else {
      for (c_face = 0; c_face < c_n_faces; c_face++) {
        c_xa[2*c_face]    =         relax_param  * c_xa[2*c_face]
                            + (1. - relax_param) * c_xa0[c_face];
        c_xa[2*c_face +1] =         relax_param  * c_xa[2*c_face +1]
                            + (1. - relax_param) * c_xa0[c_face];
      }
    }

  } /* endif interp == 1 */

  if (interp != 0 && interp != 1)
    bft_error(__FILE__, __LINE__, 0, "interp incorrectly defined.");

  /* Diagonal term */

  if (db_size[0] == 1) {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      c_da[ic] += w1[ii];
    }
  }
  else {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      for (jj = 0; jj < db_size[0]; jj++) {
        for (kk = 0; kk < db_size[0]; kk++)
          c_da[ic*db_size[3] + db_size[2]*jj + kk]
            += w1[ii*db_size[3] + db_size[2]*jj + kk];
      }
    }
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    for (kk = 0; kk < db_size[0]; kk++) {
      c_da[ic*db_size[3] + db_size[2]*kk + kk] -= c_xa[c_face*isym];
      c_da[jc*db_size[3] + db_size[2]*kk + kk] -= c_xa[(c_face +1)*isym -1];
    }
  }

  BFT_FREE(w1);

  /* Optional verification */

  if (verbosity > 3)
    _verify_coarse_quantities(fine_grid,
                              coarse_grid,
                              n_clips_min,
                              n_clips_max,
                              interp);

#if 0
  /* Experimental: force diagonal dominance on coarse grid */
  {
    cs_lnum_t c_n_cells = coarse_grid->n_cells;

    BFT_MALLOC(w1, c_n_cells_ext, cs_real_t);

#   pragma omp parallel for if(c_n_cells > CS_THR_MIN)
    for (ic = 0; ic < c_n_cells; ic++)
      w1[ic] = CS_ABS(c_da[ic]);

#   pragma omp parallel for if(c_n_cells_ext - c_n_cells > CS_THR_MIN)
    for (ic = c_n_cells; ic < c_n_cells_ext; ic++)
      w1[ic] = 0.;

    for (c_face = 0; c_face < c_n_faces; c_face++) {
      ic = c_face_cell[c_face][0];
      jc = c_face_cell[c_face][1];
      w1[ic] -= fabs(c_xa[c_face*isym]);
      w1[jc] -= fabs(c_xa[(c_face +1)*isym -1]);
    }

    for (ic = 0; ic < c_n_cells; ic++) {
      if (w1[ic] < 0.) {
        if (c_da[ic] <= 0.)
          c_da[ic] += w1[ic];
        else
          c_da[ic] -= w1[ic];
      }
    }

    BFT_FREE(w1);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Build a coarse level from a finer level for convection/diffusion case.
 *
 * parameters:
 *   fine_grid   <-- Fine grid structure
 *   coarse_grid <-> Coarse grid structure
 *   relax_param <-- P0/P1 relaxation parameter
 *   verbosity   <-- verbosity level
 *----------------------------------------------------------------------------*/

static void
_compute_coarse_quantities_conv_diff(const cs_grid_t  *fine_grid,
                                     cs_grid_t        *coarse_grid,
                                     cs_real_t         relax_param,
                                     int               verbosity)
{
  int interp;

  cs_lnum_t ic, jc, ii, jj, kk, c_face, face_id;

  cs_real_t dsigjg, dsxaij, agij;

  cs_lnum_t f_n_cells = fine_grid->n_cells;
  cs_lnum_t f_n_cells_ext = fine_grid->n_cells_ext;
  cs_lnum_t f_n_faces = fine_grid->n_faces;

  cs_lnum_t c_n_cells_ext = coarse_grid->n_cells_ext;
  cs_lnum_t c_n_faces = coarse_grid->n_faces;

  cs_lnum_t *c_coarse_cell = coarse_grid->coarse_cell;
  cs_lnum_t *c_coarse_face = coarse_grid->coarse_face;

  cs_gnum_t n_clips_min = 0, n_clips_max = 0;

  cs_real_t *f_xa0ij = fine_grid->xa0ij;

  cs_real_t *c_cell_cen = coarse_grid->_cell_cen;
  cs_real_t *c_face_normal = coarse_grid->_face_normal;

  cs_real_t *c_xa0 = coarse_grid->_xa0;
  cs_real_t *c_xa0_diff = coarse_grid->_xa0_diff;
  cs_real_t *c_xa0ij = coarse_grid->xa0ij;
  cs_real_t *c_da = coarse_grid->_da;
  cs_real_t *c_xa = coarse_grid->_xa;
  cs_real_t *c_da_conv = coarse_grid->_da_conv;
  cs_real_t *c_xa_conv = coarse_grid->_xa_conv;
  cs_real_t *c_da_diff = coarse_grid->_da_diff;
  cs_real_t *c_xa_diff = coarse_grid->_xa_diff;

  cs_real_t *w1 = NULL, *w1_conv = NULL, *w1_diff = NULL;

  const int *db_size = fine_grid->diag_block_size;

  const cs_lnum_2_t *f_face_cell = fine_grid->face_cell;
  const cs_lnum_2_t *c_face_cell = coarse_grid->face_cell;

  const cs_real_t *f_face_normal = fine_grid->face_normal;
  const cs_real_t *f_xa0 = fine_grid->xa0;
  const cs_real_t *f_da = fine_grid->da;
  const cs_real_t *f_da_conv = fine_grid->da_conv;
  const cs_real_t *f_xa_conv = fine_grid->xa_conv;
  const cs_real_t *f_xa0_diff = fine_grid->xa0_diff;
  const cs_real_t *f_da_diff = fine_grid->da_diff;
  const cs_real_t *f_xa_diff = fine_grid->xa_diff;

  BFT_MALLOC(w1,      2*f_n_cells_ext*db_size[3], cs_real_t);
  BFT_MALLOC(w1_conv, 2*f_n_cells_ext*db_size[3], cs_real_t);
  BFT_MALLOC(w1_diff, 2*f_n_cells_ext*db_size[3], cs_real_t);

  assert(fine_grid->symmetric == false);

  /* P0 restriction of matrixes, "interior" surface: */
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
      c_xa_conv[2*c_face]     += f_xa_conv[2*face_id];
      c_xa_conv[2*c_face +1]  += f_xa_conv[2*face_id +1];
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
      c_xa_conv[2*c_face]     += f_xa_conv[2*face_id +1];
      c_xa_conv[2*c_face +1]  += f_xa_conv[2*face_id];
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

  interp = 1;

  /* Initialize non differential fine grid term saved in w1 */

  if (db_size[0] == 1) {
#   pragma omp parallel for if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      w1_conv[ii] = f_da_conv[ii];
      w1_diff[ii] = f_da_diff[ii];
      w1[ii]      = f_da[ii] - f_da_conv[ii] - f_da_diff[ii];
    }
  }
  else {
#   pragma omp parallel for private(jj, kk) if(f_n_cells > CS_THR_MIN)
    for (ii = 0; ii < f_n_cells; ii++) {
      for (jj = 0; jj < db_size[0]; jj++) {
        for (kk = 0; kk < db_size[0]; kk++) {
          w1_conv[ii*db_size[3] + db_size[2]*jj + kk]
            = f_da_conv[ii*db_size[3] + db_size[2]*jj + kk];
          w1_diff[ii*db_size[3] + db_size[2]*jj + kk]
            = f_da_diff[ii*db_size[3] + db_size[2]*jj + kk];
          w1[ii*db_size[3] + db_size[2]*jj + kk]
            = f_da[ii*db_size[3] + db_size[2]*jj + kk]
              - f_da_conv[ii*db_size[3] + db_size[2]*jj + kk]
              - f_da_diff[ii*db_size[3] + db_size[2]*jj + kk];
        }
      }
    }
  }

# pragma omp parallel for if(f_n_cells_ext - f_n_cells > CS_THR_MIN)
  for (ii = f_n_cells*db_size[3]; ii < f_n_cells_ext*db_size[3]; ii++)
    w1[ii] = 0.;

  if (db_size[0] == 1) {
    for (face_id = 0; face_id < f_n_faces; face_id++) {
      ii = f_face_cell[face_id][0];
      jj = f_face_cell[face_id][1];
      w1_conv[ii] += f_xa_conv[2*face_id];
      w1_conv[jj] += f_xa_conv[2*face_id +1];
      w1_diff[ii] += f_xa_diff[face_id];
      w1_diff[jj] += f_xa_diff[face_id];
    }
  }
  else {
    for (face_id = 0; face_id < f_n_faces; face_id++) {
      ii = f_face_cell[face_id][0];
      jj = f_face_cell[face_id][1];
      for (kk = 0; kk < db_size[0]; kk++) {
        w1_conv[ii*db_size[3] + db_size[2]*kk + kk]
          += f_xa_conv[2*face_id];
        w1_conv[jj*db_size[3] + db_size[2]*kk + kk]
          += f_xa_conv[2*face_id +1];
        w1_diff[ii*db_size[3] + db_size[2]*kk + kk]
          += f_xa_diff[face_id];
        w1_diff[jj*db_size[3] + db_size[2]*kk + kk]
          += f_xa_diff[face_id];
      }
    }
  }

  /* Initialize coarse matrix storage on (c_da, c_xa) */

# pragma omp parallel for if(c_n_cells_ext > CS_THR_MIN)
  for (ic = 0; ic < c_n_cells_ext*db_size[3]; ic++) {
    c_da[ic]      = 0.;
    c_da_conv[ic] = 0.;
    c_da_diff[ic] = 0.;
  }

  /* Extradiagonal terms */
  /* (symmetric matrixes for now, even with non symmetric storage isym = 2) */

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

  if (db_size[0] == 1) {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      c_da_conv[ic] += w1_conv[ii];
      c_da_diff[ic] += w1_diff[ii];
    }
  }
  else {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      for (jj = 0; jj < db_size[0]; jj++) {
        for (kk = 0; kk < db_size[0]; kk++) {
          c_da_conv[ic*db_size[3] + db_size[2]*jj + kk]
            += w1_conv[ii*db_size[3] + db_size[2]*jj + kk];
          c_da_diff[ic*db_size[3] + db_size[2]*jj + kk]
            += w1_diff[ii*db_size[3] + db_size[2]*jj + kk];
        }
      }
    }
  }

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    ic = c_face_cell[c_face][0];
    jc = c_face_cell[c_face][1];
    c_da_conv[ic] -= c_xa_conv[2*c_face];
    c_da_conv[jc] -= c_xa_conv[2*c_face +1];
    c_da_diff[ic] -= c_xa_diff[c_face];
    c_da_diff[jc] -= c_xa_diff[c_face];
  }

  /* Convection/diffusion matrix */

  for (c_face = 0; c_face < c_n_faces; c_face++) {
    c_xa[2*c_face]    = c_xa_conv[2*c_face]    + c_xa_diff[c_face];
    c_xa[2*c_face +1] = c_xa_conv[2*c_face +1] + c_xa_diff[c_face];
  }

  if (db_size[0] == 1) {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      c_da[ic] += w1[ii];
    }
  }
  else {
    for (ii = 0; ii < f_n_cells; ii++) {
      ic = c_coarse_cell[ii] -1;
      for (jj = 0; jj < db_size[0]; jj++) {
        for (kk = 0; kk < db_size[0]; kk++)
          c_da[ic*db_size[3] + db_size[2]*jj + kk]
            += w1[ii*db_size[3] + db_size[2]*jj + kk];
      }
    }
  }

  for (ic = 0; ic < c_n_cells_ext; ic++)
    c_da[ic] += c_da_conv[ic] + c_da_diff[ic];

  BFT_FREE(w1);
  BFT_FREE(w1_conv);
  BFT_FREE(w1_diff);

  /* Optional verification */

  if (verbosity > 3)
    _verify_coarse_quantities(fine_grid,
                              coarse_grid,
                              n_clips_min,
                              n_clips_max,
                              interp);

}

/*============================================================================
 * Semi-private function definitions
 *
 * The following functions are intended to be used by the multigrid layer
 * (cs_multigrid.c), not directly by the user, so they are no more
 * documented than private static functions)
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create base grid by mapping from shared mesh values.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   n_cells               <-- Local number of cells
 *   n_cells_ext           <-- Local number of cells + ghost cells
 *   n_faces               <-- Local number of faces
 *   symmetric             <-- True if xam is symmetric, false otherwise
 *   diag_block_size       <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size <-- Block sizes for diagonal, or NULL
 *   face_cell             <-- Face -> cells connectivity
 *   halo                  <-- Halo structure associated with this level,
 *                             or NULL.
 *   cell_cen              <-- Cell center (size: 3.n_cells_ext)
 *   cell_vol              <-- Cell volume (size: n_cells_ext)
 *   face_normal           <-- Internal face normals (size: 3.n_faces)
 *   a                     <-- Associated matrix
 *   a_conv                <-- Associated matrix (convection)
 *   a_diff                <-- Associated matrix (diffusion)
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_shared(cs_lnum_t              n_cells,
                           cs_lnum_t              n_cells_ext,
                           cs_lnum_t              n_faces,
                           bool                   symmetric,
                           const int             *diag_block_size,
                           const int             *extra_diag_block_size,
                           const cs_lnum_2_t     *face_cell,
                           const cs_halo_t       *halo,
                           const cs_real_t       *cell_cen,
                           const cs_real_t       *cell_vol,
                           const cs_real_t       *face_normal,
                           const cs_matrix_t     *a,
                           const cs_matrix_t     *a_conv,
                           const cs_matrix_t     *a_diff)
{
  cs_lnum_t ii, jj, kk, face_id;

  cs_grid_t *g = NULL;

  const cs_real_t *da      = cs_matrix_get_diagonal(a);
  const cs_real_t *xa      = cs_matrix_get_extra_diagonal(a);
  const cs_real_t *da_conv = NULL, *da_diff = NULL;
  const cs_real_t *xa_conv = NULL, *xa_diff = NULL;

  /* Create empty structure and map base data */

  g = _create_grid();

  g->level = 0;
  g->conv_diff = false;
  g->symmetric = symmetric;

  if (a_conv != NULL || a_diff != NULL) {
    g->conv_diff = true;
    da_conv = cs_matrix_get_diagonal(a_conv);
    da_diff = cs_matrix_get_diagonal(a_diff);
    xa_conv = cs_matrix_get_extra_diagonal(a_conv);
    xa_diff = cs_matrix_get_extra_diagonal(a_diff);
  }

  if (diag_block_size != NULL) {
    for (ii = 0; ii < 4; ii++)
      g->diag_block_size[ii] = diag_block_size[ii];
  }
  else {
    for (ii = 0; ii < 4; ii++)
      g->diag_block_size[ii] = 1;
  }

  if (extra_diag_block_size != NULL) {
    for (ii = 0; ii < 4; ii++)
      g->extra_diag_block_size[ii] = extra_diag_block_size[ii];
  }
  else {
    for (ii = 0; ii < 4; ii++)
      g->extra_diag_block_size[ii] = 1;
  }


  g->n_cells = n_cells;
  g->n_cells_ext = n_cells_ext;
  g->n_faces = n_faces;
  g->n_g_cells = n_cells;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t _n_cells = n_cells;
    MPI_Allreduce(&_n_cells, &(g->n_g_cells), 1, CS_MPI_GNUM, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  g->face_cell = face_cell;

  g->cell_cen = cell_cen;
  g->cell_vol = cell_vol;
  g->face_normal = face_normal;

  g->halo = halo;

  /* Set shared matrix coefficients */

  g->da = da;
  g->_da = NULL;

  if (g->conv_diff) {
    g->da_conv = da_conv;
    g->_da_conv = NULL;
    g->da_diff = da_diff;
    g->_da_diff = NULL;
  }

  g->xa = xa;
  g->_xa = NULL;

  if (g->conv_diff) {
    g->xa_conv = xa_conv;
    g->_xa_conv = NULL;
    g->xa_diff = xa_diff;
    g->_xa_diff = NULL;
  }

  /* Build symmetrized extra-diagonal terms if necessary,
     or point to existing terms if already symmetric */

  if (symmetric == true) {
    g->xa0 = g->xa;
    g->_xa0 = NULL;
  }
  else if (g->conv_diff) {
    g->xa0  = g->xa;
    g->_xa0 = NULL;
    g->xa0_diff = g->xa_diff;
    g->_xa0_diff = NULL;
  }
  else {
    BFT_MALLOC(g->_xa0, n_faces, cs_real_t);
    for (face_id = 0; face_id < n_faces; face_id++)
      g->_xa0[face_id] = 0.5 * (xa[face_id*2] + xa[face_id*2+1]);
    g->xa0 = g->_xa0;
  }

  /* Compute multigrid-specific terms */

  BFT_MALLOC(g->xa0ij, n_faces*3, cs_real_t);

  const cs_real_t *restrict g_xa0 = g->xa0;
  if (g->conv_diff)
    g_xa0 = g->xa0_diff;

# pragma omp parallel for private(ii, jj, kk) if(n_faces > CS_THR_MIN)
  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = face_cell[face_id][0];
    jj = face_cell[face_id][1];
    for (kk = 0; kk < 3; kk++) {
      g->xa0ij[face_id*3 + kk] =   g_xa0[face_id]
                                 * (cell_cen[jj*3 + kk] - cell_cen[ii*3 + kk]);
    }
  }

  g->matrix_struct = NULL;
  g->matrix = a;
  g->_matrix = NULL;

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

    BFT_FREE(g->_face_cell);

    BFT_FREE(g->coarse_cell);
    BFT_FREE(g->coarse_face);

    if (g->_cell_cen != NULL)
      BFT_FREE(g->_cell_cen);
    if (g->_cell_vol != NULL)
      BFT_FREE(g->_cell_vol);
    if (g->_face_normal != NULL)
      BFT_FREE(g->_face_normal);

    if (g->_halo != NULL)
      cs_halo_destroy(&(g->_halo));

    if (g->_da != NULL)
      BFT_FREE(g->_da);
    if (g->_da_conv != NULL)
      BFT_FREE(g->_da_conv);
    if (g->_da_diff != NULL)
      BFT_FREE(g->_da_diff);
    if (g->_xa != NULL)
      BFT_FREE(g->_xa);
    if (g->_xa_conv != NULL)
      BFT_FREE(g->_xa_conv);
    if (g->_xa_diff != NULL)
      BFT_FREE(g->_xa_diff);
    if (g->_xa0 != NULL)
      BFT_FREE(g->_xa0);
    if (g->_xa0_diff != NULL)
      BFT_FREE(g->_xa0_diff);

    BFT_FREE(g->xa0ij);

    cs_matrix_destroy(&(g->_matrix));
    cs_matrix_structure_destroy(&(g->matrix_struct));

#if defined(HAVE_MPI)
    BFT_FREE(g->merge_cell_idx);
#endif

    BFT_FREE(*grid);
  }
}

/*----------------------------------------------------------------------------
 * Get grid information.
 *
 * parameters:
 *   g           <-- Grid structure
 *   level       --> Level in multigrid hierarchy (or NULL)
 *   symmetric   --> Symmetric matrix coefficients indicator (or NULL)
 *   db_size     --> Size of the diagonal block (or NULL)
 *   eb_size     --> Size of the extra diagonal block (or NULL)
 *   n_ranks     --> number of ranks with data (or NULL)
 *   n_cells     --> Number of local cells (or NULL)
 *   n_cells_ext --> Number of cells including ghosts (or NULL)
 *   n_faces     --> Number of faces (or NULL)
 *   n_g_cells   --> Number of global cells (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 bool             *symmetric,
                 int              *db_size,
                 int              *eb_size,
                 int              *n_ranks,
                 cs_lnum_t        *n_cells,
                 cs_lnum_t        *n_cells_ext,
                 cs_lnum_t        *n_faces,
                 cs_gnum_t        *n_g_cells)
{
  assert(g != NULL);

  if (level != NULL)
    *level = g->level;

  if (symmetric != NULL)
    *symmetric = g->symmetric;

  if (db_size != NULL) {
    db_size[0] = g->diag_block_size[0];
    db_size[1] = g->diag_block_size[1];
    db_size[2] = g->diag_block_size[2];
    db_size[3] = g->diag_block_size[3];
  }

 if (eb_size != NULL) {
    eb_size[0] = g->extra_diag_block_size[0];
    eb_size[1] = g->extra_diag_block_size[1];
    eb_size[2] = g->extra_diag_block_size[2];
    eb_size[3] = g->extra_diag_block_size[3];
  }

  if (n_ranks != NULL) {
#if defined(HAVE_MPI)
    *n_ranks = g->n_ranks;
#else
    *n_ranks = 1;
#endif
  }

  if (n_cells != NULL)
    *n_cells = g->n_cells;
  if (n_cells_ext != NULL)
    *n_cells_ext = g->n_cells_ext;
  if (n_faces != NULL)
    *n_faces = g->n_faces;

  if (n_g_cells != NULL)
    *n_g_cells = g->n_g_cells;
}

/*----------------------------------------------------------------------------
 * Get number of cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of cells of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cells(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_cells;
}

/*----------------------------------------------------------------------------
 * Get number of extended (local + ghost) cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   number of extended cells of grid structure
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_grid_get_n_cells_ext(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_cells_ext;
}

/*----------------------------------------------------------------------------
 * Get maximum number of extended (local + ghost) cells corresponding to
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
cs_grid_get_n_cells_max(const cs_grid_t  *g)
{
  cs_lnum_t retval = 0;

  if (g != NULL)
    retval = CS_MAX(g->n_cells_ext, g->n_cells_r[1]);

  return retval;
}

/*----------------------------------------------------------------------------
 * Get global number of cells corresponding to a grid.
 *
 * parameters:
 *   g <-- Grid structure
 *
 * returns:
 *   global number of cells of grid structure
 *----------------------------------------------------------------------------*/

cs_gnum_t
cs_grid_get_n_g_cells(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_g_cells;
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
  int grid_id;

  MPI_Comm comm = cs_glob_mpi_comm;

  assert(g != NULL);

  if (g->n_ranks != cs_glob_n_ranks) {
    grid_id = 0;
    while (_grid_ranks[grid_id] != g->n_ranks && grid_id < _n_grid_comms)
      grid_id++;
    comm = _grid_comm[grid_id];
  }

  return comm;
}

#endif

/*----------------------------------------------------------------------------
 * Create coarse grid from fine grid.
 *
 * parameters:
 *   f                    <-- Fine grid structure
 *   verbosity            <-- Verbosity level
 *   coarsening_type      <-- Coarsening traversal type:
 *                              0: algebraic with natural face traversal;
 *                              1: algebraic with face traveral by criteria;
 *                              2: algebraic with Hilbert face traversal
 *   aggregation_limit    <-- Maximum allowed fine cells per coarse cell
 *   relaxation_parameter <-- P0/P1 relaxation factor
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen(const cs_grid_t   *f,
                int                verbosity,
                int                coarsening_type,
                int                aggregation_limit,
                double             relaxation_parameter)
{
  cs_lnum_t isym = 2;
  bool conv_diff = f->conv_diff;

  /* By default, always use MSR structure, as it often seems to provide the
     best performance, and is required for the hybrid Gauss-Seidel-Jacobi
     smoothers. In multithreaded case, we also prefer to use a matrix
     structure allowing threading without a specific renumbering, as
     structures are rebuilt often (so only CSR and MSR can be considered) */

  cs_matrix_type_t coarse_matrix_type = CS_MATRIX_MSR;
  cs_matrix_variant_t *coarse_mv = NULL;

  cs_grid_t *c = NULL;

  const int *db_size = f->diag_block_size;

  assert(f != NULL);

  /* Initialization */

  c = _coarse_init(f);

  if (f->symmetric == true)
    isym = 1;

  /* Determine fine->coarse cell connectivity (aggregation) */

  _automatic_aggregation(f,
                         coarsening_type,
                         aggregation_limit,
                         relaxation_parameter,
                         verbosity,
                         c->coarse_cell);

  /* Build coarse grid connectivity */

  _coarsen(f, c);

  /* Allocate permanent arrays in coarse grid */

  BFT_MALLOC(c->_cell_cen, c->n_cells_ext*3, cs_real_t);
  c->cell_cen = c->_cell_cen;

  BFT_MALLOC(c->_cell_vol, c->n_cells_ext, cs_real_t);
  c->cell_vol = c->_cell_vol;

  BFT_MALLOC(c->_face_normal, c->n_faces*3, cs_real_t);
  c->face_normal = c->_face_normal;

  BFT_MALLOC(c->_da, c->n_cells_ext * c->diag_block_size[3], cs_real_t);
  c->da = c->_da;

  if (conv_diff) {
    BFT_MALLOC(c->_da_conv, c->n_cells_ext * c->diag_block_size[3], cs_real_t);
    c->da_conv = c->_da_conv;
    BFT_MALLOC(c->_da_diff, c->n_cells_ext * c->diag_block_size[3], cs_real_t);
    c->da_diff = c->_da_diff;
  }

  BFT_MALLOC(c->_xa, c->n_faces*isym, cs_real_t);
  c->xa = c->_xa;

  if (conv_diff) {
    BFT_MALLOC(c->_xa_conv, c->n_faces*2, cs_real_t);
    c->xa_conv = c->_xa_conv;
    BFT_MALLOC(c->_xa_diff, c->n_faces, cs_real_t);
    c->xa_diff = c->_xa_diff;
  }

  /* We could have xa0 point to xa if symmetric, but this would require
     caution in CRSTGR to avoid overwriting. */

  BFT_MALLOC(c->_xa0, c->n_faces*isym, cs_real_t);
  c->xa0 = c->_xa0;

  if (conv_diff) {
    BFT_MALLOC(c->_xa0_diff, c->n_faces, cs_real_t);
    c->xa0_diff = c->_xa0_diff;
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

  if (conv_diff)
    _compute_coarse_quantities_conv_diff(f, c, relaxation_parameter, verbosity);
  else
    _compute_coarse_quantities(f, c, relaxation_parameter, verbosity);

  /* Synchronize matrix's geometric quantities */

  if (c->halo != NULL)
    cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_da, db_size[3]);

  /* Merge grids if we are below the threshold */

#if defined(HAVE_MPI)
  if (c->n_ranks > _grid_merge_min_ranks && _grid_merge_stride > 1) {
    cs_gnum_t  _n_ranks = c->n_ranks;
    cs_gnum_t  _n_mean_g_cells = c->n_g_cells / _n_ranks;
    if (   _n_mean_g_cells < _grid_merge_mean_threshold
        || c->n_g_cells < _grid_merge_glob_threshold)
    _merge_grids(c, verbosity);
  }
#endif

  if (_grid_tune_max_level > 0) {

    cs_matrix_fill_type_t mft
      = cs_matrix_get_fill_type(f->symmetric,
                                f->diag_block_size,
                                f->extra_diag_block_size);

    if (_grid_tune_max_level > f->level) {
      int k = CS_MATRIX_N_FILL_TYPES*(f->level) + mft;
      coarse_mv = _grid_tune_variant[k];

      /* Create tuned variant upon first pass for this level and
         fill type */

      if  (   coarse_mv == NULL
           && _grid_tune_max_fill_level[mft] > f->level) {

        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("\n"
                        "Tuning for coarse matrices of level %d and type: %s\n"
                        "==========================\n"),
                      f->level + 1, cs_matrix_fill_type_name[mft]);

        int n_min_products;
        double t_measure;
        int n_fill_types = 1;
        cs_matrix_fill_type_t fill_types[1] = {mft};
        double fill_weights[1] = {1};

        cs_matrix_get_tuning_runs(&n_min_products, &t_measure);

        coarse_mv = cs_matrix_variant_tuned(t_measure,
                                            0, /* n_matrix_types, */
                                            n_fill_types,
                                            NULL, /* matrix_types, */
                                            fill_types,
                                            fill_weights,
                                            n_min_products,
                                            c->n_cells,
                                            c->n_cells_ext,
                                            c->n_faces,
                                            c->face_cell,
                                            c->halo,
                                            NULL);  /* face_numbering */

        _grid_tune_variant[k] = coarse_mv;

        if  (_grid_tune_max_fill_level[mft] == f->level + 1) {
          cs_log_printf(CS_LOG_PERFORMANCE, "\n");
          cs_log_separator(CS_LOG_PERFORMANCE);
        }
      }

    }

  }

  if (coarse_mv != NULL)
    coarse_matrix_type = cs_matrix_variant_type(coarse_mv);

  c->matrix_struct = cs_matrix_structure_create(coarse_matrix_type,
                                                true,
                                                c->n_cells,
                                                c->n_cells_ext,
                                                c->n_faces,
                                                c->face_cell,
                                                c->halo,
                                                NULL);

  if (coarse_mv != NULL)
    c->_matrix = cs_matrix_create_by_variant(c->matrix_struct, coarse_mv);
  else
    c->_matrix = cs_matrix_create(c->matrix_struct);

  cs_matrix_set_coefficients(c->_matrix,
                             c->symmetric,
                             c->diag_block_size,
                             c->extra_diag_block_size,
                             c->n_faces,
                             c->face_cell,
                             c->da,
                             c->xa);

  c->matrix = c->_matrix;

  /* Return new (coarse) grid */

  return c;
}

/*----------------------------------------------------------------------------
 * Compute coarse cell variable values from fine cell values
 *
 * parameters:
 *   f       <-- Fine grid structure
 *   c       <-- Fine grid structure
 *   f_var   <-- Variable defined on fine grid cells
 *   c_var   --> Variable defined on coarse grid cells
 *----------------------------------------------------------------------------*/

void
cs_grid_restrict_cell_var(const cs_grid_t  *f,
                          const cs_grid_t  *c,
                          const cs_real_t  *f_var,
                          cs_real_t        *c_var)
{
  cs_lnum_t ii;
  int i;

  cs_lnum_t f_n_cells = f->n_cells;
  cs_lnum_t c_n_cells_ext = c->n_cells_r[1];

  const cs_lnum_t *coarse_cell;
  const int *db_size = f->diag_block_size;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL || f_n_cells == 0);
  assert(f_var != NULL || f_n_cells == 0);
  assert(c_var != NULL || c_n_cells_ext == 0);

  /* Set coarse values */

  coarse_cell = c->coarse_cell;

# pragma omp parallel for private(i) if(c_n_cells_ext > CS_THR_MIN)
  for (ii = 0; ii < c_n_cells_ext; ii++)
    for (i = 0; i < db_size[0]; i++)
      c_var[ii*db_size[1]+i] = 0.;

  for (ii = 0; ii < f_n_cells; ii++)
    for (i = 0; i < db_size[0]; i++)
      c_var[(coarse_cell[ii]-1)*db_size[1]+i] += f_var[ii*db_size[1]+i];

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
        MPI_Recv(c_var + c->merge_cell_idx[rank_id]*db_size[1],
                 n_recv*db_size[1], CS_MPI_REAL, dist_rank, tag, comm, &status);
      }
    }
    else
      MPI_Send(c_var, c->n_cells_r[0]*db_size[1], CS_MPI_REAL,
               c->merge_sub_root, tag, comm);
  }

#endif /* defined(HAVE_MPI) */
}

/*----------------------------------------------------------------------------
 * Compute fine cell integer values from coarse cell values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_num   --> Variable defined on coarse grid cells
 *   f_num   <-- Variable defined on fine grid cells
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_cell_num(const cs_grid_t  *c,
                         const cs_grid_t  *f,
                         int              *c_num,
                         int              *f_num)
{
  cs_lnum_t ii;
  const cs_lnum_t *coarse_cell;
  const int *_c_num = c_num;

  cs_lnum_t f_n_cells = f->n_cells;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL || f_n_cells == 0);
  assert(f_num != NULL);
  assert(c_num != NULL);

#if defined(HAVE_MPI)
  _scatter_cell_num(c, c_num);
#endif

  /* Set fine values */

  coarse_cell = c->coarse_cell;

# pragma omp parallel for if(f_n_cells > CS_THR_MIN)
  for (ii = 0; ii < f_n_cells; ii++)
    f_num[ii] = _c_num[coarse_cell[ii] - 1];
}

/*----------------------------------------------------------------------------
 * Compute fine cell variable values from coarse cell values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_var   --> Variable defined on coarse grid cells
 *   f_var   <-- Variable defined on fine grid cells
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_cell_var(const cs_grid_t  *c,
                         const cs_grid_t  *f,
                         cs_real_t        *c_var,
                         cs_real_t        *f_var)
{
  cs_lnum_t ii;
  int i;
  const cs_lnum_t *coarse_cell;
  const cs_real_t *_c_var = c_var;

  const int *db_size = f->diag_block_size;

  cs_lnum_t f_n_cells = f->n_cells;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL || f_n_cells == 0);
  assert(f_var != NULL);
  assert(c_var != NULL);

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
        MPI_Send(c_var + c->merge_cell_idx[rank_id]*db_size[1],
                 n_send*db_size[1], CS_MPI_REAL, dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(c_var, c->n_cells_r[0]*db_size[1], CS_MPI_REAL,
               c->merge_sub_root, tag, comm, &status);
    }
  }

#endif /* defined(HAVE_MPI) */

  /* Set fine values */

  coarse_cell = c->coarse_cell;

# pragma omp parallel for private(i) if(f_n_cells > CS_THR_MIN)
  for (ii = 0; ii < f_n_cells; ii++)
    for (i = 0; i < db_size[0]; i++)
      f_var[ii*db_size[1]+i] = _c_var[(coarse_cell[ii]-1)*db_size[1]+i];
}

/*----------------------------------------------------------------------------
 * Project coarse grid cell numbers to base grid.
 *
 * If a global coarse grid cell number is larger than max_num, its
 * value modulo max_num is used.
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   max_num      <-- Values of c_cell_num = global_num % max_num
 *   c_cell_num   --> Global coarse cell number (modulo max_num)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_cell_num(const cs_grid_t  *g,
                         cs_lnum_t         n_base_cells,
                         int               max_num,
                         int               c_cell_num[])
{
  cs_lnum_t ii = 0;
  cs_gnum_t base_shift = 1;
  cs_gnum_t _max_num = max_num;
  cs_lnum_t n_max_cells = 0;
  cs_lnum_t *tmp_num_1 = NULL, *tmp_num_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(c_cell_num != NULL);

  /* Initialize array */

  n_max_cells = g->n_cells;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_cells > n_max_cells)
      n_max_cells = _g->n_cells;
  }

  BFT_MALLOC(tmp_num_1, n_max_cells, cs_lnum_t);

  /* Compute local base starting cell number in parallel mode */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_gnum_t local_shift = g->n_cells;
    cs_gnum_t global_shift = 0;
    MPI_Scan(&local_shift, &global_shift, 1, CS_MPI_GNUM, MPI_SUM,
             cs_glob_mpi_comm);
    base_shift = 1 + global_shift - g->n_cells;
  }
#endif

  for (ii = 0; ii < g->n_cells; ii++)
    tmp_num_1[ii] = (cs_gnum_t)(ii + base_shift) % _max_num;

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_num_2, n_max_cells, cs_lnum_t);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_cells = _g->parent->n_cells;

#if defined(HAVE_MPI)
      _scatter_cell_num(_g, tmp_num_1);
#endif /* defined(HAVE_MPI) */

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_num_2[ii] = tmp_num_1[_g->coarse_cell[ii] - 1];

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_num_1[ii] = tmp_num_2[ii];

    }

    assert(_g->level == 0);
    assert(_g->n_cells == n_base_cells);

    /* Free temporary arrays */

    BFT_FREE(tmp_num_2);
  }

  memcpy(c_cell_num, tmp_num_1, n_base_cells*sizeof(int));

  BFT_FREE(tmp_num_1);
}

/*----------------------------------------------------------------------------
 * Project coarse grid cell rank to base grid.
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   c_cell_rank  --> Global coarse cell rank projected to fine cells
 *----------------------------------------------------------------------------*/

void
cs_grid_project_cell_rank(const cs_grid_t  *g,
                          cs_lnum_t         n_base_cells,
                          int               c_cell_rank[])
{
  cs_lnum_t ii;
  cs_lnum_t n_max_cells = 0;
  int *tmp_rank_1 = NULL, *tmp_rank_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(c_cell_rank != NULL || g->n_cells == 0);

  n_max_cells = g->n_cells;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_cells > n_max_cells)
      n_max_cells = _g->n_cells;
  }

  BFT_MALLOC(tmp_rank_1, n_max_cells, int);

  for (ii = 0; ii < g->n_cells; ii++)
    tmp_rank_1[ii] = cs_glob_rank_id;

  /* Project to finer levels */

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_rank_2, n_max_cells, int);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_cells = _g->parent->n_cells;

      cs_grid_prolong_cell_num(_g,
                               _g->parent,
                               tmp_rank_1,
                               tmp_rank_2);

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_rank_1[ii] = tmp_rank_2[ii];

    }

    assert(_g->level == 0);
    assert(_g->n_cells == n_base_cells);

    /* Free temporary arrays */

    BFT_FREE(tmp_rank_2);
  }

  memcpy(c_cell_rank, tmp_rank_1, n_base_cells*sizeof(int));

  BFT_FREE(tmp_rank_1);
}

/*----------------------------------------------------------------------------
 * Project variable from coarse grid to base grid
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   c_var        <-- Cell variable on coarse grid
 *   f_var        --> Cell variable projected to fine grid
 *----------------------------------------------------------------------------*/

void
cs_grid_project_var(const cs_grid_t  *g,
                    cs_lnum_t         n_base_cells,
                    const cs_real_t   c_var[],
                    cs_real_t         f_var[])
{
  cs_lnum_t ii;
  int i;
  cs_lnum_t n_max_cells = 0;
  cs_real_t *tmp_var_1 = NULL, *tmp_var_2 = NULL;
  const cs_grid_t *_g = g;

  const int *db_size = g->diag_block_size;

  assert(g != NULL);
  assert(c_var != NULL || g->n_cells == 0);
  assert(f_var != NULL);

  n_max_cells = g->n_cells;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_cells > n_max_cells)
      n_max_cells = _g->n_cells;
  }

  BFT_MALLOC(tmp_var_1, n_max_cells*db_size[1], cs_real_t);
  memcpy(tmp_var_1, c_var, g->n_cells*db_size[1]*sizeof(cs_real_t));

  /* Project to finer levels */

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_var_2, n_max_cells*db_size[1], cs_real_t);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      cs_lnum_t n_parent_cells = _g->parent->n_cells;

      cs_grid_prolong_cell_var(_g,
                               _g->parent,
                               tmp_var_1,
                               tmp_var_2);

      for (ii = 0; ii < n_parent_cells; ii++)
        for (i = 0; i < db_size[0]; i++)
          tmp_var_1[ii*db_size[1]+i] = tmp_var_2[ii*db_size[1]+i];

    }

    assert(_g->level == 0);
    assert(_g->n_cells == n_base_cells);

    /* Free temporary arrays */

    BFT_FREE(tmp_var_2);
  }

  memcpy(f_var, tmp_var_1, n_base_cells*db_size[1]*sizeof(cs_real_t));

  BFT_FREE(tmp_var_1);
}

/*----------------------------------------------------------------------------
 * Compute diagonal dominance metric and project it to base grid
 *
 * parameters:
 *   g            <-- Grid structure
 *   n_base_cells <-- Number of cells in base grid
 *   diag_dom     --> Diagonal dominance metric (on fine grid)
 *----------------------------------------------------------------------------*/

void
cs_grid_project_diag_dom(const cs_grid_t  *g,
                         cs_lnum_t         n_base_cells,
                         cs_real_t         diag_dom[])
{
  cs_lnum_t ii, jj, face_id;
  int i, j;

  cs_real_t *dd = NULL;
  const int *db_size = g->diag_block_size;

  assert(g != NULL);
  assert(diag_dom != NULL);

  if (g->level == 0)
    dd = diag_dom;
  else
    BFT_MALLOC(dd, g->n_cells_ext*db_size[3], cs_real_t);

  /* Compute coarse diagonal dominance */

  {
    const cs_lnum_t n_cells = g->n_cells;
    const cs_lnum_t n_faces = g->n_faces;
    const cs_lnum_2_t *face_cel = g->face_cell;
    cs_real_t abs_trace;

    /* Diagonal part of matrix.vector product */

    for (ii = 0; ii < n_cells; ii++)
      for (i = 0; i < db_size[0]; i++)
        for (j = 0; j < db_size[0]; j++)
          dd[ii*db_size[3] + db_size[2]*i + j] =
            fabs(g->da[ii*db_size[3] + db_size[2]*i + j]);

    if (g->halo != NULL)
      cs_halo_sync_var_strided(g->halo, CS_HALO_STANDARD, dd, db_size[3]);

    if (g->symmetric) {
      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = face_cel[face_id][0];
        jj = face_cel[face_id][1];

        for (i = 0; i < db_size[0]; i++) {
          dd[ii*db_size[3] + db_size[2]*i + i] -= fabs(g->xa[face_id]);
          dd[jj*db_size[3] + db_size[2]*i + i] -= fabs(g->xa[face_id]);
        }
      }
    }
    else {
      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = face_cel[face_id][0];
        jj = face_cel[face_id][1];

        for (i = 0; i < db_size[0]; i++) {
          dd[ii*db_size[3] + db_size[2]*i + i] -= fabs(g->xa[face_id*2]);
          dd[jj*db_size[3] + db_size[2]*i + i] -= fabs(g->xa[face_id*2 + 1]);
        }
      }
    }

    for (ii = 0; ii < n_cells; ii++) {
      abs_trace = 0.;
      for (i = 0; i < db_size[0]; i++)
        abs_trace += g->da[ii*db_size[3] + db_size[2]*i + i];

      abs_trace = fabs(abs_trace);
      if (abs_trace > 1.e-18)
        for (i = 0; i < db_size[0]; i++)
          for (j = 0; j < db_size[0]; j++)
            dd[ii*db_size[3] + db_size[2]*i + j] /= abs_trace;
    }

  }

  /* Now project to finer levels */

  if (dd != diag_dom) {
    cs_grid_project_var(g, n_base_cells, dd, diag_dom);
    BFT_FREE(dd);
  }
}

/*----------------------------------------------------------------------------
 * Return the merge_stride if merging is active.
 *
 * returns:
 *   grid merge stride if merging is active, 1 otherwise
 *----------------------------------------------------------------------------*/

int
cs_grid_get_merge_stride(void)
{
  int retval = 1;

#if defined(HAVE_MPI)
  if (_grid_merge_min_ranks < cs_glob_n_ranks && _grid_merge_stride > 1)
    retval = _grid_merge_stride;
#endif

  return retval;
}

/*----------------------------------------------------------------------------
 * Finalize global info related to multigrid solvers
 *----------------------------------------------------------------------------*/

void
cs_grid_finalize(void)
{
#if defined(HAVE_MPI)
  _finalize_reduced_communicators();
#endif

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
             "  grid:           %p\n"
             "  level:          %d (parent: %p)\n"
             "  n_cells:        %d\n"
             "  n_cells_ext:    %d\n"
             "  n_faces:        %d\n"
             "  n_g_cells:      %d\n"
             "  n_cells_r:      [%d, %d]\n",
             (const void *)g, g->level, (const void *)(g->parent),
             (int)(g->n_cells), (int)(g->n_cells_ext),
             (int)(g->n_faces), (int)(g->n_g_cells),
             (int)(g->n_cells_r[0]), (int)(g->n_cells_r[1]));

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
      bft_printf("    %d: %d\n", i, g->merge_cell_idx[i]);
  }

#endif
  bft_printf("\n"
             "  face_cell:      %p\n"
             "  _face_cell:     %p\n"
             "  coarse_cell:    %p\n"
             "  coarse_face:    %p\n"
             "  halo:           %p\n",
             (const void *)g->face_cell, (const void *)g->_face_cell,
             (const void *)g->coarse_cell, (const void *)g->coarse_face,
             (const void *)g->halo);

  if (g->face_cell != NULL) {
    bft_printf("\n"
               "  face -> cell connectivity;\n");
    for (i = 0; i < g->n_faces; i++)
      bft_printf("    %d : %d, %d\n", (int)(i+1),
                 (int)(g->face_cell[i][0]), (int)(g->face_cell[i][1]));
  }

  if (g->coarse_cell != NULL && g->parent != NULL) {
    bft_printf("\n"
               "  coarse_cell;\n");
    for (i = 0; i < g->parent->n_cells; i++)
      bft_printf("    %d : %d\n",
                 (int)(i+1), (int)(g->coarse_cell[i]));
  }

  if (g->coarse_face != NULL && g->parent != NULL) {
    bft_printf("\n"
               "  coarse_face;\n");
    for (i = 0; i < g->parent->n_faces; i++)
      bft_printf("    %d : %d\n",
                 (int)(i+1), (int)(g->coarse_face[i]));
  }

  cs_halo_dump(g->halo, 1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query the global multigrid parameters for parallel grid merging.
 *
 * \param[out]  rank_stride           number of ranks over which merging
 *                                    takes place, or NULL
 * \param[out]  cells_mean_threshold  mean number of cells under which merging
 *                                    should be applied, or NULL
 * \param[out]  cells_glob_threshold  global number of cells under which merging
 *                                    should be applied, or NULL
 * \param[out]  min_ranks             number of active ranks under which
 *                                    no merging takes place, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_get_merge_options(int         *rank_stride,
                          int         *cells_mean_threshold,
                          cs_gnum_t   *cells_glob_threshold,
                          int         *min_ranks)
{
#if defined(HAVE_MPI)
  if (rank_stride != NULL)
    *rank_stride = _grid_merge_stride;
  if (cells_mean_threshold != NULL)
    *cells_mean_threshold = _grid_merge_mean_threshold;
  if (cells_glob_threshold != NULL)
    *cells_glob_threshold = _grid_merge_glob_threshold;
  if (min_ranks != NULL)
    *min_ranks = _grid_merge_min_ranks;
#else
  if (rank_stride != NULL)
    *rank_stride = 0;
  if (cells_mean_threshold != NULL)
    *cells_mean_threshold = 0;
  if (cells_glob_threshold != NULL)
    *cells_glob_threshold = 0;
  if (min_ranks != NULL)
    *min_ranks = 1;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set global multigrid parameters for parallel grid merging behavior.
 *
 * \param[in]  rank_stride           number of ranks over which merging
 *                                   takes place
 * \param[in]  cells_mean_threshold  mean number of cells under which merging
 *                                   should be applied
 * \param[in]  cells_glob_threshold  global number of cells under which merging
 *                                   should be applied
 * \param[in]  min_ranks             number of active ranks under which
 *                                   no merging takes place
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_set_merge_options(int         rank_stride,
                          int         cells_mean_threshold,
                          cs_gnum_t   cells_glob_threshold,
                          int         min_ranks)
{
#if defined(HAVE_MPI)
  _grid_merge_stride = rank_stride;
  if (_grid_merge_stride > cs_glob_n_ranks)
    _grid_merge_stride = cs_glob_n_ranks;
  _grid_merge_mean_threshold = cells_mean_threshold;
  _grid_merge_glob_threshold = cells_glob_threshold;
  _grid_merge_min_ranks = min_ranks;
#endif
}

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
/*!
 * \brief Force matrix variant selection for multigrid coarse meshes.
 *
 * The finest mesh (level 0) is handled by the default tuning options,
 * so only coarser meshes are considered here.
 *
 * \param[in]  fill_type  associated matrix fill type
 * \param[in]  level      level for which variant is assiged
 * \param[in]  mv         matrix variant to assign (NULL to unassign)
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_set_matrix_variant(cs_matrix_fill_type_t       fill_type,
                           int                         level,
                           const cs_matrix_variant_t  *mv)
{
  if (_grid_tune_max_level < level) {

    if (_grid_tune_max_level == 0) {
      BFT_MALLOC(_grid_tune_max_fill_level, CS_MATRIX_N_FILL_TYPES, int);
      for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++)
        _grid_tune_max_fill_level[i] = 0;
    }

    BFT_REALLOC(_grid_tune_variant,
                CS_MATRIX_N_FILL_TYPES*level, cs_matrix_variant_t *);

    for (int i = _grid_tune_max_level; i < level; i++) {
      for (int j = 0; j < CS_MATRIX_N_FILL_TYPES; j++) {
        _grid_tune_variant[CS_MATRIX_N_FILL_TYPES*i + j] = NULL;
      }
    }

    _grid_tune_max_level = level;
  }

  int k = CS_MATRIX_N_FILL_TYPES*(level-1) + fill_type;

  if (_grid_tune_variant[k] != NULL)
    cs_matrix_variant_destroy(&(_grid_tune_variant[k]));

  if (mv != NULL) {
    cs_matrix_type_t m_type = cs_matrix_variant_type(mv);
    _grid_tune_variant[k] = cs_matrix_variant_create(m_type, NULL);
    cs_matrix_variant_merge(_grid_tune_variant[k], mv, fill_type);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the current settings for multigrid parallel merging.
 */
/*----------------------------------------------------------------------------*/

void
cs_grid_log_merge_options(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    cs_log_printf(CS_LOG_SETUP,
                  _("\n"
                    "Multigrid rank merge parameters:\n"
                    "  merge rank stride:                  %d\n"
                    "  mean  coarse cells merge threshold: %d\n"
                    "  total coarse cells merge threshold: %llu\n"
                    "  minimum active ranks:               %d\n"),
                  _grid_merge_stride,
                  (int)_grid_merge_mean_threshold,
                  (unsigned long long)_grid_merge_glob_threshold,
                  _grid_merge_min_ranks);
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
