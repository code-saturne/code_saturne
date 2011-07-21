/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Grid connectivity and data used for multigrid coarsening
 * and associated matrix construction.
 *============================================================================*/

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
#include <fvm_order.h>
#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_matrix.h"
#include "cs_prototypes.h"
#include "cs_perio.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_grid.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

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
  cs_bool_t           symmetric;    /* Symmetric matrix coefficients
                                       indicator */

  fvm_lnum_t          n_cells;      /* Local number of cells */
  fvm_lnum_t          n_cells_ext;  /* Local number of participating cells
                                       (cells + ghost cells sharing a face) */
  fvm_lnum_t          n_faces;      /* Local number of faces */
  fvm_gnum_t          n_g_cells;    /* Global number of cells */

  fvm_lnum_t          n_cells_r[2]; /* Size of array used for restriction
                                       operations ({ncells, n_cells_ext} when
                                       no grid merging has taken place) */

  /* Grid hierarchy information */

  const struct _cs_grid_t  *parent; /* Pointer to parent (finer) grid */

  /* Connectivity information */

  const fvm_lnum_t   *face_cell;    /* Face -> cells connectivity (1 to n) */
  fvm_lnum_t         *_face_cell;   /* Face -> cells connectivity
                                       (private array) */

  /* Restriction from parent to current level */

  fvm_lnum_t         *coarse_cell;  /* Fine -> coarse cell connectivity
                                       (1 to n); size: parent n_cells_ext */
  fvm_lnum_t         *coarse_face;  /* Fine -> coarse face connectivity
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

  const cs_real_t  *xa;             /* Extra-diagonal (shared) */
  cs_real_t        *_xa;            /* Extra-diagonal (private) */

  const cs_real_t  *xa0;            /* Symmetrized extra-diagonal (shared) */
  cs_real_t        *_xa0;           /* Symmetrized extra-diagonal (private) */

  cs_real_t        *xa0ij;

  cs_matrix_structure_t   *matrix_struct;  /* Associated matrix structure */
  cs_matrix_t             *matrix;         /* Associated matrix */

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

  fvm_lnum_t       *merge_cell_idx;    /* start cell_id for each sub-rank
                                          when merge_sub_rank = 0
                                          (size: merge_size + 1) */

  int               n_ranks;           /* Number of active ranks */
  int               comm_id;           /* Associated communicator
                                          (owner when merge_idx != NULL) */

#endif
};

/*============================================================================
 *  Global variables
 *============================================================================*/

#if defined(HAVE_MPI)

static fvm_gnum_t  _grid_merge_threshold[2] = {300, 500};
static int         _grid_merge_min_ranks = 1;
static int         _grid_merge_stride = 4;

static int        _n_grid_comms = 0;
static int       *_grid_ranks = NULL;
static MPI_Comm  *_grid_comm = NULL;

#endif /* defined(HAVE_MPI) */

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

  g->level = 0;

  g->n_cells = 0;
  g->n_cells_ext = 0;
  g->n_faces = 0;
  g->n_g_cells = 0;

  g->n_cells_r[0] = 0;
  g->n_cells_r[1] = 0;

  g->parent = NULL;

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
  g->xa = NULL;
  g->_xa = NULL;
  g->xa0 = NULL;
  g->_xa0 = NULL;

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
  fvm_lnum_t ii;
  cs_grid_t *c = NULL;

  c = _create_grid();

  c->parent = f;

  c->level = f->level + 1;
  c->symmetric = f->symmetric;

  BFT_MALLOC(c->coarse_cell, f->n_cells_ext, fvm_lnum_t);

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
 *   coarse_face_cell --> coarse face -> cell connectivity (1 to n numbering)
 *                        size: n_coarse_facees * 2
 *----------------------------------------------------------------------------*/

static void
_coarsen_faces(const cs_grid_t    *fine,
               const fvm_lnum_t   *restrict coarse_cell,
               fvm_lnum_t          n_coarse_cells,
               fvm_lnum_t         *n_coarse_faces,
               fvm_lnum_t        **coarse_face,
               fvm_lnum_t        **coarse_face_cell)
{
  fvm_lnum_t  ii, jj, face_id, connect_size;

  fvm_lnum_t  *restrict c_cell_cell_cnt = NULL;
  fvm_lnum_t  *restrict c_cell_cell_idx = NULL;
  fvm_lnum_t  *restrict c_cell_cell_id = NULL;
  fvm_lnum_t  *restrict c_cell_cell_face = NULL;

  fvm_lnum_t  *restrict _coarse_face = NULL;
  fvm_lnum_t  *restrict _c_face_cell = NULL;

  fvm_lnum_t   c_n_faces = 0;

  const fvm_lnum_t c_n_cells = n_coarse_cells;
  const fvm_lnum_t f_n_faces = fine->n_faces;
  const fvm_lnum_t *restrict f_face_cell = fine->face_cell;

  /* Pre-allocate return values
     (coarse face->cell connectivity is over-allocated) */

  BFT_MALLOC(_coarse_face, f_n_faces, fvm_lnum_t);
  BFT_MALLOC(_c_face_cell, f_n_faces*2, fvm_lnum_t);

  for (face_id = 0; face_id < f_n_faces; _coarse_face[face_id++] = 0);

  /* Allocate index */

  BFT_MALLOC(c_cell_cell_idx, c_n_cells + 1, fvm_lnum_t);

  for (ii = 0; ii <= c_n_cells; c_cell_cell_idx[ii++] = 0);

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    ii = coarse_cell[f_face_cell[face_id*2]     - 1] - 1;
    jj = coarse_cell[f_face_cell[face_id*2 + 1] - 1] - 1;

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
             fvm_lnum_t);

  BFT_MALLOC(c_cell_cell_face,
             c_cell_cell_idx[c_n_cells],
             fvm_lnum_t);

  for (ii = 0, connect_size = c_cell_cell_idx[c_n_cells];
       ii < connect_size;
       c_cell_cell_face[ii++] = 0);

  /* Use a counter for array population, as array will usually
     not be fully populated */

  BFT_MALLOC(c_cell_cell_cnt, c_n_cells, fvm_lnum_t);

  for (ii = 0; ii < c_n_cells; c_cell_cell_cnt[ii++] = 0);

  /* Build connectivity */

  c_n_faces = 0;

  for (face_id = 0; face_id < f_n_faces; face_id++) {

    fvm_lnum_t kk, start_id, end_id;
    fvm_lnum_t sign = 1;

    ii = coarse_cell[f_face_cell[face_id*2]     - 1] - 1;
    jj = coarse_cell[f_face_cell[face_id*2 + 1] - 1] - 1;

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
        _c_face_cell[c_n_faces*2]     = ii + 1;
        _c_face_cell[c_n_faces*2 + 1] = jj + 1;
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

  BFT_REALLOC(_c_face_cell, c_n_faces*2, fvm_lnum_t);

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
                          cs_int_t          coarse_send[],
                          cs_int_t          coarse_cell[])
{
  fvm_lnum_t  i, start, length;

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
      length = halo->index[2*rank_id + 2] - halo->index[2*rank_id] ;

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

      cs_int_t *_coarse_cell
        = coarse_cell + halo->n_local_elts + halo->index[2*local_rank_id];

      start = halo->send_index[2*local_rank_id];
      length =   halo->send_index[2*local_rank_id + 2]
               - halo->send_index[2*local_rank_id];

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
  fvm_lnum_t ii, jj;
  fvm_lnum_t start_id, end_id, sub_count;

  fvm_lnum_t *start_end_id = NULL;
  fvm_lnum_t *sub_num = NULL;
  cs_int_t   *coarse_send = NULL;

  fvm_lnum_t *restrict coarse_cell = c->coarse_cell;

  cs_halo_t *c_halo = NULL;
  const cs_halo_t *f_halo = f->halo;

  const fvm_lnum_t c_n_cells = c->n_cells;

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

  for (ii = f->n_cells; ii < f->n_cells_ext; coarse_cell[ii++] = 0);

  /* Allocate and initialize counters */

  BFT_MALLOC(start_end_id, n_sections*2, fvm_lnum_t);
  BFT_MALLOC(sub_num, c_n_cells, fvm_lnum_t);
  BFT_MALLOC(coarse_send, n_f_send, cs_int_t);

  for (ii = 0; ii < c_n_cells; ii++)
    sub_num[ii] = -1;

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

  BFT_MALLOC(c_halo->send_list, c_halo->n_send_elts[0], cs_int_t);

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
  fvm_lnum_t  ii, jj, face_id;

  fvm_lnum_t  c_n_cells = 0;

  const fvm_lnum_t f_n_faces = f->n_faces;
  const fvm_lnum_t *restrict f_face_cell = f->face_cell;

  /* Sanity check */

  for (face_id = 0; face_id < f_n_faces; face_id++) {
    ii = f_face_cell[face_id*2]     - 1;
    jj = f_face_cell[face_id*2 + 1] - 1;
    if (ii == jj)
      bft_error(__FILE__, __LINE__, 0,
                _("Connectivity error:\n"
                  "Face %d has same cell %d on both sides."),
                (int)(face_id+1), (int)(ii+1));
  }

  /* Compute number of coarse cells */

  for (ii = 0; ii < f->n_cells; ii++) {
    if (c->coarse_cell[ii] > c_n_cells)
      c_n_cells = c->coarse_cell[ii];
  }

  c->n_cells = c_n_cells;
  c->n_g_cells = c_n_cells;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    fvm_gnum_t _c_n_cells = c_n_cells;
    MPI_Allreduce(&_c_n_cells, &(c->n_g_cells), 1, FVM_MPI_GNUM, MPI_SUM,
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

  c->face_cell = c->_face_cell;
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
  BFT_MALLOC(_grid_ranks, _n_grid_comms, cs_int_t);

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
_rebuild_halo_send_lists(cs_halo_t   *h,
                         fvm_lnum_t   new_src_cell_id[])
{
  /* As halos are based on interior faces, every domain is both
     sender and receiver */

  fvm_lnum_t start, length;

  int rank_id, tr_id;
  int n_sections = 1 + h->n_transforms;
  int request_count = 0;
  fvm_lnum_t *send_buf = NULL, *recv_buf = NULL;
  MPI_Status *status = NULL;
  MPI_Request *request = NULL;

  BFT_MALLOC(status, h->n_c_domains*2, MPI_Status);
  BFT_MALLOC(request, h->n_c_domains*2, MPI_Request);
  BFT_MALLOC(send_buf, h->n_c_domains*n_sections, fvm_lnum_t);
  BFT_MALLOC(recv_buf, h->n_c_domains*n_sections, fvm_lnum_t);

  /* Exchange sizes */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++)
    MPI_Irecv(recv_buf + rank_id*n_sections,
              n_sections,
              FVM_MPI_LNUM,
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
              FVM_MPI_LNUM,
              h->c_domain_rank[rank_id],
              cs_glob_rank_id,
              cs_glob_mpi_comm,
              &(request[request_count++]));
  }

  /* Wait for all exchanges */

  MPI_Waitall(request_count, request, status);
  request_count = 0;

  /* Update sizes */

  BFT_MALLOC(h->send_index, h->n_c_domains*2 + 1, fvm_lnum_t);
  h->send_index[0] = 0;
  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    h->send_index[rank_id*2 + 1]
      = h->send_index[rank_id*2] + recv_buf[rank_id*n_sections];
    h->send_index[rank_id*2 + 2] = h->send_index[rank_id*2 + 1];
  }

  /* Update send_perio_lst in case of transforms */

  if (h->n_transforms > 0) {
    BFT_MALLOC(h->send_perio_lst, h->n_c_domains*h->n_transforms*4, fvm_lnum_t);
    for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
      fvm_lnum_t n_cur_vals = recv_buf[rank_id*n_sections];
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++)
        n_cur_vals -= recv_buf[rank_id*n_sections + 1 + tr_id];
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
        fvm_lnum_t n_tr_vals = recv_buf[rank_id*n_sections + 1 + tr_id];
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

  BFT_MALLOC(h->send_list, h->n_send_elts[0], fvm_lnum_t);

  /* Receive data from distant ranks */

  for (rank_id = 0; rank_id < h->n_c_domains; rank_id++) {
    start = h->send_index[2*rank_id];
    length = h->send_index[2*rank_id + 1] - h->send_index[2*rank_id];
    MPI_Irecv(h->send_list + start,
              length,
              FVM_MPI_LNUM,
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
              FVM_MPI_LNUM,
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
 *   new_halo_cell_num --> new cell number for each halo cell
 *                         (< n_new_cells for cells that have become local,
 *                         >= n_new_cells for cells still in halo)
 *----------------------------------------------------------------------------*/

static void
_merge_halo_data(cs_halo_t   *h,
                 int          loc_rank_id,
                 fvm_lnum_t   n_new_cells,
                 fvm_lnum_t   new_src_cell_id[],
                 fvm_lnum_t   new_halo_cell_num[])
{
  int  rank_id, prev_rank_id, tr_id, prev_section_id;
  fvm_lnum_t  ii, cur_id, section_id, src_id, prev_src_id;

  int   stride = (h->n_transforms > 0) ? 3 : 2;
  int   n_c_domains_ini = h->n_c_domains;

  fvm_lnum_t   n_elts_ini = h->n_elts[0];
  fvm_lnum_t  *order = NULL, *section_idx = NULL;
  fvm_gnum_t  *tmp_num = NULL;

  const int  n_sections = h->n_transforms + 1;

  if (h->n_elts[0] < 1) {
    _empty_halo(h);
    return;
  }

  /* Order list by rank, transform, and new element number */

  BFT_MALLOC(tmp_num, n_elts_ini*stride, fvm_gnum_t);

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
      tmp_num[ii*3 + 2] = 0;
      tmp_num[ii*3 + 2] = new_src_cell_id[ii];
    }

    for (rank_id = 0; rank_id < n_c_domains_ini; rank_id++) {
      for (tr_id = 0; tr_id < h->n_transforms; tr_id++) {
        fvm_lnum_t ii_0
          = h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id];
        fvm_lnum_t ii_1
          = ii_0 + h->perio_lst[h->n_c_domains*4*tr_id + 4*rank_id + 1];
        for (ii = ii_0; ii < ii_1; ii++)
          tmp_num[ii*3 + 1] = tr_id + 1;
      }
    }
  }

  order = fvm_order_local_s(NULL, tmp_num, stride, n_elts_ini);

  /* Rebuilt lists and build renumbering */

  cur_id = order[0];
  h->n_c_domains = 0;
  h->index[0] = 0;
  h->n_elts[0] = 0;
  h->n_elts[1] = 0;

  prev_rank_id = -1;
  prev_src_id = 0;

  if (stride == 2) {

    for (ii = 0; ii < n_elts_ini; ii++) {

      cs_bool_t is_same = true;
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

        new_halo_cell_num[cur_id] = n_new_cells + h->n_elts[0];

        prev_rank_id = rank_id;
        prev_src_id = src_id;

      }
      else /* if (rank_id == loc_rank_id) */
        new_halo_cell_num[cur_id] = tmp_num[cur_id*2 + 1] + 1;
    }

  }
  else { /* if (stride == 3) */

    const fvm_lnum_t section_idx_size = n_sections * h->n_c_domains + 1;

    prev_section_id = -1;

    /* Initialize index as count */

    BFT_MALLOC(section_idx, section_idx_size, fvm_lnum_t);
    for (ii = 0; ii < section_idx_size; ii++)
      section_idx[ii] = 0;

    for (ii = 0; ii < n_elts_ini; ii++) {

      cs_bool_t is_same = true;
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
        }
        if (section_id != prev_section_id)
          is_same = false;
        if (src_id != prev_src_id)
          is_same = false;

        if (is_same == false) {
          new_src_cell_id[h->n_elts[0]] = src_id;
          h->n_elts[0] += 1;
          section_idx[h->n_c_domains*n_sections + section_id + 1] += 1;
        }

        new_halo_cell_num[cur_id] = n_new_cells + h->n_elts[0];

        prev_rank_id = rank_id;
        prev_section_id = section_id;
        prev_src_id = src_id;

      }
      else /* if (rank_id == loc_rank_id && tmp_num[cur_id*3 + 1] == 0) */
        new_halo_cell_num[cur_id] = tmp_num[cur_id*2 + 1];

    }

    /* Transform count to index */

    for (ii = 1; ii < section_idx_size; ii++)
      section_idx[ii] += section_idx[ii - 1];
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
  BFT_REALLOC(h->index, h->n_c_domains*2+1, fvm_lnum_t);
  if (h->n_transforms > 0)
    BFT_REALLOC(h->perio_lst,
                h->n_c_domains * h->n_transforms * 4,
                fvm_lnum_t);

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
 *   new_cel_num <-> new cell numbering for local cells
 *                   in: defined for local cells
 *                   out: updated also for halo cells
 *----------------------------------------------------------------------------*/

static void
_append_halos(cs_grid_t   *g,
              fvm_lnum_t  *new_cell_num)
{
  fvm_lnum_t ii, jj;
  int rank_id;
  int counts[3];

  int *recv_count = NULL;
  fvm_lnum_t *new_src_cell_id = NULL, *new_halo_cell_num = NULL;

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

  /* Reallocate arrays for receiving rank and append data */

  if (g->merge_sub_rank == 0) {

    BFT_MALLOC(new_src_cell_id, counts[2], fvm_lnum_t);
    for (ii = g->n_cells, jj = 0; ii < g->n_cells_ext; ii++, jj++)
      new_src_cell_id[jj] = new_cell_num[ii] - 1;

    BFT_REALLOC(h->c_domain_rank, counts[0], int);
    BFT_REALLOC(h->index, counts[0]*2 + 1, fvm_lnum_t);
    BFT_REALLOC(h->perio_lst, counts[0]*n_transforms*4, fvm_lnum_t);

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      int n_c_domains_r = recv_count[rank_id*3 + 0];
      fvm_lnum_t n_recv = recv_count[rank_id*3 + 2];

      fvm_lnum_t index_shift = h->index[2*h->n_c_domains];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(h->c_domain_rank + h->n_c_domains, n_c_domains_r,
               MPI_INT, dist_rank, tag, comm, &status);
      MPI_Recv(new_src_cell_id + h->n_elts[0], n_recv,
               FVM_MPI_LNUM, dist_rank, tag, comm, &status);
      MPI_Recv(h->index + h->n_c_domains*2, n_c_domains_r*2+1,
               FVM_MPI_LNUM, dist_rank, tag, comm, &status);

      for (ii = 0, jj = h->n_c_domains*2;
           ii < n_c_domains_r*2+1;
           ii++, jj++)
        h->index[jj] += index_shift;

      if (n_transforms > 0)
        MPI_Recv(h->perio_lst + h->n_c_domains*n_transforms*4,
                 n_c_domains_r*n_transforms*4,
                 FVM_MPI_LNUM, dist_rank, tag, comm, &status);

      /* Update halo sizes */

      h->n_local_elts += recv_count[rank_id*3 + 1];
      h->n_c_domains += n_c_domains_r;
      h->n_elts[0] += n_recv;
      h->n_elts[1] = h->n_elts[0];
    }

  }
  else if (g->merge_sub_size > 1) {

    BFT_MALLOC(new_src_cell_id, h->n_elts[0], fvm_lnum_t);
    for (ii = g->n_cells, jj = 0; ii < g->n_cells_ext; ii++, jj++)
      new_src_cell_id[jj] = new_cell_num[ii] - 1;

    MPI_Send(h->c_domain_rank, h->n_c_domains, MPI_INT,
             g->merge_sub_root, tag, comm);
    MPI_Send(new_src_cell_id, h->n_elts[0], FVM_MPI_LNUM,
             g->merge_sub_root, tag, comm);
    MPI_Send(h->index, h->n_c_domains*2+1, FVM_MPI_LNUM,
             g->merge_sub_root, tag, comm);

    if (n_transforms > 0)
      MPI_Send(h->perio_lst, h->n_c_domains*n_transforms*4,
               FVM_MPI_LNUM, g->merge_sub_root, tag, comm);

    _empty_halo(h);
  }

  /* Cleanup halo and set pointer for coarsening (sub_rooot) ranks*/

  if (h != NULL) {

    BFT_MALLOC(new_halo_cell_num, h->n_elts[0], fvm_lnum_t);

    _merge_halo_data(h,
                     cs_glob_rank_id,
                     counts[1],
                     new_src_cell_id,
                     new_halo_cell_num);

    _rebuild_halo_send_lists(h, new_src_cell_id);

  }

  if (new_src_cell_id != NULL)
    BFT_FREE(new_src_cell_id);

  cs_halo_update_buffers(h);

  g->halo = h;
  g->_halo = h;

  /* Finally, update halo section of cell renumbering array */

  if (g->merge_sub_rank == 0) {

    fvm_lnum_t n_send = recv_count[2];
    fvm_lnum_t send_shift = n_send;

    for (ii = 0; ii < n_send; ii++)
      new_cell_num[g->n_cells + ii] = new_halo_cell_num[ii];

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      n_send = recv_count[rank_id*3 + 2];
      MPI_Send(new_halo_cell_num + send_shift, n_send, FVM_MPI_LNUM,
               dist_rank, tag, comm);
      send_shift += n_send;
    }

    BFT_FREE(recv_count);
  }
  else if (g->merge_sub_size > 1)
    MPI_Recv(new_cell_num + g->n_cells, g->n_cells_ext - g->n_cells,
             FVM_MPI_LNUM, g->merge_sub_root, tag, comm, &status);

  BFT_FREE(new_halo_cell_num);
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
    BFT_REALLOC(g->_da, g->n_cells_ext, cs_real_t);

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      fvm_lnum_t n_recv = (  g->merge_cell_idx[rank_id+1]
                           - g->merge_cell_idx[rank_id]);
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_cell_cen + g->merge_cell_idx[rank_id]*3, n_recv*3,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->_cell_vol + g->merge_cell_idx[rank_id], n_recv,
               CS_MPI_REAL, dist_rank, tag, comm, &status);

      MPI_Recv(g->_da + g->merge_cell_idx[rank_id], n_recv,
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

    MPI_Send(g->_da, g->n_cells, CS_MPI_REAL, g->merge_sub_root, tag, comm);
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
      cs_perio_sync_coords(g->halo, CS_HALO_STANDARD, g->_cell_cen);

    cs_halo_sync_var(g->halo, CS_HALO_STANDARD, g->_cell_vol);

    cs_halo_sync_var(g->halo, CS_HALO_STANDARD, g->_da);

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
                  fvm_lnum_t   n_faces,
                  fvm_lnum_t  *face_list)
{
  int rank_id;

  fvm_lnum_t *recv_count = NULL;

  MPI_Status status;
  MPI_Comm  comm = cs_glob_mpi_comm;

  static const int tag = 'a'+'p'+'p'+'e'+'n'+'d'+'_'+'f';

  /* Exchange counters needed for concatenation */

  if (g->merge_sub_rank == 0) {
    BFT_MALLOC(recv_count, g->merge_sub_size, fvm_lnum_t);
    recv_count[0] = g->n_faces;
    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;
      MPI_Recv(recv_count + rank_id, 1, FVM_MPI_LNUM, dist_rank,
               tag, comm, &status);
    }
  }
  else
    MPI_Send(&n_faces, 1, FVM_MPI_LNUM, g->merge_sub_root, tag, comm);

  /* Reallocate arrays for receiving rank and append data */

  BFT_FREE(g->coarse_face);

  if (g->merge_sub_rank == 0) {

    fvm_lnum_t n_faces_tot = 0;

    for (rank_id = 0; rank_id < g->merge_sub_size; rank_id++)
      n_faces_tot += recv_count[rank_id];

    BFT_REALLOC(g->_face_cell, n_faces_tot*2, fvm_lnum_t);

    BFT_REALLOC(g->_face_normal, n_faces_tot*3, cs_real_t);

    if (g->symmetric == true)
      BFT_REALLOC(g->_xa, n_faces_tot, cs_real_t);
    else
      BFT_REALLOC(g->_xa, n_faces_tot*2, cs_real_t);

    BFT_REALLOC(g->_xa0, n_faces_tot, cs_real_t);
    BFT_REALLOC(g->xa0ij, n_faces_tot*3, cs_real_t);

    g->n_faces = recv_count[0];

    for (rank_id = 1; rank_id < g->merge_sub_size; rank_id++) {

      fvm_lnum_t n_recv = recv_count[rank_id];
      int dist_rank = g->merge_sub_root + g->merge_stride*rank_id;

      MPI_Recv(g->_face_cell + g->n_faces*2, n_recv*2,
               FVM_MPI_LNUM, dist_rank, tag, comm, &status);

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

    fvm_lnum_t face_id = 0;

    /* Filter face connectivity then send it */
    for (face_id = 0; face_id < n_faces; face_id++) {
      fvm_lnum_t p_face_id = face_list[face_id];
      g->_face_cell[face_id*2] = g->_face_cell[p_face_id*2];
      g->_face_cell[face_id*2 + 1] = g->_face_cell[p_face_id*2 + 1];
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

    MPI_Send(g->_face_cell, n_faces*2, FVM_MPI_LNUM,
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

  g->face_cell = g->_face_cell;
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
  fvm_lnum_t j, face_id;
  int base_rank = cs_glob_rank_id;
  fvm_lnum_t cell_shift = 0;
  fvm_lnum_t n_faces = 0;
  fvm_lnum_t *new_cell_num = NULL, *face_list = NULL;
  cs_bool_t  *halo_cell_flag = NULL;
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
      BFT_MALLOC(g->merge_cell_idx, g->merge_sub_size + 1, fvm_lnum_t);
      g->merge_cell_idx[0] = 0; g->merge_cell_idx[1] = g->n_cells;
      for (i = 1; i < g->merge_sub_size; i++) {
        fvm_lnum_t recv_val;
        int dist_rank = g->merge_sub_root + g->merge_stride*i;
        MPI_Recv(&recv_val, 1, FVM_MPI_LNUM, dist_rank, tag, comm, &status);
        g->merge_cell_idx[i + 1] = recv_val + g->merge_cell_idx[i];
      }
    }
    else {
      fvm_lnum_t send_val = g->n_cells;
      MPI_Send(&send_val, 1, FVM_MPI_LNUM, g->merge_sub_root, tag, comm);
    }

    /* Now return computed info to grids that will be merged */

    if (g->merge_sub_rank == 0) {
      for (i = 1; i < g->merge_sub_size; i++) {
        fvm_lnum_t send_val;
        send_val = g->merge_cell_idx[i];
        MPI_Send(&send_val, 1, FVM_MPI_LNUM,
                 g->merge_sub_root + g->merge_stride*i, tag, comm);
      }
    }
    else
      MPI_Recv(&cell_shift, 1, FVM_MPI_LNUM,
               g->merge_sub_root, tag, comm, &status);
  }

  /* Compute and exchange new cell numbers */

  BFT_MALLOC(new_cell_num, g->n_cells_ext, fvm_lnum_t);
  for (j = 0; j < g->n_cells; j++)
    new_cell_num[j] = cell_shift + j + 1;
  for (j = g->n_cells; j < g->n_cells; j++)
    new_cell_num[j] = 0;

  cs_halo_sync_untyped(g->halo,
                       CS_HALO_STANDARD,
                       sizeof(fvm_lnum_t),
                       new_cell_num);

  /* Now build face filter list (before halo is modified) */

  if (g->merge_sub_size > 1 && g->merge_sub_rank > 0 && g->n_faces > 0) {

    fvm_lnum_t n_ghost_cells = g->n_cells_ext - g->n_cells;

    /* Mark which faces should be merged: to avoid duplicates, a face
       connected to a cell on a lower rank in the same merge set is
       discarded, as it has already been accounted for by that rank. */

    BFT_MALLOC(face_list, g->n_faces, fvm_lnum_t);
    BFT_MALLOC(halo_cell_flag, n_ghost_cells, cs_bool_t);
    for (j = 0; j < n_ghost_cells; j++)
      halo_cell_flag[j] = false;

    for (i = 0; i < g->halo->n_c_domains; i++) {
      rank_id = g->halo->c_domain_rank[i];
      if (rank_id >= g->merge_sub_root && rank_id < cs_glob_rank_id) {
        for (t_id = 0; t_id < g->halo->n_transforms; t_id++) {
          int t_shift = 4 * g->halo->n_c_domains * t_id;
          fvm_lnum_t t_start = g->halo->perio_lst[t_shift + 4*i];
          fvm_lnum_t t_end = t_start + g->halo->perio_lst[t_shift + 4*i + 1];
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
      cs_bool_t use_face = true;
      fvm_lnum_t ii = g->face_cell[face_id*2] - g->n_cells;
      fvm_lnum_t jj = g->face_cell[face_id*2 + 1] - g->n_cells;
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

  _append_halos(g, new_cell_num);

  /* Update face ->cells connectivity */

  for (face_id = 0; face_id < g->n_faces; face_id++) {
    fvm_lnum_t ii = g->face_cell[face_id*2] - 1;
    fvm_lnum_t jj = g->face_cell[face_id*2 + 1] - 1;
    assert(ii != jj && new_cell_num[ii] != new_cell_num[jj]);
    g->_face_cell[face_id*2] = new_cell_num[ii];
    g->_face_cell[face_id*2 + 1] = new_cell_num[jj];
  }

  BFT_FREE(new_cell_num);

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
        fvm_lnum_t n_send = (  g->merge_cell_idx[rank_id+1]
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

/*============================================================================
 *  Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Set the default parameters for multigrid coarsening.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmopt, CLMOPT)
(
 const cs_int_t   *mltmmn,    /* <-- Mean number of cells under which merging
                               *     should take place */
 const cs_int_t   *mltmgl,    /* <-- Global number of cells under which
                               *     merging should take place */
 const cs_int_t   *mltmmr,    /* <-- Number of active ranks under which no
                               *     merging takes place */
 const cs_int_t   *mltmst     /* <-- Number of ranks over which merging
                               *     takes place */
)
{
  cs_grid_set_defaults(*mltmmn, *mltmgl, *mltmmr, *mltmst);
}

/*----------------------------------------------------------------------------
 * Print the default parameters for multigrid coarsening.
 *----------------------------------------------------------------------------*/

void CS_PROCF(clmimp, CLMIMP)
(
 void
)
{
  cs_grid_log_defaults();
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create base grid by mapping from shared mesh values.
 *
 * Note that as arrays given as arguments are shared by the created grid
 * (which can only access them, not modify them), the grid should be
 * destroyed before those arrays.
 *
 * parameters:
 *   n_cells     <-- Local number of cells
 *   n_cells_ext <-- Local number of cells + ghost cells
 *   n_faces     <-- Local number of faces
 *   symmetric   <-- True if xam is symmetric, false otherwise
 *   face_cell   <-- Face -> cells connectivity (1 to n)
 *   halo        <-- Halo structure associated with this level, or NULL.
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   cell_cen    <-- Cell center (size: 3.n_cells_ext)
 *   cell_vol    <-- Cell volume (size: n_cells_ext)
 *   face_normal <-- Internal face normals (size: 3.n_faces)
 *   da          <-- Matrix diagonal (size: n_cell_ext)
 *   xa          <-- Matrix extra-diagonal terms
 *                   (size: n_faces if symmetric, 2.n_faces otherwise)
 *
 * returns:
 *   base grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_create_from_shared(fvm_lnum_t             n_cells,
                           fvm_lnum_t             n_cells_ext,
                           fvm_lnum_t             n_faces,
                           cs_bool_t              symmetric,
                           const fvm_lnum_t      *face_cell,
                           const cs_halo_t       *halo,
                           const cs_numbering_t  *numbering,
                           const cs_real_t       *cell_cen,
                           const cs_real_t       *cell_vol,
                           const cs_real_t       *face_normal,
                           const cs_real_t       *da,
                           const cs_real_t       *xa)
{
  fvm_lnum_t ii, jj, kk, face_id;

  cs_grid_t *g = NULL;

  /* Create empty structure and map base data */

  g = _create_grid();

  g->level = 0;
  g->symmetric = symmetric;

  g->n_cells = n_cells;
  g->n_cells_ext = n_cells_ext;
  g->n_faces = n_faces;
  g->n_g_cells = n_cells;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    fvm_gnum_t _n_cells = n_cells;
    MPI_Allreduce(&_n_cells, &(g->n_g_cells), 1, FVM_MPI_GNUM, MPI_SUM,
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

  g->xa = xa;
  g->_xa = NULL;

  /* Build symmetrized extra-diagonal terms if necessary,
     or point to existing terms if already symmetric */

  if (symmetric == true) {
    g->xa0 = g->xa;
    g->_xa0 = NULL;
  }
  else {
    BFT_MALLOC(g->_xa0, n_faces, cs_real_t);
    for (face_id = 0; face_id < n_faces; face_id++)
      g->_xa0[face_id] = 0.5 * (xa[face_id*2] + xa[face_id*2+1]);
    g->xa0 = g->_xa0;
  }

  /* Compute multigrid-specific terms */

  BFT_MALLOC(g->xa0ij, n_faces*3, cs_real_t);

  for (face_id = 0; face_id < n_faces; face_id++) {
    ii = face_cell[face_id*2] - 1;
    jj = face_cell[face_id*2 + 1] - 1;
    for (kk = 0; kk < 3; kk++)
      g->xa0ij[face_id*3 + kk] = g->xa0[face_id] * (  cell_cen[jj*3 + kk]
                                                    - cell_cen[ii*3 + kk]);
  }

  g->matrix_struct = cs_matrix_structure_create(CS_MATRIX_NATIVE,
                                                true,
                                                n_cells,
                                                n_cells_ext,
                                                n_faces,
                                                NULL,
                                                face_cell,
                                                halo,
                                                numbering);


  g->matrix = cs_matrix_create(g->matrix_struct);

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
      g->_halo = cs_halo_destroy(g->_halo);

    if (g->_da != NULL)
      BFT_FREE(g->_da);
    if (g->_xa != NULL)
      BFT_FREE(g->_xa);
    if (g->_xa0 != NULL)
      BFT_FREE(g->_xa0);

    BFT_FREE(g->xa0ij);

    cs_matrix_destroy(&(g->matrix));
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
 *   n_ranks     --> number of ranks with data (or NULL)
 *   n_cells     --> Number of local cells (or NULL)
 *   n_cells_ext --> Number of cells including ghosts (or NULL)
 *   n_faces     --> Number of faces (or NULL)
 *   n_g_cells   --> Number of global cells (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 cs_bool_t        *symmetric,
                 int              *n_ranks,
                 fvm_lnum_t       *n_cells,
                 fvm_lnum_t       *n_cells_ext,
                 fvm_lnum_t       *n_faces,
                 fvm_gnum_t       *n_g_cells)
{
  assert(g != NULL);

  if (level != NULL)
    *level = g->level;

  if (symmetric != NULL)
    *symmetric = g->symmetric;

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

fvm_lnum_t
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

fvm_lnum_t
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

fvm_lnum_t
cs_grid_get_n_cells_max(const cs_grid_t  *g)
{
  fvm_lnum_t retval = 0;

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

fvm_gnum_t
cs_grid_get_n_g_cells(const cs_grid_t  *g)
{
  assert(g != NULL);

  return g->n_g_cells;
}

/*----------------------------------------------------------------------------
 * Get grid's associated matrix information.
 *
 * parameters:
 *   g           <-- Grid structure
 *   da          --> Diagonal matrix coefficients
 *   xa          --> Non-diagonal matrix coefficients
 *   m           --> Associated matrix structure
 *----------------------------------------------------------------------------*/

void
cs_grid_get_matrix(const cs_grid_t   *g,
                   const cs_real_t  **da,
                   const cs_real_t  **xa,
                   cs_matrix_t      **m)
{
  assert(g != NULL);

  if (da != NULL)
    *da = g->da;

  if (xa != NULL)
    *xa = g->xa;

  if (m != NULL)
    *m = g->matrix;
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
 *   f                   <-- Fine grid structure
 *   verbosity           <-- Verbosity level
 *   agglomeration_limit <-- Maximum allowed fine cells per coarse cell
 *   max_agglomeration   <-> Maximum fine cells per coarse cell
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

cs_grid_t *
cs_grid_coarsen(const cs_grid_t   *f,
                int                verbosity,
                int                agglomeration_limit,
                int               *max_agglomeration,
                double             relaxation_parameter)
{
  cs_int_t iappel, iusmgr;
  cs_int_t igr;
  cs_int_t isym = 2;
  cs_int_t niw = 0, nrw = 0;

  cs_int_t iwarnp = verbosity;
  cs_int_t f_n_cells = f->n_cells;
  cs_int_t f_n_cells_ext = f->n_cells_ext;
  cs_int_t f_n_faces = f->n_faces;
  cs_int_t c_n_cells = 0;
  cs_int_t c_n_cells_ext = 0;
  cs_int_t c_n_faces = 0;

  cs_real_t *rwc1 = NULL, *rwc2 = NULL, *rwc3 = NULL, *rwc4 = NULL;

  cs_grid_t *c = NULL;

  assert(f != NULL);

  /* Initialization */

  c = _coarse_init(f);

  if (f->symmetric == true)
    isym = 1;

  igr = c->level;

  /* Determine if user or automatic method is chosen */

  iusmgr = 0;

  CS_PROCF (ustmgr, USTMGR) (&iappel, &igr, &isym,
                             &f_n_cells, &f_n_cells_ext, &f_n_faces, &iwarnp,
                             &iusmgr, &niw, &nrw,
                             f->face_cell,
                             f->da, f->xa,
                             f->face_normal, f->cell_vol, f->cell_cen,
                             NULL, NULL, NULL);

  /* Automatic coarsening */

  if (iusmgr == 0) {

    cs_int_t nagmax = agglomeration_limit;
    cs_int_t iagmax = *max_agglomeration;
    cs_int_t *indic = NULL, *inombr = NULL;
    cs_int_t *irsfac = NULL, *indicf = NULL;
    cs_real_t *w1 = NULL, *w2 = NULL;

    /* Allocate working arrays */

    BFT_MALLOC(indic, f_n_cells_ext*3 + f_n_faces*2, cs_int_t);
    inombr = indic + f_n_cells_ext;
    irsfac = inombr + f_n_cells_ext;
    indicf = irsfac + f_n_faces;

    BFT_MALLOC(w1, f_n_cells_ext*2, cs_real_t);
    w2 = w1 + f_n_cells_ext;

    /* Determine fine->coarse cell connectivity (agglomeration) */

    CS_PROCF(autmgr, AUTMGR) (&igr, &isym, &iagmax, &nagmax,
                              &f_n_cells, &f_n_cells_ext, &f_n_faces, &iwarnp,
                              f->face_cell,
                              f->da, f->xa,
                              f->face_normal, f->cell_vol, f->cell_cen,
                              c->coarse_cell,
                              indic, inombr, irsfac, indicf, w1, w2);

    /* Free working arrays */

    BFT_FREE(indic);
    inombr = NULL;
    irsfac = NULL;
    indicf = NULL;

    BFT_FREE(w1);
    w2 = NULL;

    *max_agglomeration = iagmax;
  }

  /* User coarsening */

  else { /* if iusmgr == 1 */

    cs_int_t *iw = NULL;
    cs_real_t *rw = NULL;

    BFT_MALLOC(iw, niw, cs_int_t);
    BFT_MALLOC(rw, nrw, cs_real_t);

    iappel = 2;

    CS_PROCF (ustmgr, USTMGR) (&iappel, &igr, &isym,
                               &f_n_cells, &f_n_cells_ext, &f_n_faces, &iwarnp,
                               &iusmgr, &niw, &nrw,
                               f->face_cell,
                               f->da, f->xa,
                               f->face_normal, f->cell_vol, f->cell_cen,
                               c->coarse_cell,
                               iw, rw);

    BFT_FREE(iw);
    BFT_FREE(rw);

  }

  /* Build coarse grid connectivity */

  _coarsen(f, c);

  /* Allocate permanent arrays in coarse grid */

  BFT_MALLOC(c->_cell_cen, c->n_cells_ext*3, cs_real_t);
  c->cell_cen = c->_cell_cen;

  BFT_MALLOC(c->_cell_vol, c->n_cells_ext, cs_real_t);
  c->cell_vol = c->_cell_vol;

  BFT_MALLOC(c->_face_normal, c->n_faces*3, cs_real_t);
  c->face_normal = c->_face_normal;

  BFT_MALLOC(c->_da, c->n_cells_ext, cs_real_t);
  c->da = c->_da;

  BFT_MALLOC(c->_xa, c->n_faces*isym, cs_real_t);
  c->xa = c->_xa;

  /* We could have xa0 point to xa if symmetric, but this would require
     caution in CRSTGR to avoid overwriting. */

  BFT_MALLOC(c->_xa0, c->n_faces, cs_real_t);
  c->xa0 = c->_xa0;

  BFT_MALLOC(c->xa0ij, c->n_faces*3, cs_real_t);

  /* Matrix-related data */

  c_n_cells = c->n_cells;
  c_n_cells_ext = c->n_cells_ext;
  c_n_faces = c->n_faces;

  BFT_MALLOC(rwc1, f_n_cells_ext*4, cs_real_t);
  rwc2 = rwc1 + f_n_cells_ext;
  rwc3 = rwc2 + f_n_cells_ext;
  rwc4 = rwc3 + f_n_cells_ext;

  iappel = 1;

  CS_PROCF(crstgr, CRSTGR) (&iappel, &isym, &igr,
                            &f_n_cells, &c_n_cells,
                            &f_n_cells_ext, &c_n_cells_ext,
                            &f_n_faces, &c_n_faces,
                            &iwarnp,
                            f->face_cell, c->face_cell,
                            c->coarse_cell, c->coarse_face,
                            &relaxation_parameter,
                            f->cell_vol, f->cell_cen, f->face_normal,
                            f->xa0, f->xa0ij, f->da, f->xa,
                            c->cell_vol, c->cell_cen, c->face_normal,
                            c->_xa0, c->xa0ij, c->_da, c->_xa,
                            rwc1, rwc2, rwc3, rwc4);

  /* Synchronize grid's geometric quantities */

  if (c->halo != NULL) {

    cs_halo_sync_var_strided(c->halo, CS_HALO_STANDARD, c->_cell_cen, 3);
    if (c->halo->n_transforms > 0)
      cs_perio_sync_coords(c->halo, CS_HALO_STANDARD, c->_cell_cen);

    cs_halo_sync_var(c->halo, CS_HALO_STANDARD, c->_cell_vol);

  }

  iappel = 2;

  CS_PROCF(crstgr, CRSTGR) (&iappel, &isym, &igr,
                            &f_n_cells, &c_n_cells,
                            &f_n_cells_ext, &c_n_cells_ext,
                            &f_n_faces, &c_n_faces,
                            &iwarnp,
                            f->face_cell, c->face_cell,
                            c->coarse_cell, c->coarse_face,
                            &relaxation_parameter,
                            f->cell_vol, f->cell_cen, f->face_normal,
                            f->xa0, f->xa0ij, f->da, f->xa,
                            c->cell_vol, c->cell_cen, c->face_normal,
                            c->_xa0, c->xa0ij, c->_da, c->_xa,
                            rwc1, rwc2, rwc3, rwc4);

  BFT_FREE(rwc1);
  rwc2 = NULL;
  rwc3 = NULL;
  rwc4 = NULL;

  /* Synchronize matrix's geometric quantities */

  if (c->halo != NULL)
    cs_halo_sync_var(c->halo, CS_HALO_STANDARD, c->_da);

  /* Merge grids if we are below the threshold */

#if defined(HAVE_MPI)
  if (c->n_ranks > _grid_merge_min_ranks && _grid_merge_stride > 1) {
    fvm_gnum_t  _n_ranks = c->n_ranks;
    fvm_gnum_t  _n_mean_g_cells = c->n_g_cells / _n_ranks;
    if (   _n_mean_g_cells < _grid_merge_threshold[0]
        || c->n_g_cells < _grid_merge_threshold[1])
    _merge_grids(c, verbosity);
  }
#endif

  c->matrix_struct = cs_matrix_structure_create(CS_MATRIX_NATIVE,
                                                true,
                                                c->n_cells,
                                                c->n_cells_ext,
                                                c->n_faces,
                                                NULL,
                                                c->face_cell,
                                                c->halo,
                                                NULL);

  c->matrix = cs_matrix_create(c->matrix_struct);

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
  fvm_lnum_t ii;

  const fvm_lnum_t *coarse_cell;

  fvm_lnum_t f_n_cells = f->n_cells;
  fvm_lnum_t c_n_cells_ext = c->n_cells_r[1];

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL || f_n_cells == 0);
  assert(f_var != NULL || f_n_cells == 0);
  assert(c_var != NULL || c_n_cells_ext == 0);

  /* Set coarse values */

  coarse_cell = c->coarse_cell;

  for (ii = 0; ii < c_n_cells_ext; ii++)
    c_var[ii] = 0.0;

  for (ii = 0; ii < f_n_cells; ii++)
    c_var[coarse_cell[ii] - 1] += f_var[ii];

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
        fvm_lnum_t n_recv = (  c->merge_cell_idx[rank_id+1]
                             - c->merge_cell_idx[rank_id]);
        int dist_rank = c->merge_sub_root + c->merge_stride*rank_id;
        MPI_Recv(c_var + c->merge_cell_idx[rank_id], n_recv, CS_MPI_REAL,
                 dist_rank, tag, comm, &status);
      }
    }
    else
      MPI_Send(c_var, c->n_cells_r[0], CS_MPI_REAL,
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
  fvm_lnum_t ii;
  const fvm_lnum_t *coarse_cell;
  const int *_c_num = c_num;

  fvm_lnum_t f_n_cells = f->n_cells;

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
  fvm_lnum_t ii;
  const fvm_lnum_t *coarse_cell;
  const cs_real_t *_c_var = c_var;

  fvm_lnum_t f_n_cells = f->n_cells;

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
        fvm_lnum_t n_send = (  c->merge_cell_idx[rank_id+1]
                             - c->merge_cell_idx[rank_id]);
        int dist_rank = c->merge_sub_root + c->merge_stride*rank_id;
        MPI_Send(c_var + c->merge_cell_idx[rank_id], n_send, CS_MPI_REAL,
                 dist_rank, tag, comm);
      }
    }
    else {
      MPI_Status status;
      MPI_Recv(c_var, c->n_cells_r[0], CS_MPI_REAL,
               c->merge_sub_root, tag, comm, &status);
    }
  }

#endif /* defined(HAVE_MPI) */

  /* Set fine values */

  coarse_cell = c->coarse_cell;

  for (ii = 0; ii < f_n_cells; ii++)
    f_var[ii] = _c_var[coarse_cell[ii] - 1];
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
                         fvm_lnum_t        n_base_cells,
                         int               max_num,
                         int               c_cell_num[])
{
  fvm_lnum_t ii = 0;
  fvm_gnum_t base_shift = 1;
  fvm_gnum_t _max_num = max_num;
  fvm_lnum_t n_max_cells = 0;
  fvm_lnum_t *tmp_num_1 = NULL, *tmp_num_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(c_cell_num != NULL);

  /* Initialize array */

  n_max_cells = g->n_cells;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_cells > n_max_cells)
      n_max_cells = _g->n_cells;
  }

  BFT_MALLOC(tmp_num_1, n_max_cells, fvm_lnum_t);

  /* Compute local base starting cell number in parallel mode */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    fvm_gnum_t local_shift = g->n_cells;
    fvm_gnum_t global_shift = 0;
    MPI_Scan(&local_shift, &global_shift, 1, FVM_MPI_GNUM, MPI_SUM,
             cs_glob_mpi_comm);
    base_shift = 1 + global_shift - g->n_cells;
  }
#endif

  for (ii = 0; ii < g->n_cells; ii++)
    tmp_num_1[ii] = (fvm_gnum_t)(ii + base_shift) % _max_num;

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_num_2, n_max_cells, fvm_lnum_t);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      fvm_lnum_t n_parent_cells = _g->parent->n_cells;

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
                          fvm_lnum_t        n_base_cells,
                          int               c_cell_rank[])
{
  fvm_lnum_t ii;
  fvm_lnum_t n_max_cells = 0;
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

      fvm_lnum_t n_parent_cells = _g->parent->n_cells;

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
                    fvm_lnum_t        n_base_cells,
                    const cs_real_t   c_var[],
                    cs_real_t         f_var[])
{
  fvm_lnum_t ii;
  fvm_lnum_t n_max_cells = 0;
  cs_real_t *tmp_var_1 = NULL, *tmp_var_2 = NULL;
  const cs_grid_t *_g = g;

  assert(g != NULL);
  assert(c_var != NULL || g->n_cells == 0);
  assert(f_var != NULL);

  n_max_cells = g->n_cells;
  for (_g = g; _g != NULL; _g = _g->parent) {
    if (_g->n_cells > n_max_cells)
      n_max_cells = _g->n_cells;
  }

  BFT_MALLOC(tmp_var_1, n_max_cells, cs_real_t);
  memcpy(tmp_var_1, c_var, g->n_cells*sizeof(cs_real_t));

  /* Project to finer levels */

  if (g->level > 0) {

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_var_2, n_max_cells, cs_real_t);

    for (_g = g; _g->level > 0; _g = _g->parent) {

      fvm_lnum_t n_parent_cells = _g->parent->n_cells;

      cs_grid_prolong_cell_var(_g,
                               _g->parent,
                               tmp_var_1,
                               tmp_var_2);

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_var_1[ii] = tmp_var_2[ii];

    }

    assert(_g->level == 0);
    assert(_g->n_cells == n_base_cells);

    /* Free temporary arrays */

    BFT_FREE(tmp_var_2);
  }

  memcpy(f_var, tmp_var_1, n_base_cells*sizeof(cs_real_t));

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
                         fvm_lnum_t        n_base_cells,
                         cs_real_t         diag_dom[])
{
  fvm_lnum_t ii, jj, face_id;

  cs_real_t *dd = NULL;

  assert(g != NULL);
  assert(diag_dom != NULL);

  if (g->level == 0)
    dd = diag_dom;
  else
    BFT_MALLOC(dd, g->n_cells_ext, cs_real_t);

  /* Compute coarse diagonal dominance */

  {
    const fvm_lnum_t n_cells = g->n_cells;
    const fvm_lnum_t n_faces = g->n_faces;
    const fvm_lnum_t *face_cel = g->face_cell;

    /* Diagonal part of matrix.vector product */

    for (ii = 0; ii < n_cells; ii++)
      dd[ii] = fabs(g->da[ii]);

    if (g->halo != NULL)
      cs_halo_sync_var(g->halo, CS_HALO_STANDARD, dd);

    if (g->symmetric) {
      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = face_cel[2*face_id] -1;
        jj = face_cel[2*face_id + 1] -1;
        dd[ii] -= fabs(g->xa[face_id]);
        dd[jj] -= fabs(g->xa[face_id]);
      }
    }
    else {
      for (face_id = 0; face_id < n_faces; face_id++) {
        ii = face_cel[2*face_id] -1;
        jj = face_cel[2*face_id + 1] -1;
        dd[ii] -= fabs(g->xa[face_id*2]);
        dd[jj] -= fabs(g->xa[face_id*2 + 1]);
      }
    }

    for (ii = 0; ii < n_cells; ii++) {
      if (fabs(g->da[ii]) > 1.e-18)
        dd[ii] /= fabs(g->da[ii]);
    }

  }

  /* Now project to finer levels */

  if (dd != diag_dom) {
    cs_grid_project_var(g, n_base_cells, dd, diag_dom);
    BFT_FREE(dd);
  }
}

/*----------------------------------------------------------------------------
 * Get the default parameters for multigrid coarsening.
 *
 * parameters:
 *   merge_mean_threshold --> mean number of cells under which merging
 *                            should take place, or NULL
 *   merge_glob_threshold --> global number of cells under which merging
 *                            should take place, or NULL
 *   merge_min_ranks      --> number of active ranks under which no merging
 *                            takes place, or NULL
 *   merge_stride         --> number of ranks over which merging takes place,
 *                            or NULL
 *----------------------------------------------------------------------------*/

void
cs_grid_get_defaults(int  *merge_mean_threshold,
                     int  *merge_glob_threshold,
                     int  *merge_min_ranks,
                     int  *merge_stride)
{
#if defined(HAVE_MPI)
  if (merge_mean_threshold != NULL)
    *merge_mean_threshold = _grid_merge_threshold[0];
  if (merge_glob_threshold != NULL)
    *merge_glob_threshold = _grid_merge_threshold[1];
  if (merge_min_ranks != NULL)
    *merge_min_ranks = _grid_merge_min_ranks;
  if (merge_stride != NULL)
    *merge_stride = _grid_merge_stride;
#else
  if (merge_mean_threshold != NULL)
    *merge_mean_threshold = 0;
  if (merge_glob_threshold != NULL)
    *merge_glob_threshold = 0;
  if (merge_min_ranks != NULL)
    *merge_min_ranks = 0;
  if (merge_stride != NULL)
    *merge_stride = 0;
#endif
}

/*----------------------------------------------------------------------------
 * Set the default parameters for multigrid coarsening.
 *
 * parameters:
 *   merge_mean_threshold <-- mean number of cells under which merging
 *                            should take place
 *   merge_glob_threshold <-- global number of cells under which merging
 *                            should take place
 *   merge_min_ranks      <-- number of active ranks under which no merging
 *                            takes place
 *   merge_stride         <-- number of ranks over which merging takes place
 *----------------------------------------------------------------------------*/

void
cs_grid_set_defaults(int  merge_mean_threshold,
                     int  merge_glob_threshold,
                     int  merge_min_ranks,
                     int  merge_stride)
{
#if defined(HAVE_MPI)
  _grid_merge_threshold[0] = merge_mean_threshold;
  _grid_merge_threshold[1] = merge_glob_threshold;
  _grid_merge_min_ranks = merge_min_ranks;
  _grid_merge_stride = merge_stride;
#endif
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
 * Print the default parameters for multigrid coarsening.
 *----------------------------------------------------------------------------*/

void
cs_grid_log_defaults(void)
{
#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1)
    bft_printf(_("\n"
                 "  Multigrid rank merge parameters:\n"
                 "    mean  coarse cells merge threshold: %d\n"
                 "    total coarse cells merge threshold: %d\n"
                 "    minimum ranks merge threshold:      %d\n"
                 "    merge stride:                       %d\n"),
               _grid_merge_threshold[0], _grid_merge_threshold[1],
               _grid_merge_min_ranks, _grid_merge_stride);
#endif
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
  fvm_lnum_t  i;

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
             g, g->level, g->parent,
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
             g->face_cell, g->_face_cell, g->coarse_cell, g->coarse_face,
             g->halo);

  if (g->face_cell != NULL) {
    bft_printf("\n"
               "  face -> cell connectivity;\n");
    for (i = 0; i < g->n_faces; i++)
      bft_printf("    %d : %d, %d\n", (int)(i+1),
                 (int)(g->face_cell[i*2]), (int)(g->face_cell[i*2+1]));
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
