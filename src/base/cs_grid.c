/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Grid connectivity and data used for multgrid coarsening
 * and associated matrix construction.
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_defs.h>
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

  cs_matrix_t      *matrix;         /* Associated matrix structure */
};

/*============================================================================
 *  Global variables
 *============================================================================*/

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

  return g;
}

/*----------------------------------------------------------------------------
 * Initialize coarse grid from fine grid
 *
 * This creates au quasi-empty grid structure, with only symmetry and
 * level information initialized, and coarsening array allocated and
 * zeroed.
 *
 * After this function is called, the coarsening array mus be determined,
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

  /* Build face coarsening and coarse grid face -> cells connectivity */

  _coarsen_faces(f,
                 c->coarse_cell,
                 c->n_cells,
                 &(c->n_faces),
                 &(c->coarse_face),
                 &(c->_face_cell));

  c->face_cell = c->_face_cell;
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
    const cs_real_t *xa_sym = xa + n_faces;
    BFT_MALLOC(g->_xa0, n_faces, cs_real_t);
    for (face_id = 0; face_id < n_faces; face_id++)
      g->_xa0[face_id] = 0.5 * (xa[face_id] + xa_sym[face_id]);
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

  g->matrix = cs_matrix_create(CS_MATRIX_NATIVE,
                               symmetric,
                               true,
                               false, /* No periodicity here yet */
                               n_cells,
                               n_cells_ext,
                               n_faces,
                               NULL,
                               face_cell,
                               halo,
                               numbering);

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
 *   n_cells_ext --> Number of local cells (or NULL)
 *   n_cells_ext --> Number of cells including ghosts (or NULL)
 *   n_cells_ext --> Number of faces (or NULL)
 *   n_g_cells   --> Number of global cells (or NULL)
 *----------------------------------------------------------------------------*/

void
cs_grid_get_info(const cs_grid_t  *g,
                 int              *level,
                 cs_bool_t        *symmetric,
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
cs_grid_coarsen(const cs_grid_t  *f,
                int               verbosity,
                int               agglomeration_limit,
                int              *max_agglomeration)
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
                            f->cell_vol, f->cell_cen, f->face_normal,
                            f->xa0, f->xa0ij, f->da, f->xa,
                            c->cell_vol, c->cell_cen, c->face_normal,
                            c->_xa0, c->xa0ij, c->_da, c->_xa,
                            rwc1, rwc2, rwc3, rwc4);

  /* Synchronize matrix's geometric quantities */

  if (c->halo != NULL)
    cs_halo_sync_var(c->halo, CS_HALO_STANDARD, c->_da);

  BFT_FREE(rwc1);
  rwc2 = NULL;
  rwc3 = NULL;
  rwc4 = NULL;

  c->matrix = cs_matrix_create(CS_MATRIX_NATIVE,
                               c->symmetric,
                               true,
                               false, /* No periodicity here yet */
                               c->n_cells,
                               c->n_cells_ext,
                               c->n_faces,
                               NULL,
                               c->face_cell,
                               c->halo,
                               NULL);

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
 *
 * returns:
 *   coarse grid structure
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
  fvm_lnum_t c_n_cells_ext = c->n_cells_ext;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL);
  assert(f_var != NULL);
  assert(c_var != NULL);

  /* Set coarse values */

  coarse_cell = c->coarse_cell;

  for (ii = 0; ii < c_n_cells_ext; ii++)
    c_var[ii] = 0.0;

  for (ii = 0; ii < f_n_cells; ii++)
    c_var[coarse_cell[ii] - 1] += f_var[ii];
}

/*----------------------------------------------------------------------------
 * Compute fine cell variable values from coarse cell values
 *
 * parameters:
 *   c       <-- Fine grid structure
 *   f       <-- Fine grid structure
 *   c_var   --> Variable defined on coarse grid cells
 *   f_var   <-- Variable defined on fine grid cells
 *
 * returns:
 *   coarse grid structure
 *----------------------------------------------------------------------------*/

void
cs_grid_prolong_cell_var(const cs_grid_t  *c,
                         const cs_grid_t  *f,
                         const cs_real_t  *c_var,
                         cs_real_t        *f_var)
{
  fvm_lnum_t ii;
  const fvm_lnum_t *coarse_cell;

  fvm_lnum_t f_n_cells = f->n_cells;

  assert(f != NULL);
  assert(c != NULL);
  assert(c->coarse_cell != NULL);
  assert(f_var != NULL);
  assert(c_var != NULL);

  /* Set fine values */

  coarse_cell = c->coarse_cell;

  for (ii = 0; ii < f_n_cells; ii++)
    f_var[ii] = c_var[coarse_cell[ii] - 1];
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

  assert(g != NULL);
  assert(c_cell_num != NULL);

  /* Initialize array */

  for (ii = 0; ii < n_base_cells; ii++)
    c_cell_num[ii] = 0;

  if (g->level == 0)
    return;

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

  if (g->level == 1) {

    for (ii = 0; ii < n_base_cells; ii++)
      c_cell_num[ii]
        = ((fvm_gnum_t)(g->coarse_cell[ii]) + base_shift) % _max_num;

  }
  else { /* if g->level > 1) */

    fvm_lnum_t *tmp_num_1 = NULL, *tmp_num_2 = NULL;
    const cs_grid_t *_g = g;

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_num_1, n_base_cells, fvm_lnum_t);
    BFT_MALLOC(tmp_num_2, n_base_cells, fvm_lnum_t);

    for (ii = 0; ii < g->n_cells; ii++)
      tmp_num_1[ii] = ii;

    for (_g = g; _g->level > 1; _g = _g->parent) {

      fvm_lnum_t n_parent_cells = _g->parent->n_cells;

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_num_2[ii] = tmp_num_1[_g->coarse_cell[ii] - 1];

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_num_1[ii] = tmp_num_2[ii];

    }

    assert(_g->level == 1);
    assert(_g->parent->n_cells == n_base_cells);

    for (ii = 0; ii < n_base_cells; ii++)
      c_cell_num[ii]
        =   ((fvm_gnum_t)(tmp_num_1[_g->coarse_cell[ii] - 1]) + base_shift)
          % _max_num;

    /* Free temporary arrays */

    BFT_FREE(tmp_num_1);
    BFT_FREE(tmp_num_2);
  }

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

  assert(g != NULL);
  assert(c_var != NULL);
  assert(f_var != NULL);

  if (g->level == 0)
    memcpy(f_var, c_var, n_base_cells*sizeof(cs_real_t));

  /* Project to finer levels */

  else if (g->level == 1) {

    for (ii = 0; ii < n_base_cells; ii++)
      f_var[ii] = c_var[g->coarse_cell[ii] - 1];

  }
  else { /* if g->level > 1) */

    cs_real_t *tmp_var_1 = NULL, *tmp_var_2 = NULL;
    const cs_grid_t *_g = g;

    /* Allocate temporary arrays */

    BFT_MALLOC(tmp_var_1, n_base_cells, cs_real_t);
    BFT_MALLOC(tmp_var_2, n_base_cells, cs_real_t);

    for (ii = 0; ii < g->n_cells; ii++)
      tmp_var_1[ii] = c_var[ii];

    for (_g = g; _g->level > 1; _g = _g->parent) {

      fvm_lnum_t n_parent_cells = _g->parent->n_cells;

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_var_2[ii] = tmp_var_1[_g->coarse_cell[ii] - 1];

      for (ii = 0; ii < n_parent_cells; ii++)
        tmp_var_1[ii] = tmp_var_2[ii];

    }

    assert(_g->level == 1);
    assert(_g->parent->n_cells == n_base_cells);

    for (ii = 0; ii < n_base_cells; ii++)
      f_var[ii] = tmp_var_1[_g->coarse_cell[ii] - 1];

    /* Free temporary arrays */

    BFT_FREE(tmp_var_1);
    BFT_FREE(tmp_var_2);
  }
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
        dd[ii] -= fabs(g->xa[face_id]);
        dd[jj] -= fabs(g->xa[face_id + n_faces]);
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
