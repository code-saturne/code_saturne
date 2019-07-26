/*============================================================================
 * Default Sparse Matrix structure and Tuning.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_blas.h"
#include "cs_field.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_internal_coupling.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_numbering.h"
#include "cs_prototypes.h"
#include "cs_range_set.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"
#include "cs_matrix_tuning.h"

#include "cs_matrix_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*
  Tuned matrix structures re-used for various resolutions.

  These structures are kept throughout the whole run, to avoid paying the
  CPU overhead for their construction at each system resolution
  (at the cost of extra memory use, depending on the chosen structure).
*/

static bool _initialized = false;

/* _matrix_variant_tuned[mft] may be defined for a given fill type mft
   when variants are defined,  but are merged if possible
   upon calling cs_matrix_initialize(), so access for a given fill
   type is done using _matrix_variant_tuned[tuned_matrix_id[mft]]
   after that call */

static cs_matrix_variant_t *_matrix_variant_tuned[CS_MATRIX_N_FILL_TYPES];

static cs_matrix_structure_t *_matrix_struct_tuned[CS_MATRIX_N_FILL_TYPES];
static cs_matrix_t *_matrix_tuned[CS_MATRIX_N_FILL_TYPES];

/* _tuned_matrix_id[mft] is initialized to -1, and may be set to -2
   (or reset -to -1) using cs_matrix_set_tuning(mft, tune) to indicate
   that autotuning is requested */

static int _tuned_matrix_id[CS_MATRIX_N_FILL_TYPES];

/* MSR matrix structure, if needed */

static cs_matrix_structure_t *_matrix_struct_msr = NULL;
static cs_matrix_t *_matrix_msr = NULL;

/* Native matrix structure, if needed */

static cs_matrix_structure_t *_matrix_struct_native = NULL;
static cs_matrix_t *_matrix_native = NULL;

/* Tuning options */

static double _t_measure = 0.5;
static int _n_min_products = 50;

/* Pointer to global (block-based) numbering, if used */

static cs_lnum_t  _row_num_size = 0;
static cs_gnum_t  *_global_row_id = NULL;
static cs_gnum_t  _l_range[2] = {0, 0};

/* Pointer to internal coupling oriented matrix structures */

static cs_matrix_assembler_t  **_matrix_assembler_coupled = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize sparse matrix API.
 *----------------------------------------------------------------------------*/

static void
_initialize_api(void)
{
  if (! _initialized) {
    for (cs_matrix_fill_type_t mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++) {
      _matrix_variant_tuned[mft] = NULL;
      _matrix_struct_tuned[mft] = NULL;
      _matrix_tuned[mft] = NULL;
      _tuned_matrix_id[mft] = -1;
    }
    _matrix_struct_msr = NULL;
    _matrix_msr = NULL;
    _matrix_struct_native = NULL;
    _matrix_native = NULL;
    _initialized = true;
  }
}

/*----------------------------------------------------------------------------
 * Build a global block row numbering.
 *
 * parameters:
 *   n_rows <-- associated number of local rows
 *   halo   <-- associated halo, or NULL
 *----------------------------------------------------------------------------*/

static void
_build_block_row_g_id(cs_lnum_t         n_rows,
                      const cs_halo_t  *halo)
{
  cs_lnum_t _n_rows = n_rows;

  _row_num_size = n_rows;
  if (halo != NULL) {
    assert(n_rows == halo->n_local_elts);
    _n_rows += halo->n_elts[CS_HALO_EXTENDED];
  }

  BFT_REALLOC(_global_row_id, _n_rows, cs_gnum_t);

  cs_range_set_define(NULL,
                      halo,
                      n_rows,
                      false,
                      0, /* g_id_base */
                      _l_range,
                      _global_row_id);
}

/*----------------------------------------------------------------------------
 * Create a matrix structure using a matrix assembler
 *
 * parameters:
 *   coupling_id <-- internal coupling id, or -1
 *
 * returns:
 *   pointer to created matrix assembler
 *----------------------------------------------------------------------------*/

static cs_matrix_assembler_t *
_create_assembler(int  coupling_id)
{
  const cs_mesh_t *m = cs_glob_mesh;

  const cs_lnum_t     n_rows = m->n_cells;
  const cs_lnum_t     n_edges = m->n_i_faces;
  const cs_lnum_2_t  *edges = (const cs_lnum_2_t *)(m->i_face_cells);

  /* Global cell ids, based on range/scan */

  if (_global_row_id == NULL)
    _build_block_row_g_id(n_rows, m->halo);

  const cs_gnum_t *r_g_id = _global_row_id;
  cs_gnum_t l_range[2] = {_l_range[0], _l_range[1]};

  /* Build matrix assembler;
     assemble by blocks to amortize overhead */

  assert(edges != NULL || n_edges == 0);

  cs_matrix_assembler_t  *ma
    = cs_matrix_assembler_create(l_range, true);

#if 0 /* TODO: test and check performance with this flag */
  int ma_flags = CS_MATRIX_DISTANT_ROW_USE_COL_IDX;
#else
  int ma_flags = 0;
#endif

  cs_matrix_assembler_set_options(ma, ma_flags);

  /* First, add diagonal terms */

  cs_matrix_assembler_add_g_ids(ma, n_rows, r_g_id, r_g_id);

  /* Then add standard local off-diagonal terms */

  {
    const cs_lnum_t block_size = 800;
    cs_gnum_t g_row_id[800];
    cs_gnum_t g_col_id[800];

    cs_lnum_t jj = 0;
    for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
      cs_lnum_t i0 = edges[ii][0];
      cs_lnum_t i1 = edges[ii][1];
      if (i0 < n_rows) {
        g_row_id[jj] = r_g_id[i0];
        g_col_id[jj] = r_g_id[i1];
        jj++;
      }
      if (i1 < n_rows) {
        g_row_id[jj] = r_g_id[i1];
        g_col_id[jj] = r_g_id[i0];
        jj++;
      }
      if (jj >= block_size - 1) {
        cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
        jj = 0;
      }
    }
    if (jj > 0)
      cs_matrix_assembler_add_g_ids(ma, jj, g_row_id, g_col_id);
  }

  if (coupling_id > -1)
    cs_internal_coupling_matrix_add_ids(coupling_id, r_g_id,  ma);

  /* Now compute structure */
  cs_matrix_assembler_compute(ma);

  return ma;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Matrix (native format) vector product
 *
 * parameters:
 *   symmetric     <-- Symmetry indicator:
 *   db_size       <-- block sizes for diagonal
 *   eb_size       <-- block sizes for extra diagonal
 *   rotation_mode <-- halo update option for rotational periodicity
 *   f_id          <-- associated field id, or < 0
 *   dam           <-- Matrix diagonal
 *   xam           <-- Matrix extra-diagonal terms
 *   vx            <-- A*vx
 *   vy            <-> vy = A*vx
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_native_multiply(bool                symmetric,
                                 const int           db_size[4],
                                 const int           eb_size[4],
                                 cs_halo_rotation_t  rotation_mode,
                                 int                 f_id,
                                 const cs_real_t    *dam,
                                 const cs_real_t    *xam,
                                 cs_real_t          *vx,
                                 cs_real_t          *vy)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_matrix_t *a;

  a = cs_matrix_native(symmetric, db_size, eb_size);

  cs_matrix_set_coefficients(a,
                             symmetric,
                             db_size,
                             eb_size,
                             m->n_i_faces,
                             (const cs_lnum_2_t *)m->i_face_cells,
                             dam,
                             xam);

  cs_matrix_vector_multiply(rotation_mode,
                            a,
                            vx,
                            vy);

  /* Add extended contribution for domain coupling */

  if (f_id != -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    int coupling_id = cs_field_get_key_int(f,
                                           cs_field_key_id("coupling_entity"));

    if (coupling_id > -1)
      cs_internal_coupling_spmv_contribution(false,
                                             f,
                                             vx,
                                             vy);
  }

}

/*----------------------------------------------------------------------------
 * Initialize sparse matrix API.
 *----------------------------------------------------------------------------*/

void
cs_matrix_initialize(void)
{
  cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

  int n_tuned_types = 0;
  bool matrix_tune = false;

  assert(mesh != NULL);

  if (!_initialized)
    _initialize_api();

  /* Compute tuned variants for matrix */

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {

    cs_matrix_variant_t *mv = _matrix_variant_tuned[i];

    _matrix_variant_tuned[i] = NULL;

    if (mv == NULL) {

      if (_tuned_matrix_id[i] < -1) {

        matrix_tune = true;

        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("\n"
                        "Tuning for matrices of type: %s\n"
                        "===========================\n"),
                      cs_matrix_fill_type_name[i]);

        int n_fill_types = 1;
        cs_matrix_fill_type_t fill_types[1] = {i};
        double fill_weights[1] = {1};

        mv = cs_matrix_variant_tuned(_t_measure,
                                     0, /* n_matrix_types, */
                                     n_fill_types,
                                     NULL, /* matrix_types, */
                                     fill_types,
                                     fill_weights,
                                     mesh->n_cells,
                                     mesh->n_cells_with_ghosts,
                                     mesh->n_i_faces,
                                     (const cs_lnum_2_t *)(mesh->i_face_cells),
                                     mesh->halo,
                                     mesh->i_face_numbering);

      }

      else {

        cs_matrix_type_t m_type = CS_MATRIX_NATIVE;

        mv = cs_matrix_variant_create(m_type,
                                      mesh->i_face_numbering);

      }

    }

    /* Prepare to share matrix variants and structures if possible */

    int m_id = -1;
    cs_matrix_type_t m_type = cs_matrix_variant_type(mv);

    for (int j = 0; j < n_tuned_types; j++) {
      if (m_type == _matrix_struct_tuned[j]->type) {
        m_id = j;
        cs_matrix_variant_merge(_matrix_variant_tuned[m_id], mv, i);
        _tuned_matrix_id[i] = j;
        cs_matrix_variant_destroy(&mv);
        break;
      }
    }

    /* Build new structure otherwise */

    if (m_id < 0) {

      m_id = n_tuned_types;

      _matrix_variant_tuned[m_id] = mv;

      _tuned_matrix_id[i] = m_id;

      if (m_type == CS_MATRIX_MSR && ma != NULL)
        _matrix_struct_tuned[m_id]
          = cs_matrix_structure_create_msr_shared(true,
                                                  ma->single_faces_to_cells,
                                                  mesh->n_cells,
                                                  mesh->n_cells_with_ghosts,
                                                  ma->cell_cells_idx,
                                                  ma->cell_cells,
                                                  mesh->halo,
                                                  mesh->i_face_numbering);

      else
        _matrix_struct_tuned[m_id]
          = cs_matrix_structure_create(m_type,
                                       true,
                                       mesh->n_cells,
                                       mesh->n_cells_with_ghosts,
                                       mesh->n_i_faces,
                                       (const cs_lnum_2_t *)(mesh->i_face_cells),
                                       mesh->halo,
                                       mesh->i_face_numbering);

      _matrix_tuned[m_id]
        = cs_matrix_create_by_variant(_matrix_struct_tuned[m_id], mv);

      n_tuned_types += 1;

    }

  }

  if (matrix_tune > 0) {
    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);
  }

  /* Matrices for internal couplings */

  int n_ic = cs_internal_coupling_n_couplings();

  if (n_ic > 0) {

    BFT_MALLOC(_matrix_assembler_coupled, n_ic, cs_matrix_assembler_t *);

    for (int i = 0; i < n_ic; i++)
      _matrix_assembler_coupled[i] = _create_assembler(i);

  }
}

/*----------------------------------------------------------------------------
 * Finalize sparse matrix API.
 *----------------------------------------------------------------------------*/

void
cs_matrix_finalize(void)
{
  BFT_FREE(_global_row_id);

  for (cs_matrix_fill_type_t mft = 0; mft < CS_MATRIX_N_FILL_TYPES; mft++)
    _tuned_matrix_id[mft] = -1;

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    if (_matrix_tuned[i] != NULL)
      cs_matrix_destroy(&(_matrix_tuned[i]));
    if (_matrix_struct_tuned[i] != NULL)
      cs_matrix_structure_destroy(&(_matrix_struct_tuned[i]));
    if (_matrix_variant_tuned[i] != NULL)
      cs_matrix_variant_destroy(&(_matrix_variant_tuned[i]));
  }

  if (_matrix_msr != NULL)
    cs_matrix_destroy(&(_matrix_msr));
  if (_matrix_struct_msr != NULL)
    cs_matrix_structure_destroy(&(_matrix_struct_msr));

  if (_matrix_native != NULL)
    cs_matrix_destroy(&(_matrix_native));
  if (_matrix_struct_native != NULL)
    cs_matrix_structure_destroy(&(_matrix_struct_native));

  /* Matrices for internal couplings */

  int n_ic = cs_internal_coupling_n_couplings();
  for (int i = 0; i < n_ic; i++) {
    cs_matrix_assembler_destroy(&(_matrix_assembler_coupled[i]));
  }
  BFT_FREE(_matrix_assembler_coupled);

  /* Exit status */

  _initialized = false;
  _initialize_api();
  _initialized = false;
}

/*----------------------------------------------------------------------------
 * Update sparse matrix API in case of mesh modification.
 *----------------------------------------------------------------------------*/

void
cs_matrix_update_mesh(void)
{
  const cs_mesh_t  *mesh = cs_glob_mesh;
  const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

  if (_global_row_id != NULL)
    _build_block_row_g_id(mesh->n_cells, mesh->halo);

  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {

    if (_matrix_tuned[i] != NULL) {

      const cs_matrix_type_t m_type = (_matrix_tuned[i])->type;

      cs_matrix_destroy(&(_matrix_tuned[i]));
      cs_matrix_structure_destroy(&(_matrix_struct_tuned[i]));

      if (m_type == CS_MATRIX_MSR && ma != NULL)
        _matrix_struct_tuned[i]
          = cs_matrix_structure_create_msr_shared(true,
                                                  ma->single_faces_to_cells,
                                                  mesh->n_cells,
                                                  mesh->n_cells_with_ghosts,
                                                  ma->cell_cells_idx,
                                                  ma->cell_cells,
                                                  mesh->halo,
                                                  mesh->i_face_numbering);

      else
        _matrix_struct_tuned[i]
          = cs_matrix_structure_create(m_type,
                                       true,
                                       mesh->n_cells,
                                       mesh->n_cells_with_ghosts,
                                       mesh->n_i_faces,
                                       (const cs_lnum_2_t *)(mesh->i_face_cells),
                                       mesh->halo,
                                       mesh->i_face_numbering);

      assert(_matrix_variant_tuned[i] != NULL);

      _matrix_tuned[i]
        = cs_matrix_create_by_variant(_matrix_struct_tuned[i],
                                      _matrix_variant_tuned[i]);

    }

  }

  /* MSR might also be required separately */

  if (_matrix_msr != NULL) {

    cs_matrix_destroy(&(_matrix_msr));
    cs_matrix_structure_destroy(&(_matrix_struct_msr));

    if (ma != NULL)
      _matrix_struct_msr
        = cs_matrix_structure_create_msr_shared(true,
                                                ma->single_faces_to_cells,
                                                mesh->n_cells,
                                                mesh->n_cells_with_ghosts,
                                                ma->cell_cells_idx,
                                                ma->cell_cells,
                                                mesh->halo,
                                                mesh->i_face_numbering);
    else
      _matrix_struct_msr
        = cs_matrix_structure_create(CS_MATRIX_MSR,
                                     true,
                                     mesh->n_cells,
                                     mesh->n_cells_with_ghosts,
                                     mesh->n_i_faces,
                                     (const cs_lnum_2_t *)(mesh->i_face_cells),
                                     mesh->halo,
                                     mesh->i_face_numbering);

    _matrix_msr = cs_matrix_create(_matrix_struct_msr);

  }

  /* Same for native... */

  if (_matrix_native != NULL) {

    cs_matrix_destroy(&(_matrix_native));
    cs_matrix_structure_destroy(&(_matrix_struct_native));

    _matrix_struct_native
      = cs_matrix_structure_create(CS_MATRIX_NATIVE,
                                   true,
                                   mesh->n_cells,
                                   mesh->n_cells_with_ghosts,
                                   mesh->n_i_faces,
                                   (const cs_lnum_2_t *)(mesh->i_face_cells),
                                   mesh->halo,
                                   mesh->i_face_numbering);

    _matrix_native = cs_matrix_create(_matrix_struct_native);

  }

  /* Matrices for internal couplings */

  int n_ic = cs_internal_coupling_n_couplings();

  for (int i = 0; i < n_ic; i++) {
    cs_matrix_assembler_destroy(&(_matrix_assembler_coupled[i]));
    _matrix_assembler_coupled[i] = _create_assembler(i);
  }
}

/*----------------------------------------------------------------------------
 * Return default matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to default matrix structure adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_default(bool        symmetric,
                  const int  *diag_block_size,
                  const int  *extra_diag_block_size)
{
  cs_matrix_t *m = NULL;

  cs_matrix_fill_type_t mft = cs_matrix_get_fill_type(symmetric,
                                                      diag_block_size,
                                                      extra_diag_block_size);

  if (_tuned_matrix_id[mft] > -1)
    m = _matrix_tuned[_tuned_matrix_id[mft]];

  return m;
}

/*----------------------------------------------------------------------------
 * Return MSR matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to MSR matrix adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_msr(bool        symmetric,
              const int  *diag_block_size,
              const int  *extra_diag_block_size)
{
  cs_matrix_t *m = NULL;

  /* If default matrix for fill type is already MSR, return that */

  cs_matrix_fill_type_t mft = cs_matrix_get_fill_type(symmetric,
                                                      diag_block_size,
                                                      extra_diag_block_size);

  if (_matrix_tuned[mft] != NULL) {
    if ((_matrix_tuned[mft])->type == CS_MATRIX_MSR)
      m = cs_matrix_default(symmetric,
                            diag_block_size,
                            extra_diag_block_size);
  }

  if (m == NULL) {

    /* Create matrix if not done yet */

    if (_matrix_msr == NULL) {

      const cs_mesh_t  *mesh = cs_glob_mesh;
      const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

      if (ma != NULL)
        _matrix_struct_msr
          = cs_matrix_structure_create_msr_shared(true,
                                                  ma->single_faces_to_cells,
                                                  mesh->n_cells,
                                                  mesh->n_cells_with_ghosts,
                                                  ma->cell_cells_idx,
                                                  ma->cell_cells,
                                                  mesh->halo,
                                                  mesh->i_face_numbering);
      else
        _matrix_struct_msr
          = cs_matrix_structure_create(CS_MATRIX_MSR,
                                       true,
                                       mesh->n_cells,
                                       mesh->n_cells_with_ghosts,
                                       mesh->n_i_faces,
                                       (const cs_lnum_2_t *)(mesh->i_face_cells),
                                       mesh->halo,
                                       mesh->i_face_numbering);

      _matrix_msr = cs_matrix_create(_matrix_struct_msr);

    }

    m = _matrix_msr;

  }

  return m;
}

/*----------------------------------------------------------------------------
 * Return native matrix for a given fill type
 *
 * parameters:
 *   symmetric              <-- Indicates if matrix coefficients are symmetric
 *   diag_block_size        <-- Block sizes for diagonal, or NULL
 *   extra_diag_block_size  <-- Block sizes for extra diagonal, or NULL
 *
 * returns:
 *   pointer to native matrix adapted to fill type
 *----------------------------------------------------------------------------*/

cs_matrix_t  *
cs_matrix_native(bool        symmetric,
                 const int  *diag_block_size,
                 const int  *extra_diag_block_size)
{
  cs_matrix_t *m = NULL;

  /* If default matrix for fill type is already native, return that */

  cs_matrix_fill_type_t mft = cs_matrix_get_fill_type(symmetric,
                                                      diag_block_size,
                                                      extra_diag_block_size);

  if (_matrix_tuned[_tuned_matrix_id[mft]] != NULL) {
    if ((_matrix_tuned[_tuned_matrix_id[mft]])->type == CS_MATRIX_NATIVE)
      m = cs_matrix_default(symmetric,
                            diag_block_size,
                            extra_diag_block_size);
  }

  if (m == NULL) {

    /* Create matrix if not done yet */

    if (_matrix_native == NULL) {

      cs_mesh_t  *mesh = cs_glob_mesh;

      _matrix_struct_native
        = cs_matrix_structure_create(CS_MATRIX_NATIVE,
                                     true,
                                     mesh->n_cells,
                                     mesh->n_cells_with_ghosts,
                                     mesh->n_i_faces,
                                     (const cs_lnum_2_t *)(mesh->i_face_cells),
                                     mesh->halo,
                                     mesh->i_face_numbering);

      _matrix_native = cs_matrix_create(_matrix_struct_native);

    }

    m = _matrix_native;

  }

  return m;
}

/*----------------------------------------------------------------------------
 * Force matrix variant for a given fill type
 *
 * Information from the variant used fo this definition is copied,
 * so it may be freed after calling this function.
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *   mv         <-- Matrix variant to use for this type
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_variant(cs_matrix_fill_type_t       fill_type,
                      const cs_matrix_variant_t  *mv)
{
  if (!_initialized)
    _initialize_api();

  /* Create default variant for copy if none present */

  if (_matrix_variant_tuned[fill_type] == NULL) {
    cs_matrix_type_t m_type = cs_matrix_variant_type(mv);
    _matrix_variant_tuned[fill_type] = cs_matrix_variant_create(m_type,
                                                                NULL);
  }

  cs_matrix_variant_t *_mv = _matrix_variant_tuned[fill_type];
  cs_matrix_variant_merge(_mv, mv, fill_type);
}

/*----------------------------------------------------------------------------
 * Set matrix tuning behavior for a given fill type
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *   tune       <-- 1 to activate tuning, 0 to deactivate
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_tuning(cs_matrix_fill_type_t   fill_type,
                     int                     tune)
{
  if (!_initialized)
    _initialize_api();

  if (_tuned_matrix_id[fill_type] < 0) {
    if (tune)
      _tuned_matrix_id[fill_type] = -2;
    else
      _tuned_matrix_id[fill_type] = -1;
  }
}

/*----------------------------------------------------------------------------
 * Return matrix tuning behavior for a given fill type.
 *
 * parameters:
 *   fill type  <-- Fill type for which tuning behavior is set
 *
 * returns:
 *   1 if tuning is active, 0 otherwise
 *----------------------------------------------------------------------------*/

int
cs_matrix_get_tuning(cs_matrix_fill_type_t   fill_type)
{
  int retval = 0;

  if (!_initialized)
    _initialize_api();

  if (_tuned_matrix_id[fill_type] < -1)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------
 * Set number of matrix computation runs for tuning.
 *
 * If this function is not called, defaults are:
 *  - minimum of 10 runs
 *  - minimum of 0.5 seconds of running
 *
 * parameters:
 *   n_min_products <-- minimum number of expected SpM.V products for
 *                      coefficients assign amortization.
 *   t_measure      <-- minimum running time per measure
 *----------------------------------------------------------------------------*/

void
cs_matrix_set_tuning_runs(int     n_min_products,
                          double  t_measure)
{
  if (!_initialized)
    _initialize_api();

  _n_min_products = n_min_products;
  _t_measure = t_measure;
}

/*----------------------------------------------------------------------------
 * Get number of matrix computation runs for tuning.
 *
 * parameters:
 *   n_min_products --> minimum number of expected SpM.V products for
 *                      coefficients assign amortization.
 *   t_measure      --> minimum running time per measure, or NULL
 *----------------------------------------------------------------------------*/

void
cs_matrix_get_tuning_runs(int     *n_min_products,
                          double  *t_measure)
{
  if (!_initialized)
    _initialize_api();

  if (n_min_products != NULL)
    *n_min_products = _n_min_products;

  if (t_measure != NULL)
    *t_measure = _t_measure;
}

/*----------------------------------------------------------------------------
 * Return a (0-based) global block row numbering.
 *
 * The numbering is built if not previously present, and returned otherwise.
 *
 * Currently, the function only handles one n_rows/halo combination, and does
 * not check for consistency.
 *
 * parameters:
 *   n_rows <-- associated number of local rows
 *   halo   <-- associated halo, or NULL
 *
 * returns:
 *   pointer to requested global numbering
 *----------------------------------------------------------------------------*/

const cs_gnum_t *
cs_matrix_get_block_row_g_id(cs_lnum_t         n_rows,
                             const cs_halo_t  *halo)
{
  const cs_gnum_t  *g_row_num = _global_row_id;

  if (_global_row_id == NULL || n_rows > _row_num_size) {
    _build_block_row_g_id(n_rows, halo);
    g_row_num = _global_row_id;
  }

  return g_row_num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign coefficients to a matrix using a matrix assembler.
 *
 * \param[in]  f                      pointer to associated field
 * \param[in]  type                   matrix type
 * \param[in]  symmetric              is matrix symmetric ?
 * \param[in]  diag_block_size        block sizes for diagonal, or NULL
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal, or NULL
 * \param[in]  da                     diagonal values (NULL if zero)
 * \param[in]  xa                     extradiagonal values (NULL if zero)
 *                                    casts as:
 *                                      xa[n_edges]    if symmetric,
 *                                      xa[n_edges][2] if non symmetric
 *
 * \return  pointer to associated matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_matrix_set_coefficients_coupled(const cs_field_t  *f,
                                   cs_matrix_type_t   type,
                                   bool               symmetric,
                                   const int         *diag_block_size,
                                   const int         *extra_diag_block_size,
                                   const cs_real_t   *da,
                                   const cs_real_t   *xa)
{
  int coupling_id = cs_field_get_key_int(f,
                                         cs_field_key_id("coupling_entity"));

  assert(coupling_id > -1);  /* Only reason for this restriction is
                                storage/access to assembler */

  const cs_mesh_t *mesh = cs_glob_mesh;

  const cs_lnum_t     n_rows = mesh->n_cells;
  const cs_lnum_t     n_edges = mesh->n_i_faces;
  const cs_lnum_2_t  *edges = (const cs_lnum_2_t *)(mesh->i_face_cells);

  cs_lnum_t s0 = 2;
  cs_lnum_t s1 = 1;

  if (symmetric) {
    s0 = 1;
    s1 = 0;
  }

  cs_matrix_assembler_t  *ma = _matrix_assembler_coupled[coupling_id];

  cs_matrix_t *m = cs_matrix_create_from_assembler(type, ma);

  cs_matrix_assembler_values_t *mav
    = cs_matrix_assembler_values_init(m,
                                      diag_block_size,
                                      extra_diag_block_size);

  /* Range information already built for assembler */

  assert(n_rows == cs_matrix_get_n_rows(m));

  const cs_gnum_t *r_g_id = _global_row_id;

  /* Set coefficients */

  const cs_lnum_t block_size = 800;
  cs_gnum_t g_row_id[800];
  cs_gnum_t g_col_id[800];
  cs_real_t val[1600];

  /* Diagonal values */

  cs_matrix_assembler_values_add_g(mav, n_rows, r_g_id, r_g_id, da);

  /* Extradiagonal values based on internal faces */

  cs_lnum_t db_size = 1;
  if (diag_block_size != NULL)
    db_size = diag_block_size[0];

  cs_lnum_t eb_size = 1;
  if (extra_diag_block_size != NULL)
    eb_size = extra_diag_block_size[0];

  cs_lnum_t jj = 0;

  if (eb_size == 1) {
    for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
      cs_lnum_t i0 = edges[ii][0];
      cs_lnum_t i1 = edges[ii][1];
      if (i0 < n_rows) {
        g_row_id[jj] = r_g_id[i0];
        g_col_id[jj] = r_g_id[i1];
        val[jj] = xa[ii*s0];
        jj++;
      }
      if (i1 < n_rows) {
        g_row_id[jj] = r_g_id[i1];
        g_col_id[jj] = r_g_id[i0];
        val[jj] = xa[ii*s0+s1];
        jj++;
      }
      if (jj >= block_size - 1) {
        cs_matrix_assembler_values_add_g(mav, jj,
                                         g_row_id, g_col_id, val);
        jj = 0;
      }
    }
    cs_matrix_assembler_values_add_g(mav, jj,
                                     g_row_id, g_col_id, val);
    jj = 0;
  }
  else {
    assert(0); /* TODO handle extra-diagonal blocks of size > 1 */
  }

  /* Set extended contribution for domain coupling */

  cs_internal_coupling_matrix_add_values(f,
                                         db_size,
                                         eb_size,
                                         r_g_id,
                                         mav);

  /* Finalize assembly */

  cs_matrix_assembler_values_finalize(&mav);

  return m;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
