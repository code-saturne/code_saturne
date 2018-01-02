/*============================================================================
 * Default Sparse Matrix structure and Tuning.
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

  cs_gnum_t l_range[2];

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
                      l_range,
                      _global_row_id);
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
 *   dam           <-- Matrix diagonal
 *   xam           <-- Matrix extra-diagonal terms
 *   vx            <-- A*vx
 *   vy            <-> vy = A*vx
 *----------------------------------------------------------------------------*/

void
cs_matrix_vector_native_multiply(bool                symmetric,
                                 int                 db_size[4],
                                 int                 eb_size[4],
                                 cs_halo_rotation_t  rotation_mode,
                                 int                 f_id,
                                 const cs_real_t    *dam,
                                 const cs_real_t    *xam,
                                 cs_real_t          *vx,
                                 cs_real_t          *vy)
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_matrix_t *a;

  a = cs_matrix_native(symmetric,
                       db_size,
                       eb_size);

  cs_matrix_set_coefficients(a,
                             symmetric,
                             db_size,
                             eb_size,
                             m->n_i_faces,
                             (const cs_lnum_2_t *)m->i_face_cells,
                             dam,
                             xam);

  /* Set extended contribution for domain coupling */
  if (f_id != -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    int coupling_id = cs_field_get_key_int(f,
                                           cs_field_key_id("coupling_entity"));

    if (coupling_id > -1) {
      cs_matrix_set_extend
        (a,
         cs_internal_coupling_spmv_contribution,
         cs_matrix_preconditionning_add_coupling_contribution,
         f);
    }
  }

  cs_matrix_vector_multiply(rotation_mode,
                            a,
                            vx,
                            vy);
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
                                     _n_min_products,
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

END_C_DECLS
