/*============================================================================
 * Sparse Matrix Benchmarking of various matrix representations and operations..
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#if defined (HAVE_MKL)
#include <mkl_spblas.h>
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
#include "cs_log.h"
#include "cs_numbering.h"
#include "cs_parall.h"
#include "cs_prototypes.h"
#include "cs_timer.h"

#if defined(HAVE_HYPRE)
#include "cs_matrix_hypre.h"
#endif

#if defined(HAVE_PETSC)
#include "cs_matrix_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"

#include "cs_benchmark_matrix.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Note that most types are declared in cs_matrix_priv.h.
   only those only handled here are declared here. */

/* Structure used for timing variants */
/*------------------------------------*/

typedef struct {

  char                   name[32];            /* Variant name */

  char                   external_type[32];   /* External matrix type, or NULL */
  cs_matrix_type_t       type;                /* Matrix type */

  /* Function names, with variants:
     fill_type + exclude_diagonal_flag */

  char vector_multiply_name[CS_MATRIX_N_FILL_TYPES]
                           [CS_MATRIX_SPMV_N_TYPES][64];

  /* Measured structure creation cost, or -1 otherwise */

  double  matrix_create_cost;

  /* Measured assignment costs for each available fill type, or -1 otherwise */

  double  matrix_assign_cost[CS_MATRIX_N_FILL_TYPES];

  /* Measured operation costs for each available operation, or -1 otherwise
     fill_type + operation type + mean/variance + local/parallel */

  double     matrix_vector_cost[CS_MATRIX_N_FILL_TYPES]
                               [CS_MATRIX_SPMV_N_TYPES][2][2];
  cs_gnum_t  matrix_vector_n_ops[CS_MATRIX_N_FILL_TYPES]
                                [CS_MATRIX_SPMV_N_TYPES];

} cs_matrix_timing_variant_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

static const char *_matrix_fill_name[CS_MATRIX_N_FILL_TYPES]
  = {"scalar",
     "scalar symmetric",
     "block 3 diagonal",
     "block 6 diagonal",
     "block 3 diagonal symmetric",
     "block 3"};

static const char *_matrix_operation_name[CS_MATRIX_N_FILL_TYPES][2]
  = {{"y <- A.x",
      "y <- (A-D).x"},
     {"Symmetric y <- A.x",
      "Symmetric y <- (A-D).x"},
     {"Block 3 diagonal y <- A.x",
      "Block 3 diagonal y <- (A-D).x"},
     {"Block 6 diagonal y <- A.x",
      "Block 6 diagonal y <- (A-D).x"},
     {"Block 3 diagonal symmetric y <- A.x",
      "Block 3 diagonal symmetric y <- (A-D).x"},
     {"Block 3 y <- A.x",
      "Block 3 y <- (A-D).x"}};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize local variant matrix.
 *
 * parameters:
 *   v  <-> pointer to matrix variant
 *----------------------------------------------------------------------------*/

static void
_variant_init(cs_matrix_timing_variant_t  *v)
{
  v->matrix_create_cost = -1.;
  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    for (int j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
      v->vector_multiply_name[i][j][0] = '\0';
      for (int k = 0; k < 2; k++) {
        v->matrix_vector_cost[i][j][0][k] = -1.;
        v->matrix_vector_cost[i][j][1][k] = -1.;
      }
      v->matrix_vector_n_ops[i][j] = 0;
    v->matrix_assign_cost[i] = -1.;
    }
  }
}

/*----------------------------------------------------------------------------
 * Add variant
 *
 * parameters:
 *   name                 <-- matrix variant name
 *   external_type        <-- matrix external type name if applicable, or NULL
 *   type                 <-- matrix type
 *   n_fill_types         <-- number of fill types tuned for
 *   fill_types           <-- array of fill types tuned for
 *   op_flag              <-- bit mask for each operation type
 *   vector_multiply      <-- function pointer for A.x
 *   b_vector_multiply    <-- function pointer for block A.x
 *   bb_vector_multiply   <-- function pointer for block A.x
 *                             with block extra diag
 *   n_variants           <-> number of variants
 *   n_variants_max       <-> current maximum number of variants
 *   m_variant            <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_variant_add(const char                        *name,
             const char                        *external_type,
             cs_matrix_type_t                   type,
             int                                n_fill_types,
             cs_matrix_fill_type_t              fill_types[],
             int                                op_flag,
             const char                        *vector_multiply,
             const char                        *b_vector_multiply,
             const char                        *bb_vector_multiply,
             int                               *n_variants,
             int                               *n_variants_max,
             cs_matrix_timing_variant_t       **m_variant)
{
  cs_matrix_timing_variant_t  *v;
  int i = *n_variants;

  if (*n_variants_max == *n_variants) {
    if (*n_variants_max == 0)
      *n_variants_max = 8;
    else
      *n_variants_max *= 2;
    BFT_REALLOC(*m_variant, *n_variants_max, cs_matrix_timing_variant_t);
  }

  v = (*m_variant) + i;

  _variant_init(v);

  strcpy(v->name, name);
  memset(v->external_type, 0, 32);
  if (external_type != NULL)
    strncpy(v->external_type, external_type, 31);

  v->vector_multiply_name[0][0][0] = '\0';
  v->vector_multiply_name[0][1][0] = '\0';
  v->type = type;

  for (int j = 0; j < n_fill_types; j++) {

    cs_matrix_fill_type_t mft = fill_types[j];

    switch(mft) {

    case CS_MATRIX_SCALAR:
    case  CS_MATRIX_SCALAR_SYM:
      if (vector_multiply != NULL) {
        for (int k = 0; k < CS_MATRIX_SPMV_N_TYPES; k++) {
          if (op_flag & (1<<k))
            strncpy(v->vector_multiply_name[mft][k], vector_multiply, 63);
        }
      }
      break;

    case CS_MATRIX_BLOCK_D:
    case CS_MATRIX_BLOCK_D_66:
    case CS_MATRIX_BLOCK_D_SYM:
      if (b_vector_multiply != NULL) {
        for (int k = 0; k < CS_MATRIX_SPMV_N_TYPES; k++) {
          if (op_flag & (1<<k))
            strncpy(v->vector_multiply_name[mft][k], b_vector_multiply, 63);
        }
      }
      break;

    case CS_MATRIX_BLOCK:
      if (bb_vector_multiply != NULL) {
        for (int k = 0; k < CS_MATRIX_SPMV_N_TYPES; k++) {
          if (op_flag & (1<<k))
            strncpy(v->vector_multiply_name[mft][k], bb_vector_multiply, 63);
        }
      }
      break;

    default:
      assert(0);
      break;
    }

  }

  *n_variants += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build list of variants for testing.
 *
 * \param[in]   n_fill_types  number of fill types tuned for
 * \param[in]   fill_types    array of fill types tuned for
 * \param[in]   type_filter   true for matrix types tuned for, false for others
 * \param[in]   numbering     vectorization or thread-related numbering info,
 *                            or NULL
 * \param[out]  n_variants    number of variants
 * \param[out]  m_variant     array of matrix variants
 */
/*----------------------------------------------------------------------------*/

static void
_variant_build_list(int                             n_fill_types,
                    cs_matrix_fill_type_t           fill_types[],
                    bool                            type_filter[],
                    const cs_numbering_t           *numbering,
                    int                            *n_variants,
                    cs_matrix_timing_variant_t    **m_variant)
{
  int  n_variants_max = 0;

  int op_flag_a = 1<<CS_MATRIX_SPMV;
  int op_flag_e = 1<<CS_MATRIX_SPMV_E;
  int op_flag_ae = op_flag_a | op_flag_e;

  *n_variants = 0;
  *m_variant = NULL;

  if (type_filter[CS_MATRIX_NATIVE]) {

    _variant_add("Native, baseline",
                 NULL,
                 CS_MATRIX_NATIVE,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "baseline",
                 "baseline",
                 "baseline",
                 n_variants,
                 &n_variants_max,
                 m_variant);

    if (numbering != NULL) {

#if defined(HAVE_OPENMP)

      if (numbering->type == CS_NUMBERING_THREADS)
        _variant_add("Native, OpenMP",
                     NULL,
                     CS_MATRIX_NATIVE,
                     n_fill_types,
                     fill_types,
                     op_flag_ae,
                     "omp",
                     "omp",
                     NULL,
                     n_variants,
                     &n_variants_max,
                     m_variant);

      _variant_add("Native, OpenMP atomic",
                   NULL,
                   CS_MATRIX_NATIVE,
                   n_fill_types,
                   fill_types,
                   op_flag_ae,
                   "omp_atomic",
                   "omp_atomic",
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif

      if (numbering->type == CS_NUMBERING_VECTORIZE)
        _variant_add("Native, vectorized",
                     NULL,
                     CS_MATRIX_NATIVE,
                     n_fill_types,
                     fill_types,
                     op_flag_ae,
                     "vector",
                     NULL,
                     NULL,
                     n_variants,
                     &n_variants_max,
                     m_variant);

    }

  }

  if (type_filter[CS_MATRIX_CSR]) {

    _variant_add("CSR",
                 NULL,
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "default",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add("CSR, with MKL",
                 NULL,
                 CS_MATRIX_CSR,
                 n_fill_types,
                 fill_types,
                 op_flag_a,
                 "mkl",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

#if defined(HAVE_CUDA)

    if (cs_get_device_id() > -1)
      _variant_add("CSR, CUDA",
                   NULL,
                   CS_MATRIX_CSR,
                   n_fill_types,
                   fill_types,
                   op_flag_ae,
                   "cuda",
                   NULL,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUDA) */

#if defined(HAVE_CUSPARSE)

    if (cs_get_device_id() > -1)
      _variant_add("CSR, with cuSPARSE",
                   NULL,
                   CS_MATRIX_CSR,
                   n_fill_types,
                   fill_types,
                   op_flag_ae,
                   "cusparse",
                   NULL,
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUSPARSE) */

  }

  if (type_filter[CS_MATRIX_MSR]) {

    _variant_add("MSR",
                 NULL,
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "default",
                 "default",
                 "default",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add("MSR, with MKL",
                 NULL,
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "mkl",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

#if defined(HAVE_CUDA)

    if (cs_get_device_id() > -1)
      _variant_add("MSR, CUDA",
                   NULL,
                   CS_MATRIX_MSR,
                   n_fill_types,
                   fill_types,
                   op_flag_ae,
                   "cuda",
                   "cuda",
                   NULL,
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUDA) */

#if defined(HAVE_CUSPARSE)

    if (cs_get_device_id() > -1)
      _variant_add("MSR, with cuSPARSE",
                   NULL,
                   CS_MATRIX_MSR,
                   n_fill_types,
                   fill_types,
                   op_flag_ae,
                   "cusparse",
#if defined(HAVE_CUSPARSE_GENERIC_API)
                   "cusparse",
#else
                   NULL,
#endif
                   "cusparse",
                   n_variants,
                   &n_variants_max,
                   m_variant);

#endif /* defined(HAVE_CUSPARSE) */

    _variant_add("MSR, OpenMP scheduling",
                 NULL,
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "omp_sched",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_HYPRE)

    _variant_add("HYPRE (PARCSR)",
                 "HYPRE",
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 op_flag_a,
                 "external",
                 "external",
                 "external",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_HYPRE) */

#if defined(HAVE_PETSC)

    _variant_add("PETSc (MATAIJ)",
                 "PETSC",
                 CS_MATRIX_MSR,
                 n_fill_types,
                 fill_types,
                 op_flag_a,
                 "external",
                 "external",
                 "external",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

  }

  if (type_filter[CS_MATRIX_DIST]) {

    _variant_add("Distributed",
                 NULL,
                 CS_MATRIX_DIST,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "default",
                 "default",
                 "default",
                 n_variants,
                 &n_variants_max,
                 m_variant);

#if defined(HAVE_MKL)

    _variant_add("Distributed, with MKL",
                 NULL,
                 CS_MATRIX_DIST,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "mkl",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

#endif /* defined(HAVE_MKL) */

    _variant_add("Distributed, OpenMP scheduling",
                 NULL,
                 CS_MATRIX_DIST,
                 n_fill_types,
                 fill_types,
                 op_flag_ae,
                 "omp_sched",
                 NULL,
                 NULL,
                 n_variants,
                 &n_variants_max,
                 m_variant);

  }

  n_variants_max = *n_variants;
  BFT_REALLOC(*m_variant, *n_variants, cs_matrix_timing_variant_t);
}

/*----------------------------------------------------------------------------
 * Copy array to reference for matrix computation check.
 *
 * parameters:
 *   n_elts      <-- number values to compare
 *   y           <-- array to copare or copy
 *   yr          <-- reference array
 *
 * returns:
 *   maximum difference between values
 *----------------------------------------------------------------------------*/

static double
_matrix_check_compare(cs_lnum_t        n_elts,
                      const cs_real_t  y[],
                      cs_real_t        yr[])
{
  cs_lnum_t  ii;

  double dmax = 0.0;

  for (ii = 0; ii < n_elts; ii++) {
    double d = CS_ABS(y[ii] - yr[ii]);
    if (d > dmax)
      dmax = d;
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double dmaxg;
    MPI_Allreduce(&dmax, &dmaxg, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
    dmax = dmaxg;
  }

#endif

  return dmax;
}

/*----------------------------------------------------------------------------
 * Check local matrix.vector product operations.
 *
 * parameters:
 *   t_measure   <-- minimum time for each measure
 *   n_variants  <-- number of variants in array
 *   n_rows      <-- local number of rows
 *   n_cols_ext  <-- number of local + ghost columns
 *   n_edges     <-- local number of (undirected) graph edges
 *   edges       <-- edges (symmetric row <-> column) connectivity
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   m_variant   <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_matrix_check(int                          n_variants,
              cs_lnum_t                    n_rows,
              cs_lnum_t                    n_cols_ext,
              cs_lnum_t                    n_edges,
              const cs_lnum_2_t           *edges,
              const cs_halo_t             *halo,
              const cs_numbering_t        *numbering,
              cs_matrix_timing_variant_t  *m_variant)
{
  bool print_subtitle = false;
  cs_real_t  *da = NULL, *xa = NULL, *x = NULL, *y = NULL;
  cs_real_t  *yr0 = NULL;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;
  cs_lnum_t d_block_size = 3;
  cs_lnum_t e_block_size = 3;

  bft_printf
    ("\n"
     "Checking matrix structure and operation variants (diff/reference):\n"
     "------------------------------------------------\n");

  /* Allocate and initialize  working arrays
     (for a maximum block size of 6) */

  CS_MALLOC_HD(x, n_cols_ext*6, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(y, n_cols_ext*6, cs_real_t, cs_alloc_mode);
  BFT_MALLOC(yr0, n_cols_ext*6, cs_real_t);

  CS_MALLOC_HD(da, n_cols_ext*6*6, cs_real_t, cs_alloc_mode);
  BFT_MALLOC(xa, n_edges*2*6*6, cs_real_t);

  /* Initialize arrays */

  cs_gnum_t *cell_gnum = NULL;
  BFT_MALLOC(cell_gnum, n_cols_ext, cs_gnum_t);
  if (cs_glob_mesh->global_cell_num != NULL) {
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      cell_gnum[ii] = cs_glob_mesh->global_cell_num[ii];
  }
  else {
    for (cs_lnum_t ii = 0; ii < n_rows; ii++)
      cell_gnum[ii] = ii+1;
  }
  if (halo != NULL)
    cs_halo_sync_untyped(halo,
                         CS_HALO_STANDARD,
                         sizeof(cs_gnum_t),
                         cell_gnum);

  /* Loop on fill options */

  for (int f_id = 0; f_id < CS_MATRIX_N_FILL_TYPES; f_id++) {

    cs_lnum_t _d_block_size
      = (f_id >= CS_MATRIX_BLOCK_D) ? d_block_size : 1;
    if (f_id == CS_MATRIX_BLOCK_D_66)
      _d_block_size = 6;
    cs_lnum_t _e_block_size
      = (f_id >= CS_MATRIX_BLOCK) ? e_block_size : 1;
    const cs_lnum_t _block_mult = _d_block_size;
    const bool sym_coeffs = (   f_id == CS_MATRIX_SCALAR_SYM
                             || f_id == CS_MATRIX_BLOCK_D_SYM) ? true : false;

    /* Generate matrix coefficients */

    const cs_lnum_t stride = _d_block_size;
    const cs_lnum_t sd = _d_block_size*_d_block_size;
    const cs_lnum_t se = _e_block_size*_e_block_size;

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cols_ext; ii++) {
      cs_gnum_t jj = (cell_gnum[ii] - 1)*sd;
      for (cs_lnum_t kk = 0; kk < sd; kk++) {
        da[ii*sd+kk] = 1.0 + cos(jj*sd+kk);
      }
    }

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_edges; ii++) {
      cs_lnum_t i0 = edges[ii][0];
      cs_lnum_t i1 = edges[ii][1];
      cs_gnum_t j0 = (cell_gnum[i0] - 1)*se;
      cs_gnum_t j1 = (cell_gnum[i1] - 1)*se;
      for (cs_lnum_t kk = 0; kk < se; kk++) {
        xa[(ii*se+kk)*2]
          = 0.5*(0.45 + cos(j0*se+kk) + cos(j1*se+kk));
        xa[(ii*se+kk)*2 + 1]
          = -0.5*(0.45 + cos(j0*se+kk) + cos(j1*se+kk));
      }
    }

#   pragma omp parallel for
    for (cs_lnum_t ii = 0; ii < n_cols_ext; ii++) {
      cs_gnum_t jj = (cell_gnum[ii] - 1)*stride;
      for (cs_lnum_t kk = 0; kk < stride; kk++)
        x[ii*stride+kk] = sin(jj*stride+kk);
    }

    /* Loop on diagonal exclusion options */

    for (int op_id = 0; op_id < CS_MATRIX_SPMV_N_TYPES; op_id++) {

      print_subtitle = true;

      /* Loop on variant types */

      for (int v_id = 0; v_id < n_variants; v_id++) {

        cs_matrix_timing_variant_t *v = m_variant + v_id;

        int j = 0;
        for (j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
          if (strlen(v->vector_multiply_name[f_id][j]) >= 1)
            break;
        }
        if (j == CS_MATRIX_SPMV_N_TYPES)
          continue;

        ms = cs_matrix_structure_create(v->type,
                                        n_rows,
                                        n_cols_ext,
                                        n_edges,
                                        edges,
                                        halo,
                                        numbering);

        m = cs_matrix_create(ms);

#if defined(HAVE_HYPRE)
        if (strcmp(v->external_type, "HYPRE") == 0) {
          int device_id = cs_get_device_id();
          int use_device = (device_id < 0) ? 0 : 1;
          cs_matrix_set_type_hypre(m, use_device);
        }
#endif

#if defined(HAVE_PETSC)
        if (strcmp(v->external_type, "PETSc") == 0) {
          cs_matrix_set_type_petsc(m, NULL);
        }
#endif

        bool is_external_type = (strlen(v->external_type) == 0) ? false : true;

        cs_matrix_set_coefficients(m,
                                   sym_coeffs,
                                   _d_block_size,
                                   _e_block_size,
                                   n_edges,
                                   edges,
                                   da,
                                   xa);

        cs_matrix_vector_product_t  *vector_multiply = NULL;

        if (strlen(v->vector_multiply_name[f_id][op_id]) > 0) {

          if (is_external_type == false) {
            cs_matrix_variant_t *mv = cs_matrix_variant_create(m);
            cs_matrix_variant_set_func(mv,
                                       f_id,
                                       op_id,
                                       numbering,
                                       v->vector_multiply_name[f_id][op_id]);
            cs_matrix_variant_apply(m, mv);
            cs_matrix_variant_destroy(&mv);
          }

          vector_multiply = m->vector_multiply[f_id][op_id];

        }

        if (vector_multiply != NULL) {
          /* Set part of y to incorrect value (to check setting) */

          for (cs_lnum_t ii = 0; ii < n_cols_ext*_d_block_size; ii++)
            y[ii] = -1;

          /* Check multiplication */

          vector_multiply(m, op_id, true, x, y);
          if (v_id == 0)
            memcpy(yr0, y, n_rows*_block_mult*sizeof(cs_real_t));
          else {
            double dmax = _matrix_check_compare(n_rows*_block_mult, y, yr0);
            if (print_subtitle) {
              bft_printf("\n%s\n",
                         _matrix_operation_name[f_id][op_id]);
              print_subtitle = false;
            }
            bft_printf("  %-32s : %12.5e\n",
                       v->name,
                       dmax);
            bft_printf_flush();
          }

        }

        cs_matrix_release_coefficients(m);

        cs_matrix_destroy(&m);
        cs_matrix_structure_destroy(&ms);

      } /* end of loop on variants */

    } /* end of loop on op_id */

  } /* end of loop on fill types */

  BFT_FREE(cell_gnum);

  BFT_FREE(yr0);

  BFT_FREE(y);
  BFT_FREE(x);

  BFT_FREE(xa);
  BFT_FREE(da);
}

/*----------------------------------------------------------------------------
 * Tune local matrix.vector product operations.
 *
 * parameters:
 *   t_measure   <-- minimum time for each measure
 *   n_variants  <-- number of variants in array
 *   n_cells     <-- number of local cells
 *   n_cells_ext <-- number of cells including ghost cells (array size)
 *   n_faces     <-- local number of internal faces
 *   cell_num    <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell   <-- face -> cells connectivity
 *   halo        <-- cell halo structure
 *   numbering   <-- vectorization or thread-related numbering info, or NULL
 *   m_variant   <-> array of matrixtiming  variants
 *----------------------------------------------------------------------------*/

static void
_matrix_time_test(double                       t_measure,
                  int                          n_variants,
                  cs_lnum_t                    n_cells,
                  cs_lnum_t                    n_cells_ext,
                  cs_lnum_t                    n_faces,
                  const cs_lnum_2_t           *face_cell,
                  const cs_halo_t             *halo,
                  const cs_numbering_t        *numbering,
                  cs_matrix_timing_variant_t  *m_variant)
{
  cs_lnum_t  ii;
  int  n_runs, run_id, v_id, f_id;
  double  wt0, wt1, wtu;
  double wti, wtf;
  cs_matrix_type_t  type, type_prev;

  double test_sum = 0.0;
  cs_real_t  *da = NULL, *xa = NULL, *x = NULL, *y = NULL;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;
  cs_lnum_t d_block_size = 3, e_block_size = 3;
  cs_lnum_t d_block_stride = d_block_size*d_block_size;
  cs_lnum_t e_block_stride = e_block_size*e_block_size;

  type_prev = CS_MATRIX_N_TYPES;

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  CS_MALLOC_HD(x, n_cells_ext*d_block_size, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(y, n_cells_ext*d_block_size, cs_real_t, cs_alloc_mode);

  CS_MALLOC_HD(da, n_cells_ext*d_block_stride, cs_real_t, cs_alloc_mode);
  CS_MALLOC_HD(xa, n_faces*e_block_stride*2, cs_real_t, cs_alloc_mode);

# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_stride; ii++)
    da[ii] = 1.0;
# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_size; ii++)
    x[ii] = ii*0.1/n_cells_ext;

# pragma omp parallel for
  for (ii = 0; ii < n_faces*e_block_stride; ii++) {
    xa[ii*2] = 0.5;
    xa[ii*2 + 1] = -0.5;
  }

  /* Loop on variant types */
  /*-----------------------*/

  for (v_id = 0; v_id < n_variants; v_id++) {

    bool test_assign = false;

    cs_matrix_timing_variant_t *v = m_variant + v_id;

    type = v->type;

    if (type != type_prev) {

      test_assign = true;

      wt0 = cs_timer_wtime(), wt1 = wt0;
      run_id = 0, n_runs = (t_measure > 0) ? 16 : 1;
      while (run_id < n_runs) {
        while (run_id < n_runs) {
          if (m != NULL)
            cs_matrix_destroy(&m);
          if (ms != NULL)
            cs_matrix_structure_destroy(&ms);
          ms = cs_matrix_structure_create(type,
                                          n_cells,
                                          n_cells_ext,
                                          n_faces,
                                          face_cell,
                                          halo,
                                          numbering);
          m = cs_matrix_create(ms);
          run_id++;
        }
        wt1 = cs_timer_wtime();
        double wt_r0 = wt1 - wt0;
        cs_parall_max(1, CS_DOUBLE, &wt_r0);
        if (wt_r0 < t_measure)
          n_runs *= 2;
      }
      v->matrix_create_cost = (wt1 - wt0) / n_runs;
    }

    /* Loop on fill patterns sizes */

    for (f_id = 0; f_id < CS_MATRIX_N_FILL_TYPES; f_id++) {

      const cs_lnum_t _d_block_size
        = (f_id >= CS_MATRIX_BLOCK_D) ? d_block_size : 1;
      const cs_lnum_t _e_block_size
        = (f_id >= CS_MATRIX_BLOCK) ? e_block_size : 1;
      const bool sym_coeffs
        = (   f_id == CS_MATRIX_SCALAR_SYM
           || f_id == CS_MATRIX_BLOCK_D_SYM) ? true : false;

      /* Loop on diagonal exclusion flags */

      double t_measure_assign = -1;

      if (   strlen(v->vector_multiply_name[f_id][0]) < 1
          && strlen(v->vector_multiply_name[f_id][1]) < 1)
        continue;

      if (m != NULL)
        cs_matrix_destroy(&m);

      m = cs_matrix_create(ms);

      /* Measure overhead of setting coefficients if not already done */

      if (test_assign) {
        t_measure_assign = t_measure;
        n_runs = 16;
      }
      else
        n_runs = 1;

      wt0 = cs_timer_wtime(), wt1 = wt0;
      run_id = 0;
      while (run_id < n_runs) {
        while (run_id < n_runs) {
          cs_matrix_set_coefficients(m,
                                     sym_coeffs,
                                     _d_block_size,
                                     _e_block_size,
                                     n_faces,
                                     (const cs_lnum_2_t *)face_cell,
                                     da,
                                     xa);
          run_id++;
        }
        wt1 = cs_timer_wtime();
        double wt_r0 = wt1 - wt0;
        cs_parall_max(1, CS_DOUBLE, &wt_r0);
        if (wt_r0 < t_measure_assign)
          n_runs *= 2;
      }
      if (n_runs > 1)
        v->matrix_assign_cost[f_id] = (wt1 - wt0) / n_runs;

      /* Measure matrix.vector operations */

      bool is_external_type = (strlen(v->external_type) == 0) ? false : true;

      for (int op_id = 0; op_id < CS_MATRIX_SPMV_N_TYPES; op_id++) {

        cs_matrix_vector_product_t  *vector_multiply = NULL;

        if (is_external_type == false) {

          if (   strlen(v->external_type) == 0
              && strlen(v->vector_multiply_name[f_id][op_id]) > 0) {
            cs_matrix_variant_t *mv = cs_matrix_variant_create(m);
            cs_matrix_variant_set_func(mv,
                                       f_id,
                                       op_id,
                                       numbering,
                                       v->vector_multiply_name[f_id][op_id]);
            cs_matrix_variant_apply(m, mv);
            cs_matrix_variant_destroy(&mv);
          }

        }

        vector_multiply = m->vector_multiply[f_id][op_id];

        if (vector_multiply == NULL)
          continue;

        cs_lnum_t n_rows = cs_matrix_get_n_rows(m);
        cs_lnum_t nnz = cs_matrix_get_n_entries(m);
        const cs_lnum_t db_size = cs_matrix_get_diag_block_size(m);
        const cs_lnum_t eb_size = cs_matrix_get_extra_diag_block_size(m);

        v->matrix_vector_n_ops[f_id][op_id] = n_rows*db_size*db_size;
        if (eb_size == 1)
          v->matrix_vector_n_ops[f_id][op_id] += (nnz-n_rows)*db_size;
        else
          v->matrix_vector_n_ops[f_id][op_id] +=  (nnz-n_rows)
                                                 *(eb_size*eb_size);

        int mpi_flag_max = (cs_glob_n_ranks > 1) ? 2 : 1;

        for (int mpi_flag = 0; mpi_flag < mpi_flag_max; mpi_flag++) {

          wt0 = cs_timer_wtime(), wt1 = wt0;
          run_id = 0, n_runs = (t_measure > 0) ? 16 : 1;

          double mean = 0.0;
          double m2 = 0.0;

          while (run_id < n_runs) {
            while (run_id < n_runs) {
              if (run_id % 16)
                test_sum = 0;
              wti = cs_timer_wtime(), wtf = wti;
              if (mpi_flag > 0)
                vector_multiply(m, op_id, false, x, y);
              else {
                if (op_id == 0)
                  cs_matrix_vector_multiply(m, x, y);
                else
                  cs_matrix_vector_multiply_partial(m, op_id, x, y);
              }
              wtf = cs_timer_wtime();
              test_sum += y[n_cells-1];
              run_id++;
              double t = (wtf - wti);
              double delta = t - mean;
              double r = delta / run_id;
              double m_n = mean + r;
              m2 = (m2*(run_id-1) + delta*(t-m_n)) / run_id;
              mean += r;
            }
            wt1 = cs_timer_wtime();
            double wt_r0 = wt1 - wt0;
            cs_parall_max(1, CS_DOUBLE, &wt_r0);
            if (wt_r0 < t_measure)
              n_runs *= 2;
          }

          wtu = (wt1 - wt0) / n_runs;
          if (n_runs > 1)
            m2 = sqrt(m2 / (n_runs - 1));
          else
            m2 = 0;
          v->matrix_vector_cost[f_id][op_id][0][mpi_flag] = wtu;
          v->matrix_vector_cost[f_id][op_id][1][mpi_flag] = m2;

        } /* End of loop on mpi_flag */

      } /* end of loop on op_id */

      cs_matrix_release_coefficients(m);

    } /* end of loop on fill patterns */

    type_prev = type;

  } /* end of loop on variants */

  if (m != NULL)
    cs_matrix_destroy(&m);
  if (ms != NULL)
    cs_matrix_structure_destroy(&ms);

  BFT_FREE(x);
  BFT_FREE(y);

  BFT_FREE(da);
  BFT_FREE(xa);
}

/*----------------------------------------------------------------------------
 * Print title for statistics on matrix timing SpMV info.
 *
 * parameters:
 *   struct_flag <-- 0: assignment; 1: structure creation
 *   fill_type   <-- matrix fill type
 *----------------------------------------------------------------------------*/

static void
_matrix_time_create_assign_title(int                    struct_flag,
                                 cs_matrix_fill_type_t  fill_type)
{
  size_t i = 0;
  size_t l = 80;
  char title[81] = "";

  /* Print title */

  if (struct_flag == 0) {
    snprintf(title + i,  l-i, " matrix %s coefficients assign",
             _matrix_fill_name[fill_type]);
    title[80] = '\0';
    i = strlen(title);
    l -= i;
  }
  else
    strncat(title + i, "matrix structure creation/destruction", l);

  title[80] = '\0';

  l = cs_log_strlen(title);

  cs_log_printf(CS_LOG_PERFORMANCE, "\n%s\n", title);

  for (i = 0; i < l; i++)
    title[i] = '-';
  title[l] = '\0';

  cs_log_printf(CS_LOG_PERFORMANCE, "%s\n", title);

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    char tmp_s[4][24] =  {"", "", "", ""};

    cs_log_strpadl(tmp_s[0], "time (s)", 10, 24);
    cs_log_strpadl(tmp_s[1], "mean", 10, 24);
    cs_log_strpadl(tmp_s[2], "max", 10, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %10s %s\n"
                  "  %24s %s %s\n",
                  " ", " ", tmp_s[0],
                  " ", tmp_s[1], tmp_s[2]);
  }

#endif

  if (cs_glob_n_ranks == 1) {

    char tmp_s[24] =  {""};

    cs_log_strpadl(tmp_s, "time (s)", 10, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %s\n",
                  " ", tmp_s);

  }
}

/*----------------------------------------------------------------------------
 * Print statistics on matrix timing creation or assignment info.
 *
 * parameters:
 *   m_variant   <-- array of matrix variants
 *   variant_id  <-- variant id
 *   struct_flag <-- 0: assignment; 1: structure creation
 *   fill_type   <-- type of matrix fill
 *----------------------------------------------------------------------------*/

static void
_matrix_time_create_assign_stats(const cs_matrix_timing_variant_t  *m_variant,
                                 int                                variant_id,
                                 int                                struct_flag,
                                 cs_matrix_fill_type_t              fill_type)
{
  char title[32];

  double t_loc = -1;

  const cs_matrix_timing_variant_t  *v = m_variant + variant_id;

  cs_log_strpad(title, v->name, 24, 32);

  if (struct_flag == 0)
    t_loc = v->matrix_assign_cost[fill_type];
  else
    t_loc = v->matrix_create_cost;

  if (t_loc < 0)
    return;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {
    double t_max, t_sum = -1;
    MPI_Allreduce(&t_loc, &t_sum, 1, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&t_loc, &t_max, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %10.3e %10.3e\n",
                  title, t_sum/cs_glob_n_ranks, t_max);
  }

#endif

  if (cs_glob_n_ranks == 1)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %10.3e\n", title, t_loc);
}

/*----------------------------------------------------------------------------
 * Print title for statistics on matrix timing SpMV info.
 *
 * parameters:
 *   fill_type   <-- type of matrix fill
 *   op_type     <-- operation type
 *   mpi_flag    <-- 0: include MPI; 1: local only
 *   mpi_flag    <-- 0: with MPI; 1: local only
 *----------------------------------------------------------------------------*/

static void
_matrix_time_spmv_title(cs_matrix_fill_type_t  fill_type,
                        cs_matrix_spmv_type_t  op_type,
                        int                    mpi_flag)
{
  size_t i = 0;
  size_t l = 80;
  char title[81] = "";

  /* Print title */

  snprintf(title, 80, "%s",
           _matrix_operation_name[fill_type][op_type]);
  title[80] = '\0';
  l = cs_log_strlen(title);

  if (mpi_flag == 0)
    cs_log_printf(CS_LOG_PERFORMANCE, "\n%s\n", title);
  else
    cs_log_printf(CS_LOG_PERFORMANCE, "\n%s (local)\n", title);

  for (i = 0; i < l; i++)
    title[i] = '-';
  title[l] = '\0';

  cs_log_printf(CS_LOG_PERFORMANCE, "%s\n", title);

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    char tmp_s[12][24] =  {"", "", "", "", "", "", "", "", "", "", "", ""};

    cs_log_strpadl(tmp_s[0], "time (s)", 12, 24);
    cs_log_strpadl(tmp_s[1], "speedup", 12, 24);
    cs_log_strpadl(tmp_s[2], "std. dev.", 12, 24);
    cs_log_strpadl(tmp_s[3], "mean", 10, 24);
    cs_log_strpadl(tmp_s[4], "max", 10, 24);
    cs_log_strpadl(tmp_s[5], "mean", 6, 24);
    cs_log_strpadl(tmp_s[6], "min", 6, 24);
    cs_log_strpadl(tmp_s[7], "max", 6, 24);
    cs_log_strpadl(tmp_s[8], "mean", 9, 24);
    cs_log_strpadl(tmp_s[9], "max", 9, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %8s %s   %7s %s   %6s %s\n"
                  "  %24s %s %s | %s %s %s | %s %s\n",
                  " ", " ", tmp_s[0], " ", tmp_s[1], " ", tmp_s[2],
                  " ", tmp_s[3], tmp_s[4], tmp_s[5],
                  tmp_s[6], tmp_s[7], tmp_s[8], tmp_s[9]);
  }

#endif

  if (cs_glob_n_ranks == 1) {

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s  time (s)  speedup  std. dev.\n", " ");

  }
}

/*----------------------------------------------------------------------------
 * Print statistics on matrix timing SpMV info.
 *
 * parameters:
 *   m_variant   <-- array of matrix variants
 *   variant_id  <-- variant id
 *   fill_type   <-- type of matrix fill
 *   op_type     <-- SpMV operation type.
 *   mpi_flag    <-- 0: with MPI; 1: local only
 *----------------------------------------------------------------------------*/

static void
_matrix_time_spmv_stats(const cs_matrix_timing_variant_t  *m_variant,
                        int                                variant_id,
                        cs_matrix_fill_type_t              fill_type,
                        cs_matrix_spmv_type_t              op_type,
                        int                                mpi_flag)
{
  char title[32];

  double v_loc[3] = {-1, -1, 0};

  const cs_matrix_timing_variant_t  *r = m_variant;
  const cs_matrix_timing_variant_t  *v = m_variant + variant_id;

  cs_log_strpad(title, v->name, 24, 32);

  /* Get timing info */

  v_loc[0] = v->matrix_vector_cost[fill_type][op_type][0][mpi_flag];
  v_loc[1] = r->matrix_vector_cost[fill_type][op_type][0][mpi_flag];
  v_loc[2] = v->matrix_vector_cost[fill_type][op_type][1][mpi_flag];

  if (v_loc[0] < 0)
    return;

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double v_max[3], speedup_min, v_sum[3];

    MPI_Allreduce(v_loc, v_sum, 3, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);

    if (v_loc[0] > 1e-12) /* Always above 10-9 with real hardware */
      v_loc[1] =   r->matrix_vector_cost[fill_type][op_type][0][mpi_flag]
                 / v_loc[0];
    else
      v_loc[1] = 1;
    MPI_Allreduce(v_loc + 1, &speedup_min, 1, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);

    MPI_Allreduce(v_loc, v_max, 3, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    cs_real_t t_mean = v_sum[0]/cs_glob_n_ranks;
    cs_real_t speedup_mean = v_sum[1] / v_sum[0]; /* weighted by v_loc[0] */
    cs_real_t stddev_mean = v_sum[2]/cs_glob_n_ranks;
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %10.3e %10.3e | %6.3f %6.3f %6.3f | %9.2e %9.2e\n",
                  title,
                  t_mean, v_max[0], speedup_mean, speedup_min, v_max[1],
                  stddev_mean, v_max[2]);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    cs_real_t speedup = v_loc[1] / v_loc[0];
    cs_real_t stddev = v_loc[2] / v_loc[0];
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %10.3e  %6.3f  %9.2e\n",
                  title,
                  v_loc[0], speedup, stddev);
  }
}

/*----------------------------------------------------------------------------
 * Print title for FLOPS statistics on matrix timing SpMv performance info.
 *----------------------------------------------------------------------------*/

static void
_matrix_time_spmv_title_ops(void)
{
  /* Skip line */

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    char tmp_s[12][24] =  {"", "", "", "", "", "", "", "", "", "", "", ""};

    cs_log_strpadl(tmp_s[0], "GFLOPS", 12, 24);
    cs_log_strpadl(tmp_s[1], "mean", 10, 24);
    cs_log_strpadl(tmp_s[2], "min", 10, 24);
    cs_log_strpadl(tmp_s[3], "max", 10, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %19s %s\n"
                  "  %24s %s %s %s\n",
                  " ", " ", tmp_s[0],
                  " ", tmp_s[1], tmp_s[2], tmp_s[3]);
  }

#endif

  if (cs_glob_n_ranks == 1) {

    char tmp_s[1][24] =  {""};

    cs_log_strpadl(tmp_s[0], "GFLOPS", 10, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %s\n",
                  " ", tmp_s[0]);

  }
}

/*----------------------------------------------------------------------------
 * Print statistics on matrix timing SpMv performance info.
 *
 * parameters:
 *   m_variant   <-- array of matrix variants
 *   variant_id  <-- variant id
 *   fill_type   <-- type of matrix fill
 *   op_type     <-- SpMV operation type
 *   mpi_flag    <-- 0: with MPI; 1: local only
 *----------------------------------------------------------------------------*/

static void
_matrix_time_spmv_stats_ops(const cs_matrix_timing_variant_t  *m_variant,
                            int                                variant_id,
                            cs_matrix_fill_type_t              fill_type,
                            cs_matrix_spmv_type_t              op_type,
                            int                                mpi_flag)
{
  char title[32];

  const cs_matrix_timing_variant_t  *v = m_variant + variant_id;

  cs_log_strpad(title, v->name, 24, 32);

  /* Get timing info */

  double f_loc = 0;

  if (v->matrix_vector_cost[fill_type][op_type][0][mpi_flag] > 1e-12)
    f_loc =   v->matrix_vector_n_ops[fill_type][op_type]
            / v->matrix_vector_cost[fill_type][op_type][0][mpi_flag];

  f_loc /= 1e9;

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double f_min, f_max, f_mean;

    MPI_Allreduce(&f_loc, &f_mean, 1, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);
    MPI_Allreduce(&f_loc, &f_min, 1, MPI_DOUBLE, MPI_MIN, cs_glob_mpi_comm);
    MPI_Allreduce(&f_loc, &f_max, 1, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    f_mean /= cs_glob_n_ranks;

    if (f_min > 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    "  %s %10.3e %10.3e %10.3e\n",
                    title,
                    f_mean, f_min, f_max);
  }

#endif

  if (cs_glob_n_ranks == 1 && f_loc > 0) {
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %10.3e\n",
                  title, f_loc);
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Time matrix operations.
 *
 * parameters:
 *   t_measure      <-- minimum time for each measure
 *   n_types        <-- number of matrix types timed, or 0
 *   n_fill_types   <-- number of fill types timed, or 0
 *   types          <-- array of matrix types timed, or NULL
 *   fill_types     <-- array of fill types timed, or NULL
 *   n_cells        <-- number of local cells
 *   n_cells_ext    <-- number of cells including ghost cells (array size)
 *   n_faces        <-- local number of internal faces
 *   cell_num       <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell      <-- face -> cells connectivity
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *----------------------------------------------------------------------------*/

void
cs_benchmark_matrix(double                 t_measure,
                    int                    n_types,
                    int                    n_fill_types,
                    cs_matrix_type_t       types[],
                    cs_matrix_fill_type_t  fill_types[],
                    cs_lnum_t              n_cells,
                    cs_lnum_t              n_cells_ext,
                    cs_lnum_t              n_faces,
                    const cs_lnum_2_t     *face_cell,
                    const cs_halo_t       *halo,
                    const cs_numbering_t  *numbering)
{
  int  t_id, f_id, v_id, ed_flag;

  bool                   type_filter[CS_MATRIX_N_BUILTIN_TYPES] = {true,
                                                                   true,
                                                                   true,
                                                                   true};

  int                    _n_fill_types_default = 4;
  cs_matrix_fill_type_t  _fill_types_default[] = {CS_MATRIX_SCALAR,
                                                  CS_MATRIX_SCALAR_SYM,
                                                  CS_MATRIX_BLOCK_D,
                                                  CS_MATRIX_BLOCK};

  int                    _n_types = n_types;
  int                    _n_fill_types = n_fill_types;
  cs_matrix_fill_type_t  *_fill_types = fill_types;

  int  n_variants = 0;
  cs_matrix_timing_variant_t  *m_variant = NULL;

  cs_timer_t  t0, t1;

  t0 = cs_timer_time();

  /* Use defaults if required */

  if (_n_types > 0) {
    for (t_id = 0; t_id < CS_MATRIX_N_BUILTIN_TYPES; t_id++)
      type_filter[t_id] = false;
    for (t_id = 0; t_id < n_types; t_id++)
      type_filter[types[t_id]] = true;
  }

  if (_n_fill_types < 1) {
    _n_fill_types =  _n_fill_types_default;
    _fill_types = _fill_types_default;
  }

  /* Build variants array */
  /*----------------------*/

  _variant_build_list(_n_fill_types,
                      _fill_types,
                      type_filter,
                      numbering,
                      &n_variants,
                      &m_variant);

  /* Check results consistency */

  _matrix_check(n_variants,
                n_cells,
                n_cells_ext,
                n_faces,
                face_cell,
                halo,
                numbering,
                m_variant);

  /* Run tests on variants */

  _matrix_time_test(t_measure,
                    n_variants,
                    n_cells,
                    n_cells_ext,
                    n_faces,
                    face_cell,
                    halo,
                    numbering,
                    m_variant);

  /* Print info on variants */

  _matrix_time_create_assign_title(1, 0);
  for (v_id = 0; v_id < n_variants; v_id++)
    _matrix_time_create_assign_stats(m_variant, v_id, 1, CS_MATRIX_SCALAR);

  for (f_id = 0; f_id < _n_fill_types; f_id++) {
    cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
    _matrix_time_create_assign_title(0, fill_type);
    for (v_id = 0; v_id < n_variants; v_id++)
      _matrix_time_create_assign_stats(m_variant,
                                       v_id,
                                       0,
                                       fill_type);
  }

  int mpi_flag_max = (cs_glob_n_ranks > 1) ? 2 : 1;

  for (int mpi_flag = 0; mpi_flag < mpi_flag_max; mpi_flag++) {

    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
      for (ed_flag = 0; ed_flag < 2; ed_flag++) {
        _matrix_time_spmv_title(fill_type, ed_flag, mpi_flag);
        for (v_id = 0; v_id < n_variants; v_id++)
          _matrix_time_spmv_stats(m_variant,
                                  v_id,
                                  fill_type,
                                  ed_flag,
                                  mpi_flag);
        _matrix_time_spmv_title_ops();
        for (v_id = 0; v_id < n_variants; v_id++)
          _matrix_time_spmv_stats_ops(m_variant,
                                      v_id,
                                      fill_type,
                                      ed_flag,
                                      mpi_flag);
      }
    }

  }

  BFT_FREE(m_variant);

  t1 = cs_timer_time();

  cs_timer_counter_t  t_tune = cs_timer_diff(&t0, &t1);

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Elapsed timing time: %12.3f s.\n",
                t_tune.nsec*1e-9);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
