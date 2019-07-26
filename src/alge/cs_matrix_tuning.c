/*============================================================================
 * Sparse Matrix Representation and Operations.
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

/*
 * Notes:
 *
 * The aim of these structures and associated functions is multiple:
 *
 * - Provide an "opaque" matrix object for linear solvers, allowing possible
 *   choice of the matrix type based on run-time tuning at code initialization
 *   (depending on matrix size, architecture, and compiler, the most efficient
 *   structure for matrix.vector products may vary).
 *
 * - Provide at least a CSR matrix structure in addition to the "native"
 *   matrix structure, as this may allow us to leverage existing librairies.
 *
 * - Provide a C interface, also so as to be able to interface more easily
 *   with external libraries.
 *
 * The structures used here could easily be extended to block matrixes,
 * using for example the same structure information with 3x3 blocks which
 * could arise from coupled velocity components. This would imply that the
 * corresponding vectors be interlaced (or an interlaced copy be used
 * for recurring operations such as sparse linear system resolution),
 * for better memory locality, and possible loop unrolling.
 */

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
#include "cs_prototypes.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_matrix.h"
#include "cs_matrix_priv.h"

#include "cs_matrix_tuning.h"

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

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for matrix types */

static const char *_matrix_operation_name[CS_MATRIX_N_FILL_TYPES][2]
  = {{N_("y <- A.x"),
      N_("y <- (A-D).x")},
     {N_("Symmetric y <- A.x"),
      N_("Symmetric y <- (A-D).x")},
     {N_("Block 3 diagonal y <- A.x"),
      N_("Block 3 diagonal y <- (A-D).x")},
     {N_("Block 6 diagonal y <- A.x"),
      N_("Block 6 diagonal y <- (A-D).x")},
     {N_("Block 3 diagonal symmetric y <- A.x"),
      N_("Block 3 diagonal symmetric y <- (A-D).x")},
     {N_("Block 3 y <- A.x"),
      N_("Block 3 y <- (A-D).x")}};

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
 *   m_variant   <-> array of matrix variants
 *----------------------------------------------------------------------------*/

static void
_matrix_tune_test(double                 t_measure,
                  int                    n_variants,
                  cs_lnum_t              n_cells,
                  cs_lnum_t              n_cells_ext,
                  cs_lnum_t              n_faces,
                  const cs_lnum_2_t     *face_cell,
                  const cs_halo_t       *halo,
                  const cs_numbering_t  *numbering,
                  cs_matrix_variant_t   *m_variant)
{
  cs_lnum_t  ii;
  int  n_runs, run_id, v_id, f_id, ed_flag;
  double  wt0, wt1, wtu;
  double wti, wtf;
  cs_matrix_type_t  type, type_prev;

  double test_sum = 0.0;
  cs_real_t  *da = NULL, *xa = NULL, *x = NULL, *y = NULL;
  cs_matrix_structure_t *ms = NULL;
  cs_matrix_t *m = NULL;
  int d_block_size[4] = {3, 3, 3, 9};
  int ed_block_size[4] = {3, 3, 3, 9};

  type_prev = CS_MATRIX_N_TYPES;

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  if (CS_MEM_ALIGN > 0) {
    BFT_MEMALIGN(x, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n_cells_ext*d_block_size[1], cs_real_t);
  }
  else {
    BFT_MALLOC(x, n_cells_ext*d_block_size[1], cs_real_t);
    BFT_MALLOC(y, n_cells_ext*d_block_size[1], cs_real_t);
  }

  BFT_MALLOC(da, n_cells_ext*d_block_size[3], cs_real_t);
  BFT_MALLOC(xa, n_faces*ed_block_size[3]*2, cs_real_t);

# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_size[3]; ii++)
    da[ii] = 1.0;
# pragma omp parallel for
  for (ii = 0; ii < n_cells_ext*d_block_size[1]; ii++)
    x[ii] = ii*0.1/n_cells_ext;

# pragma omp parallel for
  for (ii = 0; ii < n_faces*ed_block_size[3]; ii++) {
    xa[ii*2] = 0.5;
    xa[ii*2 + 1] = -0.5;
  }

  /* Loop on variant types */
  /*-----------------------*/

  for (v_id = 0; v_id < n_variants; v_id++) {

    cs_matrix_variant_t *v = m_variant + v_id;

    type = v->type;

    if (type != type_prev) {

      if (m != NULL)
        cs_matrix_destroy(&m);
      if (ms != NULL)
        cs_matrix_structure_destroy(&ms);
      ms = cs_matrix_structure_create(type,
                                      true,
                                      n_cells,
                                      n_cells_ext,
                                      n_faces,
                                      face_cell,
                                      halo,
                                      numbering);
      m = cs_matrix_create(ms);

    }

    /* Loop on fill patterns sizes */

    for (f_id = 0; f_id < CS_MATRIX_N_FILL_TYPES; f_id++) {

      const int *_d_block_size
        = (f_id >= CS_MATRIX_BLOCK_D) ? d_block_size : NULL;
      const int *_ed_block_size
        = (f_id >= CS_MATRIX_BLOCK) ? ed_block_size : NULL;
      const bool sym_coeffs
        = (   f_id == CS_MATRIX_SCALAR_SYM
           || f_id == CS_MATRIX_BLOCK_D_SYM) ? true : false;

      /* Loop on diagonal exclusion flags */

      if (   v->vector_multiply[f_id][0] == NULL
          && v->vector_multiply[f_id][1] == NULL)
        continue;

      /* Measure overhead of setting coefficients if not already done */

      cs_matrix_set_coefficients(m,
                                 sym_coeffs,
                                 _d_block_size,
                                 _ed_block_size,
                                 n_faces,
                                 (const cs_lnum_2_t *)face_cell,
                                 da,
                                 xa);

      /* Measure matrix.vector operations */

      for (ed_flag = 0; ed_flag < 2; ed_flag++) {

        cs_matrix_vector_product_t
          *vector_multiply = v->vector_multiply[f_id][ed_flag];

        if (vector_multiply == NULL)
          continue;

        if (vector_multiply != NULL) {
          wt0 = cs_timer_wtime(), wt1 = wt0;
          run_id = 0, n_runs = 8;

          double mean = 0.0;
          double m2 = 0.0;
          double delta = 0.0;

          while (run_id < n_runs) {
            while (run_id < n_runs) {
              if (run_id % 8)
                test_sum = 0;
              wti = cs_timer_wtime(), wtf = wti;
              vector_multiply(ed_flag, m, x, y);
              wtf = cs_timer_wtime();
              test_sum += y[n_cells-1];
              run_id++;
              delta = (wtf - wti) - mean;
              mean = mean + delta/run_id;
              m2 = m2 + delta * (wtf - wti - mean);
            }
            wt1 = cs_timer_wtime();
            if (wt1 - wt0 < t_measure)
              n_runs *= 2;
          }
          wtu = (wt1 - wt0) / n_runs;
          m2 = sqrt(m2 / (n_runs - 1));
          v->matrix_vector_cost[f_id][ed_flag][0] = wtu;
          v->matrix_vector_cost[f_id][ed_flag][1] = m2;
        }
      } /* end of loop on ed_flag */

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
 * Print title for statistics on matrix tuning SpMv info.
 *
 * parameters:
 *   fill_type   <-- type of matrix fill
 *   ed_flag     <-- 0: include diagonal; 1: exclude diagonal
 *----------------------------------------------------------------------------*/

static void
_matrix_tune_spmv_title(cs_matrix_fill_type_t  fill_type,
                        int                    ed_flag)
{
  size_t i = 0;
  size_t l = 80;
  char title[81] = "";

  /* Print title */

  snprintf(title, 80, "%s",
           _(_matrix_operation_name[fill_type][ed_flag]));
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

    char tmp_s[12][24] =  {"", "", "", "", "", "", "", "", "", "", "", ""};

    cs_log_strpadl(tmp_s[0], _("time (s)"), 16, 24);
    cs_log_strpadl(tmp_s[1], _("speedup"), 16, 24);
    cs_log_strpadl(tmp_s[2], _("std. dev."), 16, 24);
    cs_log_strpadl(tmp_s[3], _("mean"), 12, 24);
    cs_log_strpadl(tmp_s[4], _("max"), 12, 24);
    cs_log_strpadl(tmp_s[5], _("mean"), 8, 24);
    cs_log_strpadl(tmp_s[6], _("min"), 8, 24);
    cs_log_strpadl(tmp_s[7], _("max"), 8, 24);
    cs_log_strpadl(tmp_s[8], _("mean"), 8, 24);
    cs_log_strpadl(tmp_s[9], _("max"), 8, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %8s %s %9s %s  %s\n"
                  "  %24s %s %s %s %s %s %s %s\n",
                  " ", " ", tmp_s[0], " ", tmp_s[1], tmp_s[2],
                  " ", tmp_s[3], tmp_s[4], tmp_s[5],
                  tmp_s[6], tmp_s[7], tmp_s[8], tmp_s[9]);
  }

#endif

  if (cs_glob_n_ranks == 1) {

    char tmp_s[3][24] =  {"", "", ""};

    cs_log_strpadl(tmp_s[0], _("time (s)"), 12, 24);
    cs_log_strpadl(tmp_s[1], _("speedup"), 8, 24);
    cs_log_strpadl(tmp_s[2], _("std. dev."), 8, 24);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %24s %s %s %s\n",
                  " ", tmp_s[0], tmp_s[1], tmp_s[2]);

  }
}

/*----------------------------------------------------------------------------
 * Print statistics on matrix tuning SpMv info.
 *
 * parameters:
 *   m_variant   <-- array of matrix variants
 *   variant_id  <-- variant id
 *   fill_type   <-- type of matrix fill
 *   ed_flag     <-- 0: include diagonal; 1: exclude diagonal
 *----------------------------------------------------------------------------*/

static void
_matrix_tune_spmv_stats(const cs_matrix_variant_t  *m_variant,
                        int                         variant_id,
                        cs_matrix_fill_type_t       fill_type,
                        int                         ed_flag)
{
  char title[32];

  double v_loc[3] = {-1, -1, 0};

  const cs_matrix_variant_t  *r = m_variant;
  const cs_matrix_variant_t  *v = m_variant + variant_id;

  cs_log_strpad(title, v->name, 24, 32);

  /* Get timing info */

  v_loc[0] = v->matrix_vector_cost[fill_type][ed_flag][0];
  v_loc[1] = r->matrix_vector_cost[fill_type][ed_flag][0];
  v_loc[2] = v->matrix_vector_cost[fill_type][ed_flag][1];

  if (v_loc[0] < 0)
    return;

  /* Compute local ratios */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double v_max[3], speedup_min, v_sum[3];

    MPI_Allreduce(v_loc, v_sum, 3, MPI_DOUBLE, MPI_SUM, cs_glob_mpi_comm);

    if (v_loc[0] > 1e-9)
      v_loc[1] = r->matrix_vector_cost[fill_type][ed_flag][0] / v_loc[0];
    else
      v_loc[1] = 1;
    MPI_Allreduce(v_loc + 1, &speedup_min, 1, MPI_DOUBLE, MPI_MIN,
                  cs_glob_mpi_comm);

    if (v_loc[0] > 1e-9)
      v_loc[2] = v->matrix_vector_cost[fill_type][ed_flag][1] / v_loc[0];
    else {
      v_loc[1] = 0;
      v_loc[2] = 0;
    }
    MPI_Allreduce(v_loc, v_max, 3, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    cs_real_t t_mean = v_sum[0]/cs_glob_n_ranks;
    cs_real_t speedup_mean = v_sum[1] / v_sum[0]; /* weighted by v_loc[0] */
    cs_real_t stddev_mean = v_sum[2]/v_sum[0];
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %12.5e %12.5e %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                  title,
                  t_mean, v_max[0], speedup_mean, speedup_min, v_max[1],
                  stddev_mean, v_max[2]);
  }

#endif

  if (cs_glob_n_ranks == 1) {
    cs_real_t speedup = v_loc[1] / v_loc[0];
    cs_real_t stddev = v_loc[2] / v_loc[0];
    cs_log_printf(CS_LOG_PERFORMANCE,
                  "  %s %12.5e %8.4f %8.4f\n",
                  title,
                  v_loc[0], speedup, stddev);
  }
}

/*----------------------------------------------------------------------------
 * Initialize local variant matrix.
 *
 * parameters:
 *   v  <-> pointer to matrix variant
 *----------------------------------------------------------------------------*/

static void
_variant_init(cs_matrix_variant_t  *v)
{
  for (int i = 0; i < CS_MATRIX_N_FILL_TYPES; i++) {
    for (int j = 0; j < 2; j++) {
      v->vector_multiply[i][j] = NULL;
      v->matrix_vector_cost[i][j][0] = -1.;
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Tune local matrix.vector product operations.
 *
 * To avoid multiplying structures for multiple matrix fill-ins,
 * an array of tuning types may be provided, and weights may be
 * associated to each type based on the expected usage of each fill-in
 * type. If n_fill_types is set to 0, these arrays are ignored, and their
 * following default is used:
 *
 *   CS_MATRIX_SCALAR      0.5
 *   CS_MATRIX_SCALAR_SYM  0.25
 *   CS_MATRIX_33_BLOCK_D  0.25
 *
 * parameters:
 *   t_measure      <-- minimum time for each measure
 *   n_types        <-- number of matrix types tuned for, or 0
 *   n_fill_types   <-- number of fill types tuned for, or 0
 *   types          <-- array of matrix types tuned for, or NULL
 *   fill_types     <-- array of fill types tuned for, or NULL
 *   fill_weights   <-- weight of fill types tuned for, or NULL
 *   n_cells        <-- number of local cells
 *   n_cells_ext    <-- number of cells including ghost cells (array size)
 *   n_faces        <-- local number of internal faces
 *   cell_num       <-- Optional global cell numbers (1 to n), or NULL
 *   face_cell      <-- face -> cells connectivity
 *   halo           <-- cell halo structure
 *   numbering      <-- vectorization or thread-related numbering info, or NULL
 *
 * returns:
 *   pointer to tuning results structure
 *----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_tuned(double                 t_measure,
                        int                    n_types,
                        int                    n_fill_types,
                        cs_matrix_type_t       types[],
                        cs_matrix_fill_type_t  fill_types[],
                        double                 fill_weights[],
                        cs_lnum_t              n_cells,
                        cs_lnum_t              n_cells_ext,
                        cs_lnum_t              n_faces,
                        const cs_lnum_2_t     *face_cell,
                        const cs_halo_t       *halo,
                        const cs_numbering_t  *numbering)
{
  int  t_id, t_id_max, f_id, v_id, ed_flag;

  double speedup, max_speedup;
  double t_speedup[CS_MATRIX_N_TYPES][CS_MATRIX_N_FILL_TYPES];
  int cur_select[CS_MATRIX_N_FILL_TYPES][2];

  bool                   type_filter[CS_MATRIX_N_TYPES] = {true,
                                                           true,
                                                           true,
                                                           true};

  int                    _n_fill_types_default = 3;
  cs_matrix_fill_type_t  _fill_types_default[] = {CS_MATRIX_SCALAR,
                                                  CS_MATRIX_SCALAR_SYM,
                                                  CS_MATRIX_BLOCK_D};
  double                 _fill_weights_default[] = {0.5, 0.25, 0.25};

  int                    _n_types = n_types;
  int                    _n_fill_types = n_fill_types;
  cs_matrix_fill_type_t  *_fill_types = fill_types;
  double                 *_fill_weights = fill_weights;
  double                  tot_weight = 0.;

  int  n_variants = 0;
  cs_matrix_variant_t  *m_variant = NULL, *v = NULL;

  cs_matrix_variant_t  *r = NULL;

  cs_timer_t           t0, t1;

  t0 = cs_timer_time();

  /* Use defaults if required */

  if (_n_types > 0) {
    for (t_id = 0; t_id < CS_MATRIX_N_TYPES; t_id++)
      type_filter[t_id] = false;
    for (t_id = 0; t_id < n_types; t_id++)
      type_filter[types[t_id]] = true;
  }

  if (_n_fill_types < 1) {
    _n_fill_types =  _n_fill_types_default;
    _fill_types = _fill_types_default;
    _fill_weights = _fill_weights_default;
  }

  /* Base flags on weights */

  for (t_id = 0; t_id < CS_MATRIX_N_TYPES; t_id++) {
    for (cs_matrix_fill_type_t fill_type = 0;
         fill_type < CS_MATRIX_N_FILL_TYPES;
         fill_type++) {
      t_speedup[t_id][fill_type] = -1;
    }
  }

  /* Build variants array */
  /*----------------------*/

  cs_matrix_variant_build_list(_n_fill_types,
                               _fill_types,
                               type_filter,
                               numbering,
                               &n_variants,
                               &m_variant);

  /* Run tests on variants */

  _matrix_tune_test(t_measure,
                    n_variants,
                    n_cells,
                    n_cells_ext,
                    n_faces,
                    face_cell,
                    halo,
                    numbering,
                    m_variant);

  /* Print info on variants */

  for (f_id = 0; f_id < _n_fill_types; f_id++) {
    cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
    tot_weight += _fill_weights[f_id];
    for (ed_flag = 0; ed_flag < 2; ed_flag++) {
      _matrix_tune_spmv_title(fill_type, ed_flag);
      for (v_id = 0; v_id < n_variants; v_id++)
        _matrix_tune_spmv_stats(m_variant,
                                v_id,
                                fill_type,
                                ed_flag);
    }
  }

  /* Select type of matrix with best possible performance */

  for (v_id = 0; v_id < n_variants; v_id++) {
    v = m_variant + v_id;
    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
      speedup =   m_variant->matrix_vector_cost[fill_type][0][0]
                / v->matrix_vector_cost[fill_type][0][0];
      if (t_speedup[v->type][fill_type] < speedup)
        t_speedup[v->type][fill_type] = speedup;
    }
  }

  max_speedup = 0;
  t_id_max = 0;

  for (t_id = 0; t_id < CS_MATRIX_N_TYPES; t_id++) {
    speedup = 0.;
    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
      speedup += fill_weights[f_id] * t_speedup[t_id][fill_type];
    }
    speedup /= tot_weight;
    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
      if (t_speedup[t_id][fill_type] < 0)
        speedup = -1;
    }
    if (speedup > max_speedup) {
      max_speedup = speedup;
      t_id_max = t_id;
    }
  }

  /* Now that the best type is chosen, build the variant */

  BFT_MALLOC(r, 1, cs_matrix_variant_t);

  _variant_init(r);

  strncpy(r->name, cs_matrix_type_name[t_id_max], 31);
  r->type = t_id_max;

  for (f_id = 0; f_id < _n_fill_types; f_id++) {
    for (ed_flag = 0; ed_flag < 2; ed_flag++)
      cur_select[f_id][ed_flag] = -1;
  }

  for (v_id = 0; v_id < n_variants; v_id++) {

    v = m_variant + v_id;
    if (v->type != r->type)
      continue;

    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      for (ed_flag = 1; ed_flag >= 0; ed_flag--) { /* full matrix priority */

        cs_matrix_fill_type_t  fill_type = _fill_types[f_id];

        if (v->matrix_vector_cost[fill_type][ed_flag][0] > 0) {
          if (   v->matrix_vector_cost[fill_type][ed_flag][0]
               < r->matrix_vector_cost[fill_type][ed_flag][0]
              || r->matrix_vector_cost[fill_type][ed_flag][0] < 0) {
            r->vector_multiply[fill_type][ed_flag]
              = v->vector_multiply[fill_type][ed_flag];
            for (int o_id = 0; o_id < 1; o_id++)
              r->matrix_vector_cost[fill_type][ed_flag][o_id]
                = v->matrix_vector_cost[fill_type][ed_flag][o_id];
            cur_select[f_id][ed_flag] = v_id;
          }
        }

      }
    }

  } /* End of loop on variants */

  /* print info on selected variants */

  cs_log_printf(CS_LOG_PERFORMANCE,
                "\n"
                "Selected matrix operation implementations:\n"
                "------------------------------------------\n");

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int step = _n_fill_types*2;
    int *select_loc, *select_sum;

    BFT_MALLOC(select_sum, n_variants*step, int);
    BFT_MALLOC(select_loc, n_variants*step, int);

    for (v_id = 0; v_id < n_variants; v_id++) {
      for (int sub_id = 0; sub_id < step; sub_id++)
        select_loc[v_id*step + sub_id] = 0;
    }
    for (int sub_id = 0; sub_id < step; sub_id++) {
      f_id = sub_id / 2;
      ed_flag = sub_id %2;
      if (cur_select[f_id][ed_flag] > -1)
        select_loc[cur_select[f_id][ed_flag]*step + sub_id] = 1;
    }

    MPI_Allreduce(select_loc, select_sum, n_variants*step, MPI_INT, MPI_SUM,
                  cs_glob_mpi_comm);

    BFT_FREE(select_loc);

    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      for (ed_flag = 0; ed_flag < 2; ed_flag++) {

        int count_tot = 0;

        for (v_id = 0; v_id < n_variants; v_id++)
          count_tot += (select_sum[v_id*step + f_id*2 + ed_flag]);

        if (count_tot > 0) {
          cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
          cs_log_printf(CS_LOG_PERFORMANCE,
                        _("\n  -%s:\n"),
                        _(_matrix_operation_name[fill_type][ed_flag]));
          for (v_id = 0; v_id < n_variants; v_id++) {
            int scount = select_sum[v_id*step + f_id*2 + ed_flag];
            if (scount > 0) {
              char title[36] =  {""};
              v = m_variant + v_id;
              cs_log_strpad(title, _(v->name), 32, 36);
              cs_log_printf(CS_LOG_PERFORMANCE,
                            _("    %s : %d ranks\n"), title, scount);
            }
          }
        }

      }
    }

    BFT_FREE(select_sum);

  } /* if (cs_glob_n_ranks > 1) */

#endif

  if (cs_glob_n_ranks == 1) {

    cs_log_printf(CS_LOG_PERFORMANCE, "\n");

    for (f_id = 0; f_id < _n_fill_types; f_id++) {
      cs_matrix_fill_type_t  fill_type = _fill_types[f_id];
      for (ed_flag = 0; ed_flag < 2; ed_flag++) {
        v_id = cur_select[f_id][ed_flag];
        if (v_id > -1) {
          v = m_variant + v_id;
          cs_log_printf(CS_LOG_PERFORMANCE,
                        _("  %-44s : %s\n"),
                        _(_matrix_operation_name[fill_type][ed_flag]),
                        _(v->name));
        }

      }
    }

  } /* if (cs_glob_n_ranks == 1) */

  BFT_FREE(m_variant);

  t1 = cs_timer_time();

  cs_timer_counter_t  t_tune = cs_timer_diff(&t0, &t1);

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n"
                  "Elapsed tuning time: %12.3f s.\n"),
                t_tune.wall_nsec*1e-9);

  return r;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
