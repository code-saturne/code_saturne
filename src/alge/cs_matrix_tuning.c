/*============================================================================
 * Sparse Matrix Representation and Operations.
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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Tune local matrix.vector product operations.
 *
 * parameters:
 *   m           <-- matrix to tune
 *   n_measure   <-- minimum number of measures
 *   t_measure   <-- minimum time for each measure
 *   n_variants  <-- number of variants in array
 *   m_variant   <-- array of matrix variants
 *   spmv_cost   --> SpMV cost
 *----------------------------------------------------------------------------*/

static void
_matrix_tune_test(const cs_matrix_t     *m,
                  int                    n_measure,
                  double                 t_measure,
                  int                    n_variants,
                  cs_matrix_variant_t   *m_variant,
                  double                 spmv_cost[])
{
  int  n_runs, run_id, v_id;
  double  wt0, wt1, wtu;

  double test_sum = 0.0;
  cs_real_t  *x = NULL, *y = NULL;

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  cs_lnum_t n_cols = cs_matrix_get_n_columns(m);
  cs_lnum_t b_size = cs_matrix_get_diag_block_size(m);

  cs_lnum_t n = n_cols*b_size;

  CS_MALLOC_HD(x, n, cs_real_t, m->alloc_mode);
  CS_MALLOC_HD(y, n, cs_real_t, m->alloc_mode);

# pragma omp parallel for  if(n > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n; i++) {
    x[i] = 1.0;
    y[i] = 0.0;
  }

  /* Loop on variant types */
  /*-----------------------*/

  for (v_id = 0; v_id < n_variants; v_id++) {

    const cs_matrix_variant_t *v = m_variant + v_id;

    /* Measure matrix.vector operations */

    for (int op_type = 0; op_type < CS_MATRIX_SPMV_N_TYPES; op_type++) {

      cs_matrix_vector_product_t
        *vector_multiply = v->vector_multiply[op_type];

      if (vector_multiply != NULL) {

        cs_matrix_t m_t;
        memcpy(&m_t, m, sizeof(cs_matrix_t));

        m_t.vector_multiply[m->fill_type][op_type] = vector_multiply;

        wt0 = cs_timer_wtime(), wt1 = wt0;
        run_id = 0, n_runs = (n_measure > 0) ? n_measure : 1;

        while (run_id < n_runs) {
          while (run_id < n_runs) {
            if (run_id % 8)
              test_sum = 0;
            if (op_type == 0)
              cs_matrix_vector_multiply(&m_t, x, y);
            else
              cs_matrix_vector_multiply_partial(&m_t, op_type, x, y);
            test_sum += y[n-1];
            run_id++;
          }
          wt1 = cs_timer_wtime();
          double wt_r0 = wt1 - wt0;

#if defined(HAVE_MPI)

          if (cs_glob_n_ranks > 1) {
            double _wt_r0 = wt_r0;
            MPI_Allreduce(&_wt_r0, &wt_r0, 1, MPI_DOUBLE, MPI_MAX,
                          cs_glob_mpi_comm);
          }

#endif /* defined(HAVE_MPI) */

          if (wt_r0 < t_measure)
            n_runs *= 2;
        }
        wtu = (wt1 - wt0) / n_runs;
        spmv_cost[v_id*CS_MATRIX_SPMV_N_TYPES + op_type] = wtu;

        if (m_t.destroy_adaptor != NULL)
          m_t.destroy_adaptor(&m_t);

      }
      else
        spmv_cost[v_id*CS_MATRIX_SPMV_N_TYPES + op_type] = -1;

    } /* end of loop on ed_flag */

  } /* end of loop on variants */

  CS_FREE_HD(x);
  CS_FREE_HD(y);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select spmv variant with best performance
 *
 * The first variant of the list is modified to select the function
 * with best performance.
 *
 * \param[in]   m             associated matrix
 * \param[in]   verbosity     verbosity level
 * \param[in]   fill_type     fill type tuned for
 * \param[in]   n_variants    number of variants
 * \param[in]   n_r_variants  number of return variants (1 on host only,
 *                            3 with device)
 * \param[in]   m_variant     array of input matrix variants
 * \param[out]  r_variant     array of output matrix variants
 * \param[in]   spmv_cost     costs for each variant
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_tune_spmv_select(const cs_matrix_t    *m,
                         int                   verbosity,
                         int                   n_variants,
                         int                   n_r_variants,
                         cs_matrix_variant_t  *m_variant,
                         cs_matrix_variant_t  *r_variant,
                         double                spmv_cost[])
{
  /* Use maximum value for comparisons */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int     n = n_variants*2;
    double *cost_local;
    BFT_MALLOC(cost_local, n, double);
    for (int i = 0; i < n; i++)
      cost_local[i] = spmv_cost[i];

    MPI_Allreduce(cost_local, spmv_cost, n, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    BFT_FREE(cost_local);

  }

#endif

  int min_c[3][2] = {{-1, -1}, {-1, -1}, {-1, -1}};

  for (int i = 0; i < n_variants; i++) {
    const cs_matrix_variant_t *mv = m_variant+i;
    for (int j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
      if (spmv_cost[i*CS_MATRIX_SPMV_N_TYPES + j] > 0) {
        int k = 1;
        if (mv->vector_multiply_xy_hd[j] == 'd')
          k = 2;
        if (min_c[k][j] < 0)
          min_c[k][j] = i;
        else if (  spmv_cost[i*CS_MATRIX_SPMV_N_TYPES + j]
                 < spmv_cost[min_c[k][j]*CS_MATRIX_SPMV_N_TYPES + j])
          min_c[k][j] = i;
      }
    }
  }

  for (int j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
    if (min_c[1][j] > -1)
      min_c[0][j] = min_c[1][j];
    if (min_c[2][j] > -1) {
      if (min_c[0][j] > -1) {
        if (  spmv_cost[min_c[2][j]*CS_MATRIX_SPMV_N_TYPES + j]
            < spmv_cost[min_c[0][j]*CS_MATRIX_SPMV_N_TYPES + j])
          min_c[0][j] = min_c[2][j];
      }
      else
        min_c[0][j] = min_c[2][j];
    }
  }

  for (int k = 0; k < n_r_variants; k++) {
    cs_matrix_variant_t *o_variant = r_variant + k;
    o_variant->fill_type = m->fill_type;

    for (int j = 0; j < CS_MATRIX_SPMV_N_TYPES; j++) {
      if (min_c[k][j] > -1) {
        const cs_matrix_variant_t *mv_s = m_variant+min_c[k][j];
        strcpy(o_variant->name[j], mv_s->name[j]);
        o_variant->vector_multiply[j] = mv_s->vector_multiply[j];
        o_variant->vector_multiply_xy_hd[j] = mv_s->vector_multiply_xy_hd[j];
      }
    }
  }

  if (verbosity > 0) {
    const char* hd_type[] = {"", "host ", "device "};
    for (int k = 0; k < n_r_variants; k++) {
      cs_matrix_variant_t *o_variant = r_variant + k;
      cs_log_printf
        (CS_LOG_PERFORMANCE,
         _("\n"
           "Selected %sSpMV variant for matrix of type %s and fill %s:\n"),
         hd_type[k],
         _(cs_matrix_get_type_name(m)),
         _(cs_matrix_fill_type_name[m->fill_type]));
      if (min_c[k][0] > -1)
        cs_log_printf
          (CS_LOG_PERFORMANCE,
           _("  %32s for y <= A.x       (speedup: %6.2f)\n"),
           o_variant->name[0],
           spmv_cost[0]/spmv_cost[min_c[k][0]*CS_MATRIX_SPMV_N_TYPES]);
      if (min_c[k][1] > -1)
        cs_log_printf
          (CS_LOG_PERFORMANCE,
           _("  %32s for y <= (A-D).x   (speedup: %6.2f)\n"),
           o_variant->name[1],
           spmv_cost[1]/spmv_cost[min_c[k][1]*CS_MATRIX_SPMV_N_TYPES+1]);
    }
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a matrix variant tuned matrix.vector product operations.
 *
 * The variant may later be applied to matrices of the same type and fill type.
 *
 * In presence of supported accelerated devices, an array of 3 variants
 * is returned; the second one applies to the host only, the third one
 * to the device only.
 *
 * \param[in]  m           associated matrix
 * \param[in]  verbosity   verbosity level
 * \param[in]  n_measure   minimum number of measuring runs
 * \param[in]  t_measure   minimum measure time
 *
 * \returns  pointer to tuning results structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_variant_t *
cs_matrix_variant_tuned(const cs_matrix_t  *m,
                        int                 verbosity,
                        int                 n_measure,
                        double              t_measure)
{
  int  n_variants = 0, n_r_variants = 1;
  cs_matrix_variant_t  *r_variant = NULL;
  cs_matrix_variant_t  *m_variant = NULL;

  if (cs_get_device_id() > -1)
    n_r_variants = 3;

  BFT_MALLOC(r_variant, n_r_variants, cs_matrix_variant_t);

  cs_matrix_variant_build_list(m, &n_variants, &m_variant);

  if (n_variants > 1) {

    if (verbosity > 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("\n"
                      "Tuning for matrices of type %s and fill %s\n"
                      "===========================\n"),
                    cs_matrix_get_type_name(m),
                    cs_matrix_fill_type_name[m->fill_type]);

    double *spmv_cost;
    BFT_MALLOC(spmv_cost, n_variants*CS_MATRIX_SPMV_N_TYPES, double);

    _matrix_tune_test(m,
                      n_measure,
                      t_measure,
                      n_variants,
                      m_variant,
                      spmv_cost);

    _matrix_tune_spmv_select(m,
                             verbosity,
                             n_variants,
                             n_r_variants,
                             m_variant,
                             r_variant,
                             spmv_cost);

    BFT_FREE(spmv_cost);

    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

  }

  /* Single-variant case */

  else
    memcpy(r_variant, m_variant, sizeof(cs_matrix_variant_t));

  BFT_FREE(m_variant);

  return r_variant;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
