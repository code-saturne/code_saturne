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
  int  n_runs, run_id, v_id, ed_flag;
  double  wt0, wt1, wtu;

  double test_sum = 0.0;
  cs_real_t  *x = NULL, *y = NULL;

  /* Allocate and initialize  working arrays */
  /*-----------------------------------------*/

  cs_lnum_t n_cols = cs_matrix_get_n_columns(m);
  cs_lnum_t b_size = cs_matrix_get_diag_block_size(m)[0];

  cs_lnum_t n = n_cols*b_size;

  if (CS_MEM_ALIGN > 0) {
    BFT_MEMALIGN(x, CS_MEM_ALIGN, n, cs_real_t);
    BFT_MEMALIGN(y, CS_MEM_ALIGN, n, cs_real_t);
  }
  else {
    BFT_MALLOC(x, n, cs_real_t);
    BFT_MALLOC(y, n, cs_real_t);
  }

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

    for (ed_flag = 0; ed_flag < 2; ed_flag++) {

      cs_matrix_vector_product_t
        *vector_multiply = v->vector_multiply[ed_flag];

      if (vector_multiply == NULL)
        continue;

      if (vector_multiply != NULL) {

        cs_matrix_t m_t;
        memcpy(&m_t, m, sizeof(cs_matrix_t));

        m_t.vector_multiply[m->fill_type][ed_flag] = vector_multiply;

        wt0 = cs_timer_wtime(), wt1 = wt0;
        run_id = 0, n_runs = (n_measure > 0) ? n_measure : 1;

        while (run_id < n_runs) {
          while (run_id < n_runs) {
            if (run_id % 8)
              test_sum = 0;
            if (ed_flag == 0)
              cs_matrix_vector_multiply(CS_HALO_ROTATION_COPY, &m_t, x, y);
            else
              cs_matrix_exdiag_vector_multiply(CS_HALO_ROTATION_COPY, &m_t, x, y);
            test_sum += y[n-1];
            run_id++;
          }
          wt1 = cs_timer_wtime();
          double wt_r0 = wt1 - wt0;
          if (cs_glob_n_ranks > 1) {
            double _wt_r0 = wt_r0;
            MPI_Allreduce(&_wt_r0, &wt_r0, 1, MPI_DOUBLE, MPI_MAX,
                          cs_glob_mpi_comm);
          }
          if (wt_r0 < t_measure)
            n_runs *= 2;
        }
        wtu = (wt1 - wt0) / n_runs;
        spmv_cost[v_id*2 + ed_flag] = wtu;
      }

    } /* end of loop on ed_flag */

  } /* end of loop on variants */

  BFT_FREE(x);
  BFT_FREE(y);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select spmv variant with best performance
 *
 * The first variant of the list is modified to select the function
 * with best performance.
 *
 * \param[in]       m             associated matrix
 * \param[in]       verbosity     verbosity level
 * \param[in]       fill_type     fill type tuned for
 * \param[in]       n_variants    number of variants
 * \param[in, out]  m_variant     array of matrix variants
 * \param[in]       spmv_cost     costs for each variant
 */
/*----------------------------------------------------------------------------*/

static void
_matrix_tune_spmv_select(const cs_matrix_t    *m,
                         int                   verbosity,
                         int                   n_variants,
                         cs_matrix_variant_t  *m_variant,
                         double                spmv_cost[])
{
  /* Use maximum value for comparisons */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    int     n = n_variants*2;
    double *cost_local;
    BFT_MALLOC(cost_local, n, double);
    for (int i = 0; i < n_variants; i++)
      cost_local[i] = spmv_cost[i];

    MPI_Allreduce(cost_local, spmv_cost, n, MPI_DOUBLE, MPI_MAX, cs_glob_mpi_comm);

    BFT_FREE(cost_local);

  }

#endif

  int min_c[2] = {0, 0};

  for (int i = 1; i < n_variants; i++) {
    for (int j = 0; j < 2; j++) {
      if (spmv_cost[i*2 + j] < spmv_cost[min_c[j]*2 + j])
        min_c[j] = i;
    }
  }

  for (int j = 0; j < 2; j++) {
    if (spmv_cost[min_c[j]*2 + j] < spmv_cost[j]) {
      const cs_matrix_variant_t *mv_s = m_variant+min_c[j];
      strcpy(m_variant->name[j], mv_s->name[j]);
      m_variant->vector_multiply[j] = mv_s->vector_multiply[j];
    }
  }

  if (verbosity > 0) {
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n"
                    "Selected SpMV variant for matrix of type %s and fill %s:\n"
                    "  %32s for y <= A.x       (speedup: %6.2f)\n"
                    "  %32s for y <= (A-D).x   (speedup: %6.2f)\n"),
                  _(cs_matrix_type_name[m->type]),
                  _(cs_matrix_fill_type_name[m->fill_type]),
                  m_variant[0].name[0], spmv_cost[0]/spmv_cost[min_c[0]*2],
                  m_variant[0].name[1], spmv_cost[1]/spmv_cost[min_c[1]*2+1]);
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
  int  n_variants = 0;
  cs_matrix_variant_t  *m_variant = NULL;

  cs_matrix_variant_build_list(m, &n_variants, &m_variant);

  if (n_variants > 1) {

    if (verbosity > 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("\n"
                      "Tuning for matrices of type %s and fill %s\n"
                      "===========================\n"),
                    cs_matrix_type_name[m->type],
                    cs_matrix_fill_type_name[m->fill_type]);

    double *spmv_cost;
    BFT_MALLOC(spmv_cost, n_variants*2, double);

    _matrix_tune_test(m,
                      n_measure,
                      t_measure,
                      n_variants,
                      m_variant,
                      spmv_cost);

    _matrix_tune_spmv_select(m,
                             verbosity,
                             n_variants,
                             m_variant,
                             spmv_cost);

    BFT_FREE(spmv_cost);

    cs_log_printf(CS_LOG_PERFORMANCE, "\n");
    cs_log_separator(CS_LOG_PERFORMANCE);

    n_variants = 1;
    BFT_REALLOC(m_variant, 1, cs_matrix_variant_t);

  }

  return m_variant;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
