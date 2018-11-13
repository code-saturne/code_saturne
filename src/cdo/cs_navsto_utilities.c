/*============================================================================
 * Routines common to all the velocity-pressure couplings
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_sles.h"
#include "cs_equation_priv.h"
#include "cs_log.h"
#include "cs_matrix.h"
#include "cs_dbg.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_parall.h"
#include "cs_range_set.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_utilities.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_NAVSTO_UTILITIES_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the grad-div part to the local matrix (i.e. for the current
 *         cell)
 *
 * \param[in]      n_fc       local number of faces for the current cell
 * \param[in]      zeta       scalar coefficient for the grad-div operator
 * \param[in]      div        divergence
 * \param[in, out] mat        local system matrix to update
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_grad_div(short int          n_fc,
                       const cs_real_t    zeta,
                       const cs_real_t    div[],
                       cs_sdm_t          *mat)
{
  cs_sdm_t  *b = NULL;

  /* Avoid dealing with cell DoFs which are not impacted */
  for (short int bi = 0; bi < n_fc; bi++) {

    const cs_real_t  *divi = div + 3*bi;
    const cs_real_t  zt_di[3] = {zeta*divi[0], zeta*divi[1], zeta*divi[2]};

    /* Begin with the diagonal block */
    b = cs_sdm_get_block(mat, bi, bi);
    assert(b->n_rows == b->n_cols && b->n_rows == 3);
    for (short int l = 0; l < 3; l++) {
      cs_real_t *m_l = b->val + 3*l;
      for (short int m = 0; m < 3; m++)
        m_l[m] += zt_di[l] * divi[m];
    }

    /* Continue with the extra-diag. blocks */
    for (short int bj = bi+1; bj < n_fc; bj++) {

      b = cs_sdm_get_block(mat, bi, bj);
      assert(b->n_rows == b->n_cols && b->n_rows == 3);
      cs_real_t *mij  = b->val;
      b = cs_sdm_get_block(mat, bj, bi);
      assert(b->n_rows == b->n_cols && b->n_rows == 3);
      cs_real_t *mji  = b->val;

      const cs_real_t *divj = div + 3*bj;

      for (short int l = 0; l < 3; l++) {

        /* Diagonal: 3*l+l = 4*l */
        const cs_real_t  gd_coef_ll = zt_di[l]*divj[l];
        mij[4*l] += gd_coef_ll;
        mji[4*l] += gd_coef_ll;

        /* Extra-diagonal: Use the symmetry of the grad-div */
        for (short int m = l+1; m < 3; m++) {
          const short int  lm = 3*l+m, ml = 3*m+l;
          const cs_real_t  gd_coef_lm = zt_di[l]*divj[m];
          mij[lm] += gd_coef_lm;
          mji[ml] += gd_coef_lm;
          const cs_real_t  gd_coef_ml = zt_di[m]*divj[l];
          mij[ml] += gd_coef_ml;
          mji[lm] += gd_coef_ml;
        }
      }

    } /* Loop on column blocks: bj */
  } /* Loop on row blocks: bi */

}

/*----------------------------------------------------------------------------*/
 /*!
 * \brief  Solve the linear system (adaptation of cs_equation_solve())
 *
 * \param[in, out]  eq    pointer to the momentum cs_equation_t structure
 * \param[in, out]  x     pointer to temporary velocity on faces
 * \param[in, out]  b     pointer to auxiliary rhs, should be NULL
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_navsto_solve(cs_equation_t   *eq,
                cs_real_t       *x,
                cs_real_t       *b)
{
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);

  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_itsol_t  itsol_info = eq->param->itsol_info;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   eq->matrix,
                                                   CS_HALO_ROTATION_IGNORE,
                                                   itsol_info.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   b,
                                                   x,
                                                   0,      /* aux. size */
                                                   NULL);  /* aux. buffers */

  if (eq->param->sles_verbosity > 0) {

    const cs_lnum_t  size = eq->n_sles_gather_elts;
    const cs_lnum_t  *row_index, *col_id;
    const cs_real_t  *d_val, *x_val;

    cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, &d_val, &x_val);

    cs_gnum_t  nnz = row_index[size];
    if (cs_glob_n_ranks > 1)
      cs_parall_counter(&nnz, 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  "  <%s/sles_cvg> code  %-d n_iters  %d residual  % -8.4e"
                  " nnz %lu\n",
                  eq->param->name, code, n_iters, residual, nnz);

#if defined(DEBUG) && !defined(NDEBUG) && CS_NAVSTO_UTILITIES_DBG
    cs_dbg_dump_linear_system(eq->param->name, size, CS_NAVSTO_UTILITIES_DBG,
                              x, b,
                              row_index, col_id, x_val, d_val);
#endif
  }

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(eq->rset, CS_REAL_TYPE, 1, /* type and stride */
                         x, x);                     /* in/out */

    cs_range_set_scatter(eq->rset, CS_REAL_TYPE, 1, /* type and stride */
                         b, eq->rhs);               /* in/out */

  }

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  return  n_iters; /* Number of inner iterations for solving the linear system
                      at this step */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
