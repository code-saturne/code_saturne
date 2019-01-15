/*============================================================================
 * Routines shared among all face-based schemes for the discretization of the
 * Navier-Stokes system
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_cdo_bc.h"
#include "cs_equation_bc.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_post.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_navsto.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_navsto.c
 *
 * \brief Shared routines among all face-based schemes for building and solving
 *        Stokes and Navier-Stokes problem
 *
 */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_NAVSTO_DBG      0

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
 * \brief  Compute the divergence of a cell using the \ref cs_cdo_quantities_t
 *         structure
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_dof        values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdofb_navsto_cell_divergence(const cs_lnum_t               c_id,
                                const cs_cdo_quantities_t    *quant,
                                const cs_adjacency_t         *c2f,
                                const cs_real_t              *f_dof)
{
  cs_real_t  div = 0.0;

  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

    const cs_lnum_t  f_id = c2f->ids[f];
    const cs_real_t  *_val = f_dof + 3*f_id;

    if (f_id < quant->n_i_faces) {
      const cs_real_t *_nuf = quant->i_face_normal + 3*f_id;

      div += c2f->sgn[f]*quant->i_face_surf[f_id]*
        cs_math_3_dot_product(_val, _nuf) / cs_math_3_norm(_nuf);

    }
    else {

      const cs_lnum_t  bf_id = f_id - quant->n_i_faces;
      const cs_real_t  *_nuf = quant->b_face_normal + 3*bf_id;

      div += c2f->sgn[f]*quant->b_face_surf[bf_id]*
        cs_math_3_dot_product(_val, _nuf) / cs_math_3_norm(_nuf);

    } /* Boundary face */

  } /* Loop on cell faces */

  div /= quant->cell_vol[c_id];

  return div;
}

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
cs_cdofb_navsto_add_grad_div(short int          n_fc,
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
 * \brief  Initialize the pressure values
 *
 * \param[in]       nsp     pointer to a \ref cs_navsto_param_t structure
 * \param[in]       quant   pointer to a \ref cs_cdo_quantities_t structure
 * \param[in]       ts      pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr      pointer to the pressure \ref cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_pressure(const cs_navsto_param_t     *nsp,
                              const cs_cdo_quantities_t   *quant,
                              const cs_time_step_t        *ts,
                              cs_field_t                  *pr)
{
  /* Sanity checks */
  assert(nsp != NULL);

  /* Initial conditions for the pressure */
  if (nsp->n_pressure_ic_defs == 0)
    return; /* Nothing to do */

  assert(nsp->pressure_ic_defs != NULL);

  const cs_real_t  t_cur = ts->t_cur;
  const cs_flag_t  dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;

  cs_real_t  *values = pr->val;

  for (int def_id = 0; def_id < nsp->n_pressure_ic_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    cs_xdef_t  *def = nsp->pressure_ic_defs[def_id];

    /* Initialize face-based array */
    switch (def->type) {

      /* Evaluating the integrals: the averages will be taken care of at the
       * end when ensuring zero-mean valuedness */
    case CS_XDEF_BY_VALUE:
      cs_evaluate_density_by_value(dof_flag, def, values);
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;

        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          /* Forcing BARY so that it is equivalent to DeRham (JB?)*/
          cs_xdef_set_quadrature(def, CS_QUADRATURE_BARY);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          /* Restoring the original */
          cs_xdef_set_quadrature(def, nsp->qtype);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Incompatible reduction for the field %s.\n"),
                    __func__, pr->name);

        }  /* Switch on possible reduction types */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Incompatible way to initialize the field %s.\n"),
                __func__, pr->name);
      break;

    }  /* Switch on possible type of definition */

  }  /* Loop on definitions */

  /* We should ensure that the mean of the pressure is zero. Thus we compute
   * it and subtract it from every value.
   *
   * NOTES:
   *  - It could be useful to stored this average somewhere
   *  - The procedure is not optimized (we can avoid setting the average if
   *    it's a value), but it is the only way to allow multiple definitions
   *    and definitions that do not cover all the domain. Moreover, we need
   *    information (e.g. cs_cdo_quantities_t) which we do not know here
   */
  cs_cdofb_navsto_set_zero_mean_pressure(quant, values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a zero-mean
 *         average
 *
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_set_zero_mean_pressure(const cs_cdo_quantities_t  *quant,
                                       cs_real_t                   values[])
{
  /* We should ensure that the mean of the pressure is zero. Thus we compute
   * it and subtract it from every value. */
  /* NOTES:
   *  - It could be useful to stored this average somewhere
   *  - The procedure is not optimized (we can avoid setting the average if
   *    it's a value), but it is the only way to allow multiple definitions
   *    and definitions that do not cover all the domain. */

  const cs_lnum_t  n_cells = quant->n_cells;
 /*
  * The algorithm used for summing is l3superblock60, based on the article:
  * "Reducing Floating Point Error in Dot Product Using the Superblock Family
  * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
  * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
  * 2008 Society for Industrial and Applied Mathematics
  */

  const cs_real_t  intgr = cs_sum(n_cells, values);
  assert(quant->vol_tot > 0.);
  const cs_real_t  g_avg = intgr / quant->vol_tot;

  const cs_real_t *cv = quant->cell_vol;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    values[c_id] = values[c_id] / cv[c_id] - g_avg;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
