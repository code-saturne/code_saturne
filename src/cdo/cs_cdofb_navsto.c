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
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sdm.h"
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

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute \f$ \int_{fb} \nabla (u) \cdot \nu_{fb} v \f$ where \p fb
 *         is a boundary face (Co+St algorithm)
 *
 * \param[in]       fb        index of the boundary face on the local numbering
 * \param[in]       beta      value of coefficient in front of the STAB part
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       kappa_f   diffusion property against face vector for all
 *                            faces
 * \param[in, out]  ntrgrd    pointer to a local matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_normal_flux_reco(short int                  fb,
                  const double               beta,
                  const cs_cell_mesh_t      *cm,
                  const cs_real_3_t         *kappa_f,
                  cs_sdm_t                  *ntrgrd)
{
  /* Sanity check */
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_HFQ));
  assert(cm->f_sgn[fb] == 1);  /* +1 because it's a boundary face */

  const short int  nfc = cm->n_fc;
  const cs_quant_t  pfbq = cm->face[fb];
  const cs_nvec3_t  debq = cm->dedge[fb];

  /* |fb|^2 * nu_{fb}^T.kappa.nu_{fb} */
  const cs_real_t  fb_k_fb = pfbq.meas * _dp3(kappa_f[fb], pfbq.unitv);
  const cs_real_t  beta_fbkfb_o_pfc = beta * fb_k_fb / cm->pfc[fb];
  const cs_real_t  ov_vol = 1./cm->vol_c;

  cs_real_t  *ntrgrd_fb = ntrgrd->val + fb * (nfc + 1);
  cs_real_t  row_sum = 0.0;
  for (short int f = 0; f < nfc; f++) {

    const cs_real_t  if_ov = ov_vol * cm->f_sgn[f];
    const cs_real_t  f_k_fb = pfbq.meas * _dp3(kappa_f[f], pfbq.unitv);
    const cs_quant_t  pfq = cm->face[f];

    cs_real_t  stab = -pfq.meas*debq.meas * _dp3(debq.unitv, pfq.unitv);
    if (f == fb) stab += cm->vol_c;

    const cs_real_t  int_gradf_dot_f = if_ov *
      ( f_k_fb                        /* Cons */
        + beta_fbkfb_o_pfc * stab);   /* Stab */
    ntrgrd_fb[f] -= int_gradf_dot_f;  /* Minus because -du/dn */
    row_sum      += int_gradf_dot_f;

  } /* Loop on f */

  /* Cell column */
  ntrgrd_fb[nfc] += row_sum;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the members of the cs_cdofb_navsto_builder_t structure
 *
 * \param[in]      t_eval     time at which one evaluates the pressure BC
 * \param[in]      nsp        set of parameters to define the NavSto system
 * \param[in]      cm         cellwise view of the mesh
 * \param[in]      csys       cellwise view of the algebraic system
 * \param[in]      pr_bc      set of definitions for the presuure BCs
 * \param[in]      bf_type    type of boundaries for all boundary faces
 * \param[in, out] nsb        builder to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_define_builder(cs_real_t                    t_eval,
                               const cs_navsto_param_t     *nsp,
                               const cs_cell_mesh_t        *cm,
                               const cs_cell_sys_t         *csys,
                               const cs_cdo_bc_face_t      *pr_bc,
                               const cs_boundary_type_t    *bf_type,
                               cs_cdofb_navsto_builder_t   *nsb)
{
  assert(cm != NULL && csys != NULL && nsp != NULL); /* sanity checks */

  const short int n_fc = cm->n_fc;

  /* Build the divergence operator:
   *        D(\hat{u}) = \frac{1}{|c|} \sum_{f_c} \iota_{fc} u_f.f
   * But in the linear system what appears is after integration
   *        [[ -div(u), q ]]_{P_c} = -|c| div(u)_c q_c
   * Thus, the volume in the divergence drops
   */

  for (short int f = 0; f < n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_real_t  sgn_f = -cm->f_sgn[f] * pfq.meas;

    cs_real_t  *_div_f = nsb->div_op + 3*f;
    _div_f[0] = sgn_f * pfq.unitv[0];
    _div_f[1] = sgn_f * pfq.unitv[1];
    _div_f[2] = sgn_f * pfq.unitv[2];

  } /* Loop on cell faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_NAVSTO_DBG > 2
  if (cs_dbg_cw_test(NULL, cm, csys)) {
    const cs_real_t  sovc = -1. / cm->vol_c;
#   pragma omp critical
    {
      cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
      for (short int f = 0; f < n_fc; f++)
        cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e, %- .4e, %- .4e\n",
                      f, nsb->div_op[3*f]*sovc, nsb->div_op[3*f+1]*sovc,
                      nsb->div_op[3*f+2]*sovc);
    } /* Critical section */
  }
#endif

  /* Build local arrays related to the boundary conditions */
  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering and the boundary face id in
       the mesh numbering */
    const short int  f = csys->_f_ids[i];
    const cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;

    /* Set the type of boundary */
    nsb->bf_type[i] = bf_type[bf_id];

    /* Set the pressure BC if required */
    if (nsb->bf_type[i] & CS_BOUNDARY_IMPOSED_P) {

      assert(nsb->bf_type[i] & (CS_BOUNDARY_INLET | CS_BOUNDARY_OUTLET));

      /* Add a Dirichlet for the pressure field */
      const short int  def_id = pr_bc->def_ids[bf_id];
      const cs_xdef_t  *def = nsp->pressure_bc_defs[def_id];
      assert(pr_bc != NULL);

      switch(def->type) {
      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;
          nsb->pressure_bc_val[i] = constant_val[0];
        }
        break;

      case CS_XDEF_BY_ARRAY:
        {
          cs_xdef_array_input_t  *a_in = (cs_xdef_array_input_t *)def->input;
          assert(a_in->stride == 1);
          assert(cs_flag_test(a_in->loc, cs_flag_primal_face));
          nsb->pressure_bc_val[i] = a_in->values[bf_id];
        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        switch(nsp->dof_reduction_mode) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_cw_eval_at_xyz_by_analytic(cm, 1, cm->face[f].center,
                                             t_eval,
                                             def->input,
                                             nsb->pressure_bc_val + i);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_cw_eval_scalar_face_avg_by_analytic(cm, f, t_eval,
                                                      def->input,
                                                      def->qtype,
                                                      nsb->pressure_bc_val + i);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Invalid type of reduction.\n"
                      " Stop computing the Dirichlet value.\n"), __func__);

        } /* switch on reduction */
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Invalid type of definition.\n"
                    " Stop computing the Dirichlet value.\n"), __func__);
        break;

      } /* def->type */

    }
    else
      nsb->pressure_bc_val[i] = 0.;

  } /* Loop on boundary faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence of a cell using the \ref cs_cdo_quantities_t
 *         structure
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdofb_navsto_cell_divergence(const cs_lnum_t               c_id,
                                const cs_cdo_quantities_t    *quant,
                                const cs_adjacency_t         *c2f,
                                const cs_real_t              *f_vals)
{
  cs_real_t  div = 0.0;

  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

    const cs_lnum_t  f_id = c2f->ids[f];
    const cs_real_t  *_val = f_vals + 3*f_id;

    if (f_id < quant->n_i_faces)
      div += c2f->sgn[f]*cs_math_3_dot_product(_val,
                                               quant->i_face_normal + 3*f_id);
    else {

      const cs_lnum_t  bf_id = f_id - quant->n_i_faces;

      div += c2f->sgn[f]* cs_math_3_dot_product(_val,
                                                quant->b_face_normal + 3*bf_id);

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

  /* We should ensure that the mean value of the pressure is zero. Thus we
   * compute it and subtract it from every value.
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
 * \brief  Initialize the pressure values when the pressure is defined at
 *         faces
 *
 * \param[in]       nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]       connect   pointer to a \ref cs_cdo_connect_t structure
 * \param[in]       ts        pointer to a \ref cs_time_step_t structure
 * \param[in, out]  pr_f      pointer to the pressure values at faces
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_init_face_pressure(const cs_navsto_param_t     *nsp,
                                   const cs_cdo_connect_t      *connect,
                                   const cs_time_step_t        *ts,
                                   cs_real_t                   *pr_f)
{
  /* Sanity checks */
  assert(nsp != NULL && pr_f != NULL);

  /* Initial conditions for the pressure */
  if (nsp->n_pressure_ic_defs == 0)
    return; /* Nothing to do */

  assert(nsp->pressure_ic_defs != NULL);

  cs_lnum_t  *def2f_ids = (cs_lnum_t *)cs_equation_get_tmpbuf();
  cs_lnum_t  *def2f_idx = NULL;
  BFT_MALLOC(def2f_idx, nsp->n_pressure_ic_defs + 1, cs_lnum_t);

  cs_equation_sync_definitions_at_faces(connect,
                                        nsp->n_pressure_ic_defs,
                                        nsp->pressure_ic_defs,
                                        def2f_idx,
                                        def2f_ids);

  const cs_real_t  t_cur = ts->t_cur;

  for (int def_id = 0; def_id < nsp->n_pressure_ic_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    cs_xdef_t  *def = nsp->pressure_ic_defs[def_id];
    const cs_lnum_t  n_f_selected = def2f_idx[def_id+1] - def2f_idx[def_id];
    const cs_lnum_t  *selected_lst = def2f_ids + def2f_idx[def_id];

    /* Initialize face-based array */
    switch (def->type) {

      /* Evaluating the integrals: the averages will be taken care of at the
       * end when ensuring zero-mean valuedness */
    case CS_XDEF_BY_VALUE:
      cs_evaluate_potential_at_faces_by_value(def,
                                              n_f_selected,
                                              selected_lst,
                                              pr_f);
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;

        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          cs_evaluate_potential_at_faces_by_analytic(def,
                                                     t_cur,
                                                     n_f_selected,
                                                     selected_lst,
                                                     pr_f);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_average_on_faces_by_analytic(def,
                                                   t_cur,
                                                   n_f_selected,
                                                   selected_lst,
                                                   pr_f);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    _(" %s: Incompatible reduction for the pressure field\n"),
                    __func__);

        }  /* Switch on possible reduction types */

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Incompatible way to initialize the pressure field.\n"),
                __func__);
      break;

    }  /* Switch on possible type of definition */

  }  /* Loop on definitions */

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

  cs_real_t  intgr = cs_weighted_sum(n_cells, quant->cell_vol, values);

  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &intgr);

  assert(quant->vol_tot > 0.);
  const cs_real_t  g_avg = intgr / quant->vol_tot;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    values[c_id] = values[c_id] - g_avg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform extra-operation related to Fb schemes when solving
 *         Navier-Stokes.
 *         - Compute the mass flux accross the boundaries.
 *
 * \param[in]  nsp        pointer to a \ref cs_navsto_param_t struct.
 * \param[in]  quant      pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]  connect    pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  adv_field  pointer to a \ref cs_adv_field_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_extra_op(const cs_navsto_param_t     *nsp,
                         const cs_cdo_quantities_t   *quant,
                         const cs_cdo_connect_t      *connect,
                         const cs_adv_field_t        *adv_field)
{
  CS_UNUSED(connect);

  const cs_boundary_t  *boundaries = nsp->boundaries;

  /* Retrieve the boundary velocity flux (mass flux) */
  cs_field_t  *nflx
    = cs_advection_field_get_field(adv_field, CS_MESH_LOCATION_BOUNDARY_FACES);

  /* 1. Compute for each boundary the integrated flux */
  bool  *belong_to_default = NULL;
  BFT_MALLOC(belong_to_default, quant->n_b_faces, bool);
# pragma omp parallel for if  (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_b_faces; i++)
    belong_to_default[i] = true;

  cs_real_t  *boundary_fluxes = NULL;
  BFT_MALLOC(boundary_fluxes, boundaries->n_boundaries + 1, cs_real_t);
  memset(boundary_fluxes, 0, (boundaries->n_boundaries + 1)*sizeof(cs_real_t));

  for (int b_id = 0; b_id < boundaries->n_boundaries; b_id++) {

    const cs_zone_t  *z = cs_boundary_zone_by_id(boundaries->zone_ids[b_id]);

    for (cs_lnum_t i = 0; i < z->n_elts; i++) {
      const cs_lnum_t  bf_id = z->elt_ids[i];
      belong_to_default[bf_id] = false;
      boundary_fluxes[b_id] += nflx->val[bf_id];
    }

  } /* Loop on domain boundaries */

  /* Update the flux through the default boundary */
  for (cs_lnum_t i = 0; i < quant->n_b_faces; i++) {
    if (belong_to_default[i])
      boundary_fluxes[boundaries->n_boundaries] += nflx->val[i];
  }

  /* Parallel synchronization if needed */
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(boundaries->n_boundaries + 1, CS_REAL_TYPE, boundary_fluxes);

  /* Output result */
  cs_log_printf(CS_LOG_DEFAULT,
                "--- Balance of the mass flux across the boundaries:\n");

  for (int b_id = 0; b_id < boundaries->n_boundaries; b_id++) {

    const cs_zone_t  *z = cs_boundary_zone_by_id(boundaries->zone_ids[b_id]);

    char descr[33];
    cs_boundary_get_type_descr(boundaries, boundaries->types[b_id], 33, descr);

    cs_log_printf(CS_LOG_DEFAULT, "-b- %-22s |%-32s |% -8.6e\n",
                  descr, z->name, boundary_fluxes[b_id]);

  } /* Loop on boundaries */

  /* Default boundary */
  switch (boundaries->default_type) {
  case CS_BOUNDARY_SYMMETRY:
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-22s |%-32s |% -8.6e\n",
                  "symmetry", "default boundary",
                  boundary_fluxes[boundaries->n_boundaries]);
    break;
  case CS_BOUNDARY_WALL:
    cs_log_printf(CS_LOG_DEFAULT, "-b- %-22s |%-32s |% -8.6e\n",
                  "wall", "default boundary",
                  boundary_fluxes[boundaries->n_boundaries]);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid type of default boundary.\n"
                " A valid choice is either \"CS_BOUNDARY_WALL\" or"
                " \"CS_BOUNDARY_SYMMETRY\"."), __func__);
  }

  /* Free temporary buffers */
  BFT_FREE(belong_to_default);
  BFT_FREE(boundary_fluxes);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         an algebraic technique.
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_alge(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(cm);

  double  *x_dir = cb->values;
  double  *ax_dir = cb->values + 3;
  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;
  assert(bd != NULL);
  assert(bd->n_row_blocks == cm->n_fc || bd->n_row_blocks == cm->n_fc + 1);

  /* Build x_dir */
  bool  is_non_homogeneous = true;

  memset(cb->values, 0, 6*sizeof(double));

  for (int k = 0; k < 3; k++) {
    if (csys->dof_flag[3*f+k] & CS_CDO_BC_DIRICHLET) {
      x_dir[k] = csys->dir_values[3*f+k];
      is_non_homogeneous = true;
    }
  }

  if (is_non_homogeneous) {

    for (int bi = 0; bi < bd->n_row_blocks; bi++) {

      if (bi == f)
        continue;

      cs_real_t  *_rhs = csys->rhs + 3*bi;
      cs_sdm_t  *mIF = cs_sdm_get_block(m, bi, f);

      cs_sdm_square_matvec(mIF, x_dir, ax_dir);
      for (int k = 0; k < 3; k++) _rhs[k] -= ax_dir[k];

    }

  } /* Non-homogeneous Dirichlet BC */

  /* Set RHS to the Dirichlet value for the related face */
  for (int k = 0; k < 3; k++)
    csys->rhs[3*f+k] = x_dir[k];

  /* Second pass: Replace the Dirichlet block by a diagonal block and fill with
   * zero the remaining row and column */
  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    if (bi != f) {

      /* Reset block (I,F) which is a 3x3 block */
      cs_sdm_t  *mIF = cs_sdm_get_block(m, bi, f);
      memset(mIF->val, 0, 9*sizeof(double));

    }
    else { /* bi == f */

      /* Reset block (I==F,J) which is a 3x3 block */
      for (int bj = 0; bj < bd->n_col_blocks; bj++) {
        cs_sdm_t  *mFJ = cs_sdm_get_block(m, f, bj);
        memset(mFJ->val, 0, 9*sizeof(double));
      }

      cs_sdm_t  *mFF = cs_sdm_get_block(m, f, f);
      assert((mFF->n_cols == 3) && (mFF->n_rows == 3));
      for (int k = 0; k < 3; k++)
        mFF->val[4*k] = 1; /* 4 == mFF->n_rows + 1 */

    }

  } /* Block bi */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a penalization technique (with a large coefficient).
 *         One assumes that static condensation has been performed and that
 *         the velocity-block has size 3*n_fc
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_pena(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  CS_UNUSED(cb);
  CS_UNUSED(cm);
  assert(cm != NULL && csys != NULL);

  cs_sdm_t  *m = csys->mat;
  assert(m->block_desc != NULL);
  assert(m->block_desc->n_row_blocks == cm->n_fc);

  const cs_flag_t  *_flag = csys->dof_flag + 3*f;
  const cs_real_t  *_dir_val = csys->dir_values + 3*f;

  bool  is_non_homogeneous = true;
  for (int k = 0; k < 3; k++) {
    if (_flag[k] & CS_CDO_BC_DIRICHLET)
      is_non_homogeneous = true;
  }

  /* Penalize diagonal entry (and its rhs if needed) */
  cs_sdm_t  *mFF = cs_sdm_get_block(m, f, f);
  assert((mFF->n_rows == 3) &&  (mFF->n_cols == 3));

  if (is_non_homogeneous) {

    cs_real_t  *_rhs = csys->rhs + 3*f;
    for (int k = 0; k < 3; k++) {
      mFF->val[4*k] += eqp->strong_pena_bc_coeff; /* 4 == mFF->n_rows + 1 */
      _rhs[k] += _dir_val[k] * eqp->strong_pena_bc_coeff;
    }

  }
  else {

    for (int k = 0; k < 3; k++)
      mFF->val[4*k] += eqp->strong_pena_bc_coeff; /* 4 == mFF->n_rows + 1 */

  } /* Homogeneous BC */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_weak(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  const cs_param_hodge_t  hodgep = eqp->diffusion_hodge;

  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);
  assert(cs_equation_param_has_diffusion(eqp));
  assert(hodgep.is_iso == true);

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */
  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int i = 0; i < cm->n_fc; i++) {
    const cs_real_t  coef = cm->face[i].meas*cb->dpty_val;
    for (short int k = 0; k < 3; k++)
      kappa_f[i][k] = coef * cm->face[i].unitv[k];
  }

  /* 1) Build the bc_op matrix (scalar-valued version) */

  /* Initialize the matrix related to the flux reconstruction operator */
  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix */
  _normal_flux_reco(f, hodgep.coef, cm, (const cs_real_t (*)[3])kappa_f, bc_op);

  /* 2) Update the bc_op matrix and the RHS with the Dirichlet values. */

  /* coeff * \meas{f} / h_f  */
  const cs_real_t pcoef = eqp->weak_pena_bc_coeff * sqrt(cm->face[f].meas);

  bc_op->val[f*(n_dofs + 1)] += pcoef; /* Diagonal term */

  for (short int k = 0; k < 3; k++)
    csys->rhs[3*f + k] += pcoef * csys->dir_values[3*f + k];

  /* 3) Update the local system matrix */
  for (int bi = 0; bi < n_dofs; bi++) { /* n_(scalar)_dofs == n_blocks */
    for (int bj = 0; bj < n_dofs; bj++) {

      /* Retrieve the 3x3 matrix */
      cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
      assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

      const cs_real_t  _val = bc_op->val[n_dofs*bi + bj];
      /* Update diagonal terms only */
      bij->val[0] += _val;
      bij->val[4] += _val;
      bij->val[8] += _val;

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         a weak penalization technique (symmetrized Nitsche).
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has size 3*(n_fc + 1)
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_wsym(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  const cs_param_hodge_t  hodgep = eqp->diffusion_hodge;

  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);
  assert(cs_equation_param_has_diffusion(eqp));
  assert(hodgep.is_iso == true);

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */
  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int i = 0; i < cm->n_fc; i++) {
    const cs_real_t  coef = cm->face[i].meas*cb->dpty_val;
    for (short int k = 0; k < 3; k++)
      kappa_f[i][k] = coef * cm->face[i].unitv[k];
  }

  /* 1) Build the bc_op matrix (scalar-valued version) */

  /* Initialize the matrix related to the flux reconstruction operator */
  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc, *bc_op_t = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix */
  _normal_flux_reco(f, hodgep.coef, cm, (const cs_real_t (*)[3])kappa_f, bc_op);

  /* 2) Update the bc_op matrixand the RHS with the Dirichlet values */

  /* Update bc_op = bc_op + transp and transp = transpose(bc_op)
     cb->loc plays the role of the flux operator */
  cs_sdm_square_add_transpose(bc_op, bc_op_t);

  /* Update the RHS with bc_op_t*dir_face */
  for (short int k = 0; k < 3; k++) {

    /* Only this value is not zero (simplify the matvec operation) */
    const cs_real_t  dir_f = csys->dir_values[3*f+k];
    for (short int i = 0; i < n_dofs; i++)
      csys->rhs[3*i+k] += bc_op_t->val[i*n_dofs+f] * dir_f;

  } /* Loop on components */

  /* 3) Update the bc_op matrix and the RHS with the penalization */

  /* coeff * \meas{f} / h_f  */
  const cs_real_t pcoef = eqp->weak_pena_bc_coeff * sqrt(cm->face[f].meas);

  bc_op->val[f*(n_dofs + 1)] += pcoef; /* Diagonal term */

  for (short int k = 0; k < 3; k++)
    csys->rhs[3*f + k] += pcoef * csys->dir_values[3*f + k];

  /* 4) Update the local system matrix */
  for (int bi = 0; bi < n_dofs; bi++) { /* n_(scalar)_dofs == n_blocks */
    for (int bj = 0; bj < n_dofs; bj++) {

      /* Retrieve the 3x3 matrix */
      cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
      assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

      const cs_real_t  _val = bc_op->val[n_dofs*bi + bj];
      /* Update diagonal terms only */
      bij->val[0] += _val;
      bij->val[4] += _val;
      bij->val[8] += _val;

    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a symmetric boundary (treated as a sliding BCs on
 *         the three velocity components.
 *         A weak penalization technique (symmetrized Nitsche) is used.
 *         One assumes that static condensation has not been performed yet and
 *         that the velocity-block has (n_fc + 1) blocks of size 3x3.
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_symmetry(short int                       f,
                  const cs_equation_param_t      *eqp,
                  const cs_cell_mesh_t           *cm,
                  cs_cell_builder_t              *cb,
                  cs_cell_sys_t                  *csys)
{
  const cs_param_hodge_t  hodgep = eqp->diffusion_hodge;

  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);
  assert(hodgep.is_iso == true); /* if not the case something else TODO ? */
  assert(cs_equation_param_has_diffusion(eqp));

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */
  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int i = 0; i < cm->n_fc; i++) {
    const cs_real_t  coef = cm->face[i].meas*cb->dpty_val;
    for (short int k = 0; k < 3; k++)
      kappa_f[i][k] = coef * cm->face[i].unitv[k];
  }

  /* Initialize the matrix related this flux reconstruction operator */
  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->hdg;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix */
  _normal_flux_reco(f, hodgep.coef, cm, (const cs_real_t (*)[3])kappa_f, bc_op);

  /* 2) Update the bc_op matrix and nothing done to the RHS since a sliding
     means homogeneous Dirichlet values on the normal component and hommogeneous
     Neumann on the tangential flux */

  const cs_quant_t  pfq = cm->face[f];
  const cs_real_t  *nf = pfq.unitv;
  const cs_real_33_t  nf_nf = { {nf[0]*nf[0], nf[0]*nf[1], nf[0]*nf[2]},
                                {nf[1]*nf[0], nf[1]*nf[1], nf[1]*nf[2]},
                                {nf[2]*nf[0], nf[2]*nf[1], nf[2]*nf[2]} };

  /* chi * \meas{f} / h_f  */
  const cs_real_t  pcoef = eqp->weak_pena_bc_coeff * sqrt(pfq.meas);

  /* Handle the diagonal block: Retrieve the 3x3 matrix */
  cs_sdm_t  *bFF = cs_sdm_get_block(csys->mat, f, f);
  assert(bFF->n_rows == bFF->n_cols && bFF->n_rows == 3);

  const cs_real_t  _val = pcoef + 2*bc_op->val[f*(n_dofs+1)];
  for (short int k = 0; k < 3; k++) {
    bFF->val[3*k  ] += nf_nf[0][k] * _val;
    bFF->val[3*k+1] += nf_nf[1][k] * _val;
    bFF->val[3*k+2] += nf_nf[2][k] * _val;
  }

  for (short int xj = 0; xj < n_dofs; xj++) {

    if (xj == f)
      continue;

    /* It should be done both for face- and cell-defined DoFs */
    /* Retrieve the 3x3 matrix */
    cs_sdm_t  *bFJ = cs_sdm_get_block(csys->mat, f, xj);
    assert(bFJ->n_rows == bFJ->n_cols && bFJ->n_rows == 3);
    cs_sdm_t  *bJF = cs_sdm_get_block(csys->mat, xj, f);
    assert(bJF->n_rows == bJF->n_cols && bJF->n_rows == 3);

    const cs_real_t  op_fj = bc_op->val[n_dofs*f  + xj];
    const cs_real_t  op_jf = bc_op->val[n_dofs*xj + f];
    const cs_real_t  _val_fj_jf = op_fj + op_jf;

    for (int k = 0; k < 3; k++) {

      bFJ->val[3*k  ] += nf_nf[0][k] * _val_fj_jf;
      bFJ->val[3*k+1] += nf_nf[1][k] * _val_fj_jf;
      bFJ->val[3*k+2] += nf_nf[2][k] * _val_fj_jf;

      /* nf_nf is symmetric */
      bJF->val[3*k  ] += nf_nf[0][k] * _val_fj_jf;
      bJF->val[3*k+1] += nf_nf[1][k] * _val_fj_jf;
      bJF->val[3*k+2] += nf_nf[2][k] * _val_fj_jf;

    }

  } /* Loop on xj */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a wall BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_fixed_wall(short int                       f,
                    const cs_equation_param_t      *eqp,
                    const cs_cell_mesh_t           *cm,
                    cs_cell_builder_t              *cb,
                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(cb);
  assert(cm != NULL && csys != NULL);  /* Sanity checks */

  const cs_quant_t  pfq = cm->face[f];
  const cs_real_t  *ni = pfq.unitv;
  const cs_real_t  ni_ni[9] = { ni[0]*ni[0], ni[0]*ni[1], ni[0]*ni[2],
                                ni[1]*ni[0], ni[1]*ni[1], ni[1]*ni[2],
                                ni[2]*ni[0], ni[2]*ni[1], ni[2]*ni[2]};

  /* chi * \meas{f} / h_f  */
  const cs_real_t  pcoef = eqp->weak_pena_bc_coeff * sqrt(pfq.meas);

  cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, f, f);
  assert(bii->n_rows == bii->n_cols && bii->n_rows == 3);

  for (short int k = 0; k < 9; k++)
    bii->val[k] += pcoef * ni_ni[k];
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
