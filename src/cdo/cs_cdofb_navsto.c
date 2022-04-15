/*============================================================================
 * Functions shared among all face-based schemes for the discretization of the
 * Navier-Stokes system
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
#include "cs_cdo_blas.h"
#include "cs_cdo_toolbox.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_reco.h"
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
 * \brief Shared functions among all face-based schemes for building and
 *        solving Stokes and Navier-Stokes problem
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

static cs_cdofb_navsto_boussinesq_type_t  cs_cdofb_navsto_boussinesq_type =
  CS_CDOFB_NAVSTO_BOUSSINESQ_FACE_DOF;

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
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_PFC));
  assert(cm->f_sgn[fb] == 1);  /* +1 because it's a boundary face */

  const short int  nfc = cm->n_fc;
  const double  inv_volc = 1./cm->vol_c;
  const cs_quant_t  pfbq = cm->face[fb];
  const cs_nvec3_t  debq = cm->dedge[fb];

  /* |fb|^2 * nu_{fb}^T.kappa.nu_{fb} */

  const double  stab_scaling =
    beta * pfbq.meas * _dp3(kappa_f[fb], pfbq.unitv) / cm->pvol_f[fb];

  /* Get the 'fb' row */

  cs_real_t  *ntrgrd_fb = ntrgrd->val + fb * (nfc + 1);
  double  row_sum = 0.0;

  for (short int f = 0; f < nfc; f++) {

    const cs_quant_t  pfq = cm->face[f];

    /* consistent part */

    const double  consist_scaling = cm->f_sgn[f] * pfq.meas * inv_volc;
    const double  consist_part = consist_scaling * _dp3(kappa_f[fb], pfq.unitv);

    /* stabilization part */

    double  stab_part = -consist_scaling*debq.meas*_dp3(debq.unitv, pfq.unitv);
    if (f == fb) stab_part += 1;
    stab_part *= stab_scaling;

    const double  fb_f_part = consist_part + stab_part;

    ntrgrd_fb[f] -= fb_f_part;  /* Minus because -du/dn */
    row_sum      += fb_f_part;

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
 * \brief  Set the way to compute the Boussinesq approximation
 *
 * \param[in] type     type of algorithm to use
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_set_boussinesq_algo(cs_cdofb_navsto_boussinesq_type_t   type)
{
  cs_cdofb_navsto_boussinesq_type = type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and allocate a local NavSto builder when Fb schemes are used
 *
 * \param[in] nsp         set of parameters to define the NavSto system
 * \param[in] connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a cs_cdofb_navsto_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdofb_navsto_builder_t
cs_cdofb_navsto_create_builder(const cs_navsto_param_t  *nsp,
                               const cs_cdo_connect_t   *connect)
{
  cs_cdofb_navsto_builder_t  nsb = {.rho_c = 1.,
                                    .div_op = NULL,
                                    .bf_type = NULL,
                                    .pressure_bc_val = NULL};

  assert(nsp != NULL);
  nsb.rho_c = nsp->mass_density->ref_value;

  if (connect == NULL)
    return nsb;

  BFT_MALLOC(nsb.div_op, 3*connect->n_max_fbyc, cs_real_t);
  BFT_MALLOC(nsb.bf_type, connect->n_max_fbyc, cs_boundary_type_t);
  BFT_MALLOC(nsb.pressure_bc_val, connect->n_max_fbyc, cs_real_t);

  return nsb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy the given cs_cdofb_navsto_builder_t structure
 *
 * \param[in, out] nsb   pointer to the cs_cdofb_navsto_builder_t to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_free_builder(cs_cdofb_navsto_builder_t   *nsb)
{
  if (nsb != NULL) {
    BFT_FREE(nsb->div_op);
    BFT_FREE(nsb->bf_type);
    BFT_FREE(nsb->pressure_bc_val);
  }
}

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

  /* Update the value of the mass density for the current cell if needed */
  /* TODO: Case of a uniform but not constant in time */

  if (!cs_property_is_uniform(nsp->mass_density))
    nsb->rho_c = cs_property_value_in_cell(cm, nsp->mass_density, t_eval);

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
      for (short int f = 0; f < cm->n_fc; f++)
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
          const cs_real_t  *constant_val = (cs_real_t *)def->context;
          nsb->pressure_bc_val[i] = constant_val[0];
        }
        break;

      case CS_XDEF_BY_ARRAY:
        {
          cs_xdef_array_context_t  *c = def->context;
          assert(c->stride == 1);
          assert(cs_flag_test(c->loc, cs_flag_primal_face));
          nsb->pressure_bc_val[i] = c->values[bf_id];
        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        switch(nsp->dof_reduction_mode) {

        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_cw_eval_at_xyz_by_analytic(cm, 1, cm->face[f].center,
                                             t_eval,
                                             def->context,
                                             nsb->pressure_bc_val + i);
          break;

        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_cw_eval_scalar_face_avg_by_analytic(cm, f, t_eval,
                                                      def->context,
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
 * \brief  Compute the mass flux playing the role of the advection field in
 *         the Navier-Stokes equations
 *         One considers the mass flux across primal faces which relies on the
 *         velocity vector defined on each face.
 *
 * \param[in]      nsp         set of parameters to define the NavSto system
 * \param[in]      quant       set of additional geometrical quantities
 * \param[in]      face_vel    velocity vectors for each face
 * \param[in, out] mass_flux   array of mass flux values to update (allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_mass_flux(const cs_navsto_param_t     *nsp,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *face_vel,
                          cs_real_t                   *mass_flux)
{
  if (mass_flux == NULL)
    return;

  assert(face_vel != NULL);
  assert(nsp->space_scheme == CS_SPACE_SCHEME_CDOFB);
  assert(cs_property_is_uniform(nsp->mass_density));
  assert(nsp->mass_density->n_definitions == 1);

  const cs_real_t  rho_val = nsp->mass_density->ref_value;

  /* Define the mass flux */

# pragma omp parallel for if (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++) {

    const cs_real_t  *fq = cs_quant_get_face_vector_area(f_id, quant);
    mass_flux[f_id] = rho_val*cs_math_3_dot_product(face_vel + 3*f_id, fq);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence in a cell of a vector-valued array defined at
 *         faces (values are defined both at interior and border faces).
 *         Variant based on the usage of \ref cs_cdo_quantities_t structure.
 *
 * \param[in]     c_id         cell id
 * \param[in]     quant        pointer to a \ref cs_cdo_quantities_t
 * \param[in]     c2f          pointer to cell-to-face \ref cs_adjacency_t
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

double
cs_cdofb_navsto_cell_divergence(const cs_lnum_t               c_id,
                                const cs_cdo_quantities_t    *quant,
                                const cs_adjacency_t         *c2f,
                                const cs_real_t              *f_vals)
{
  const cs_lnum_t  thd = 3 * quant->n_i_faces;

  double  div = 0.0;
  for (cs_lnum_t f = c2f->idx[c_id]; f < c2f->idx[c_id+1]; f++) {

    const cs_lnum_t  shift = 3*c2f->ids[f];
    const cs_real_t  *_val = f_vals + shift;
    const cs_real_t  *fnorm = (shift < thd) ?
      quant->i_face_normal + shift : quant->b_face_normal + (shift - thd);

    div += c2f->sgn[f]*cs_math_3_dot_product(_val, fnorm);

  } /* Loop on cell faces */

  div /= quant->cell_vol[c_id];

  return div;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an estimation of the pressure at faces
 *
 * \param[in]       mesh       pointer to a cs_mesh_t structure
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]       time_step  pointer to a cs_time_step_t structure
 * \param[in]       nsp        pointer to a \ref cs_navsto_param_t struct.
 * \param[in]       p_cell     value of the pressure inside each cell
 * \param[in, out]  p_face     value of the pressure at each face
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_compute_face_pressure(const cs_mesh_t             *mesh,
                                      const cs_cdo_connect_t      *connect,
                                      const cs_cdo_quantities_t   *quant,
                                      const cs_time_step_t        *ts,
                                      const cs_navsto_param_t     *nsp,
                                      const cs_real_t             *p_cell,
                                      cs_real_t                   *p_face)
{
  /* Interior faces */

  for (cs_lnum_t  i = 0; i < mesh->n_i_faces; i++) {

    cs_lnum_t  iid = mesh->i_face_cells[i][0];
    cs_lnum_t  jid = mesh->i_face_cells[i][1];

    p_face[i] = 0.5*(p_cell[iid] + p_cell[jid]);

  }

  /* Border faces
   * Use the knowledge from the BCs if possible, otherwise assume a homogeneous
   * Neumann (i.e. p_face = p_cell). */

  cs_real_t  *p_bface = p_face + mesh->n_i_faces;

  for (cs_lnum_t  i = 0; i < mesh->n_b_faces; i++)
    p_bface[i] = p_cell[mesh->b_face_cells[i]];

  for (int def_id = 0; def_id < nsp->n_pressure_bc_defs; def_id++) {

    cs_xdef_t  *pbc_def = nsp->pressure_bc_defs[def_id];
    const cs_zone_t  *z = cs_boundary_zone_by_id(pbc_def->z_id);

    assert(pbc_def->meta & CS_CDO_BC_DIRICHLET);

    switch(pbc_def->type) {

    case CS_XDEF_BY_VALUE:
      cs_xdef_eval_scalar_by_val(z->n_elts, z->elt_ids,
                                 false, /* dense output */
                                 mesh,
                                 connect,
                                 quant,
                                 ts->t_cur,
                                 pbc_def->context,
                                 p_bface);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Type of definition not handle.\n"
                  " Stop computing the pressure BC value.\n"), __func__);
      break;

    } /* def->type */

  } /* Loop on pressure BCs */
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
      /* When constant mean-value or the value at cell center is the same */
      cs_evaluate_potential_at_cells_by_value(def, values);
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        const cs_param_dof_reduction_t  red = nsp->dof_reduction_mode;

        switch (red) {
        case CS_PARAM_REDUCTION_DERHAM:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_density_by_analytic(dof_flag, def, t_cur, values);
          break;
        case CS_PARAM_REDUCTION_AVERAGE:
          cs_xdef_set_quadrature(def, nsp->qtype);
          cs_evaluate_average_on_cells_by_analytic(def, t_cur, values);
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

  cs_cdofb_navsto_rescale_pressure_to_ref(nsp, quant, values);
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
  CS_UNUSED(connect);
  assert(nsp != NULL && pr_f != NULL);

  /* Initial conditions for the pressure */

  if (nsp->n_pressure_ic_defs == 0)
    return; /* Nothing to do */

  assert(nsp->pressure_ic_defs != NULL);

  cs_lnum_t  *def2f_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
  cs_lnum_t  *def2f_idx = NULL;

  BFT_MALLOC(def2f_idx, nsp->n_pressure_ic_defs + 1, cs_lnum_t);

  cs_cdo_sync_vol_def_at_faces(nsp->n_pressure_ic_defs, nsp->pressure_ic_defs,
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

  BFT_FREE(def2f_idx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the pressure field in order to get a field with a mean-value
 *         equal to the reference value
 *
 * \param[in]       nsp       pointer to a cs_navsto_param_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  values    pressure field values
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_rescale_pressure_to_ref(const cs_navsto_param_t    *nsp,
                                        const cs_cdo_quantities_t  *quant,
                                        cs_real_t                   values[])
{
  const cs_lnum_t  n_cells = quant->n_cells;

  /*
   * The algorithm used for summing is l3superblock60, based on the article:
   * "Reducing Floating Point Error in Dot Product Using the Superblock Family
   * of Algorithms" by Anthony M. Castaldo, R. Clint Whaley, and Anthony
   * T. Chronopoulos, SIAM J. SCI. COMPUT., Vol. 31, No. 2, pp. 1156--1174
   * 2008 Society for Industrial and Applied Mathematics
   */

  cs_real_t  intgr = cs_weighted_sum(n_cells, quant->cell_vol, values);

  cs_parall_sum(1, CS_REAL_TYPE, &intgr);

  assert(quant->vol_tot > 0.);
  const cs_real_t  g_avg = intgr / quant->vol_tot;
  const cs_real_t  p_shift = nsp->reference_pressure - g_avg;

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    values[c_id] = values[c_id] + p_shift;
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
 *         Navier-Stokes. Computation of the following quantities according to
 *         post-processing flags beeing activated.
 *         - The mass flux accross the boundaries.
 *         - The global mass in the computational domain
 *         - The norm of the velocity divergence
 *         - the cellwise mass flux balance
 *         - the kinetic energy
 *         - the velocity gradient
 *         - the pressure gradient
 *         - the vorticity
 *         - the helicity
 *         - the enstrophy
 *         - the stream function
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t struct.
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      quant         pointer to a \ref cs_cdo_quantities_t struct.
 * \param[in]      connect       pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]      ts            pointer to a \ref cs_time_step_t struct.
 * \param[in,out]  time_plotter  pointer to a \ref cs_time_plot_t struct.
 * \param[in]      adv_field     pointer to a \ref cs_adv_field_t struct.
 * \param[in]      mass_flux     scalar-valued mass flux for each face
 * \param[in]      p_cell        scalar-valued pressure in each cell
 * \param[in]      u_cell        vector-valued velocity in each cell
 * \param[in]      u_face        vector-valued velocity on each face
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_extra_op(const cs_navsto_param_t     *nsp,
                         const cs_mesh_t             *mesh,
                         const cs_cdo_quantities_t   *quant,
                         const cs_cdo_connect_t      *connect,
                         const cs_time_step_t        *ts,
                         cs_time_plot_t              *time_plotter,
                         const cs_adv_field_t        *adv_field,
                         const cs_real_t             *mass_flux,
                         const cs_real_t             *p_cell,
                         const cs_real_t             *u_cell,
                         const cs_real_t             *u_face)
{
  CS_UNUSED(adv_field);

  const cs_boundary_t  *boundaries = nsp->boundaries;
  const cs_real_t  *bmass_flux = mass_flux + quant->n_i_faces;

  /* 1. Compute for each boundary the integrated mass flux to perform mass
   *    balance
   */

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
      boundary_fluxes[b_id] += bmass_flux[bf_id];
    }

  } /* Loop on domain boundaries */

  /* Update the flux through the default boundary */

  cs_lnum_t  default_case_count = 0;
  for (cs_lnum_t i = 0; i < quant->n_b_faces; i++) {
    if (belong_to_default[i]) {
      default_case_count += 1;
      boundary_fluxes[boundaries->n_boundaries] += bmass_flux[i];
    }
  }

  /* Parallel synchronization if needed */

  cs_parall_sum(boundaries->n_boundaries + 1, CS_REAL_TYPE, boundary_fluxes);
  cs_parall_counter_max(&default_case_count, 1);

  /* Output result */

  cs_log_printf(CS_LOG_DEFAULT,
                "\n- Balance of the mass flux across the boundaries:\n");

  char descr[32];
  for (int b_id = 0; b_id < boundaries->n_boundaries; b_id++) {

    const cs_zone_t  *z = cs_boundary_zone_by_id(boundaries->zone_ids[b_id]);

    cs_boundary_get_type_descr(boundaries, boundaries->types[b_id], 32, descr);

    cs_log_printf(CS_LOG_DEFAULT, "b %-32s | %-32s |% -8.6e\n",
                  descr, z->name, boundary_fluxes[b_id]);

  } /* Loop on boundaries */

  /* Default boundary (if something to do) */

  if (default_case_count > 0) {

    cs_boundary_get_type_descr(boundaries, boundaries->default_type, 32, descr);
    cs_log_printf(CS_LOG_DEFAULT, "b %-32s | %-32s |% -8.6e\n",
                  descr, "default boundary",
                  boundary_fluxes[boundaries->n_boundaries]);

  }

  /* Free temporary buffers */

  BFT_FREE(belong_to_default);
  BFT_FREE(boundary_fluxes);

  /* Predefined post-processing */
  /* ========================== */

  /* There are five values if all flags are activated for the monitoring plot */

  int  n_cols = 0;
  cs_real_t  col_vals[5] = {0, 0, 0, 0, 0};

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE) {

    double  div_norm2 = 0.;
    cs_field_t  *vel_div = cs_field_by_name("velocity_divergence");
    assert(vel_div != NULL);

    if (nsp->coupling != CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY)
      cs_field_current_to_previous(vel_div);

    /* Only the face velocity is used */

#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN) \
    reduction(+:div_norm2)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double  div_c = cs_cdofb_navsto_cell_divergence(c_id,
                                                      quant,
                                                      connect->c2f,
                                                      u_face);

      vel_div->val[c_id] = div_c;
      div_norm2 += quant->cell_vol[c_id] * div_c * div_c;

    } /* Loop on cells */

    cs_parall_sum(1, CS_DOUBLE, &div_norm2);
    col_vals[n_cols++] = sqrt(div_norm2);

  } /* Velocity divergence */

  if (nsp->post_flag & CS_NAVSTO_POST_MASS_DENSITY) {

    double  mass_integral = 0.;
    cs_field_t  *rho = cs_field_by_name("mass_density");
    assert(rho != NULL);

    cs_field_current_to_previous(rho);

#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN) \
    reduction(+:mass_integral)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double  boussi_coef = 1;

      for (int i = 0; i < nsp->n_boussinesq_terms; i++) {

        cs_navsto_param_boussinesq_t  *bp = nsp->boussinesq_param + i;
        boussi_coef += -bp->beta*(bp->var[c_id] - bp->var0);

      } /* Loop on Boussinesq terms */

      double  rho_c = nsp->mass_density->ref_value * boussi_coef;
      rho->val[c_id] = rho_c;
      mass_integral += quant->cell_vol[c_id] * rho_c;

    } /* Loop on cells */

    cs_parall_sum(1, CS_DOUBLE, &mass_integral);
    col_vals[n_cols++] = mass_integral;

  } /* Mass density */

  if (nsp->post_flag & CS_NAVSTO_POST_CELL_MASS_FLUX_BALANCE) {

    cs_field_t  *mf_balance = cs_field_by_name("mass_flux_balance");
    assert(mf_balance != NULL);

    cs_field_current_to_previous(mf_balance);

    const cs_adjacency_t  *c2f = connect->c2f;

#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      double  balance = 0.;
      for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {
        balance += c2f->sgn[j]*mass_flux[c2f->ids[j]];

      }
      mf_balance->val[c_id] = balance;

    } /* Loop on cells */

  } /* Cell mass flux balance */

  if (nsp->post_flag & CS_NAVSTO_POST_PRESSURE_GRADIENT) {

    cs_field_t  *pr_grd = cs_field_by_name("pressure_gradient");
    assert(pr_grd != NULL);

    cs_field_current_to_previous(pr_grd);

    /* Compute a face pressure */

    cs_real_t  *p_face = NULL;
    BFT_MALLOC(p_face, quant->n_faces, cs_real_t);

    cs_cdofb_navsto_compute_face_pressure(mesh,
                                          connect,
                                          quant,
                                          ts,
                                          nsp,
                                          p_cell,
                                          p_face);

#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      cs_reco_grad_cell_from_fb_dofs(c_id, connect, quant,
                                     p_cell, p_face,
                                     pr_grd->val + 3*c_id);

    BFT_FREE(p_face);

  } /* Pressure gradient */

  if (nsp->post_flag & CS_NAVSTO_POST_KINETIC_ENERGY) {

    double  k_integral = 0.;
    cs_field_t  *kinetic_energy = cs_field_by_name("kinetic_energy");
    assert(kinetic_energy != NULL);

    cs_field_current_to_previous(kinetic_energy);

    if (cs_property_is_uniform(nsp->mass_density)) {

      /* This can be any cell but one assumes that there is at least one cell by
         MPI rank */

      const cs_real_t  rho = cs_property_get_cell_value(0, /* cell_id */
                                                        ts->t_cur,
                                                        nsp->mass_density);

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)  \
      reduction(+:k_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        double kc = 0.5*rho*cs_math_3_square_norm(u_cell + 3*c_id);

        kinetic_energy->val[c_id] = kc;
        k_integral += quant->cell_vol[c_id]*kc;

      }

    }
    else { /* Mass density is not uniform in space */

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN) \
      reduction(+:k_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        cs_real_t  rho_c = cs_property_get_cell_value(c_id,
                                                      ts->t_cur,
                                                      nsp->mass_density);
        double  kc = 0.5*rho_c*cs_math_3_square_norm(u_cell + 3*c_id);

        kinetic_energy->val[c_id] = kc;
        k_integral += quant->cell_vol[c_id] * kc;

      }

    }

    cs_parall_sum(1, CS_DOUBLE, &k_integral); /* Sync. parallel computations */
    col_vals[n_cols++] = k_integral;

  } /* Kinetic energy */

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_GRADIENT) {

    cs_field_t  *velocity_gradient = cs_field_by_name("velocity_gradient");
    assert(velocity_gradient != NULL);

    cs_field_current_to_previous(velocity_gradient);

#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                        u_cell, u_face,
                                        velocity_gradient->val + 9*c_id);

  } /* Velocity gradient */

  cs_flag_t  mask_velgrd[3] = { CS_NAVSTO_POST_VORTICITY,
                                CS_NAVSTO_POST_HELICITY,
                                CS_NAVSTO_POST_ENSTROPHY };

  if (cs_flag_at_least(nsp->post_flag, 3, mask_velgrd)) {

    cs_field_t  *vorticity = cs_field_by_name("vorticity");
    assert(vorticity != NULL);
    cs_field_current_to_previous(vorticity);

    double  e_integral = 0.;
    cs_field_t  *enstrophy = cs_field_by_name_try("enstrophy");
    if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY)
      cs_field_current_to_previous(enstrophy);

    double  h_integral = 0.;
    cs_field_t  *helicity = cs_field_by_name_try("helicity");
    if (nsp->post_flag & CS_NAVSTO_POST_HELICITY)
      cs_field_current_to_previous(helicity);

    cs_field_t  *velocity_gradient = cs_field_by_name_try("velocity_gradient");

    if (velocity_gradient == NULL) {

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN) \
      reduction(+:e_integral) reduction(+:h_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        cs_real_t  grd_uc[9];

        /* Compute the velocity gradient */

        cs_reco_grad_33_cell_from_fb_dofs(c_id, connect, quant,
                                          u_cell, u_face, grd_uc);

        /* Compute the cell vorticity */

        cs_real_t  *w = vorticity->val + 3*c_id;
        w[0] = grd_uc[7] - grd_uc[5];
        w[1] = grd_uc[2] - grd_uc[6];
        w[2] = grd_uc[3] - grd_uc[1];

        if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY) {

          double  ec = cs_math_3_square_norm(w);
          enstrophy->val[c_id] = ec;
          e_integral += quant->cell_vol[c_id] * ec;

        }

        if (nsp->post_flag & CS_NAVSTO_POST_HELICITY) {

          double  hc = cs_math_3_dot_product(u_cell + 3*c_id, w);
          helicity->val[c_id] = hc;
          h_integral += quant->cell_vol[c_id] * hc;

        }

      } /* Loop on cells */

    }
    else {

#     pragma omp parallel for if (quant->n_cells > CS_THR_MIN)  \
      reduction(+:e_integral) reduction(+:h_integral)
      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        cs_real_t  *grd_uc = velocity_gradient->val + 9*c_id;

        /* Compute the cell vorticity */

        cs_real_t  *w = vorticity->val + 3*c_id;
        w[0] = grd_uc[7] - grd_uc[5];
        w[1] = grd_uc[2] - grd_uc[6];
        w[2] = grd_uc[3] - grd_uc[1];

        if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY) {

          double  ec = cs_math_3_square_norm(w);
          enstrophy->val[c_id] = ec;
          e_integral += quant->cell_vol[c_id] * ec;

        }

        if (nsp->post_flag & CS_NAVSTO_POST_HELICITY) {

          double  hc = cs_math_3_dot_product(u_cell + 3*c_id, w);
          helicity->val[c_id] = hc;
          h_integral += quant->cell_vol[c_id] * hc;

        }

      } /* Loop on cells */

    } /* velocity gradient has been computed previously */

    if (nsp->post_flag & CS_NAVSTO_POST_ENSTROPHY) {

      cs_parall_sum(1, CS_DOUBLE, &e_integral);
      col_vals[n_cols++] = e_integral;

    }

    if (nsp->post_flag & CS_NAVSTO_POST_HELICITY) {

      cs_parall_sum(1, CS_DOUBLE, &h_integral);
      col_vals[n_cols++] = h_integral;

    }

  } /* vorticity, helicity or enstrophy computations */

  if (cs_glob_rank_id < 1 && time_plotter != NULL)
    cs_time_plot_vals_write(time_plotter,
                            ts->nt_cur,
                            ts->t_cur,
                            n_cols,
                            col_vals);

  /* Stream function */
  /* --------------- */

  if (nsp->post_flag & CS_NAVSTO_POST_STREAM_FUNCTION) {

    cs_equation_t  *eq = cs_equation_by_name(CS_NAVSTO_STREAM_EQNAME);
    assert(eq != NULL);
    cs_equation_solve_steady_state(mesh, eq);

    cs_equation_param_t  *eqp = cs_equation_get_param(eq);
    if (eqp->n_bc_defs == 0) {

      /* Since this is an equation solved with only homogeneous Neumann BCs, one
       * substracts the mean value to get a unique solution */

      cs_real_t  mean_value;
      cs_equation_integrate_variable(connect, quant, eq, &mean_value);
      mean_value /= quant->vol_tot;

      cs_real_t  *psi_v = cs_equation_get_vertex_values(eq, false);
      for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
        psi_v[i] -= mean_value;

    } /* If homogeneous Neumann everywhere */

  } /* Computation of the stream function is requested */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account a Dirichlet BCs on the three velocity components.
 *         For instance for a velocity inlet boundary or a wall
 *         Handle the velocity-block in the global algebraic system in case of
 *         an algebraic technique.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_alge(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              const cs_property_data_t       *pty,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(cm);
  CS_UNUSED(pty);

  double  *x_dir = cb->values;
  double  *ax_dir = cb->values + 3;
  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;
  assert(bd != NULL);
  assert(bd->n_row_blocks == cm->n_fc || bd->n_row_blocks == cm->n_fc + 1);

  /* Build x_dir */

  bool  is_non_homogeneous = false; /* Assume homogeneous by default */

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
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       f         face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_pena(short int                       f,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              const cs_property_data_t       *pty,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  CS_UNUSED(cb);
  CS_UNUSED(cm);
  CS_UNUSED(pty);

  assert(csys != NULL);

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
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_weak(short int                       fb,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              const cs_property_data_t       *pty,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  assert(cm != NULL && cb != NULL && csys != NULL && pty != NULL);
  assert(pty->is_iso == true);

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int f = 0; f < cm->n_fc; f++) {
    const double  coef = cm->face[f].meas*pty->value;
    for (short int k = 0; k < 3; k++)
      kappa_f[f][k] = coef * cm->face[f].unitv[k];
  }

  /* 1) Build the bc_op matrix (scalar-valued version) */

  /* Initialize the matrix related to the flux reconstruction operator */

  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix */

  _normal_flux_reco(fb, eqp->diffusion_hodgep.coef, cm,
                    (const cs_real_t (*)[3])kappa_f,
                    bc_op);

  /* 2) Update the bc_op matrix and the RHS with the Dirichlet values. */

  /* coeff * \meas{f} / h_f  */

  const cs_real_t  pcoef = eqp->weak_pena_bc_coeff * sqrt(cm->face[fb].meas);

  bc_op->val[fb*(n_dofs + 1)] += pcoef; /* Diagonal term */

  for (short int k = 0; k < 3; k++)
    csys->rhs[3*fb + k] += pcoef * csys->dir_values[3*fb + k];

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
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_block_dirichlet_wsym(short int                       fb,
                              const cs_equation_param_t      *eqp,
                              const cs_cell_mesh_t           *cm,
                              const cs_property_data_t       *pty,
                              cs_cell_builder_t              *cb,
                              cs_cell_sys_t                  *csys)
{
  assert(cm != NULL && cb != NULL && csys != NULL && pty != NULL);
  assert(cs_equation_param_has_diffusion(eqp));
  assert(pty->is_iso == true);

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int f = 0; f < cm->n_fc; f++) {
    const double  coef = cm->face[f].meas*pty->value;
    for (int k = 0; k < 3; k++)
      kappa_f[f][k] = coef * cm->face[f].unitv[k];
  }

  /* 1) Build the bc_op matrix (scalar-valued version) */

  /* Initialize the matrix related to the flux reconstruction operator */

  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc, *bc_op_t = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix. Only the line associated to
     fb has non-zero values */

  _normal_flux_reco(fb, eqp->diffusion_hodgep.coef, cm,
                    (const cs_real_t (*)[3])kappa_f,
                    bc_op);

  /* 2) Update the bc_op matrix and the RHS with the Dirichlet values */

  /* Update bc_op = bc_op + transp and transp = transpose(bc_op)
     cb->loc plays the role of the flux operator */

  cs_sdm_square_add_transpose(bc_op, bc_op_t);

  /* Update the RHS with bc_op_t * dir_fb */

  for (int k = 0; k < 3; k++) {

    /* Only the fb column has non-zero values in bc_op_t
       One thus simplifies the matrix-vector operation */

    const cs_real_t  dir_fb = csys->dir_values[3*fb+k];
    for (short int i = 0; i < n_dofs; i++)
      csys->rhs[3*i+k] += bc_op_t->val[i*n_dofs+fb] * dir_fb;

  } /* Loop on components */

  /* 3) Update the bc_op matrix and the RHS with the penalization */

  /* coeff * \meas{f} / h_f  \equiv coeff * h_f */

  const double  pcoef = eqp->weak_pena_bc_coeff * sqrt(cm->face[fb].meas);

  bc_op->val[fb*(n_dofs + 1)] += pcoef; /* Diagonal term */

  for (short int k = 0; k < 3; k++)
    csys->rhs[3*fb + k] += pcoef * csys->dir_values[3*fb + k];

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
 * \brief  Take into account a boundary defined as 'symmetry' (treated as a
 *         sliding BCs on the three velocity components.)
 *         A weak penalization technique (symmetrized Nitsche) is used.  One
 *         assumes that static condensation has not been performed yet and that
 *         the velocity-block has (n_fc + 1) blocks of size 3x3.
 *         This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_symmetry(short int                       fb,
                  const cs_equation_param_t      *eqp,
                  const cs_cell_mesh_t           *cm,
                  const cs_property_data_t       *pty,
                  cs_cell_builder_t              *cb,
                  cs_cell_sys_t                  *csys)
{
  assert(cm != NULL && cb != NULL && csys != NULL && pty != NULL);
  assert(pty->is_iso == true); /* if not the case something else TODO ? */
  assert(cs_equation_param_has_diffusion(eqp));

  /* 0) Pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  for (short int f = 0; f < cm->n_fc; f++) {
    const double  coef = cm->face[f].meas*pty->value;
    for (int k = 0; k < 3; k++)
      kappa_f[f][k] = coef * cm->face[f].unitv[k];
  }

  /* 1) Build the bc_op matrix (scalar-valued version) */

  /* Initialize the matrix related this flux reconstruction operator */

  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */

  cs_sdm_t *bc_op = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* Compute \int_f du/dn v and update the matrix. Only the line associated to
     fb has non-zero values */

  _normal_flux_reco(fb, eqp->diffusion_hodgep.coef, cm,
                    (const cs_real_t (*)[3])kappa_f,
                    bc_op);

  /* 2) Update the bc_op matrix and nothing is added to the RHS since a sliding
     means homogeneous Dirichlet values on the normal component and hommogeneous
     Neumann on the tangential flux */

  const cs_quant_t  pfq = cm->face[fb];
  const cs_real_t  *nf = pfq.unitv;
  const cs_real_33_t  nf_nf = { {nf[0]*nf[0], nf[0]*nf[1], nf[0]*nf[2]},
                                {nf[1]*nf[0], nf[1]*nf[1], nf[1]*nf[2]},
                                {nf[2]*nf[0], nf[2]*nf[1], nf[2]*nf[2]} };

  /* coeff * \meas{f} / h_f  \equiv coeff * h_f */

  const double  pcoef = eqp->weak_pena_bc_coeff * sqrt(pfq.meas);

  /* Handle the diagonal block: Retrieve the 3x3 matrix */

  cs_sdm_t  *bFF = cs_sdm_get_block(csys->mat, fb, fb);
  assert(bFF->n_rows == bFF->n_cols && bFF->n_rows == 3);

  const cs_real_t  _val = pcoef + 2*bc_op->val[fb*(n_dofs+1)];
  for (short int k = 0; k < 3; k++) {
    bFF->val[3*k  ] += nf_nf[0][k] * _val;
    bFF->val[3*k+1] += nf_nf[1][k] * _val;
    bFF->val[3*k+2] += nf_nf[2][k] * _val;
  }

  for (short int xj = 0; xj < n_dofs; xj++) {

    if (xj == fb)
      continue;

    /* It should be done both for face- and cell-defined DoFs */
    /* Retrieve the 3x3 matrix */

    cs_sdm_t  *bFJ = cs_sdm_get_block(csys->mat, fb, xj);
    assert(bFJ->n_rows == bFJ->n_cols && bFJ->n_rows == 3);
    cs_sdm_t  *bJF = cs_sdm_get_block(csys->mat, xj, fb);
    assert(bJF->n_rows == bJF->n_cols && bJF->n_rows == 3);

    const cs_real_t  op_fj = bc_op->val[n_dofs*fb  + xj];
    const cs_real_t  op_jf = bc_op->val[n_dofs*xj + fb];
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
 *          This prototype matches the function pointer cs_cdo_apply_boundary_t
 *
 * \param[in]       fb        face id in the cell mesh numbering
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_fixed_wall(short int                       fb,
                    const cs_equation_param_t      *eqp,
                    const cs_cell_mesh_t           *cm,
                    const cs_property_data_t       *pty,
                    cs_cell_builder_t              *cb,
                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(cb);
  CS_UNUSED(pty);

  assert(cm != NULL && csys != NULL);  /* Sanity checks */

  const cs_quant_t  pfq = cm->face[fb];
  const cs_real_t  *ni = pfq.unitv;
  const cs_real_t  ni_ni[9] = { ni[0]*ni[0], ni[0]*ni[1], ni[0]*ni[2],
                                ni[1]*ni[0], ni[1]*ni[1], ni[1]*ni[2],
                                ni[2]*ni[0], ni[2]*ni[1], ni[2]*ni[2] };

  /* coeff * \meas{fb} / h_fb \equiv coeff * h_fb */
  const cs_real_t  pcoef = eqp->weak_pena_bc_coeff * sqrt(pfq.meas);

  cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, fb, fb);
  assert(bii->n_rows == bii->n_cols && bii->n_rows == 3);

  for (short int k = 0; k < 9; k++)
    bii->val[k] += pcoef * ni_ni[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one has to do one more non-linear iteration.
 *         Test if performed on the relative norm on the increment between
 *         two iterations
 *
 * \param[in]      nl_algo_type   type of non-linear algorithm
 * \param[in]      pre_iterate    previous state of the mass flux iterate
 * \param[in]      cur_iterate    current state of the mass flux iterate
 * \param[in, out] algo           pointer to a cs_iter_algo_t structure
 *
 * \return the convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_cdofb_navsto_nl_algo_cvg(cs_param_nl_algo_t        nl_algo_type,
                            const cs_real_t          *pre_iterate,
                            cs_real_t                *cur_iterate,
                            cs_iter_algo_t           *algo)
{
  assert(algo != NULL);

  if (nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON && algo->n_algo_iter > 0) {

    cs_iter_algo_param_aa_t  aap = cs_iter_algo_get_anderson_param(algo);

    switch (aap.dp_type) {

    case CS_PARAM_DOTPROD_EUCLIDEAN:
      cs_iter_algo_aa_update(algo,
                             cur_iterate,
                             pre_iterate,
                             cs_cdo_blas_dotprod_face,
                             cs_cdo_blas_square_norm_face);
      break;

    case CS_PARAM_DOTPROD_CDO:
      cs_iter_algo_aa_update(algo,
                             cur_iterate,
                             pre_iterate,
                             cs_cdo_blas_dotprod_pfsf,
                             cs_cdo_blas_square_norm_pfsf);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid case.\n", __func__);

    }

  } /* Anderson acceleration */

  /* Update the residual values. Compute the norm of the difference between the
     two mass fluxes (the current one and the previous one) */

  algo->prev_res = algo->res;
  algo->res = cs_cdo_blas_square_norm_pfsf_diff(pre_iterate, cur_iterate);
  assert(algo->res > -DBL_MIN);
  algo->res = sqrt(algo->res);

  if (algo->n_algo_iter < 1) /* Store the first residual to detect a
                                possible divergence of the algorithm */
    algo->res0 = algo->res;

  /* Update the convergence members */

  cs_iter_algo_update_cvg(algo);

  if (algo->param.verbosity > 0) {

    if (algo->n_algo_iter == 1)
      cs_log_printf(CS_LOG_DEFAULT,
                    "## %10s.It    Algo.Res   Inner  Cumul  Tolerance\n",
                    cs_param_get_nl_algo_label(nl_algo_type));
    cs_log_printf(CS_LOG_DEFAULT,
                  "## %10s.It%02d   %5.3e  %5d  %5d  %6.4e\n",
                  cs_param_get_nl_algo_label(nl_algo_type),
                  algo->n_algo_iter, algo->res, algo->last_inner_iter,
                  algo->n_inner_iter, algo->tol);

  } /* verbosity > 0 */

  return algo->cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the function pointer computing the source term in the momentum
 *         equation related to the gravity effect (hydrostatic pressure or the
 *         Boussinesq approximation)
 *
 * \param[in]  nsp          set of parameters for the Navier-Stokes system
 * \param[out] p_func       way to compute the gravity effect
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_set_gravity_func(const cs_navsto_param_t      *nsp,
                                 cs_cdofb_navsto_source_t    **p_func)
{
  if (nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ) {

    switch (cs_cdofb_navsto_boussinesq_type) {

    case CS_CDOFB_NAVSTO_BOUSSINESQ_FACE_DOF:
      *p_func = cs_cdofb_navsto_boussinesq_at_face;
      break;

    case CS_CDOFB_NAVSTO_BOUSSINESQ_CELL_DOF:
      *p_func = cs_cdofb_navsto_boussinesq_at_cell;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid type of algorithm to compute the Boussinesq"
                " approximation.\n", __func__);
    }

  }
  else if (nsp->model_flag & CS_NAVSTO_MODEL_GRAVITY_EFFECTS)
    *p_func = cs_cdofb_navsto_gravity_term;

  else
    *p_func = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the gravity effects.
 *         Compute and add the source term to the local RHS.
 *         This is a special treatment since of face DoFs are involved
 *         contrary to the standard case where only the cell DoFs is involved.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_gravity_term(const cs_navsto_param_t           *nsp,
                             const cs_cell_mesh_t              *cm,
                             const cs_cdofb_navsto_builder_t   *nsb,
                             cs_cell_sys_t                     *csys)
{
  assert(nsp->model_flag & CS_NAVSTO_MODEL_GRAVITY_EFFECTS);

  const cs_real_t  *gravity_vector = nsp->phys_constants->gravity;
  const cs_real_t  cell_coef = _dp3(gravity_vector, cm->xc);

  for (int f = 0; f < cm->n_fc; f++) {

    const cs_real_t  *_div_f = nsb->div_op + 3*f;
    const cs_quant_t  pfq = cm->face[f];
    const cs_real_t  face_coef = _dp3(gravity_vector, pfq.center);

    /* div_op is built such that _div_f[k] = -i_{f,c} * |f| * n_f[k] */

    for (int k = 0; k < 3; k++)
      csys->rhs[3*f+k] += _div_f[k] * nsb->rho_c * (cell_coef - face_coef);

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This is the standard case where the face DoFs are used for the
 *         constant part rho0 . g[] and only the cell DoFs are involved for the
 *         remaining part (the Boussinesq approximation).
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_boussinesq_at_cell(const cs_navsto_param_t           *nsp,
                                   const cs_cell_mesh_t              *cm,
                                   const cs_cdofb_navsto_builder_t   *nsb,
                                   cs_cell_sys_t                     *csys)
{
  CS_UNUSED(nsb);
  assert(nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ);

  /* Boussinesq term: rho0 * g[] * ( 1 - beta * (var[c] - var0) ).  The
   * remaining part rho0 * g[] * ( -beta * (var - var_c) has a zero mean-value
   * if one considers the reconstruction var = var_c + grad(vard)|_c * ( x -
   * x_c) which has a mean value equal to var_c */

  const cs_real_t  rho0 = nsp->mass_density->ref_value;
  const cs_real_t  *gravity_vector = nsp->phys_constants->gravity;

  cs_real_t  rho0g[3] = { rho0 * gravity_vector[0],
                          rho0 * gravity_vector[1],
                          rho0 * gravity_vector[2] };

  /* Constant part: rho_ref * g[] => Should be in balance with the pressure
     gradient in order to retrieve the hydrostatic case */

  for (int f = 0; f < cm->n_fc; f++) {
    const cs_real_t  *_div_f = nsb->div_op + 3*f;
    for (int k = 0; k < 3; k++)
      csys->rhs[3*f+k] += rho0g[k] * _div_f[k] * cm->xc[k];
  }

  /* Volume part  */

  double  boussi_coef = 0;

  for (int i = 0; i < nsp->n_boussinesq_terms; i++) {

    cs_navsto_param_boussinesq_t  *bp = nsp->boussinesq_param + i;
    boussi_coef += -bp->beta*(bp->var[cm->c_id] - bp->var0);

  }

  for (int k = 0; k < 3; k++)
    csys->rhs[3*cm->n_fc+k] += rho0g[k] * boussi_coef * cm->vol_c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Take into account the buoyancy force with the Boussinesq approx.
 *         Compute and add the source term to the local RHS.
 *         This way to compute the Boussinesq approximation applies only to
 *         face DoFs. This should enable to keep a stable (no velocity) in
 *         case of a stratified configuration.
 *
 * \param[in]      nsp     set of parameters to handle the Navier-Stokes system
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      nsb     pointer to a builder structure for the NavSto system
 * \param[in, out] csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_boussinesq_at_face(const cs_navsto_param_t           *nsp,
                                   const cs_cell_mesh_t              *cm,
                                   const cs_cdofb_navsto_builder_t   *nsb,
                                   cs_cell_sys_t                     *csys)
{
  CS_UNUSED(nsb);
  assert(nsp->model_flag & CS_NAVSTO_MODEL_BOUSSINESQ);

  /* Boussinesq term: rho0 * g[] * ( 1 - beta * (var[c] - var0) ). */

  /* 1. Compute the mass density for this cell taking into account the
     Boussinesq approximation */

  double  boussi_coef = 1;

  for (int i = 0; i < nsp->n_boussinesq_terms; i++) {

    cs_navsto_param_boussinesq_t  *bp = nsp->boussinesq_param + i;
    boussi_coef += -bp->beta*(bp->var[cm->c_id] - bp->var0);

  } /* Loop on Boussinesq terms */

  const cs_real_t  rho_c = nsp->mass_density->ref_value * boussi_coef;
  const cs_real_t  *gravity_vector = nsp->phys_constants->gravity;
  const double  cell_coef = _dp3(gravity_vector, cm->xc);

  for (int f = 0; f < cm->n_fc; f++) {

    /* div_op is built such that _div_f[k] = -i_{f,c} * |f| * n_f[k] */

    const cs_real_t  *_div_f = nsb->div_op + 3*f;
    const cs_quant_t  pfq = cm->face[f];
    const double  face_coef = _dp3(gravity_vector, pfq.center);
    const double rhs_coef = rho_c * (cell_coef - face_coef);

    for (int k = 0; k < 3; k++)
      csys->rhs[3*f+k] += rhs_coef * _div_f[k];

  } /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the source term for computing the stream function.
 *         This relies on the prototype associated to the generic function
 *         pointer \ref cs_dof_function_t
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        result of the function. Must be allocated.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_navsto_stream_source_term(cs_lnum_t            n_elts,
                                   const cs_lnum_t     *elt_ids,
                                   bool                 dense_output,
                                   void                *input,
                                   cs_real_t           *retval)
{
  assert(input != NULL);
  assert(retval != NULL);

  /* input is a pointer to the vorticity field */
  const cs_real_t  *w = (cs_real_t *)input;

  for (cs_lnum_t i = 0; i < n_elts; i++) {

    cs_lnum_t  id = (elt_ids == NULL) ? i : elt_ids[i];
    cs_lnum_t  r_id = dense_output ? i : id;

    retval[r_id] = w[3*id+2];   /* Extract the z component */

  }
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
