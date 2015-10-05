/*============================================================================
 * Build discrete convection operators for CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <limits.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_evaluate.h"

#include "cs_cdo_convection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdo_convection.c

  \brief Build discrete convection operators for CDO schemes

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CONVECTION_OP_DBG 0

struct _convection_builder_t {

  short int  *loc_v_ids;   /* Indirection array between a local number (in each
                              cell) and the (initial) vertex numbering */

  cs_real_t  *fluxes;      /* flux of the advection field across each dual
                              face (size: number of edges in a cell) */

  cs_real_t  *criter;      /* criterion used to evaluate how to upwind
                              (size: number of edges in a cell) */

  cs_toolbox_locmat_t   *loc;  // Local matrix for the convection operator

};

/*============================================================================
 * Local variables
 *============================================================================*/

static cs_real_t  one_third = 1./3;

/*============================================================================
 * Private constant variables
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value of the weighting function related to upwinding
 *
 * \param[in]  criterion  dot product between advection and normal vectors or
 *                        estimation of a local Peclet number
 * \param[in]  adv_info   set of options for the computation
 *
 * \return the weight value
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_upwind_weight(cs_real_t                   criterion,
               const cs_param_advection_t  adv_info)
{
  cs_real_t  weight = -1;

  switch (adv_info.weight_algo) {

  case CS_PARAM_ADVECTION_WEIGHT_ALGO_UPWIND:
    if (criterion > 0)
      weight = 1;
    else if (criterion < 0)
      weight = 0;
    else
      weight = 0.5;
    break;

  case CS_PARAM_ADVECTION_WEIGHT_ALGO_SAMARSKII:
    if (criterion < 0)
      weight = 1./(2 - criterion);
    else
      weight = (1 + criterion)/(2 + criterion);
    break;

  case CS_PARAM_ADVECTION_WEIGHT_ALGO_SG: // Sharfetter-Gummel
    if (criterion < 0)
      weight = 0.5*exp(criterion);
    else
      weight = 1 - 0.5*exp(-criterion);
    break;

  case CS_PARAM_ADVECTION_WEIGHT_ALGO_CENTERED:
    weight = 0.5;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible type of algorithm to compute the weight of"
              " upwind.");

  } // Switch on type of algo

  return weight;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the advective flux across dual faces lying inside the
 *          computational domain (only when the advective field is not uniform)
 *
 * \param[in]  quant      pointer to the cdo quantities structure
 * \param[in]  a_info     set of options to manage the advection term
 * \param[in]  tcur       value of the current time
 * \param[in]  qe         quantities related to an edge
 * \param[in]  qdf        quantities related to a dual face
 *
 * \return the value of the flux
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_compute_adv_flux(const cs_cdo_quantities_t   *quant,
                  const cs_param_advection_t   a_info,
                  double                       tcur,
                  const cs_real_3_t            xc,
                  const cs_quant_t             qe,
                  const cs_dface_t             qdf)
{
  int  j, k;
  cs_real_t  w;
  cs_real_3_t  beta_g, gpts[3], xg;

  cs_real_t  beta_flx = 0;

  /* Compute the flux of beta accros the dual face */
  switch (a_info.quad_type) {

  case CS_QUADRATURE_BARY:
    for (j = 0; j < 2; j++) {

      cs_lnum_t  f_id = qdf.parent_id[j];
      cs_quant_t  qf = quant->face[f_id];

      for (k = 0; k < 3; k++)
        xg[k] = one_third * (xc[k] + qe.center[k] + qf.center[k]);

      cs_evaluate_adv_field(a_info.adv_id, tcur, xg, &beta_g);
      beta_flx += qdf.meas[j] * _dp3(beta_g, &(qdf.unitv[3*j]));

    } // Loop on the two triangles composing the dual face inside a cell
    break;

  case CS_QUADRATURE_HIGHER:
    for (j = 0; j < 2; j++) {

      cs_lnum_t  f_id = qdf.parent_id[j];
      cs_quant_t  qf = quant->face[f_id];

      cs_quadrature_tria_3pts(qe.center, qf.center, xc, qdf.meas[j],
                              gpts, &w);

      cs_real_t  add = 0;
      for (int p = 0; p < 3; p++) {
        cs_evaluate_adv_field(a_info.adv_id, tcur, gpts[p], &beta_g);
        add += _dp3(beta_g, &(qdf.unitv[3*j]));
      }
      beta_flx += add * w;

    } // Loop on the two triangles composing the dual face inside a cell
    break;

  case CS_QUADRATURE_HIGHEST:
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of quadrature for computing the advected flux.");
  }

  return beta_flx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the builder structure and the local matrix related to
 *          the convection operator when diffusion is also activated
 *
 * \param[in]      c_id      cell id
 * \param[in]      connect   pointer to the connectivity structure
 * \param[in]      quant     pointer to the cdo quantities structure
 * \param[in]      a_info    set of options for the advection term
 * \param[in]      d_info    set of options for the diffusion term
 * \param[in]      tcur      value of the current time
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_with_diffusion(cs_lnum_t                    c_id,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     const cs_param_advection_t   a_info,
                     const cs_param_hodge_t       d_info,
                     double                       tcur,
                     cs_convection_builder_t     *builder)
{
  cs_lnum_t  i, k, id, n;
  cs_real_t  _crit, beta_flx, inv_val, inv_diff;
  cs_real_3_t  beta_c, xc, xg, beta_g, matnu;
  cs_real_33_t  matpty;
  cs_qvect_t  dface, beta;

  cs_toolbox_locmat_t  *loc = builder->loc;

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_connect_index_t  *c2e = connect->c2e;

  /* Keep the link between local cell numbering and vertex numbering */
  for (i = c2v->idx[c_id], n = 0; i < c2v->idx[c_id+1]; i++, n++) {
    id = c2v->ids[i];
    loc->ids[n] = id;
    builder->loc_v_ids[id] = n;
  }
  loc->n_ent = n;

  /* Initialize local matrix structure */
  for (i = 0; i < n*n; i++)
    loc->mat[i] = 0;

  bool  adv_uniform = cs_param_adv_field_is_uniform(a_info.adv_id);

  /* Retrieve cell center */
  for (k = 0; k < 3; k++)
    xc[k] = quant->cell_centers[3*c_id+k];

  /* Get the value of the material property at the cell center */
  cs_evaluate_pty(d_info.pty_id, tcur, xc, false, &matpty);

  /* Compute the flux induced by the advection field along with the criterion
     to evaluate how to upwind */

  if (adv_uniform) {

    cs_evaluate_adv_field(a_info.adv_id, tcur, xc, &beta_c);
    cs_qvect(beta_c, &beta);

    /* Loop on cell edges */
    for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

      cs_lnum_t  e_id = c2e->ids[i];
      cs_dface_t  qdf = quant->dface[i];

      /* Compute dual face vector split into a unit vector and its measure */
      dface.meas = qdf.meas[0] + qdf.meas[1];
      inv_val = 1/dface.meas;
      for (k = 0; k < 3; k++)
        dface.unitv[k] = inv_val*qdf.vect[k];

      _crit = beta.meas * _dp3(beta.unitv, dface.unitv);
      _mv3(matpty, dface.unitv, &matnu);
      inv_diff = 1/_dp3(dface.unitv, matnu);
      builder->criter[id] = quant->edge[e_id].meas * _crit * inv_diff;
      builder->fluxes[id] = dface.meas *_crit;

    } // Loop on cell edges

  }
  else { // advection is not uniform

    /* Loop on cell edges */
    for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

      cs_lnum_t  e_id = c2e->ids[i];
      cs_quant_t  qe = quant->edge[e_id];
      cs_dface_t  qdf = quant->dface[i];

      /* Compute dual face vector split into a unit vector and its measure */
      dface.meas = qdf.meas[0] + qdf.meas[1];
      inv_val = 1/dface.meas;
      for (k = 0; k < 3; k++)
        dface.unitv[k] = inv_val*qdf.vect[k];

      beta_flx = _compute_adv_flux(quant, a_info, tcur, xc, qe, qdf);

      /* Compute the criterion for evaluating the weight of upwinding */
      switch (a_info.weight_criterion) {

      case CS_PARAM_ADVECTION_WEIGHT_FLUX:
        _crit = beta_flx * inv_val;
        break;

      case CS_PARAM_ADVECTION_WEIGHT_XEXC:
        for (k = 0; k < 3; k++)
          xg[k] = 0.5 * (qe.center[k] + xc[k]);
        cs_evaluate_adv_field(a_info.adv_id, tcur, xg, &beta_g);
        _crit = _dp3(beta_g, dface.unitv);
        break;

      default:
        _crit = -1; // To avoid a warning at the compilation stage
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid choice of algorithm for computing the weight"
                  " criterion.");
      }

      _mv3(matpty, dface.unitv, &matnu);
      inv_diff = 1/_dp3(dface.unitv, matnu);
      builder->criter[id] = qe.meas * _crit * inv_diff;
      builder->fluxes[id] = beta_flx;

    } // Loop on cell edges

  } // advection uniform or not ?

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the builder structure and the local matrix related to
 *          the convection operator
 *
 * \param[in]      c_id      cell id
 * \param[in]      connect   pointer to the connectivity structure
 * \param[in]      quant     pointer to the cdo quantities structure
 * \param[in]      a_info    set of options for the advection term
 * \param[in]      tcur      value of the current time
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_init(cs_lnum_t                    c_id,
      const cs_cdo_connect_t      *connect,
      const cs_cdo_quantities_t   *quant,
      const cs_param_advection_t   a_info,
      double                       tcur,
      cs_convection_builder_t     *builder)
{
  cs_lnum_t  i, k, id, n;
  cs_real_t  _crit, beta_flx, invval;
  cs_real_3_t  beta_c, xc, xg, beta_g;
  cs_qvect_t  dface, beta;

  cs_toolbox_locmat_t  *loc = builder->loc;

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_connect_index_t  *c2e = connect->c2e;

  /* Keep the link between local cell numbering and vertex numbering */
  for (i = c2v->idx[c_id], n = 0; i < c2v->idx[c_id+1]; i++, n++) {
    id = c2v->ids[i];
    loc->ids[n] = id;
    builder->loc_v_ids[id] = n;
  }
  loc->n_ent = n;

  /* Initialize local matrix structure */
  for (i = 0; i < n*n; i++)
    loc->mat[i] = 0;

  /* Compute the flux induced by the advection field along with the criterion
     to evaluate how to upwind */
  bool  adv_uniform = cs_param_adv_field_is_uniform(a_info.adv_id);

  /* Retrieve cell center */
  for (k = 0; k < 3; k++)
    xc[k] = quant->cell_centers[3*c_id+k];

  if (adv_uniform) {

    cs_evaluate_adv_field(a_info.adv_id, tcur, xc, &beta_c);
    cs_qvect(beta_c, &beta);

    /* Loop on cell edges */
    for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

      cs_dface_t  qdf = quant->dface[i];

      /* Compute dual face vector split into a unit vector and its measure */
      dface.meas = qdf.meas[0] + qdf.meas[1];
      invval = 1/dface.meas;
      for (k = 0; k < 3; k++)
        dface.unitv[k] = invval*qdf.vect[k];

      _crit = beta.meas * _dp3(beta.unitv, dface.unitv);
      builder->criter[id] = _crit;
      builder->fluxes[id] = dface.meas *_crit;

    } // Loop on cell edges

  }
  else { // advection is not uniform

    /* Loop on cell edges */
    for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

      cs_lnum_t  e_id = c2e->ids[i];
      cs_quant_t  qe = quant->edge[e_id];
      cs_dface_t  qdf = quant->dface[i];

      /* Compute dual face vector split into a unit vector and its measure */
      dface.meas = qdf.meas[0] + qdf.meas[1];
      invval = 1/dface.meas;
      for (k = 0; k < 3; k++)
        dface.unitv[k] = invval*qdf.vect[k];

      beta_flx = _compute_adv_flux(quant, a_info, tcur, xc, qe, qdf);

      /* Compute the criterion for evaluating the weight of upwinding */
      switch (a_info.weight_criterion) {

      case CS_PARAM_ADVECTION_WEIGHT_FLUX:
        _crit = beta_flx * invval;
        break;

      case CS_PARAM_ADVECTION_WEIGHT_XEXC:
        for (k = 0; k < 3; k++)
          xg[k] = 0.5 * (qe.center[k] + xc[k]);
        cs_evaluate_adv_field(a_info.adv_id, tcur, xg, &beta_g);
        _crit = _dp3(beta_g, dface.unitv);
        break;

      default:
        _crit = 0; // To avoid a warning at the compilation stage
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid choice of algorithm for computing the weight"
                  " criterion.");
      }

      builder->criter[id] = _crit;
      builder->fluxes[id] = beta_flx;

    } // Loop on cell edges

  } // advection uniform or not ?

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local EpCd convection operator
 *
 * \param[in]      c_id      cell id
 * \param[in]      c2e       cell -> edges connectivity structure
 * \param[in]      e2v       edge -> vertices connectivity structure
 * \param[in]      a_info    set of options for the advection term
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_epcd(cs_lnum_t                    c_id,
                  const cs_connect_index_t    *c2e,
                  const cs_sla_matrix_t       *e2v,
                  const cs_param_advection_t   a_info,
                  cs_convection_builder_t     *builder)
{
  cs_lnum_t  i, id;

  /* Loop on cell edges */
  for (i = c2e->idx[c_id], id  = 0; i < c2e->idx[c_id+1]; i++, id++) {

    cs_real_t  beta_flx = builder->fluxes[id];

    if (fabs(beta_flx) > 0) {

      cs_lnum_t  e_id = c2e->ids[i];
      cs_lnum_t  e_shft = e2v->idx[e_id];
      cs_lnum_t  v1_id = e2v->col_id[e_shft];
      short int  sgn_v1 = e2v->sgn[e_shft];
      cs_lnum_t  v2_id = e2v->col_id[e_shft+1];
      short int  sgn_v2 = e2v->sgn[e_shft+1];  // fd(e),cd(v) = -v,e

      cs_real_t  criter = builder->criter[id];

      // Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e
      cs_real_t  weight_v1 = _upwind_weight(-sgn_v1 * criter, a_info);

      /* Update local convection matrix */
      short int  _v1 = builder->loc_v_ids[v1_id];
      short int  _v2 = builder->loc_v_ids[v2_id];
      int  shft_v1 = _v1*builder->loc->n_ent, shft_v2 = _v2*builder->loc->n_ent;
      cs_real_t  c1mw = beta_flx * (1 - weight_v1), cw = beta_flx * weight_v1;

      /* Sanity check */
      assert(_v1 != -1 && _v2 != -1);

      builder->loc->mat[shft_v1 + _v1] += sgn_v1 * c1mw;
      builder->loc->mat[shft_v2 + _v2] += sgn_v2 * cw;
      builder->loc->mat[shft_v1 + _v2] =  sgn_v2 * c1mw;
      builder->loc->mat[shft_v2 + _v1] =  sgn_v1 * cw;

    } // convective flux is greater than zero

  } // Loop on cell edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local VpCd convection operator
 *
 * \param[in]      c_id      cell id
 * \param[in]      c2e       cell -> edges connectivity structure
 * \param[in]      e2v       edge -> vertices connectivity structure
 * \param[in]      a_info    set of options for the advection term
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_vpfd(cs_lnum_t                    c_id,
                  const cs_connect_index_t    *c2e,
                  const cs_sla_matrix_t       *e2v,
                  const cs_param_advection_t   a_info,
                  cs_convection_builder_t     *builder)
{
  cs_lnum_t  i, id;

  /* Loop on cell edges */
  for (i = c2e->idx[c_id], id  = 0; i < c2e->idx[c_id+1]; i++, id++) {

    cs_real_t  beta_flx = builder->fluxes[id];

    if (fabs(beta_flx) > 0) {

      cs_lnum_t  e_id = c2e->ids[i];
      cs_lnum_t  e_shft = e2v->idx[e_id];
      cs_lnum_t  v1_id = e2v->col_id[e_shft];
      short int  sgn_v1 = e2v->sgn[e_shft];
      cs_lnum_t  v2_id = e2v->col_id[e_shft+1];
      short int  sgn_v2 = e2v->sgn[e_shft+1];  // fd(e),cd(v) = -v,e

      cs_real_t  criter = builder->criter[id];

      // Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e
      cs_real_t  weight_v1 = _upwind_weight(-sgn_v1 * criter, a_info);
      cs_real_t  weight_v2 = 1 - weight_v1;

      /* Update local convection matrix */
      short int  _v1 = builder->loc_v_ids[v1_id];
      short int  _v2 = builder->loc_v_ids[v2_id];
      int  shft_v1 = _v1*builder->loc->n_ent, shft_v2 = _v2*builder->loc->n_ent;
      cs_real_t  cw1 = beta_flx * weight_v1, cw2 = beta_flx * weight_v2;

      /* Sanity check */
      assert(_v1 != -1 && _v2 != -1);

      builder->loc->mat[shft_v1 + _v1] += -sgn_v1 * cw1;
      builder->loc->mat[shft_v2 + _v2] += -sgn_v2 * cw2;
      builder->loc->mat[shft_v1 + _v2] =  -sgn_v1 * cw2;
      builder->loc->mat[shft_v2 + _v1] =  -sgn_v2 * cw1;

    } // convective flux is greater than zero

  } // Loop on cell edges

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection flux accross the dual face df(e) lying
 *          inside the cell c and associated to the edge e.
 *          This function is associated to vertex-based discretization.
 *
 * \param[in]  quant    pointer to the cdo quantities structure
 * \param[in]  a_info   set of options for the advection term
 * \param[in]  tcur     value of the current time
 * \param[in]  xc       center of the cell c
 * \param[in]  qe       quantities related to edge e in E_c
 * \param[in]  qdf      quantities to the dual face df(e)
 *
 * \return the value of the convective flux accross the triangle
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_convection_vbflux_compute(const cs_cdo_quantities_t   *quant,
                             const cs_param_advection_t   a_info,
                             double                       tcur,
                             const cs_real_3_t            xc,
                             const cs_quant_t             qe,
                             const cs_dface_t             qdf)
{
  cs_real_t  beta_flx = 0;

  beta_flx = _compute_adv_flux(quant, a_info, tcur, xc, qe, qdf);

  return beta_flx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure for the convection operator
 *
 * \param[in]      connect  pointer to the connectivity structure
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

cs_convection_builder_t *
cs_convection_builder_init(const cs_cdo_connect_t      *connect)
{
  cs_lnum_t  i;

  cs_convection_builder_t  *b = NULL;

  cs_lnum_t  n_vertices = connect->v_info->n_ent;
  cs_lnum_t  n_ec = connect->n_max_ebyc;

  BFT_MALLOC(b, 1, cs_convection_builder_t);

  BFT_MALLOC(b->loc_v_ids, n_vertices, short int);
  for (i = 0; i < n_vertices; i++)
    b->loc_v_ids[i] = -1;

  BFT_MALLOC(b->fluxes, n_ec, cs_real_t);
  BFT_MALLOC(b->criter, n_ec, cs_real_t);
  for (i = 0; i < n_ec; i++) {
    b->fluxes[i] = 0;
    b->criter[i] = 0;
  }

  b->loc = cs_toolbox_locmat_create(connect->n_max_vbyc);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a builder structure for the convection operator
 *
 * \param[in, out] b   pointer to a cs_convection_builder_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_convection_builder_t *
cs_convection_builder_free(cs_convection_builder_t  *b)
{
  if (b == NULL)
    return b;

  BFT_FREE(b->loc_v_ids);
  BFT_FREE(b->fluxes);
  BFT_FREE(b->criter);

  b->loc = cs_toolbox_locmat_free(b->loc);

  BFT_FREE(b);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator for pure convection
 *
 * \param[in]      connect  pointer to the connectivity structure
 * \param[in]      quant    pointer to the cdo quantities structure
 * \param[in]      tcur     value of the current time
 * \param[in]      a_info   set of options for the advection term
 * \param[in, out] builder  pointer to a builder structure
 * \param[in, out] matrix   pointer to the matrix structure of the system
 */
/*----------------------------------------------------------------------------*/

void
cs_convection(const cs_cdo_connect_t      *connect,
              const cs_cdo_quantities_t   *quant,
              double                       tcur,
              const cs_param_advection_t   a_info,
              cs_convection_builder_t     *builder,
              cs_sla_matrix_t             *matrix)
{
  cs_lnum_t  c_id;

  if (a_info.form == CS_PARAM_ADVECTION_FORM_NONCONS) {

    for (c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Initialize the builder structure and compute geometrical quantities */
      _init(c_id, connect, quant, a_info, tcur, builder);

      /* Build the local convection operator */
      _build_local_epcd(c_id, connect->c2e, connect->e2v, a_info, builder);

      /* Update system matrix from this local contribution */
      cs_sla_assemble_msr(builder->loc, matrix);

    } /* Loop on cells */

  }
  else if (a_info.form == CS_PARAM_ADVECTION_FORM_CONSERV) {

    for (c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Initialize the builder structure and compute geometrical quantities */
      _init(c_id, connect, quant, a_info, tcur, builder);

      /* Build the local convection operator */
      _build_local_vpfd(c_id, connect->c2e, connect->e2v, a_info, builder);

      /* Update system matrix from this local contribution */
      cs_sla_assemble_msr(builder->loc, matrix);

    } /* Loop on cells */

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of advection operation.\n"
              " Choices are the following: conservative or not");

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator when diffusion is activated
 *
 * \param[in]      connect  pointer to the connectivity structure
 * \param[in]      quant    pointer to the cdo quantities structure
 * \param[in]      tcur     value of the current time
 * \param[in]      a_info   set of options for the advection term
 * \param[in]      d_info   set of options for the diffusion term
 * \param[in, out] builder  pointer to a builder structure
 * \param[in, out] matrix   pointer to the matrix structure of the system
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_with_diffusion(const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             double                       tcur,
                             const cs_param_advection_t   a_info,
                             const cs_param_hodge_t       d_info,
                             cs_convection_builder_t     *builder,
                             cs_sla_matrix_t             *matrix)
{
  cs_lnum_t   c_id;

  if (a_info.form == CS_PARAM_ADVECTION_FORM_NONCONS) {

    for (c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Initialize the builder structure and compute geometrical quantities */
      _init_with_diffusion(c_id, connect, quant, a_info, d_info, tcur, builder);

      /* Build the local convection operator */
      _build_local_epcd(c_id, connect->c2e, connect->e2v, a_info, builder);

      /* Update system matrix from this local contribution */
      cs_sla_assemble_msr(builder->loc, matrix);

    } /* End of loop on cells */

  }
  else if (a_info.form == CS_PARAM_ADVECTION_FORM_CONSERV) {

    for (c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Initialize the builder structure and compute geometrical quantities */
      _init_with_diffusion(c_id, connect, quant, a_info, d_info, tcur, builder);

      /* Build the local convection operator */
      _build_local_vpfd(c_id, connect->c2e, connect->e2v, a_info, builder);

      /* Update system matrix from this local contribution */
      cs_sla_assemble_msr(builder->loc, matrix);

    } /* End of loop on cells */

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of advection operation.\n"
              " Choices are the following: conservative or not");

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell in a given direction
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in]      a_info   set of options for the advection term
 * \param[in]      d_info   set of options for the diffusion term
 * \param[in]      dir_vect  direction in which we estimate the Peclet number
 * \param[in]      tcur      value of the current time
 * \param[in, out] peclet    pointer to the pointer of real numbers to fill
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_get_peclet_cell(const cs_cdo_quantities_t   *cdoq,
                              const cs_param_advection_t   a_info,
                              const cs_param_hodge_t       d_info,
                              const cs_real_3_t            dir_vect,
                              cs_real_t                    tcur,
                              cs_real_t                   *p_peclet[])
{
  cs_real_33_t  matpty;
  cs_real_3_t  beta_c, ptydir;

  cs_real_3_t  xc = {0, 0, 0};
  cs_real_t  *peclet = *p_peclet;

  bool  pty_uniform = cs_param_pty_is_uniform(d_info.pty_id);
  bool  adv_uniform = cs_param_adv_field_is_uniform(a_info.adv_id);

  if (peclet == NULL)
    BFT_MALLOC(peclet, cdoq->n_cells, cs_real_t);

  if (pty_uniform)
    cs_evaluate_pty(d_info.pty_id, tcur, xc, false, &matpty);
  if (adv_uniform)
    cs_evaluate_adv_field(a_info.adv_id, tcur, xc, &beta_c);

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    if (!pty_uniform || !adv_uniform) { /* Retrieve cell center */
      for (int k = 0; k < 3; k++)
        xc[k] = cdoq->cell_centers[3*c_id+k];

      /* Get the value of the material property at the cell center */
      if (!pty_uniform)
        cs_evaluate_pty(d_info.pty_id, tcur, xc, false, &matpty);

      if (!adv_uniform)
        cs_evaluate_adv_field(a_info.adv_id, tcur, xc, &beta_c);

    } // not uniform

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], one_third);
    cs_real_t  dp = _dp3(beta_c, dir_vect);

    _mv3(matpty, dir_vect, &ptydir);

    cs_real_t  inv_denum = 1/(_dp3(dir_vect, ptydir));

    peclet[c_id] = hc * dp * inv_denum;

  } // Loop on cells

  /* Return pointer */
  *p_peclet = peclet;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value in each cell of the upwinding coefficient given
 *          a related Peclet number
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_convection_get_upwind_coef_cell(const cs_cdo_quantities_t   *cdoq,
                                   const cs_param_advection_t   a_info,
                                   cs_real_t                    coefval[])
{
  /* Sanity check */
  assert(coefval != NULL);

  for (cs_lnum_t  c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_t  coef = _upwind_weight(coefval[c_id], a_info);

    coefval[c_id] = coef;

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
