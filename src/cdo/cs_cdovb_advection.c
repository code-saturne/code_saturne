/*============================================================================
 * Build discrete convection operators for CDO vertex-based schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_cdovb_advection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdovb_advection.c

  \brief Build discrete advection operators for CDO vertex-based schemes

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_ADVECTION_DBG 0

struct _cs_cdovb_adv_t {

  /* Settings for the advection operator */
  cs_param_advection_t    a_info;
  const cs_adv_field_t   *adv;   // shared pointer

  bool                    with_diffusion;

  cs_real_t  *fluxes;  /* flux of the advection field across each dual
                          face (size: number of edges in a cell) */
  cs_real_t  *criter;  /* criterion used to evaluate how to upwind
                          (size: number of edges in a cell) */

  cs_locmat_t   *loc;  /* Local matrix for the convection operator */

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

  case CS_PARAM_ADVECTION_WEIGHT_ALGO_D10G5: // Specific (JB)
    {
      cs_real_t  x = fabs(criterion);

      if (x > 0)
        weight = 0.5 + 0.5*( x*(x + 1) / ( x*(x + 5) + 10) );
      else
        weight = 0.5 - 0.5*( x*(x + 1) / ( x*(x + 5) + 10) );
    }
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
 * \brief   Initialize the builder structure and the local matrix related to
 *          the convection operator when diffusion is also activated
 *
 * \param[in]      c_id      cell id
 * \param[in]      connect   pointer to the connectivity structure
 * \param[in]      quant     pointer to the cdo quantities structure
 * \param[in]      matpty    tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_with_diffusion(cs_lnum_t                    c_id,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     const cs_real_33_t           matpty,
                     cs_cdovb_adv_t              *b)
{
  cs_lnum_t  i, k, id;
  cs_real_3_t  matnu;
  cs_nvec3_t  dface;

  const cs_connect_index_t  *c2e = connect->c2e;

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_advection_field_get_flux_dfaces(c_id, b->a_info, b->adv, b->fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

    const cs_lnum_t  e_id = c2e->ids[i];
    const cs_dface_t  qdf = quant->dface[i];

    /* Compute dual face vector split into a unit vector and its measure */
    const double  dfmeas = qdf.sface[0].meas + qdf.sface[1].meas;
    const double  inv_dfmeas = 1/dfmeas;
    const double  mean_flux = inv_dfmeas * b->fluxes[id];

    for (k = 0; k < 3; k++)
      dface.unitv[k] = inv_dfmeas*qdf.vect[k];

    _mv3((const cs_real_t (*)[3])matpty, dface.unitv, matnu);

    b->criter[id] = quant->edge[e_id].meas * mean_flux /_dp3(dface.unitv, matnu);

  } // Loop on cell edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the builder structure and the local matrix related to
 *          the convection operator
 *
 * \param[in]      c_id      cell id
 * \param[in]      connect   pointer to the connectivity structure
 * \param[in]      quant     pointer to the cdo quantities structure
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_init(cs_lnum_t                    c_id,
      const cs_cdo_connect_t      *connect,
      const cs_cdo_quantities_t   *quant,
      cs_cdovb_adv_t              *b)
{
  cs_lnum_t  i, id;

  const cs_connect_index_t  *c2e = connect->c2e;

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_advection_field_get_flux_dfaces(c_id, b->a_info, b->adv, b->fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  for (i = c2e->idx[c_id], id = 0; i < c2e->idx[c_id+1]; i++, id++) {

    cs_dface_t  qdf = quant->dface[i];

    /* Compute dual face vector split into a unit vector and its measure */
    const double  inv_dfmeas = 1/(qdf.sface[0].meas + qdf.sface[1].meas);

    b->criter[id] = inv_dfmeas * b->fluxes[id];

  } // Loop on cell edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local EpCd convection operator
 *
 * \param[in]      c_id      cell id
 * \param[in]      c2e       cell -> edges connectivity structure
 * \param[in]      e2v       edge -> vertices connectivity structure
 * \param[in]      loc_ids   store the value of the local id for each entity
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_epcd(cs_lnum_t                    c_id,
                  const cs_connect_index_t    *c2e,
                  const cs_sla_matrix_t       *e2v,
                  const cs_lnum_t             *loc_ids,
                  cs_cdovb_adv_t              *builder)
{
  cs_lnum_t  i, id;

  cs_locmat_t  *loc = builder->loc;

  const cs_param_advection_t  a_info = builder->a_info;

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
      short int  _v1 = loc_ids[v1_id];
      short int  _v2 = loc_ids[v2_id];
      int  shft_v1 = _v1*loc->n_ent, shft_v2 = _v2*loc->n_ent;
      cs_real_t  c1mw = beta_flx * (1 - weight_v1), cw = beta_flx * weight_v1;

      /* Sanity check */
      assert(_v1 != -1 && _v2 != -1);

      loc->mat[shft_v1 + _v1] += sgn_v1 * c1mw;
      loc->mat[shft_v2 + _v2] += sgn_v2 * cw;
      loc->mat[shft_v1 + _v2] =  sgn_v2 * c1mw;
      loc->mat[shft_v2 + _v1] =  sgn_v1 * cw;

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
 * \param[in]      loc_ids   store the value of the local id for each entity
 * \param[in, out] builder   pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_vpfd(cs_lnum_t                   c_id,
                  const cs_connect_index_t   *c2e,
                  const cs_sla_matrix_t      *e2v,
                  const cs_lnum_t            *loc_ids,
                  cs_cdovb_adv_t             *builder)
{
  cs_lnum_t  i, id;

  cs_locmat_t  *loc = builder->loc;

  const cs_param_advection_t  a_info = builder->a_info;

  /* Loop on cell edges */
  for (i = c2e->idx[c_id], id  = 0; i < c2e->idx[c_id+1]; i++, id++) {

    cs_real_t  beta_flx = builder->fluxes[id];

    if (fabs(beta_flx) > 0) { /* Flux induced by the advection field > 0 */

      const cs_lnum_t  e_id = c2e->ids[i];
      const cs_lnum_t  e_shft = e2v->idx[e_id];
      const cs_lnum_t  v1_id = e2v->col_id[e_shft];
      const cs_lnum_t  v2_id = e2v->col_id[e_shft+1];
      const short int  sgn_v1 = e2v->sgn[e_shft]; // fd(e),cd(v) = -v,e
      const short int  sgn_v2 = -sgn_v1;

      const cs_real_t  criter = builder->criter[id];

      /* Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e */
      const cs_real_t  weight_v1 = _upwind_weight(-sgn_v1 * criter, a_info);
      const cs_real_t  weight_v2 = 1 - weight_v1;
      const cs_real_t  cw1 = beta_flx * weight_v1, cw2 = beta_flx * weight_v2;

      /* Update local convection matrix */
      const short int  _v1 = loc_ids[v1_id];
      const short int  _v2 = loc_ids[v2_id];
      const int  shft_v1 = _v1*loc->n_ent, shft_v2 = _v2*loc->n_ent;

      /* Sanity check */
      assert(_v1 != -1 && _v2 != -1);

      /* Vertex _v1 */
      loc->mat[shft_v1 + _v1] += -sgn_v1 * cw1;
      loc->mat[shft_v1 + _v2] =  -sgn_v1 * cw2;
      /* Vertex _v2 */
      loc->mat[shft_v2 + _v2] += -sgn_v2 * cw2;
      loc->mat[shft_v2 + _v1] =  -sgn_v2 * cw1;

    } // convective flux is greater than zero

  } // Loop on cell edges

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure for the convection operator
 *
 * \param[in]  connect       pointer to the connectivity structure
 * \param[in]  adv_field     pointer to a cs_adv_field_t structure
 * \param[in]  a_info        set of options for the advection term
 * \param[in]  do_diffusion  true is diffusion is activated
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_adv_t *
cs_cdovb_advection_builder_init(const cs_cdo_connect_t      *connect,
                                const cs_adv_field_t        *adv,
                                const cs_param_advection_t   a_info,
                                bool                         do_diffusion)
{
  cs_cdovb_adv_t  *b = NULL;
  cs_lnum_t  n_max_ec = connect->n_max_ebyc;

  BFT_MALLOC(b, 1, cs_cdovb_adv_t);

  b->adv = adv; // share the pointer to an advection field structure

  /* Copy a cs_param_convection_t structure */
  b->a_info.formulation = a_info.formulation;
  b->a_info.weight_algo = a_info.weight_algo;
  b->a_info.weight_criterion = a_info.weight_criterion;
  b->a_info.quad_type = a_info.quad_type;

  b->with_diffusion = do_diffusion;

  /* Allocate and initialize buffers */
  BFT_MALLOC(b->fluxes, n_max_ec, cs_real_t);
  BFT_MALLOC(b->criter, n_max_ec, cs_real_t);
  for (int i = 0; i < n_max_ec; i++) {
    b->fluxes[i] = 0;
    b->criter[i] = 0;
  }

  b->loc = cs_locmat_create(connect->n_max_vbyc);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a builder structure for the convection operator
 *
 * \param[in, out] b   pointer to a cs_cdovb_adv_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_adv_t *
cs_cdovb_advection_builder_free(cs_cdovb_adv_t  *b)
{
  if (b == NULL)
    return b;

  BFT_FREE(b->fluxes);
  BFT_FREE(b->criter);

  b->loc = cs_locmat_free(b->loc);

  BFT_FREE(b);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator for pure convection
 *
 * \param[in]      c_id       cell id
 * \param[in]      connect    pointer to the connectivity structure
 * \param[in]      quant      pointer to the cdo quantities structure
 * \param[in]      loc_ids    store the local entity ids for this cell
 * \param[in]      diffmat    tensor related to the diffusion property
 * \param[in, out] builder    pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_advection_build_local(cs_lnum_t                    c_id,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_lnum_t             *loc_ids,
                               const cs_real_33_t           diffmat,
                               cs_cdovb_adv_t              *builder)
{
  cs_lnum_t  i;
  short int  n_cell_vertices = 0;

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_param_advection_t  a_info = builder->a_info;

  /* Initialize the builder structure and compute geometrical quantities */
  /* Keep the link between local cell numbering and vertex numbering */
  if (builder->with_diffusion)
    _init_with_diffusion(c_id, connect, quant, diffmat, builder);
  else
    _init(c_id, connect, quant, builder);

  /* Initialize local matrix structure */
  for (i = c2v->idx[c_id]; i < c2v->idx[c_id+1]; i++)
    builder->loc->ids[n_cell_vertices++] = c2v->ids[i];
  builder->loc->n_ent = n_cell_vertices;
  for (i = 0; i < n_cell_vertices*n_cell_vertices; i++)
    builder->loc->mat[i] = 0;

  /* Build the local convection operator */
  switch (a_info.formulation) {

  case CS_PARAM_ADVECTION_FORM_NONCONS:
    _build_local_epcd(c_id, connect->c2e, connect->e2v, loc_ids,  builder);
    break;

  case CS_PARAM_ADVECTION_FORM_CONSERV:
    _build_local_vpfd(c_id, connect->c2e, connect->e2v, loc_ids, builder);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of advection operation.\n"
              " Choices are the following: conservative or not");
    break;

  } // Switch on the formulation

  return builder->loc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator for pure convection
 *
 * \param[in]      connect       pointer to the connectivity structure
 * \param[in]      quant         pointer to the cdo quantities structure
 * \param[in]      dir_vals      values of the Dirichlet boundary condition
 * \param[in, out] builder       pointer to a convection builder structure
 * \param[in, out] rhs_contrib   array storing the rhs contribution
 * \param[in, out] diag_contrib  array storing the diagonal contribution
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_add_bc(const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *dir_vals,
                          cs_cdovb_adv_t              *builder,
                          cs_real_t                    rhs_contrib[],
                          cs_real_t                    diag_contrib[])
{
  cs_lnum_t  i, f_id;
  cs_nvec3_t  advf;

  const cs_adv_field_t  *adv = builder->adv;
  const cs_param_advection_t  a_info = builder->a_info;
  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  /* Loop on border faces.
     Add diagonal term for vertices attached to a boundary face where
     the advection field points inward */
  if (cs_advection_field_is_cellwise(adv)) {

    for (f_id = quant->n_i_faces; f_id < quant->n_faces; f_id++) {

      const cs_quant_t  qf = quant->face[f_id];

      /* Sanity check (this is a border face) */
      assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

      const cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];

      /* Retrieve the value of the advection field in the current cell */
      cs_advection_field_get_cell_vector(c_id, adv, &advf);

      const double  dp = _dp3(advf.unitv, qf.unitv);
      if (fabs(dp) > cs_get_zero_threshold()) {

        /* Loop on border face edges */
        for (i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

          const cs_lnum_t  e_id = f2e->col_id[i];
          const cs_lnum_t  e_shft = e2v->idx[e_id];
          const cs_lnum_t  v1_id = e2v->col_id[e_shft];
          const cs_lnum_t  v2_id = e2v->col_id[e_shft+1];
          const cs_real_t  *xv1 = xyz + 3*v1_id;
          const cs_real_t  *xv2 = xyz + 3*v2_id;
          const double  surf = 0.5*cs_surftri(xv1, xv2, qf.center);
          const double  _flx = dp * advf.meas * surf;

          if (dp < 0) { // advection field is inward w.r.t. the face normal

            rhs_contrib[v1_id] -= _flx * dir_vals[v1_id];
            rhs_contrib[v2_id] -= _flx * dir_vals[v2_id];

            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS) {
              diag_contrib[v1_id] -= _flx;
              diag_contrib[v2_id] -= _flx;
            }

          }
          else { // advection is oriented outward

            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV) {
              diag_contrib[v1_id] += _flx;
              diag_contrib[v2_id] += _flx;
            }

          }

        } // Loop on face edges

      } // abs(dp) > 0

    } // Loop on border faces

  }
  else { // Advection field is not uniform

    for (f_id = quant->n_i_faces; f_id < quant->n_faces; f_id++) {

      /* Sanity check (this is a border face) */
      assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

      /* Loop on border face edges */
      for (i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

        const cs_lnum_t  e_id = f2e->col_id[i];
        const cs_lnum_t  e_shft = e2v->idx[e_id];
        const cs_lnum_t  v1_id = e2v->col_id[e_shft];
        const cs_lnum_t  v2_id = e2v->col_id[e_shft+1];

        const double  flux_v1 =
          cs_advection_field_get_flux_svef(v1_id, e_id, f_id, a_info, adv);
        const double  flux_v2 =
          cs_advection_field_get_flux_svef(v2_id, e_id, f_id, a_info, adv);

        if (flux_v1 < 0) { // advection field is inward w.r.t. the face normal

          rhs_contrib[v1_id] -= flux_v1 * dir_vals[v1_id];
          if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS)
            diag_contrib[v1_id] -= flux_v1;

        }
        else  // advection is oriented outward
          if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV)
            diag_contrib[v1_id] += flux_v1;

        if (flux_v2 < 0) { // advection field is inward w.r.t. the face normal

          rhs_contrib[v2_id] -= flux_v2 * dir_vals[v2_id];
          if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS)
            diag_contrib[v2_id] -= flux_v2;

        }
        else  // advection is oriented outward
          if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV)
            diag_contrib[v2_id] += flux_v2;

      } // Loop on face edges

    } // Loop on border faces

  } // Advection field uniform or not

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell in a given direction
 *
 * \param[in]      cdoq           pointer to the cdo quantities structure
 * \param[in]      adv            pointer to the advection field struct.
 * \param[in]      diff_property  pointer to the diffusion property struct.
 * \param[in]      dir_vect       direction for estimating the Peclet number
 * \param[in, out] peclet         pointer to the pointer of real numbers to fill
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_get_peclet_cell(const cs_cdo_quantities_t   *cdoq,
                                   const cs_adv_field_t        *adv,
                                   const cs_property_t         *diff_property,
                                   const cs_real_3_t            dir_vect,
                                   cs_real_t                   *p_peclet[])
{
  cs_real_t  ptymat[3][3];
  cs_real_3_t  ptydir;
  cs_nvec3_t  adv_field;

  cs_real_t  *peclet = *p_peclet;
  bool  pty_uniform = cs_property_is_uniform(diff_property);

  if (peclet == NULL)
    BFT_MALLOC(peclet, cdoq->n_cells, cs_real_t);

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    /* Get the value of the material property at the cell center */
    if (!pty_uniform || c_id == 0)
      cs_property_get_cell_tensor(c_id, diff_property, false, ptymat);

    cs_advection_field_get_cell_vector(c_id, adv, &adv_field);

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], one_third);
    cs_real_t  dp = adv_field.meas * _dp3(adv_field.unitv, dir_vect);

    _mv3((const cs_real_t (*)[3])ptymat, dir_vect, ptydir);

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
 * \param[in]      a_info    set of options for the advection term
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_get_upwind_coef_cell(const cs_cdo_quantities_t   *cdoq,
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
