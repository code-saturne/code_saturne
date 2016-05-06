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

#include "cs_math.h"

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
 * Local Macro definitions
 *============================================================================*/

#define CS_CDOVB_ADVECTION_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

struct _cs_cdovb_adv_t {

  /* Settings for the advection operator */
  cs_param_advection_t    a_info;

  const cs_adv_field_t   *adv;   // shared pointer

  bool           with_diffusion;

  cs_real_t     *fluxes;    /* flux of the advection field across each dual
                               face (size: number of edges in a cell) */
  cs_real_t     *criter;    /* criterion used to evaluate how to upwind
                               (size: number of edges in a cell) */
  cs_locmat_t   *loc;       /* Local matrix for the convection operator */
  cs_real_t     *tmp_rhs;   /* local buffer (size: n_max_vbyc) */
};

/*============================================================================
 * Local variables
 *============================================================================*/

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

inline static double
_upwind_weight(double                      criterion,
               const cs_param_advection_t  adv_info)
{
  double  weight = -1;

  switch (adv_info.scheme) {

  case CS_PARAM_ADVECTION_SCHEME_UPWIND:
    if (criterion > 0)
      weight = 1;
    else if (criterion < 0)
      weight = 0;
    else
      weight = 0.5;
    break;

  case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
    if (criterion < 0)
      weight = 1./(2 - criterion);
    else
      weight = (1 + criterion)/(2 + criterion);
    break;

  case CS_PARAM_ADVECTION_SCHEME_SG: // Sharfetter-Gummel
    if (criterion < 0)
      weight = 0.5*exp(criterion);
    else
      weight = 1 - 0.5*exp(-criterion);
    break;

  case CS_PARAM_ADVECTION_SCHEME_CENTERED:
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
 * \param[in]      lm        pointer to a cs_cdo_locmesh_t structure
 * \param[in]      matpty    tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_with_diffusion(const cs_cdo_locmesh_t      *lm,
                     const cs_real_33_t           matpty,
                     cs_cdovb_adv_t              *b)
{
  cs_real_3_t  matnu;

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  for (short int e = 0; e < lm->n_ec; e++) {

    const cs_nvec3_t  dfq = lm->dface[e];
    const double  mean_flux = b->fluxes[e]/dfq.meas;

    cs_math_33_3_product((const cs_real_t (*)[3])matpty, dfq.unitv, matnu);

    cs_real_t  diff_contrib = _dp3(dfq.unitv, matnu);
    if (diff_contrib > cs_math_zero_threshold)
      b->criter[e] = lm->edge[e].meas * mean_flux / diff_contrib;
    else
      b->criter[e] = mean_flux * cs_math_big_r; // dominated by convection

  } // Loop on cell edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal edges and dual
 *          cells. (Non-conservative formulation)
 *
 * \param[in]      lm        pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_epcd(const cs_cdo_locmesh_t    *lm,
                  cs_cdovb_adv_t            *b)
{
  cs_locmat_t  *m = b->loc;

  const cs_param_advection_t  a_info = b->a_info;

  if (a_info.scheme == CS_PARAM_ADVECTION_SCHEME_CENTERED) {

    /* Weight is always equal to 0.5
       Loop on cell edges */
    for (short int e = 0; e < lm->n_ec; e++) {

      const cs_real_t  wflx = 0.5*b->fluxes[e];

      if (fabs(wflx) > 0) {

        short int  shft = 2*e;
        short int  sgn_v1 = lm->e2v_sgn[shft]; // fd(e),cd(v) = -v,e

        /* Update local convection matrix */
        short int  v1 = lm->e2v_ids[shft];
        short int  v2 = lm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1);        // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;
        const double contrib = sgn_v1 * wflx;

        m1[v1] +=  contrib;
        m1[v2] =  -contrib;
        m2[v2] += -contrib;
        m2[v1] =   contrib;

      } // convective flux is greater than zero

    } // Loop on cell edges

  }
  else {

    /* Loop on cell edges */
    for (short int e = 0; e < lm->n_ec; e++) {

      const cs_real_t  beta_flx = b->fluxes[e];

      if (fabs(beta_flx) > 0) {

        short int  shft = 2*e;
        short int  sgn_v1 = lm->e2v_sgn[shft];

        /* Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e */
        const double  wv1 = _upwind_weight(-sgn_v1 * b->criter[e], a_info);
        const double  c1mw = sgn_v1 * beta_flx * (1 - wv1);
        const double  cw = sgn_v1 * beta_flx * wv1;

        /* Update local convection matrix */
        short int  v1 = lm->e2v_ids[shft];
        short int  v2 = lm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1); // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;

        m1[v1] +=  c1mw;
        m1[v2] =  -c1mw; // sgn_v2 = -sgn_v1
        m2[v2] += -cw;   // sgn_v2 = -sgn_v1
        m2[v1] =   cw;

      } // convective flux is greater than zero

    } // Loop on cell edges

  } // Convection scheme is not centered

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Conservative formulation)
 *
 * \param[in]      lm        pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_vpfd(const cs_cdo_locmesh_t    *lm,
                  cs_cdovb_adv_t            *b)
{
  cs_locmat_t  *m = b->loc;

  const cs_param_advection_t  a_info = b->a_info;

  if (a_info.scheme == CS_PARAM_ADVECTION_SCHEME_CENTERED) {

    /* Weight is always equal to 0.5
       Loop on cell edges */
    for (short int e = 0; e < lm->n_ec; e++) {

      const cs_real_t  wflx = 0.5*b->fluxes[e];

      if (fabs(wflx) > 0) {

        short int  shft = 2*e;
        short int  sgn_v1 = lm->e2v_sgn[shft]; // fd(e),cd(v) = -v,e
        const double contrib = sgn_v1 * wflx;

        /* Update local convection matrix */
        short int  v1 = lm->e2v_ids[shft];
        short int  v2 = lm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1);        // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;

        m1[v1] += -contrib;
        m1[v2] =  -contrib;
        m2[v2] +=  contrib; // sgn_v2 = -sgn_v1
        m2[v1] =   contrib; // sgn_v2 = -sgn_v1

      } // convective flux is greater than zero

    } // Loop on cell edges

  }
  else {

    /* Loop on cell edges */
    for (short int e = 0; e < lm->n_ec; e++) {

      const cs_real_t  beta_flx = b->fluxes[e];

      if (fabs(beta_flx) > 0) {

        short int  shft = 2*e;
        short int  sgn_v1 = lm->e2v_sgn[shft];

        /* Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e */
        const double  wv1 = _upwind_weight(-sgn_v1 * b->criter[e], a_info);
        const double  cw1 = sgn_v1 * beta_flx * wv1;
        const double  cw2 = sgn_v1 * beta_flx * (1 - wv1);

        /* Update local convection matrix */
        short int  v1 = lm->e2v_ids[shft];
        short int  v2 = lm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1); // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;

        m1[v1] += -cw1;
        m1[v2] =  -cw2;
        m2[v2] +=  cw2;
        m2[v1] =   cw1;

      } // convective flux is greater than zero

    } // Loop on cell edges

  } // Convection scheme is not centered

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
  cs_lnum_t  n_max_ec = connect->n_max_ebyc, n_max_vc = connect->n_max_vbyc;

  BFT_MALLOC(b, 1, cs_cdovb_adv_t);

  b->adv = adv; // share the pointer to an advection field structure

  /* Copy a cs_param_convection_t structure */
  b->a_info.formulation = a_info.formulation;
  b->a_info.scheme = a_info.scheme;
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

  BFT_MALLOC(b->tmp_rhs, n_max_vc, cs_real_t);
  for (int i = 0; i < n_max_vc; i++)
    b->tmp_rhs[i] = 0;

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
  BFT_FREE(b->tmp_rhs);

  b->loc = cs_locmat_free(b->loc);

  BFT_FREE(b);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell
 *
 * \param[in]      lm        pointer to a cs_cdo_locmesh_t structure
 * \param[in]      diffmat   tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_advection_build_local(const cs_cdo_locmesh_t      *lm,
                               const cs_real_33_t           diffmat,
                               cs_cdovb_adv_t              *b)
{
  /* Initialize local matrix structure */
  for (short int v = 0; v < lm->n_vc; v++)
    b->loc->ids[v] = lm->v_ids[v];
  b->loc->n_ent = lm->n_vc;
  for (short int v = 0; v < lm->n_vc*lm->n_vc; v++)
    b->loc->val[v] = 0;

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_advection_field_get_flux_dfaces(lm->c_id, b->a_info, b->adv, b->fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  if (b->a_info.scheme != CS_PARAM_ADVECTION_SCHEME_CENTERED) {
    if (b->with_diffusion)
      _init_with_diffusion(lm, diffmat, b);
    else
      for (short int e = 0; e < lm->n_ec; e++)
        b->criter[e] = 1/lm->dface[e].meas * b->fluxes[e];
  }

  /* Build the local convection operator */
  switch (b->a_info.formulation) {

  case CS_PARAM_ADVECTION_FORM_NONCONS:
    _build_local_epcd(lm,  b);
    break;

  case CS_PARAM_ADVECTION_FORM_CONSERV:
    _build_local_vpfd(lm, b);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of advection operation.\n"
              " Choices are the following: conservative or not");
    break;

  } // Switch on the formulation

  return b->loc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      quant         pointer to the cdo quantities structure
 * \param[in]      lm            pointer to a cs_cdo_locmesh_t struct.
 * \param[in, out] advb          pointer to a convection builder structure
 * \param[in, out] ls            cell-wise structure sotring the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_add_bc(const cs_cdo_quantities_t   *quant,
                          cs_cdo_locmesh_t            *lm,
                          cs_cdovb_adv_t              *advb,
                          cs_cdo_locsys_t             *ls)
{
  cs_nvec3_t  adv_vec;
  cs_locmat_t  *m = advb->loc; /* Use this local dense matrix as an array
                                  taking into account the diagonal contrib. */

  const cs_adv_field_t  *adv_field = advb->adv;
  const cs_param_advection_t  a_info = advb->a_info;
  const cs_real_t  *xyz = quant->vtx_coord;

  /* Reset local temporay RHS and diagonal contributions */
  for (short int v = 0; v < lm->n_vc; v++) {
    m->val[v] = 0;
    advb->tmp_rhs[v] = 0;
  }

  /* Loop on border faces.
     Add diagonal term for vertices attached to a boundary face where
     the advection field points inward */
  if (cs_advection_field_is_cellwise(adv_field)) {

    for (short int f = 0; f < lm->n_fc; f++) {

      cs_lnum_t  f_id = lm->f_ids[f];
      if (f_id >= quant->n_i_faces) { // Border face

      const cs_quant_t  pfq = lm->face[f];

      /* Retrieve the value of the advection field in the current cell */
      cs_advection_field_get_cell_vector(lm->c_id, adv_field, &adv_vec);

        const double  dp = _dp3(adv_vec.unitv, pfq.unitv);
        if (fabs(dp) > cs_math_zero_threshold) {

          /* Loop on border face edges */
          for (int i = lm->f2e_idx[f]; i < lm->f2e_idx[f+1]; i++) {

            const short int  e = lm->f2e_ids[i];
            const short int  v1 = lm->e2v_ids[2*e];
            const short int  v2 = lm->e2v_ids[2*e+1];
            const double  surf = 0.5 * cs_math_surftri(xyz + 3*lm->v_ids[v1],
                                                       xyz + 3*lm->v_ids[v2],
                                                       pfq.center);
            const double  flx = dp * adv_vec.meas * surf;

            if (dp < 0) { // advection field is inward w.r.t. the face normal

              advb->tmp_rhs[v1] -= flx * ls->dir_bc[v1];
              advb->tmp_rhs[v2] -= flx * ls->dir_bc[v2];

              if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS) {
                m->val[v1] -= flx;
                m->val[v2] -= flx;
              }

            }
            else { // advection is oriented outward

              if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV) {
                m->val[v1] += flx;
                m->val[v2] += flx;
              }

            }

          } // Loop on face edges

        } // abs(dp) > 0

      } // If border face
    } // Loop on cell faces

  }
  else { // Advection field is not uniform

    for (short int f = 0; f < lm->n_fc; f++) {

      cs_lnum_t  f_id = lm->f_ids[f];
      if (f_id >= quant->n_i_faces) { // Border face

        /* Loop on border face edges */
        for (int i = lm->f2e_idx[f]; i < lm->f2e_idx[f+1]; i++) {

          const short int  e = lm->f2e_ids[i];
          const cs_lnum_t  e_id = lm->e_ids[e];
          const short int  v1 = lm->e2v_ids[2*e];
          const short int  v2 = lm->e2v_ids[2*e+1];
          const double  flx1 = cs_advection_field_get_flux_svef(lm->v_ids[v1],
                                                                e_id, f_id,
                                                                a_info,
                                                                adv_field);

          if (flx1 < 0) { // advection field is inward w.r.t. the face normal

            advb->tmp_rhs[v1] -= flx1 * ls->dir_bc[v1];
            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS)
              m->val[v1] -= flx1;

          }
          else  // advection is oriented outward
            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV)
              m->val[v1] += flx1;

          const double  flx2 = cs_advection_field_get_flux_svef(lm->v_ids[v2],
                                                                e_id, f_id,
                                                                a_info,
                                                                adv_field);

          if (flx2 < 0) { // advection field is inward w.r.t. the face normal

            advb->tmp_rhs[v2] -= flx2 * ls->dir_bc[v2];
            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS)
              m->val[v2] -= flx2;

          }
          else  // advection is oriented outward
            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV)
              m->val[v2] += flx2;

        } // Loop on face edges

      } // Loop on border faces

    } // Loop on cell faces

  } // Advection field uniform or not

  /* Update the diagonal and the RHS of the local system matrix */
  for (short int v = 0; v < lm->n_vc; v++)  {
    ls->mat->val[v*lm->n_vc + v] += m->val[v];
    ls->rhs[v] += advb->tmp_rhs[v];
  }

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

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);
    cs_real_t  dp = adv_field.meas * _dp3(adv_field.unitv, dir_vect);

    cs_math_33_3_product((const cs_real_t (*)[3])ptymat, dir_vect, ptydir);

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

#undef _dp3

END_C_DECLS
