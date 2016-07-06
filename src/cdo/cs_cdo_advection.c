/*============================================================================
 * Build discrete convection operators for CDO schemes
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

#include "cs_cdo_advection.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdo_advection.c

  \brief Build discrete advection operators for CDO vertex-based schemes

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDO_ADVECTION_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

struct _cs_cdo_adv_t {

  cs_lnum_t      n_i_faces;  // In order to detect border faces
  bool           with_diffusion;

  /* Temporary buffers to store useful quantities needed during the definition
     of the advection operator */
  cs_real_3_t   *tmp_vect;
  cs_real_t     *tmp_scal;

  cs_locmat_t   *loc;       /* Local matrix for the convection operator */

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
 * \brief   Define the local convection operator between primal edges and dual
 *          cells. (Non-conservative formulation)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      a_info   set of parameters related to the advection operator
 * \param[in, out] b        pointer to a cs_cdo_adv_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_epcd(const cs_cell_mesh_t        *cm,
                  const cs_param_advection_t   a_info,
                  cs_cdo_adv_t                *b)
{
  assert(a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS); // Sanity check

  cs_locmat_t  *m = b->loc;

  const cs_real_t  *fluxes = b->tmp_scal;
  const cs_real_t  *upwcoef = b->tmp_scal + cm->n_ec;

  if (a_info.scheme == CS_PARAM_ADVECTION_SCHEME_CENTERED) {

    /* Weight is always equal to 0.5
       Loop on cell edges */
    for (short int e = 0; e < cm->n_ec; e++) {

      const short int  shft = 2*e;
      const cs_real_t  wflx = 0.5 * fluxes[e] * cm->e2v_sgn[shft]; // sgn_v1

      if (fabs(wflx) > 0) {

        /* Update local convection matrix */
        short int  v1 = cm->e2v_ids[shft];
        short int  v2 = cm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1);        // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;

        // Use the fact that fd(e),cd(v) = -v,e and sgn_v1 = -sgn_v2
        m1[v1] +=  wflx;
        m1[v2] =  -wflx;
        m2[v2] += -wflx;
        m2[v1] =   wflx;

      } // convective flux is greater than zero

    } // Loop on cell edges

  }
  else {

    /* Loop on cell edges */
    for (short int e = 0; e < cm->n_ec; e++) {

      const short int  shft = 2*e;
      const short int  sgn_v1 = cm->e2v_sgn[shft];
      const cs_real_t  beta_flx = fluxes[e] * sgn_v1;

      if (fabs(beta_flx) > 0) {

        /* Compute the upwind coefficient knowing that fd(e),cd(v) = -v,e */
        const double  wv1 = _upwind_weight(-sgn_v1 * upwcoef[e], a_info);
        const double  c1mw = beta_flx * (1 - wv1);
        const double  cw = beta_flx * wv1;

        /* Update local convection matrix */
        short int  v1 = cm->e2v_ids[shft];
        short int  v2 = cm->e2v_ids[shft+1];
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
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      a_info   set of parameters related to the advection operator
 * \param[in, out] b        pointer to a cs_cdo_adv_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_vpfd(const cs_cell_mesh_t         *cm,
                  const cs_param_advection_t    a_info,
                  cs_cdo_adv_t                 *b)
{
  assert(a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV); // Sanity check

  cs_locmat_t  *m = b->loc;

  const cs_real_t  *fluxes = b->tmp_scal;
  const cs_real_t  *upwcoef = b->tmp_scal + cm->n_ec;

  if (a_info.scheme == CS_PARAM_ADVECTION_SCHEME_CENTERED) {

    /* Weight is always equal to 0.5
       Loop on cell edges */
    for (short int e = 0; e < cm->n_ec; e++) {

      const short int  shft = 2*e;
      const cs_real_t  wflx = 0.5*fluxes[e]*cm->e2v_sgn[shft];

      if (fabs(wflx) > 0) {

        /* Update local convection matrix */
        short int  v1 = cm->e2v_ids[shft];
        short int  v2 = cm->e2v_ids[shft+1];
        assert(v1 != -1 && v2 != -1);        // Sanity check
        double  *m1 = m->val + v1*m->n_ent, *m2 = m->val + v2*m->n_ent;

        m1[v1] += -wflx;
        m1[v2] =  -wflx;
        m2[v2] +=  wflx; // sgn_v2 = -sgn_v1
        m2[v1] =   wflx; // sgn_v2 = -sgn_v1

      } // convective flux is greater than zero

    } // Loop on cell edges

  }
  else {

    /* Loop on cell edges */
    for (short int e = 0; e < cm->n_ec; e++) {

      const cs_real_t  beta_flx = fluxes[e];

      if (fabs(beta_flx) > 0) {

        short int  shft = 2*e;
        short int  sgn_v1 = cm->e2v_sgn[shft];

        /* Compute the updwind coefficient knowing that fd(e),cd(v) = -v,e */
        const double  wv1 = _upwind_weight(-sgn_v1 * upwcoef[e], a_info);
        const double  cw1 = sgn_v1 * beta_flx * wv1;
        const double  cw2 = sgn_v1 * beta_flx * (1 - wv1);

        /* Update local convection matrix */
        short int  v1 = cm->e2v_ids[shft];
        short int  v2 = cm->e2v_ids[shft+1];
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
 * \param[in]  eqp           pointer to a cs_equation_param_t structure
 * \param[in]  do_diffusion  true is diffusion is activated
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_adv_t *
cs_cdo_advection_builder_init(const cs_cdo_connect_t      *connect,
                              const cs_equation_param_t   *eqp,
                              bool                         do_diffusion)
{
  int  n_sysc = 0;
  cs_cdo_adv_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdo_adv_t);

  b->n_i_faces = connect->f_info->n_i_elts;
  b->with_diffusion = do_diffusion;
  b->tmp_vect = NULL;
  b->tmp_scal = NULL;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      const cs_lnum_t  n_max_ec = connect->n_max_ebyc;
      const cs_lnum_t  n_max_vc = connect->n_max_vbyc;

      n_sysc = n_max_vc;

      int  s_size = 2*n_max_ec;
      BFT_MALLOC(b->tmp_scal, s_size, cs_real_t);
      for (int i = 0; i < s_size; i++)
        b->tmp_scal[i] = 0;

    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      const cs_lnum_t  n_max_vc = connect->n_max_vbyc;

      n_sysc = n_max_vc + 1;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid numerical scheme for advection."));
    break;

  } // switch on space_scheme

  b->loc = cs_locmat_create(n_sysc);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Destroy a builder structure for the convection operator
 *
 * \param[in, out] b   pointer to a cs_cdo_adv_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_adv_t *
cs_cdo_advection_builder_free(cs_cdo_adv_t  *b)
{
  if (b == NULL)
    return b;

  BFT_FREE(b->tmp_scal);
  BFT_FREE(b->tmp_vect);

  b->loc = cs_locmat_free(b->loc);

  BFT_FREE(b);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      diffmat   tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_advection_build_local(const cs_cell_mesh_t       *cm,
                               const cs_equation_param_t  *eqp,
                               const cs_real_33_t          diffmat,
                               cs_cdo_adv_t               *b)
{
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB); // Sanity check

  /* Initialize local matrix structure */
  b->loc->n_ent = cm->n_vc;

  for (short int v = 0; v < cm->n_vc; v++)
    b->loc->ids[v] = cm->v_ids[v];

  for (short int v = 0; v < cm->n_vc*cm->n_vc; v++)
    b->loc->val[v] = 0;

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = b->tmp_scal;

  cs_advection_field_get_flux_dfaces(cm->c_id,
                                     eqp->advection_info,
                                     eqp->advection_field,
                                     fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  cs_real_t  *upwcoef = b->tmp_scal + cm->n_ec;

  if (eqp->advection_info.scheme != CS_PARAM_ADVECTION_SCHEME_CENTERED) {

    if (b->with_diffusion) {

      cs_real_3_t  matnu;

      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_nvec3_t  dfq = cm->dface[e];
        const double  mean_flux = fluxes[e]/dfq.meas;

        cs_math_33_3_product((const cs_real_t (*)[3])diffmat, dfq.unitv, matnu);

        cs_real_t  diff_contrib = _dp3(dfq.unitv, matnu);
        if (diff_contrib > cs_math_zero_threshold)
          upwcoef[e] = cm->edge[e].meas * mean_flux / diff_contrib;
        else
          upwcoef[e] = mean_flux * cs_math_big_r; // dominated by convection

      } // Loop on cell edges

    }
    else
      for (short int e = 0; e < cm->n_ec; e++)
        upwcoef[e] = 1/cm->dface[e].meas * fluxes[e];

  }

  /* Build the local convection operator */
  switch (eqp->advection_info.formulation) {

  case CS_PARAM_ADVECTION_FORM_NONCONS:
    _build_local_epcd(cm, eqp->advection_info, b);
    break;

  case CS_PARAM_ADVECTION_FORM_CONSERV:
    _build_local_vpfd(cm, eqp->advection_info, b);
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
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      diffmat   tensor related to the diffusion property
 * \param[in, out] b         pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovcb_advection_build_local(const cs_cell_mesh_t       *cm,
                                const cs_equation_param_t  *eqp,
                                const cs_real_33_t          diffmat,
                                cs_cdo_adv_t               *b)
{
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB); // Sanity check

  /* Initialize local matrix structure */
  int  n_sysc = cm->n_vc + 1;

  for (short int v = 0; v < cm->n_vc; v++)
    b->loc->ids[v] = cm->v_ids[v];
  b->loc->ids[cm->n_vc] = cm->c_id;
  b->loc->n_ent = n_sysc;

  for (short int v = 0; v < n_sysc*n_sysc; v++)
    b->loc->val[v] = 0;

  return b->loc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] b         pointer to a convection builder structure
 * \param[in, out] ls        cell-wise structure sotring the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_advection_add_bc(const cs_cell_mesh_t       *cm,
                          const cs_equation_param_t  *eqp,
                          cs_cdo_adv_t               *b,
                          cs_cdo_locsys_t            *ls)
{
  cs_nvec3_t  adv_vec;

  cs_real_t  *tmp_rhs = b->tmp_scal;
  cs_locmat_t  *m = b->loc; /* Use this local dense matrix as an array
                               taking into account the diagonal contrib. */

  const cs_adv_field_t  *adv_field = eqp->advection_field;
  const cs_param_advection_t  a_info = eqp->advection_info;

  /* Reset local temporay RHS and diagonal contributions */
  for (short int v = 0; v < cm->n_vc; v++) {
    m->val[v] = 0;
    tmp_rhs[v] = 0;
  }

  /* Loop on border faces.
     Add diagonal term for vertices attached to a boundary face where
     the advection field points inward */
  if (cs_advection_field_is_cellwise(adv_field)) {

    for (short int f = 0; f < cm->n_fc; f++) {

      if (cm->f_ids[f] >= b->n_i_faces) { // Border face

        const cs_quant_t  pfq = cm->face[f];

        /* Retrieve the value of the advection field in the current cell */
        cs_advection_field_get_cell_vector(cm->c_id, adv_field, &adv_vec);

        const double  dp = _dp3(adv_vec.unitv, pfq.unitv);
        if (fabs(dp) > cs_math_zero_threshold) {

          /* Loop on border face edges */
          for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

            const short int  eshft = 2*cm->f2e_ids[i];
            const short int  v1 = cm->e2v_ids[eshft];
            const short int  v2 = cm->e2v_ids[eshft+1];
            const double  surf = 0.5 * cs_math_surftri(cm->xv + 3*v1,
                                                       cm->xv + 3*v2,
                                                       pfq.center);
            const double  flx = dp * adv_vec.meas * surf;

            if (dp < 0) { // advection field is inward w.r.t. the face normal

              tmp_rhs[v1] -= flx * ls->dir_bc[v1];
              tmp_rhs[v2] -= flx * ls->dir_bc[v2];

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

    for (short int f = 0; f < cm->n_fc; f++) {

      if (cm->f_ids[f] >= b->n_i_faces) { // Border face

        /* Loop on border face edges */
        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

          const short int  e = cm->f2e_ids[i];
          const short int  v1 = cm->e2v_ids[2*e];
          const short int  v2 = cm->e2v_ids[2*e+1];
          const double  flx_tef = cs_advection_field_get_flux_tef(adv_field,
                                                                  a_info,
                                                                  cm,
                                                                  f, e, v1, v2);
          /* Assume that flx_tef has the same level of contribution for the
             two vertices */
          const double  flx = 0.5 * flx_tef;


          if (flx < 0) { // advection field is inward w.r.t. the face normal

            tmp_rhs[v1] -= flx * ls->dir_bc[v1];
            tmp_rhs[v2] -= flx * ls->dir_bc[v2];

            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS) {
              m->val[v1] -= flx;
              m->val[v2] -= flx;
            }

          }
          else { // advection is oriented outward

            if (a_info.formulation == CS_PARAM_ADVECTION_FORM_CONSERV) {
              m->val[v1] += flx;
              m->val[v1] += flx;

            }

          } // outward

        } // Loop on face edges

      } // Loop on border faces

    } // Loop on cell faces

  } // Advection field uniform or not

  /* Update the diagonal and the RHS of the local system matrix */
  for (short int v = 0; v < cm->n_vc; v++)  {
    ls->mat->val[v*cm->n_vc + v] += m->val[v];
    ls->rhs[v] += tmp_rhs[v];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator with CDO
 *          V+C schemes
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] b         pointer to a convection builder structure
 * \param[in, out] ls        cell-wise structure sotring the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_advection_add_bc(const cs_cell_mesh_t       *cm,
                           const cs_equation_param_t  *eqp,
                           cs_cdo_adv_t               *b,
                           cs_cdo_locsys_t            *ls)
{
  // TODO

  return;
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
cs_cdo_advection_get_upwind_coef_cell(const cs_cdo_quantities_t   *cdoq,
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
