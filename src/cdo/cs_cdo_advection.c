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

#include "cs_cdo_scheme_geometry.h"
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
  double        *tmp_scal;

  cs_locmat_t   *loc;       /* Local matrix for the convection operator */
  cs_locmat_t   *f_loc;     /* Local matrix reduced to face vertices (Only
                               used for V+C algo. */

};

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

// Advanced developper parameters (stabilization coefficient)
static double  cs_cip_stab_coef = 1e-2;

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
 * \brief  Retrieve the face edge id from a given cell edge id.
 *         Compute the weights related to each vertex of face from pre-computed
 *         quantities
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      e        edge id in the cell numbering
 * \param[in]      f        face id in the cell numbering
 * \param[in]      tef      value of the area of the triangles tef
 * \param[in, out] wvf      weights of vertices for the current face
 *
 * \return the face edge id
 */
/*----------------------------------------------------------------------------*/

inline static short int
_set_fquant(const cs_cell_mesh_t    *cm,
            short int                e,
            short int                f,
            const double            *tef,
            double                  *wvf)
{
  short int ef = -1;

  for (short int v = 0; v < cm->n_vc; v++)
    wvf[v] = 0;

  for (short int ee = cm->f2e_idx[f]; ee < cm->f2e_idx[f+1]; ee++) {

    short int  eab = cm->f2e_ids[ee];
    short int  va = cm->e2v_ids[2*eab];
    short int  vb = cm->e2v_ids[2*eab+1];

    wvf[va] += tef[ee];
    wvf[vb] += tef[ee];

    if (eab == e)  ef = ee;

  } // Loop on edges of f1

  const double  invf = 0.5/cm->face[f].meas;
  for (short int v = 0; v < cm->n_vc; v++)
    wvf[v] *= invf;

  assert(ef != -1);
  return ef;
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the consistent part of the convection operator attached to
 *          a cell with a CDO vertex+cell-based scheme when the advection is
 *          considered as constant inside a cell
 *
 * \param[in]      adv_cell  constant vector field inside the current cell
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] b         pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_cellwise_consistent_part(const cs_nvec3_t            adv_cell,
                              const cs_cell_mesh_t       *cm,
                              const cs_face_mesh_t       *fm,
                              cs_cdo_adv_t               *b)
{
  cs_real_3_t  grd_v1, grd_v2, grd_c;

  double  *af = b->f_loc->val;

  /* Temporary buffers */
  double  *bgc_save = b->tmp_scal;            // size = n_fc
  double  *tef_save = b->tmp_scal + cm->n_fc; // size = 2*n_ec

  const int  n_sysf = fm->n_vf + 1;
  const short int  fshift = cm->f2e_idx[fm->f_id];

  /* Useful quantities are stored in b->tmp_scal and b->tmp_vect */
  double  *l_vc = tef_save + 2*cm->n_ec;
  cs_real_3_t  *bgvf = b->tmp_vect + fshift;
  cs_real_3_t  *u_vc = b->tmp_vect + 2*cm->n_ec;

  const cs_nvec3_t  deq = fm->dedge;
  const cs_quant_t  pfq = fm->face;
  const double  hf_coef = cs_math_onethird * _dp3(pfq.unitv,deq.unitv)*deq.meas;
  const double  pfc_vol = hf_coef * pfq.meas;

  /* Set the consistent part for already known part.
     Consistent part (c,c) contribution: sum_(f \in F_c) |pfc| bgc
     Consistent part (i,c) contribution: sum_(f \in F_c) 0.75 * wif * |pfc| bgc
  */

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdc(fm, grd_c);

  const double  bgc = _dp3(grd_c, adv_cell.unitv);
  const double  pfc_bgc = adv_cell.meas * pfc_vol * bgc;

  bgc_save[fm->f_id] = bgc; // Store it for a future use

  af[n_sysf*fm->n_vf + fm->n_vf] = 0.25 * pfc_bgc;         // (c,c)
  for (short int v = 0; v < fm->n_vf; v++)
    af[n_sysf*v + fm->n_vf] = 0.75 * fm->wvf[v] * pfc_bgc; // (i,c)

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute the gradient of Lagrange functions related to the partition
     into tetrahedron of p_(fc) */
  for (short int e = 0; e < fm->n_ef; e++) {

    const double  pef_coef = 0.25 * hf_coef * fm->tef[e] * adv_cell.meas;
    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    const double  bgv1 = _dp3(grd_v1, adv_cell.unitv);
    const double  bgv2 = _dp3(grd_v2, adv_cell.unitv);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity
       grd_f = -( grd_c + grd_v1 + grd_v2) */
    const double  bgf = -(bgc + bgv1 + bgv2);

    for (short int vi = 0; vi < fm->n_vf; vi++) {

      double  *afi = af + n_sysf*vi;
      double  lvci_part = fm->wvf[vi];
      if (vi == v1 || vi == v2)
        lvci_part += 1;
      const double  lvci = pef_coef * lvci_part;

      for (short int vj = 0; vj < fm->n_vf; vj++) {

        double  glj = fm->wvf[vj]*bgf;
        if (vj == v1)
          glj += bgv1;
        else if (vj == v2)
          glj += bgv2;

        afi[vj] += glj * lvci; // consistent part (i,j) face mat.

      } // Loop on vj in Vf

    } // Loop on vi in Vf

    /* (c,j) entries */
    double  *afc = af + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++) {

      double  glj = fm->wvf[vj]*bgf;
      if (vj == v1)
        glj += bgv1;
      else if (vj == v2)
        glj += bgv2;

      afc[vj] += glj * pef_coef; // consistent part (c,j) face mat.

    } // Loop on vj in Vf

    /* Store bgv1, bgv2, bgf */
    bgvf[e][0] = bgv1;
    bgvf[e][1] = bgv2;
    bgvf[e][2] = bgf;

  } // Loop on face edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the consistent part of the convection operator attached to
 *          a cell with a CDO vertex+cell-based scheme when the advection field
 *          is not considered constant inside the current cell
 *
 * \param[in]      adv_field  pointer to a cs_adv_field_t structure
 * \param[in]      adv_cell   constant vector field inside the current cell
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      fm         pointer to a cs_face_mesh_t structure
 * \param[in, out] b          pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_consistent_part(const cs_adv_field_t     *adv_field,
                     const cs_nvec3_t          adv_cell,
                     const cs_cell_mesh_t     *cm,
                     const cs_face_mesh_t     *fm,
                     cs_cdo_adv_t             *b)
{
  cs_real_3_t  grd_v1, grd_v2, grd_c, xg;
  cs_nvec3_t  bnv;
  double  *af = b->f_loc->val;

  /* Temporary buffers */
  double  *bgc_save = b->tmp_scal;            // size = n_fc
  double  *tef_save = b->tmp_scal + cm->n_fc; // size = 2*n_ec

  const int  n_sysf = fm->n_vf + 1;
  const short int  fshift = cm->f2e_idx[fm->f_id];

  /* Useful quantities are stored in b->tmp_scal and b->tmp_vect */
  double  *l_vc = tef_save + 2*cm->n_ec;
  cs_real_3_t  *bgvf = b->tmp_vect + fshift;
  cs_real_3_t  *u_vc = b->tmp_vect + 2*cm->n_ec;

  const cs_nvec3_t  deq = fm->dedge;
  const cs_quant_t  pfq = fm->face;
  const double  hf_coef = cs_math_onethird * _dp3(pfq.unitv,deq.unitv)*deq.meas;

  /* Set the consistent part for already known part.
     Consistent part (c,c) contribution: sum_(f \in F_c) |pfc| bgc
     Consistent part (i,c) contribution: sum_(f \in F_c) 0.75 * wif * |pfc| bgc
  */

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdc(fm, grd_c);

  const double  bgcc = _dp3(grd_c, adv_cell.unitv);

  bgc_save[fm->f_id] = bgcc; // Store it for a future use in the stanilization

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute the gradient of Lagrange functions related to the partition
     into tetrahedron of p_(fc) */
  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    for (int k = 0; k < 3; k++)
      xg[k] = 0.25 * (fm->xv[3*v1+k]+fm->xv[3*v2+k]+fm->xc[k]+pfq.center[k]);
    cs_advection_field_get_at_xyz(adv_field, xg, &bnv);

    const double  bgc = _dp3(grd_c, bnv.unitv);
    const double  pef_coef = 0.25 * bnv.meas * hf_coef * fm->tef[e];

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    const double  bgv1 = _dp3(grd_v1, bnv.unitv);
    const double  bgv2 = _dp3(grd_v2, bnv.unitv);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity
       grd_f = -( grd_c + grd_v1 + grd_v2) */
    const double  bgf = -(bgc + bgv1 + bgv2);

    for (short int vi = 0; vi < fm->n_vf; vi++) {

      double  *afi = af + n_sysf*vi;
      double  lvci_part = fm->wvf[vi];
      if (vi == v1 || vi == v2)
        lvci_part += 1;
      const double  lvci = pef_coef * lvci_part;

      for (short int vj = 0; vj < fm->n_vf; vj++) {

        double  glj = fm->wvf[vj]*bgf;
        if (vj == v1)
          glj += bgv1;
        else if (vj == v2)
          glj += bgv2;

        afi[vj] += glj * lvci; // consistent part (i,j) face mat.

      } // Loop on vj in Vf

      /* (i, c) entries */
      afi[fm->n_vf] += lvci * bgc;

    } // Loop on vi in Vf

    /* (c,j) entries */
    double  *afc = af + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++) {

      double  glj = fm->wvf[vj]*bgf;
      if (vj == v1)
        glj += bgv1;
      else if (vj == v2)
        glj += bgv2;

      afc[vj] += glj * pef_coef; // consistent part (c,j) face mat.

    } // Loop on vj in Vf

    /* (c,c) entries */
    afc[fm->n_vf] += pef_coef * bgc;

    /* Store bgv1, bgv2, bgf */
    bgvf[e][0] = _dp3(grd_v1, adv_cell.unitv);
    bgvf[e][1] = _dp3(grd_v2, adv_cell.unitv);
    /* This formula is a consequence of the Partition of the Unity
       grd_f = -( grd_c + grd_v1 + grd_v2) */
    bgvf[e][2] =  -(bgcc + bgvf[e][0] + bgvf[e][1]);

  } // Loop on face edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stabilization part of the convection operator attached
 *          to a cell with a CDO vertex+cell-based scheme (inside pfc)
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      fm         pointer to a cs_face_mesh_t structure
 * \param[in]      stab_coef  default value of the stabilization coefficient
 * \param[in, out] b          pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_stabilization_part1(const cs_cell_mesh_t        *cm,
                         const cs_face_mesh_t        *fm,
                         const double                 stab_coef,
                         cs_cdo_adv_t                *b)
{
  const short int  n_sysf = fm->n_vf + 1;
  const short int  fshift = cm->f2e_idx[fm->f_id];

  cs_real_3_t  *bgvf = b->tmp_vect + fshift;
  double  *af = b->f_loc->val;

  /* First part of the stabilization (inside each face) */
  for (short int e1 = 0; e1 < fm->n_ef; e1++) {

    short int  e2 = e1 - 1;
    if (e1 == 0)
      e2 = fm->n_ef - 1; // Last edge face

    short int  v1e1 = fm->e2v_ids[2*e1];
    short int  v2e1 = fm->e2v_ids[2*e1+1];
    short int  v1e2 = fm->e2v_ids[2*e2];    // v1_prev
    short int  v2e2 = fm->e2v_ids[2*e2+1];  // v2_prev

    const double  jump_bgf = bgvf[e1][2] - bgvf[e2][2];

    /* Area of the triangle which is the interface between the two pefc */
    double  svfc = stab_coef;
    if (v1e1 == v1e2 || v1e1 == v2e2)
      svfc *= cs_math_surftri(fm->xv + 3*v1e1, fm->face.center, fm->xc);
    else
      svfc *= cs_math_surftri(fm->xv + 3*v2e1, fm->face.center, fm->xc);

    /* Update convection matrix related to the current face
       (symmetric contributions) */
    for (short int vi = 0; vi < fm->n_vf; vi++) {

      double  *afi = af + n_sysf*vi;
      double  jump_i = fm->wvf[vi] * jump_bgf;

      if (vi == v1e1)
        jump_i += bgvf[e1][0];
      else if (vi == v2e1)
        jump_i += bgvf[e1][1];

      if (vi == v1e2)
        jump_i -= bgvf[e2][0];
      else if (vi == v2e2)
        jump_i -= bgvf[e2][1];

      const double  coef_i = svfc * jump_i;

      afi[vi] += coef_i * jump_i; // Stab. (i,i) face mat.

      for (short int vj = vi + 1; vj < fm->n_vf; vj++) {

        double  jump_j = fm->wvf[vj] * jump_bgf;

        if (vj == v1e1)
          jump_j += bgvf[e1][0];
        else if (vj == v2e1)
          jump_j += bgvf[e1][1];

        if (vj == v1e2)
          jump_j -= bgvf[e2][0];
        else if (vj == v2e2)
          jump_j -= bgvf[e2][1];

        const double  coef_ij = coef_i * jump_j;

        afi[vj] += coef_ij;            // Stab. (i,j) face mat.
        af[vj*n_sysf + vi] += coef_ij; // Stab. (j,i) face mat. symmetric

      } // Loop on vj in Vf
    } // Loop on vi in Vf

  } // Loop on edge faces

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stabilization part of the convection operator attached
 *          to a cell with a CDO vertex+cell-based scheme (between pfc)
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      stab_coef  value of the stabilization coefficient
 * \param[in, out] b          pointer to a convection builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_stabilization_part2(const cs_cell_mesh_t        *cm,
                         const double                 stab_coef,
                         cs_cdo_adv_t                *b)
{
  const short int  n_sysc  = cm->n_vc + 1;

  double  *a = b->loc->val;

  /* Temporary buffers used to store pre-computed data */
  double  *bgc_save = b->tmp_scal;            // size = n_fc
  double  *tef_save = b->tmp_scal + cm->n_fc; // size = 2*n_ec
  double  *wvf1 = tef_save + 2*cm->n_ec;
  double  *wvf2 = wvf1 + cm->n_vc;
  cs_real_3_t  *bgvf_save = b->tmp_vect;      // size = 2*n_ec

  for (short int e = 0; e < cm->n_ec; e++) {

    const double  tec = stab_coef * cs_compute_tec(e, cm);
    const short int  eshift = 2*e;
    const short int  v1 = cm->e2v_ids[eshift];
    const short int  v2 = cm->e2v_ids[eshift+1];
    const short int  _v = (v1 < v2) ? 0 : 1;
    const short int  f1 = cm->e2f_ids[eshift];
    const short int  f2 = cm->e2f_ids[eshift+1];
    const double  jump_c = bgc_save[f2] - bgc_save[f1];

    // (c, c) contrib
    double  *ac = a + cm->n_vc*n_sysc;
    ac[cm->n_vc] += tec * jump_c * jump_c;

    const short int ef1 = _set_fquant(cm, e, f1, tef_save, wvf1);
    const short int ef2 = _set_fquant(cm, e, f2, tef_save, wvf2);

    for (short int vi = 0; vi < cm->n_vc; vi++) {
      if (wvf2[vi] + wvf1[vi] > 0) { // vi belongs at least to f1 or f2

        double  *ai = a + vi*n_sysc;
        double  jump_i =
          wvf2[vi]*bgvf_save[ef2][2] - wvf1[vi]*bgvf_save[ef1][2];

        if (vi == v1)
          jump_i += bgvf_save[ef2][_v] - bgvf_save[ef1][_v];
        else if (vi == v2)
          jump_i += bgvf_save[ef2][1-_v] - bgvf_save[ef1][1-_v];

        const double  coef_i = tec * jump_i;

        // (i, i) contrib
        ai[vi] += jump_i * coef_i;

        for (short int vj = vi + 1; vj < cm->n_vc; vj++) {
          if (wvf2[vj] + wvf1[vj] > 0) { // vj belongs at least to f1 or f2

            double  jump_j =
              wvf2[vj]*bgvf_save[ef2][2] - wvf1[vj]*bgvf_save[ef1][2];
            if (vj == v1)
              jump_j += bgvf_save[ef2][_v] - bgvf_save[ef1][_v];
            else if (vj == v2)
              jump_j += bgvf_save[ef2][1-_v] - bgvf_save[ef1][1-_v];

            const double coef_ij = coef_i * jump_j;
            ai[vj] += coef_ij;
            a[vj*n_sysc + vi] += coef_ij;  // symmetric

          } // vj belongs to f1 or f2
        } // vj loop

        // (i, c) contrib
        const double coef_ic = coef_i * jump_c;
        ai[cm->n_vc] += coef_ic;
        ac[vi] += coef_ic;                 // symmetric

      } // vi belongs to f1 or f2
    } // vi loop

  } // End of loop on cell edges

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the stabilization coefficient used in CIP scheme
 *
 * \param[in]  new_value     value of the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_set_cip_coef(double     new_value)
{
  cs_cip_stab_coef = new_value;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the value of the stabilization coefficient used in CIP scheme
 *
 * \return   the value the stabilization coefficient
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_advection_get_cip_coef(void)
{
  return cs_cip_stab_coef;
}

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
  int  n_sysc = 0, s_size = 0, v_size = 0;
  cs_cdo_adv_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdo_adv_t);

  b->n_i_faces = connect->f_info->n_i_elts;
  b->with_diffusion = do_diffusion;

  b->f_loc = NULL;
  b->tmp_vect = NULL;
  b->tmp_scal = NULL;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    n_sysc = connect->n_max_vbyc;
    s_size = 2*connect->n_max_ebyc;
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    n_sysc = connect->n_max_vbyc + 1;
    s_size = 2*(connect->n_max_ebyc+connect->n_max_vbyc) + connect->n_max_fbyc;
    v_size = 2*connect->n_max_ebyc + connect->n_max_vbyf;

    b->f_loc = cs_locmat_create(connect->n_max_vbyf + 1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid numerical scheme for advection."));
    break;

  } // switch on space_scheme

  b->loc = cs_locmat_create(n_sysc);

  if (s_size > 0) { // Temporary scalar buffers

    BFT_MALLOC(b->tmp_scal, s_size, double);
    for (int i = 0; i < s_size; i++)
      b->tmp_scal[i] = 0;

  }

  if (v_size > 0) { // Temporary vector buffers

    BFT_MALLOC(b->tmp_vect, v_size, cs_real_3_t);
    for (int i = 0; i < v_size; i++)
      b->tmp_vect[i][0] = b->tmp_vect[i][1] = b->tmp_vect[i][2]= 0;

  }

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
  b->f_loc = cs_locmat_free(b->f_loc);

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
cs_cdovb_advection_build(const cs_cell_mesh_t       *cm,
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
 * \param[in, out] b         pointer to a convection builder structure
 *
 * \return a pointer to a local dense matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovcb_advection_build(const cs_cell_mesh_t       *cm,
                          const cs_equation_param_t  *eqp,
                          cs_cdo_adv_t               *b)
{
  cs_nvec3_t  adv_cell;

  cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(0); // Modify for OpenMP
  double  *af = b->f_loc->val;
  double  *a = b->loc->val;

  /* Temporary buffers */
  double  *tef_save = b->tmp_scal + cm->n_fc; // size = 2*n_ec

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB);
  assert(eqp->advection_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS);
  assert(eqp->advection_info.scheme == CS_PARAM_ADVECTION_SCHEME_CIP);

  const bool  is_cellwise =
    cs_advection_field_is_cellwise(eqp->advection_field);

  /* Initialize local matrix structure */
  const int  n_sysc = cm->n_vc + 1;
  for (short int v = 0; v < n_sysc*n_sysc; v++) a[v] = 0;
  for (short int v = 0; v < cm->n_vc; v++)
    b->loc->ids[v] = cm->v_ids[v];
  b->loc->ids[cm->n_vc] = cm->c_id;
  b->loc->n_ent = n_sysc;

  /* Use a cellwise constant approximation of the advection field */
  cs_advection_field_get_cell_vector(cm->c_id,
                                     eqp->advection_field,
                                     &adv_cell);

  if (adv_cell.meas < cs_math_get_machine_epsilon())
    return b->loc;

  /* Stabilization coefficent * |beta_c| */
  double  stab_coef = cs_cip_stab_coef * adv_cell.meas;

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_face_mesh_build_from_cell_mesh(cm, f, fm);

    /* Initialize local face matrix. */
    const int  n_sysf = fm->n_vf + 1;
    for (short int i = 0; i < n_sysf*n_sysf; i++) af[i] = 0.;

    /* Store tef areas for a future usage (Second part of the stabilization) */
    const short int  fshift = cm->f2e_idx[f];
    double  *tef = tef_save + fshift;
    for (short int e = 0; e < fm->n_ef; e++) tef[e] = fm->tef[e];

    /* Initialize and update the face matrix inside */
    if (is_cellwise)
      _vcb_cellwise_consistent_part(adv_cell, cm, fm, b);
    else
      _vcb_consistent_part(eqp->advection_field, adv_cell, cm, fm, b);

    /* Build the first part of the stabilization. (inside pfc) */
    _vcb_stabilization_part1(cm, fm, stab_coef, b);

    /* Reorder v1 and v2 to insure coherency for edges shared between
       two faces */
    cs_real_3_t  *bgvf = b->tmp_vect + fshift;
    for (short int e = 0; e < fm->n_ef; e++) {
      if (fm->v_ids[fm->e2v_ids[2*e]] > fm->v_ids[fm->e2v_ids[2*e+1]]) {
        double  save = bgvf[e][0];
        bgvf[e][0] = bgvf[e][1];
        bgvf[e][1] = save;
      }
    }

    /* Add the face matrix to the cell matrix */
    for (short int vi = 0; vi < fm->n_vf; vi++) {

      double  *aci = a + n_sysc*fm->v_ids[vi], *afi = af + n_sysf*vi;
      for (short int vj = 0; vj < fm->n_vf; vj++)
        aci[fm->v_ids[vj]] += afi[vj];  // (i,j) face --> cell
      aci[cm->n_vc] += afi[fm->n_vf];   // (i,c) face --> cell

    }

    double  *acc = a + n_sysc*cm->n_vc, *afc = af + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++)
      acc[fm->v_ids[vj]] += afc[vj];     // (c,j) face --> cell
    acc[cm->n_vc] += afc[fm->n_vf];

  } /* Loop on cell faces */

  /* Build the second part of the stabilization. */
  _vcb_stabilization_part2(cm, stab_coef, b);

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
  cs_nvec3_t  adv_vec;

  cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(0); // Modify for OpenMP

  const int  n_sysc = cm->n_vc + 1;
  const cs_adv_field_t  *adv_field = eqp->advection_field;

  /* Auxiliary buffers to build this operator */
  double  *rhsf = b->tmp_scal;
  double  *dirf = b->tmp_scal + fm->n_max_vbyf;

  double  *mf = b->f_loc->val; /* Local dense matrix on face vertices */
  cs_locmat_t  *m = ls->mat;   /* Local dense matrix (cell + vertices)
                                  storing the cellwise view of the system */

  /* Sanity check */
  assert(eqp->advection_info.formulation == CS_PARAM_ADVECTION_FORM_NONCONS);
  assert(eqp->advection_info.scheme == CS_PARAM_ADVECTION_SCHEME_CIP);

  if (cs_advection_field_is_cellwise(adv_field)) {

    /* Retrieve the value of the advection field in the current cell */
    cs_advection_field_get_cell_vector(cm->c_id, adv_field, &adv_vec);

    for (short int f = 0; f < cm->n_fc; f++) {

      if (cm->f_ids[f] >= b->n_i_faces) { // Border face

        const cs_quant_t  pfq = cm->face[f];
        const double  dp = _dp3(adv_vec.unitv, pfq.unitv);
        const double  beta_nf = adv_vec.meas * 0.5 * (fabs(dp) - dp);

        if (beta_nf > 0) {

          cs_face_mesh_build_from_cell_mesh(cm, f, fm);
          b->f_loc->n_ent = fm->n_vf;

          /* Initialize local face matrix. */
          for (short int i = 0; i < fm->n_vf*fm->n_vf; i++) mf[i] = 0.;

          /* Compute a surfacic Hodge-like operator */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

            double  *mi = mf + vfi*fm->n_vf;

            /* Default contribution */
            const double  default_coef = 0.5 * fm->wvf[vfi] * pfq.meas;
            for (short int vfj = 0; vfj < fm->n_vf; vfj++)
              mi[vfj] = default_coef * fm->wvf[vfj];

            /* Specific diagonal contribution */
            mi[vfi] += 2 * default_coef * cs_math_onethird;

          } // Loop on face vertices

          /* Specific extra-diag contribution */
          for (short int e = 0; e < fm->n_ef; e++) {

            const short int  vf1 = fm->e2v_ids[2*e];
            const short int  vf2 = fm->e2v_ids[2*e+1];
            const double  wtef = cs_math_onetwelve * fm->tef[e];

            mf[vf1*fm->n_vf + vf2] += wtef;
            mf[vf2*fm->n_vf + vf1] += wtef;

          } /* Loop on face edges */

          /* Compute the Dirichlet contribution to RHS */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++)
            dirf[vfi] = beta_nf * ls->dir_bc[fm->v_ids[vfi]];
          cs_locmat_matvec(b->f_loc, dirf, rhsf);

          /* Update RHS */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++)
            ls->rhs[fm->v_ids[vfi]] += rhsf[vfi];

          /* Update cellwise matrix */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

            double  *mi = m->val + fm->v_ids[vfi]*n_sysc;
            double  *mfi = mf + vfi*fm->n_vf;

            for (short int vfj = 0; vfj < fm->n_vf; vfj++)
              mi[fm->v_ids[vfj]] += beta_nf * mfi[vfj];

          } // Loop on face vertices

        } // beta_nf > 0

      } // Border face

    } // Loop on cell faces

  } // Cellwise constant advection field
  else if (cs_advection_field_get_deftype(adv_field) ==
           CS_PARAM_DEF_BY_ANALYTIC_FUNCTION) {

    cs_nvec3_t  beta;
    cs_real_3_t  xg;

    for (short int f = 0; f < cm->n_fc; f++) {

      if (cm->f_ids[f] >= b->n_i_faces) { // Border face

        const cs_quant_t  pfq = cm->face[f];

        cs_face_mesh_build_from_cell_mesh(cm, f, fm);
        b->f_loc->n_ent = fm->n_vf;

        /* Initialize local face matrix. */
        for (short int i = 0; i < fm->n_vf*fm->n_vf; i++) mf[i] = 0.;

        for (short int e = 0; e < fm->n_ef; e++) {

          const short int  vf1 = fm->e2v_ids[2*e];
          const short int  vf2 = fm->e2v_ids[2*e+1];

          for (int k = 0; k < 3; k++)
            xg[k] = cs_math_onethird *
              (pfq.center[k] + fm->xv[3*vf1 + k] + fm->xv[3*vf2 + k]);

          cs_advection_field_get_at_xyz(adv_field, xg, &beta);

          const double  dp = _dp3(beta.unitv, pfq.unitv);
          const double  beta_nef = beta.meas * 0.5 * (fabs(dp) - dp);

          if (beta_nef > 0) {

            const double  wtef = cs_math_onetwelve * fm->tef[e] * beta_nef;

            /* Compute a surfacic Hodge-like operator */
            for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

              double  *mi = mf + vfi*fm->n_vf;

              const bool i_in_e = (vfi == vf1 || vfi == vf2) ? true : false;

              double coef_ii = 2*fm->wvf[vfi]*fm->wvf[vfi];
              if (i_in_e)
                coef_ii += 2*(1 + fm->wvf[vfi]);

              mi[vfi] += coef_ii * wtef;

              for (short int vfj = vfi + 1; vfj < fm->n_vf; vfj++) {

                const bool j_in_e = (vfj == vf1 || vfj == vf2) ? true : false;

                double  coef_ij = 2*fm->wvf[vfi]*fm->wvf[vfj];
                if (i_in_e)
                  coef_ij += fm->wvf[vfj];
                if (j_in_e)
                  coef_ij += fm->wvf[vfi];
                if (i_in_e && j_in_e)
                  coef_ij += 1;
                coef_ij *= wtef; // tef/12 * betanm

                mi[vfj] += coef_ij;
                mf[vfj*fm->n_vf + vfi] += coef_ij;

              } // Loop on face vertices (j)

            } // Loop on face vertices (i)

          } /* beta_nef > 0 */

        } /* Loop on face edges */

        /* Compute the Dirichlet contribution to RHS */
        for (short int v = 0; v < fm->n_vf; v++)
          dirf[v] = ls->dir_bc[fm->v_ids[v]];
        cs_locmat_matvec(b->f_loc, dirf, rhsf);

        /* Update cellwise matrix and RHS */
        for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

          const short int  vci = fm->v_ids[vfi];
          double  *mi = m->val + vci*n_sysc;
          double  *mfi = mf + vfi*fm->n_vf;

          ls->rhs[vci] += rhsf[vfi];
          for (short int vfj = 0; vfj < fm->n_vf; vfj++)
            mi[fm->v_ids[vfj]] += mfi[vfj];

        } // Loop on face vertices

      } // Border face

    } // Loop on cell faces

  }
  else {

    for (short int f = 0; f < cm->n_fc; f++) {

      if (cm->f_ids[f] >= b->n_i_faces) { // Border face

        const cs_quant_t  pfq = cm->face[f];

        /* Retrieve the value of the advection field in the current cell */
        cs_advection_field_get_at_xyz(adv_field, pfq.center, &adv_vec);

        cs_face_mesh_build_from_cell_mesh(cm, f, fm);

        const double  dp = _dp3(adv_vec.unitv, pfq.unitv);
        const double  beta_nf = adv_vec.meas * 0.5 * (fabs(dp) - dp);

        if (beta_nf > 0) {

          b->f_loc->n_ent = fm->n_vf;

          /* Initialize local face matrix. */
          for (short int i = 0; i < fm->n_vf*fm->n_vf; i++) mf[i] = 0.;

          /* Compute a surfacic Hodge-like operator */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

            double  *mi = mf + vfi*fm->n_vf;

            /* Default contribution */
            const double  default_coef = 0.5 * fm->wvf[vfi] * pfq.meas;
            for (short int vfj = 0; vfj < fm->n_vf; vfj++)
              mi[vfj] = default_coef * fm->wvf[vfj];

            /* Specific diagonal contribution */
            mi[vfi] += 2 * default_coef * cs_math_onethird;

          } // Loop on face vertices

          /* Specific extra-diag contribution */
          for (short int e = 0; e < fm->n_ef; e++) {

            const short int  vf1 = fm->e2v_ids[2*e];
            const short int  vf2 = fm->e2v_ids[2*e+1];
            const double  wtef = cs_math_onetwelve * fm->tef[e];

            mf[vf1*fm->n_vf + vf2] += wtef;
            mf[vf2*fm->n_vf + vf1] += wtef;

          } /* Loop on face edges */

          /* Compute the Dirichlet contribution to RHS */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++)
            dirf[vfi] = beta_nf * ls->dir_bc[fm->v_ids[vfi]];
          cs_locmat_matvec(b->f_loc, dirf, rhsf);

          /* Update cellwise matrix and RHS */
          for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

            const short int  vci = fm->v_ids[vfi];
            double  *mi = m->val + vci*n_sysc;
            double  *mfi = mf + vfi*fm->n_vf;

            ls->rhs[vci] += rhsf[vfi];
            for (short int vfj = 0; vfj < fm->n_vf; vfj++)
              mi[fm->v_ids[vfj]] += beta_nf * mfi[vfj];

          } // Loop on face vertices

        } // beta_nf > 0

      } // Border face

    } // Loop on cell faces

  } // Advection field is not uniform on this cell

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
