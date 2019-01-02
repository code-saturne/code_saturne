/*============================================================================
 * Build discrete convection operators for CDO schemes
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
#include <limits.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_bc.h"
#include "cs_hodge.h"
#include "cs_math.h"
#include "cs_scheme_geometry.h"

#if defined(DEBUG) && !defined(NDEBUG) /* For debugging purpose */
#include "cs_dbg.h"
#include "cs_log.h"
#endif

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value of the weighting function related to upwinding
 *
 * \param[in]  criterion  dot product between advection and normal vectors or
 *                        estimation of a local Peclet number
 *
 * \return the weight value
 */
/*----------------------------------------------------------------------------*/

typedef double
(_upwind_weight_t)(double   criterion);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the system by taking into account the boundary conditions
 *
 * \param[in]      csys        pointer to a cellwise view of the system
 * \param[in]      bf          local face number id (0..n_bc_faces)
 * \param[in]      flx         advective flux across the triangle "tef"
 * \param[in]      v1          first vertex to consider in the cell numbering
 * \param[in]      v2          second vertex to consider in the cell numbering
 * \param[in, out] rhs         rhs of the local system
 * \param[in, out] diag        diagonal of the local system
 */
/*----------------------------------------------------------------------------*/

typedef void
(_update_vb_system_with_bc_t)(const cs_cell_sys_t         *csys,
                              const short int              bf,
                              const double                 flx,
                              const short int              v1,
                              const short int              v2,
                              double                      *rhs,
                              double                      *diag);

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

/* Advanced developper parameters (stabilization coefficient) */
static double  cs_cip_stab_coef = 1e-2;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the weight of upwinding for an upwind scheme
 *
 * \param[in]  criterion  dot product between advection and normal vectors or
 *                        estimation of a local Peclet number
 *
 * \return the weight value
 */
/*----------------------------------------------------------------------------*/

inline static double
_get_upwind_weight(double  criterion)
{
  if (criterion > cs_math_zero_threshold)
    return  1;
  else if (criterion < -cs_math_zero_threshold)
    return  0;
  else
    return  0.5;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the weight of upwinding for an Samarskii scheme
 *
 * \param[in]  criterion  dot product between advection and normal vectors or
 *                        estimation of a local Peclet number
 *
 * \return the weight value
 */
/*----------------------------------------------------------------------------*/

inline static double
_get_samarskii_weight(double  criterion)
{
  double  weight = 1.0;
  if (criterion < 0)
    weight = 1./(2 - criterion);
  else
    weight = (1 + criterion)/(2 + criterion);

  return weight;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the weight of upwinding for an Sharfetter-Gummel scheme
 *
 * \param[in]  criterion  dot product between advection and normal vectors or
 *                        estimation of a local Peclet number
 *
 * \return the weight value
 */
/*----------------------------------------------------------------------------*/

inline static double
_get_sg_weight(double  criterion)
{
  double  weight = 1.0;

  if (criterion < 0)
    weight = 0.5*exp(criterion);
  else
    weight = 1 - 0.5*exp(-criterion);

  return weight;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the function related to the kind of upwind weighting
 *
 * \param[in]  scheme    schme used for discretizing the advection term
 *
 * \return  a function pointer
 */
/*----------------------------------------------------------------------------*/

inline static _upwind_weight_t *
_assign_weight_func(const cs_param_advection_scheme_t    scheme)
{
  switch (scheme) {

  case CS_PARAM_ADVECTION_SCHEME_UPWIND:
    return _get_upwind_weight;

  case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
    return _get_samarskii_weight;

  case CS_PARAM_ADVECTION_SCHEME_SG:  /* Sharfetter-Gummel */
    return _get_sg_weight;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible type of algorithm to compute the weight of"
              " upwind.");

  }  /* Switch on the type of function to return */

  return NULL;  /* Avoid warning */
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

  }  /* Loop on edges of f1 */

  const double  invf = 0.5/cm->face[f].meas;
  for (short int v = 0; v < cm->n_vc; v++)
    wvf[v] *= invf;

  assert(ef != -1);
  return ef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Assemble a face matrix into a cell matrix
 *
 * \param[in]      cm    pointer to a cs_cell_mesh_t structure
 * \param[in]      fm    pointer to a cs_face_mesh_t structure
 * \param[in]      af    pointer to a cs_sdm_t structure related to a face
 * \param[in, out] ac    pointer to a cs_sdm_t structure related to a cell
 */
/*----------------------------------------------------------------------------*/

static void
_assemble_face(const cs_cell_mesh_t     *cm,
               const cs_face_mesh_t     *fm,
               const cs_sdm_t           *af,
               cs_sdm_t                 *ac)
{
  /* Add the face matrix to the cell matrix */
  for (short int vi = 0; vi < fm->n_vf; vi++) {

    const double *afi = af->val + af->n_rows*vi;

    double  *aci = ac->val + ac->n_rows*fm->v_ids[vi];
    for (short int vj = 0; vj < fm->n_vf; vj++)
      aci[fm->v_ids[vj]] += afi[vj];   /* (i,j) face --> cell */
    aci[cm->n_vc] += afi[fm->n_vf];    /* (i,c) face --> cell */

  }

  const double  *afc = af->val + af->n_rows*fm->n_vf;

  double  *acc = ac->val + ac->n_rows*cm->n_vc;
  for (short int vj = 0; vj < fm->n_vf; vj++)
    acc[fm->v_ids[vj]] += afc[vj];      /* (c,j) face --> cell */
  acc[cm->n_vc] += afc[fm->n_vf];       /* (c,c) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Conservative formulation and Upwind flux)
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      get_weight  pointer to a function for evaluating upw. weight
 * \param[in]      fluxes      array of fluxes on dual faces
 * \param[in]      upwcoef     array of coefficient to the weight the upwinding
 * \param[in, out] adv         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_epcd_upw(const cs_cell_mesh_t      *cm,
                     _upwind_weight_t          *get_weight,
                     const cs_real_t            fluxes[],
                     const cs_real_t            upwcoef[],
                     cs_sdm_t                  *adv)
{
  /* Loop on cell edges */
  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  sgn_v1 = cm->e2v_sgn[e];
    const cs_real_t  beta_flx = fluxes[e] * sgn_v1;

    if (fabs(beta_flx) > 0) {

      /* Compute the upwind coefficient knowing that fd(e),cd(v) = -v,e */
      const double  wv1 = get_weight(-sgn_v1 * upwcoef[e]);
      const double  c1mw = beta_flx * (1 - wv1);
      const double  cw = beta_flx * wv1;

      /* Update local convection matrix.
         Remember that sgn_v2 = -sgn_v1 */
      short int  v1 = cm->e2v_ids[2*e];
      short int  v2 = cm->e2v_ids[2*e+1];
      assert(v1 != -1 && v2 != -1);  /* Sanity check */

      /* Update for the vertex v1 */
      double  *m1 = adv->val + v1*adv->n_rows;
      m1[v1] +=  c1mw;           /* diagonal term */
      m1[v2] =  -c1mw;           /* extra-diag term */


      /* Update for the vertex v2 */
      double  *m2 = adv->val + v2*adv->n_rows;
      m2[v2] += -cw;           /* diagonal term */
      m2[v1] =   cw;           /* extra-diag term */

    } /* convective flux is greater than zero */

  } /* Loop on cell edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Non-conservative formulation and centered flux)
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes      array of fluxes on dual faces
 * \param[in, out] adv         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_epcd_cen(const cs_cell_mesh_t          *cm,
                     const cs_real_t                fluxes[],
                     cs_sdm_t                      *adv)
{
  /* Weight is always equal to 0.5 Loop on cell edges */
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_real_t  wflx = 0.5 * fluxes[e] * cm->e2v_sgn[e]; // sgn_v1

    if (fabs(wflx) > 0) {

      /* Update local convection matrix */
      short int  v1 = cm->e2v_ids[2*e];
      short int  v2 = cm->e2v_ids[2*e+1];
      assert(v1 != -1 && v2 != -1); /* Sanity check */

      /* Update for the vertex v1 */
      double  *adv1 = adv->val + v1*adv->n_rows;
      adv1[v1] +=  wflx;           /* diagonal term */
      adv1[v2] =  -wflx;           /* extra-diag term */

      /* Update for the vertex v2
         Use the fact that fd(e),cd(v) = -v,e and sgn_v1 = -sgn_v2 */
      double  *adv2 = adv->val + v2*adv->n_rows;
      adv2[v2] += -wflx;           /* diagonal term */
      adv2[v1] =   wflx;           /* extra-diag term */

    } /* convective flux is greater than zero */

  } /* Loop on cell edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Conservative formulation and Upwind flux)
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      get_weight  pointer to a function for evaluating upw. weight
 * \param[in]      fluxes      array of fluxes on dual faces
 * \param[in]      upwcoef     array of coefficient to the weight the upwinding
 * \param[in, out] adv         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_vpfd_upw(const cs_cell_mesh_t          *cm,
                     _upwind_weight_t              *get_weight,
                     const cs_real_t                fluxes[],
                     const cs_real_t                upwcoef[],
                     cs_sdm_t                      *adv)
{
  /* Loop on cell edges */
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_real_t  beta_flx = fluxes[e];

    if (fabs(beta_flx) > 0) {

      short int  sgn_v1 = cm->e2v_sgn[e];

      /* Compute the upwind coefficient knowing that fd(e),cd(v) = -v,e */
      const double  wv1 = get_weight(-sgn_v1 * upwcoef[e]); /* in [0,1] */
      const double  cw1 = sgn_v1 * beta_flx * wv1;
      const double  cw2 = sgn_v1 * beta_flx * (1 - wv1);

      /* Update local convection matrix */
      short int  v1 = cm->e2v_ids[2*e];
      short int  v2 = cm->e2v_ids[2*e+1];
      assert(v1 != -1 && v2 != -1);  /* Sanity check */

      /* Update for the vertex v1 */
      double  *m1 = adv->val + v1*adv->n_rows;
      m1[v1] += -cw1;           /* diagonal term */
      m1[v2] =  -cw2;           /* extra-diag term */

      /* Update for the vertex v2 */
      double  *m2 = adv->val + v2*adv->n_rows;
      m2[v2] +=  cw2;           /* diagonal term */
      m2[v1] =   cw1;           /* extra-diag term */

    } /* convective flux is greater than zero */

  } /* Loop on cell edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Conservative formulation and centered flux)
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes      array of fluxes on dual faces
 * \param[in, out] adv         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_vpfd_cen(const cs_cell_mesh_t          *cm,
                     const cs_real_t                fluxes[],
                     cs_sdm_t                      *adv)
{
  /* Weight is always equal to 0.5. Loop on cell edges */
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_real_t  wflx = 0.5*fluxes[e]*cm->e2v_sgn[e];

    if (fabs(wflx) > 0) { /* Update local convection matrix */

      short int  v1 = cm->e2v_ids[2*e];
      short int  v2 = cm->e2v_ids[2*e+1];
      assert(v1 != -1 && v2 != -1);  /* Sanity check */

      /* Update for the vertex v1 */
      double  *adv1 = adv->val + v1*adv->n_rows;
      adv1[v1] += -wflx;           /* diagonal term */
      adv1[v2] =  -wflx;           /* extra-diag term */

      /* Update for the vertex v2 (sgn_v2 = -sgn_v1) */
      double  *adv2 = adv->val + v2*adv->n_rows;
      adv2[v2] +=  wflx;           /* diagonal term */
      adv2[v1] =   wflx;           /* extra-diag term */

    } /* convective flux is greater than zero */

  } /* Loop on cell edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local convection operator between primal vertices and
 *          dual faces. (Conservative formulation and mixed centered/upwind)
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      upwp     portion of upwinding to apply
 * \param[in]      fluxes   array of fluxes on dual faces
 * \param[in, out] adv      pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_cell_vpfd_mcu(const cs_cell_mesh_t    *cm,
                     const cs_real_t          upwp,
                     const cs_real_t          fluxes[],
                     cs_sdm_t                *adv)
{
  /* Loop on cell edges */
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_real_t  beta_flx = fluxes[e];

    if (fabs(beta_flx) > 0) {

      short int  sgn_v1 = cm->e2v_sgn[e];

      /* Compute the upwind coefficient knowing that fd(e),cd(v) = -v,e */
      const double  wup1 = _get_upwind_weight(-sgn_v1*beta_flx); /* 0 or 1 */
      /* (1-upwp)*0.5 --> centered scheme contribution
         upwp * wup1  --> upwind scheme contribution */
      const double  wv1 = wup1 * upwp + 0.5*(1 - upwp);
      const double  cw1 = sgn_v1 * beta_flx * wv1;
      const double  cw2 = sgn_v1 * beta_flx * (1 - wv1);

      /* Update local convection matrix */
      short int  v1 = cm->e2v_ids[2*e];
      short int  v2 = cm->e2v_ids[2*e+1];
      assert(v1 != -1 && v2 != -1);  /* Sanity check */

      /* Update for the vertex v1 */
      double  *m1 = adv->val + v1*adv->n_rows;
      m1[v1] += -cw1;           /* diagonal term */
      m1[v2] =  -cw2;           /* extra-diag term */

      /* Update for the vertex v2 */
      double  *m2 = adv->val + v2*adv->n_rows;
      m2[v2] +=  cw2;           /* diagonal term */
      m2[v1] =   cw1;           /* extra-diag term */

    } /* convective flux is greater than zero */

  } /* Loop on cell edges */

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the cellwise system with the contributions arising from
 *         the boundary conditions
 *
 * \param[in]      beta_nf   advection field coefficient
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      matf      pointer to a local dense matrix related to a face
 * \param[in, out] csys      pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_vcb_system_with_bc(const cs_real_t              beta_nf,
                           const cs_face_mesh_t        *fm,
                           const cs_sdm_t              *matf,
                           cs_cell_sys_t               *csys)
{
  double  _dirf[10], _rhsf[10];
  double  *dirf = NULL, *rhsf = NULL;

  if (csys->has_dirichlet) {

    if (fm->n_vf > 10) {
      BFT_MALLOC(dirf, 2*fm->n_vf, double);
      rhsf = dirf + fm->n_vf;
    }
    else {
      dirf = _dirf;
      rhsf = _rhsf;
    }

    /* Compute the Dirichlet contribution to RHS */
    for (short int vfi = 0; vfi < fm->n_vf; vfi++)
      dirf[vfi] = beta_nf * csys->dir_values[fm->v_ids[vfi]];
    cs_sdm_square_matvec(matf, dirf, rhsf);

    /* Update RHS */
    for (short int vfi = 0; vfi < fm->n_vf; vfi++)
      csys->rhs[fm->v_ids[vfi]] += rhsf[vfi];

  }  /* There is at least one dirichlet */

  /* Update cellwise matrix */
  const int  n_cell_dofs = csys->mat->n_rows;
  double *matc = csys->mat->val;

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    const double  *mfi = matf->val + vfi*fm->n_vf;
    double  *mi = matc + fm->v_ids[vfi]*n_cell_dofs;

    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {
      const short int  vj = fm->v_ids[vfj];
      mi[vj] += beta_nf * mfi[vfj];
    }

  }  /* Loop on face vertices */

  if (dirf != _dirf)
    BFT_FREE(dirf);
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
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_cellwise_consistent_part(const cs_nvec3_t            adv_cell,
                              const cs_cell_mesh_t       *cm,
                              const cs_face_mesh_t       *fm,
                              cs_cell_builder_t          *cb)
{
  cs_real_3_t  grd_v1, grd_v2, grd_c;

  const int  n_sysf = fm->n_vf + 1;
  const short int  fshift = cm->f2e_idx[fm->f_id];

  /* Temporary buffers storing useful quantities to build the consistent part */
  cs_sdm_t  *af = cb->aux;

  /* tef_save is not set here but its length is 2*n_ec */
  double  *bgc_save = cb->values; // size = n_fc
  double  *l_vc = cb->values + cm->n_fc + 2*cm->n_ec;  // size = n_vc

  cs_real_3_t  *bgvf = cb->vectors + fshift;     // size 2*n_ec
  cs_real_3_t  *u_vc = cb->vectors + 2*cm->n_ec; // size n_vc

  /* Useful quantities are stored in cb->values and cb->vectors */
  const cs_nvec3_t  deq = fm->dedge;
  const cs_quant_t  pfq = fm->face;
  const double  hf_coef = cs_math_1ov3 * cm->hfc[fm->f_id];
  const double  pfc_vol = hf_coef * pfq.meas;

  /* Set the consistent part for already known part.
     Consistent part (c,c) contribution: sum_(f \in F_c) |pfc| bgc
     Consistent part (i,c) contribution: sum_(f \in F_c) 0.75 * wif * |pfc| bgc
  */

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdfc(fm->f_sgn, pfq, deq, grd_c);

  const double  bgc = _dp3(grd_c, adv_cell.unitv);
  const double  pfc_bgc = adv_cell.meas * pfc_vol * bgc;

  bgc_save[fm->f_id] = bgc;  /* Store it for a future use */

  af->val[n_sysf*fm->n_vf + fm->n_vf] = 0.25 * pfc_bgc;         // (c,c)
  for (short int v = 0; v < fm->n_vf; v++)
    af->val[n_sysf*v + fm->n_vf] = 0.75 * fm->wvf[v] * pfc_bgc; // (i,c)

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

      double  *afi = af->val + n_sysf*vi;
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

        afi[vj] += glj * lvci;  /* consistent part (i,j) face mat. */

      }  /* Loop on vj in Vf */

    }  /* Loop on vi in Vf */

    /* (c,j) entries */
    double  *afc = af->val + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++) {

      double  glj = fm->wvf[vj]*bgf;
      if (vj == v1)
        glj += bgv1;
      else if (vj == v2)
        glj += bgv2;

      afc[vj] += glj * pef_coef;  /* consistent part (c,j) face mat. */

    }  /* Loop on vj in Vf */

    /* Store bgv1, bgv2, bgf */
    bgvf[e][0] = bgv1;
    bgvf[e][1] = bgv2;
    bgvf[e][2] = bgf;

  }  /* Loop on face edges */

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
 * \param[in]      t_eval     time at which one evaluates the advection field
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_consistent_part(const cs_adv_field_t     *adv_field,
                     const cs_nvec3_t          adv_cell,
                     const cs_cell_mesh_t     *cm,
                     const cs_face_mesh_t     *fm,
                     cs_real_t                 t_eval,
                     cs_cell_builder_t        *cb)
{
  cs_real_3_t  grd_v1, grd_v2, grd_c, xg;
  cs_nvec3_t  bnv;

  const short int  fshift = cm->f2e_idx[fm->f_id];
  const int  n_sysf = fm->n_vf + 1;
  const cs_nvec3_t  deq = fm->dedge;
  const cs_quant_t  pfq = fm->face;
  const double  hf_coef = cs_math_1ov3 * cm->hfc[fm->f_id];

  /* Useful quantities are stored in cb->values and cb->vectors */
  cs_sdm_t  *af = cb->aux;
  double  *bgc_save = cb->values;                      // size n_fc
  double  *l_vc = cb->values + cm->n_fc + 2*cm->n_ec;  // size n_vc
  cs_real_3_t  *bgvf = cb->vectors + fshift;
  cs_real_3_t  *u_vc = cb->vectors + 2*cm->n_ec;

  /* Set the consistent part for already known part.
     Consistent part (c,c) contribution: sum_(f \in F_c) |pfc| bgc
     Consistent part (i,c) contribution: sum_(f \in F_c) 0.75 * wif * |pfc| bgc
  */

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  const cs_real_t ohf = -fm->f_sgn/cm->hfc[fm->f_id];
  for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];

  const double  bgcc = _dp3(grd_c, adv_cell.unitv);

  bgc_save[fm->f_id] = bgcc;  /* Store it for using it in the stabilization */

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
    cs_advection_field_cw_eval_at_xyz(adv_field, cm, xg, t_eval, &bnv);

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

      double  *afi = af->val + n_sysf*vi;
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

        afi[vj] += glj * lvci;  /* consistent part (i,j) face mat. */

      }  /* Loop on vj in Vf */

      /* (i, c) entries */
      afi[fm->n_vf] += lvci * bgc;

    }  /* Loop on vi in Vf */

    /* (c,j) entries */
    double  *afc = af->val + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++) {

      double  glj = fm->wvf[vj]*bgf;
      if (vj == v1)
        glj += bgv1;
      else if (vj == v2)
        glj += bgv2;

      afc[vj] += glj * pef_coef;  /* consistent part (c,j) face mat. */

    }  /* Loop on vj in Vf */

    /* (c,c) entries */
    afc[fm->n_vf] += pef_coef * bgc;

    /* Store bgv1, bgv2, bgf */
    bgvf[e][0] = _dp3(grd_v1, adv_cell.unitv);
    bgvf[e][1] = _dp3(grd_v2, adv_cell.unitv);
    /* This formula is a consequence of the Partition of the Unity
       grd_f = -( grd_c + grd_v1 + grd_v2) */
    bgvf[e][2] =  -(bgcc + bgvf[e][0] + bgvf[e][1]);

  }  /* Loop on face edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stabilization part of the convection operator attached
 *          to a cell with a CDO vertex+cell-based scheme (inside pfc)
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      fm         pointer to a cs_face_mesh_t structure
 * \param[in]      stab_coef  default value of the stabilization coefficient
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_stabilization_part1(const cs_cell_mesh_t     *cm,
                         const cs_face_mesh_t     *fm,
                         const double              stab_coef,
                         cs_cell_builder_t        *cb)
{
  const short int  n_sysf = fm->n_vf + 1;
  const short int  fshift = cm->f2e_idx[fm->f_id];

  cs_real_3_t  *bgvf = cb->vectors + fshift;
  cs_sdm_t  *af = cb->aux;

  /* First part of the stabilization (inside each face) */
  for (short int e1 = 0; e1 < fm->n_ef; e1++) {

    short int  e2 = e1 - 1;
    if (e1 == 0)
      e2 = fm->n_ef - 1;  /* Last edge face */

    short int  v1e1 = fm->e2v_ids[2*e1];
    short int  v2e1 = fm->e2v_ids[2*e1+1];
    short int  v1e2 = fm->e2v_ids[2*e2];     /* v1_prev */
    short int  v2e2 = fm->e2v_ids[2*e2+1];   /* v2_prev */

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

      double  *afi = af->val + n_sysf*vi;
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

      afi[vi] += coef_i * jump_i;  /* Stab. (i,i) face mat. */

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

        afi[vj] += coef_ij;                 /* Stab. (i,j) face mat. */
        af->val[vj*n_sysf + vi] += coef_ij; /* Stab. (j,i) face mat. symm. */

      }  /* Loop on vj in Vf */
    }  /* Loop on vi in Vf */

  }  /* Loop on edge faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stabilization part of the convection operator attached
 *          to a cell with a CDO vertex+cell-based scheme (between pfc)
 *
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      stab_coef   value of the stabilization coefficient
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_stabilization_part2(const cs_cell_mesh_t     *cm,
                         const double              stab_coef,
                         cs_cell_builder_t        *cb)
{
  const short int  n_sysc  = cm->n_vc + 1;

  cs_sdm_t  *a = cb->loc;

  /* Temporary buffers used to store pre-computed data */
  double  *bgc_save, *tef_save, *wvf1, *wvf2;  /* scalar-valued buffers */
  int  m_shft = 0;

  bgc_save = cb->values, m_shft = cm->n_fc;             // size = n_fc
  tef_save = cb->values + m_shft, m_shft += 2*cm->n_ec; // size = 2*n_ec
  wvf1 = cb->values + m_shft, m_shft += cm->n_vc;       // size = n_vc
  wvf2 = cb->values + m_shft, m_shft += cm->n_vc;       // size = n_vc

  cs_real_3_t  *bgvf_save = cb->vectors;      /* size = 2*n_ec */

  for (short int e = 0; e < cm->n_ec; e++) {

    const double  tec = stab_coef * cs_compute_area_from_quant(cm->edge[e],
                                                               cm->xc);
    const short int  eshift = 2*e;
    const short int  v1 = cm->e2v_ids[eshift];
    const short int  v2 = cm->e2v_ids[eshift+1];
    const short int  _v = (v1 < v2) ? 0 : 1;
    const short int  f1 = cm->e2f_ids[eshift];
    const short int  f2 = cm->e2f_ids[eshift+1];
    const double  jump_c = bgc_save[f2] - bgc_save[f1];

    /* (c, c) contrib */
    double  *ac = a->val + cm->n_vc*n_sysc;
    ac[cm->n_vc] += tec * jump_c * jump_c;

    const short int ef1 = _set_fquant(cm, e, f1, tef_save, wvf1);
    const short int ef2 = _set_fquant(cm, e, f2, tef_save, wvf2);

    for (short int vi = 0; vi < cm->n_vc; vi++) {
      if (wvf2[vi] + wvf1[vi] > 0) {  /* vi belongs at least to f1 or f2 */

        double  *ai = a->val + vi*n_sysc;
        double  jump_i =
          wvf2[vi]*bgvf_save[ef2][2] - wvf1[vi]*bgvf_save[ef1][2];

        if (vi == v1)
          jump_i += bgvf_save[ef2][_v] - bgvf_save[ef1][_v];
        else if (vi == v2)
          jump_i += bgvf_save[ef2][1-_v] - bgvf_save[ef1][1-_v];

        const double  coef_i = tec * jump_i;

        /* (i, i) contrib */
        ai[vi] += jump_i * coef_i;

        for (short int vj = vi + 1; vj < cm->n_vc; vj++) {
          if (wvf2[vj] + wvf1[vj] > 0) {  /* vj belongs at least to f1 or f2 */

            double  jump_j =
              wvf2[vj]*bgvf_save[ef2][2] - wvf1[vj]*bgvf_save[ef1][2];
            if (vj == v1)
              jump_j += bgvf_save[ef2][_v] - bgvf_save[ef1][_v];
            else if (vj == v2)
              jump_j += bgvf_save[ef2][1-_v] - bgvf_save[ef1][1-_v];

            const double coef_ij = coef_i * jump_j;
            ai[vj] += coef_ij;
            a->val[vj*n_sysc + vi] += coef_ij;   /* symmetric */

          }  /* vj belongs to f1 or f2 */
        }  /* vj loop */

        /* (i, c) contrib */
        const double coef_ic = coef_i * jump_c;
        ai[cm->n_vc] += coef_ic;
        ac[vi] += coef_ic;  /* symmetric */

      }  /* vi belongs to f1 or f2 */
    }  /* vi loop */

  }  /* End of loop on cell edges */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
 * \brief   Build the cellwise advection operator for CDO-Fb schemes
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval      time at which one evaluates the advection field
 * \param[in]      build_func  pointer to the function building the system
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_advection_build(const cs_equation_param_t   *eqp,
                         const cs_cell_mesh_t        *cm,
                         cs_real_t                    t_eval,
                         cs_cdofb_advection_t        *build_func,
                         cs_cell_builder_t           *cb)
{
  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOFB);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ));

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, adv);

  /* Compute the flux across the primal faces. Store in cb->adv_fluxes */
  cs_advection_field_cw_face_flux(cm, eqp->adv_field, t_eval, cb->adv_fluxes);

  /* Define the local operator for advection */
  build_func(cm, cb->adv_fluxes, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection fluxes");
    cs_log_printf(CS_LOG_DEFAULT, "\n beta_fluxes>>");
    for (int f = 0; f < cm->n_fc; f++)
      cs_log_printf(CS_LOG_DEFAULT, "f%d;% -5.3e|",
                    cm->f_ids[f], cb->adv_fluxes[f]);
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (without diffusion)
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc(const cs_equation_param_t   *eqp,
                       const cs_cell_mesh_t        *cm,
                       cs_cell_builder_t           *cb,
                       cs_cell_sys_t               *csys)
{
  CS_UNUSED(eqp);

  /* Sanity checks */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));
  assert(csys->n_dofs == cm->n_fc + 1);

  const cs_real_t  *fluxes = cb->adv_fluxes;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    const short int  f = csys->_f_ids[i];
    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];

    cs_real_t  *f_row = csys->mat->val + f*csys->n_dofs;

    if (fabs(beta_flx) > cs_math_zero_threshold) {

      const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);
      const cs_real_t  beta_plus =  0.5*(fabs(beta_flx) + beta_flx);

      /* Outward flux */
      f_row[f] += beta_plus;

      /* Weak enforcement of the Dirichlet BCs. Update RHS for faces attached
         to a boundary face */
      csys->rhs[f] += beta_minus * csys->dir_values[f];

    }
    else { /* If non diffusion, add this term to avoid a singularity in the
              linear system. Set the value at the current face as the mean value
              of the two adjacent cells */

      f_row[cm->n_fc]  = -1.0;
      f_row[f]        +=  1.0;

    }

  } /* Loop on border faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add the contribution of the boundary conditions to the local system
 *          in CDO-Fb schemes (with a diffusion term activated)
 *
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_bc_wdi(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_cell_builder_t           *cb,
                           cs_cell_sys_t               *csys)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));
  assert(csys->n_dofs == cm->n_fc + 1);

  const cs_real_t  *fluxes = cb->adv_fluxes;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    const short int  f = csys->_f_ids[i];
    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];

    cs_real_t  *f_row = csys->mat->val + f*csys->n_dofs;

    const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);
    const cs_real_t  beta_plus =  0.5*(fabs(beta_flx) + beta_flx);

    /* Outward flux */
    if (eqp->adv_formulation == CS_PARAM_ADVECTION_FORM_CONSERV)
      f_row[f] += beta_minus;
    else
      f_row[f] += beta_plus;

    /* Weak enforcement of the Dirichlet BCs. Update RHS for faces attached
       to a boundary face */
    csys->rhs[f] += beta_minus * csys->dir_values[f];

  } /* Loop on border faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the conservative formulation
 *         - upwind scheme
 *         - no diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] adv       pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwcsv(const cs_cell_mesh_t      *cm,
                           const cs_real_t            fluxes[],
                           cs_sdm_t                  *adv)
{
  const short int  c = cm->n_fc;  /* current cell's location in the matrix */

  /* Access the row containing current cell */
  double  *c_row = adv->val + c*adv->n_rows;

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];
    const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);

    /* access the row containing the current face */
    double  *f_row = adv->val + f*adv->n_rows;

    if (fabs(beta_flx) > cs_math_zero_threshold) {

      /* Consistent part */
      f_row[c] -= beta_flx;
      c_row[c] += beta_flx;

      /* Stabilization part */
      f_row[f] += beta_minus;
      f_row[c] -= beta_minus;
      c_row[f] -= beta_minus;
      c_row[c] += beta_minus;

    }
    else { /* If non diffusion, add this term to avoid a singularity in the
              linear system. Set the value at the current face as the mean value
              of the two adjacent cells */

      f_row[c]  = -1.0;
      f_row[f] +=  1.0;

    }

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the conservative formulation
 *         - upwind scheme
 *         - diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] adv       pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwcsv_di(const cs_cell_mesh_t      *cm,
                              const cs_real_t            fluxes[],
                              cs_sdm_t                  *adv)
{
  const short int  c = cm->n_fc;  /* current cell's location in the matrix */

  /* Access the row containing current cell */
  double  *c_row = adv->val + c*adv->n_rows;

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];
    const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);

    /* access the row containing the current face */
    double  *f_row = adv->val + f*adv->n_rows;

    if (fabs(beta_flx) > cs_math_zero_threshold) {

      /* Consistent part */
      f_row[c] -= beta_flx;
      c_row[c] += beta_flx;

      /* Stabilization part */
      f_row[f] += beta_minus;
      f_row[c] -= beta_minus;
      c_row[f] -= beta_minus;
      c_row[c] += beta_minus;

    }

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the non-conservative formulation
 *         - upwind scheme
 *         - no diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] cb        pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwnoc(const cs_cell_mesh_t      *cm,
                           const cs_real_t            fluxes[],
                           cs_sdm_t                  *adv)
{
  const short int  c = cm->n_fc;  /* current cell's location in the matrix */

  /* Access the row containing current cell */
  double  *c_row = adv->val + c*adv->n_rows;

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];
    const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);

    /* access the row containing the current face */
    double  *f_row = adv->val + f*adv->n_rows;

    if (fabs(beta_flx) > cs_math_zero_threshold) {

      /* Consistent part */
      f_row[c] -= beta_flx;

      /* Stabilization part */
      f_row[f] += beta_minus;
      f_row[c] -= beta_minus;
      c_row[f] -= beta_minus;
      c_row[c] += beta_minus;

    }
    else { /* If non diffusion, add this term to avoid a singularity in the
              linear system. Set the value at the current face as the mean value
              of the two adjacent cells */

      f_row[c]  = -1.0;
      f_row[f] +=  1.0;

    }

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the convection operator attached to a cell with a CDO
 *         face-based scheme in the non-conservative formulation
 *         - upwind scheme
 *         - diffusion is present
 *         Rely on the article: Di Pietro, Droniou, Ern (2015)
 *         A discontinuous-skeletal method for advection-diffusion-reaction on
 *         general meshes
 *         The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      fluxes    array of computed fluxes across cell faces
 * \param[in, out] cb        pointer to a local matrix to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_fb_upwnoc_di(const cs_cell_mesh_t      *cm,
                              const cs_real_t            fluxes[],
                              cs_sdm_t                  *adv)
{
  const short int  c = cm->n_fc;  /* current cell's location in the matrix */

  /* Access the row containing current cell */
  double  *c_row = adv->val + c*adv->n_rows;

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_real_t  beta_flx = cm->f_sgn[f]*fluxes[f];
    const cs_real_t  beta_minus = 0.5*(fabs(beta_flx) - beta_flx);

    /* access the row containing the current face */
    double  *f_row = adv->val + f*adv->n_rows;

    if (fabs(beta_flx) > cs_math_zero_threshold) {

      /* Consistent part */
      f_row[c] -= beta_flx;

      /* Stabilization part */
      f_row[f] += beta_minus;
      f_row[c] -= beta_minus;
      c_row[f] -= beta_minus;
      c_row[c] += beta_minus;

    }

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when diffusion is activated and an upwind
 *          scheme and a conservative formulation is used
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv_di(const cs_equation_param_t   *eqp,
                              const cs_cell_mesh_t        *cm,
                              cs_real_t                    t_eval,
                              cs_face_mesh_t              *fm,
                              cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);
  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_PEQ |
                      CS_CDO_LOCAL_DFQ));

  const cs_param_advection_scheme_t  adv_scheme = eqp->adv_scheme;

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  cs_real_t  *upwcoef = cb->values + cm->n_ec;
  cs_real_3_t  matnu;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_nvec3_t  dfq = cm->dface[e];

    cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat,
                         dfq.unitv, matnu);

    /* Define a coefficient close to a Peclet number */
    const cs_real_t  diff_contrib = _dp3(dfq.unitv, matnu);
    const double  mean_flux = fluxes[e]/dfq.meas;

    if (diff_contrib > cs_math_zero_threshold)
      upwcoef[e] = cm->edge[e].meas * mean_flux / diff_contrib;
    else
      upwcoef[e] = mean_flux * cs_math_big_r;  /* dominated by convection */

  }  /* Loop on cell edges */

  /* Set the function to compute the weight of upwinding */
  _upwind_weight_t  *get_weight = _assign_weight_func(adv_scheme);

  /* Define the local operator for advection */
  _build_cell_vpfd_upw(cm, get_weight, fluxes, upwcoef, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion and an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwcsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  const cs_param_advection_scheme_t  adv_scheme = eqp->adv_scheme;

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  cs_real_t  *upwcoef = cb->values + cm->n_ec;
  for (short int e = 0; e < cm->n_ec; e++)
    upwcoef[e] = fluxes[e]/cm->dface[e].meas;

  /* Set the function to compute the weight of upwinding */
  _upwind_weight_t  *get_weight = _assign_weight_func(adv_scheme);

  /* Define the local operator for advection */
  _build_cell_vpfd_upw(cm, get_weight, fluxes, upwcoef, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cencsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);

  /* Sanity check */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);  /* Sanity check */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV));

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Define the local operator for advection */
  _build_cell_vpfd_cen(cm, fluxes, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a mixed centered/upwind scheme with
 *          a conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_mcucsv(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);

  /* Sanity check */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);  /* Sanity check */
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV));

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Define the local operator for advection */
  _build_cell_vpfd_mcu(cm, eqp->upwind_portion, fluxes, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when diffusion is activated and an upwind
 *          scheme and a conservative formulation is used
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc_di(const cs_equation_param_t   *eqp,
                              const cs_cell_mesh_t        *cm,
                              cs_real_t                    t_eval,
                              cs_face_mesh_t              *fm,
                              cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ));

  const cs_param_advection_scheme_t  adv_scheme = eqp->adv_scheme;

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  cs_real_t  *upwcoef = cb->values + cm->n_ec;
  cs_real_3_t  matnu;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_nvec3_t  dfq = cm->dface[e];

    cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat, dfq.unitv, matnu);

    const cs_real_t  diff_contrib = _dp3(dfq.unitv, matnu);
    const double  mean_flux = fluxes[e]/dfq.meas;

    if (diff_contrib > cs_math_zero_threshold)
      upwcoef[e] = cm->edge[e].meas * mean_flux / diff_contrib;
    else
      upwcoef[e] = mean_flux * cs_math_big_r;  /* dominated by convection */

  }  /* Loop on cell edges */

  /* Set the function to compute the weight of upwinding */
  _upwind_weight_t  *get_weight = _assign_weight_func(adv_scheme);

  /* Define the local operator for advection */
  _build_cell_epcd_upw(cm, get_weight, fluxes, upwcoef,        adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme without diffusion when an upwind scheme and a
 *          conservative formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_upwnoc(const cs_equation_param_t   *eqp,
                           const cs_cell_mesh_t        *cm,
                           cs_real_t                    t_eval,
                           cs_face_mesh_t              *fm,
                           cs_cell_builder_t           *cb)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);  /* Sanity check */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EV));

  const cs_param_advection_scheme_t  adv_scheme = eqp->adv_scheme;

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Compute the criterion attached to each edge of the cell which is used
     to evaluate how to upwind */
  cs_real_t  *upwcoef = cb->values + cm->n_ec;
  for (short int e = 0; e < cm->n_ec; e++)
    upwcoef[e] = fluxes[e]/cm->dface[e].meas;

  /* Set the function to compute the weight of upwinding */
  _upwind_weight_t  *get_weight = _assign_weight_func(adv_scheme);

  /* Define the local operator for advection */
  _build_cell_epcd_upw(cm, get_weight, fluxes, upwcoef, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex-based scheme when a centered scheme and a non-conservative
 *          formulation is used.
 *          The local matrix related to this operator is stored in cb->loc
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_cennoc(const cs_equation_param_t    *eqp,
                           const cs_cell_mesh_t         *cm,
                           cs_real_t                     t_eval,
                           cs_face_mesh_t               *fm,
                           cs_cell_builder_t            *cb)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV));

  /* Initialize the local matrix structure */
  cs_sdm_t  *adv = cb->loc;
  cs_sdm_square_init(cm->n_vc, adv);

  /* Compute the flux across the dual face attached to each edge of the cell */
  cs_real_t  *fluxes = cb->values;  /* size n_ec */
  cs_advection_field_cw_dface_flux(cm, eqp->adv_field, t_eval, fluxes);

  /* Define the local operator for advection */
  _build_cell_epcd_cen(cm, fluxes, adv);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme when the advection field is cellwise
 *          constant
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_cw_cst(const cs_equation_param_t   *eqp,
                            const cs_cell_mesh_t        *cm,
                            cs_real_t                    t_eval,
                            cs_face_mesh_t              *fm,
                            cs_cell_builder_t           *cb)
{
  CS_UNUSED(t_eval);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB);
  assert(eqp->adv_formulation == CS_PARAM_ADVECTION_FORM_NONCONS);
  assert(eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_CIP);
  assert(cs_advection_field_is_cellwise(eqp->adv_field));
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_EF  | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  const int  n_sysc = cm->n_vc + 1;

  /* Initialize local matrix structure */
  cs_sdm_t  *a = cb->loc;
  cs_sdm_square_init(n_sysc, a);

  /* Use a cellwise constant approximation of the advection field */
  cs_nvec3_t  adv_cell;
  cs_advection_field_get_cell_vector(cm->c_id, eqp->adv_field, &adv_cell);

  if (adv_cell.meas < cs_math_get_machine_epsilon())
    return;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 2
  const int n_vc = cm->n_vc;
  cs_sdm_t  *cons = cs_sdm_square_create(n_sysc);

  cons->n_rows = n_sysc;
  for (short int v = 0; v < n_sysc*n_sysc; v++) cons->val[v] = 0;
#endif

  /* Stabilization coefficent * |beta_c| */
  const double  stab_coef = cs_cip_stab_coef * adv_cell.meas;

  /* Temporary buffers
     bgc  stored in cb->values (size n_fc)
     tef  stored in cb->values + cm->n_fc (size 2*n_ec)

     bgvf stored in cb->vectors (size: 2*n_ec)
  */
  cs_sdm_t  *af = cb->aux;
  double  *tef_save = cb->values + cm->n_fc;  /* size = 2*n_ec */

  for (short int f = 0; f < cm->n_fc; f++) {  /* Loop on cell faces */

    /* Build a facewise view of the mesh */
    cs_face_mesh_build_from_cell_mesh(cm, f, fm);

    const int  n_sysf = fm->n_vf + 1;

    /* Initialize the local face matrix */
    af->n_rows = n_sysf;
    for (short int i = 0; i < n_sysf*n_sysf; i++) af->val[i] = 0.;

    /* Store tef areas for a future usage (Second part of the stabilization) */
    const short int  fshift = cm->f2e_idx[f];
    double  *tef = tef_save + fshift;
    for (short int e = 0; e < fm->n_ef; e++) tef[e] = fm->tef[e];

    /* Update the face matrix inside with the consistent */
    _vcb_cellwise_consistent_part(adv_cell, cm, fm, cb);

    /* Build the first part of the stabilization. (inside pfc) */
    _vcb_stabilization_part1(cm, fm, stab_coef, cb);

    /* Reorder v1 and v2 to insure coherency for edges shared between
       two faces */
    cs_real_3_t  *bgvf = cb->vectors + fshift;
    for (short int e = 0; e < fm->n_ef; e++) {
      if (fm->v_ids[fm->e2v_ids[2*e]] > fm->v_ids[fm->e2v_ids[2*e+1]]) {
        double  save = bgvf[e][0];
        bgvf[e][0] = bgvf[e][1];
        bgvf[e][1] = save;
      }
    }

    /* Add the face matrix to the cell matrix */
    for (short int vi = 0; vi < fm->n_vf; vi++) {

      double  *aci = a->val + n_sysc*fm->v_ids[vi];
      const double *afi = af->val + n_sysf*vi;
      for (short int vj = 0; vj < fm->n_vf; vj++)
        aci[fm->v_ids[vj]] += afi[vj];   /* (i,j) face --> cell */
      aci[cm->n_vc] += afi[fm->n_vf];    /* (i,c) face --> cell */

    }

    double  *acc = a->val + n_sysc*cm->n_vc;
    const double  *afc = af->val + n_sysf*fm->n_vf;
    for (short int vj = 0; vj < fm->n_vf; vj++)
      acc[fm->v_ids[vj]] += afc[vj];      /* (c,j) face --> cell */
    acc[cm->n_vc] += afc[fm->n_vf];

  } /* Loop on cell faces */

  /* Build the second part of the stabilization. */
  _vcb_stabilization_part2(cm, stab_coef, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix (CW version)");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
  }
#if CS_CDO_ADVECTION_DBG > 2
  cons = cs_sdm_free(cons);
#endif
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the convection operator attached to a cell with a CDO
 *          vertex+cell-based scheme
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval    time at which one evaluates the advection field
 * \param[in, out] fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb(const cs_equation_param_t   *eqp,
                     const cs_cell_mesh_t        *cm,
                     cs_real_t                    t_eval,
                     cs_face_mesh_t              *fm,
                     cs_cell_builder_t           *cb)
{
  CS_UNUSED(t_eval);

  /* Sanity checks */
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB);
  assert(eqp->adv_formulation == CS_PARAM_ADVECTION_FORM_NONCONS);
  assert(eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_CIP);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_EF  | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  const int  n_sysc = cm->n_vc + 1;

  /* Initialize local matrix structure */
  cs_sdm_t  *a = cb->loc;
  cs_sdm_square_init(n_sysc, a);

  /* Use a cellwise constant approximation of the advection field */
  cs_nvec3_t  adv_cell;
  cs_advection_field_get_cell_vector(cm->c_id, eqp->adv_field, &adv_cell);

  if (adv_cell.meas < cs_math_get_machine_epsilon())
    return;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 2
  const int n_vc = cm->n_vc;
  cs_sdm_t  *cons = cs_sdm_square_create(n_sysc);

  cons->n_rows = n_sysc;
  for (short int v = 0; v < n_sysc*n_sysc; v++) cons->val[v] = 0;
#endif

  /* Stabilization coefficent * |beta_c| */
  const double  stab_coef = cs_cip_stab_coef * adv_cell.meas;

  /* Temporary buffers:
     bgc  stored in cb->values (size n_fc)
     tef  stored in cb->values + cm->n_fc (size 2*n_ec)

     bgvf stored in cb->vectors (size: 2*n_ec)
  */
  cs_sdm_t  *af = cb->aux;
  double  *tef_save = cb->values + cm->n_fc;  /* size = 2*n_ec */

  for (short int f = 0; f < cm->n_fc; f++) {  /* Loop on cell faces */

    /* Build a facewise view of the mesh */
    cs_face_mesh_build_from_cell_mesh(cm, f, fm);

    const int  n_sysf = fm->n_vf + 1;

    /* Initialize the local face matrix */
    cs_sdm_square_init(n_sysf, af);

    /* Store tef areas for a future usage (Second part of the stabilization) */
    const short int  fshift = cm->f2e_idx[f];
    double  *tef = tef_save + fshift;
    for (short int e = 0; e < fm->n_ef; e++) tef[e] = fm->tef[e];

    /* Initialize and update the face matrix inside (build bgvf inside) */
    _vcb_consistent_part(eqp->adv_field, adv_cell, cm, fm, t_eval, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 2
    _assemble_face(cm, fm, af, cons);
#endif

    /* Build the first part of the stabilization: that inside pfc */
    _vcb_stabilization_part1(cm, fm, stab_coef, cb);

    /* Reorder v1 and v2 to insure coherency for edges shared between
       two faces */
    cs_real_3_t  *bgvf = cb->vectors + fshift;
    for (short int e = 0; e < fm->n_ef; e++) {
      if (fm->v_ids[fm->e2v_ids[2*e]] > fm->v_ids[fm->e2v_ids[2*e+1]]) {
        double  save = bgvf[e][0];
        bgvf[e][0] = bgvf[e][1];
        bgvf[e][1] = save;
      }
    }

    /* Add the face matrix to the cell matrix */
    _assemble_face(cm, fm, af, a);

  } /* Loop on cell faces */

  /* Build the second part of the stabilization. */
  _vcb_stabilization_part2(cm, stab_coef, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, NULL)) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Cell advection matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, cb->loc);
#if CS_CDO_ADVECTION_DBG > 2
  cs_log_printf(CS_LOG_DEFAULT, "\n>> Advection matrix (CONSISTENCY PART)");
  cs_sdm_dump(cm->c_id, NULL, NULL, cons);
  cons = cs_sdm_free(cons);
#endif
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a convection builder structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vb_bc(const cs_cell_mesh_t       *cm,
                       const cs_equation_param_t  *eqp,
                       cs_real_t                   t_eval,
                       cs_face_mesh_t             *fm,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FV));

  cs_real_t  *tmp_rhs = cb->values;
  cs_real_t  *mat_diag = cb->values + cm->n_vc;
  cs_real_t  *v_nflx = cb->values + 2*cm->n_vc;

  const cs_adv_field_t  *adv = eqp->adv_field;

  /* Reset local temporay RHS and diagonal contributions */
  for (short int v = 0; v < cm->n_vc; v++) mat_diag[v] = tmp_rhs[v] = 0;

  /* Add diagonal term for vertices attached to a boundary face where
     the advection field points inward. */
  for (short int i = 0; i < csys->n_bc_faces; i++) {  /* Loop on border faces */

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    cs_advection_field_cw_boundary_f2v_flux(cm, adv, f, t_eval, v_nflx);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, " %s: f:%d --> bndy_flux:", __func__, f);
      for (short int v = 0; v < cm->n_vc; v++)
        cs_log_printf(CS_LOG_DEFAULT, " v%d:%e", v, v_nflx[v]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }
#endif

    if (eqp->adv_formulation == CS_PARAM_ADVECTION_FORM_CONSERV) {

      for (int iv = cm->f2v_idx[f]; iv < cm->f2v_idx[f+1]; iv++) {

        const short int  v_id = cm->f2v_ids[iv];

        if (v_nflx[v_id] < 0) {
          /* advection field is inward w.r.t. the face normal */
          if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET)
            /* Homogoneous Dirichlet don't contribute. Other Bcs are invalid */
            tmp_rhs[v_id] -= v_nflx[v_id] * csys->dir_values[v_id];
        }
        else  /* advection is oriented outward */
          mat_diag[v_id] += v_nflx[v_id];

      } /* Loop on face vertices */

    }
    else { /* Non-conservative formulation */

      for (int iv = cm->f2v_idx[f]; iv < cm->f2v_idx[f+1]; iv++) {

        const short int  v_id = cm->f2v_ids[iv];

        if (v_nflx[v_id] < 0) {

          /* advection field is inward w.r.t. the face normal */
          if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET)
            /* Homogoneous Dirichlet don't contribute. Other Bcs are invalid */
            tmp_rhs[v_id] -= v_nflx[v_id] * csys->dir_values[v_id];

          mat_diag[v_id] -= v_nflx[v_id];

        }

      } /* Loop on face vertices */

    } /* Type of formulation */

  } /* Loop on border faces */

  /* Update the diagonal and the RHS of the local system matrix */
  for (short int v = 0; v < cm->n_vc; v++)  {
    csys->mat->val[v*cm->n_vc + v] += mat_diag[v];
    csys->rhs[v] += tmp_rhs[v];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the BC contribution for the convection operator with CDO
 *          V+C schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      eqp     pointer to a cs_equation_param_t structure
 * \param[in]      t_eval  time at which one evaluates the advection field
 * \param[in, out] fm      pointer to a cs_face_mesh_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 * \param[in, out] csys    cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_vcb_bc(const cs_cell_mesh_t        *cm,
                        const cs_equation_param_t   *eqp,
                        cs_real_t                    t_eval,
                        cs_face_mesh_t              *fm,
                        cs_cell_builder_t           *cb,
                        cs_cell_sys_t               *csys)
{
  if (csys == NULL)
    return;

  /* Sanity checks */
  assert(cm != NULL && eqp != NULL);
  assert(csys->mat->n_rows == cm->n_vc + 1);
  assert(eqp->adv_formulation == CS_PARAM_ADVECTION_FORM_NONCONS);
  assert(eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_CIP);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  const cs_adv_field_t  *adv = eqp->adv_field;

  /* Loop on border faces */
  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];
    const cs_real_t  nflx = cs_advection_field_cw_boundary_face_flux(t_eval, f,
                                                                     cm, adv);
    const cs_real_t  beta_nflx = 0.5 * (fabs(nflx) - nflx);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_ADVECTION_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_log_printf(CS_LOG_DEFAULT, " %s: f:%d --> bndy_flux: %e\n",
                    __func__, f, nflx);
#endif

    if (beta_nflx > 0) {

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      cs_hodge_compute_wbs_surfacic(fm, cb->aux);

      _update_vcb_system_with_bc(beta_nflx/fm->face.meas, fm, cb->aux, csys);

    }  /* At least for one face vertex, beta_nf > 0 */

  }  /* Loop on border faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the value of the upwinding coefficient in each cell knowing
 *          the related Peclet number
 *
 * \param[in]      cdoq      pointer to the cdo quantities structure
 * \param[in]      scheme    type of scheme used for the advection term
 * \param[in, out] coefval   pointer to the pointer of real numbers to fill
 *                           in: Peclet number in each cell
 *                           out: value of the upwind coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_advection_cell_upwind_coef(const cs_cdo_quantities_t    *cdoq,
                                  cs_param_advection_scheme_t   scheme,
                                  cs_real_t                     coefval[])
{
  /* Sanity check */
  assert(coefval != NULL);

  /* Set the function to compute the weight of upwinding */
  _upwind_weight_t  *get_weight = _assign_weight_func(scheme);

  for (cs_lnum_t  c_id = 0; c_id < cdoq->n_cells; c_id++)
    coefval[c_id] = get_weight(coefval[c_id]);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
