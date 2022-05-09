/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO schemes
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
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_property.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_diffusion.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdo_diffusion.c

  \brief  Build discrete stiffness matrices and handled boundary conditions for
          diffusion term in CDO vertex-based and vertex+cell schemes.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDO_DIFFUSION_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */

#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Local variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Pre-compute the product between a property and the face vector areas
 *
 * \param[in]       pty       pointer to a \ref cs_property_data_t structure
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  kappa_f   property against face vector for all faces
 */
/*----------------------------------------------------------------------------*/

static inline void
_compute_kappa_f(const cs_property_data_t  *pty,
                 const cs_cell_mesh_t      *cm,
                 cs_real_3_t               *kappa_f)
{
  if (pty->is_unity) {
    for (short int f = 0; f < cm->n_fc; f++) {
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = cm->face[f].meas*cm->face[f].unitv[k];
    }
  }
  else if (pty->is_iso) {
    for (short int f = 0; f < cm->n_fc; f++) {
      const cs_real_t  coef = cm->face[f].meas*pty->value;
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = coef * cm->face[f].unitv[k];
    }
  }
  else {
    for (short int f = 0; f < cm->n_fc; f++) {
      cs_math_33_3_product(pty->tensor, cm->face[f].unitv, kappa_f[f]);
      for (short int k = 0; k < 3; k++) kappa_f[f][k] *= cm->face[f].meas;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute \f$ \int_{fb} \nabla (u) \cdot \nu_{fb} v \f$ where \p fb
 *         is a boundary face (Co+St algorithm)
 *
 * \param[in]       fb       index of the boundary face on the local numbering
 * \param[in]       cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       hodgep   pointer to a \ref cs_hodge_param_t structure
 * \param[in]       kappa_f  property against face vector for all faces
 * \param[in, out]  ntrgrd   pointer to a local matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdofb_normal_flux_reco(short int                  fb,
                        const cs_cell_mesh_t      *cm,
                        const cs_hodge_param_t    *hodgep,
                        const cs_real_3_t         *kappa_f,
                        cs_sdm_t                  *ntrgrd)
{
  assert(hodgep->type == CS_HODGE_TYPE_EDFP);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFC));
  assert(cm->f_sgn[fb] == 1);  /* +1 because it's a boundary face */

  const short int  n_fc = cm->n_fc;
  const cs_quant_t  pfbq = cm->face[fb];
  const cs_nvec3_t  debq = cm->dedge[fb];
  const double  inv_volc = 1./cm->vol_c;

  /* beta * |fb|^2 * nu_{fb}^T.kappa.nu_{fb} / |pvol_fb| */

  const double  stab_scaling =
    hodgep->coef * pfbq.meas * _dp3(kappa_f[fb], pfbq.unitv) / cm->pvol_f[fb];

  /* Get the 'fb' row */

  cs_real_t  *ntrgrd_fb = ntrgrd->val + fb * (n_fc + 1);
  double  row_sum = 0.0;

  for (short int f = 0; f < n_fc; f++) {

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

  ntrgrd_fb[n_fc] += row_sum;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the cellwise consistent gradient for Vb schemes
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in, out] grd_cell  cellwise gradient vector attached to each vertex
 *                           degree of freedom
 */
/*----------------------------------------------------------------------------*/

static void
_svb_cellwise_grd(const cs_cell_mesh_t   *cm,
                  cs_real_3_t             grd_cell[])
{
  /* Reset array */

  for (short int v = 0; v < cm->n_vc; v++)
    grd_cell[v][0] = grd_cell[v][1] = grd_cell[v][2] = 0;

  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *_v = cm->e2v_ids + 2*e; /* v1 = _v[0], v2 = _v[1] */
    const short int  sgn_v1 = cm->e2v_sgn[e]; /* sgn_v2 = - sgn_v1 */
    const cs_nvec3_t  dfq = cm->dface[e];

    grd_cell[_v[0]][0] += sgn_v1 * dfq.meas * dfq.unitv[0];
    grd_cell[_v[0]][1] += sgn_v1 * dfq.meas * dfq.unitv[1];
    grd_cell[_v[0]][2] += sgn_v1 * dfq.meas * dfq.unitv[2];
    grd_cell[_v[1]][0] -= sgn_v1 * dfq.meas * dfq.unitv[0];
    grd_cell[_v[1]][1] -= sgn_v1 * dfq.meas * dfq.unitv[1];
    grd_cell[_v[1]][2] -= sgn_v1 * dfq.meas * dfq.unitv[2];

  } /* Loop on cell edges */

  const cs_real_t  over_vol_c = 1/cm->vol_c;
  for (short int v = 0; v < cm->n_vc; v++) {
    grd_cell[v][0] *= over_vol_c;
    grd_cell[v][1] *= over_vol_c;
    grd_cell[v][2] *= over_vol_c;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal diffusive flux reconstruction for Vb schemes
 *          using the COST/VORONOI/BUBBLE algorithm. This normal flux is
 *          computed in each t_{ek,f}: the triangle with base e_k and apex xf
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      dbeta     d times the value of the stabilization coef.
 * \param[in]      ek        edge id to handle (related to p_{ek,c} and t_{ek,f}
 * \param[in]      pekq      edge quantities for ek
 * \param[in]      dfkq      dual face quantities for ek
 * \param[in]      grd_cell  cellwise gradient vector attached to each vertex
 *                           degree of freedom
 * \param[in]      mnu       diff. property times the face normal of t_{ek,f}
 * \param[in, out] nflux     reconstructed normal flux related to each vertex
 */
/*----------------------------------------------------------------------------*/

static void
_nflux_reco_svb_ocs_in_pec(const cs_cell_mesh_t   *cm,
                           const double            dbeta,
                           const short int         ek,
                           const cs_quant_t        pekq,
                           const cs_nvec3_t        dfkq,
                           const cs_real_3_t       grd_cell[],
                           const cs_real_3_t       mnu,
                           cs_real_t               nflux[])
{
  /* Reset the normal flux for each vertex of the cell */

  for (short int v = 0; v < cm->n_vc; v++) nflux[v] = 0;

  /* Compute the gradient reconstruction for the potential arrays of the
   * canonical basis (size = number of cell vertices)
   * L_{E_c}(GRAD(p)) on t_{e_k,f} for each vertex of the cell
   * The restriction to t_{e_k,f} is identical to the restriction to p_{e_k,f}
   *
   * L_{E_c}(GRAD(p))|t_{e_k,f}
   *   =  Gc(p)                                              ──> Consist. part
   *    + beta*fd_ek(c)/(p_{ek,c} * (GRAD(p)|ek - ek.Gc(p) ) ──> Stabili. part
   */

  const cs_real_t  stab_ek = dbeta/(_dp3(pekq.unitv, dfkq.unitv));
  const short int  *_vk = cm->e2v_ids + 2*ek; /* vk0 = _v[0], vk1 = _v[1] */
  const short int  sgn_vk0 = cm->e2v_sgn[ek]; /* sgn_vk1 = - sgn_vk0 */

  /* Consistent + Stabilization part */

  for (short int v = 0; v < cm->n_vc; v++) {

    const cs_real_t  ekgc = _dp3(pekq.unitv, grd_cell[v]);

    /* Initialize the reconstruction of the gradient with the consistent part */

    cs_real_3_t  le_grd = {grd_cell[v][0], grd_cell[v][1], grd_cell[v][2]};

    if (v == _vk[0]) { /* v = vk0 belonging to edge ek */

      const cs_real_t  stab_coef = stab_ek * (sgn_vk0/pekq.meas - ekgc);
      for (int p = 0; p < 3; p++)
        le_grd[p] += stab_coef * dfkq.unitv[p];  /* Stabilization part */

    }
    else if (v == _vk[1]) { /* v = vk1 belonging to edge ek */

      const cs_real_t  stab_coef = -stab_ek * (sgn_vk0/pekq.meas + ekgc);
      for (int p = 0; p < 3; p++)
        le_grd[p] += stab_coef * dfkq.unitv[p];  /* Stabilization part */

    }
    else { /* Other cell vertices */

      const cs_real_t  stab_coef = -stab_ek * ekgc;
      for (int p = 0; p < 3; p++)
        le_grd[p] += stab_coef * dfkq.unitv[p];  /* Stabilization part */

    }

    nflux[v] = -_dp3(mnu, le_grd); /* Minus because -du/dn */

  } /* Loop on cell vertices */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          COST/BUBLLE or Voronoï algo. is used for reconstructing a flux
 *          from the degrees of freedom
 *
 * \param[in]      fm      pointer to a cs_face_mesh_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu     property tensor times the face normal
 * \param[in]      dbeta   d times the  value of the stabilization coef.
 * \param[in, out] cb      pointer to a cell builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op. i.e.
 *                         the flux operator
 */
/*----------------------------------------------------------------------------*/

static void
_vb_ocs_normal_flux_op(const short int           f,
                       const cs_cell_mesh_t     *cm,
                       const cs_real_3_t         mnu,
                       double                    dbeta,
                       cs_cell_builder_t        *cb,
                       cs_sdm_t                 *ntrgrd)
{
  cs_real_t  *nflux = cb->values;       /* size = cm->n_vc */
  cs_real_3_t  *grd_cell = cb->vectors; /* size = cm->n_vc */

  /* Compute the cellwise-constant gradient for each array of the canonical
   * basis of the potential DoFs (size = cm->n_vc)
   */

  _svb_cellwise_grd(cm, grd_cell);

  /* Loop on border face edges */

  for (int fe_idx = cm->f2e_idx[f]; fe_idx < cm->f2e_idx[f+1]; fe_idx++) {

    const short int  ek = cm->f2e_ids[fe_idx];
    const cs_real_t  tekf = cm->tef[fe_idx];
    const cs_quant_t  pekq = cm->edge[ek];
    const cs_nvec3_t  dfkq = cm->dface[ek];
    const short int  *_vk = cm->e2v_ids + 2*ek; /* v1 = _v[0], v2 = _v[1] */

    /* Compute the reconstructed normal flux */

    _nflux_reco_svb_ocs_in_pec(cm, dbeta, ek, pekq, dfkq,
           (const cs_real_3_t *)grd_cell,
                                mnu,
                                nflux);

    for (short int v = 0; v < cm->n_vc; v++) {

      const double  contrib_v = 0.5 * nflux[v] * tekf;

      ntrgrd->val[_vk[0]*cm->n_vc + v] += contrib_v;
      ntrgrd->val[_vk[1]*cm->n_vc + v] += contrib_v;

    } /* Loop on cell vertices */

  } /* Loop on face edges */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 2
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id, cm->f_ids[f]);
  cs_sdm_dump(cm->c_id, NULL, NULL, ntrgrd);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the transposed of the normal trace operator and the full
 *          boundary operator to enforce weakly a Robin BC for a given border
 *          face when a COST algo. is used for reconstructing the degrees of
 *          freedom
 *
 * \param[in]      fm         pointer to a cs_face_mesh_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu        property tensor times the face normal
 * \param[in]      beta       value of the stabilization coef. related to reco.
 * \param[in, out] cb         pointer to a cell builder structure
 * \param[in, out] bc_op      local matrix storing du/dn.dv/dn
 * \param[in, out] ntrgrd_tr  local matrix related to the transposed version of
 *                            the normal trace op. (i.e. the flux operator)
 */
/*----------------------------------------------------------------------------*/

static void
_vb_cost_full_flux_op(const short int           f,
                      const cs_cell_mesh_t     *cm,
                      const cs_real_3_t         mnu,
                      double                    beta,
                      cs_cell_builder_t        *cb,
                      cs_sdm_t                 *bc_op,
                      cs_sdm_t                 *ntrgrd_tr)
{
  cs_real_t  *nflux = cb->values;         /* size = cm->n_vc */
  cs_real_3_t  *grd_cell = cb->vectors;   /* size = cm->n_vc */

  /* Initialize the local operators */

  cs_sdm_square_init(cm->n_vc, bc_op);
  cs_sdm_square_init(cm->n_vc, ntrgrd_tr);

  /* Compute the cellwise-constant gradient for each array of the canonical
   * basis of the potential DoFs (size = cm->n_vc) */

  _svb_cellwise_grd(cm, grd_cell);

  /* Loop on border face edges */

  for (int fe_idx = cm->f2e_idx[f]; fe_idx < cm->f2e_idx[f+1]; fe_idx++) {

    const short int  ek = cm->f2e_ids[fe_idx];
    const cs_real_t  tekf = cm->tef[fe_idx];
    const cs_quant_t  pekq = cm->edge[ek];
    const cs_nvec3_t  dfkq = cm->dface[ek];
    const short int  *_vk = cm->e2v_ids + 2*ek; /* v1 = _v[0], v2 = _v[1] */

    /* Compute the reconstructed normal flux */

    _nflux_reco_svb_ocs_in_pec(cm, beta, ek, pekq, dfkq,
          (const cs_real_3_t *)grd_cell,
                               mnu,
                               nflux);

    for (short int vi = 0; vi < cm->n_vc; vi++) {

      const double  contrib_vi = nflux[vi] * tekf;

      ntrgrd_tr->val[vi*cm->n_vc + _vk[0]] += contrib_vi;
      ntrgrd_tr->val[vi*cm->n_vc + _vk[1]] += contrib_vi;

      for (short int vj = 0; vj < cm->n_vc; vj++)
        bc_op->val[vi*cm->n_vc + vj] += contrib_vi * nflux[vj];

    } /* Loop on cell vertices (vi) */

  } /* Loop on face edges */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 2
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Full flux.Op (NGRD.NGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id, cm->f_ids[f]);
  cs_log_printf(CS_LOG_DEFAULT, " grd_cell:\n>> ");
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < cm->n_vc; i++)
      cs_log_printf(CS_LOG_DEFAULT, "% .4e ", grd_cell[i][k]);
    if (k < 2)
      cs_log_printf(CS_LOG_DEFAULT, "\n>> ");
    else
      cs_log_printf(CS_LOG_DEFAULT, "\n");
  }
  cs_sdm_dump(cm->c_id, NULL, NULL, bc_op);
  cs_sdm_dump(cm->c_id, NULL, NULL, ntrgrd_tr);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

static void
_vb_wbs_normal_flux_op(const cs_face_mesh_t     *fm,
                       const cs_cell_mesh_t     *cm,
                       const cs_real_3_t         pty_nuf,
                       cs_cell_builder_t        *cb,
                       cs_sdm_t                 *ntrgrd)
{
  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in cb->values and cb->vectors */

  cs_real_t  *l_vc = cb->values;
  cs_real_3_t  *mng_ef = cb->vectors;
  cs_real_3_t  *u_vc = cb->vectors + fm->n_vf;

  /* Initialize the local operator */

  cs_sdm_square_init(cm->n_vc, ntrgrd);

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */

  cs_compute_grdfc_fw(fm, grd_c);

  /* f_coef = (pty_tensor*nu_f).grd_c */

  const cs_real_t  f_coef = fm->face.meas * _dp3(pty_nuf, grd_c);

  /* Compute xc --> xv length and unit vector for all face vertices */

  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute a weight for each vertex of the current face */

  double  sum_ef = 0.;
  for (short int e = 0; e < fm->n_ef; e++) {

    /* Gradient of the Lagrange function related to v1 and v2 */

    cs_compute_grd_ve(fm->e2v_ids[2*e],
                      fm->e2v_ids[2*e+1],
                      fm->dedge,
                      (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the unity */

    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = fm->tef[e] * cs_math_1ov3;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

    sum_ef += mng_ef[e][2];

  } /* End of loop on face edges */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    const short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_vi = ntrgrd->val + vi*cm->n_vc;

    /* Default contribution for this line */

    const double  default_coef = f_coef * fm->wvf[vfi];
    for (short int vj = 0; vj < cm->n_vc; vj++)
      ntrgrd_vi[vj] = default_coef * cm->wvc[vj];  /* two contributions */

    /* Block Vf x Vf */

    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

      short int vj = fm->v_ids[vfj];
      ntrgrd_vi[vj] += sum_ef * fm->wvf[vfi] * fm->wvf[vfj];

      double  entry_ij = 0.;
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->e2v_ids[2*e];
        const short int  v2 = fm->e2v_ids[2*e+1];

        if (vfj == v1)
          entry_ij += mng_ef[e][0] * fm->wvf[vfi];
        else if (vfj == v2)
          entry_ij += mng_ef[e][1] * fm->wvf[vfi];

        if (vfi == v1 || vfi == v2) { /* i in e */
          entry_ij += fm->wvf[vfj] * mng_ef[e][2];
          if (vfj == v1) /* j is also in e */
            entry_ij += mng_ef[e][0];
          if (vfj == v2)
            entry_ij += mng_ef[e][1];
        }

      }  /* Loop on face edges */

      ntrgrd_vi[vj] += entry_ij;

    }  /* Loop on face vertices (vj) */

  }  /* Loop on face vertices (vi) */

  /* The flux is -1 times the reconstruction of the gradient */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    const short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_vi = ntrgrd->val + vi*cm->n_vc;
    for (short int vj = 0; vj < cm->n_vc; vj++)
      ntrgrd_vi[vj] *= -1;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 2
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id,cm->f_ids[fm->f_id]);
  cs_sdm_dump(cm->c_id,  NULL, NULL, ntrgrd);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *          Specific to CDO-V+C schemes
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

static void
_vcb_wbs_normal_flux_op(const cs_face_mesh_t     *fm,
                        const cs_cell_mesh_t     *cm,
                        const cs_real_3_t         pty_nuf,
                        cs_cell_builder_t        *cb,
                        cs_sdm_t                 *ntrgrd)
{
  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in cb->values and cb->tmp-vect */

  cs_real_t  *l_vc = cb->values;
  cs_real_3_t  *mng_ef = cb->vectors;
  cs_real_3_t  *u_vc = cb->vectors + fm->n_vf;

  /* Initialize the local operator */

  cs_sdm_square_init(cm->n_vc + 1, ntrgrd);

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */

  cs_compute_grdfc_fw(fm, grd_c);

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); /* (pty_tensor*nu_f).grd_c */

  /* Compute xc --> xv length and unit vector for all face vertices */

  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute a weight for each vertex of the current face */

  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* Gradient of the Lagrange function related to v1 and v2 */

    cs_compute_grd_ve(v1, v2, fm->dedge, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the unity */

    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = fm->tef[e] * cs_math_1ov3;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

  } /* End of loop on face edges */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_vi = ntrgrd->val + vi*ntrgrd->n_rows;

    /* Contribution to the cell column */

    ntrgrd_vi[cm->n_vc] = fm->wvf[vfi] * fm->face.meas * mng_cf;

    /* Block Vf x Vf */

    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

      double  entry_ij = 0.;
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->e2v_ids[2*e];
        const short int  v2 = fm->e2v_ids[2*e+1];

        double  coef_i = fm->wvf[vfi];
        if (vfi == v1 || vfi == v2)
          coef_i += 1;

        double  coef_j = fm->wvf[vfj] * mng_ef[e][2];
        if (vfj == v1)
          coef_j += mng_ef[e][0];
        else if (vfj == v2)
          coef_j += mng_ef[e][1];

        entry_ij += coef_i * coef_j;

      }  /* Loop on face edges */

      ntrgrd_vi[fm->v_ids[vfj]] += entry_ij;

    }  /* Loop on face vertices (vj) */

  }  /* Loop on face vertices (vi) */

  /* The flux is -1 times the reconstruction of the gradient */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    const short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_vi = ntrgrd->val + vi*ntrgrd->n_rows;
    for (short int vj = 0; vj < ntrgrd->n_rows; vj++)
      ntrgrd_vi[vj] *= -1;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 2
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id,cm->f_ids[fm->f_id]);
  cs_sdm_dump(cm->c_id,  NULL, NULL, ntrgrd);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Weak enforcement of Dirichlet BCs using a Nitsche technique.
 *          Scalar-valued case for CDO-VB or CDO-VCb schemes
 *
 * \param[in]       pcoef     value of the penalization coefficient
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  ntrgrd    pointer to a local matrix structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

static inline void
_svb_nitsche(const double              pcoef,
             const cs_face_mesh_t     *fm,
             cs_sdm_t                 *ntrgrd,
             cs_cell_sys_t            *csys)
{
  assert(pcoef >= 0);
  assert(csys->mat->n_rows == ntrgrd->n_rows);

  /* Update local RHS and system */

  for (short int v = 0; v < fm->n_vf; v++) {

    /* Set the penalty diagonal coefficient */

    const double  pcoef_v = pcoef * fm->wvf[v];
    const short int  vi = fm->v_ids[v];

    ntrgrd->val[vi*(1 + ntrgrd->n_rows)] += pcoef_v;
    csys->rhs[vi] += pcoef_v * csys->dir_values[vi];

  }  /* Dirichlet or homogeneous Dirichlet */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_face_mesh_t                  *fm,
                                cs_hodge_t                      *hodge,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys)
{
  CS_UNUSED(hodge); /* Prototype common to cs_cdo_enforce_bc_t */
  CS_UNUSED(fm);    /*  Hence the unused parameters */
  CS_UNUSED(cm);
  CS_UNUSED(cb);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  /* Penalize diagonal entry (and its rhs if needed) */

  for (short int i = 0; i < csys->n_dofs; i++) {

    if (csys->dof_flag[i] & CS_CDO_BC_HMG_DIRICHLET) {
      csys->mat->val[i + csys->n_dofs*i] += eqp->strong_pena_bc_coeff;
    }
    else if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET) {
      csys->mat->val[i + csys->n_dofs*i] += eqp->strong_pena_bc_coeff;
      csys->rhs[i] += csys->dir_values[i] * eqp->strong_pena_bc_coeff;
    }

  } /* Loop on degrees of freedom */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value.
 *          Case of a cellwise system defined by block.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_face_mesh_t                  *fm,
                                      cs_hodge_t                      *hodge,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys)
{
  CS_UNUSED(hodge);   /* Prototype common to cs_cdo_enforce_bc_t */
  CS_UNUSED(fm);      /*  Hence the unused parameters */
  CS_UNUSED(cm);
  CS_UNUSED(cb);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;
  assert(bd != NULL);

  /* Penalize diagonal entry (and its rhs if needed) */

  int  shift = 0;
  for (short int bi = 0; bi < bd->n_row_blocks; bi++) {

    cs_sdm_t  *mII = cs_sdm_get_block(m, bi, bi);
    assert(mII->n_rows == mII->n_cols);
    cs_real_t  *_rhs = csys->rhs + shift;
    const cs_flag_t  *_flag = csys->dof_flag + shift;
    const cs_real_t  *_dir_val = csys->dir_values + shift;

    for (int i = 0; i < mII->n_rows; i++) {

      if (_flag[i] & CS_CDO_BC_HMG_DIRICHLET) {
        mII->val[i + mII->n_rows*i] += eqp->strong_pena_bc_coeff;
      }
      else if (_flag[i] & CS_CDO_BC_DIRICHLET) {
        mII->val[i + mII->n_rows*i] += eqp->strong_pena_bc_coeff;
        _rhs[i] += _dir_val[i] * eqp->strong_pena_bc_coeff;
      }

    }

    shift += mII->n_rows;

  } /* Block bi */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by keeping the DoFs related to
 *          Dirichlet BCs in the algebraic system (i.e. a weak enforcement)
 *          The corresponding DoFs are algebraically "removed" of the system
 *
 *          |      |     |     |      |     |     |  |     |          |
 *          | Aii  | Aid |     | Aii  |  0  |     |bi|     |bi-Aid.bd |
 *          |------------| --> |------------| and |--| --> |----------|
 *          |      |     |     |      |     |     |  |     |          |
 *          | Adi  | Add |     |  0   |  Id |     |bd|     |    xd    |
 *
 *          where xd is the value of the Dirichlet BC
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_face_mesh_t                  *fm,
                                cs_hodge_t                      *hodge,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys)
{
  CS_UNUSED(eqp);   /* Prototype common to cs_cdo_enforce_bc_t */
  CS_UNUSED(fm);    /* Hence the unused parameters */
  CS_UNUSED(cm);
  CS_UNUSED(hodge);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  double  *x_dir = cb->values;
  double  *ax_dir = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  /* Build x_dir */

  for (short int i = 0; i < csys->n_dofs; i++)
    if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET) /* Only non-homogeneous */
      x_dir[i] = csys->dir_values[i];

  /* Contribution of the Dirichlet conditions */

  cs_sdm_matvec(csys->mat, x_dir, ax_dir);

  /* Second pass: Replace the Dirichlet block by a diagonal block */

  for (short int i = 0; i < csys->n_dofs; i++) {

    if (cs_cdo_bc_is_dirichlet(csys->dof_flag[i])) { /* All Dirichlet:
                                                        homogeneous or not */
      /* Reset row i */

      memset(csys->mat->val + csys->n_dofs*i, 0, csys->n_dofs*sizeof(double));

      /* Reset column i */

      for (short int j = 0; j < csys->n_dofs; j++)
        csys->mat->val[i + csys->n_dofs*j] = 0;
      csys->mat->val[i*(1 + csys->n_dofs)] = 1;

      /* Set the RHS */

      csys->rhs[i] = csys->dir_values[i];

    } /* DoF associated to a Dirichlet BC */
    else
      csys->rhs[i] -= ax_dir[i];  /* Update RHS */

  } /* Loop on degrees of freedom */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by keeping the DoFs related to
 *          Dirichlet BCs in the algebraic system (i.e. a weak enforcement)
 *          The corresponding DoFs are algebraically "removed" of the system
 *          Block version.
 *
 *          |      |     |     |      |     |     |  |     |          |
 *          | Aii  | Aid |     | Aii  |  0  |     |bi|     |bi-Aid.xd |
 *          |------------| --> |------------| and |--| --> |----------|
 *          |      |     |     |      |     |     |  |     |          |
 *          | Adi  | Add |     |  0   |  Id |     |bd|     |    xd    |
 *
 *          where xd is the value of the Dirichlet BC
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_face_mesh_t                  *fm,
                                      cs_hodge_t                      *hodge,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys)
{
  CS_UNUSED(eqp);   /* Prototype common to cs_cdo_enforce_bc_t */
  CS_UNUSED(fm);    /* Hence the unused parameters */
  CS_UNUSED(cm);
  CS_UNUSED(hodge);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  double  *x_dir = cb->values;
  double  *ax_dir = cb->values + csys->n_dofs;
  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;
  assert(bd != NULL);

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  /* Build x_dir */

  for (short int i = 0; i < csys->n_dofs; i++)
    if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET) /* Only non-homogeneous */
      x_dir[i] = csys->dir_values[i];

  /* Contribution of the Dirichlet conditions */

  cs_sdm_block_matvec(csys->mat, x_dir, ax_dir);

  /* Second pass: Replace the Dirichlet block by a diagonal block */

  int  s = 0;
  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    cs_real_t  *_rhs = csys->rhs + s;
    cs_real_t  *_dir = csys->dir_values + s;
    cs_flag_t  *_flg = csys->dof_flag + s;
    cs_sdm_t  *mII = cs_sdm_get_block(m, bi, bi);
    assert(mII->n_rows == mII->n_cols);

    /* Is the current block associated to a Dirichlet BC ? */

    int  n_dir = 0;
    for (int i = 0; i < mII->n_rows; i++)
      if (cs_cdo_bc_is_dirichlet(_flg[i]))
        n_dir += 1;

    if (n_dir > 0) {

      if (n_dir != mII->n_rows)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Partial Dirichlet block not yet implemented", __func__);

      /* Reset block row bi and block colum bi */

      for (int bj = 0; bj < bd->n_col_blocks; bj++) {

        if (bj != bi) {

          cs_sdm_t  *mIJ = cs_sdm_get_block(m, bi, bj);
          cs_sdm_t  *mJI = cs_sdm_get_block(m, bj, bi);

          memset(mIJ->val, 0, mIJ->n_rows*mIJ->n_cols*sizeof(double));
          memset(mJI->val, 0, mJI->n_rows*mJI->n_cols*sizeof(double));

        }
        else { /* mII block --> diagonal */

          memset(mII->val, 0, mII->n_rows*mII->n_rows*sizeof(double));

          for (int i = 0; i < mII->n_rows; i++) {
            mII->val[i + mII->n_rows*i] = 1;
            _rhs[i] = _dir[i]; /* Set the RHS */
          }

        }

      } /* Loop on row/columnn blocks */

    } /* DoF associated to a Dirichlet BC */
    else {

      for (int i = 0; i < mII->n_rows; i++)
        _rhs[i] -= ax_dir[s+i];  /* Update RHS */

    } /* Not a Dirichlet block */

    s += mII->n_rows;

  } /* Block bi */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique.
 *          Case of scalar-valued CDO Face-based schemes
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL);

  const cs_property_data_t  *pdata = hodge->pty_data;
  const double chi =
    eqp->weak_pena_bc_coeff * fabs(pdata->eigen_ratio)*pdata->eigen_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  _compute_kappa_f(pdata, cm, kappa_f);

  /* Initialize the matrix related to this flux reconstruction operator */

  const short int n_dofs = cm->n_fc + 1;
  cs_sdm_t  *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute \int_f du/dn v and update the matrix */

      _cdofb_normal_flux_reco(f, cm, hodge->param,
                              (const cs_real_t (*)[3])kappa_f,
                              bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Second pass: Update the cell system with the bc_op matrix and the Dirichlet
     values. Avoid a truncation error if the arbitrary coefficient of the
     Nitsche algorithm is large
  */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* chi * \meas{f} / h_f  */

      const cs_real_t pcoef = chi * sqrt(cm->face[f].meas);

      /* Diagonal term */

      bc_op->val[f*(n_dofs + 1)] += pcoef;

      /* rhs */

      csys->rhs[f] += pcoef * csys->dir_values[f];

    } /* If Dirichlet */

  } /* Loop on boundary faces */

  /* Update the local system matrix */

  cs_sdm_add(csys->mat, bc_op);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique.
 *          Case of vector-valued CDO Face-based schemes
 *          The idea is to compute the scalar version and dispatch it three
 *          times, one for each Cartesian components
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && csys != NULL && hodge);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  _compute_kappa_f(pty, cm, kappa_f);

  /* Initialize the matrix related this flux reconstruction operator */

  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute \int_f du/dn v and update the matrix */

      _cdofb_normal_flux_reco(f, cm, hodge->param,
                              (const cs_real_t (*)[3])kappa_f,
                              bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Second pass: Update the cell system with the bc_op matrix and the Dirichlet
   * values. Avoid a truncation error if the arbitrary coefficient of the
   * Nitsche algorithm is large
   */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* chi * \meas{f} / h_f  */

      const cs_real_t pcoef = chi * sqrt(cm->face[f].meas);

      /* Diagonal term */

      bc_op->val[f*(n_dofs + 1)] += pcoef;

      /* rhs */

      for (short int k = 0; k < 3; k++)
        csys->rhs[3*f + k] += pcoef * csys->dir_values[3*f + k];

    } /* If Dirichlet */

  } /* Loop on boundary faces */

  assert(pty->is_iso == true); /* if not the case something else TODO ? */

  /* Update the local system matrix */

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
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment
 *          Case of scalar-valued CDO Face-based schemes
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && csys != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  _compute_kappa_f(pty, cm, kappa_f);

  const short int n_dofs = cm->n_fc + 1, n_f = cm->n_fc;
  cs_sdm_t  *bc_op = cb->loc, *bc_op_t = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute \int_f du/dn v and update the matrix */

      _cdofb_normal_flux_reco(f, cm, hodge->param,
                              (const cs_real_t (*)[3])kappa_f,
                              bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Symmetrize the Nitsche contribution */

  cs_real_t *dir_val = cb->values, *u0_trgradv = cb->values + n_dofs;

  /* Putting the face DoFs of the BC, into a face- and cell-DoFs array */

  memcpy(dir_val, csys->dir_values, n_f*sizeof(cs_real_t));
  dir_val[n_f] = 0.;

  /* Update bc_op = bc_op + transp and transp = transpose(bc_op) cb->loc
     plays the role of the flux operator */

  cs_sdm_square_add_transpose(bc_op, bc_op_t);
  cs_sdm_square_matvec(bc_op_t, dir_val, u0_trgradv);

  /* Second pass: Update the cell system with the bc_op matrix and the
   * Dirichlet values. Avoid a truncation error if the arbitrary coefficient of
   * the Nitsche algorithm is large */

  for (short int i = 0; i < n_dofs; i++) /* Cell too! */
    csys->rhs[i] += u0_trgradv[i];

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* chi * \meas{f} / h_f  */

      const cs_real_t pcoef = chi * sqrt(cm->face[f].meas);

      /* Diagonal term */

      bc_op->val[f*(n_dofs + 1)] += pcoef;

      /* rhs */

      csys->rhs[f] += pcoef * csys->dir_values[f];

    } /* If Dirichlet */

  } /* Loop on boundary faces */

  /* Update the local system matrix */

  cs_sdm_add(csys->mat, bc_op);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *          The idea is to compute the scalar version and dispatch it three
 *          times, one for each Cartesian components
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && csys != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  _compute_kappa_f(pty, cm, kappa_f);

  /* Initialize the matrix related this flux reconstruction operator */

  const short int  n_dofs = cm->n_fc + 1; /* n_blocks or n_scalar_dofs */
  cs_sdm_t *bc_op = cb->loc, *bc_op_t = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute \int_f du/dn v and update the matrix */

      _cdofb_normal_flux_reco(f, cm, hodge->param,
                              (const cs_real_t (*)[3])kappa_f,
                              bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Second pass: add the bc_op matrix, add the BC */

  /* Update bc_op = bc_op + transp and bc_op_t = transpose(bc_op) */

  cs_sdm_square_add_transpose(bc_op, bc_op_t);

  /* Fill a face- and cell-DoFs array with the face DoFs values associated to
   * Dirichlet BC. These Dirichlet BC values are interlaced. One first
   * de-interlaces this array to perform the local matrix-vector product:
   * bc_op_t dir_val = u0_trgradv */

  cs_real_t *dir_val = cb->values, *u0_trgradv = cb->values + n_dofs;

  dir_val[cm->n_fc] = 0.;     /* cell DoF is not attached to a Dirichlet */

  for (int k = 0; k < 3; k++) {

    /* Build dir_val */

    for (short int f = 0; f < cm->n_fc; f++)
      dir_val[f] = csys->dir_values[3*f+k];

    cs_sdm_square_matvec(bc_op_t, dir_val, u0_trgradv);

    for (short int i = 0; i < n_dofs; i++) /* Cell too! */
      csys->rhs[3*i+k] += u0_trgradv[i];

  } /* Loop on components */

  /* Third pass: Update the cell system with the bc_op matrix and the Dirichlet
   * values. Avoid a truncation error if the arbitrary coefficient of the
   * Nitsche algorithm is large */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* chi * \meas{f} / h_f  */

      const cs_real_t pcoef = chi * sqrt(cm->face[f].meas);

      /* Diagonal term */

      bc_op->val[f*(n_dofs + 1)] += pcoef;

      /* rhs */

      for (short int k = 0; k < 3; k++)
        csys->rhs[3*f + k] += pcoef * csys->dir_values[3*f + k];

    } /* If Dirichlet */

  } /* Loop on boundary faces */

  assert(pty->is_iso == true); /* if not the case something else TODO ? */

  /* Update the local system matrix */

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
 * \brief   Take into account sliding BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment.
 *          Case of vector-valued CDO Face-based schemes
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vfb_wsym_sliding(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_face_mesh_t                 *fm,
                                  cs_hodge_t                     *hodge,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);
  assert(csys != NULL);

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_sliding == false)
    return;  /* Nothing to do */

  const cs_property_data_t  *pty = hodge->pty_data;

  assert(cm != NULL && cb != NULL);
  assert(pty->is_iso == true); /* if not the case something else TODO ? */

  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;
  const short int  n_f = cm->n_fc;
  const short int  n_dofs = n_f + 1; /* n_blocks or n_scalar_dofs */

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */

  cs_real_3_t  *kappa_f = cb->vectors;
  _compute_kappa_f(pty, cm, kappa_f);

  /* Initialize the matrix related this flux reconstruction operator */

  cs_sdm_t *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_SLIDING) {

      /* Compute \int_f du/dn v and update the matrix */

      _cdofb_normal_flux_reco(f, cm, hodge->param,
                              (const cs_real_t (*)[3])kappa_f,
                              bc_op);

    } /* If sliding */

  } /* Loop boundary faces */

  /* Second pass: add the bc_op matrix, add the BC */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  fi = csys->_f_ids[i];

    if (csys->bf_flag[fi] & CS_CDO_BC_SLIDING) {

      const cs_quant_t  pfq = cm->face[fi];
      const cs_real_t  *ni = pfq.unitv;
      const cs_real_t  ni_ni[9] =
        { ni[0]*ni[0], ni[0]*ni[1], ni[0]*ni[2],
          ni[1]*ni[0], ni[1]*ni[1], ni[1]*ni[2],
          ni[2]*ni[0], ni[2]*ni[1], ni[2]*ni[2]};

      /* chi * \meas{f} / h_f  */

      const cs_real_t  pcoef = chi * sqrt(pfq.meas);

      for (short int xj = 0; xj < n_dofs; xj++) {

        /* It should be done both for face- and cell-defined DoFs */

        if ( xj == fi ) {

          /* Retrieve the 3x3 matrix */

          cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, fi, fi);
          assert(bii->n_rows == bii->n_cols && bii->n_rows == 3);

          const cs_real_t  _val = pcoef + 2 * bc_op->val[n_dofs*fi + fi];

          for (short int k = 0; k < 9; k++)
            bii->val[k] += ni_ni[k] * _val;

        }
        else { /* xj != fi */

          /* Retrieve the 3x3 matrix */

          cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, fi, xj);
          assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);
          cs_sdm_t  *bji = cs_sdm_get_block(csys->mat, xj, fi);
          assert(bji->n_rows == bji->n_cols && bji->n_rows == 3);

          const cs_real_t  _val = bc_op->val[n_dofs*fi + xj];

          for (short int k = 0; k < 9; k++) {
            bij->val[k] += ni_ni[k] * _val;
            /* ni_ni is symmetric */
            bji->val[k] += ni_ni[k] * _val;
          }

        } /* End if */

      } /* Loop on xj */

    } /* If sliding */

  } /* Loop boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Robin BCs.
 *          Case of scalar-valued CDO-Fb schemes with a CO+ST algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_cost_robin(const cs_equation_param_t      *eqp,
                                const cs_cell_mesh_t           *cm,
                                cs_face_mesh_t                 *fm,
                                cs_hodge_t                     *hodge,
                                cs_cell_builder_t              *cb,
                                cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(fm);
  CS_UNUSED(hodge);
  CS_UNUSED(cb);

  assert(cm != NULL && csys != NULL);

  /* Enforcement of the Robin BCs */

  if (csys->has_robin == false)
    return;  /* Nothing to do */

  /* Robin BC expression: K du/dn + alpha*(u - u0) = g */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_ROBIN) {

      /* Robin BC expression: K du/dn + alpha*(u - u0) = g */
      /* ------------------------------------------------- */

      const double  alpha = csys->rob_values[3*f];
      const double  u0 = csys->rob_values[3*f+1];
      const double  g = csys->rob_values[3*f+2];
      const cs_real_t  f_meas = cm->face[f].meas;

      /* Update the RHS and the local system */

      csys->rhs[f] += (alpha*u0 + g)*f_meas;

      cs_real_t  *row_f_val = csys->mat->val + f*csys->n_dofs;
      row_f_val[f] += alpha*f_meas;

    }  /* Robin face */

  } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, csys))
    cs_cell_sys_dump(">> Cell %d, scalar Fb: After Robin", csys);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Robin BCs.
 *          Case of scalar-valued CDO-Vb schemes with a CO+ST algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_cost_robin(const cs_equation_param_t      *eqp,
                                const cs_cell_mesh_t           *cm,
                                cs_face_mesh_t                 *fm,
                                cs_hodge_t                     *hodge,
                                cs_cell_builder_t              *cb,
                                cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(hodge);
  assert(cm != NULL && cb != NULL && csys != NULL);

  /* Enforcement of the Robin BCs */

  if (csys->has_robin == false)
    return;  /* Nothing to do */

  /* Robin BC expression: K du/dn + alpha*(u - u0) = beta */

  cs_sdm_t  *bc_op = cb->loc;

  /* Reset local operator */

  cs_sdm_square_init(cm->n_vc, bc_op);

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_ROBIN) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Robin BC expression: K du/dn + alpha*(u - u0) = beta */
      /* ---------------------------------------------------- */

      const double  alpha = csys->rob_values[3*f];
      const double  u0 = csys->rob_values[3*f+1];
      const double  beta = csys->rob_values[3*f+2];
      const double  g = beta + alpha*u0;

      /* Update the RHS and the local system */

      for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

        const double  f_cap_dualv = fm->face.meas * fm->wvf[vfi];
        const short int  vi = fm->v_ids[vfi];

        csys->rhs[vi] += g*f_cap_dualv;
        bc_op->val[vi*(1+cm->n_vc)] += alpha*f_cap_dualv;

      }

    }  /* Robin face */
  } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, csys)) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Cell %d Vb COST Robin bc matrix",
                  cm->c_id);
    cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, bc_op);
  }
#endif

  /* Add contribution to the linear system */

  cs_sdm_add(csys->mat, bc_op);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account generic BCs by a weak enforcement using Nitsche
 *          technique. According to the settings one can apply Neumann BCs if
 *          alpha = 0, Dirichlet BCs if alpha >> 1 or Robin BCs
 *          Case of scalar-valued CDO-Vb schemes with a CO+ST algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_cost_generic(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_face_mesh_t                 *fm,
                                  cs_hodge_t                     *hodge,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  assert(csys != NULL);

  /* Enforcement of the Robin BCs */

  if (csys->has_robin == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const double  cs_gamma_generic = 1e-2;
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->need_tensor);

  /* Robin BC expression: K du/dn + alpha*(u - u0) = g
   * Store x = u0 + g/alpha */

  cs_real_t  *x = cb->values;
  cs_real_t  *ax = cb->values + cm->n_vc;
  cs_sdm_t  *bc_op = cb->loc;
  cs_sdm_t  *ntrgrd_tr = cb->aux;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_ROBIN) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_cost_full_flux_op(f, cm, pty_nuf, hodgep->coef, cb, bc_op, ntrgrd_tr);

      /* Robin BC expression: K du/dn + alpha*(u - u0) = g */
      /* ------------------------------------------------- */

      const double  alpha = csys->rob_values[3*f];
      const double  u0 = csys->rob_values[3*f+1];
      const double  g = csys->rob_values[3*f+2];
      assert(fabs(alpha) > 0);
      const double  epsilon = 1./alpha;

      const double  hf = sqrt(cm->face[f].meas);
      const double  pena_coef = 1./(epsilon + cs_gamma_generic*hf);
      const double  flux_coef = cs_gamma_generic*hf*pena_coef;

      memset(x, 0, sizeof(cs_real_t)*cm->n_vc);
      for (short int v = 0; v < fm->n_vf; v++) {
        const short int  vi = fm->v_ids[v];
        x[vi] = u0 + epsilon*g;
      }

      /* Update the RHS */

      cs_sdm_square_matvec(ntrgrd_tr, x, ax);
      for (short int v = 0; v < fm->n_vf; v++) {
        const double  pcoef_v = pena_coef * fm->face.meas * fm->wvf[v];
        const short int  vi = fm->v_ids[v];
        csys->rhs[vi] += pcoef_v*x[vi] - flux_coef*ax[vi];
      }

      /* Update the local system matrix:
       *
       * Initially bc_op = du/dn.dv/dn
       * 1) bc_op *= -epsilon*gamma*hf/(epsilon + gamma*hf)
       * 2) bc_op += (u,v) * 1/(epsilon + gamma*hf)
       * 3) bc_op += (-gamma*hf)/(epsilon + gamma*hf) [ntrgrd + ntrgrd_tr]
       */

      /* Step 1 */

      for (int ij = 0; ij < cm->n_vc*cm->n_vc; ij++)
        bc_op->val[ij] *= -epsilon*flux_coef;

      /* Step 2 */

      for (short int v = 0; v < fm->n_vf; v++) {
        const double  pcoef_v = pena_coef * fm->face.meas * fm->wvf[v];
        const short int  vi = fm->v_ids[v];
        bc_op->val[vi*(1 + bc_op->n_rows)] += pcoef_v;
      }

      /* Step 3 */

      cs_sdm_square_2symm(ntrgrd_tr); /* ntrgrd_tr is now a symmetric matrix */
      cs_sdm_add_mult(bc_op, -flux_coef, ntrgrd_tr);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell Vb COST generic bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, bc_op);
      }
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, bc_op);

    }  /* Robin face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of scalar-valued CDO-Vb schemes with an orthogonal
 *          splitting between the consistency/stabilization parts (OCS)
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_ocs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_hodge_t                     *hodge,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const cs_hodge_param_t  *hodgep = hodge->param;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;
  const cs_real_t  dbeta =
    (hodgep->algo == CS_HODGE_ALGO_BUBBLE) ? hodgep->coef : 3*hodgep->coef;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  /* Initialize the local operator */

  cs_sdm_square_init(cm->n_vc, ntrgrd);

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_ocs_normal_flux_op(f, cm, pty_nuf, dbeta, cb, ntrgrd);

      /* Update the RHS and the local system matrix */

      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);

    }  /* Dirichlet face */
  } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, csys)) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Cell Vb.COST Weak bc matrix");
    cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
  }
#endif

  /* Add contribution to the linear system */

  cs_sdm_add(csys->mat, ntrgrd);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. A Dirichlet is set for the three components of the
 *          vector. Case of vector-valued CDO-Vb schemes with an orthogonal
 *          splitting between the consistency/stabilization parts (OCS)
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vvb_ocs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_hodge_t                     *hodge,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_real_t  dbeta =
    (hodgep->algo == CS_HODGE_ALGO_BUBBLE) ? hodgep->coef : 3*hodgep->coef;

  assert(pty->need_tensor);
  cs_sdm_t  *ntrgrd = cb->loc;

  /* Initialize the local operator */

  cs_sdm_square_init(cm->n_vc, ntrgrd);

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_ocs_normal_flux_op(f, cm, pty_nuf, dbeta, cb, ntrgrd);

      /* Update the RHS and the local system matrix */

      const double  pcoef = chi/sqrt(cm->face[f].meas);

      for (short int v = 0; v < fm->n_vf; v++) {

        /* Set the penalty diagonal coefficient */

        const short int  vi = fm->v_ids[v];
        const double  pcoef_v = pcoef * fm->wvf[v];

        ntrgrd->val[vi*(1 + ntrgrd->n_rows)] += pcoef_v;

        /* Update the RHS */

        for (int k = 0; k < 3; k++)
          csys->rhs[3*vi+k] += pcoef_v * csys->dir_values[3*vi+k];

      }  /* Dirichlet or homogeneous Dirichlet */

    }  /* Dirichlet face */
  } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
  if (cs_dbg_cw_test(eqp, cm, csys)) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Cell Vb.COST Weak bc matrix");
    cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
  }
#endif

  /* Add contribution to the linear system */

  assert(pty->is_iso == true); /* if not the case something else TODO ? */

  for (int bvi = 0; bvi < cm->n_vc; bvi++) {
    for (int bvj = 0; bvj < cm->n_vc; bvj++) {

      /* Retrieve the 3x3 matrix */

      cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bvi, bvj);
      assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

      const cs_real_t  _val = ntrgrd->val[cm->n_vc*bvi + bvj];

      /* Update diagonal terms only */

      bij->val[0] += _val;
      bij->val[4] += _val;
      bij->val[8] += _val;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account a sliding BCs.
 *          Case of vector-valued CDO-Vb schemes with a OCS algorithm.
 *          Orthogonal splitting between Consistency/Stabilization parts.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vvb_ocs_sliding(const cs_equation_param_t      *eqp,
                                 const cs_cell_mesh_t           *cm,
                                 cs_face_mesh_t                 *fm,
                                 cs_hodge_t                     *hodge,
                                 cs_cell_builder_t              *cb,
                                 cs_cell_sys_t                  *csys)
{
  /* Enforcement of the sliding BCs */

  if (csys->has_sliding == false)
    return;  /* Nothing to do */

  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_real_t  dbeta =
    (hodgep->algo == CS_HODGE_ALGO_BUBBLE) ? hodgep->coef : 3*hodgep->coef;
  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->is_iso == true); /* if not the case something else TODO ? */

  cs_sdm_t  *ntrgrd = cb->loc;

  /* Initialize the local operator (same as the scalar-valued one) */

  cs_sdm_square_init(cm->n_vc, ntrgrd);

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_SLIDING) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      const cs_real_t  nf[3] = {fm->face.unitv[0],
                                fm->face.unitv[1],
                                fm->face.unitv[2]};

      /* Compute the product: matpty*face unit normal */

      cs_real_t  pty_nuf[3] = {pty->value * nf[0],
                               pty->value * nf[1],
                               pty->value * nf[2]};

      /* Compute the flux operator related to the trace on the current face
       * of the normal gradient
       * cb->values (size = cm->n_vc) is used inside */

      _vb_ocs_normal_flux_op(f, cm, pty_nuf, dbeta, cb, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT, ">> Cell Vb.COST Sliding bc matrix");
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif

      /* Add contribution to the linear system of the penalized part */

      const double  pcoef = eqp->weak_pena_bc_coeff/sqrt(cm->face[f].meas);

      for (int vfi = 0; vfi < fm->n_vf; vfi++) {

        const short int  bvi = fm->v_ids[vfi];

        /* Retrieve the diagonal 3x3 block matrix */

        cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, bvi, bvi);

        /* Penalized part (u.n)(v.n) and tangential flux part */

        const cs_real_t  vcoef = pcoef * fm->wvf[vfi];
        const cs_real_t  nii = ntrgrd->val[cm->n_vc*bvi + bvi];
        const cs_real_t  bval = vcoef + 2*nii;

        for (int k = 0; k < 3; k++) {
          bii->val[3*k  ] += bval * nf[0]*nf[k];
          bii->val[3*k+1] += bval * nf[1]*nf[k];
          bii->val[3*k+2] += bval * nf[2]*nf[k];
        }

      } /* Loop on face vertices */

      /* Update the system matrix on the extra-diagonal blocks */

      for (int bvi = 0; bvi < cm->n_vc; bvi++) {
        for (int bvj = 0; bvj < cm->n_vc; bvj++) {

          if (bvi == bvj)
            continue;

          /* Retrieve the 3x3 matrix */

          cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bvi, bvj);
          assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

          const cs_real_t  nij = ntrgrd->val[cm->n_vc*bvi + bvj];
          const cs_real_t  nji = ntrgrd->val[cm->n_vc*bvj + bvi];

          for (int k = 0; k < 3; k++) {
            bij->val[3*k  ] += nf[0]*nf[k] * (nij + nji);
            bij->val[3*k+1] += nf[1]*nf[k] * (nij + nji);
            bij->val[3*k+2] += nf[2]*nf[k] * (nij + nji);
          }

        } /* vfj */
      } /* vfi */

    }  /* Sliding face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-Vb schemes with a
 *          COST/Bubble or Voronoi algorithm. One assumes an Orthogonal
 *          splitting between Consistency/Stabilization parts (OCS).
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_ocs_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_hodge_t                     *hodge,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Initialize the local operator */

      cs_sdm_square_init(cm->n_vc, ntrgrd);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_ocs_normal_flux_op(f, cm, pty_nuf, hodgep->coef, cb, ntrgrd);

      /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */

      cs_sdm_t  *ntrgrd_tr = cb->aux;
      cs_sdm_square_add_transpose(ntrgrd, ntrgrd_tr);

     /* Update RHS according to the add of ntrgrd_tr (cb->aux) */

      cs_sdm_square_matvec(ntrgrd_tr, csys->dir_values, cb->values);
      for (short int v = 0; v < csys->n_dofs; v++)
        csys->rhs[v] += cb->values[v];

      /* Update the RHS and the local system matrix */

      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell Vb COST WeakSym bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif
    }  /* Dirichlet face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Robin BCs.
 *          Case of scalar-valued CDO-Vb schemes with a WBS algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_robin(const cs_equation_param_t      *eqp,
                               const cs_cell_mesh_t           *cm,
                               cs_face_mesh_t                 *fm,
                               cs_hodge_t                     *hodge,
                               cs_cell_builder_t              *cb,
                               cs_cell_sys_t                  *csys)
{
  CS_UNUSED(eqp);
  CS_UNUSED(hodge);

  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Robin BCs */

  if (csys->has_robin == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL);

  /* Robin BC expression: K du/dn + alpha*(u - u0) = g
   * Store x = alpha*u0 + g */

  cs_real_t  *x = cb->values;
  cs_sdm_t  *bc_op = cb->loc;
  cs_sdm_t  *hloc = cb->aux; /* 2D Hodge operator on a boundary face */

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_ROBIN) {

      /* Reset local operator */

      cs_sdm_square_init(csys->n_dofs, bc_op);

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      cs_hodge_compute_wbs_surfacic(fm, hloc);  /* hloc is of size n_vf */

      /* Robin BC expression: K du/dn + alpha*(u - u0) = g */
      /* ------------------------------------------------- */

      const double  alpha = csys->rob_values[3*f];
      const double  u0 = csys->rob_values[3*f+1];
      const double  g = csys->rob_values[3*f+2];

      memset(x, 0, sizeof(cs_real_t)*cm->n_vc);
      for (short int v = 0; v < fm->n_vf; v++) {
        const short int  vi = fm->v_ids[v];
        x[vi] = alpha*u0 + g;
      }

      /* Update the RHS and the local system */

      for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

        const short int  vi = fm->v_ids[vfi];
        const cs_real_t  *hfi = hloc->val + vfi*fm->n_vf;
        cs_real_t  *opi = bc_op->val + vi*bc_op->n_rows;

        for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

          const short int  vj = fm->v_ids[vfj];
          csys->rhs[vi] += hfi[vfj]*x[vj];
          opi[vj] += alpha * hfi[vfj];

        }

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell Vb WBS Robin bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, bc_op);
      }
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, bc_op);

    }  /* Robin face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of CDO-Vb schemes with a WBS algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_weak_dirichlet(const cs_equation_param_t      *eqp,
                                        const cs_cell_mesh_t           *cm,
                                        cs_face_mesh_t                 *fm,
                                        cs_hodge_t                     *hodge,
                                        cs_cell_builder_t              *cb,
                                        cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_wbs_normal_flux_op(fm, cm, pty_nuf, cb, ntrgrd);

      /* Update the RHS and the local system matrix */

#if 1 /* Default choice */
      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);
#else  /* This option seems less robust w.r.t the linear algebra */
      _wbs_nitsche(chi/sqrt(cm->face[f].meas), fm, ntrgrd, cb, csys);
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell Vb WBS Weak bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif
    }  /* Dirichlet face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-Vb schemes with a
 *          WBS algorithm
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_wbs_wsym_dirichlet(const cs_equation_param_t     *eqp,
                                        const cs_cell_mesh_t          *cm,
                                        cs_face_mesh_t                *fm,
                                        cs_hodge_t                    *hodge,
                                        cs_cell_builder_t             *cb,
                                        cs_cell_sys_t                 *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Initialize the local operator */

      cs_sdm_square_init(cm->n_vc, ntrgrd);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vb_wbs_normal_flux_op(fm, cm, pty_nuf, cb, ntrgrd);

      /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) cb->loc
         plays the role of the flux operator */

      cs_sdm_t  *ntrgrd_tr = cb->aux;
      cs_sdm_square_add_transpose(ntrgrd, ntrgrd_tr);

      /* Update RHS according to the add of ntrgrd_tr (cb->aux) */

      cs_sdm_square_matvec(ntrgrd_tr, csys->dir_values, cb->values);
      for (short int v = 0; v < csys->n_dofs; v++)
        csys->rhs[v] += cb->values[v];

      /* Update the RHS and the local system matrix */

#if 1 /* Default choice */
      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);
#else  /* This option seems less robust w.r.t the linear algebra */
      _wbs_nitsche(chi/sqrt(cm->face[f].meas), fm, ntrgrd, cb, csys);
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell Vb WBS WeakSym bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif
    }  /* Dirichlet face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique. Case of CDO-VCb schemes with a WBS algorithm.
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcb_weak_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vcb_wbs_normal_flux_op(fm, cm, pty_nuf, cb, ntrgrd);

      /* Update the RHS and the local system matrix */

#if 1 /* Default choice */
      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);
#else  /* This option seems less robust w.r.t the linear algebra */
      _wbs_nitsche(chi/sqrt(cm->face[f].meas), fm, ntrgrd, cb, csys);
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell VCb Weak bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif
    }  /* Dirichlet face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment. Case of CDO-VCb schemes with
 *          a WBS algorithm
 *          Predefined prototype to match the function pointer
 *          cs_cdo_enforce_bc_t
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcb_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                    const cs_cell_mesh_t           *cm,
                                    cs_face_mesh_t                 *fm,
                                    cs_hodge_t                     *hodge,
                                    cs_cell_builder_t              *cb,
                                    cs_cell_sys_t                  *csys)
{
  assert(csys != NULL);  /* Sanity checks */

  /* Enforcement of the Dirichlet BCs */

  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  assert(cm != NULL && cb != NULL && hodge != NULL);

  const cs_property_data_t  *pty = hodge->pty_data;
  const double  chi =
    eqp->weak_pena_bc_coeff * fabs(pty->eigen_ratio)*pty->eigen_max;

  assert(pty->need_tensor);

  cs_sdm_t  *ntrgrd = cb->loc;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */

    const short int  f = csys->_f_ids[i];

    if (cs_cdo_bc_is_dirichlet(csys->bf_flag[f])) {

      /* Compute the face-view of the mesh */

      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */

      cs_real_3_t  pty_nuf;
      cs_math_33_3_product(pty->tensor, fm->face.unitv, pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */

      _vcb_wbs_normal_flux_op(fm, cm, pty_nuf, cb, ntrgrd);

      /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) cb->loc
         plays the role of the flux operator */

      cs_sdm_t  *ntrgrd_tr = cb->aux;
      cs_sdm_square_add_transpose(ntrgrd, ntrgrd_tr);

      /* Update RHS according to the add of ntrgrd_tr (cb->aux) */

      cs_sdm_square_matvec(ntrgrd_tr, csys->dir_values, cb->values);
      for (short int v = 0; v < csys->n_dofs; v++)
        csys->rhs[v] += cb->values[v];

      /* Update the RHS and the local system matrix */

#if 1 /* Default choice */
      _svb_nitsche(chi/sqrt(fm->face.meas), fm, ntrgrd, csys);
#else  /* This option seems less robust w.r.t the linear algebra */
      _wbs_nitsche(chi/sqrt(cm->face[f].meas), fm, ntrgrd, cb, csys);
#endif

      /* Add contribution to the linear system */

      cs_sdm_add(csys->mat, ntrgrd);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      ">> Cell VCb WeakSym bc matrix (f_id: %d)", fm->f_id);
        cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
      }
#endif
    }  /* Dirichlet face */
  } /* Loop on boundary faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across primal faces for a given cell.
 *          Use the same consistent approximation as in the discrete Hodge op.
 *          for this computation.
 *          This function is dedicated to scalar-valued face-based schemes.
 *                       Flux = -Consistent(Hdg) * GRAD(pot)
 *          Predefined prototype to match the function pointer
 *          cs_cdo_diffusion_cw_flux_t
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in]      hodge   pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across primal faces
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_get_face_flux(const cs_cell_mesh_t      *cm,
                                   const double              *pot,
                                   const cs_hodge_t          *hodge,
                                   cs_cell_builder_t         *cb,
                                   double                    *flx)
{
  if (flx == NULL)
    return;

  /* Cellwise DoFs related to the discrete gradient (size: n_fc)
   * flux(f,c) = - Hodge * gradient */

  double  *grd_f = cb->values;
  for (int f = 0; f < cm->n_fc; f++)
    grd_f[f] = cm->f_sgn[f]*(pot[cm->n_fc] - pot[f]);

  /* Store the local fluxes. flux = -Hdg * gradient (grd_f = -gradient)
   * hodge->matrix has been computed just before the call to this function */

  cs_sdm_square_matvec(hodge->matrix, grd_f, flx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell.
 *          Use the same consistent approximation as in the discrete Hodge op.
 *          for this computation.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Consistent(Hdg) * GRAD(pot)
 *          Predefined prototype to match the function pointer
 *          cs_cdo_diffusion_cw_flux_t
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in]      hodge   pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_get_dfbyc_flux(const cs_cell_mesh_t      *cm,
                                    const double              *pot,
                                    const cs_hodge_t          *hodge,
                                    cs_cell_builder_t         *cb,
                                    double                    *flx)
{
  if (flx == NULL)
    return;

  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_EV));

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */

  double  *gec = cb->values;
  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;

    /* sgn_v2 = -sgn_v1; flux = - HDG * GRAD(P) */

    gec[e] = cm->e2v_sgn[e]*(pot[v[1]] - pot[v[0]]);

  }  /* Loop on cell edges */

  /* Store the local fluxes. flux = -Hdg * grd_c(pdi_c)
   * hodge->matrix has been computed just before the call to this function */

  cs_sdm_square_matvec(hodge->matrix, gec, flx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the constant approximation of the diffusive flux inside a
 *          (primal) cell. Use the same consistent approximation as in the
 *          discrete Hodge op. for this computation. This function is dedicated
 *          to vertex-based schemes.
 *          Flux = -Hdg * GRAD(pot)
 *          Predefined prototype to match the function pointer
 *          cs_cdo_diffusion_cw_flux_t
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in]      hodge   pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux inside the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_get_cell_flux(const cs_cell_mesh_t      *cm,
                                   const double              *pot,
                                   const cs_hodge_t          *hodge,
                                   cs_cell_builder_t         *cb,
                                   double                    *flx)
{
  CS_UNUSED(cb);

  if (flx == NULL)
    return;

  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->need_tensor);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_EV | CS_FLAG_COMP_DFQ));

  cs_real_t  grd[3] = {0., 0., 0.};

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */

  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;

    /* sgn_v1 = -sgn_v0; flux = - Kc * GRAD(P) */

    const double  ge = cm->e2v_sgn[e]*(pot[v[1]] - pot[v[0]]);
    const double  contrib = ge * cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      grd[k] += contrib * cm->dface[e].unitv[k];

  }  /* Loop on cell edges */

  cs_math_33_3_product(pty->tensor, grd, flx);

  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) flx[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the normal flux for a face assuming only the knowledge of
 *         the potential at cell vertices. Valid for algorithm relying on a
 *         spliting of the consistency/stabilization part as in OCS (also
 *         called CO+ST) or Bubble algorithm. This is used for reconstructing
 *         the normal flux from the degrees of freedom. The contribution for
 *         each vertex of the face is then computed.
 *
 * \param[in]      f       face id in the cell mesh
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     array of values of the potential (all the mesh)
 * \param[in]      hodge   pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb      auxiliary structure for building the flux
 * \param[in, out] flux    array of values to set (size: n_vc)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_svb_vbyf_flux(short int                   f,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               const cs_hodge_t           *hodge,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux)
{
  if (flux == NULL)
    return;

  const cs_property_data_t  *pty = hodge->pty_data;
  const cs_hodge_param_t  *hodgep = hodge->param;

  assert(pty->need_tensor);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_EV |
                       CS_FLAG_COMP_DFQ | CS_FLAG_COMP_FE));

  const cs_real_t  beta =
    (hodgep->algo == CS_HODGE_ALGO_BUBBLE) ? hodgep->coef : 3*hodgep->coef;
  const cs_quant_t  pfq = cm->face[f];

  /* Reset the fluxes */

  memset(flux, 0, cm->n_vc*sizeof(cs_real_t));

  /* Compute the product: matpty*face unit normal */

  cs_real_t  pty_nuf[3] = {0, 0, 0};;
  cs_math_33_3_product(pty->tensor, pfq.unitv, pty_nuf);

  /* Cellwise constant and consistent gradient */

  cs_real_t  grd_cc[3] = {0, 0, 0};
  cs_real_t  *g = cb->values;

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */

  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;

    /* sgn_v1 = -sgn_v0; GRAD(P) */

    g[e] = cm->e2v_sgn[e]*(pot[v[0]] - pot[v[1]]);

    const double  ge_coef = g[e] * cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      grd_cc[k] += ge_coef * cm->dface[e].unitv[k];

  } /* Loop on cell edges */

  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) grd_cc[k] *= invvol;

  /* Add the stabilisation part which is constant on p_{e,c} --> t_ef if one
     restricts to the face f */

  for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {

    cs_real_t  grd_tef[3] = {0, 0, 0};

    const short int  e = cm->f2e_ids[ie];
    const short int  *v = cm->e2v_ids + 2*e;
    const cs_quant_t  peq = cm->edge[e];
    const cs_nvec3_t  dfq = cm->dface[e];
    const cs_real_t  pec_coef = beta/(peq.meas*_dp3(peq.unitv, dfq.unitv));
    const cs_real_t  delta = g[e] - peq.meas*_dp3(peq.unitv, grd_cc);
    const cs_real_t  stab_coef = pec_coef * delta;

    for (int k = 0; k < 3; k++)
      grd_tef[k] = grd_cc[k] + stab_coef * dfq.unitv[k];

    const cs_real_t  tef = (cs_eflag_test(cm->flag, CS_FLAG_COMP_FEQ)) ?
      cm->tef[ie] : cs_compute_area_from_quant(peq, pfq.center);

    const double  _flx = -0.5 * tef * _dp3(grd_tef, pty_nuf);

    flux[v[0]] += _flx;
    flux[v[1]] += _flx;

  } /* Loop on face edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *          Predefined prototype to match the function pointer
 *          cs_cdo_diffusion_cw_flux_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in]      hodge    pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_dfbyc_flux(const cs_cell_mesh_t   *cm,
                                    const cs_real_t        *pot,
                                    const cs_hodge_t       *hodge,
                                    cs_cell_builder_t      *cb,
                                    cs_real_t              *flx)
{
  /* Temporary buffers */

  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];
  const cs_property_data_t  *pty = hodge->pty_data;

  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_SEF));
  assert(pty->need_tensor);

  /* Reset local fluxes */

  for (short int e = 0; e < cm->n_ec; e++) flx[e] = 0.;

  /* Store segments xv --> xc for this cell */

  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_real_3_t  grd_c, grd_v1, grd_v2, grd_pef, mgrd;

    /* Compute for the current face:
       - the area of each triangle defined by a base e and an apex f
       - the gradient of the Lagrange function related xc in p_{f,c} */

    cs_compute_grdfc_cw(f, cm, grd_c);

    /* Compute the reconstructed value of the potential at p_f */

    double  p_f = 0.;
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];

      p_f += cm->tef[i]*(  p_v[cm->e2v_ids[2*e]]       /* p_v1 */
                         + p_v[cm->e2v_ids[2*e+1]] );  /* p_v2 */
    }
    p_f *= 0.5/cm->face[f].meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */

    const cs_nvec3_t  deq = cm->dedge[f];
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  ee = 2*e;
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      /* Gradient of the lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c} */

      const double  dp1f = p_v[v1] - p_f, dp2f = p_v[v2] - p_f;
      for (int k = 0; k < 3; k++)
        grd_pef[k] = dp_cf*grd_c[k]  + dp1f*grd_v1[k] + dp2f*grd_v2[k];

      cs_math_33_3_product(pty->tensor, grd_pef, mgrd);

      flx[e] -= cm->sefc[i].meas * _dp3(cm->sefc[i].unitv, mgrd);

    }  /* Loop on face edges */

  }  /* Loop on cell faces */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux inside a given primal cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *          Predefined prototype to match the function pointer
 *          cs_cdo_diffusion_cw_flux_t
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in]      hodge    pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux vector inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_cell_flux(const cs_cell_mesh_t   *cm,
                                   const cs_real_t        *pot,
                                   const cs_hodge_t       *hodge,
                                   cs_cell_builder_t      *cb,
                                   cs_real_t              *flx)
{
  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->need_tensor);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                       CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV  | CS_FLAG_COMP_HFQ));

  cs_real_t  cgrd[3] = {0, 0, 0};

  /* Compute the mean-value of the cell gradient */

  cs_reco_cw_cgrd_wbs_from_pvc(cm, pot, cb, cgrd);

  cs_math_33_3_product(pty->tensor, cgrd, flx);
  for (int k = 0; k < 3; k++) flx[k] *= -1;  /* Flux = - tensor * grd */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices (and at cell center).
 *          WBS algorithm is used for reconstructing the normal flux from the
 *          degrees of freedom.
 *
 * \param[in]      f          face id in the cell mesh
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      pot        array of values of the potential (all the mesh)
 * \param[in]      hodge      pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb         auxiliary structure dedicated to diffusion
 * \param[in, out] vf_flux    array of values to set (size: n_vc)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_vbyf_flux(short int                   f,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               const cs_hodge_t           *hodge,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux)
{
  if (flux == NULL)
    return;

  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->need_tensor);
  assert(hodge->param->algo == CS_HODGE_ALGO_WBS);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_FV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ |
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_FE));

  /* Reset the fluxes */

  memset(flux, 0, cm->n_vc*sizeof(cs_real_t));

  /* Compute the product: matpty*face unit normal */

  cs_real_t  mnuf[3] = {0, 0, 0};
  cs_math_33_3_product(pty->tensor, cm->face[f].unitv, mnuf);

  /* Compute xc --> xv length and unit vector for all face vertices */

  double  *l_vc = cb->values;
  cs_real_3_t  *u_vc = cb->vectors;
  for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++) {
    short int  v = cm->f2v_ids[i];
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);
  }

  /* Compute for the current face, the gradient of the Lagrange function
     related to xc in p_{f,c} */

  cs_real_t  grd_c[3], grd_pef[3], grd_v0[3], grd_v1[3];
  cs_compute_grdfc_cw(f, cm, grd_c);

  const cs_real_t  *p_v = pot;
  const cs_real_t  p_f = cs_reco_cw_scalar_pv_at_face_center(f, cm, p_v);

  /* Compute p_c - p_f (where p_c is the reconstructed values at the
     cell center */

  const double  dp_cf = pot[cm->n_vc] - p_f;

  /* Loop on face edges to scan p_{ef,c} subvolumes */

  for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {

    const short int  *v = cm->e2v_ids + 2*cm->f2e_ids[ie];

    /* Compute the gradient of the Lagrange function related xv0, xv1
       in each p_{ef,c} for e in E_f */

    cs_compute_grd_ve(v[0], v[1], cm->dedge[f],
                      (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v0, grd_v1);

    /* Gradient of the lagrange function related to a face.
       grd_f = -(grd_c + grd_v1 + grd_v2)
       This formula is a consequence of the Partition of the Unity.
       This yields the following formula for grd(Lv^conf)|_p_{ef,c} */

    const double  dp_v0f = p_v[v[0]] - p_f;
    const double  dp_v1f = p_v[v[1]] - p_f;
    for (int k = 0; k < 3; k++)
      grd_pef[k] =  dp_cf*grd_c[k] + dp_v0f*grd_v0[k] + dp_v1f*grd_v1[k];

    /* tef: area of the triangle defined by the base e and the apex f */

    const double  _flx = -0.5*cm->tef[ie] * _dp3(mnuf, grd_pef);

    flux[v[0]] += _flx;
    flux[v[1]] += _flx;

  }  /* Loop on face edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across a face (based on a subdivision
 *          into tetrahedra of the volume p_{f,c})
 *
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       pty_tens  3x3 matrix related to the diffusion property
 * \param[in]       p_v       array of values attached to face vertices
 * \param[in]       p_f       value attached to the face
 * \param[in]       p_c       value attached to the cell
 * \param[in, out]  cb        auxiliary structure dedicated to diffusion
 *
 * \return the value of the diffusive flux across the current face
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_diffusion_wbs_face_flux(const cs_face_mesh_t      *fm,
                               const cs_real_t            pty_tens[3][3],
                               const double              *p_v,
                               const double               p_f,
                               const double               p_c,
                               cs_cell_builder_t         *cb)
{
  cs_real_t  grd_c[3], grd_pef[3], grd_v1[3], grd_v2[3], mnuf[3] = {0, 0, 0};
  double  f_flux = 0.;

  /* Retrieve temporary buffers */

  double  *l_vc = cb->values;
  cs_real_3_t  *u_vc = cb->vectors;

  cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, fm->face.unitv,
                       mnuf);

  /* Compute xc --> xv length and unit vector for all face vertices */

  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute for the current face, the gradient of the Lagrange function
     related to xc in p_{f,c} */

  cs_compute_grdfc_fw(fm, grd_c);

  /* Compute p_c - p_f (where p_c is the reconstructed values at the
     cell center */

  const double  dp_cf = p_c - p_f;

  /* Loop on face edges to scan p_{ef,c} subvolumes */

  for (int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* Compute the gradient of the Lagrange function related xv1, xv2
       in each p_{ef,c} for e in E_f */

    cs_compute_grd_ve(v1, v2, fm->dedge, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the lagrange function related to a face.
       grd_f = -(grd_c + grd_v1 + grd_v2)
       This formula is a consequence of the Partition of the Unity.
       This yields the following formula for grd(Lv^conf)|_p_{ef,c} */

    const double  dp1f = p_v[v1] - p_f, dp2f = p_v[v2] - p_f;
    for (int k = 0; k < 3; k++)
      grd_pef[k] = dp_cf*grd_c[k] + dp1f*grd_v1[k] + dp2f*grd_v2[k];

    /* Area of the triangle defined by the base e and the apex f */

    f_flux -= fm->tef[e] * _dp3(mnuf, grd_pef);

  }  /* Loop on face edges */

  return f_flux;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal diffusive flux for a face assuming only the
 *          knowledge of the potential at faces and cell.
 *          CO+ST algorithm is used for reconstructing the normal flux from
 *          the degrees of freedom.
 *
 * \param[in]      f       face id in the cell mesh
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     array of values of the potential (all the mesh)
 * \param[in]      hodge   pointer to a \ref cs_hodge_t structure
 * \param[in, out] cb      auxiliary structure dedicated to diffusion
 * \param[out]     flux    pointer to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_sfb_cost_flux(short int                   f,
                               const cs_cell_mesh_t       *cm,
                               const cs_real_t            *pot,
                               const cs_hodge_t           *hodge,
                               cs_cell_builder_t          *cb,
                               cs_real_t                  *flux)
{
  if (flux == NULL)
    return;

  const cs_property_data_t  *pty = hodge->pty_data;

  assert(pty->need_tensor);
  assert(hodge->param->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  const cs_real_t  beta = hodge->param->coef;
  const cs_quant_t  pfq = cm->face[f];

  /* Compute the product: matpty*face unit normal */

  cs_real_t  pty_nuf[3] = {0, 0, 0};
  cs_math_33_3_product(pty->tensor, pfq.unitv, pty_nuf);

  /* Define a cellwise constant and consistent gradient */
  /* -------------------------------------------------- */

  cs_real_t  grd_cc[3] = {0, 0, 0};
  cs_real_t  *g = cb->values;

  /* Cellwise DoFs related to the discrete gradient (size: n_fc) */

  for (short int ff = 0; ff < cm->n_fc; ff++) {

    /* Delta of the potential along the dual edge (x_f - x_c) */

    g[ff] = cm->f_sgn[ff]*(pot[ff] - pot[cm->n_fc]);

    const double  gcoef = g[ff] * cm->face[ff].meas;
    for (int k = 0; k < 3; k++)
      grd_cc[k] += gcoef * cm->face[ff].unitv[k];

  } /* Loop on cell faces */

  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) grd_cc[k] *= invvol;

  /* Add the stabilisation part which is constant on p_{f,c} */

  const cs_nvec3_t  *deq = cm->dedge + f;
  const cs_real_t  pfc_coef = 3*beta/(_dp3(pfq.unitv, deq->unitv));
  const cs_real_t  delta = g[f] - deq->meas*_dp3(deq->unitv, grd_cc);
  const cs_real_t  stab_coef = pfc_coef * delta;

  cs_real_t  grd_reco[3] = {0., 0., 0.};
  for (int k = 0; k < 3; k++)
    grd_reco[k] = grd_cc[k] + stab_coef * pfq.unitv[k];

  *flux = -pfq.meas * _dp3(grd_reco, pty_nuf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing a piecewise constant gradient from the degrees of
 *          freedom.
 *
 * \param[in]      f            face id in the cell mesh
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_tensor  property tensor times the face normal
 * \param[in]      pot_values   array of values of the potential (all the mesh)
 * \param[in, out] fluxes       values of the fluxes related to each vertex
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_p0_face_flux(const short int           f,
                                const cs_cell_mesh_t     *cm,
                                const cs_real_t           diff_tensor[3][3],
                                const cs_real_t          *pot_values,
                                cs_real_t                *fluxes)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ));

  cs_real_3_t  mnuf;
  cs_real_3_t  gc = {0, 0, 0};

  cs_math_33_3_product(diff_tensor, cm->face[f].unitv, mnuf);

  cs_reco_dfbyc_in_cell(cm, pot_values, gc);

  for (short int v = 0; v < cm->n_vc; v++)  fluxes[v] = 0;

  const cs_real_t  flux_coef = 0.5 * _dp3(gc, mnuf);
  for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {
    const cs_real_t  _flx = flux_coef * cm->tef[i];
    const int  eshft = 2*cm->f2e_ids[i];

    fluxes[cm->e2v_ids[eshft]]   -= _flx;
    fluxes[cm->e2v_ids[eshft+1]] -= _flx;
  }
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
