/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*============================================================================
 * Local variables
 *============================================================================*/

// Advanced developper parameters (weakly enforced the boundary conditions)
static const double  cs_nitsche_pena_coef = 500;
static const double  cs_big_pena_coef = 1e13;

/*============================================================================
 * Private constant variables
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Using the local (cellwise) "normal trace gradient" matrix takes
 *          into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique
 *
 * \param[in]       pcoef     value of the penalization coefficient
 * \param[in]       h_info    cs_param_hodge_t structure for diffusion
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       cbc       pointer to a cs_cell_bc_t structure
 * \param[in, out]  ntrgrd    pointer to a local matrix structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_nitsche(const double              pcoef,
                 const cs_param_hodge_t    h_info,
                 const cs_face_mesh_t     *fm,
                 const cs_cell_bc_t       *cbc,
                 cs_sdm_t                 *ntrgrd,
                 cs_cell_builder_t        *cb,
                 cs_cell_sys_t            *csys)
{
  /* Sanity checks */
  assert(pcoef > 0);
  assert(csys->mat->n_rows == ntrgrd->n_rows);

  const short int  n_csys = csys->mat->n_rows;

  switch (h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
  case CS_PARAM_HODGE_ALGO_AUTO:

    for (short int v = 0; v < fm->n_vf; v++) {

      const short int  vi = fm->v_ids[v];

      // Set the penalty diagonal coefficient
      const double  pcoef_v = pcoef * fm->wvf[v];

      ntrgrd->val[vi + vi*ntrgrd->n_rows] += pcoef_v;
      csys->rhs[vi] += pcoef_v * cbc->dir_values[vi];

    } // Dirichlet or homogeneous Dirichlet
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    { /* Build the related border Hodge operator */
      cs_sdm_t  *hloc = cb->aux;

      cs_hodge_compute_wbs_surfacic(fm, hloc); // hloc is of size n_vf

      /* Add the border Hodge op. to the normal trace op.
         Update RHS whith H*p^{dir} */
      for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

        const double  *hi = hloc->val + vfi*fm->n_vf;
        const short int  vi = fm->v_ids[vfi];

        double  *ntrg_i = ntrgrd->val + vi*n_csys;

        for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

          const double  pcoef_ij = pcoef * hi[vfj];
          const short int  vj = fm->v_ids[vfj];

          ntrg_i[vj] += pcoef_ij;
          csys->rhs[vi] += pcoef_ij * cbc->dir_values[vj];

        } // Loop on face vertices vj

      } // Loop on face vertices vi

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet.");

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Local weak bc matrix (f_id: %d)", fm->f_id);
  cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, ntrgrd);
#endif

  /* Add contribution to the linear system */
  cs_sdm_add(csys->mat, ntrgrd);

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *          Specific to CDO-V+C schemes
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in]      beta      not useful here (prototype of function pointer)
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_diffusion_flux_op(const cs_face_mesh_t     *fm,
                            const cs_cell_mesh_t     *cm,
                            const cs_real_3_t         pty_nuf,
                            double                    beta,
                            cs_cell_builder_t        *cb,
                            cs_sdm_t                 *ntrgrd)
{
  CS_UNUSED(beta);

  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in cb->values and cb->tmp-vect */
  cs_real_t  *l_vc = cb->values;
  cs_real_3_t  *mng_ef = cb->vectors;
  cs_real_3_t  *u_vc = cb->vectors + fm->n_vf;

  const cs_quant_t  pfq = fm->face;
  const cs_nvec3_t  deq = fm->dedge;

  /* Initialize the local operator */
  cs_sdm_square_init(cm->n_vc + 1, ntrgrd);

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdfc(fm->f_sgn, pfq, deq, grd_c);

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); // (pty_tensor * nu_f).grd_c

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute a weight for each vertex of the current face */
  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity */
    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = fm->tef[e] * cs_math_onethird;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

  } /* End of loop on face edges */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_i = ntrgrd->val + vi*ntrgrd->n_rows;

    /* Contribution to the cell column */
    ntrgrd_i[cm->n_vc] = fm->wvf[vfi] * pfq.meas * mng_cf;

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

      } // Loop on face edges

      ntrgrd_i[fm->v_ids[vfj]] += entry_ij;

    } // Loop on face vertices (vj)

  } // Loop on face vertices (vi)

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in]      beta      not useful here (prototype of function pointer)
 * \param[in, out] cb        pointer to a cell builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op. i.e.
 *                           the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_wbs_flux_op(const cs_face_mesh_t     *fm,
                               const cs_cell_mesh_t     *cm,
                               const cs_real_3_t         pty_nuf,
                               double                    beta,
                               cs_cell_builder_t        *cb,
                               cs_sdm_t                 *ntrgrd)
{
  CS_UNUSED(beta);

  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in cb->values and cb->vectors */
  cs_real_t  *l_vc = cb->values;
  cs_real_3_t  *mng_ef = cb->vectors;
  cs_real_3_t  *u_vc = cb->vectors + fm->n_vf;

  /* Initialize the local operator */
  cs_sdm_square_init(cm->n_vc, ntrgrd);

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  const cs_quant_t  pfq = fm->face;
  const cs_nvec3_t  deq = fm->dedge;

  cs_compute_grdfc(fm->f_sgn, pfq, deq, grd_c);

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); // (pty_tensor * nu_f).grd_c

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute a weight for each vertex of the current face */
  double  sum_ef = 0.;
  for (short int e = 0; e < fm->n_ef; e++) {

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(fm->e2v_ids[2*e],
                      fm->e2v_ids[2*e+1],
                      deq,
                      (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity */
    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = fm->tef[e] * cs_math_onethird;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

    sum_ef += mng_ef[e][2];

  } /* End of loop on face edges */

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_i = ntrgrd->val + vi*cm->n_vc;

    /* Default contribution for this line */
    const double  default_coef = pfq.meas * fm->wvf[vfi] * mng_cf;
    for (short int vj = 0; vj < cm->n_vc; vj++)
      ntrgrd_i[vj] = default_coef * cm->wvc[vj]; // two contributions

    /* Block Vf x Vf */
    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

      short int vj = fm->v_ids[vfj];
      ntrgrd_i[vj] += sum_ef * fm->wvf[vfi] * fm->wvf[vfj];

      double  entry_ij = 0.;
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->e2v_ids[2*e];
        const short int  v2 = fm->e2v_ids[2*e+1];

        if (vfj == v1)
          entry_ij += mng_ef[e][0] * fm->wvf[vfi];
        if (vfj == v2)
          entry_ij += mng_ef[e][1] * fm->wvf[vfi];

        if (vfi == v1 || vfi == v2) { /* i in e */
          entry_ij += fm->wvf[vfj] * mng_ef[e][2];
          if (vfj == v1) /* j is also in e */
            entry_ij += mng_ef[e][0];
          if (vfj == v2)
            entry_ij += mng_ef[e][1];
        }

      } // Loop on face edges

      ntrgrd_i[vj] += entry_ij;

    } // Loop on face vertices (vj)

  } // Loop on face vertices (vi)

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id,cm->f_ids[fm->f_id]);
  cs_sdm_dump(cm->c_id,  NULL, NULL, ntrgrd);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          COST algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm      pointer to a cs_face_mesh_t structure
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu     property tensor times the face normal
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] cb      pointer to a cell builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op. i.e.
 *                         the flux operator
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_cost_flux_op(const cs_face_mesh_t     *fm,
                                const cs_cell_mesh_t     *cm,
                                const cs_real_3_t         mnu,
                                double                    beta,
                                cs_cell_builder_t        *cb,
                                cs_sdm_t                 *ntrgrd)
{
  cs_real_3_t  lek;
  cs_real_3_t  *le_grd = cb->vectors;

  const cs_real_t  over_vol_c = 1/cm->vol_c;

  /* Initialize the local operator */
  cs_sdm_square_init(cm->n_vc, ntrgrd);

  /* Loop on border face edges */
  for (short int e = 0; e < fm->n_ef; e++) {

    const cs_quant_t  peq_i = fm->edge[e];
    const short int  vi1 = fm->v_ids[fm->e2v_ids[2*e]];
    const short int  vi2 = fm->v_ids[fm->e2v_ids[2*e+1]];
    const short int  ei = fm->e_ids[e]; // edge id in the cell mesh
    const cs_nvec3_t  dfq_i = cm->dface[ei];
    const double  dp = _dp3(peq_i.unitv, dfq_i.unitv);
    const double  tmp_val = peq_i.meas * dfq_i.meas * dp; // 3*pec_vol
    const double  beta_pec_vol = 3. * beta/tmp_val;
    const double  coef_ei = 3. * beta/dp;

    /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
    for (short int v = 0; v < cm->n_vc; v++)
      le_grd[v][0] = le_grd[v][1] = le_grd[v][2] = 0;

    /* Term related to the flux reconstruction:
       Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell */
    for (short int ek = 0; ek < cm->n_ec; ek++) {

      const short int  vj1 = cm->e2v_ids[2*ek], vj2 = cm->e2v_ids[2*ek+1];
      const short int  sgn_vj1 = cm->e2v_sgn[ek]; // sgn_vj2 = - sgn_vj1
      const cs_nvec3_t  dfq_k = cm->dface[ek];

      /* Compute l_(ek,c)|p(ei,f,c) */
      const double  eik_part = coef_ei * _dp3(dfq_k.unitv, peq_i.unitv);
      const double  coef_mult = dfq_k.meas * over_vol_c;

      for (int k = 0; k < 3; k++)
        lek[k] = coef_mult * (dfq_k.unitv[k] - eik_part * dfq_i.unitv[k]);

      if (ek == ei)
        for (int k = 0; k < 3; k++)
          lek[k] += dfq_k.meas * dfq_k.unitv[k] * beta_pec_vol;

      for (int k = 0; k < 3; k++) {
        le_grd[vj1][k] += sgn_vj1 * lek[k];
        le_grd[vj2][k] -= sgn_vj1 * lek[k];
      }

    } // Loop on cell edges

    for (short int v = 0; v < cm->n_vc; v++) {

      const double  contrib = _dp3(mnu, le_grd[v]) * fm->tef[e];
      ntrgrd->val[vi1*cm->n_vc + v] += contrib;
      ntrgrd->val[vi2*cm->n_vc + v] += contrib;

    }

  } // Loop on face edges

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)",
                cm->c_id, cm->f_ids[fm->f_id]);
  cs_sdm_dump(cm->c_id, NULL, NULL, ntrgrd);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique
 *
 * \param[in]       h_info    cs_param_hodge_t structure for diffusion
 * \param[in]       cbc       pointer to a cs_cell_bc_t structure
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_weak_dirichlet(const cs_param_hodge_t          h_info,
                                  const cs_cell_bc_t             *cbc,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  /* Sanity checks */
  assert(cbc != NULL && cm != NULL);
  assert(cb != NULL && csys != NULL);

  /* For VCB schemes cbc->n_dofs = cm->n_vc and csys->mat->n_rows = cm->n_vc + 1
     For VB schemes  cbc->n_dofs = csys->mat->n_rows = cm->n_vc */

  /* Enforcement of the Dirichlet BCs */
  if (cbc->n_dirichlet == 0)
    return; // Nothing to do

  const double chi = cs_nitsche_pena_coef * fabs(cb->eig_ratio) * cb->eig_max;

  for (short int i = 0; i < cbc->n_bc_faces; i++) {
    if (cbc->face_flag[i] & CS_CDO_BC_DIRICHLET ||
        cbc->face_flag[i] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute the face-view of the mesh */
      cs_face_mesh_build_from_cell_mesh(cm, cbc->bf_ids[i], fm);

      /* Compute the product: matpty*face unit normal */
      cs_real_3_t  pty_nuf;
      cs_math_33_3_product((const cs_real_t (*)[3])cb->pty_mat,
                           fm->face.unitv,
                           pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */
      flux_op(fm, cm, pty_nuf, h_info.coef, cb, cb->loc);

      /* Update the RHS and the local system matrix */
      _enforce_nitsche(chi/sqrt(fm->face.meas), // Penalization coeff.
                       h_info,
                       fm, cbc,
                       cb->loc, cb, csys);

    } // Dirichlet face
  } /* Loop on boundary faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment
 *
 * \param[in]       h_info    cs_param_hodge_t structure for diffusion
 * \param[in]       cbc       pointer to a cs_cell_bc_t structure
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_wsym_dirichlet(const cs_param_hodge_t           h_info,
                                  const cs_cell_bc_t              *cbc,
                                  const cs_cell_mesh_t            *cm,
                                  cs_cdo_diffusion_flux_trace_t   *flux_op,
                                  cs_face_mesh_t                  *fm,
                                  cs_cell_builder_t               *cb,
                                  cs_cell_sys_t                   *csys)
{
  /* Sanity checks */
  assert(cbc != NULL && cm != NULL);
  assert(cb != NULL && csys != NULL);
  /* For VCB schemes cbc->n_dofs = cm->n_vc and csys->mat->n_rows = cm->n_vc + 1
     For VB schemes  cbc->n_dofs = csys->mat->n_rows = cm->n_vc */

  /* Enforcement of the Dirichlet BCs */
  if (cbc->n_dirichlet == 0)
    return; // Nothing to do

  const double chi = cs_nitsche_pena_coef * cb->eig_ratio * cb->eig_max;

  for (short int i = 0; i < cbc->n_bc_faces; i++) {
    if (cbc->face_flag[i] & CS_CDO_BC_DIRICHLET ||
        cbc->face_flag[i] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute the face-view of the mesh */
      cs_face_mesh_build_from_cell_mesh(cm, cbc->bf_ids[i], fm);

      /* Compute the product: matpty*face unit normal */
      cs_real_3_t  pty_nuf;
      cs_math_33_3_product((const cs_real_t (*)[3])cb->pty_mat,
                           fm->face.unitv,
                           pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */
      flux_op(fm, cm, pty_nuf, h_info.coef, cb, cb->loc);

      /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) cb->loc
         plays the role of the flux operator */
      cs_sdm_square_add_transpose(cb->loc, cb->aux);

      /* Update RHS according to the add of transp (cb->aux) */
      cs_sdm_square_matvec(cb->aux, cbc->dir_values, cb->values);
      for (short int v = 0; v < cbc->n_dofs; v++)
        csys->rhs[v] += cb->values[v];

      /* Update the RHS and the local system matrix */
      _enforce_nitsche(chi/sqrt(fm->face.meas), // Penalization coeff.
                       h_info,
                       fm, cbc,
                       cb->loc, cb, csys);

    } /* This boundary face is attached to a Dirichlet BC */
  } /* Loop on boundary faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value
 *
 * \param[in]       h_info    cs_param_hodge_t structure for diffusion
 * \param[in]       cbc       pointer to a cs_cell_bc_t structure
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_dirichlet(const cs_param_hodge_t           h_info,
                                const cs_cell_bc_t              *cbc,
                                const cs_cell_mesh_t            *cm,
                                cs_cdo_diffusion_flux_trace_t   *flux_op,
                                cs_face_mesh_t                  *fm,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys)
{
  CS_UNUSED(h_info); // Prototype common to cs_cdo_diffusion_enforce_dir_t
  CS_UNUSED(fm);
  CS_UNUSED(cm);
  CS_UNUSED(cb);
  CS_UNUSED(flux_op);

  /* Sanity checks */
  assert(cbc != NULL && cm != NULL);
  assert(csys != NULL);
  /* For VCB schemes cbc->n_dofs = cm->n_vc and csys->mat->n_rows = cm->n_vc + 1
     For VB schemes  cbc->n_dofs = csys->mat->n_rows = cm->n_vc
     For FB schemes  cbc->n_dofs = csys->mat->n_rows = cm->n_fc */

  /* Enforcement of the Dirichlet BCs */
  if (cbc->n_dirichlet == 0)
    return; // Nothing to do

  const short int n_dofs = csys->mat->n_rows;

  // Penalize diagonal entry (and its rhs if needed)
  for (short int i = 0; i < cbc->n_dofs; i++) {

    if (cbc->dof_flag[i] & CS_CDO_BC_DIRICHLET) {
      csys->mat->val[i + n_dofs*i] += cs_big_pena_coef;
      csys->rhs[i] += cbc->dir_values[i] * cs_big_pena_coef;
    }
    else if (cbc->dof_flag[i] & CS_CDO_BC_HMG_DIRICHLET)
      csys->mat->val[i + n_dofs*i] += cs_big_pena_coef;

  } /* Loop on degrees of freedom */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          Use the COST algo. for computing the discrete Hodge op.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Hdg * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcost_get_dfbyc_flux(const cs_cell_mesh_t      *cm,
                                      const double              *pot,
                                      cs_cell_builder_t         *cb,
                                      double                    *flx)
{
  /* Sanity checks */
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_EV));

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */
  double  *gec = cb->values;
  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;
    // sgn_v2 = -sgn_v1; flux = - HDG * GRAD(P)
    gec[e] = cm->e2v_sgn[e]*(pot[v[1]] - pot[v[0]]);

  } // Loop on cell edges

  /* Store the local fluxes. flux = -Hdg * grd_c(pdi_c)
     cb->hdg has been computed just before the call to this function */
  cs_sdm_square_matvec(cb->hdg, gec, flx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux inside a (primal) cell
 *          Use the COST algo. for computing the discrete Hodge op.
 *          This function is dedicated to vertex-based schemes.
 *                       Flux = -Hdg * GRAD(pot)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      pot     values of the potential fields at specific locations
 * \param[in, out] cb      auxiliary structure for computing the flux
 * \param[in, out] flx     values of the flux across specific entities
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_vcost_get_pc_flux(const cs_cell_mesh_t      *cm,
                                   const double              *pot,
                                   cs_cell_builder_t         *cb,
                                   double                    *flx)
{
  /* Sanity checks */
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  cs_real_t  grd[3] = {0., 0., 0.};

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */
  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;

    // sgn_v1 = -sgn_v0; flux = - Kc * GRAD(P)
    const double  ge = cm->e2v_sgn[e]*(pot[v[1]] - pot[v[0]]);
    const double  contrib = ge * cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      grd[k] += contrib * cm->dface[e].unitv[k];

  } // Loop on cell edges

  cs_math_33_3_product((const cs_real_t (*)[3])cb->pty_mat, grd, flx);
  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) flx[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_dfbyc_flux(const cs_cell_mesh_t   *cm,
                                    const cs_real_t        *pot,
                                    cs_cell_builder_t      *cb,
                                    cs_real_t              *flx)
{
  /* Sanity checks */
  assert(cs_test_flag(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_EFQ));

  cs_real_3_t  grd_c, grd_v1, grd_v2, grd_pef, mgrd;

  /* Temporary buffers */
  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];

  /* Reset local fluxes */
  for (short int e = 0; e < cm->n_ec; e++) flx[e] = 0.;

  /* Store segments xv --> xc for this cell */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute for the current face:
       - the area of each triangle defined by a base e and an apex f
       - the gradient of the Lagrange function related xc in p_{f,c} */
    cs_compute_grdfc(cm->f_sgn[f], pfq, deq, grd_c);

    /* Compute the reconstructed value of the potential at p_f */
    double  p_f = 0.;

    /* Loop on face edges */
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];

      p_f += cm->tef[i]*(  p_v[cm->e2v_ids[2*e]]      // p_v1
                         + p_v[cm->e2v_ids[2*e+1]] ); // p_v2
    }
    p_f *= 0.5/pfq.meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */
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
      for (int k = 0; k < 3; k++)
        grd_pef[k] = dp_cf          *grd_c[k]  +
                     (p_v[v1] - p_f)*grd_v1[k] +
                     (p_v[v2] - p_f)*grd_v2[k];

      cs_math_33_3_product((const cs_real_t (*)[3])cb->pty_mat, grd_pef, mgrd);

      if (f == cm->e2f_ids[ee])
        flx[e] -= cm->sefc[ee].meas * _dp3(cm->sefc[ee].unitv, mgrd);
      else {
        assert(f == cm->e2f_ids[ee+1]);
        flx[e] -= cm->sefc[ee+1].meas * _dp3(cm->sefc[ee+1].unitv, mgrd);
      }

    } // Loop on face edges

  } // Loop on cell faces

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux inside a given primal cell
 *          Use the WBS algo. for approximating the gradient
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] flx      flux vector inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_wbs_get_pc_flux(const cs_cell_mesh_t   *cm,
                                 const cs_real_t        *pot,
                                 cs_cell_builder_t      *cb,
                                 cs_real_t              *flx)
{
  /* Sanity checks */
  assert(cs_test_flag(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  cs_real_3_t  cgrd;

  /* Compute the mean-value of the cell gradient */
  cs_reco_cw_cgrd_wbs_from_pvc(cm, pot, cb, cgrd);

  cs_math_33_3_product((const cs_real_t (*)[3])cb->pty_mat, cgrd, flx);
  for (int k = 0; k < 3; k++) flx[k] *= -1; // Flux = - tensor * grd
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
cs_cdo_diffusion_face_flux(const cs_face_mesh_t      *fm,
                           const cs_real_t            pty_tens[3][3],
                           const double              *p_v,
                           const double               p_f,
                           const double               p_c,
                           cs_cell_builder_t         *cb)
{
  cs_real_3_t  grd_c, grd_v1, grd_v2, grd_pef, mnuf;

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
     related xc in p_{f,c} */
  cs_compute_grdfc(fm->f_sgn, fm->face, fm->dedge, grd_c);

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
    for (int k = 0; k < 3; k++)
      grd_pef[k] =  dp_cf          *grd_c[k]  +
                   (p_v[v1] - p_f) *grd_v1[k] +
                   (p_v[v2] - p_f) *grd_v2[k];

    /* Area of the triangle defined by the base e and the apex f */
    f_flux -= fm->tef[e] * _dp3(mnuf, grd_pef);

  } // Loop on face edges

  return f_flux;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
