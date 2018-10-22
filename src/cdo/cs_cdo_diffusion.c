/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 * \brief   Using the local (cellwise) "normal trace gradient" matrix takes
 *          into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique
 *
 * \param[in]       pcoef     value of the penalization coefficient
 * \param[in]       h_info    cs_param_hodge_t structure for diffusion
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  ntrgrd    pointer to a local matrix structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_nitsche(const double              pcoef,
                 const cs_param_hodge_t    h_info,
                 const cs_face_mesh_t     *fm,
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

      /* Set the penalty diagonal coefficient */
      const double  pcoef_v = pcoef * fm->wvf[v];

      ntrgrd->val[vi + vi*ntrgrd->n_rows] += pcoef_v;
      csys->rhs[vi] += pcoef_v * csys->dir_values[vi];

    }  /* Dirichlet or homogeneous Dirichlet */
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    { /* Build the related border Hodge operator */
      cs_sdm_t  *hloc = cb->aux;

      cs_hodge_compute_wbs_surfacic(fm, hloc);  /* hloc is of size n_vf */

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
          csys->rhs[vi] += pcoef_ij * csys->dir_values[vj];

        }  /* Loop on face vertices vj */

      }  /* Loop on face vertices vi */

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute \f$ \int_{fb} \nabla (u) \cdot \nu_{fb} v \f$ where \p fb
 *         is a boundary faces (Co+St algorithm)
 *
 * \param[in]       fb       index of the boundary face on the local numbering
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in]       h_info    \ref cs_param_hodge_t structure for diffusion
 * \param[in]       kappa_f   diffusion property against face vector for all
 *                            faces
 * \param[in, out]  ntrgrd    pointer to a local matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdofb_normal_flux_reco(short int                  fb,
                        const cs_cell_mesh_t      *cm,
                        const cs_cell_builder_t   *cb,
                        const cs_param_hodge_t     h_info,
                        const cs_real_3_t         *kappa_f,
                        cs_sdm_t                  *ntrgrd)
{
  /* Sanity check */
  assert(cb != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_HFQ));
  assert(cm->f_sgn[fb] == 1);  /* +1 because it's a boundary face */

  const short int  nfc = cm->n_fc;
  const cs_quant_t  pfbq = cm->face[fb];
  const cs_nvec3_t  debq = cm->dedge[fb];

  /* |fb|^2 * nu_{fb}^T.kappa.nu_{fb} */
  const cs_real_t  fb_k_fb = pfbq.meas * _dp3(kappa_f[fb], pfbq.unitv);
  const cs_real_t  beta_fbkfb_o_pfc = h_info.coef * fb_k_fb / cm->pfc[fb];

  cs_real_t  *ntrgrd_fb = ntrgrd->val + fb * (nfc + 1);
  cs_real_t  row_sum = 0.0;
  for (short int f = 0; f < nfc; f++) {

    const cs_real_t  if_ov = cm->f_sgn[f] / cm->vol_c;
    const cs_real_t  f_k_fb = pfbq.meas * _dp3(kappa_f[f], pfbq.unitv);
    const cs_quant_t  pfq = cm->face[f];
    cs_real_t  stab = - pfq.meas*debq.meas * _dp3(debq.unitv, pfq.unitv);
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
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique - Face-based version
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_diffusion_weak_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);
  assert(cs_equation_param_has_diffusion(eqp));

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const double chi =
    eqp->bc_penalization_coeff * fabs(cb->eig_ratio)*cb->eig_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */
  cs_real_3_t  *kappa_f = cb->vectors;
  if (h_info.is_unity) {
    for (short int f = 0; f < cm->n_fc; f++) {
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = cm->face[f].meas*cm->face[f].unitv[k];
    }
  }
  else if (h_info.is_iso) {
    for (short int f = 0; f < cm->n_fc; f++) {
      const cs_real_t  coef = cm->face[f].meas*cb->dpty_val;
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = coef * cm->face[f].unitv[k];
    }
  }
  else {
    for (short int f = 0; f < cm->n_fc; f++) {
      const cs_real_t  coef = cm->face[f].meas*cb->dpty_val;
      cs_math_33_3_product((const cs_real_3_t *)cb->dpty_mat, cm->face[f].unitv,
                           kappa_f[f]);
      for (short int k = 0; k < 3; k++) kappa_f[f][k] *= cm->face[f].meas;
    }
  }

  /* Initialize the matrix related this flux reconstruction operator */
  const short int n_dofs = cm->n_fc + 1;
  cs_sdm_t *bc_op = cb->loc;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */
  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute \int_f du/dn v and update the matrix */
      _cdofb_normal_flux_reco(f, cm, cb, h_info, kappa_f, bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Second pass: add the bc_op matrix, add the BC */
  cs_real_t *ntrgrd_v = bc_op->val;

  /* !!! ATTENTION !!!
   * Two passes in order to avoid truncation error if the arbitrary coefficient
   * of the Nitsche algorithm is large
   */
  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

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
 *          technique plus a symmetric treatment - Face-based version
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_diffusion_wsym_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  CS_UNUSED(fm);

  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);
  assert(cs_equation_param_has_diffusion(eqp));

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const double chi =
    eqp->bc_penalization_coeff * fabs(cb->eig_ratio)*cb->eig_max;

  /* First step: pre-compute the product between diffusion property and the
     face vector areas */
  cs_real_3_t  *kappa_f = cb->vectors;
  if (h_info.is_unity) {
    for (short int f = 0; f < cm->n_fc; f++) {
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = cm->face[f].meas*cm->face[f].unitv[k];
    }
  }
  else if (h_info.is_iso) {
    for (short int f = 0; f < cm->n_fc; f++) {
      const cs_real_t  coef = cm->face[f].meas*cb->dpty_val;
      for (short int k = 0; k < 3; k++)
        kappa_f[f][k] = coef * cm->face[f].unitv[k];
    }
  }
  else {
    for (short int f = 0; f < cm->n_fc; f++) {
      const cs_real_t  coef = cm->face[f].meas*cb->dpty_val;
      cs_math_33_3_product((const cs_real_3_t *)cb->dpty_mat, cm->face[f].unitv,
                           kappa_f[f]);
      for (short int k = 0; k < 3; k++) kappa_f[f][k] *= cm->face[f].meas;
    }
  }

  const short int n_dofs = cm->n_fc + 1, n_f = cm->n_fc;
  cs_sdm_t  *bc_op = cb->loc, *bc_op_t = cb->aux;
  cs_sdm_square_init(n_dofs, bc_op);

  /* First pass: build the bc_op matrix */
  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute \int_f du/dn v and update the matrix */
      _cdofb_normal_flux_reco(f, cm, cb, h_info, kappa_f, bc_op);

    } /* If Dirichlet */

  } /* Loop boundary faces */

  /* Second pass: add the bc_op matrix, add the BC */
  /* !!! ATTENTION !!!
   * Two passes in order to avoid truncation error if the arbitrary coefficient
   * of the Nitsche algo is large
   */
  cs_real_t *dir_val = cb->values, *u0_trgradv = cb->values + n_dofs;

  /* Putting the face DoFs of the BC, into a face- and cell-DoFs array */
  memcpy(dir_val, csys->dir_values, n_f*sizeof(cs_real_t));
  dir_val[n_f] = 0.;

  /* Update bc_op = bc_op + transp and transp = transpose(bc_op) cb->loc
     plays the role of the flux operator */
  cs_sdm_square_add_transpose(bc_op, bc_op_t);
  cs_sdm_square_matvec(bc_op_t, dir_val, u0_trgradv);

  for (short int i = 0; i < n_dofs; i++) /* Cell too! */
    csys->rhs[i] += u0_trgradv[i];

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

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

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); /* (pty_tensor*nu_f).grd_c */

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

      }  /* Loop on face edges */

      ntrgrd_i[fm->v_ids[vfj]] += entry_ij;

    }  /* Loop on face vertices (vj) */

  }  /* Loop on face vertices (vi) */

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

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); /* (pty_tensor*nu_f).grd_c */

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
      ntrgrd_i[vj] = default_coef * cm->wvc[vj];  /* two contributions */

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

      }  /* Loop on face edges */

      ntrgrd_i[vj] += entry_ij;

    }  /* Loop on face vertices (vj) */

  }  /* Loop on face vertices (vi) */

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
 * \param[in]      beta    value of the stabilization coef. related to reco.
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
    const short int  ei = fm->e_ids[e];  /* edge id in the cell mesh */
    const cs_nvec3_t  dfq_i = cm->dface[ei];
    const double  dp = _dp3(peq_i.unitv, dfq_i.unitv);
    const double  tmp_val = peq_i.meas * dfq_i.meas * dp;  /* 3*pec_vol */
    const double  beta_pec_vol = 3. * beta/tmp_val;
    const double  coef_ei = 3. * beta/dp;

    /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
    for (short int v = 0; v < cm->n_vc; v++)
      le_grd[v][0] = le_grd[v][1] = le_grd[v][2] = 0;

    /* Term related to the flux reconstruction:
       Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell */
    for (short int ek = 0; ek < cm->n_ec; ek++) {

      const short int  vj1 = cm->e2v_ids[2*ek], vj2 = cm->e2v_ids[2*ek+1];
      const short int  sgn_vj1 = cm->e2v_sgn[ek];  /* sgn_vj2 = - sgn_vj1 */
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

    }  /* Loop on cell edges */

    for (short int v = 0; v < cm->n_vc; v++) {

      const double  contrib = _dp3(mnu, le_grd[v]) * fm->tef[e];
      ntrgrd->val[vi1*cm->n_vc + v] += contrib;
      ntrgrd->val[vi2*cm->n_vc + v] += contrib;

    }

  }  /* Loop on face edges */

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
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cellwise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_weak_dirichlet(const cs_equation_param_t      *eqp,
                                  const cs_cell_mesh_t           *cm,
                                  cs_cdo_diffusion_flux_trace_t  *flux_op,
                                  cs_face_mesh_t                 *fm,
                                  cs_cell_builder_t              *cb,
                                  cs_cell_sys_t                  *csys)
{
  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const double chi = eqp->bc_penalization_coeff*fabs(cb->eig_ratio)*cb->eig_max;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];

    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute the face-view of the mesh */
      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */
      cs_real_3_t  pty_nuf;
      cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat,
                           fm->face.unitv,
                           pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */
      flux_op(fm, cm, pty_nuf, h_info.coef, cb, cb->loc);

      /* Update the RHS and the local system matrix */
      _enforce_nitsche(chi/sqrt(fm->face.meas),  /* Penalization coeff. */
                       h_info,
                       fm,
                       cb->loc, cb, csys);

    }  /* Dirichlet face */
  } /* Loop on boundary faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique plus a symmetric treatment
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a \ref cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a \ref cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_wsym_dirichlet(const cs_equation_param_t       *eqp,
                                  const cs_cell_mesh_t            *cm,
                                  cs_cdo_diffusion_flux_trace_t   *flux_op,
                                  cs_face_mesh_t                  *fm,
                                  cs_cell_builder_t               *cb,
                                  cs_cell_sys_t                   *csys)
{
  /* Sanity checks */
  assert(cm != NULL && cb != NULL && csys != NULL);

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const double chi = eqp->bc_penalization_coeff*fabs(cb->eig_ratio)*cb->eig_max;

  for (short int i = 0; i < csys->n_bc_faces; i++) {

    /* Get the boundary face in the cell numbering */
    const short int  f = csys->_f_ids[i];
    if (csys->bf_flag[f] & CS_CDO_BC_DIRICHLET ||
        csys->bf_flag[f] & CS_CDO_BC_HMG_DIRICHLET) {

      /* Compute the face-view of the mesh */
      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute the product: matpty*face unit normal */
      cs_real_3_t  pty_nuf;
      cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat,
                           fm->face.unitv,
                           pty_nuf);

      /* Compute the flux operator related to the trace on the current face
         of the normal gradient */
      flux_op(fm, cm, pty_nuf, h_info.coef, cb, cb->loc);

      /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) cb->loc
         plays the role of the flux operator */
      cs_sdm_square_add_transpose(cb->loc, cb->aux);

      /* Update RHS according to the add of transp (cb->aux) */
      cs_sdm_square_matvec(cb->aux, csys->dir_values, cb->values);
      for (short int v = 0; v < csys->n_dofs; v++)
        csys->rhs[v] += cb->values[v];

      /* Update the RHS and the local system matrix */
      _enforce_nitsche(chi/sqrt(fm->face.meas),  /* Penalization coeff. */
                       h_info,
                       fm,
                       cb->loc, cb, csys);

    } /* This boundary face is attached to a Dirichlet BC */
  } /* Loop on boundary faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_cdo_diffusion_flux_trace_t   *flux_op,
                                cs_face_mesh_t                  *fm,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys)
{
  /* Prototype common to cs_cdo_diffusion_enforce_dir_t.
     Hence the unused parameters */
  CS_UNUSED(fm);
  CS_UNUSED(cm);
  CS_UNUSED(cb);
  CS_UNUSED(flux_op);

  /* Sanity checks */
  assert(cm != NULL && csys != NULL);

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  /* Penalize diagonal entry (and its rhs if needed) */
  for (short int i = 0; i < csys->n_dofs; i++) {

    if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET) {
      csys->mat->val[i + csys->n_dofs*i] += eqp->bc_penalization_coeff;
      csys->rhs[i] += csys->dir_values[i] * eqp->bc_penalization_coeff;
    }
    else if (csys->dof_flag[i] & CS_CDO_BC_HMG_DIRICHLET)
      csys->mat->val[i + csys->n_dofs*i] += eqp->bc_penalization_coeff;

  } /* Loop on degrees of freedom */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account Dirichlet BCs by a weak enforcement by a
 *          penalization technique with a huge value.
 *          Case of a cellwise system defined by block.
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_pena_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_cdo_diffusion_flux_trace_t   *flux_op,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys)
{
  /* Prototype common to cs_cdo_diffusion_enforce_dir_t
     Hence the unused parameters */
  CS_UNUSED(fm);
  CS_UNUSED(cm);
  CS_UNUSED(cb);
  CS_UNUSED(flux_op);

  /* Sanity checks */
  assert(cm != NULL && csys != NULL);

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

      if (_flag[i] & CS_CDO_BC_DIRICHLET) {
        mII->val[i + mII->n_rows*i] += eqp->bc_penalization_coeff;
        _rhs[i] += _dir_val[i] * eqp->bc_penalization_coeff;
      }
      else if (_flag[i] & CS_CDO_BC_HMG_DIRICHLET)
        mII->val[i + mII->n_rows*i] += eqp->bc_penalization_coeff;

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
 * where xd is the value of the Dirichlet BC
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_dirichlet(const cs_equation_param_t       *eqp,
                                const cs_cell_mesh_t            *cm,
                                cs_cdo_diffusion_flux_trace_t   *flux_op,
                                cs_face_mesh_t                  *fm,
                                cs_cell_builder_t               *cb,
                                cs_cell_sys_t                   *csys)
{
  /* Prototype common to cs_cdo_diffusion_enforce_dir_t
     Hence the unused parameters */
  CS_UNUSED(eqp);
  CS_UNUSED(fm);
  CS_UNUSED(cm);
  CS_UNUSED(flux_op);

  /* Enforcement of the Dirichlet BCs */
  if (csys->has_dirichlet == false)
    return;  /* Nothing to do */

  double  *x_dir = cb->values;
  double  *ax_dir = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  /* Build x_dir */
  for (short int i = 0; i < csys->n_dofs; i++)
    if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
      x_dir[i] = csys->dir_values[i];

  /* Contribution of the Dirichlet conditions */
  cs_sdm_matvec(csys->mat, x_dir, ax_dir);

  /* Second pass: Replace the Dirichlet block by a diagonal block */
  for (short int i = 0; i < csys->n_dofs; i++) {

    if (csys->dof_flag[i] & (CS_CDO_BC_DIRICHLET | CS_CDO_BC_HMG_DIRICHLET)) {

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
 * where xd is the value of the Dirichlet BC
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       flux_op   function pointer to the flux trace operator
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_alge_block_dirichlet(const cs_equation_param_t       *eqp,
                                      const cs_cell_mesh_t            *cm,
                                      cs_cdo_diffusion_flux_trace_t   *flux_op,
                                      cs_face_mesh_t                  *fm,
                                      cs_cell_builder_t               *cb,
                                      cs_cell_sys_t                   *csys)
{
  /* Prototype common to cs_cdo_diffusion_enforce_dir_t
     Hence the unused parameters */
  CS_UNUSED(eqp);
  CS_UNUSED(fm);
  CS_UNUSED(cm);
  CS_UNUSED(flux_op);

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
    if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
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
      if (_flg[i] & (CS_CDO_BC_DIRICHLET | CS_CDO_BC_HMG_DIRICHLET))
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
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_EV));

  /* Cellwise DoFs related to the discrete gradient (size: n_ec) */
  double  *gec = cb->values;
  for (short int e = 0; e < cm->n_ec; e++) {

    const short int  *v = cm->e2v_ids + 2*e;
    /* sgn_v2 = -sgn_v1; flux = - HDG * GRAD(P) */
    gec[e] = cm->e2v_sgn[e]*(pot[v[1]] - pot[v[0]]);

  }  /* Loop on cell edges */

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
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

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

  cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat, grd, flx);
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
  assert(cs_flag_test(cm->flag,
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

      p_f += cm->tef[i]*(  p_v[cm->e2v_ids[2*e]]       /* p_v1 */
                         + p_v[cm->e2v_ids[2*e+1]] );  /* p_v2 */
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

      cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat, grd_pef, mgrd);

      if (f == cm->e2f_ids[ee])
        flx[e] -= cm->sefc[ee].meas * _dp3(cm->sefc[ee].unitv, mgrd);
      else {
        assert(f == cm->e2f_ids[ee+1]);
        flx[e] -= cm->sefc[ee+1].meas * _dp3(cm->sefc[ee+1].unitv, mgrd);
      }

    }  /* Loop on face edges */

  }  /* Loop on cell faces */

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
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  cs_real_3_t  cgrd;

  /* Compute the mean-value of the cell gradient */
  cs_reco_cw_cgrd_wbs_from_pvc(cm, pot, cb, cgrd);

  cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat, cgrd, flx);
  for (int k = 0; k < 3; k++) flx[k] *= -1;  /* Flux = - tensor * grd */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing a piecewise constant gradient from the degrees of
 *          freedom.
 *
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      diff_tensor  property tensor times the face normal
 * \param[in]      pot_values   array of values of the potential (all the mesh)
 * \param[in]      f            face id in the cell mesh
 * \param[in]      t_eval       time at which one evaluates the advection field
 * \param[in, out] fluxes       values of the fluxes related to each vertex
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_diffusion_face_p0_flux(const cs_cell_mesh_t     *cm,
                                const cs_real_3_t        *diff_tensor,
                                const cs_real_t          *pot_values,
                                short int                 f,
                                cs_real_t                 t_eval,
                                cs_real_t                *fluxes)
{
  CS_UNUSED(t_eval);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ));

  cs_real_3_t  mnuf;
  cs_real_3_t  gc = {0, 0, 0};

  cs_math_33_3_product((const cs_real_t (*)[3])diff_tensor, cm->face[f].unitv,
                       mnuf);

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
cs_cdo_diffusion_face_wbs_flux(const cs_face_mesh_t      *fm,
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

  }  /* Loop on face edges */

  return f_flux;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal flux for a face assuming only the knowledge
 *          of the potential at cell vertices. COST algorithm is used for
 *          reconstructing the degrees of freedom.
 *          This routine shares similarities with
 *          \ref cs_cdovb_diffusion_cost_flux_op
 *
 * \param[in]   f             face id in the cell mesh
 * \param[in]   cm            pointer to a cs_cell_mesh_t structure
 * \param[in]   diff_tensor   property tensor times the face normal
 * \param[in]   pot_values    array of values of the potential (all the mesh)
 * \param[in]   beta          value of the stabilization coef. related to reco.
 * \param[in, out]  cb        auxiliary structure dedicated to diffusion
 *
 * \return the diffusive flux across the face f
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdovb_diffusion_face_cost_flux(short int                 f,
                                  const cs_cell_mesh_t     *cm,
                                  const cs_real_3_t        *diff_tensor,
                                  const cs_real_t          *pot_values,
                                  double                    beta,
                                  cs_cell_builder_t        *cb)
{
  CS_UNUSED(f);
  CS_UNUSED(cm);
  CS_UNUSED(diff_tensor);
  CS_UNUSED(pot_values);
  CS_UNUSED(beta);
  CS_UNUSED(cb);

  cs_real_t  flux = 0.;

/*   cs_real_3_t  lek; */
/*   cs_real_3_t  *le_grd = cb->vectors; */

/*   const cs_real_t  over_vol_c = 1/cm->vol_c; */


/*   /\* Loop on border face edges *\/ */
/*   for (short int e = 0; e < fm->n_ef; e++) { */

/*     const cs_quant_t  peq_i = fm->edge[e]; */
/*     const short int  vi1 = fm->v_ids[fm->e2v_ids[2*e]]; */
/*     const short int  vi2 = fm->v_ids[fm->e2v_ids[2*e+1]]; */
/*     const short int  ei = fm->e_ids[e]; // edge id in the cell mesh */
/*     const cs_nvec3_t  dfq_i = cm->dface[ei]; */
/*     const double  dp = _dp3(peq_i.unitv, dfq_i.unitv); */
/*     const double  tmp_val = peq_i.meas * dfq_i.meas * dp; // 3*pec_vol */
/*     const double  beta_pec_vol = 3. * beta/tmp_val; */
/*     const double  coef_ei = 3. * beta/dp; */

/*     /\* Reset L_Ec(GRAD(p_j)) for each vertex of the cell *\/ */
/*     for (short int v = 0; v < cm->n_vc; v++) */
/*       le_grd[v][0] = le_grd[v][1] = le_grd[v][2] = 0; */

/*     /\* Term related to the flux reconstruction: */
/*        Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell *\/ */
/*     for (short int ek = 0; ek < cm->n_ec; ek++) { */

/*       const short int  vj1 = cm->e2v_ids[2*ek], vj2 = cm->e2v_ids[2*ek+1]; */
/*       const short int  sgn_vj1 = cm->e2v_sgn[ek]; // sgn_vj2 = - sgn_vj1 */
/*       const cs_nvec3_t  dfq_k = cm->dface[ek]; */

/*       /\* Compute l_(ek,c)|p(ei,f,c) *\/ */
/*       const double  eik_part = coef_ei * _dp3(dfq_k.unitv, peq_i.unitv); */
/*       const double  coef_mult = dfq_k.meas * over_vol_c; */

/*       for (int k = 0; k < 3; k++) */
/*         lek[k] = coef_mult * (dfq_k.unitv[k] - eik_part * dfq_i.unitv[k]); */

/*       if (ek == ei) */
/*         for (int k = 0; k < 3; k++) */
/*           lek[k] += dfq_k.meas * dfq_k.unitv[k] * beta_pec_vol; */

/*       for (int k = 0; k < 3; k++) { */
/*         le_grd[vj1][k] += sgn_vj1 * lek[k]; */
/*         le_grd[vj2][k] -= sgn_vj1 * lek[k]; */
/*       } */

/*     } // Loop on cell edges */

/*     for (short int v = 0; v < cm->n_vc; v++) { */

/*       const double  contrib = _dp3(mnu, le_grd[v]) * fm->tef[e]; */
/*       ntrgrd->val[vi1*cm->n_vc + v] += contrib; */
/*       ntrgrd->val[vi2*cm->n_vc + v] += contrib; */

/*     } */

/*   } // Loop on face edges */

/* #if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 1 */
/*   cs_log_printf(CS_LOG_DEFAULT, */
/*                 ">> Flux.Op (NTRGRD) matrix (c_id: %d,f_id: %d)", */
/*                 cm->c_id, cm->f_ids[fm->f_id]); */
/*   cs_sdm_dump(cm->c_id, NULL, NULL, ntrgrd); */
/* #endif */

  return flux;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
