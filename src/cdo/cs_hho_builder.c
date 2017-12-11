/*============================================================================
 * Low-level functions and structures used to build the algebraic system with
 * a cellwise process when Hybrid High Order schemes are set for the space
 * discretization
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

#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "cs_log.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hho_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product
#define _mv3  cs_math_33_3_product

#define CS_HHO_BUILDER_DBG  1

/*============================================================================
 * Private variables
 *============================================================================*/

static const double  cs_hho_builder_pena_coef = 1e13;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Fill the volume-related part of the matrices when building the
 *         reconstruction operator
 *
 * \param[in, out]  stiffness   pointer to the stiffness matrix
 * \param[in, out]  rhs_c_t     pointer to the right-hand side (matrix)
 * \param[in, out]  hhob        pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_fill_vol_reco_op(cs_sdm_t             *stiffness,
                  cs_sdm_t             *rhs_c_t,
                  cs_hho_builder_t     *hhob)
{
  const short int  gs = hhob->grad_basis->size - 1;
  const short int  cs = hhob->cell_basis->size;

  /* Stiffness matrix is symmetric. Only the right upper is build */
  cs_real_t  *mg = stiffness->val;
  for (short int i = 0; i < gs; i++) {
    const cs_real_t  *mg_i = mg + i*gs;
    for (short int j = i + 1; j < gs; j++)
      mg[j*gs + i] = mg_i[j];
  }

  /* rhs_c shares a part of m_g (the stiffness matrix). Add this contribution */
  cs_real_t *r_j = rhs_c_t->val + gs; // Skip the first row (only 0)
  for (short int j = 0; j < cs - 1; j++, r_j += gs) {
    const cs_real_t  *mg_j = mg + j*gs;
    for (short int i = 0; i < gs; i++)
      r_j[i] += mg_j[i];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the face polynomial space of a function
 *         defined by an analytical expression depending on the location and
 *         the current time
 *
 * \param[in]       anai     pointer to an analytical definition
 * \param[in]       fbf      pointer to a structure for face basis functions
 * \param[in]       xv1      first vertex
 * \param[in]       xv2      second vertex
 * \param[in]       xv3      third vertex
 * \param[in]       surf     area of the triangle
 * \param[in, out]  cb       pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array    array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tria_reduction(const cs_xdef_analytic_input_t  *anai,
                    const cs_basis_func_t           *fbf,
                    const cs_real_3_t                xv1,
                    const cs_real_3_t                xv2,
                    const cs_real_3_t                xv3,
                    const double                     surf,
                    cs_cell_builder_t               *cb,
                    cs_real_t                        array[])
{
  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *ana_eval = cb->values + 7;
  cs_real_t  *phi_eval = cb->values + 14;

  /* Compute Gauss points and related weights */
  cs_quadrature_tria_7pts(xv1, xv2, xv3, surf, gpts, gw);

  /* Evaluate the analytical function at the Gauss points */
  anai->func(cs_glob_time_step->t_cur, 7, NULL, (const cs_real_t *)gpts, true,
             anai->input, ana_eval);

  for (short int gp = 0; gp < 7; gp++) {

    fbf->eval_all_at_point(fbf, gpts[gp], phi_eval);

    const cs_real_t  w = gw[gp] * ana_eval[gp];
    for (short int i = 0; i < fbf->size; i++)
      array[i] += w * phi_eval[i];

  }  /* End of loop on Gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the cell polynomial space of a function
 *         defined by an analytical expression depending on the location and
 *         the current time
 *
 * \param[in]       anai     pointer to an analytical definition
 * \param[in]       cbf      pointer to a structure for face basis functions
 * \param[in]       xv1      first vertex
 * \param[in]       xv2      second vertex
 * \param[in]       xv3      third vertex
 * \param[in]       xv4      third vertex
 * \param[in]       vol      volume of the tetrahedron
 * \param[in, out]  cb       pointer to a cs_cell_builder_structure_t
 * \param[in, out]  array    array storing values to compute
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_tetra_reduction(const cs_xdef_analytic_input_t  *anai,
                     const cs_basis_func_t           *cbf,
                     const cs_real_3_t                xv1,
                     const cs_real_3_t                xv2,
                     const cs_real_3_t                xv3,
                     const cs_real_3_t                xv4,
                     const double                     vol,
                     cs_cell_builder_t               *cb,
                     cs_real_t                        array[])
{
  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *ana_eval = cb->values + 15;
  cs_real_t  *phi_eval = cb->values + 30;

  /* Compute Gauss points and related weights */
  cs_quadrature_tet_15pts(xv1, xv2, xv3, xv4, vol, gpts, gw);

  /* Evaluate the analytical function at the Gauss points */
  anai->func(cs_glob_time_step->t_cur, 15, NULL, (const cs_real_t *)gpts, true,
             anai->input, ana_eval);

  for (short int gp = 0; gp < 15; gp++) {

    cbf->eval_all_at_point(cbf, gpts[gp], phi_eval);

    const cs_real_t  w = gw[gp] * ana_eval[gp];
    for (short int i = 0; i < cbf->size; i++)
      array[i] += w * phi_eval[i];

  }  /* End of loop on Gauss points */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine for computing volumetric integrals over a tetrahedron which
 *         are contributions to the local stiffness matrix on the gradient
 *         basis and to the right-hand side
 *
 * \param[in]      xv1        first vertex
 * \param[in]      xv2        second vertex
 * \param[in]      xv3        third vertex
 * \param[in]      surf       area of the triangle
 * \param[in]      fbf        pointer to the related set of face basis functions
 * \param[in]      kappa_nfc  permeability tensor times the related face normal
 * \param[in, out] gpts       coordinates of the Gauss points
 * \param[in, out] rc         right-hand side matrix to compute (cell part)
 * \param[in, out] rf         right-hand side matrix to compute (face part)
 * \param[in, out] gpts       coordinates of the Gauss points
 * \param[in, out] kappa_nfc  coordinates of the Gauss points
 * \param[in, out] cb         pointer to a cs_cell_builder_structure_t
 * \param[in, out] hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_tria_to_reco_op(const cs_real_3_t         xv1,
                     const cs_real_3_t         xv2,
                     const cs_real_3_t         xv3,
                     const double              surf,
                     const cs_basis_func_t    *fbf,
                     const cs_real_t          *kappa_nfc,
                     cs_real_3_t              *gpts,
                     cs_sdm_t                 *rc,
                     cs_sdm_t                 *rf,
                     cs_cell_builder_t        *cb,
                     cs_hho_builder_t         *hhob)
{
  const cs_basis_func_t  *cbf = hhob->cell_basis;
  const cs_basis_func_t  *gbf = hhob->grad_basis;

  const short int  fsize = fbf->size;
  const short int  csize = cbf->size;
  const short int  gsize = gbf->size - 1;

  cs_real_t  *gw = cb->values;
  cs_real_t  *f_phi = cb->values + fbf->n_gpts_tria;
  cs_real_t  *c_phi = cb->values + fbf->n_gpts_tria + fsize;
  cs_real_t  *g_phi = cb->values + fbf->n_gpts_tria + fsize + csize;

  /* Compute Gauss points and related weights */
  fbf->quadrature_tria(xv1, xv2, xv3, surf, gpts, gw);

  /* Build the stiffness matrix M_g by adding the contribution of the
     current tetrahedron.
     M_g is symmetric (fill the right upper part first) */
  for (short int gp = 0; gp < fbf->n_gpts_tria; gp++) {

    gbf->eval_all_at_point(gbf, gpts[gp], g_phi);
    cbf->eval_all_at_point(cbf, gpts[gp], c_phi);
    fbf->eval_all_at_point(fbf, gpts[gp], f_phi);

    for (short int i = 0; i < gsize; i++) {

      const cs_real_t  coef = gw[gp] * _dp3(kappa_nfc, g_phi + 3*(i+1));

      /* cell part */
      for (short int j = 0; j < csize; j++)
        rc->val[j*gsize + i] -= coef * c_phi[j];

      /* face part */
      for (short int j = 0; j < fsize; j++)
        rf->val[j*gsize + i] += coef * f_phi[j];

    } /* End of loop on matrix columns */

  }  /* End of loop on Gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine for computing volumetric integrals over a tetrahedron which
 *         are contributions to the local stiffness matrix on the gradient
 *         basis and to the right-hand side
 *
 * \param[in]      xv1         first vertex
 * \param[in]      xv2         second vertex
 * \param[in]      xv3         third vertex
 * \param[in]      xv4         fourth vertex
 * \param[in ]     vol         volume of the terahedron
 * \param[in, out] stiffness   stiffness matrix to compute
 * \param[in, out] gpts        coordinates of the Gauss points
 * \param[in, out] cb          pointer to a cs_cell_builder_structure_t
 * \param[in, out] hhob        pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_tetra_to_reco_op(const cs_real_3_t         xv1,
                      const cs_real_3_t         xv2,
                      const cs_real_3_t         xv3,
                      const cs_real_3_t         xv4,
                      const double              vol,
                      cs_sdm_t                 *stiffness,
                      cs_real_3_t              *gpts,
                      cs_cell_builder_t        *cb,
                      cs_hho_builder_t         *hhob)
{
  const cs_basis_func_t  *gbf = hhob->grad_basis;
  const short int gs = gbf->size - 1;

  cs_real_3_t Kgrad_phi_i;
  cs_real_t  *gw = cb->values;
  cs_real_t  *g_phi = cb->values + gbf->n_gpts_tetra;

  /* Compute Gauss points and related weights */
  gbf->quadrature_tetra(xv1, xv2, xv3, xv4, vol, gpts, gw);

  /* Build the stiffness matrix M_g by adding the contribution of the
     current tetrahedron.
     M_g is symmetric (fill the right upper part first) */
  for (short int gp = 0; gp < gbf->n_gpts_tetra; gp++) {

    const cs_real_t  w = gw[gp];

    gbf->eval_all_at_point(gbf, gpts[gp], g_phi);

    for (short int i = 0; i < gs; i++) {

      _mv3((const cs_real_t (*)[3])cb->pty_mat, g_phi + 3*(i+1), Kgrad_phi_i);

      cs_real_t  *mg_i = stiffness->val + i*gs;
      for (short int j = i; j < gs; j++)
        mg_i[j] += w * _dp3(Kgrad_phi_i, g_phi + 3*(j+1));

    } /* End of loop on matrix rows */

  }  /* End of loop on Gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine for computing volumetric integrals over a tetrahedron which
 *         are contributions to the M_cg matrix
 *
 * \param[in]      xv1        first vertex
 * \param[in]      xv2        second vertex
 * \param[in]      xv3        third vertex
 * \param[in]      xv4        fourth vertex
 * \param[in]      vol        volume of the tetrahedron
 * \param[in]      cbf_kp1    pointer to a cell basis function order:k+1
 * \param[in, out] m_cg       matrix to compute
 * \param[in, out] cb         pointer to a cs_cell_builder_structure_t
 * \param[in, out] hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_contrib_mcg(const cs_real_3_t         xv1,
                 const cs_real_3_t         xv2,
                 const cs_real_3_t         xv3,
                 const cs_real_3_t         xv4,
                 const double              vol,
                 const cs_basis_func_t    *cbf_kp1,
                 cs_sdm_t                 *m_cg,
                 cs_cell_builder_t        *cb,
                 cs_hho_builder_t         *hhob)
{
  const int n_gpts = cbf_kp1->n_gpts_tetra;
  const short int  cs_kp1 = cbf_kp1->size;
  const short int  cs = hhob->cell_basis->size;
  const short int  gs = hhob->grad_basis->size - 1;

  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *c_phi = cb->values + n_gpts;

  /* Compute Gauss points and related weights */
  cbf_kp1->quadrature_tetra(xv1, xv2, xv3, xv4, vol, gpts, gw);

  for (short int gp = 0; gp < cbf_kp1->n_gpts_tetra; gp++) {

    cbf_kp1->eval_all_at_point(cbf_kp1, gpts[gp], c_phi);

    const cs_real_t  w = gw[gp];
    for (short int i = 0; i < cs; i++) {

      const cs_real_t  coef_i = w * c_phi[i];
      cs_real_t  *m_cg_i = m_cg->val + i*gs;
      for (short int j = cs; j < cs_kp1; j++)
        m_cg_i[j-1] += coef_i * c_phi[j];

    } /* End of loop on matrix rows */

  }  /* End of loop on Gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Routine for computing surfacic integrals over a triangle which
 *         are contributions to the Mf_cg matrix
 *
 * \param[in]      xv1        first vertex
 * \param[in]      xv2        second vertex
 * \param[in]      xv3        third vertex
 * \param[in]      surf       surface of the triangle to consider
 * \param[in]      fbf        pointer to a set of face basis functions
 * \param[in]      cbf_kp1    pointer to a set of cell basis functions (k+1)
 * \param[in, out] mf_cg      matrix to compute
 * \param[in, out] cb         pointer to a cs_cell_builder_structure_t
 * \param[in, out] hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_add_contrib_mf_cg(const cs_real_3_t         xv1,
                   const cs_real_3_t         xv2,
                   const cs_real_3_t         xv3,
                   const double              surf,
                   const cs_basis_func_t    *fbf,
                   const cs_basis_func_t    *cbf_kp1,
                   cs_sdm_t                 *mf_cg,
                   cs_cell_builder_t        *cb,
                   cs_hho_builder_t         *hhob)
{
  const int n_gpts = fbf->n_gpts_tria;
  const short int  cs = hhob->cell_basis->size;
  const short int  fs = fbf->size;
  const short int  cs_kp1 = cbf_kp1->size;

  cs_real_3_t  *gpts = cb->vectors;
  cs_real_t  *gw = cb->values;
  cs_real_t  *f_phi = cb->values + n_gpts;
  cs_real_t  *c_phi = f_phi + fbf->size;

  /* Compute Gauss points and related weights */
  fbf->quadrature_tria(xv1, xv2, xv3, surf, gpts, gw);

  for (short int gp = 0; gp < fbf->n_gpts_tria; gp++) {

    cbf_kp1->eval_all_at_point(cbf_kp1, gpts[gp], c_phi);
    fbf->eval_all_at_point(fbf, gpts[gp], f_phi);

    const cs_real_t  w = gw[gp];

    /* First block (column) c_phi0 */
    cs_sdm_t  *m0 = cs_sdm_get_block(mf_cg, 0, 0);
    const cs_real_t  coef0 = w * c_phi[0];
    for (short int j = 0; j < fs; j++)
      m0->val[j] += coef0 * f_phi[j];

    /* Second block (c_phi_k \ 0) */
    cs_sdm_t  *mk = cs_sdm_get_block(mf_cg, 0, 1);

    /* Third block (c_phi_k+1 \ k) */
    cs_sdm_t  *mkp1 = cs_sdm_get_block(mf_cg, 0, 2);

    for (short int i = 0; i < fs; i++) {
      const cs_real_t  coef_i = w * f_phi[i];

      cs_real_t  *mk_i = mk->val + i*(cs-1);
      for (short int j = 1; j < cs; j++)
        mk_i[j-1] += coef_i * c_phi[j];

      cs_real_t  *mkp1_i = mkp1->val + i*(cs_kp1-cs);
      short int shift = 0;
      for (short int j = cs; j < cs_kp1; j++, shift++)
        mkp1_i[shift] += coef_i * c_phi[j];

    } /* Loop on rows */

  }  /* End of loop on Gauss points */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusion operator. The gradient reconstruction operator
 *         has to be built just before this call (cb->aux stores the rhs)
 *
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in]       cbf_kp1  pointer to a set of cell basis functions order=k+1
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 *
 * \return a pointer to a cs_sdm_t structure storing m_cg
 */
/*----------------------------------------------------------------------------*/

static cs_sdm_t *
_compute_mcg(const cs_cell_mesh_t    *cm,
             cs_basis_func_t         *cbf_kp1,
             cs_cell_builder_t       *cb,
             cs_hho_builder_t        *hhob)
{
  const cs_basis_func_t  *cbf = hhob->cell_basis;
  const short int  cs = cbf->size;
  const short int  cs_kp1 = cbf_kp1->size;
  const short int  gs = hhob->grad_basis->size - 1;

  assert(gs + 1 == cs_kp1);

  cs_sdm_t  *mcg = cb->hdg;

  /* Set values to zero */
  cs_sdm_init(cs, gs, mcg);

  /* First copy the square block cs x cs from Mcc (projector)
     Skip the first column */
  for (short int i = 0; i < cs; i++) {
    const cs_real_t  *pcc_i = cbf->projector->val + i*cs;
    cs_real_t  *mcg_i = mcg->val + i*gs;
    for (short int j = 1; j < cs; j++)
      mcg_i[j-1] = pcc_i[j];
  }

  /* Switching on cell-type: optimised version for tetra */
  switch (cm->type) {

  case FVM_CELL_TETRA:
    _add_contrib_mcg(cm->xv, cm->xv + 3, cm->xv + 6, cm->xv + 9,
                     cm->vol_c, cbf_kp1, mcg, cb, hhob);
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_onethird * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;

      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          _add_contrib_mcg(cm->xv + 3*v0, cm->xv + 3*v1, cm->xv + 3*v2, cm->xc,
                           hf_coef*pfq.meas, cbf_kp1, mcg, cb, hhob);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _add_contrib_mcg(xv0, xv1, pfq.center, cm->xc, hf_coef*tef[e],
                             cbf_kp1, mcg, cb, hhob);
          }
        }
        break;

      } /* End of switch */

    }   /* End of loop on faces */

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */

  return mcg;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate  a cs_hho_builder_t structure
 *
 * \param[in] order             order of the polynomial basis function
 * \param[in] n_fc              max. number of faces in a cell
 *
 * \return a pointer to a new allocated cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hho_builder_t *
cs_hho_builder_create(int     order,
                      int     n_fc)
{
  cs_hho_builder_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_hho_builder_t);

  /* Retrieve options for building basis functions */
  cs_flag_t  face_basis_flag, cell_basis_flag;
  cs_basis_func_get_hho_flag(&face_basis_flag, &cell_basis_flag);

  b->n_face_basis = 0;
  b->n_max_face_basis = n_fc;
  BFT_MALLOC(b->face_basis, n_fc, cs_basis_func_t *);
  for (int i = 0; i < n_fc; i++)
    b->face_basis[i] = cs_basis_func_create(face_basis_flag, order, 2);

  b->cell_basis = cs_basis_func_create(cell_basis_flag, order, 3);
  b->grad_basis = cs_basis_func_grad_create(b->cell_basis);

  const short int fbs = b->face_basis[0]->size;
  const short int cbs = b->cell_basis->size;
  const short int gbs = b->grad_basis->size - 1;

  /* Reconstruction operator for the gradient of the potential.
     Store the transposed matrix for a more convenient manipulation */
  short int  *block_size = NULL;
  BFT_MALLOC(block_size, n_fc + 1, short int);
  for (short int f = 0; f < n_fc; f++) {
    block_size[f] = fbs;
  }
  block_size[n_fc] = cbs;
  assert(gbs >= cbs);

  b->grad_reco_op = cs_sdm_block_create(n_fc + 1, 1, block_size, &gbs);
  b->tmp = cs_sdm_block_create(1 + n_fc, 1, block_size, &fbs);
  b->bf_t = cs_sdm_block_create(1 + n_fc, 1, block_size, &fbs);
  b->jstab = cs_sdm_block_create(1 + n_fc, 1 + n_fc, block_size, block_size);

  BFT_FREE(block_size);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_hho_builder_t structure
 *
 * \param[in, out] p_builder  pointer of pointer on a cs_hho_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_free(cs_hho_builder_t  **p_builder)
{
  if (p_builder == NULL)
    return;

  cs_hho_builder_t  *b = *p_builder;

  /* Free all basis */
  b->grad_basis = cs_basis_func_free(b->grad_basis);
  b->cell_basis = cs_basis_func_free(b->cell_basis);
  for (int i = 0; i < b->n_max_face_basis; i++)
    b->face_basis[i] = cs_basis_func_free(b->face_basis[i]);
  BFT_FREE(b->face_basis);

  /* Free matrices */
  b->grad_reco_op = cs_sdm_free(b->grad_reco_op);
  b->tmp = cs_sdm_free(b->tmp);
  b->bf_t = cs_sdm_free(b->bf_t);
  b->jstab = cs_sdm_free(b->jstab);

  BFT_FREE(b);

  *p_builder = NULL;
  p_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set-up the basis functions related to a cell, its gradient and to
 *         the faces of this cell.
 *         Compute cell and face projection and its related modified Cholesky
 *         factorization.
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_cellwise_setup(const cs_cell_mesh_t    *cm,
                              cs_cell_builder_t       *cb,
                              cs_hho_builder_t        *hhob)
{
  if (hhob == NULL)
    return;

  /* Sanity checks */
  assert(cm != NULL);

  hhob->n_face_basis = cm->n_fc;

  /* Setup cell basis functions */
  hhob->cell_basis->setup(hhob->cell_basis, cm, 0, cm->xc, cb);

  /* Compute M_{cc} and its modified Cholesky decomposition */
  hhob->cell_basis->compute_projector(hhob->cell_basis, cm, 0);
  hhob->cell_basis->compute_factorization(hhob->cell_basis);

  /* Setup gradient basis functions (for cell) */
  cs_basis_func_copy_setup(hhob->cell_basis, hhob->grad_basis);

  /* Setup face basis functions */
  for (short int f = 0; f < cm->n_fc; f++) {

    hhob->face_basis[f]->setup(hhob->face_basis[f],
                               cm, f, cm->face[f].center, cb);

    /* Compute M_{ff} and its modified Cholesky decomposition */
    hhob->face_basis[f]->compute_projector(hhob->face_basis[f], cm, f);
    hhob->face_basis[f]->compute_factorization(hhob->face_basis[f]);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the gradient operator stemming from the relation
 *            stiffness * grad_op = rhs
 *         where stiffness is a square matrix of size grd_size
 *               rhs is matrix of size (n_fc*f_size + c_size) * grd_size
 *         Hence, grad_op a matrix grd_size * (n_fc*f_size + c_size)
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_compute_grad_reco(const cs_cell_mesh_t    *cm,
                                 cs_cell_builder_t       *cb,
                                 cs_hho_builder_t        *hhob)
{
  if (hhob == NULL)
    return;
  assert(hhob->cell_basis != NULL && hhob->grad_basis != NULL &&
         hhob->face_basis != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  const short int  gs = hhob->grad_basis->size - 1;
  cs_sdm_t  *stiffness = cb->hdg, *rhs_t = cb->aux;
  cs_sdm_square_init(gs, stiffness);

  /* Initialize the matrix related to the gradient reconstruction operator.
     The right-hand side is a matrix since we are building an operator
     We define its transposed version to have a more convenient access to its
     columns */
  short int  tots = 0;
  for (short int f = 0; f < cm->n_fc; f++) {
    cb->ids[f] = hhob->face_basis[f]->size;
    tots += hhob->face_basis[f]->size;
  }
  tots += hhob->cell_basis->size;
  cb->ids[cm->n_fc] = hhob->cell_basis->size;

  cs_sdm_block_init(rhs_t, cm->n_fc + 1, 1, cb->ids, &gs);
  cs_sdm_block_init(hhob->grad_reco_op, cm->n_fc + 1, 1, cb->ids, &gs);

  /* Retrieve the rhs block related to "cell" and to this face.
     rhs is a matrix. We store its transposed. */
  cs_sdm_t  *rc = cs_sdm_get_block(rhs_t, cm->n_fc, 0);

  /* Pre-compute tensor x outward normal */
  cs_real_3_t  *kappa_nfc = cb->vectors;
  cs_real_3_t  *gpts = cb->vectors + cm->n_fc;
  for (short int f = 0; f < cm->n_fc; f++) {
    _mv3((const cs_real_t (*)[3])cb->pty_mat, cm->face[f].unitv, kappa_nfc[f]);
    kappa_nfc[f][0] *= cm->f_sgn[f];
    kappa_nfc[f][1] *= cm->f_sgn[f];
    kappa_nfc[f][2] *= cm->f_sgn[f];
  }

  /* Switching on cell-type: optimised version for tetra */
  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);
      _add_tetra_to_reco_op(cm->xv, cm->xv+3, cm->xv+6, cm->xv+9, cm->vol_c,
                            stiffness, gpts, cb, hhob);
      _fill_vol_reco_op(stiffness, rc, hhob);

      short int  v0, v1, v2;
      for (short int f = 0; f < cm->n_fc; ++f) {

        const cs_basis_func_t  *fbf = hhob->face_basis[f];
        const cs_real_t  *knfc = (const cs_real_t *)(kappa_nfc[f]);
        const cs_quant_t  pfq = cm->face[f];
        const short int  *f2e_ids = cm->f2e_ids + cm->f2e_idx[f];

        /* Retrieve the rhs block related to this face.
           rhs is a matrix. We store its transposed. */
        cs_sdm_t  *rf = cs_sdm_get_block(rhs_t, f, 0);

        cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

        _add_tria_to_reco_op(cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
                             fbf, knfc, gpts,
                             rc, rf, cb, hhob);
      }
    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_basis_func_t  *fbf = hhob->face_basis[f];
      const cs_real_t  *knfc = (const cs_real_t *)(kappa_nfc[f]);
      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_onethird * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;

      /* Retrieve the rhs block related to this face.
         rhs is a matrix. We store its transposed. */
      cs_sdm_t  *rf = cs_sdm_get_block(rhs_t, f, 0);

      assert(n_vf > 2);
      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          const double  *xv0 = cm->xv + 3*v0;
          const double  *xv1 = cm->xv + 3*v1;
          const double  *xv2 = cm->xv + 3*v2;
          _add_tetra_to_reco_op(xv0, xv1, xv2, cm->xc, hf_coef * pfq.meas,
                                stiffness, gpts, cb, hhob);

          _add_tria_to_reco_op(xv0, xv1, xv2, pfq.meas,
                               fbf, knfc, gpts,
                               rc, rf, cb, hhob);

        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _add_tetra_to_reco_op(xv0, xv1, pfq.center, cm->xc, hf_coef*tef[e],
                                  stiffness, gpts, cb, hhob);

            _add_tria_to_reco_op(xv0, xv1, pfq.center, tef[e],
                                 fbf, knfc, gpts,
                                 rc, rf, cb, hhob);

          }
        }
        break;

      } /* End of switch */

    }   /* End of loop on faces */

    _fill_vol_reco_op(stiffness, rc, hhob);

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */

  /* Compute the modified Cholesky factorization of the stiffness matrix */
  cs_sdm_t  *rhs = NULL, *gop = NULL;
  cs_real_t  *facto = cb->values;

  if (hhob->face_basis[0]->poly_order == 0) {

    cs_sdm_33_ldlt_compute(stiffness, facto);
    assert(gs == 3);

    /* Retrieve the rhs block related to "cell" and to this face.
       rhs is a matrix. We store its transposed. */
    for (short int f = 0; f < cm->n_fc; ++f) {
      rhs = cs_sdm_get_block(rhs_t, f, 0);
      gop = cs_sdm_get_block(hhob->grad_reco_op, f, 0);
      assert(rhs->n_cols == 3 && rhs->n_rows == 1);
      assert(gop->n_cols == 3 && gop->n_rows == 1);
      cs_sdm_33_ldlt_solve(facto, rhs->val, gop->val);
    }

    /* Cell contributions */
    rhs = cs_sdm_get_block(rhs_t, cm->n_fc, 0);
    gop = cs_sdm_get_block(hhob->grad_reco_op, cm->n_fc, 0);
    cs_sdm_33_ldlt_solve(facto, rhs->val, gop->val);

  }
  else {

    short int  fs, cs, facto_shift;
    if (hhob->face_basis[0]->poly_order == 1) {
      fs = 3;
      cs = 4;
      facto_shift = 45;
      assert(gs == 9);
    }
    else if (hhob->face_basis[0]->poly_order == 2) {
      fs = 6;
      cs = 10;
      facto_shift = 190;
      assert(gs == 19);
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "Polynomial order is limited to 2 up to now.");

    cs_sdm_ldlt_compute(stiffness, facto, facto + facto_shift);

    /* Retrieve the rhs block related to "cell" and to this face.
       rhs is a matrix. We store its transposed. */
    for (short int f = 0; f < cm->n_fc; ++f) {
      rhs = cs_sdm_get_block(rhs_t, f, 0);
      gop = cs_sdm_get_block(hhob->grad_reco_op, f, 0);
      assert(rhs->n_cols == gs && rhs->n_rows == fs);
      assert(gop->n_cols == gs && gop->n_rows == fs);
      for (int i = 0; i < fs; i++)
        cs_sdm_ldlt_solve(gs, facto, rhs->val + i*gs, gop->val + i*gs);
    }

    /* Cell contributions */
    rhs = cs_sdm_get_block(rhs_t, cm->n_fc, 0);
    gop = cs_sdm_get_block(hhob->grad_reco_op, cm->n_fc, 0);
    assert(rhs->n_cols == gs && rhs->n_rows == cs);
    assert(gop->n_cols == gs && gop->n_rows == cs);
    for (int i = 0; i < cs; i++)
      cs_sdm_ldlt_solve(gs, facto, rhs->val + i*gs, gop->val + i*gs);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusion operator. The gradient reconstruction operator
 *         has to be built just before this call (cb->aux stores the rhs)
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       cb         pointer to a cell builder_t structure
 * \param[in, out]  hhob       pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_diffusion(const cs_cell_mesh_t    *cm,
                         cs_cell_builder_t       *cb,
                         cs_hho_builder_t        *hhob)
{
  if (hhob == NULL)
    return;
  assert(hhob->cell_basis != NULL && hhob->grad_basis != NULL &&
         hhob->face_basis != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  const cs_basis_func_t  *cbf = hhob->cell_basis;
  const short int  cs = cbf->size;
  const short int  gs = hhob->grad_basis->size - 1;
  const short int  fs = hhob->face_basis[0]->size;

  cs_sdm_t  *rhs_t = cb->aux;
  cs_sdm_t  *gop_t = hhob->grad_reco_op;

  /* Initialize the matrix related to the diffusion term */
  short int  tots = 0;
  for (short int f = 0; f < cm->n_fc; f++) {
    assert(fs == hhob->face_basis[f]->size);
    cb->ids[f] = fs;
    tots += fs;
  }
  cb->ids[cm->n_fc] = cs;
  tots += cs;

  cs_sdm_block_init(cb->loc, cm->n_fc + 1, cm->n_fc + 1, cb->ids, cb->ids);

  /* Compute the consistent part */
  cs_sdm_block_multiply_rowrow_sym(rhs_t, gop_t, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_BUILDER_DBG > 0
  if (cm->c_id == 0) {
    printf(" Diffusion matrix: rhs_t\n");
    cs_sdm_block_fprintf(NULL, NULL, 1e-16, rhs_t);

    printf(" Diffusion matrix: grad_op_t\n");
    cs_sdm_block_fprintf(NULL, NULL, 1e-16, gop_t);

    printf(" Diffusion matrix: consistent part\n");
    cs_sdm_block_fprintf(NULL, NULL, 1e-16, cb->loc);
  }
#endif

  /* Compute the stabilization part.
     Several steps are needed :
     1. Define a set of basis functions for P^(k+1)_c
     2. Compute additional matrices related to cell and gradient functions
     3. Compose/Multiply cell-related matrices
     4. Define for each face a contribution to the stabilization term
  */

  cs_basis_func_t  *cbf_kp1 = cs_basis_func_create(cbf->flag,
                                                   cbf->poly_order + 1,
                                                   cbf->dim);
  cs_basis_func_copy_setup(cbf, cbf_kp1);

  const short int  cs_kp1 = cbf_kp1->size;

  /* M_cg is stored in cb->hdg (temporary matrix) */
  cs_sdm_t  *m_cg = _compute_mcg(cm, cbf_kp1, cb, hhob);

  /* First part of the stabilization. This part only relies on cell quantities.

     Compute Id_c - M_cc^-1 * M_cg * GradOp (transposed of GradOp is stored)
     The resulting matrix is denoted by M_ccgg. Its transposed matrix is
     computed.
     This is done in one pass thanks to the block structure.

     The projection onto the cell basis function is performed thanks to a
     modified Cholesky decomposition.
   */

  cs_sdm_t  *m_ccgg = cb->aux;  /* Constant for this cell */
  cs_sdm_block_init(m_ccgg, cm->n_fc + 1, 1, cb->ids, &cs);

  cs_real_t  *array = NULL;
  cs_real_t  _array[10];
  if (cs > 10)
    BFT_MALLOC(array, cs, cs_real_t);
  else
    array = _array;

  assert(hhob->grad_reco_op->block_desc->n_row_blocks == cm->n_fc + 1);

  /* Blocks related to faces (b < n_fc) or cell (b == n_fc) */
  for (short int b = 0; b < cm->n_fc + 1; b++) {

    const cs_sdm_t  *gb = cs_sdm_get_block(hhob->grad_reco_op, b, 0);

    /* Block mb = fs (or cs) * cs to compute */
    cs_sdm_t  *mb = cs_sdm_get_block(m_ccgg, b, 0);

    assert(mb->n_rows == gb->n_rows); /* sanity check */

    for (short int j = 0; j < mb->n_rows; j++) { /* Works column by column */

      /* Extract column array (easy since we have stored the transposed) */
      const cs_real_t  *const gb_j = gb->val + j*gs;

      /* Compute the "j"th column array of M_cg * GradOp  */
      for (short int i = 0; i < cs; i++) {

        /* Equivalent to a part of the matrix-matrix computation */
        const cs_real_t  *const m_cg_i = m_cg->val + i*gs;
        array[i] = 0;
        for (short int k = 0; k < gs; k++)
          array[i] += gb_j[k] * m_cg_i[k];

      }

      /* Define the row "j": M_cc^-1 then * -1
         The transposed of M_ccgg is stored */
      cs_real_t  *mb_j = mb->val + j*cs;
      cbf->project(cbf, array, mb_j);
      for (short int i = 0; i < cs; i++)
        mb_j[i] *= -1;

      if (b == cm->n_fc) /* cell block */
        mb_j[j] += 1;

    } /* Build the row j of the matrix in block b */

  } /* Loop on blocks (face or cell) */

  /* Add the second contribution for the stabilization part which is related
     to each face
     This part results from M_ff^-1 *{M_fg*GradOp + M_fc*M_ccgg} - Id_f
     where transposed of GradOp and M_ccgg are stored.
  */

  cs_sdm_t  *bf_t_mff = hhob->tmp;

  assert(m_ccgg->block_desc->n_row_blocks == cm->n_fc + 1);

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_basis_func_t  *fbf = hhob->face_basis[f];

    /* --1-- Build M_fc and M_fg (stored in cb->hdg) */
    cs_sdm_t  *mf_cg = cb->hdg;
    const short int  m_sizes[3] = {1, cs-1, cs_kp1 - cs};
    cs_sdm_block_init(mf_cg, 1, 3, &fs, m_sizes);

    const cs_quant_t  pfq = cm->face[f];
    const int  start = cm->f2e_idx[f];
    const int  end = cm->f2e_idx[f+1];
    const short int n_vf = end - start; // #vertices (=#edges)
    const short int *f2e_ids = cm->f2e_ids + start;

    switch(n_vf){

    case 3: /* triangle (optimized version, no subdivision) */
      {
        short int  v0, v1, v2;
        cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

        _add_contrib_mf_cg(cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
                           fbf, cbf_kp1, mf_cg, cb, hhob);
      }
      break;

    default:
      {
        const double  *tef = cm->tef + start;

        for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

          // Edge-related variables
          const short int e0  = f2e_ids[e];
          const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
          const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

          _add_contrib_mf_cg(xv0, xv1, pfq.center, tef[e],
                             fbf, cbf_kp1, mf_cg, cb, hhob);
        }
      }
      break;

    } /* End of switch */

    /* Blocks of Mf_cg: c_phi0, (c_phi_k\0) and (c_phi_k+1\k) respectively */
    const cs_sdm_t  *const m0 = cs_sdm_get_block(mf_cg, 0, 0);
    const cs_sdm_t  *const mk = cs_sdm_get_block(mf_cg, 0, 1);
    const cs_sdm_t  *const mkp1 = cs_sdm_get_block(mf_cg, 0, 2);

    /* Reset transposed of Bf and temporary matrix */
    cs_sdm_block_init(hhob->bf_t, cm->n_fc + 1, 1, cb->ids, &fs);
    cs_sdm_block_init(bf_t_mff,  cm->n_fc + 1, 1, cb->ids, &fs);

    /* --2-- Compute blocks related to faces (b < n_fc) or cell (b == n_fc) */
    for (short int b = 0; b < cm->n_fc + 1; b++) {

      const cs_sdm_t  *mb = cs_sdm_get_block(m_ccgg, b, 0);
      const cs_sdm_t  *gb = cs_sdm_get_block(hhob->grad_reco_op, b, 0);

      /* Block bf_f = fs (or cs) * fs to compute */
      cs_sdm_t  *bfb = cs_sdm_get_block(hhob->bf_t, b, 0);

      assert(bfb->n_rows == gb->n_rows); /* sanity check */
      assert(mb->n_rows == gb->n_rows); /* sanity check */

      for (short int j = 0; j < mb->n_rows; j++) { /* Works column by column */

        /* Extract column array (easy since we have stored the transposed) */
        const cs_real_t  *const mb_j = mb->val + j*cs;
        const cs_real_t  *const gb_j1 = gb->val + j*gs;
        const cs_real_t  *const gb_j2 = gb_j1 + m_sizes[1];

        /* Compute the "j"th column array of M_fg * GradOp
           Compute the "j"th column array of M_fc * M_ccgg */
        for (short int i = 0; i < fs; i++) {

          /* Optimized way to perform a matrix-matrix computation */
          array[i] = m0->val[i] * mb_j[0];
          const cs_real_t  *const mk_i = mk->val + i*m_sizes[1];
          for (short int k = 0; k < m_sizes[1]; k++)
            array[i] += (gb_j1[k] + mb_j[k+1]) * mk_i[k];
          const cs_real_t  *const mkp1_i = mkp1->val + i*m_sizes[2];
          for (short int k = 0; k < m_sizes[2]; k++)
            array[i] += gb_j2[k] * mkp1_i[k];

        }

        /* Define the row "j": M_ff^-1. The transposed of Bf is stored */
        cs_real_t  *bfb_j = bfb->val + j*fs;
        fbf->project(fbf, array, bfb_j);
        if (b == f) /* current face block */
          bfb_j[j] -= 1;

      } /* Build the row j of the matrix in block b */

      /* --3-- The current block of Bf_t is now defined. Perform Bf_t * Mff.
               Mff is symmetric by construction --> use the rowrow version
      */
      cs_sdm_t  *tmp = cs_sdm_get_block(bf_t_mff, b, 0);

      cs_sdm_multiply_rowrow(bfb, fbf->projector, tmp);

    } /* Loop on blocks (face or cell) for the current face */

    /* Store the cumulated contribution over faces in jstab */
    cs_sdm_block_init(hhob->jstab, cm->n_fc+1, cm->n_fc+1, cb->ids, cb->ids);

    /* Update the upper right part of the matrix storing the stabilization
       The lower left part is set by symmetry */
    assert(bf_t_mff->block_desc->n_col_blocks == 1);
    for (short int bi = 0; bi < cm->n_fc + 1; bi++) {
      const cs_sdm_t  *bmff_i0 = cs_sdm_get_block(bf_t_mff, bi, 0);

      for (short int bj = bi; bj < cm->n_fc + 1; bj++) {
        const cs_sdm_t  *bf_j0 = cs_sdm_get_block(hhob->bf_t, bj, 0);

        cs_sdm_t  *js_ij = cs_sdm_get_block(hhob->jstab, bi, bj);
        cs_sdm_multiply_rowrow(bmff_i0, bf_j0, js_ij);

        if (bj > bi) {  /* The stabilization part is symmetric */
          cs_sdm_t  *js_ji = cs_sdm_get_block(hhob->jstab, bj, bi);
          cs_sdm_transpose_and_update(js_ij, js_ji);
        }

      } /* Loop on col blocks */
    } /* Loop on row blocks */

    /* Determine the stabilization coefficient for this face */
    cs_real_3_t  k_nf;
    _mv3((const cs_real_t (*)[3])cb->pty_mat, cm->face[f].unitv, k_nf);
    const cs_real_t  f_coef = _dp3(k_nf, cm->face[f].unitv) / cm->f_diam[f];

    cs_sdm_block_add_mult(cb->loc, f_coef, hhob->jstab);

  } /* Loop on cell faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_BUILDER_DBG > 0
  if (cm->c_id == 0) {
    printf(" Diffusion matrix: full version\n");
    cs_sdm_block_fprintf(NULL, NULL, 1e-16, cb->loc);
  }
#endif

  /* Free temporary buffers and structures */
  cbf_kp1 = cs_basis_func_free(cbf_kp1);
  if (array != _array)
    BFT_FREE(array);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the reduction onto the polynomial spaces (cell and faces)
 *         of a function defined by an analytical expression depending on the
 *         location and the current time
 *         red array has to be allocated before calling this function
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  red      vector containing the reduction
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_reduction_from_analytic(const cs_xdef_t         *def,
                                       const cs_cell_mesh_t    *cm,
                                       cs_cell_builder_t       *cb,
                                       cs_hho_builder_t        *hhob,
                                       cs_real_t                red[])
{
  /* Sanity checks */
  if (hhob == NULL || def == NULL)
    return;
  assert(hhob->cell_basis != NULL && hhob->face_basis != NULL);
  if (red == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s : array storing the reduction has to be allocated.\n",
              __func__);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

  const cs_basis_func_t  *cbf = hhob->cell_basis;
  assert(cbf->facto != NULL && cbf->project != NULL);

  int  shift = 0;
  cs_real_t  *c_rhs = cb->values + 30 + cbf->size;
  cs_real_t  *f_rhs = c_rhs + cbf->size;

  /* Initialize cell related array */
  memset(c_rhs, 0, cbf->size*sizeof(cs_real_t));

  /* Switch according to the cell type: optimised version for tetra */
  switch (cm->type) {

  case FVM_CELL_TETRA:
    {
      assert(cm->n_fc == 4 && cm->n_vc == 4);

      _add_tetra_reduction(anai, cbf,
                           cm->xv, cm->xv+3, cm->xv+6, cm->xv+9, cm->vol_c,
                           cb, c_rhs);

      short int  v0, v1, v2;
      for (short int f = 0; f < cm->n_fc; ++f) {

        const cs_quant_t  pfq = cm->face[f];
        const short int  *f2e_ids = cm->f2e_ids + cm->f2e_idx[f];
        const cs_basis_func_t  *fbf = hhob->face_basis[f];

        assert(fbf->facto != NULL);
        assert(cbf->size >= fbf->size);

        memset(f_rhs, 0, fbf->size*sizeof(cs_real_t));

        cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

        _add_tria_reduction(anai, fbf,
                            cm->xv+3*v0, cm->xv+3*v1, cm->xv+3*v2, pfq.meas,
                            cb, f_rhs);

        /* Modified Cholesky decomposition to compute DoF */
        fbf->project(fbf, f_rhs, red + shift);
        shift += fbf->size;

      } /* Loop on cell faces */

    }
    break;

  case FVM_CELL_PYRAM:
  case FVM_CELL_PRISM:
  case FVM_CELL_HEXA:
  case FVM_CELL_POLY:
  {
    for (short int f = 0; f < cm->n_fc; ++f) {

      const cs_quant_t  pfq = cm->face[f];
      const double  hf_coef = cs_math_onethird * cm->hfc[f];
      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;
      const cs_basis_func_t  *fbf = hhob->face_basis[f];

      memset(f_rhs, 0, fbf->size*sizeof(cs_real_t));

      assert(fbf->facto != NULL);
      assert(cbf->size >= fbf->size);
      assert(n_vf > 2);
      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          const double  *xv0 = cm->xv + 3*v0;
          const double  *xv1 = cm->xv + 3*v1;
          const double  *xv2 = cm->xv + 3*v2;

          _add_tria_reduction(anai, fbf, xv0, xv1, xv2, pfq.meas, cb, f_rhs);

          _add_tetra_reduction(anai, cbf,
                               xv0, xv1, xv2, cm->xc, hf_coef * pfq.meas,
                               cb, c_rhs);

        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _add_tetra_reduction(anai, cbf, xv0, xv1, pfq.center, cm->xc,
                                 hf_coef*tef[e], cb, c_rhs);

            _add_tria_reduction(anai, fbf, xv0, xv1, pfq.center, tef[e],
                                cb, f_rhs);

          }
        }
        break;

      } /* End of switch */

      /* Modified Cholesky decomposition to compute DoF */
      fbf->project(fbf, f_rhs, red + shift);
      shift += fbf->size;

    } /* End of loop on faces */

  }
  break;

  default:
    bft_error(__FILE__, __LINE__, 0,  _(" Unknown cell-type.\n"));
    break;

  } /* End of switch on the cell-type */

  /* Modified Cholesky decomposition to compute DoF */
  cbf->project(cbf, c_rhs, red + shift);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the projection of the Dirichlet boundary conditions onto
 *         the polynomial spaces on faces
 *
 * \param[in]       def      pointer to a cs_xdef_t structure
 * \param[in]       f        local face id in the cellwise view of the mesh
 * \param[in]       cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out]  cb       pointer to a cell builder_t structure
 * \param[in, out]  hhob     pointer to a cs_hho_builder_t structure
 * \param[in, out]  res      vector containing the result
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_builder_compute_dirichlet(const cs_xdef_t         *def,
                                 short int                f,
                                 const cs_cell_mesh_t    *cm,
                                 cs_cell_builder_t       *cb,
                                 cs_hho_builder_t        *hhob,
                                 cs_real_t                res[])
{
  /* Sanity checks */
  if (hhob == NULL || def == NULL)
    return;

  assert(hhob->face_basis != NULL);

  const cs_quant_t  pfq = cm->face[f];
  const cs_basis_func_t  *fbf = hhob->face_basis[f];

  /* See _add_tria_reduction to understand the following shift */
  cs_real_t  *rhs = cb->values + 14 + fbf->size;

  assert(fbf != NULL);
  assert(fbf->facto != NULL);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

  memset(res, 0, fbf->size*sizeof(cs_real_t));
  memset(rhs, 0, fbf->size*sizeof(cs_real_t));

  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->input;

      /* The bc is constant thus its projection is a multiple of the
         constant basis function */
      cs_real_t  phi0;

      fbf->eval_at_point(fbf, pfq.center, 0, 1, &phi0);

      res[0] = constant_val[0] / phi0;
      for (short int i = 1; i < fbf->size; i++)
        res[i] = 0.;

    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

      const int  start = cm->f2e_idx[f];
      const int  end = cm->f2e_idx[f+1];
      const short int n_vf = end - start; // #vertices (=#edges)
      const short int *f2e_ids = cm->f2e_ids + start;

      assert(n_vf > 2);
      switch(n_vf){

      case 3: /* triangle (optimized version, no subdivision) */
        {
          short int  v0, v1, v2;
          cs_cell_mesh_get_next_3_vertices(f2e_ids, cm->e2v_ids, &v0, &v1, &v2);

          const double  *xv0 = cm->xv + 3*v0;
          const double  *xv1 = cm->xv + 3*v1;
          const double  *xv2 = cm->xv + 3*v2;

          _add_tria_reduction(anai, fbf, xv0, xv1, xv2, pfq.meas, cb, rhs);
        }
        break;

      default:
        {
          const double  *tef = cm->tef + start;

          for (short int e = 0; e < n_vf; e++) { /* Loop on face edges */

            // Edge-related variables
            const short int e0  = f2e_ids[e];
            const double  *xv0 = cm->xv + 3*cm->e2v_ids[2*e0];
            const double  *xv1 = cm->xv + 3*cm->e2v_ids[2*e0+1];

            _add_tria_reduction(anai, fbf, xv0, xv1, pfq.center, tef[e],
                                cb, rhs);

          }
        }
        break;

      } /* End of switch */

      /* Modified Cholesky decomposition to compute DoF */
      fbf->project(fbf, rhs, res);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Stop execution.\n"
                " Invalid type of definition.\n"), __func__);

  } // switch def_type

}

/*----------------------------------------------------------------------------*/

#undef _dp3
#undef _mv3

END_C_DECLS
