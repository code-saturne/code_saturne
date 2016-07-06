/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based and vertex+cell-based schemes
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
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo_scheme_geometry.h"
#include "cs_math.h"
#include "cs_property.h"

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

#define CS_CDO_DIFFUSION_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/* Stiffness matrix builder */
struct _cs_cdo_diff_t {

  cs_space_scheme_t      scheme;  // space discretization scheme
  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs

  /* Data related to the discrete Hodge operator attached to the
     diffusion property */
  bool                  is_uniform;  // Diffusion tensor is uniform ?
  cs_param_hodge_t      h_info;      // Parameters defining a discrete Hodge op.
  cs_hodge_builder_t   *hb;          // Helper for building a discrete Hodge op.

  /* Temporary buffers */
  cs_real_3_t    *tmp_vect;  // set of local vectors
  cs_real_t      *tmp_real;  // set of local arrays of double

  /* Specific members for weakly enforced BCs */
  double          eig_ratio; // ratio of the eigenvalues of the diffusion tensor
  double          eig_max;   // max. value among eigenvalues
  cs_locmat_t    *transp;    // Specific to the symmetric treatment

  /* Local matrix (stiffness or normal trace gradient) */
  cs_locmat_t    *loc;

};

/*============================================================================
 * Local variables
 *============================================================================*/

// Advanced developper parameters (weakly enforced the boundary conditions)
static const double  cs_weak_nitsche_pena_coef = 500;
static const double  cs_wbs_tef_weight = 1/24.;

/*============================================================================
 * Private constant variables
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stiffness matrix on this cell using a Whitney
 *          Barycentric Subdivision (WBS) algo.
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      cm          pointer to a cs_cell_mesh_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_locmat_t *
_compute_wbs_stiffness(const cs_cdo_quantities_t   *quant,
                       const cs_cell_mesh_t        *cm,
                       const cs_real_3_t           *tensor,
                       cs_cdo_diff_t               *diff)
{
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *uvc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + cm->n_max_vbyc;
  cs_real_t  *lvc = diff->tmp_real;
  cs_real_t  *wvf = diff->tmp_real + cm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*cm->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*cm->c_id;

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(xc, xyz + 3*cm->v_ids[v], lvc + v, uvc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute the gradient of the lagrange function related to a cell
       in each p_{f,c} and the weights for each vertex related to this face */
    cs_compute_fwbs_q2(f, cm, grd_c, wvf, pefc_vol);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = cm->f2e_idx[f], jj = 0; i < cm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 and v2 */
      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])uvc, lvc,
                        grd_v1, grd_v2);

      /* Gradient of the lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < sloc->n_ent; si++) {

        for (int k = 0; k < 3; k++)
          glv[si][k] = cm->wvc[si]*grd_c[k];

        if (wvf[si] > 0) // Face contrib.
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) // Vertex 1 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) // Vertex 2 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } // Loop on cell vertices

      /* Build the upper right part */
      for (int si = 0; si < sloc->n_ent; si++) {

        cs_math_33_3_product(tensor, glv[si], matg);

        /* Diagonal contribution */
        double  *mi = sloc->val + si*sloc->n_ent;
        mi[si] += subvol * _dp3(matg, glv[si]);

        /* Loop on vertices v_j (j > i) */
        for (int sj = si+1; sj < sloc->n_ent; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } // Loop on cell faces

  /* Matrix is symmetric by construction */
  for (int si = 0; si < sloc->n_ent; si++) {
    double  *mi = sloc->val + si*sloc->n_ent;
    for (int sj = si+1; sj < sloc->n_ent; sj++)
      sloc->val[sj*sloc->n_ent+si] = mi[sj];
  }

  return sloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stiffness matrix on the current cell using a Whitney
 *          Barycentric Subdivision (WBS) algo and a V+C CDO scheme.
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      cm          pointer to a cs_cell_mesh_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_locmat_t *
_compute_vcb_stiffness(const cs_cdo_quantities_t   *quant,
                       const cs_cell_mesh_t        *cm,
                       const cs_real_3_t           *tensor,
                       cs_cdo_diff_t               *diff)
{
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg, matg_c;

  cs_real_3_t  *uvc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + cm->n_max_vbyc;
  cs_real_t  *lvc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + cm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*cm->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const int  msize = cm->n_vc + 1;
  const int  cc = msize*cm->n_vc + cm->n_vc;
  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*cm->c_id;

  assert(sloc->n_ent == msize);

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(xc, xyz + 3*cm->v_ids[v], lvc + v, uvc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute for the current face:
       - the gradient of the Lagrange function related xc in p_{f,c}
       - weights related to vertices
       - subvolume p_{ef,c} related to edges
    */
    const double  pfc_vol = cs_compute_fwbs_q3(f, cm, grd_c, wvf, pefc_vol);

    /* Compute the contribution to the entry A(c,c) */
    cs_math_33_3_product(tensor, grd_c, matg_c);
    sloc->val[cc] += pfc_vol * _dp3(grd_c, matg_c);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = cm->f2e_idx[f], jj = 0; i < cm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 and v2 */
      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])uvc, lvc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < cm->n_vc; si++) {

        for (int k = 0; k < 3; k++)
          glv[si][k] = 0;

        if (wvf[si] > 0) // Face contrib.
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) // Vertex 1 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) // Vertex 2 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } // Loop on cell vertices

      /* Build the upper right part (v-v and v-c)
         Be careful: sloc->n_ent = cm->n_vc + 1 */
      for (int si = 0; si < cm->n_vc; si++) {

        double  *mi = sloc->val + si*sloc->n_ent;

        /* Add v-c contribution */
        mi[cm->n_vc] += subvol * _dp3(matg_c, glv[si]);

        /* Add v-v contribution */
        cs_math_33_3_product(tensor, glv[si], matg);

        /* Loop on vertices v_j (j >= i) */
        for (int sj = si; sj < cm->n_vc; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } // Loop on cell faces

  /* Matrix is symmetric by construction */
  for (int si = 0; si < sloc->n_ent; si++) {
    double  *mi = sloc->val + si*sloc->n_ent;
    for (int sj = si+1; sj < sloc->n_ent; sj++)
      sloc->val[sj*sloc->n_ent+si] = mi[sj];
  }

  return sloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          COST algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      f       face id in the cellwise numbering
 * \param[in]      mn      property tensor times the face normal
 * \param[in]      xf      face center
 * \param[in]      xyz     vertex coordinates
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] diff    pointer to a builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_cost_algo(short int                 f,
                  const cs_real_3_t         mn,
                  const cs_real_t          *xf,
                  const cs_real_t          *xyz,
                  const cs_cell_mesh_t     *cm,
                  double                    beta,
                  cs_cdo_diff_t            *diff,
                  cs_locmat_t              *ntrgrd)
{
  cs_real_3_t  lek;
  cs_real_3_t  *leg = diff->tmp_vect;
  cs_real_t  *v_coef = diff->tmp_real;

  const cs_real_t  over_vol_c = 1/cm->vol_c;

  for (short int v = 0; v < cm->n_vc; v++)
    v_coef[v] = 0;

  /* Loop on border face edges */
  for (int jf = cm->f2e_idx[f]; jf < cm->f2e_idx[f+1]; jf++) {

    short int  ei = cm->f2e_ids[jf];
    short int  vi1 = cm->e2v_ids[2*ei];
    short int  vi2 = cm->e2v_ids[2*ei+1];
    cs_quant_t  peq_i = cm->edge[ei];
    cs_nvec3_t  dfq_i = cm->dface[ei];

    const double  dp = _dp3(peq_i.unitv, dfq_i.unitv);
    const double  tmp_val = peq_i.meas * dfq_i.meas * dp; // 3*pec_vol
    const double  beta_pec_vol = 3. * beta/tmp_val;
    const double  coef_ei =  3. * beta/dp;

    /* Area of the triangle t_{e,f} defined by x_vi1, x_vi2 and x_f */
    const double  surf = cs_math_surftri(xyz + 3*cm->v_ids[vi1],
                                         xyz + 3*cm->v_ids[vi2],
                                         xf);

    /* Penalization term */
    v_coef[vi1] += 0.5*surf;
    v_coef[vi2] += 0.5*surf;

    /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
    for (short int v = 0; v < cm->n_vc; v++)
      leg[v][0] = leg[v][1] = leg[v][2] = 0;

    /* Term related to the flux reconstruction:
       Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell */
    for (short int ek = 0; ek < cm->n_ec; ek++) {

      const short int  shft = 2*ek;
      const short int  vj1 = cm->e2v_ids[shft];
      const short int  vj2 = cm->e2v_ids[shft+1];
      const short int  sgn_vj1 = cm->e2v_sgn[shft];
      const cs_nvec3_t  dfq_k = cm->dface[ek];

      /* Compute l_(ek,c)|p(_ei,f,c) */
      const double  eik_part = coef_ei * _dp3(dfq_k.unitv, peq_i.unitv);
      const double  coef_mult = dfq_k.meas * over_vol_c;

      for (int k = 0; k < 3; k++)
        lek[k] = coef_mult * (dfq_k.unitv[k] - eik_part * dfq_i.unitv[k]);

      if (ek == ei)
        for (int k = 0; k < 3; k++)
          lek[k] += dfq_k.meas * dfq_k.unitv[k] * beta_pec_vol;

      for (int k = 0; k < 3; k++) {
        leg[vj1][k] += sgn_vj1 * lek[k];
        leg[vj2][k] -= sgn_vj1 * lek[k];
      }

    } // Loop on cell edges

    for (short int v = 0; v < cm->n_vc; v++) {

      const double  contrib = _dp3(mn, leg[v]) * surf;
      ntrgrd->val[vi1*cm->n_vc + v] += contrib;
      ntrgrd->val[vi2*cm->n_vc + v] += contrib;

    }

  } // border face edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute useful quantities for building operators when WBS algo.
 *          is requested
 *          wvf, mng, tef/12
 *
 * \param[in]      f       face id in the cellwise numbering
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      xyz     vertex coordinates
 * \param[in]      xc      cell center
 * \param[in]      mn      property tensor times the face normal
 * \param[in, out] grd_c   gradient for the Lagrange function related to x_c
 * \param[in, out] diff    pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_quant_wbs_algo(short int                 f,
                    const cs_cell_mesh_t     *cm,
                    const cs_real_t          *xyz,
                    const cs_real_t          *xc,
                    const cs_real_3_t         mn,
                    cs_real_3_t               grd_c,
                    cs_cdo_diff_t            *diff)
{
  double  len;
  cs_real_3_t  un, cp, grd_f, grd_v1, grd_v2;

  /* Useful quantities are stored in diff->tmp_real and diff->tmp-vect */
  cs_real_3_t  *mng = diff->tmp_vect;
  cs_real_3_t  *uvc = diff->tmp_vect + cm->n_vc;
  cs_real_t  *wvf = diff->tmp_real;
  cs_real_t  *lvc = diff->tmp_real + cm->n_vc;
  cs_real_t  *wtef = diff->tmp_real + 2*cm->n_vc;

  const cs_nvec3_t  deq = cm->dedge[f];
  const cs_quant_t  pfq = cm->face[f];
  const double  f_coef = 0.25/pfq.meas;

  /* Compute the gradient of the Lagrange function related to xc
     which is constant inside p_{f,c} */
  const double  hf = _dp3(pfq.unitv, deq.unitv) * deq.meas;
  const cs_real_t  ohf = -cm->f_sgn[f]/hf;
  for (int k = 0; k < 3; k++)
    grd_c[k] = ohf * pfq.unitv[k];

  for (short int v = 0; v < cm->n_vc; v++) {
    cs_math_3_length_unitv(xc, xyz + 3*cm->v_ids[v], lvc + v, uvc[v]);
    wvf[v] = 0;
  }

  /* Compute a weight for each vertex of the current face */
  for (int ii = 0, i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++, ii++) {

    const short int  e = cm->f2e_ids[i];
    const cs_quant_t  peq = cm->edge[e];
    const short int  v1 = cm->e2v_ids[2*e];
    const short int  v2 = cm->e2v_ids[2*e+1];

    cs_math_3_length_unitv(peq.center, pfq.center, &len, un);
    cs_math_3_cross_product(un, peq.unitv, cp);

    const double  tef2 = len * peq.meas * cs_math_3_norm(cp); // 2*|t_{e,f}|
    const double  wvf_contrib = tef2 * f_coef;
    const double  tef_coef = tef2 * cs_wbs_tef_weight;  // |t_{e,f}| / 12

    wvf[v1] += wvf_contrib;
    wvf[v2] += wvf_contrib;
    wtef[ii] = tef_coef;

    /* Gradient of the lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])uvc, lvc,
                      grd_v1, grd_v2);

    /* Gradient of the lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity */
    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    mng[ii][0] = _dp3(mn, grd_v1) * tef_coef;
    mng[ii][1] = _dp3(mn, grd_v2) * tef_coef;
    mng[ii][2] = _dp3(mn, grd_f) * tef_coef;

  } /* End of loop on face edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          COST algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      f       face id in the cellwise numbering
 * \param[in]      mngc      property tensor times the face normal
 * \param[in]      xf      face center
 * \param[in]      xyz     vertex coordinates
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] diff    pointer to a builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_wbs_algo(short int                 f,
                 double                    f_meas,
                 double                    mngc,
                 const cs_cell_mesh_t     *cm,
                 cs_cdo_diff_t            *diff,
                 cs_locmat_t              *ntrgrd)
{
  /* Useful quantities are stored in diff->tmp_real and diff->tmp-vect */
  cs_real_3_t  *mng = diff->tmp_vect;
  cs_real_t  *wvf = diff->tmp_real;

  double  sum_ef = 0.;
  for (int jf = cm->f2e_idx[f], ii = 0; jf < cm->f2e_idx[f+1]; jf++, ii++)
    sum_ef += mng[ii][2];

  for (short int vi = 0; vi < cm->n_vc; vi++) {
    if (wvf[vi] > 0) {

      double  *ntrgrd_i = ntrgrd->val + vi*cm->n_vc;

      for (short int vj = 0; vj < cm->n_vc; vj++) {

        /* Default contribution */
        double contrib = 0.25 * f_meas * wvf[vi] * cm->wvc[vj] * mngc; // 1ab

        if (wvf[vj] > 0) {

          contrib += wvf[vj] * wvf[vi] * sum_ef; // 2b

          for (int jf = cm->f2e_idx[f], ii = 0; jf < cm->f2e_idx[f+1];
               jf++, ii++) {

            const short int  e = cm->f2e_ids[jf];
            const short int  v1 = cm->e2v_ids[2*e];
            const short int  v2 = cm->e2v_ids[2*e+1];

            if (vi == v1 || vi == v2) {
              contrib += wvf[vj] * mng[ii][2]; // 2a
              if (vj == v1)
                contrib += mng[ii][0]; // 3a
              if (vj == v2)
                contrib += mng[ii][1]; // 3a
            }
            if (vj == v1)
              contrib += wvf[vi] * mng[ii][0]; // 3b
            if (vj == v2)
              contrib += wvf[vi] * mng[ii][1]; // 3b

          } // Loop on face edges

        } // vj belongs to f

        ntrgrd_i[vj] = contrib;

      } // Loop on cell vertices (vj)

    } // vi belongs to f
  } // Loop on cell vertices (vi)

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure used to build the stiffness matrix
 *
 * \param[in] connect       pointer to a cs_cdo_connect_t structure
 * \param[in] space_scheme  scheme used for discretizing in space
 * \param[in] is_uniform    diffusion tensor is uniform ? (true or false)
 * \param[in] h_info        cs_param_hodge_t structure
 * \param[in] bc_enforce    type of boundary enforcement for Dirichlet values
 *
 * \return a pointer to a new allocated cs_cdo_diff_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_diff_t *
cs_cdo_diffusion_builder_init(const cs_cdo_connect_t       *connect,
                              cs_space_scheme_t             space_scheme,
                              bool                          is_uniform,
                              const cs_param_hodge_t        h_info,
                              const cs_param_bc_enforce_t   bc_enforce)
{
  cs_cdo_diff_t  *diff = NULL;

  BFT_MALLOC(diff, 1, cs_cdo_diff_t);

  /* Store generic for a straightforward access */
  diff->is_uniform = is_uniform;
  diff->enforce = bc_enforce;
  diff->scheme = space_scheme;

  /* Copy the data related to a discrete Hodge operator */
  diff->h_info.type    = h_info.type;
  diff->h_info.inv_pty = h_info.inv_pty;
  diff->h_info.algo    = h_info.algo;
  diff->h_info.coef    = h_info.coef;

  bool  wnit = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) ? true : false;
  bool  wsym = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;
  bool  hwbs = (h_info.algo == CS_PARAM_HODGE_ALGO_WBS) ? true : false;

  int   dof_size = connect->n_max_vbyc; // CS_SPACE_SCHEME_CDOVB
  if (space_scheme == CS_SPACE_SCHEME_CDOVCB)
    dof_size += 1;

  diff->hb = NULL;
  diff->tmp_real = NULL;
  diff->tmp_vect = NULL;

  int  v_size = 0;
  int  s_size = connect->n_max_ebyc;
  if (wnit || wsym) { // Weak enforcement of Dirichlet BCs
    v_size = connect->n_max_ebyc + connect->n_max_vbyc;
    s_size = connect->n_max_ebyc + connect->n_max_vbyc;
  }

  if (hwbs) {
    v_size = CS_MAX(v_size, 2*connect->n_max_ebyc);
    s_size = CS_MAX(s_size,
                    2*connect->n_max_vbyc + connect->n_max_vbyf);
  }
  else {

    /* Define a builder for the related discrete Hodge operator */
    diff->hb = cs_hodge_builder_init(connect, h_info);

  }

  if (v_size > 0) BFT_MALLOC(diff->tmp_vect, v_size, cs_real_3_t);
  if (s_size > 0) BFT_MALLOC(diff->tmp_real, s_size, cs_real_t);

  /* Specific members for weakly enforced BCs */
  diff->eig_ratio = -1;
  diff->eig_max = -1;
  if (wsym || (wnit && hwbs))
    diff->transp = cs_locmat_create(dof_size);

  /* Allocate the local stiffness matrix */
  diff->loc = cs_locmat_create(dof_size);

  return diff;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdo_diff_t structure
 *
 * \param[in, out ] diff   pointer to a cs_cdo_diff_t struc.
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

cs_cdo_diff_t *
cs_cdo_diffusion_builder_free(cs_cdo_diff_t   *diff)
{
  if (diff == NULL)
    return diff;

  cs_param_bc_enforce_t  bc_enforce = diff->enforce;
  bool  wnit = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) ? true : false;
  bool  wsym = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;
  bool  hwbs = (diff->h_info.algo == CS_PARAM_HODGE_ALGO_WBS) ? true : false;

  BFT_FREE(diff->tmp_vect);
  BFT_FREE(diff->tmp_real);

  if (!hwbs)
    diff->hb = cs_hodge_builder_free(diff->hb);

  if (wsym || (wnit && hwbs))
    diff->transp = cs_locmat_free(diff->transp);

  /* Local stiffness matrix */
  diff->loc = cs_locmat_free(diff->loc);

  BFT_FREE(diff);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the related Hodge builder structure
 *
 * \param[in]  diff   pointer to a cs_cdo_diff_t structure
 *
 * \return  a pointer to a cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_cdo_diffusion_get_hodge_builder(cs_cdo_diff_t   *diff)
{
  if (diff == NULL)
    return NULL;

  return diff->hb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get temporary buffers attached to a cs_cdo_diff_t structure
 *
 * \param[in]       diff     pointer to a cs_cdo_diff_t structure
 * \param[in, out]  tmp_vec  pointer to a buffer of cs_real_3_t
 * \param[in, out]  tmp_sca  pointer to a buffer of cs_real_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_get_tmp_buffers(const cs_cdo_diff_t   *diff,
                                 cs_real_3_t          **tmp_vec,
                                 cs_real_t            **tmp_sca)
{
  *tmp_vec = NULL;
  *tmp_sca = NULL;

  if (diff == NULL)
    return;

  *tmp_vec = diff->tmp_vect;
  *tmp_sca = diff->tmp_real;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) stiffness matrix
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      cm          cell-wise connectivity and quantitites
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdo_diffusion_build_local(const cs_cdo_quantities_t   *quant,
                             const cs_cell_mesh_t        *cm,
                             const cs_real_3_t           *tensor,
                             cs_cdo_diff_t               *diff)
{
  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  /* Initialize the local matrix */
  sloc->n_ent = cm->n_vc;
  for (int i = 0; i < cm->n_vc; i++)
    sloc->ids[i] = cm->v_ids[i];
  if (diff->scheme == CS_SPACE_SCHEME_CDOVCB) {
    sloc->n_ent += 1;
    sloc->ids[cm->n_vc] = cm->c_id;
  }

  for (int i = 0; i < sloc->n_ent*sloc->n_ent; i++)
    sloc->val[i] = 0;

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:

    /* Sanity check */
    assert(diff->scheme != CS_SPACE_SCHEME_CDOVCB);

    /* Set the diffusion tensor if needed */
    if (cs_hodge_builder_get_setting_flag(diff->hb) == false ||
        diff->is_uniform == false)
      cs_hodge_builder_set_tensor(diff->hb, tensor);

    /* Build a local discrete Hodge op. and return a local dense matrix */
    cs_hodge_build_local_stiffness(cm, diff->hb, sloc);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    if (diff->scheme == CS_SPACE_SCHEME_CDOVB)
      sloc = _compute_wbs_stiffness(quant, cm, tensor, diff);
    else if (diff->scheme == CS_SPACE_SCHEME_CDOVCB)
      sloc = _compute_vcb_stiffness(quant, cm, tensor, diff);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid space scheme for building the stiffness matrix.\n"
                  " Please check your settings."));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid algorithm for building the local stiffness matrix.");

  } // End of switch

  return sloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) "normal trace gradient" matrix taking
 *          into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique (symmetrized or not)
 *
 * \param[in]       f_id      face id (a border face attached to a Dir. BC)
 * \param[in]       quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]       cm        pointer to a cs_cell_mesh_t struct.
 * \param[in]       matpty    3x3 matrix related to the diffusion property
 * \param[in, out]  diff      auxiliary structure used to build the diff. term
 * \param[in, out]  ls        cell-wise structure storing the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_weak_bc(cs_lnum_t                    f_id,
                         const cs_cdo_quantities_t   *quant,
                         const cs_cell_mesh_t        *cm,
                         const cs_real_t              matpty[3][3],
                         cs_cdo_diff_t               *diff,
                         cs_cdo_locsys_t             *ls)
{
  /* Sanity check */
  assert(diff != NULL);
  assert(diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM ||
         diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE);

  const cs_param_hodge_t  h_info = diff->h_info;
  const cs_param_hodge_algo_t  h_algo = h_info.algo;

  /* Set the diffusion tensor */
  if (diff->eig_ratio < 0 || diff->is_uniform == false)
    cs_math_33_eigen((const cs_real_t (*)[3])matpty,
                     &(diff->eig_ratio),
                     &(diff->eig_max));

  /* Set the local face id */
  short int f = -1; // not set
  for (short int ii = 0; ii < cm->n_fc; ii++) {
    if (cm->f_ids[ii] == f_id) {
      f = ii;
      break;
    }
  }
  assert(f != -1);

  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*cm->c_id;
  const cs_quant_t  pfq = cm->face[f];
  const double  f_coef = cs_weak_nitsche_pena_coef * pow(pfq.meas, -0.5) *
    diff->eig_ratio * diff->eig_max;

  assert(f_coef > 0); // Sanity check

  /* Initialize the local quantities */
  cs_locmat_t  *ntrgrd = diff->loc;

  ntrgrd->n_ent = cm->n_vc;
  for (short int v = 0; v < cm->n_vc; v++)
    ntrgrd->ids[v] = cm->v_ids[v];

  for (int i = 0; i < cm->n_vc*cm->n_vc; i++)
    ntrgrd->val[i] = 0;

  /* Compute the product: matpty*face unit normal */
  cs_real_3_t  mn;
  cs_math_33_3_product((const cs_real_t (*)[3])matpty, pfq.unitv, mn);

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    _ntrgrd_cost_algo(f, mn, pfq.center, xyz, cm, h_info.coef, diff, ntrgrd);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    {
      cs_real_3_t  grd_c;

      /* Compute useful quantities for the WBS algo. (stored in diff->tmp_*) */
      _get_quant_wbs_algo(f, cm, xyz, xc, mn, grd_c, diff);

      cs_real_t  mng_c = _dp3(mn, grd_c); // (pty_tensor * nu_f) . grd_c

      _ntrgrd_wbs_algo(f, pfq.meas, mng_c, cm, diff, ntrgrd);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  } // End of switch

  if (diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    cs_real_t  *mv = diff->tmp_real + cm->n_vc;

    /* ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
    cs_locmat_add_transpose(ntrgrd, diff->transp);

    /* Update RHS according to the add of transp */
    cs_locmat_matvec(diff->transp, ls->dir_bc, mv);
    for (short int v = 0; v < ntrgrd->n_ent; v++)
      ls->rhs[v] += mv[v];

  }

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    { // v_coef is computed inside _ntrgrd_cost_algo
      const cs_real_t  *v_coef = diff->tmp_real;

      for (short int v = 0; v < ntrgrd->n_ent; v++) {
        if (v_coef[v] > 0) {

          const double p_coef = f_coef * v_coef[v];

          // Set the penalty diagonal coefficient
          ntrgrd->val[v + v*ntrgrd->n_ent] += p_coef;
          ls->rhs[v] += p_coef * ls->dir_bc[v];

        }
      }
    }
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    {
      double  *wvf = diff->tmp_real;
      double  *wtef = diff->tmp_real + 2*cm->n_vc;
      cs_locmat_t  *hloc = diff->transp;

      /* Build the border Hodge operator */
      for (int i = 0; i < cm->n_vc*cm->n_vc; i++)
        hloc->val[i] = 0;

      for (short int vi = 0; vi < cm->n_vc; vi++) {
        if (wvf[vi] > 0) {

          /* Diagonal contribution */
          double  *hi = hloc->val + vi*cm->n_vc;

          hi[vi] = pfq.meas * wvf[vi] * (wvf[vi]*0.5 + cs_math_onethird);

          /* Extra-diagonal contribution */
          for (short int vj = vi + 1; vj < cm->n_vc; vj++) {
            if (wvf[vj] > 0) {

              hi[vj] = 0.5 * wvf[vi] * wvf[vj] * pfq.meas;

              for (int jf = cm->f2e_idx[f], ii = 0; jf < cm->f2e_idx[f+1];
                   jf++, ii++) {

                const short int  e = cm->f2e_ids[jf];
                const short int  v1 = cm->e2v_ids[2*e];
                const short int  v2 = cm->e2v_ids[2*e+1];

                if (v1 == vi && v2 == vj)
                  hi[vj] += wtef[ii];
                else if (v1 == vj && v2 == vi)
                  hi[vj] += wtef[ii];

              } /* Loop on face edges */

              /* Matrix is symmetric (used only the upper part) */
              // hloc->val[vj*cm->n_vc + vi] = hi[vj];

            } // vj belongs to f
          } // Loop on cell vertices vj

        } // vi belongs to f
      } // Loop on cell vertices vi

      /* Add the border Hodge op. to the normal trace op.
         Update RHS whith H*p^{dir} */
      for (short int vi = 0; vi < cm->n_vc; vi++) {
        if (wvf[vi] > 0) {

          double  *ntrg_i = ntrgrd->val + vi*cm->n_vc;
          const double  *hi = hloc->val + vi*cm->n_vc;
          const double pii_coef = f_coef * hi[vi];

          // Set the penalty diagonal coefficient
          ntrg_i[vi] += pii_coef;
          ls->rhs[vi] += pii_coef * ls->dir_bc[vi];

          for (short int vj = vi+1; vj < cm->n_vc; vj++) {
            if (wvf[vj] > 0) {

              const double  pij_coef = f_coef * hi[vj];

              ntrg_i[vj] += pij_coef;
              ntrgrd->val[vi + vj*cm->n_vc] += pij_coef;
              ls->rhs[vi] += pij_coef * ls->dir_bc[vj];
              ls->rhs[vj] += pij_coef * ls->dir_bc[vi];

            } // vj belongs to f
          } // Loop on cell vertices vj

        } // vi belongs to f
      } // Loop on cell vertices vi

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_DIFFUSION_DBG > 1
  bft_printf(">> Local weak bc matrix (f_id: %d)", f_id);
  cs_locmat_dump(cm->c_id, ntrgrd);
#endif

  /* Add contribution to the linear system */
  cs_locmat_add(ls->mat, ntrgrd);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the diffusive flux across dual faces for a given cell
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]       cm        pointer to a cs_face_mesh_t structure
 * \param[in]       dface     pointer to the dual faces related to cell edges
 * \param[in]       pty_tens  3x3 matrix related to the diffusion property
 * \param[in]       p_v       array of values attached to face vertices
 * \param[in]       p_c       value attached to the cell
 * \param[in, out]  diff      auxiliary structure dedicated to diffusion
 * \param[in, out]  c_flux    flux across dual faces inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_cellwise_flux(const cs_cell_mesh_t      *cm,
                               const cs_dface_t          *dface,
                               const cs_real_t            pty_tens[3][3],
                               const double              *p_v,
                               const double               p_c,
                               cs_cdo_diff_t             *diff,
                               double                    *c_flux)
{
  cs_real_3_t  grd_c, grd_v1, grd_v2, grd_pef, mgrd;

  cs_real_3_t  *u_vc = diff->tmp_vect;
  double  *l_vc = diff->tmp_real;
  double  *tef = diff->tmp_real + cm->n_vc;

  /* Reset local fluxes */
  for (short int e = 0; e < cm->n_ec; e++)
    c_flux[e] = 0.;

  /* Store segments xv --> xc for this cell */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    cs_nvec3_t  deq = cm->dedge[f];

    /* Compute for the current face:
       - the area of each triangle defined by a base e and an apex f
       - the gradient of the Lagrange function related xc in p_{f,c}
    */
    cs_compute_tef_grdc(f, cm, tef, grd_c);

    /* Compute the reconstructed value of the potential at p_f */
    double  p_f = 0.;
    for (int i = cm->f2e_idx[f], ii = 0; i < cm->f2e_idx[f+1]; i++, ii++) {

      const short int  eshft = 2*cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[eshft];
      const short int  v2 = cm->e2v_ids[eshft+1];

      p_f += tef[ii]*(p_v[v1] + p_v[v2]);
    }
    p_f *= 0.5/cm->face[f].meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      const double  dp_v1f = p_v[v1] - p_f;
      const double  dp_v2f = p_v[v2] - p_f;

      /* Gradient of the lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c}
      */
      for (int k = 0; k < 3; k++)
        grd_pef[k] = dp_cf*grd_c[k] + dp_v1f*grd_v1[k] + dp_v2f*grd_v2[k];

      cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, grd_pef, mgrd);

      const cs_dface_t  dfq = dface[e];
      if (cm->f_ids[f] == dfq.parent_id[0])
        c_flux[e] -= dfq.sface[0].meas * _dp3(dfq.sface[0].unitv, mgrd);
      else {
        assert(cm->f_ids[f] == dfq.parent_id[1]);
        c_flux[e] -= dfq.sface[1].meas * _dp3(dfq.sface[1].unitv, mgrd);
      }

    } // Loop on face edges

  } // Loop on cell face

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
 * \param[in, out]  diff      auxiliary structure dedicated to diffusion
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
                           cs_cdo_diff_t             *diff)
{
  cs_real_3_t  grd_c, grd_v1, grd_v2, grd_pef, mnuf;

  double  f_flux = 0.;

  /* Retrieve temporary buffers */
  double  *l_vc = diff->tmp_real;
  cs_real_3_t  *u_vc = diff->tmp_vect;

  cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, fm->face.unitv,
                       mnuf);

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Compute for the current face, the gradient of the Lagrange function
     related xc in p_{f,c} */
  cs_compute_grdc(fm, grd_c);

  /* Compute p_c - p_f (where p_c is the reconstructed values at the
     cell center */
  const double  dp_cf = p_c - p_f;
  const cs_nvec3_t  deq = fm->dedge;

  /* Loop on face edges to scan p_{ef,c} subvolumes */
  for (int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* Compute the gradient of the Lagrange function related xv1, xv2
       in each p_{ef,c} for e in E_f */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    const double  dp_v1f = p_v[v1] - p_f;
    const double  dp_v2f = p_v[v2] - p_f;

    /* Gradient of the lagrange function related to a face.
       grd_f = -(grd_c + grd_v1 + grd_v2)
       This formula is a consequence of the Partition of the Unity.
       This yields the following formula for grd(Lv^conf)|_p_{ef,c}
    */
    for (int k = 0; k < 3; k++)
      grd_pef[k] = dp_cf*grd_c[k] + dp_v1f*grd_v1[k] + dp_v2f*grd_v2[k];

    /* Area of the triangle defined by the base e and the apex f */
    f_flux -= cs_compute_tef(e, fm) * _dp3(mnuf, grd_pef);

  } // Loop on face edges

  return f_flux;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
