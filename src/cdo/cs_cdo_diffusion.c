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
  cs_real_3_t  *glv = diff->tmp_vect + cm->n_vc;
  cs_real_t  *lvc = diff->tmp_real;
  cs_real_t  *wvf = diff->tmp_real + cm->n_vc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*cm->n_vc;

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
  cs_real_3_t  *glv = diff->tmp_vect + cm->n_vc;
  cs_real_t  *lvc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + cm->n_vc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*cm->n_vc;

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
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      mnu     property tensor times the face normal
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] diff    pointer to a builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_cost_algo(short int                 f,
                  const cs_cell_mesh_t     *cm,
                  const cs_real_3_t         mnu,
                  double                    beta,
                  cs_cdo_diff_t            *diff,
                  cs_locmat_t              *ntrgrd)
{
  cs_real_3_t  lek;
  cs_real_3_t  *le_grd = diff->tmp_vect;
  cs_real_t  *v_area = diff->tmp_real;

  const cs_quant_t  pfq = cm->face[f];
  const cs_real_t  over_vol_c = 1/cm->vol_c;

  for (short int v = 0; v < cm->n_vc; v++)
    v_area[v] = 0;

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

    /* Half the area of the triangle t_{e,f} defined by x_vi1, x_vi2 and x_f */
    const double  tef = cs_math_surftri(cm->xv + 3*vi1,
                                        cm->xv + 3*vi2,
                                        pfq.center);

    /* Penalization term proportional to the area attached to a border vertex */
    v_area[vi1] += 0.5*tef;
    v_area[vi2] += 0.5*tef;

    /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
    for (short int v = 0; v < cm->n_vc; v++)
      le_grd[v][0] = le_grd[v][1] = le_grd[v][2] = 0;

    /* Term related to the flux reconstruction:
       Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell */
    for (short int ek = 0; ek < cm->n_ec; ek++) {

      const short int  shft = 2*ek;
      const short int  vj1 = cm->e2v_ids[shft];
      const short int  vj2 = cm->e2v_ids[shft+1];
      const short int  sgn_vj1 = cm->e2v_sgn[shft]; // sgn_vj2 = - sgn_vj1
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

      const double  contrib = _dp3(mnu, le_grd[v]) * tef;
      ntrgrd->val[vi1*cm->n_vc + v] += contrib;
      ntrgrd->val[vi2*cm->n_vc + v] += contrib;

    }

  } // border face edges

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the normal trace operator for a given border face when a
 *          WBS algo. is used for reconstructing the degrees of freedom
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in, out] diff      pointer to a builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_wbs_algo(const cs_face_mesh_t     *fm,
                 const cs_cell_mesh_t     *cm,
                 const cs_real_3_t         pty_nuf,
                 cs_cdo_diff_t            *diff,
                 cs_locmat_t              *ntrgrd)
{
  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in diff->tmp_real and diff->tmp_vect */
  cs_real_t  *wvf = diff->tmp_real;
  cs_real_t  *wtef = diff->tmp_real + fm->n_vf;
  cs_real_t  *l_vc = diff->tmp_real + 2*fm->n_vf;
  cs_real_3_t  *mng_ef = diff->tmp_vect;
  cs_real_3_t  *u_vc = diff->tmp_vect + fm->n_vf;

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdc(fm, grd_c);

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); // (pty_tensor * nu_f) . grd_c

  /* First pass: compute useful quantities to build the operator */
  for (short int v = 0; v < fm->n_vf; v++)
    wvf[v] = 0;

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  const cs_quant_t  pfq = fm->face;
  const cs_nvec3_t  deq = fm->dedge;

  /* Compute a weight for each vertex of the current face */
  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const double  tef = cs_compute_tef(e, fm);

    wvf[v1] += tef;
    wvf[v2] += tef;
    wtef[e] = tef * cs_math_onetwelve;

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity */
    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = tef * cs_math_onethird;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

  } /* End of loop on face edges */

  const double  f_coef = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++)
    wvf[v] *= f_coef;

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_i = ntrgrd->val + vi*cm->n_vc;

    /* Default contribution for this line */
    const double  default_coef = pfq.meas * wvf[vfi] * mng_cf;
    for (short int vj = 0; vj < cm->n_vc; vj++)
      ntrgrd_i[vj] = default_coef * cm->wvc[vj];

    /* Block Vf x Vf */
    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

      double  entry_ij = 0.;
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->e2v_ids[2*e];
        const short int  v2 = fm->e2v_ids[2*e+1];

        double  coef_i = wvf[vfi];
        if (vfi == v1 || vfi == v2)
          coef_i += 1;

        double  coef_j = wvf[vfj] * mng_ef[e][2];
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
 *          WBS algo. and a V+C scheme is considered
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t structure
 * \param[in]      pty_nuf   property tensor times the face normal
 * \param[in, out] diff      pointer to a builder structure
 * \param[in, out] ntrgrd    local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_vcb_algo(const cs_face_mesh_t     *fm,
                 const cs_cell_mesh_t     *cm,
                 const cs_real_3_t         pty_nuf,
                 cs_cdo_diff_t            *diff,
                 cs_locmat_t              *ntrgrd)
{
  cs_real_3_t  grd_f, grd_v1, grd_v2, grd_c;

  /* Useful quantities are stored in diff->tmp_real and diff->tmp-vect */
  cs_real_t  *wvf = diff->tmp_real;
  cs_real_t  *wtef = diff->tmp_real + fm->n_vf;
  cs_real_t  *l_vc = diff->tmp_real + 2*fm->n_vf;
  cs_real_3_t  *mng_ef = diff->tmp_vect;
  cs_real_3_t  *u_vc = diff->tmp_vect + fm->n_vf;

  /* Compute the gradient of the Lagrange function related to xc which is
     constant inside p_{f,c} */
  cs_compute_grdc(fm, grd_c);

  const cs_real_t  mng_cf = _dp3(pty_nuf, grd_c); // (pty_tensor * nu_f) . grd_c

  /* First pass: compute useful quantities to build the operator */
  for (short int v = 0; v < fm->n_vf; v++)
    wvf[v] = 0;

  /* Compute xc --> xv length and unit vector for all face vertices */
  for (short int v = 0; v < fm->n_vf; v++)
    cs_math_3_length_unitv(fm->xc, fm->xv + 3*v, l_vc + v, u_vc[v]);

  const cs_quant_t  pfq = fm->face;
  const cs_nvec3_t  deq = fm->dedge;

  /* Compute a weight for each vertex of the current face */
  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const double  tef = cs_compute_tef(e, fm);

    wvf[v1] += tef;
    wvf[v2] += tef;
    wtef[e] = tef * cs_math_onetwelve;

    /* Gradient of the Lagrange function related to v1 and v2 */
    cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                      grd_v1, grd_v2);

    /* Gradient of the Lagrange function related to a face.
       This formula is a consequence of the Partition of the Unity */
    for (int k = 0; k < 3; k++)
      grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

    const double  tef_coef = tef * cs_math_onethird;
    mng_ef[e][0] = _dp3(pty_nuf, grd_v1) * tef_coef;
    mng_ef[e][1] = _dp3(pty_nuf, grd_v2) * tef_coef;
    mng_ef[e][2] = _dp3(pty_nuf, grd_f)  * tef_coef;

  } /* End of loop on face edges */

  const double  f_coef = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++)
    wvf[v] *= f_coef;

  const int n_csys = ntrgrd->n_ent;
  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    short int  vi = fm->v_ids[vfi];
    double  *ntrgrd_i = ntrgrd->val + vi*n_csys;

    /* Contribution to the cell column */
    ntrgrd_i[cm->n_vc] = wvf[vfi] * pfq.meas * mng_cf;

    /* Block Vf x Vf */
    for (short int vfj = 0; vfj < fm->n_vf; vfj++) {

      double  entry_ij = 0.;
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->e2v_ids[2*e];
        const short int  v2 = fm->e2v_ids[2*e+1];

        double  coef_i = wvf[vfi];
        if (vfi == v1 || vfi == v2)
          coef_i += 1;

        double  coef_j = wvf[vfj] * mng_ef[e][2];
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

  int  v_size = CS_MAX(2*connect->n_max_vbyc, connect->n_max_ebyc);
  int  s_size = 3*connect->n_max_vbyc;

  BFT_MALLOC(diff->tmp_vect, v_size, cs_real_3_t);
  BFT_MALLOC(diff->tmp_real, s_size, cs_real_t);

  /* Define a builder for the related discrete Hodge operator */
  diff->hb = NULL;
  if (!hwbs)
    diff->hb = cs_hodge_builder_init(connect, h_info);

  /* Specific members for weakly enforced BCs */
  diff->eig_ratio = -1;
  diff->eig_max = -1;

  int   dof_size = connect->n_max_vbyc; // CS_SPACE_SCHEME_CDOVB
  if (space_scheme == CS_SPACE_SCHEME_CDOVCB)
    dof_size += 1;
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
 * \brief   Define a cell --> Dirichlet boundary faces connectivity
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      dir_face     pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] c2bcbf_idx   pointer to the index to build
 * \param[in, out] c2bcbf_ids   pointer to the list of ids to build
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_build_c2bcbf(const cs_cdo_connect_t    *connect,
                              const cs_cdo_bc_list_t    *dir_face,
                              cs_lnum_t                 *p_c2bcbf_idx[],
                              cs_lnum_t                 *p_c2bcbf_ids[])
{
  cs_lnum_t  *c2bcbf_idx = NULL, *c2bcbf_ids = NULL;

  const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
  const cs_lnum_t  n_cells = connect->c_info->n_elts;
  const cs_sla_matrix_t  *f2c = connect->f2c;

  /* Allocation and initialization */
  BFT_MALLOC(c2bcbf_idx, n_cells + 1, cs_lnum_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells + 1; c_id++)
    c2bcbf_idx[c_id] = 0;

  /* First pass: Loop on Dirichlet faces to build index */
  for (cs_lnum_t i = 0; i < dir_face->n_elts; i++) {

    cs_lnum_t  f_id = dir_face->elt_ids[i] + n_i_faces;
    cs_lnum_t  c_id = f2c->col_id[f2c->idx[f_id]];

    assert(f2c->idx[f_id+1] - f2c->idx[f_id] == 1); // check if border
    c2bcbf_idx[c_id+1] += 1;

  }

  for (cs_lnum_t i = 0; i < n_cells; i++)
    c2bcbf_idx[i+1] += c2bcbf_idx[i];

  /* Second pass: Loop on Dirichlet faces to build list of ids */
  BFT_MALLOC(c2bcbf_ids, c2bcbf_idx[n_cells], cs_lnum_t);

  short int  *count = NULL;
  BFT_MALLOC(count, n_cells, short int);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++)
    count[i] = 0;

  for (cs_lnum_t i = 0; i < dir_face->n_elts; i++) {

    cs_lnum_t  f_id = dir_face->elt_ids[i] + n_i_faces;
    cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];
    cs_lnum_t  shft = c2bcbf_idx[c_id] + count[c_id];

    c2bcbf_ids[shft] = f_id;
    count[c_id] += 1;

  }

  BFT_FREE(count);

  /* Return pointers */
  *p_c2bcbf_idx = c2bcbf_idx;
  *p_c2bcbf_ids = c2bcbf_ids;
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
 * \param[in]       cm        pointer to a cs_cell_mesh_t struct.
 * \param[in]       matpty    3x3 matrix related to the diffusion property
 * \param[in, out]  diff      auxiliary structure used to build the diff. term
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_weak_bc(cs_lnum_t                    f_id,
                         const cs_cell_mesh_t        *cm,
                         const cs_real_t              matpty[3][3],
                         cs_cdo_diff_t               *diff,
                         cs_cdo_locsys_t             *csys)
{
  /* Sanity check */
  assert(diff != NULL);
  assert(diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM ||
         diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE);

  /* Initialize the local quantities */
  cs_locmat_t  *ntrgrd = diff->loc;
  cs_face_mesh_t  *fm = NULL;

  const short int  n_csys = csys->mat->n_ent;
  ntrgrd->n_ent = n_csys;
  for (short int v = 0; v < n_csys; v++)
    ntrgrd->ids[v] = csys->mat->ids[v];
  for (int i = 0; i < n_csys*n_csys; i++)
    ntrgrd->val[i] = 0;

  const cs_param_hodge_t  h_info = diff->h_info;
  const cs_param_hodge_algo_t  h_algo = h_info.algo;

  /* Set the local face id */
  short int f = -1; // not set
  for (short int ii = 0; ii < cm->n_fc; ii++) {
    if (cm->f_ids[ii] == f_id) {
      f = ii;
      break;
    }
  }
  assert(f != -1);

  /* Compute the product: matpty*face unit normal */
  const cs_quant_t  pfq = cm->face[f];
  cs_real_3_t  pty_nuf;
  cs_math_33_3_product((const cs_real_t (*)[3])matpty, pfq.unitv, pty_nuf);

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    /* Compute v_area which is stored in diff->tmp_* for a future use */
    _ntrgrd_cost_algo(f, cm, pty_nuf, h_info.coef, diff, ntrgrd);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    {
      // TODO: Modify this line for a full openMP implementation
      fm = cs_cdo_local_get_face_mesh(0);
      cs_face_mesh_build_from_cell_mesh(cm, f, fm);

      /* Compute useful quantities for the WBS algo. (stored in diff->tmp_*)
         and set the normal trace operator */
      if (diff->scheme == CS_SPACE_SCHEME_CDOVB)
        _ntrgrd_wbs_algo(fm, cm, pty_nuf, diff, ntrgrd);
      else if (diff->scheme == CS_SPACE_SCHEME_CDOVCB)
        _ntrgrd_vcb_algo(fm, cm, pty_nuf, diff, ntrgrd);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  } // End of switch

  if (diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    cs_real_t  *mv = NULL;
    if (h_algo == CS_PARAM_HODGE_ALGO_COST)
      mv = diff->tmp_real + n_csys;
    else
      mv = diff->tmp_real + 2*fm->n_vf;

    /* Update ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
    cs_locmat_add_transpose(ntrgrd, diff->transp);

    /* Update RHS according to the add of transp */
    cs_locmat_matvec(diff->transp, csys->dir_bc, mv);
    for (short int v = 0; v < n_csys; v++)
      csys->rhs[v] += mv[v];

  }

  /* Set the diffusion tensor */
  if (diff->eig_ratio < 0 || diff->is_uniform == false)
    cs_math_33_eigen((const cs_real_t (*)[3])matpty,
                     &(diff->eig_ratio),
                     &(diff->eig_max));

  /* Compute the value of the penalization coefficient */
  const double  f_coef = cs_weak_nitsche_pena_coef * pow(pfq.meas, -0.5) *
    diff->eig_ratio * diff->eig_max;
  assert(f_coef > 0); // Sanity check

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    { // v_area is computed inside subroutine _ntrgrd_cost_algo
      const cs_real_t  *v_area = diff->tmp_real;

      for (short int v = 0; v < ntrgrd->n_ent; v++) {
        if (v_area[v] > 0) { // v belongs to f

          // Value of the surfacic diagonal Hodge operator
          const double p_coef = f_coef * v_area[v];

          // Set the penalty diagonal coefficient
          ntrgrd->val[v + v*ntrgrd->n_ent] += p_coef;
          csys->rhs[v] += p_coef * csys->dir_bc[v];

        }
      }
    }
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    {
      assert(fm != NULL);

      // wvf and tef are computed inside subroutine _ntrgrd_*_algo
      double  *wvf = diff->tmp_real;
      double  *wtef = diff->tmp_real + fm->n_vf;

      cs_locmat_t  *hloc = diff->transp;

      /* Build the border Hodge operator */
      for (int i = 0; i < n_csys*n_csys; i++)
        hloc->val[i] = 0;

      for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

        const short int vi = fm->v_ids[vfi];
        double  *hi = hloc->val + vi*n_csys;

        /* Default contribution */
        const double  default_coef = 0.5 * wvf[vfi] * pfq.meas;
        for (short int vfj = 0; vfj < fm->n_vf; vfj++)
          hi[fm->v_ids[vfj]] = default_coef * wvf[vfj];

        /* Specific diagonal contribution */
        hi[vi] += 2 * default_coef * cs_math_onethird;

      } // Loop on face vertices

      /* Specific extra-diag contribution */
      for (short int e = 0; e < fm->n_ef; e++) {

        const short int  v1 = fm->v_ids[fm->e2v_ids[2*e]];
        const short int  v2 = fm->v_ids[fm->e2v_ids[2*e+1]];

        hloc->val[v1*n_csys + v2] += wtef[e];
        hloc->val[v2*n_csys + v1] += wtef[e];

      } /* Loop on face edges */

      /* Add the border Hodge op. to the normal trace op.
         Update RHS whith H*p^{dir} */
      for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

        const short int  vi = fm->v_ids[vfi];
        const double  *hi = hloc->val + vi*n_csys;
        const double pii_coef = f_coef * hi[vi];

        double  *ntrg_i = ntrgrd->val + vi*n_csys;

        // Set the penalty diagonal coefficient
        ntrg_i[vi] += pii_coef;
        csys->rhs[vi] += pii_coef * csys->dir_bc[vi];

        for (short int vfj = vfi+1; vfj < fm->n_vf; vfj++) {

          const short int  vj = fm->v_ids[vfj];
          const double  pij_coef = f_coef * hi[vj];

          ntrg_i[vj] += pij_coef;
          ntrgrd->val[vi + vj*n_csys] += pij_coef;
          csys->rhs[vi] += pij_coef * csys->dir_bc[vj];
          csys->rhs[vj] += pij_coef * csys->dir_bc[vi];

        } // Loop on face vertices vj

      } // Loop on face vertices vi

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
  cs_locmat_add(csys->mat, ntrgrd);
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
