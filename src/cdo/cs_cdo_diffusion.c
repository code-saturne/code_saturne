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
 * \brief   Compute the gradient of a Lagrange hat function related to a cell c
 *          in a p_{f, c} subvolume of a cell where f is a face of c
 *
 * \param[in]       ifc        incidence number i_{f,c}
 * \param[in]       pfq        primal face quantities
 * \param[in]       deq        dual edge quantities
 * \param[in, out]  grd_func   gradient of the Lagrange function related to c
 */
/*----------------------------------------------------------------------------*/

inline static void
_grd_c_lagrange_pfc(short int             ifc,
                    const cs_quant_t      pfq,
                    const cs_nvec3_t      deq,
                    cs_real_3_t           grd_func)
{
  cs_real_t  ohf = -ifc/(deq.meas * fabs(_dp3(pfq.unitv, deq.unitv)));

  for (int k = 0; k < 3; k++)
    grd_func[k] = ohf * pfq.unitv[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the gradient of a Lagrange hat function related to primal
 *          vertices in a p_{e,f, c} subvolume of a cell c where e is an edge
 *          belonging to the face f with vertices v and v'
 *
 * \param[in]       ubase     unit vector from xc to xv' (opposite)
 * \param[in]       udir      unit vector from xc to xv
 * \param[in]       len_vc    lenght of the segment [xc, xv]
 * \param[in]       deq       dual edge quantities
 * \param[in, out]  grd_func  gradient of the Lagrange shape function
 */
/*----------------------------------------------------------------------------*/

inline static void
_grd_v_lagrange_pefc(const cs_real_3_t     ubase,
                     const cs_real_3_t     udir,
                     cs_real_t             len_vc,
                     const cs_nvec3_t      deq,
                     cs_real_3_t           grd_func)
{
  cs_real_3_t  unormal;

  /* Normal direction to the plane in opposition to xv */
  cs_math_3_cross_product(ubase, deq.unitv, unormal);

  /* Height from this plane to the vertex */
  double  sgn = _dp3(udir, unormal);
  double  height = len_vc * sgn;
  assert(fabs(height) > cs_math_get_machine_epsilon()); /* Sanity check */
  double  over_height = 1/height;

  for (int k = 0; k < 3; k++)
    grd_func[k] = unormal[k]*over_height;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stiffness matrix on this cell using a Whitney
 *          Barycentric Subdivision (WBS) algo.
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      lm          pointer to a cs_cdo_locmesh_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_locmat_t *
_compute_wbs_stiffness(const cs_cdo_quantities_t   *quant,
                       const cs_cdo_locmesh_t      *lm,
                       const cs_real_3_t           *tensor,
                       cs_cdo_diff_t               *diff)
{
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + lm->n_max_vbyc;
  cs_real_t  *len_vc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + lm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*lm->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < lm->n_vc; v++)
    cs_math_3_length_unitv(xc, xyz + 3*lm->v_ids[v], len_vc + v, unit_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < lm->n_fc; f++) {

    const cs_nvec3_t  deq = lm->dedge[f];

    /* Compute the gradient of the lagrange function related to a cell
       in each p_{f,c} and the weights for each vertex related to this face */
    cs_compute_fwbs_q2(f, lm, grd_c, wvf, pefc_vol);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = lm->f2e_idx[f], jj = 0; i < lm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = lm->f2e_ids[i];
      const short int  v1 = lm->e2v_ids[2*e];
      const short int  v2 = lm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 */
      _grd_v_lagrange_pefc(unit_vc[v2], unit_vc[v1], len_vc[v1], deq, grd_v1);

      /* Gradient of the lagrange function related to a v2 */
      _grd_v_lagrange_pefc(unit_vc[v1], unit_vc[v2], len_vc[v2], deq, grd_v2);

      /* Gradient of the lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < sloc->n_ent; si++) {

        for (int k = 0; k < 3; k++)
          glv[si][k] = lm->wvc[si]*grd_c[k];

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
 * \param[in]      lm          pointer to a cs_cdo_locmesh_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_locmat_t *
_compute_vcb_stiffness(const cs_cdo_quantities_t   *quant,
                       const cs_cdo_locmesh_t      *lm,
                       const cs_real_3_t           *tensor,
                       cs_cdo_diff_t               *diff)
{
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg, matg_c;

  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + lm->n_max_vbyc;
  cs_real_t  *len_vc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + lm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*lm->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const int  msize = lm->n_vc + 1;
  const int  cc = msize*lm->n_vc + lm->n_vc;
  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;

  assert(sloc->n_ent == msize);

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < lm->n_vc; v++)
    cs_math_3_length_unitv(xc, xyz + 3*lm->v_ids[v], len_vc + v, unit_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < lm->n_fc; f++) {

    const cs_nvec3_t  deq = lm->dedge[f];

    /* Compute for the current face:
       - the gradient of the Lagrange function related xc in p_{f,c}
       - weights related to vertices
       - subvolume p_{ef,c} related to edges
    */
    const double  pfc_vol = cs_compute_fwbs_q3(f, lm, grd_c, wvf, pefc_vol);

    /* Compute the contribution to the entry A(c,c) */
    cs_math_33_3_product(tensor, grd_c, matg_c);
    sloc->val[cc] += pfc_vol * _dp3(grd_c, matg_c);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = lm->f2e_idx[f], jj = 0; i < lm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = lm->f2e_ids[i];
      const short int  v1 = lm->e2v_ids[2*e];
      const short int  v2 = lm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 */
      _grd_v_lagrange_pefc(unit_vc[v2], unit_vc[v1], len_vc[v1], deq, grd_v1);

      /* Gradient of the lagrange function related to a v2 */
      _grd_v_lagrange_pefc(unit_vc[v1], unit_vc[v2], len_vc[v2], deq, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < lm->n_vc; si++) {

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
         Be careful: sloc->n_ent = lm->n_vc + 1 */
      for (int si = 0; si < lm->n_vc; si++) {

        double  *mi = sloc->val + si*sloc->n_ent;

        /* Add v-c contribution */
        mi[lm->n_vc] += subvol * _dp3(matg_c, glv[si]);

        /* Add v-v contribution */
        cs_math_33_3_product(tensor, glv[si], matg);

        /* Loop on vertices v_j (j >= i) */
        for (int sj = si; sj < lm->n_vc; sj++)
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
 * \param[in]      lm      pointer to a cs_cdo_locmesh_t structure
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
                  const cs_cdo_locmesh_t   *lm,
                  double                    beta,
                  cs_cdo_diff_t            *diff,
                  cs_locmat_t              *ntrgrd)
{
  cs_real_3_t  lek;
  cs_real_3_t  *leg = diff->tmp_vect;
  cs_real_t  *v_coef = diff->tmp_real;

  const cs_real_t  over_vol_c = 1/lm->vol_c;

  for (short int v = 0; v < lm->n_vc; v++)
    v_coef[v] = 0;

  /* Loop on border face edges */
  for (int jf = lm->f2e_idx[f]; jf < lm->f2e_idx[f+1]; jf++) {

    short int  ei = lm->f2e_ids[jf];
    short int  vi1 = lm->e2v_ids[2*ei];
    short int  vi2 = lm->e2v_ids[2*ei+1];
    cs_quant_t  peq_i = lm->edge[ei];
    cs_nvec3_t  dfq_i = lm->dface[ei];

    const double  dp = _dp3(peq_i.unitv, dfq_i.unitv);
    const double  tmp_val = peq_i.meas * dfq_i.meas * dp; // 3*pec_vol
    const double  beta_pec_vol = 3. * beta/tmp_val;
    const double  coef_ei =  3. * beta/dp;

    /* Area of the triangle t_{e,f} defined by x_vi1, x_vi2 and x_f */
    const double  surf = cs_math_surftri(xyz + 3*lm->v_ids[vi1],
                                         xyz + 3*lm->v_ids[vi2],
                                         xf);

    /* Penalization term */
    v_coef[vi1] += 0.5*surf;
    v_coef[vi2] += 0.5*surf;

    /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
    for (short int v = 0; v < lm->n_vc; v++)
      leg[v][0] = leg[v][1] = leg[v][2] = 0;

    /* Term related to the flux reconstruction:
       Compute L_Ec(GRAD(p_j)) on t_{e,f} for each vertex j of the cell */
    for (short int ek = 0; ek < lm->n_ec; ek++) {

      const short int  shft = 2*ek;
      const short int  vj1 = lm->e2v_ids[shft];
      const short int  vj2 = lm->e2v_ids[shft+1];
      const short int  sgn_vj1 = lm->e2v_sgn[shft];
      const cs_nvec3_t  dfq_k = lm->dface[ek];

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

    for (short int v = 0; v < lm->n_vc; v++) {

      const double  contrib = _dp3(mn, leg[v]) * surf;
      ntrgrd->val[vi1*lm->n_vc + v] += contrib;
      ntrgrd->val[vi2*lm->n_vc + v] += contrib;

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
 * \param[in]      lm      pointer to a cs_cdo_locmesh_t structure
 * \param[in]      xyz     vertex coordinates
 * \param[in]      xc      cell center
 * \param[in]      mn      property tensor times the face normal
 * \param[in]      grd_c   gradient for the Lagrange function related to x_c
 * \param[in, out] diff    pointer to a builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_quant_wbs_algo(short int                 f,
                    const cs_cdo_locmesh_t   *lm,
                    const cs_real_t          *xyz,
                    const cs_real_t          *xc,
                    const cs_real_3_t         mn,
                    const cs_real_3_t         grd_c,
                    cs_cdo_diff_t            *diff)
{
  double  len;
  cs_real_3_t  un, cp, grd_f, grd_v1, grd_v2;

  /* Useful quantities are stored in diff->tmp_real and diff->tmp-vect */
  cs_real_3_t  *mng = diff->tmp_vect;
  cs_real_3_t  *unit_vc = diff->tmp_vect + lm->n_vc;
  cs_real_t  *wvf = diff->tmp_real;
  cs_real_t  *len_vc = diff->tmp_real + lm->n_vc;
  cs_real_t  *wtef = diff->tmp_real + 2*lm->n_vc;

  const cs_nvec3_t  deq = lm->dedge[f];
  const cs_quant_t  pfq = lm->face[f];
  const double  f_coef = 0.25/pfq.meas;

  for (short int v = 0; v < lm->n_vc; v++) {
    cs_math_3_length_unitv(xc, xyz + 3*lm->v_ids[v], len_vc + v, unit_vc[v]);
    wvf[v] = 0;
  }

  /* Compute a weight for each vertex of the current face */
  for (int ii = 0, i = lm->f2e_idx[f]; i < lm->f2e_idx[f+1]; i++, ii++) {

    const short int  e = lm->f2e_ids[i];
    const cs_quant_t  peq = lm->edge[e];
    const short int  v1 = lm->e2v_ids[2*e];
    const short int  v2 = lm->e2v_ids[2*e+1];

    cs_math_3_length_unitv(peq.center, pfq.center, &len, un);
    cs_math_3_cross_product(un, peq.unitv, cp);

    const double  tef2 = len * peq.meas * cs_math_3_norm(cp); // 2*|t_{e,f}|
    const double  wvf_contrib = tef2 * f_coef;
    const double  tef_coef = tef2 * cs_wbs_tef_weight;  // |t_{e,f}| / 12

    wvf[v1] += wvf_contrib;
    wvf[v2] += wvf_contrib;
    wtef[ii] = tef_coef;

    /* Gradient of the lagrange function related to v1 */
    _grd_v_lagrange_pefc(unit_vc[v2], unit_vc[v1], len_vc[v1], deq, grd_v1);

    /* Gradient of the lagrange function related to a v2 */
    _grd_v_lagrange_pefc(unit_vc[v1], unit_vc[v2], len_vc[v2], deq, grd_v2);

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
 * \param[in]      lm      pointer to a cs_cdo_locmesh_t structure
 * \param[in]      beta    value of the stabilizarion coef. related to reco.
 * \param[in, out] diff    pointer to a builder structure
 * \param[in, out] ntrgrd  local matrix related to the normal trace op.
 */
/*----------------------------------------------------------------------------*/

static void
_ntrgrd_wbs_algo(short int                 f,
                 double                    f_meas,
                 double                    mngc,
                 const cs_cdo_locmesh_t   *lm,
                 cs_cdo_diff_t            *diff,
                 cs_locmat_t              *ntrgrd)
{
  /* Useful quantities are stored in diff->tmp_real and diff->tmp-vect */
  cs_real_3_t  *mng = diff->tmp_vect;
  cs_real_t  *wvf = diff->tmp_real;

  double  sum_ef = 0.;
  for (int jf = lm->f2e_idx[f], ii = 0; jf < lm->f2e_idx[f+1]; jf++, ii++)
    sum_ef += mng[ii][2];

  for (short int vi = 0; vi < lm->n_vc; vi++) {
    if (wvf[vi] > 0) {

      double  *ntrgrd_i = ntrgrd->val + vi*lm->n_vc;

      for (short int vj = 0; vj < lm->n_vc; vj++) {

        /* Default contribution */
        double contrib = 0.25 * f_meas * wvf[vi] * lm->wvc[vj] * mngc; // 1ab

        if (wvf[vj] > 0) {

          contrib += wvf[vj] * wvf[vi] * sum_ef; // 2b

          for (int jf = lm->f2e_idx[f], ii = 0; jf < lm->f2e_idx[f+1];
               jf++, ii++) {

            const short int  e = lm->f2e_ids[jf];
            const short int  v1 = lm->e2v_ids[2*e];
            const short int  v2 = lm->e2v_ids[2*e+1];

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
  int  v_size = 0, s_size = 0;
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

  if (wnit || wsym) { // Weak enforcement of Dirichlet BCs
    v_size = connect->n_max_ebyc + connect->n_max_vbyc;
    s_size = connect->n_max_ebyc + connect->n_max_vbyc;
  }

  if (hwbs) {
    v_size = CS_MAX(v_size, 2*connect->n_max_vbyc);
    s_size = CS_MAX(s_size,
                    2*connect->n_max_vbyc + connect->f2e->info.stencil_max);
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

  if (hwbs || wnit || wsym) {
    BFT_FREE(diff->tmp_vect);
    BFT_FREE(diff->tmp_real);
  }

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
 * \brief   Define the local (cellwise) stiffness matrix
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      lm          cell-wise connectivity and quantitites
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdo_diffusion_build_local(const cs_cdo_quantities_t   *quant,
                             const cs_cdo_locmesh_t      *lm,
                             const cs_real_3_t           *tensor,
                             cs_cdo_diff_t               *diff)
{
  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  /* Initialize the local matrix */
  sloc->n_ent = lm->n_vc;
  for (int i = 0; i < lm->n_vc; i++)
    sloc->ids[i] = lm->v_ids[i];
  if (diff->scheme == CS_SPACE_SCHEME_CDOVCB) {
    sloc->n_ent += 1;
    sloc->ids[lm->n_vc] = lm->c_id;
  }

  for (int i = 0; i < sloc->n_ent*sloc->n_ent; i++)
    sloc->val[i] = 0;

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:

    /* Sanity check */
    assert(diff->scheme != CS_SPACE_SCHEME_CDOVCB);

    /* Set the diffusion tensor */
    if (lm->c_id == 0 || diff->is_uniform == false)
      cs_hodge_builder_set_tensor(diff->hb, tensor);

    /* Build a local discrete Hodge op. and return a local dense matrix */
    cs_hodge_build_local_stiffness(lm, diff->hb, sloc);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    if (diff->scheme == CS_SPACE_SCHEME_CDOVB)
      sloc = _compute_wbs_stiffness(quant, lm, tensor, diff);
    else if (diff->scheme == CS_SPACE_SCHEME_CDOVCB)
      sloc = _compute_vcb_stiffness(quant, lm, tensor, diff);
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
 * \brief   Compute the gradient of the conforming reconstruction in each
 *          p_{ef,c} tetrahedron
 *
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      lm          cell-wise connectivity and quantitites
 * \param[in]      pdi         cellwise values of the discrete potential
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 * \param[in, out] grd_lv_conf gradient of the conforming reconstruction
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_get_grd_lvconf(const cs_cdo_quantities_t   *quant,
                                const cs_cdo_locmesh_t      *lm,
                                const double                *pdi,
                                cs_cdo_diff_t               *diff,
                                double                      *grd_lv_conf)
{
  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;

  /* Sanity check */
  assert(h_algo == CS_PARAM_HODGE_ALGO_WBS);

  cs_real_3_t  grd_c, grd_v1, grd_v2;
  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_t  *len_vc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + lm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*lm->n_max_vbyc;

  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;

  /* Define the length and unit vector of the segment x_c --> x_v
     Compute the reconstruction at the cell center
   */
  double  p_c = 0;
  for (short int v = 0; v < lm->n_vc; v++) {
    cs_math_3_length_unitv(xc, xyz + 3*lm->v_ids[v], len_vc + v, unit_vc[v]);
    p_c += lm->wvc[v]*pdi[v];
  }

  /* Loop on cell faces */
  for (short int f = 0; f < lm->n_fc; f++) {

    cs_nvec3_t  deq = lm->dedge[f];

    /* Compute for the current face:
       - the gradient of the Lagrange function related xc in p_{f,c}
       - weights related to vertices
       - subvolume p_{ef,c} related to edges
    */
    cs_compute_fwbs_q2(f, lm, grd_c, wvf, pefc_vol);

    double  p_f = 0.;
    for (short int v = 0; v < lm->n_vc; v++)
      p_f += wvf[v]*pdi[v];

    const double delta_cf = p_c - p_f;

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = lm->f2e_idx[f], jj = 0; i < lm->f2e_idx[f+1]; i++, jj++) {

      const short int  e = lm->f2e_ids[i];
      const short int  v1 = lm->e2v_ids[2*e];
      const short int  v2 = lm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 */
      _grd_v_lagrange_pefc(unit_vc[v2], unit_vc[v1], len_vc[v1], deq, grd_v1);

      /* Gradient of the lagrange function related to a v2 */
      _grd_v_lagrange_pefc(unit_vc[v1], unit_vc[v2], len_vc[v2], deq, grd_v2);

      /* Gradient of the lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c}
      */
      for (int k = 0; k < 3; k++)
        grd_lv_conf[3*i+k] = delta_cf * grd_c[k] +
          (pdi[v1] - p_f) * grd_v1[k] + (pdi[v2] - p_f) * grd_v2[k];

    } // Loop on face edges

  } // Loop on cell face

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) "normal trace gradient" matrix taking
 *          into account Dirichlet BCs by a weak enforcement using Nitsche
 *          technique (symmetrized or not)
 *
 * \param[in]       f_id      face id (a border face attached to a Dir. BC)
 * \param[in]       quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]       lm        pointer to a cs_cdo_locmesh_t struct.
 * \param[in]       matpty    3x3 matrix related to the diffusion property
 * \param[in, out]  diff      auxiliary structure used to build the diff. term
 * \param[in, out]  ls        cell-wise structure sotring the local system
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_diffusion_weak_bc(cs_lnum_t                    f_id,
                         const cs_cdo_quantities_t   *quant,
                         const cs_cdo_locmesh_t      *lm,
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
  for (short int ii = 0; ii < lm->n_fc; ii++) {
    if (lm->f_ids[ii] == f_id) {
      f = ii;
      break;
    }
  }
  assert(f != -1);

  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;
  const cs_quant_t  pfq = lm->face[f];
  const double  f_coef = cs_weak_nitsche_pena_coef * pow(pfq.meas, -0.5) *
    diff->eig_ratio * diff->eig_max;

  assert(f_coef > 0); // Sanity check

  /* Initialize the local quantities */
  cs_locmat_t  *ntrgrd = diff->loc;

  ntrgrd->n_ent = lm->n_vc;
  for (short int v = 0; v < lm->n_vc; v++)
    ntrgrd->ids[v] = lm->v_ids[v];

  for (int i = 0; i < lm->n_vc*lm->n_vc; i++)
    ntrgrd->val[i] = 0;

  /* Compute the product: matpty*face unit normal */
  cs_real_3_t  mn;
  cs_math_33_3_product((const cs_real_t (*)[3])matpty, pfq.unitv, mn);

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    _ntrgrd_cost_algo(f, mn, pfq.center, xyz, lm, h_info.coef, diff, ntrgrd);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    {
      cs_real_3_t  grd_c;

      /* Compute useful quantities for the WBS algo. (stored in diff->tmp_*)
         1) Gradient of the lagrange function related to a cell in p_{f,c} */
      _grd_c_lagrange_pfc(lm->f_sgn[f], pfq, lm->dedge[f], grd_c);

      cs_real_t  mng_c = _dp3(mn, grd_c); // (pty_tensor * nu_f) . grd_c

      _get_quant_wbs_algo(f, lm, xyz, xc, mn, grd_c, diff);

      _ntrgrd_wbs_algo(f, pfq.meas, mng_c, lm, diff, ntrgrd);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  } // End of switch

  if (diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    cs_real_t  *mv = diff->tmp_real + lm->n_vc;

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
      double  *wtef = diff->tmp_real + 2*lm->n_vc;
      cs_locmat_t  *hloc = diff->transp;

      /* Build the border Hodge operator */
      for (int i = 0; i < lm->n_vc*lm->n_vc; i++)
        hloc->val[i] = 0;

      for (short int vi = 0; vi < lm->n_vc; vi++) {
        if (wvf[vi] > 0) {

          /* Diagonal contribution */
          double  *hi = hloc->val + vi*lm->n_vc;

          hi[vi] = pfq.meas * wvf[vi] * (wvf[vi]*0.5 + cs_math_onethird);

          /* Extra-diagonal contribution */
          for (short int vj = vi + 1; vj < lm->n_vc; vj++) {
            if (wvf[vj] > 0) {

              hi[vj] = 0.5 * wvf[vi] * wvf[vj] * pfq.meas;

              for (int jf = lm->f2e_idx[f], ii = 0; jf < lm->f2e_idx[f+1];
                   jf++, ii++) {

                const short int  e = lm->f2e_ids[jf];
                const short int  v1 = lm->e2v_ids[2*e];
                const short int  v2 = lm->e2v_ids[2*e+1];

                if (v1 == vi && v2 == vj)
                  hi[vj] += wtef[ii];
                else if (v1 == vj && v2 == vi)
                  hi[vj] += wtef[ii];

              } /* Loop on face edges */

              /* Matrix is symmetric (used only the upper part) */
              // hloc->val[vj*lm->n_vc + vi] = hi[vj];

            } // vj belongs to f
          } // Loop on cell vertices vj

        } // vi belongs to f
      } // Loop on cell vertices vi

      /* Add the border Hodge op. to the normal trace op.
         Update RHS whith H*p^{dir} */
      for (short int vi = 0; vi < lm->n_vc; vi++) {
        if (wvf[vi] > 0) {

          double  *ntrg_i = ntrgrd->val + vi*lm->n_vc;
          const double  *hi = hloc->val + vi*lm->n_vc;
          const double pii_coef = f_coef * hi[vi];

          // Set the penalty diagonal coefficient
          ntrg_i[vi] += pii_coef;
          ls->rhs[vi] += pii_coef * ls->dir_bc[vi];

          for (short int vj = vi+1; vj < lm->n_vc; vj++) {
            if (wvf[vj] > 0) {

              const double  pij_coef = f_coef * hi[vj];

              ntrg_i[vj] += pij_coef;
              ntrgrd->val[vi + vj*lm->n_vc] += pij_coef;
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
  cs_locmat_dump(lm->c_id, ntrgrd);
#endif

  /* Add contribution to the linear system */
  cs_locmat_add(ls->mat, ntrgrd);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
