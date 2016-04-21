/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based schemes
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

#include "cs_math.h"
#include "cs_property.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_diffusion.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cdovb_diffusion.c

  \brief  Build discrete stiffness matrices and handled boundary conditions
          diffusion term in CDO vertex-based schemes.

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_CDOVB_DIFFUSION_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/* Stiffness matrix builder */
struct _cs_cdovb_diff_t {

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
static const cs_real_t  cs_weak_nitsche_pena_coef = 500;

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
 * \brief  Compute for each face a weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Set also the local and local numbering of the vertices of this face
 *
 * \param[in]      _f         id of the face in the cell-wise numbering
 * \param[in]      pfq        geometric quantities related to the current face
 * \param[in]      deq        geometric quantities related to its dual edge
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 */
/*----------------------------------------------------------------------------*/

static void
_compute_wvf_pefcvol(short int                   _f,
                     const cs_quant_t            pfq,
                     const cs_nvec3_t            deq,
                     const cs_cdo_locmesh_t     *lm,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol)
{
  double  len;
  cs_real_3_t  un, cp;

  const double  f_coef = 0.25/pfq.meas;
  const double  h_coef = 1/6.*_dp3(pfq.unitv, deq.unitv)*deq.meas;

  assert(h_coef > 0);

  /* Reset weights */
  for (int ii = 0; ii < lm->n_vc; ii++) wvf[ii] = 0;

  /* Compute a weight for each vertex of the current face */
  for (int ii = 0, i = lm->f2e_idx[_f]; i < lm->f2e_idx[_f+1]; i++, ii++) {

    const short int  _e = lm->f2e_ids[i];
    const cs_quant_t  peq = lm->edge[_e];
    const short int  _v1 = lm->e2v_ids[2*_e];
    const short int  _v2 = lm->e2v_ids[2*_e+1];

    cs_math_3_length_unitv(peq.center, pfq.center, &len, un);
    cs_math_3_cross_product(un, peq.unitv, cp);

    double  area = len * peq.meas * cs_math_3_norm(cp); // two times the area
    double  contrib = area * f_coef;

    pefc_vol[ii] = area * h_coef;
    wvf[_v1] += contrib;
    wvf[_v2] += contrib;

  } /* End of loop on face edges */

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
                       cs_cdovb_diff_t             *diff)
{
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + lm->n_max_vbyc;
  cs_real_t  *len_vc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + lm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*lm->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < lm->n_vc; v++) {

    const cs_real_t  *xv = quant->vtx_coord + 3*lm->v_ids[v];

    cs_math_3_length_unitv(xc, xv, len_vc + v, unit_vc[v]);

  } // Loop on cell vertices

  /* Loop on cell faces */
  for (short int f = 0; f < lm->n_fc; f++) {

    cs_quant_t  pfq = lm->face[f];
    cs_nvec3_t  deq = lm->dedge[f];

    /* Gradient of the lagrange function related to a cell in each p_{f,c} */
    _grd_c_lagrange_pfc(lm->f_sgn[f], pfq, deq, grd_c);

    /* Weights related to each vertex attached to this face */
    _compute_wvf_pefcvol(f, pfq, deq, lm, wvf, pefc_vol);

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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure used to build the stiffness matrix
 *
 * \param[in] connect      pointer to a cs_cdo_connect_t struct.
 * \param[in] is_uniform   diffusion tensor is uniform ? (true or false)
 * \param[in] h_info       cs_param_hodge_t struct.
 * \param[in] bc_enforce   type of boundary enforcement for Dirichlet values
 *
 * \return a pointer to a new allocated cs_cdovb_diffusion_builder_t struc.
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_diff_t *
cs_cdovb_diffusion_builder_init(const cs_cdo_connect_t       *connect,
                                bool                          is_uniform,
                                const cs_param_hodge_t        h_info,
                                const cs_param_bc_enforce_t   bc_enforce)
{
  int  v_size = 0, s_size = 0;
  cs_cdovb_diff_t  *diff = NULL;

  BFT_MALLOC(diff, 1, cs_cdovb_diff_t);

  diff->is_uniform = is_uniform;

  /* Copy the data related to a discrete Hodge operator */
  diff->h_info.type    = h_info.type;
  diff->h_info.inv_pty = h_info.inv_pty;
  diff->h_info.algo    = h_info.algo;
  diff->h_info.coef    = h_info.coef;

  diff->enforce = bc_enforce;
  bool  wnit = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) ? true : false;
  bool  wsym = (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;
  bool  hwbs = (h_info.algo == CS_PARAM_HODGE_ALGO_WBS) ? true : false;

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
  if (wsym)
    diff->transp = cs_locmat_create(connect->n_max_vbyc);

  /* Allocate the local stiffness matrix */
  diff->loc = cs_locmat_create(connect->n_max_vbyc);

  return diff;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_cdovb_diff_t structure
 *
 * \param[in, out ] diff   pointer to a cs_cdovb_diff_t struc.
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_diff_t *
cs_cdovb_diffusion_builder_free(cs_cdovb_diff_t   *diff)
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

  if (wsym)
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
 * \param[in]  diff   pointer to a cs_cdovb_diff_t structure
 *
 * \return  a pointer to a cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_cdovb_diffusion_get_hodge_builder(cs_cdovb_diff_t   *diff)
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
cs_cdovb_diffusion_build_local(const cs_cdo_quantities_t   *quant,
                               const cs_cdo_locmesh_t      *lm,
                               const cs_real_3_t           *tensor,
                               cs_cdovb_diff_t             *diff)
{
  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  /* Initialize the local matrix */
  sloc->n_ent = lm->n_vc;
  for (int i = 0; i < sloc->n_ent; i++)
    sloc->ids[i] = lm->v_ids[i];
  for (int i = 0; i < sloc->n_ent*sloc->n_ent; i++)
    sloc->val[i] = 0;

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:

    /* Set the diffusion tensor */
    if (lm->c_id == 0 || diff->is_uniform == false)
      cs_hodge_builder_set_tensor(diff->hb, tensor);

    /* Build a local discrete Hodge op. and return a local dense matrix */
    cs_hodge_build_local_stiffness(lm, diff->hb, sloc);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    sloc = _compute_wbs_stiffness(quant, lm, tensor, diff);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to build the local stiffness matrix.");

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
cs_cdovb_diffusion_get_grd_lvconf(const cs_cdo_quantities_t   *quant,
                                  const cs_cdo_locmesh_t      *lm,
                                  const double                *pdi,
                                  cs_cdovb_diff_t             *diff,
                                  double                      *grd_lv_conf)
{
  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;

  /* Sanity check */
  assert(h_algo == CS_PARAM_HODGE_ALGO_WBS);

  cs_real_3_t  grd_c, grd_v1, grd_v2;

  double  p_c = 0;
  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_t  *len_vc  = diff->tmp_real;
  cs_real_t  *wvf     = diff->tmp_real + lm->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 2*lm->n_max_vbyc;

  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;

  /* Define the length and unit vector of the segment x_c --> x_v
     Compute the reconstruction at the cell center
   */
  for (short int v = 0; v < lm->n_vc; v++) {

    const cs_real_t  *xv = quant->vtx_coord + 3*lm->v_ids[v];

    cs_math_3_length_unitv(xc, xv, len_vc + v, unit_vc[v]);
    p_c += lm->wvc[v]*pdi[v];

  } // Loop on cell vertices

  /* Loop on cell faces */
  for (short int f = 0; f < lm->n_fc; f++) {

    cs_quant_t  pfq = lm->face[f];
    cs_nvec3_t  deq = lm->dedge[f];

    /* Gradient of the lagrange function related to a cell in each p_{f,c} */
    _grd_c_lagrange_pfc(lm->f_sgn[f], pfq, deq, grd_c);

    /* Weights related to each vertex attached to this face */
    _compute_wvf_pefcvol(f, pfq, deq, lm, wvf, pefc_vol);

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
cs_cdovb_diffusion_weak_bc(cs_lnum_t                    f_id,
                           const cs_cdo_quantities_t   *quant,
                           cs_cdo_locmesh_t            *lm,
                           const cs_real_t              matpty[3][3],
                           cs_cdovb_diff_t             *diff,
                           cs_cdo_locsys_t             *ls)
{
  cs_real_3_t  mn;

  const double  *xyz = quant->vtx_coord;
  const cs_param_hodge_t  h_info = diff->h_info;
  const cs_param_hodge_algo_t  h_algo = h_info.algo;

  /* Set the diffusion tensor */
  if (diff->eig_ratio < 0 || diff->is_uniform == false)
    cs_math_33_eigen((const cs_real_t (*)[3])matpty,
                     &(diff->eig_ratio),
                     &(diff->eig_max));

  /* Set the local face id */
  short int _f = -1; // not set
  for (int ii = 0; ii < lm->n_fc; ii++) {
    if (lm->f_ids[ii] == f_id) {
      _f = ii;
      break;
    }
  }
  assert(_f != -1);

  const cs_quant_t  pfq = lm->face[_f];
  const double  f_coef = pow(pfq.meas, -0.5) * diff->eig_ratio * diff->eig_max;

  assert(f_coef > 0); // Sanity check

  cs_real_3_t  *leg = diff->tmp_vect;
  cs_real_t  *vcoef = diff->tmp_real;
  cs_locmat_t  *ntrgrd = diff->loc;

  ntrgrd->n_ent = lm->n_vc;
  for (int _v = 0; _v < lm->n_vc; _v++) {
    ntrgrd->ids[_v] = lm->v_ids[_v];
    vcoef[_v] = 0;
    ls->rhs[_v] = 0;
  }

  /* Initialize the local matrix */
  for (int i = 0; i < lm->n_vc*lm->n_vc; i++)
    ntrgrd->val[i] = 0;

  /* Compute the product: matpty*face unit normal */
  cs_math_33_3_product((const cs_real_t (*)[3])matpty, pfq.unitv, mn);

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    {
      cs_real_3_t  lek;

      const double  beta = h_info.coef;
      const cs_real_t  over_vol_c = 1/lm->vol_c;

      /* Loop on border face edges */
      for (int jf = lm->f2e_idx[_f]; jf < lm->f2e_idx[_f+1]; jf++) {

        short int  _ei = lm->f2e_ids[jf];
        short int  _vi1 = lm->e2v_ids[2*_ei];
        short int  _vi2 = lm->e2v_ids[2*_ei+1];
        cs_quant_t  peq_i = lm->edge[_ei];
        cs_nvec3_t  dfq_i = lm->dface[_ei];

        const double  dp = _dp3(peq_i.unitv, dfq_i.unitv);
        const double  tmp_val = peq_i.meas * dfq_i.meas * dp; // 3*pec_vol
        const double  beta_pec_vol = 3.*beta/tmp_val;
        const double  coef_ei =  3 * beta/dp;
        const double  surf = cs_math_surftri(xyz + 3*lm->v_ids[_vi1],
                                             peq_i.center,
                                             pfq.center);

        /* Penalization term */
        vcoef[_vi1] += surf*f_coef;
        vcoef[_vi2] += surf*f_coef;

        /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
        for (int _v = 0; _v < lm->n_vc; _v++)
          leg[_v][0] = leg[_v][1] = leg[_v][2] = 0;

        /* Term related to the flux reconstruction:
           Compute L_Ec(GRAD(p_j)) on t_{_e,f} for each vertex j of the cell */
        for (short int _ek = 0; _ek < lm->n_ec; _ek++) {

          short int  shft = 2*_ek;
          short int  _vj1 = lm->e2v_ids[shft];
          short int  _vj2 = lm->e2v_ids[shft+1];
          short int  sgn_vj1 = lm->e2v_sgn[shft];
          short int  sgn_vj2 = lm->e2v_sgn[shft+1];
          cs_nvec3_t  dfq_k = lm->dface[_ek];

          /* Compute l_(ek,c)|p(_ei,f,c) */
          const double  eik_part = coef_ei * _dp3(dfq_k.unitv, peq_i.unitv);
          const double  coef_mult = dfq_k.meas * over_vol_c;

          for (int k = 0; k < 3; k++)
            lek[k] = coef_mult * (dfq_k.unitv[k] - eik_part * dfq_i.unitv[k]);

          if (_ek == _ei)
            for (int k = 0; k < 3; k++)
              lek[k] += dfq_k.meas * dfq_k.unitv[k] * beta_pec_vol;

          for (int k = 0; k < 3; k++) {
            leg[_vj1][k] += sgn_vj1 * lek[k];
            leg[_vj2][k] += sgn_vj2 * lek[k];
          }

        } // Loop on cell edges

        for (short int v = 0; v < lm->n_vc; v++) {

          const double  contrib = _dp3(mn, leg[v]) * surf;
          ntrgrd->val[_vi1*lm->n_vc + v] += contrib;
          ntrgrd->val[_vi2*lm->n_vc + v] += contrib;

        }

      } // border face edges

    }
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    bft_error(__FILE__, __LINE__, 0,
              " Nitsche enforcement of boundary conditions is not yet"
              " compatible with WBS algorithm.");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  } // End of switch

  switch (diff->enforce) {
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    {
      cs_real_t  *_mvec = diff->tmp_real + lm->n_vc;

      /* ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
      cs_locmat_add_transpose(ntrgrd, diff->transp);

      /* Modify RHS according to the add of transp */
      cs_locmat_matvec(diff->transp, ls->dir_bc, _mvec);
      for (int j = 0; j < ntrgrd->n_ent; j++)
        ls->rhs[j] += _mvec[j];

      // No break to continue with Nitsche treatment
    }

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:

    for (int j = 0; j < ntrgrd->n_ent; j++) {
      if (vcoef[j] > 0) {

        const double pena_coef = cs_weak_nitsche_pena_coef * vcoef[j];

        ntrgrd->val[j*ntrgrd->n_ent + j] += pena_coef; // pena. diag. entry
        ls->rhs[j] += pena_coef * ls->dir_bc[j];
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " Invalid type of BC enforcement");

  }

#if CS_CDOVB_DIFFUSION_DBG > 1
  bft_printf(">> Local weak bc matrix (f_id: %d)", f_id);
  cs_locmat_dump(lm->c_id, ntrgrd);
#endif

  /* Add contribution to the linear system */
  cs_locmat_add(ls->mat, ntrgrd);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
