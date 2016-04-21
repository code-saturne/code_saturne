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
#include <limits.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_hodge.h"
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
  bool                  is_uniform;  /* Diffusion tensor is uniform ? */
  cs_param_hodge_t      h_info;      /* Set of parameters related to a discrete
                                        Hodge operator */
  cs_hodge_builder_t   *hb;

  /* Compact way to stored the edge --> vertices connect. related to a cell */
  int           n_bits;         // number of bits in a block
  int           n_blocks;       // number of blocks in a mask
  cs_flag_t    *emsk;           // list of masks to store the connectivity

  /* Temporary buffers */
  cs_real_3_t  *tmp_vect;  // set of local vectors
  cs_real_t    *tmp_real;  // set of local arrays of double
  short int    *tmp_ids;   // set of local ids related to a cell

  /* Local matrix (stiffness or normal trace gradient) */
  cs_locmat_t   *loc;

};

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Encode E_c cap E_v sets for all vertices of a cell using a mask
 *          of bits
 *
 * \param[in]       c_id      cell id
 * \param[in]       connect   point to a cs_cdo_connect_t struct.
 * \param[in]       vtag      tag on vertices to store local ids
 * \param[in, out]  diff      pointer to a cs_cdovb_diff_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_encode_edge_masks(cs_lnum_t                  c_id,
                   const cs_cdo_connect_t    *connect,
                   const cs_lnum_t            vtag[],
                   cs_cdovb_diff_t           *diff)
{
  short int  ie;
  cs_lnum_t  i;

  cs_flag_t  *masks = diff->emsk;

  const int  n_bits = diff->n_bits;
  const cs_connect_index_t  *c2e = connect->c2e;

  /* Reset masks */
  for (i = 0; i < diff->n_blocks*connect->n_max_vbyc; i++)
    masks[i] = 0;

  /* Loop on cell edges */
  for (ie = 0, i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; ie++, i++) {

    cs_lnum_t  eshft = 2*c2e->ids[i];
    short int  vi = vtag[connect->e2v->col_id[eshft]];
    short int  vj = vtag[connect->e2v->col_id[eshft+1]];

    /* Sanity checks */
    assert(vi > -1 && vi < connect->n_max_vbyc);
    assert(vj > -1 && vj < connect->n_max_vbyc);

    /* Encode this edge in the set E_v1,c and E_v2,c */
    short int  b = ie/n_bits, s = ie%n_bits;
    masks[vi*diff->n_blocks + b] |= (1 << s);
    masks[vj*diff->n_blocks + b] |= (1 << s);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute and store the unit vector tangential to x_c --> x_v.
 *          Store also the local weight related to each vertex wvc
 *
 * \param[in]      c_id        current cell id
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 */
/*----------------------------------------------------------------------------*/

static void
_define_wbs_builder(cs_lnum_t                    c_id,
                    const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant,
                    cs_cdovb_diff_t             *diff)
{
  cs_lnum_t  i, _i;

  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_t  *wvc = diff->tmp_real;
  cs_real_t  *len_vc = diff->tmp_real + connect->n_max_vbyc;

  const cs_real_t  over_cell_vol = 1/quant->cell_vol[c_id];
  const cs_real_t  *xc = quant->cell_centers + 3*c_id;
  const cs_connect_index_t  *c2v = connect->c2v;

  /* Loop on cell vertices */
  for (i = c2v->idx[c_id], _i = 0; i < c2v->idx[c_id+1]; i++, _i++) {

    cs_lnum_t  v_id = c2v->ids[i];
    const cs_real_t  *xv = quant->vtx_coord + 3*v_id;

    wvc[_i] = over_cell_vol * quant->dcell_vol[i];

    /* Define the segment x_c --> x_v */
    cs_math_3_length_unitv(xc, xv, len_vc + _i, unit_vc[_i]);

  } // Loop on cell vertices

}

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

static void
_grad_lagrange_cell_pfc(short int             ifc,
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

static void
_grad_lagrange_vtx_pefc(const cs_real_3_t     ubase,
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
 * \param[in]      f_id       id of the face
 * \param[in]      pfq        geometric quantities related to face f_id
 * \param[in]      deq        geometric quantities related to the dual edge
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantites_t structure
 * \param[in]      loc_ids    local vertex numbering on a cell
 * \param[in, out] wvf        pointer to an array storing the weight/vertex
 * \param[in, out] pefc_vol   pointer to an array storing the volume of pefc
 */
/*----------------------------------------------------------------------------*/

static void
_compute_wvf_pefcvol(cs_lnum_t                   f_id,
                     const cs_quant_t            pfq,
                     const cs_nvec3_t            deq,
                     const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     const cs_lnum_t            *loc_ids,
                     cs_real_t                  *wvf,
                     cs_real_t                  *pefc_vol)
{
  cs_lnum_t  i, j;
  double  len;
  cs_real_3_t  un, cp;

  const double  f_coef = 0.25/pfq.meas;
  const double  h_coef = 1/6.*_dp3(pfq.unitv, deq.unitv)*deq.meas;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  assert(h_coef > 0);

  /* Reset weights */
  for (i = 0; i < connect->n_max_vbyc; i++) wvf[i] = 0;

  /* Compute a weight for each vertex of the current face */
  for (j = 0, i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++, j++) {

    const cs_lnum_t  e_id = f2e->col_id[i];
    const cs_quant_t  eq = quant->edge[e_id];

    const cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
    const cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];
    const short int  _v1 = loc_ids[v1_id];
    const short int  _v2 = loc_ids[v2_id];

    cs_math_3_length_unitv(eq.center, pfq.center, &len, un);
    cs_math_3_cross_product(un, eq.unitv, cp);

    double  area = len * eq.meas * cs_math_3_norm(cp); // two times the area
    double  contrib = area * f_coef;

    pefc_vol[j] = area * h_coef;
    wvf[_v1] += contrib;
    wvf[_v2] += contrib;

  } /* End of loop on face edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the stiffness matrix on this cell using a Whitney
 *          Barycentric Subdivision (WBS) algo.
 *
 * \param[in]      c_id        current cell id
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      loc_ids     local vertex numbering in a cell
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_locmat_t *
_compute_wbs_stiffness(cs_lnum_t                    c_id,
                       const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_lnum_t             *loc_ids,
                       const cs_real_3_t           *tensor,
                       cs_cdovb_diff_t             *diff)
{
  cs_lnum_t  i, j, _j, k, ii, jj, ishft;
  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *unit_vc = diff->tmp_vect;
  cs_real_3_t  *glv = diff->tmp_vect + connect->n_max_vbyc;
  cs_real_t  *wvc     = diff->tmp_real;
  cs_real_t  *len_vc  = diff->tmp_real +   connect->n_max_vbyc;
  cs_real_t  *wvf     = diff->tmp_real + 2*connect->n_max_vbyc;
  cs_real_t  *pefc_vol = diff->tmp_real + 3*connect->n_max_vbyc;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  /* Loop on cell faces */
  for (i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

    short int  ifc = c2f->sgn[i];
    cs_lnum_t  f_id = c2f->col_id[i];
    cs_quant_t  pfq = quant->face[f_id];
    cs_nvec3_t  deq = quant->dedge[i];

    /* Gradient of the lagrange function related to a cell in each p_{f,c} */
    _grad_lagrange_cell_pfc(ifc, pfq, deq, grd_c);

    /* Weights related to each vertex attached to this face */
    _compute_wvf_pefcvol(f_id, pfq, deq, connect, quant, loc_ids, wvf, pefc_vol);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (j = f2e->idx[f_id], _j = 0; j < f2e->idx[f_id+1]; j++, _j++) {

      const cs_lnum_t  e_id = f2e->col_id[j];
      const cs_lnum_t  eshf = 2*e_id;
      const double  subvol = pefc_vol[_j];

      /* First vertex */
      const cs_lnum_t  v1_id = e2v->col_id[eshf];
      const short int  _v1 = loc_ids[v1_id];

      /* Second vertex */
      const cs_lnum_t  v2_id = e2v->col_id[eshf+1];
      const short int  _v2 = loc_ids[v2_id];

      /* Gradient of the lagrange function related to v1 */
      _grad_lagrange_vtx_pefc(unit_vc[_v2], unit_vc[_v1], len_vc[_v1], deq,
                              grd_v1);

      /* Gradient of the lagrange function related to a v2 */
      _grad_lagrange_vtx_pefc(unit_vc[_v1], unit_vc[_v2], len_vc[_v2], deq,
                              grd_v2);

      /* Gradient of the lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (ii = 0; ii < sloc->n_ent; ii++) {

        for (k = 0; k < 3; k++)
          glv[ii][k] = wvc[ii]*grd_c[k];

        if (wvf[ii] > 0) // Face contrib.
          for (k = 0; k < 3; k++)
            glv[ii][k] += wvf[ii]*grd_f[k];

        if (ii == _v1) // Vertex 1 contrib
          for (k = 0; k < 3; k++)
            glv[ii][k] += grd_v1[k];

        if (ii == _v2) // Vertex 2 contrib
          for (k = 0; k < 3; k++)
            glv[ii][k] += grd_v2[k];

      } // Loop on cell vertices

      /* Build the upper right part */
      for (ii = 0; ii < sloc->n_ent; ii++) {

        ishft = ii*sloc->n_ent;
        cs_math_33_3_product(tensor, glv[ii], matg);

        /* Diagonal contribution */
        sloc->mat[ishft+ii] += subvol * _dp3(matg, glv[ii]);

        for (jj = ii+1; jj < sloc->n_ent; jj++) {

          sloc->mat[ishft+jj] += subvol * _dp3(matg, glv[jj]);

        } /* Loop on vertices v_j (j > i) */

      } /* Loop on vertices v_i */

    }

  } // Loop on cell faces

  /* Matrix is symmetric by construction */
  for (ii = 0; ii < sloc->n_ent; ii++) {
    ishft = ii*sloc->n_ent;
    for (jj = ii+1; jj < sloc->n_ent; jj++)
      sloc->mat[jj*sloc->n_ent+ii] = sloc->mat[ishft+jj];
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
  cs_lnum_t  i;

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

  int  vect_size = 0;
  int  vol_size = 0;

  diff->tmp_ids = NULL;
  diff->tmp_real = NULL;
  diff->tmp_vect = NULL;
  diff->emsk = NULL;
  diff->hb = NULL;
  diff->n_bits = 0;
  diff->n_blocks = 0;

  if (wnit || wsym) {

    vect_size = connect->n_max_ebyc + connect->n_max_vbyc;
    vol_size = connect->n_max_ebyc;

    cs_lnum_t  n_ent = CS_MAX(connect->e_info->n_elts, connect->f_info->n_elts);
    BFT_MALLOC(diff->tmp_ids, n_ent, short int);
    for (i = 0; i < n_ent; i++)
      diff->tmp_ids[i] = -1;

  } // Weakly enforcement of Dirichlet BCs

  if (hwbs) {

    vect_size = CS_MAX(vect_size, 2*connect->n_max_vbyc);
    vol_size = CS_MAX(vol_size,
                      3*connect->n_max_vbyc + connect->f2e->info.stencil_max);

  }
  else {

    /* Define a builder for the related discrete Hodge operator */
    diff->hb = cs_hodge_builder_init(connect, h_info);

    /* Initialize mask features */
    diff->n_bits = sizeof(cs_flag_t)*CHAR_BIT;
    diff->n_blocks = connect->n_max_ebyc/diff->n_bits;
    if (connect->n_max_ebyc % diff->n_bits != 0)
      diff->n_blocks += 1;

    BFT_MALLOC(diff->emsk, diff->n_blocks*connect->n_max_vbyc, cs_flag_t);
    for (i = 0; i < diff->n_blocks*connect->n_max_vbyc; i++)
      diff->emsk[i] = 0;

  }

  if (vect_size > 0) BFT_MALLOC(diff->tmp_vect, vect_size, cs_real_3_t);
  if (vol_size > 0) BFT_MALLOC(diff->tmp_real, vol_size, cs_real_t);

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
    if (wnit || wsym)
      BFT_FREE(diff->tmp_ids);
  }

  if (!hwbs) {
    BFT_FREE(diff->emsk);
    diff->hb = cs_hodge_builder_free(diff->hb);
  }

  /* Local stiffness matrix */
  diff->loc = cs_locmat_free(diff->loc);

  BFT_FREE(diff);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) stiffness matrix
 *
 * \param[in]      c_id        current cell id
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      vtag        pointer to a cs_cdovb_scaleq_t struct.
 * \param[in]      tensor      3x3 matrix attached to the diffusion property
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_diffusion_build_local(cs_lnum_t                    c_id,
                               const cs_cdo_connect_t      *connect,
                               const cs_cdo_quantities_t   *quant,
                               const cs_lnum_t             *vtag,
                               const cs_real_3_t           *tensor,
                               cs_cdovb_diff_t             *diff)
{
  cs_lnum_t  i, n_ent;

  cs_locmat_t  *sloc = diff->loc; // Local stiffness matrix to build

  const cs_param_hodge_algo_t  h_algo = diff->h_info.algo;
  const cs_connect_index_t  *c2v = connect->c2v;

  /* Initialize the local matrix */
  for (i = c2v->idx[c_id], n_ent = 0; i < c2v->idx[c_id+1]; i++, n_ent++)
    sloc->ids[n_ent] = c2v->ids[i];
  sloc->n_ent = n_ent;
  for (i = 0; i < n_ent*n_ent; i++)
    sloc->mat[i] = 0;

  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    {
      short int  ek, el, vi, sgn_ik, sgn_jl;

      /* Set the diffusion tensor */
      if (c_id == 0 || diff->is_uniform == false)
        cs_hodge_builder_set_tensor(diff->hb, tensor);

      /* Build a local discrete Hodge op. and return a local dense matrix */
      const cs_locmat_t  *hloc = cs_hodge_build_local(c_id,
                                                      connect, quant,
                                                      diff->hb);

      /* Encode edge sets related to each vertex belonging to this cell */
      _encode_edge_masks(c_id, connect, vtag, diff);

      /* Build local stiffness matrix */
      for (vi = 0; vi < sloc->n_ent; vi++) {

        const short int pos_i = vi*sloc->n_ent;

        for (ek = 0; ek < hloc->n_ent; ek++) {  /* Find edges attached to vi */

          const short int b = ek/diff->n_bits, r = ek % diff->n_bits;

          if (diff->emsk[vi*diff->n_blocks+b] & (1 << r)) { // ek in E_i

            const cs_lnum_t  ek_shft = 2*hloc->ids[ek];
            const int  pos_ek = ek*hloc->n_ent;

            if (connect->e2v->col_id[ek_shft] == sloc->ids[vi])
              sgn_ik = connect->e2v->sgn[ek_shft];
            else {
              assert(connect->e2v->col_id[ek_shft+1] == sloc->ids[vi]);
              sgn_ik = connect->e2v->sgn[ek_shft+1];
            }

            for (el = 0; el < hloc->n_ent; el++) { /* Loop on cell edges */

              cs_lnum_t  el_shft = 2*hloc->ids[el];
              double  val = hloc->mat[pos_ek+el]*sgn_ik;
              cs_lnum_t  v1_id = connect->e2v->col_id[el_shft];
              short int  vj1 = vtag[v1_id];

              sgn_jl = connect->e2v->sgn[el_shft];
              sloc->mat[pos_i+vj1] += val*sgn_jl;

              cs_lnum_t  v2_id = connect->e2v->col_id[el_shft+1];
              short int  vj2 = vtag[v2_id];

              sgn_jl = connect->e2v->sgn[el_shft+1];
              sloc->mat[pos_i+vj2] += val*sgn_jl;

            } /* Loop on cell edges */

          } /* ek in E_i */

        } /* Find edges attached to _vi */

      } /* Loop on vertices vi */

    }
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    _define_wbs_builder(c_id, connect, quant, diff);
    sloc = _compute_wbs_stiffness(c_id, connect, quant, vtag, tensor, diff);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to build the local stiffness matrix.");

  } // End of switch

  return sloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the local (cellwise) "normal trace gradient" matrix
 *          This local matrix is used in Nitsche method to weakly penalized
 *          Dirichlet boundary conditions.
 *
 * \param[in]      c_id        cell id
 * \param[in]      f_id        face id (a border face attached to a Dir. BC)
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      matpty      3x3 matrix related to the diffusion property
 * \param[in]      eig_ratio   eigenvalue_max/eigenvalue_min
 * \param[in]      eig_max     eigenvalue with maximal value
 * \param[in, out] loc_v_ids   store local vertex ids
 * \param[in, out] v_coef      store local contribution on each border vertex
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 *
 * \return a pointer to a local "normal trace gradient" matrix
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_cdovb_diffusion_ntrgrd_build(cs_lnum_t                    c_id,
                                cs_lnum_t                    f_id,
                                const cs_cdo_connect_t      *connect,
                                const cs_cdo_quantities_t   *quant,
                                const cs_real_t              matpty[3][3],
                                cs_real_t                    eig_ratio,
                                cs_real_t                    eig_max,
                                cs_lnum_t                   *loc_v_ids,
                                cs_real_t                   *v_coef,
                                cs_cdovb_diff_t             *diff)
{
  cs_lnum_t  j, je, jv, k, _v, _e, v_id;
  cs_real_3_t  mn;

  cs_locmat_t  *ntrgrd = diff->loc; // Local matrix to build
  cs_real_3_t  *dfv = diff->tmp_vect;
  cs_real_3_t  *leg = diff->tmp_vect + connect->n_max_ebyc;
  short int  *loc_e_ids = diff->tmp_ids;

  const cs_param_hodge_t  h_info = diff->h_info;
  const cs_param_hodge_algo_t  h_algo = h_info.algo;
  const cs_quant_t  pfq = quant->face[f_id];
  const double  f_coef = pow(pfq.meas, -0.5) * eig_ratio * eig_max;
  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_connect_index_t  *c2v = connect->c2v;

  /*  Initialize buffers related to vertices.
      Define an id local to this cell for each vertex */
  const cs_lnum_t  start_v = c2v->idx[c_id];
  const cs_lnum_t  end_v = c2v->idx[c_id+1];
  const int  n_vbyc = end_v - start_v;
  const cs_lnum_t  start_e = c2e->idx[c_id];
  const cs_lnum_t  end_e = c2e->idx[c_id+1];

  for (jv = start_v, _v = 0; jv < end_v; jv++, _v++) {
    v_id = c2v->ids[jv];
    loc_v_ids[v_id] = _v;
    ntrgrd->ids[_v] = v_id;
  }
  ntrgrd->n_ent = n_vbyc;

  /* Sanity check */
  assert(ntrgrd->n_ent <= connect->n_max_vbyc);

  /* Initialize the local matrix */
  for (int  i = 0; i < n_vbyc*n_vbyc; i++)
    ntrgrd->mat[i] = 0;

  /* Compute the product: matpty*face unit normal */
  cs_math_33_3_product((const cs_real_t (*)[3])matpty, pfq.unitv, mn);

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
  case CS_PARAM_HODGE_ALGO_VORONOI:
    {
      short int _ek;
      cs_real_3_t  lek;

      const double  beta = h_info.coef;
      const cs_real_t  over_c_vol = 1/quant->cell_vol[c_id];
      const cs_sla_matrix_t  *e2v = connect->e2v;
      const cs_sla_matrix_t  *f2e = connect->f2e;

      /* Store the local id and useful quantities for each edge */
      for (je = start_e, _e = 0; je < end_e; je++, _e++) {

        cs_lnum_t  e_id = c2e->ids[je];
        cs_dface_t  dfq = quant->dface[je];

        loc_e_ids[e_id] = _e;
        for (k = 0; k < 3; k++)
          dfv[_e][k] = dfq.vect[k];

      } // Loop on cell edges

      /* Loop on border face edges */
      for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        cs_lnum_t  e_id = f2e->col_id[j];
        cs_quant_t  peq = quant->edge[e_id];
        cs_lnum_t  e_shft = e2v->idx[e_id];
        cs_lnum_t  v1_id = e2v->col_id[e_shft];
        cs_lnum_t  v2_id = e2v->col_id[e_shft+1];

        short int  _v1 = loc_v_ids[v1_id];
        short int  _v2 = loc_v_ids[v2_id];

        _e = loc_e_ids[e_id];

        /* Sanity checks */
        assert(_e != -1 && _v1 != -1 && _v2 != -1);

        const double  dp = _dp3(peq.unitv, dfv[_e]);
        const double  tmp_val = peq.meas*dp;
        const double  beta_pec_vol = 3.*beta/tmp_val;
        const double  ecoef =  3*beta/dp;
        const double  surf = cs_math_surftri(quant->vtx_coord + 3*v1_id,
                                             peq.center,
                                             pfq.center);

        /* Penalization term */
        v_coef[v1_id] += surf*f_coef;
        v_coef[v2_id] += surf*f_coef;

        /* Reset L_Ec(GRAD(p_j)) for each vertex of the cell */
        for (int i = 0; i < n_vbyc; i++)
          leg[i][0] = leg[i][1] = leg[i][2] = 0;

        /* Term related to the flux reconstruction:
           Compute L_Ec(GRAD(p_j)) on t_{_e,f} for each vertex j of the cell */
        for (je = start_e, _ek = 0; je < end_e; je++, _ek++) {

          cs_lnum_t  ek_id = c2e->ids[je];
          cs_lnum_t  ek_shft = e2v->idx[ek_id];
          cs_lnum_t  vj1_id = e2v->col_id[ek_shft];
          short int  sgn_j1 = e2v->sgn[ek_shft];
          cs_lnum_t  vj2_id = e2v->col_id[ek_shft+1];
          short int  sgn_j2 = e2v->sgn[ek_shft+1];

          short int  _vj1 = loc_v_ids[vj1_id];
          short int  _vj2 = loc_v_ids[vj2_id];

          /* Compute l_(ek,c)|p(_ef,c) */
          const double  eek_part = ecoef * _dp3(dfv[_ek], peq.unitv);

          for (k = 0; k < 3; k++)
            lek[k] = over_c_vol*(dfv[_ek][k] - eek_part*dfv[_e][k]);
          if (_ek == _e)
            for (k = 0; k < 3; k++)
              lek[k] += dfv[_ek][k]*beta_pec_vol;

          for (k = 0; k < 3; k++) {
            leg[_vj1][k] += sgn_j1*lek[k];
            leg[_vj2][k] += sgn_j2*lek[k];
          }

        } // Loop on cell edges

        for (int _vj = 0; _vj < n_vbyc; _vj++) {
          const double  contrib = _dp3(mn, leg[_vj])*surf;
          ntrgrd->mat[_v1*n_vbyc+_vj] += contrib;
          ntrgrd->mat[_v2*n_vbyc+_vj] += contrib;
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

#if CS_CDOVB_DIFFUSION_DBG > 1
  bft_printf(">> Local weak bc matrix (f_id: %d)", f_id);
  cs_locmat_dump(c_id, ntrgrd);
#endif

  return ntrgrd;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
