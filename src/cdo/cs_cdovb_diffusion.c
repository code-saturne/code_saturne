/*============================================================================
 * Build discrete stiffness matrices and handled boundary conditions for the
 * diffusion term in CDO vertex-based schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hodge.h"

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
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDO_DIFFUSION_DBG 0

/* Stiffness matrix builder */
struct _cs_cdovb_diff_t {

  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs

  /* Data related to the discrete Hodge operator attached to the
     diffusion property */
  cs_param_hodge_t      h_info;
  cs_hodge_builder_t   *hb;

  /* Compact way to stored the edge --> vertices connect. related to a cell */
  int           n_bits;         // number of bits in a block
  int           n_blocks;       // number of blocks in a mask
  cs_flag_t    *emsk;           // list of masks to store the connectivity

  /* Members related to the weak enforcement of BCs */
  cs_real_3_t  *local_vect;  // set of local vectors
  cs_real_t    *local_vol;   // set of local volumes
  short int    *local_ids;   // set of local ids related to a cell

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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a builder structure used to build the stiffness matrix
 *
 * \param[in] connect     pointer to a cs_cdo_connect_t struct.
 * \param[in] time_step   pointer to a time step structure
 * \param[in] h_info      cs_param_hodge_t struct.
 * \param[in] bc_enforce  type of boundary enforcement for Dirichlet values
 *
 * \return a pointer to a new allocated cs_cdovb_diffusion_builder_t struc.
 */
/*----------------------------------------------------------------------------*/

cs_cdovb_diff_t *
cs_cdovb_diffusion_builder_init(const cs_cdo_connect_t       *connect,
                                const cs_time_step_t         *time_step,
                                const cs_param_hodge_t        h_info,
                                const cs_param_bc_enforce_t   bc_enforce)
{
  cs_lnum_t  i;

  cs_cdovb_diff_t  *diff = NULL;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);

  BFT_MALLOC(diff, 1, cs_cdovb_diff_t);

  /* Copy the data related to a discrete Hodge operator */
  diff->h_info.type = h_info.type;
  diff->h_info.pty_id = h_info.pty_id;
  diff->h_info.inv_pty = h_info.inv_pty;
  diff->h_info.algo = h_info.algo;
  diff->h_info.coef = h_info.coef;

  /* Initialize mask features */
  diff->n_bits = sizeof(cs_flag_t)*CHAR_BIT;
  diff->n_blocks = connect->n_max_ebyc/diff->n_bits;
  if (connect->n_max_ebyc % diff->n_bits != 0)
    diff->n_blocks += 1;

  BFT_MALLOC(diff->emsk, diff->n_blocks*connect->n_max_vbyc, cs_flag_t);
  for (i = 0; i < diff->n_blocks*connect->n_max_vbyc; i++)
    diff->emsk[i] = 0;

  diff->enforce = bc_enforce;
  if (bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      bc_enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    
    BFT_MALLOC(diff->local_vect, 2*connect->n_max_ebyc, cs_real_3_t);
    BFT_MALLOC(diff->local_vol, connect->n_max_ebyc, cs_real_t);

    cs_lnum_t  n_ent = CS_MAX(connect->e_info->n_ent, connect->f_info->n_ent);
    BFT_MALLOC(diff->local_ids, n_ent, short int);
    for (i = 0; i < n_ent; i++)
      diff->local_ids[i] = -1;

  } // Weakly enforcement of Dirichlet BCs

  /* Define a builder for the related discrete Hodge operator */
  diff->hb = cs_hodge_builder_init(connect, time_step, h_info);

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

  BFT_FREE(diff->emsk);

  if (diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      diff->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    BFT_FREE(diff->local_vect);
    BFT_FREE(diff->local_vol);
    BFT_FREE(diff->local_ids);

  }

  diff->hb = cs_hodge_builder_free(diff->hb);
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

      /* Build a local discrete Hodge operator and return a local dense matrix */
      const cs_locmat_t  *hloc = cs_hodge_build_local(c_id,
                                                      connect, quant,
                                                      diff->hb);

      /* Encode edge sets related to each vertex belonging to the current cell */
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
    // TODO
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
  cs_locmat_t  *ntrgrd = diff->loc; // Local matrix to build
  cs_real_3_t  *dfv = diff->local_vect;
  cs_real_3_t  *pev = diff->local_vect + ntrgrd->n_ent;
  cs_real_t  *over_pec_vol = diff->local_vol;
  short int  *loc_e_ids = diff->local_ids;

  const cs_param_hodge_t  h_info = diff->h_info;
  const cs_param_hodge_algo_t  h_algo = h_info.algo;

  /* Initialize the local matrix */
  for (int  i = 0; i < ntrgrd->n_ent*ntrgrd->n_ent; i++)
    ntrgrd->mat[i] = 0;

  /* Build the local "normal trace gradient" according to the choice of
     algorithm use to build the discrete Hodge operator */
  switch (h_algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    {
      cs_lnum_t   ie, j, k, _id, v_id;
      cs_real_3_t  mn, reco_val;

      const cs_sla_matrix_t  *e2v = connect->e2v;
      const cs_sla_matrix_t  *f2e = connect->f2e;
      const cs_connect_index_t  *c2e = connect->c2e;
      const cs_connect_index_t  *c2v = connect->c2v;
      const double  beta = h_info.coef;
      const cs_quant_t  qf = quant->face[f_id];
      const cs_real_t  over_c_vol = 1/quant->cell_vol[c_id];

      double  f_coef = pow(qf.meas, -0.5) * eig_ratio * eig_max;

      /* Compute the product: matpty*face unit normal */
      _mv3((const cs_real_t (*)[3])matpty, qf.unitv, mn);

      /*  Initialize buffers related to vertices.
          Define an id local to this cell for each vertex */
      for (j = c2v->idx[c_id], _id = 0; j < c2v->idx[c_id+1]; j++, _id++) {
        v_id = c2v->ids[j];
        loc_v_ids[v_id] = _id;
        ntrgrd->ids[_id] = v_id;
      }
      ntrgrd->n_ent = _id;

      /* Store the local id and useful quantities for each edge */
      for (j = c2e->idx[c_id], _id = 0; j < c2e->idx[c_id+1]; j++, _id++) {

        cs_lnum_t  e_id = c2e->ids[j];
        cs_quant_t qpe = quant->edge[e_id];
        cs_dface_t qdf = quant->dface[j];

        loc_e_ids[e_id] = _id;
        for (k = 0; k < 3; k++) {
          pev[_id][k] = qpe.unitv[k]*qpe.meas;
          dfv[_id][k] = qdf.vect[k];
        }
        over_pec_vol[_id] = 3./_dp3(pev[_id], dfv[_id]);

      } // Loop on cell edges

      /* Loop on border face edges */
      for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        cs_lnum_t  e_id = f2e->col_id[j];
        cs_quant_t  qe = quant->edge[e_id];
        cs_lnum_t  e_shft = e2v->idx[e_id];
        cs_lnum_t  v1_id = e2v->col_id[e_shft];
        cs_lnum_t  v2_id = e2v->col_id[e_shft+1];

        short int  _e = loc_e_ids[e_id];
        short int  _v1 = loc_v_ids[v1_id];
        short int  _v2 = loc_v_ids[v2_id];

        /* Sanity checks */
        assert(_e != -1 && _v1 != -1 && _v2 != -1);

        double surf = cs_surftri(&(quant->vtx_coord[3*v1_id]),
                                 qe.center,
                                 qf.center);

        v_coef[v1_id] += surf*f_coef;
        v_coef[v2_id] += surf*f_coef;

        for (ie = c2e->idx[c_id]; ie < c2e->idx[c_id+1]; ie++) {

          cs_lnum_t  ek_id = c2e->ids[ie];
          cs_lnum_t  ek_shft = e2v->idx[ek_id];
          cs_lnum_t  vj1_id = e2v->col_id[ek_shft];
          short int  sgn_j1 = e2v->sgn[ek_shft];
          cs_lnum_t  vj2_id = e2v->col_id[ek_shft+1];
          short int  sgn_j2 = e2v->sgn[ek_shft+1];

          short int  _ek = loc_e_ids[ek_id];
          short int  _vj1 = loc_v_ids[vj1_id];
          short int  _vj2 = loc_v_ids[vj2_id];

          /* Compute L_Ec(GRAD(p_j)) on t_{ek,f} */
          if (_ek == _e)
            for (k = 0; k < 3; k++)
              reco_val[k] = dfv[_ek][k]*beta*over_pec_vol[_ek];
          else
            for (k = 0; k < 3; k++)
              reco_val[k] = 0.0;

          cs_real_t  dp = beta*over_pec_vol[_e]*_dp3(dfv[_ek], pev[_e]);

          for (k = 0; k < 3; k++)
            reco_val[k] += over_c_vol*(dfv[_ek][k] - dp*dfv[_e][k]);

          cs_real_t  contrib = _dp3(mn, reco_val)*surf;

          /* Update the coefficients of the local normal trace operator */
          ntrgrd->mat[_v1*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
          ntrgrd->mat[_v1*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;
          ntrgrd->mat[_v2*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
          ntrgrd->mat[_v2*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;

        } // Loop on cell edges

      } // border face edges


    }
    break;

  /* case CS_PARAM_HODGE_ALGO_WBS: */
  /*   // TODO */
  /*   break; */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "Invalid type of algorithm to weakly enforce Dirichlet BCs.");

  } // End of switch

  return ntrgrd;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
