/*============================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local quantities related to a cell.
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

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_locsys_t structure
 *
 * \param[in]   n_max_ent    max number of entries
 *
 * \return a pointer to a new allocated cs_cdo_locsys_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_locsys_t *
cs_cdo_locsys_create(int    n_max_ent)
{
  cs_cdo_locsys_t  *ls = NULL;

  BFT_MALLOC(ls, 1, cs_cdo_locsys_t);

  ls->mat = cs_locmat_create(n_max_ent);

  if (n_max_ent > 0) {
    BFT_MALLOC(ls->rhs, n_max_ent, double);
    BFT_MALLOC(ls->dir_bc, n_max_ent, double);
  }
  else {
    ls->rhs = NULL;
    ls->dir_bc = NULL;
  }

  return ls;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_locsys_t structure
 *
 * \param[in, out]  p_ls   pointer of pointer to a cs_cdo_locsys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_locsys_free(cs_cdo_locsys_t     **p_ls)
{
  cs_cdo_locsys_t  *ls = *p_ls;

  if (ls == NULL)
    return;

  ls->mat = cs_locmat_free(ls->mat);
  BFT_FREE(ls->rhs);
  BFT_FREE(ls->dir_bc);

  BFT_FREE(ls);
  *p_ls= NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_locmesh_t structure
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cdo_locmesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_locmesh_t *
cs_cdo_locmesh_create(const cs_cdo_connect_t     *connect)
{
  cs_cdo_locmesh_t  *lm = NULL;

  if (connect == NULL)
    return lm;

  BFT_MALLOC(lm, 1, cs_cdo_locmesh_t);

  lm->n_max_vbyc = connect->n_max_vbyc;
  lm->n_max_ebyc = connect->n_max_ebyc;
  lm->n_max_fbyc = connect->n_max_fbyc;

  lm->flag = 0;
  lm->n_vc = 0;
  lm->n_ec = 0;
  lm->n_fc = 0;

  /* Vertex information */
  BFT_MALLOC(lm->v_ids, lm->n_max_vbyc, cs_lnum_t);
  BFT_MALLOC(lm->wvc, lm->n_max_vbyc, double);

  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  BFT_MALLOC(lm->vtag, n_vertices, short int);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    lm->vtag[i] = -1;

  /* Edge information */
  BFT_MALLOC(lm->e_ids, lm->n_max_ebyc, cs_lnum_t);
  BFT_MALLOC(lm->edge, lm->n_max_ebyc, cs_quant_t);
  BFT_MALLOC(lm->dface, lm->n_max_ebyc, cs_nvec3_t);

  const cs_lnum_t  n_edges = connect->e_info->n_elts;
  BFT_MALLOC(lm->etag, n_edges, short int);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_edges; i++)
    lm->etag[i] = -1;

  /* Face information */
  BFT_MALLOC(lm->f_ids, lm->n_max_fbyc, cs_lnum_t);
  BFT_MALLOC(lm->f_sgn, lm->n_max_fbyc, short int);
  BFT_MALLOC(lm->face, lm->n_max_fbyc, cs_quant_t);
  BFT_MALLOC(lm->dedge, lm->n_max_fbyc, cs_nvec3_t);

  /* edge --> vertices connectivity */
  BFT_MALLOC(lm->e2v_ids, 2*lm->n_max_ebyc, short int);
  BFT_MALLOC(lm->e2v_sgn, 2*lm->n_max_ebyc, short int);

  /* face --> edges connectivity */
  BFT_MALLOC(lm->f2e_idx, lm->n_max_fbyc + 1, short int);
  int  n_max_ebyf = connect->f2e->info.stencil_max;
  BFT_MALLOC(lm->f2e_ids, lm->n_max_fbyc*n_max_ebyf, short int);

  /* edge --> face connectivity */
  BFT_MALLOC(lm->e2f_ids, 2*lm->n_max_ebyc, short int);
  BFT_MALLOC(lm->e2f_sgn, 2*lm->n_max_ebyc, short int);

  return lm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_locmesh_t structure
 *
 * \param[in, out]  p_lm   pointer of pointer to a cs_cdo_locmesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_locmesh_free(cs_cdo_locmesh_t     **p_lm)
{
  cs_cdo_locmesh_t  *lm = *p_lm;

  if (lm == NULL)
    return;

  BFT_FREE(lm->vtag);
  BFT_FREE(lm->v_ids);
  BFT_FREE(lm->wvc);

  BFT_FREE(lm->e_ids);
  BFT_FREE(lm->etag);
  BFT_FREE(lm->edge);
  BFT_FREE(lm->dface);

  BFT_FREE(lm->f_ids);
  BFT_FREE(lm->f_sgn);
  BFT_FREE(lm->face);
  BFT_FREE(lm->dedge);

  BFT_FREE(lm->e2v_ids);
  BFT_FREE(lm->e2v_sgn);

  BFT_FREE(lm->f2e_idx);
  BFT_FREE(lm->f2e_ids);

  BFT_FREE(lm->e2f_ids);
  BFT_FREE(lm->e2f_sgn);

  BFT_FREE(lm);
  *p_lm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_cdo_locmesh_t structure for a given cell id. According
 *         to the requested level, some quantities may not be defined;
 *
 * \param[in]       c_id      cell id
 * \param[in]       flag      indicate which members are really defined
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  lm        pointer to a cs_cdo_locmesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_locmesh_build(cs_lnum_t                    c_id,
                     cs_flag_t                    flag,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     cs_cdo_locmesh_t            *lm)
{
  if (lm == NULL)
    return;

  lm->flag = flag;
  lm->c_id = c_id;
  lm->vol_c = quant->cell_vol[c_id];

  if (flag == 0)
    return;

  if (flag & CS_CDO_LOCAL_V) {

    const double  invvol = 1/lm->vol_c;
    const cs_connect_index_t  *c2v = connect->c2v;
    const cs_lnum_t  vstart = c2v->idx[c_id];
    const cs_lnum_t  vend = c2v->idx[c_id+1];

    lm->n_vc = vend - vstart;
    for (cs_lnum_t i = vstart, ii = 0; i < vend; i++, ii++) {

      cs_lnum_t  v_id = c2v->ids[i];
      lm->vtag[v_id] = ii;
      lm->v_ids[ii] = v_id;
      lm->wvc[ii] = invvol * quant->dcell_vol[i];

    }

  } // vertex information

  if (flag & CS_CDO_LOCAL_E) {

    const cs_connect_index_t  *c2e = connect->c2e;
    const cs_lnum_t  estart = c2e->idx[c_id];
    const cs_lnum_t  eend = c2e->idx[c_id+1];

    lm->n_ec = eend - estart;
    for (cs_lnum_t i = estart, ii = 0; i < eend; i++, ii++) {

      const cs_lnum_t  e_id = c2e->ids[i];
      const cs_quant_t  qe = quant->edge[e_id];
      const cs_dface_t  qdf = quant->dface[i];

      lm->e_ids[ii] = e_id;
      lm->etag[e_id] = ii;

      const double  meas = qdf.sface[0].meas + qdf.sface[1].meas;
      const double inv_meas = 1/meas;

      lm->dface[ii].meas = meas;
      lm->edge[ii].meas = qe.meas;
      for (int k = 0; k < 3; k++) {
        lm->edge[ii].center[k] = qe.center[k];
        lm->edge[ii].unitv[k] = qe.unitv[k];
        lm->dface[ii].unitv[k] = inv_meas*qdf.vect[k];
      }

    } // Loop on cell edges

    if (flag & CS_CDO_LOCAL_EV) {

      for (cs_lnum_t i = estart, ii = 0; i < eend; i++, ii++) {

        const cs_lnum_t  e_id = lm->e_ids[ii];
        const cs_lnum_t  eshft = 2*e_id, _eshft = 2*ii;
        const cs_lnum_t  *ids = connect->e2v->col_id + eshft;
        const short int  *sgn = connect->e2v->sgn + eshft;

        lm->e2v_ids[_eshft] = lm->vtag[ids[0]];
        lm->e2v_ids[_eshft+1] = lm->vtag[ids[1]];
        lm->e2v_sgn[_eshft] = sgn[0];
        lm->e2v_sgn[_eshft+1] = sgn[1];

      } // Loop on cell edges

    } // edge --> vertices connectivity

  } // edge information

  if (flag & CS_CDO_LOCAL_F) {

    const cs_sla_matrix_t  *c2f = connect->c2f;

    cs_lnum_t  fstart = c2f->idx[c_id];
    cs_lnum_t  fend = c2f->idx[c_id+1];

    lm->n_fc = fend - fstart;
    for (cs_lnum_t i = fstart, ii = 0; i < fend; i++, ii++) {

      const cs_lnum_t  f_id = c2f->col_id[i];

      lm->f_ids[ii] = f_id;
      lm->f_sgn[ii] = c2f->sgn[i];

      const cs_quant_t  pfq = quant->face[f_id];
      const cs_nvec3_t  deq = quant->dedge[i];

      /* Related quantities */
      lm->face[ii].meas = pfq.meas;
      lm->dedge[ii].meas = deq.meas;
      for (int k = 0; k < 3; k++) {
        lm->face[ii].center[k] = pfq.center[k];
        lm->face[ii].unitv[k] = pfq.unitv[k];
        lm->dedge[ii].unitv[k] = deq.unitv[k];
      }

    } // Loop on cell faces

    if (flag & CS_CDO_LOCAL_FE) {

      const cs_sla_matrix_t  *f2e = connect->f2e;

      assert(flag & CS_CDO_LOCAL_E);
      lm->f2e_idx[0] = 0;
      for (int ii = 0; ii < lm->n_fc; ii++) {

        cs_lnum_t  f_id = lm->f_ids[ii];
        cs_lnum_t  festart = f2e->idx[f_id];
        cs_lnum_t  feend = f2e->idx[f_id+1];
        short int  _shft = lm->f2e_idx[ii];

        lm->f2e_idx[ii+1] = _shft + feend - festart;
        for (int j = 0; j < feend - festart; j++) {
          cs_lnum_t  e_id = f2e->col_id[festart + j];
          lm->f2e_ids[_shft + j] = lm->etag[e_id];
        }

      } // Loop on cell faces

    } // face --> edges connectivity

  } // face information

  if (flag & CS_CDO_LOCAL_EF) {

    assert((flag & CS_CDO_LOCAL_E) && (flag & CS_CDO_LOCAL_F));

    short int  *ecount = NULL;
    BFT_MALLOC(ecount, lm->n_ec, short int);
    for (int ii = 0; ii < lm->n_ec; ii++)
      ecount[ii] = 0;

    const cs_sla_matrix_t  *f2e = connect->f2e;

    for (short int f = 0; f < lm->n_fc; f++) {

      cs_lnum_t  f_id = lm->f_ids[f];

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        cs_lnum_t  e_id = f2e->col_id[j];
        short int  e = lm->etag[e_id];
        short int  shft = 2*e + ecount[e];

        assert(ecount[e] < 2);
        lm->e2f_ids[shft] = f;
        lm->e2f_sgn[shft] = f2e->sgn[j];
        ecount[e] += 1;

      } /* Loop on face edges */

    } // Loop on cell faces

    BFT_FREE(ecount);

  } // edge --> faces connectivity

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
