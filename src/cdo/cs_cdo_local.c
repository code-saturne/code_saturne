/* ===========================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local mesh structure related to a cell or to a couple cell/face
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

#include "cs_math.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

cs_cell_mesh_t   **cs_cdo_local_cell_meshes = NULL;
cs_face_mesh_t   **cs_cdo_local_face_meshes = NULL;

static int  cs_cdo_local_n_structures = 0;

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
 * \brief  Allocate global structures related to cs_cell_mesh_t and
 *         cs_face_mesh_t structures
 *
 * \param[in]   connect   pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_initialize(const cs_cdo_connect_t     *connect)
{
  /* Sanity check */
  assert(cs_glob_n_threads > 0);

  int  size = cs_glob_n_threads;
  cs_cdo_local_n_structures = size;
  BFT_MALLOC(cs_cdo_local_cell_meshes, size, cs_cell_mesh_t *);
  BFT_MALLOC(cs_cdo_local_face_meshes, size, cs_face_mesh_t *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdo_local_cell_meshes[t_id] = cs_cell_mesh_create(connect);
    cs_cdo_local_face_meshes[t_id] = cs_face_mesh_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdo_local_cell_meshes[0] = cs_cell_mesh_create(connect);
  cs_cdo_local_face_meshes[0] = cs_face_mesh_create(connect);
#endif /* openMP */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to cs_cell_mesh_t and cs_face_mesh_t
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_finalize(void)
{
  if (cs_cdo_local_n_structures < 1)
    return;

  assert(cs_cdo_local_n_structures == cs_glob_n_threads);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_mesh_free(&(cs_cdo_local_cell_meshes[t_id]));
    cs_face_mesh_free(&(cs_cdo_local_face_meshes[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_mesh_free(&(cs_cdo_local_cell_meshes[0]));
  cs_face_mesh_free(&(cs_cdo_local_face_meshes[0]));
#endif /* openMP */

  BFT_FREE(cs_cdo_local_cell_meshes);
  BFT_FREE(cs_cdo_local_face_meshes);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_cell_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_cell_mesh_t struct.
 *
 * \return a pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cdo_local_get_cell_mesh(int    mesh_id)
{
  if (mesh_id < 0 || mesh_id >= cs_glob_n_threads)
    return NULL;

  return cs_cdo_local_cell_meshes[mesh_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_face_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_face_mesh_t struct.
 *
 * \return a pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_cdo_local_get_face_mesh(int    mesh_id)
{
  if (mesh_id < 0 || mesh_id >= cs_glob_n_threads)
    return NULL;

  return cs_cdo_local_face_meshes[mesh_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cell_mesh_t structure
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cell_mesh_create(const cs_cdo_connect_t     *connect)
{
  cs_cell_mesh_t  *cm = NULL;

  if (connect == NULL)
    return cm;

  BFT_MALLOC(cm, 1, cs_cell_mesh_t);

  cm->n_max_vbyc = connect->n_max_vbyc;
  cm->n_max_ebyc = connect->n_max_ebyc;
  cm->n_max_fbyc = connect->n_max_fbyc;

  cm->flag = 0;
  cm->n_vc = 0;
  cm->n_ec = 0;
  cm->n_fc = 0;

  cm->c_id = -1;
  cm->xc = NULL;
  cm->vol_c = 0.;

  /* Vertex information */
  BFT_MALLOC(cm->v_ids, cm->n_max_vbyc, cs_lnum_t);
  BFT_MALLOC(cm->wvc, cm->n_max_vbyc, double);
  BFT_MALLOC(cm->xv, 3*cm->n_max_vbyc, double);

  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  BFT_MALLOC(cm->vtag, n_vertices, short int);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    cm->vtag[i] = -1;

  /* Edge information */
  BFT_MALLOC(cm->e_ids, cm->n_max_ebyc, cs_lnum_t);
  BFT_MALLOC(cm->edge, cm->n_max_ebyc, cs_quant_t);
  BFT_MALLOC(cm->dface, cm->n_max_ebyc, cs_nvec3_t);

  const cs_lnum_t  n_edges = connect->e_info->n_elts;
  BFT_MALLOC(cm->etag, n_edges, short int);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_edges; i++)
    cm->etag[i] = -1;

  /* Face information */
  BFT_MALLOC(cm->f_ids, cm->n_max_fbyc, cs_lnum_t);
  BFT_MALLOC(cm->f_sgn, cm->n_max_fbyc, short int);
  BFT_MALLOC(cm->face, cm->n_max_fbyc, cs_quant_t);
  BFT_MALLOC(cm->dedge, cm->n_max_fbyc, cs_nvec3_t);

  /* edge --> vertices connectivity */
  BFT_MALLOC(cm->e2v_ids, 2*cm->n_max_ebyc, short int);
  BFT_MALLOC(cm->e2v_sgn, 2*cm->n_max_ebyc, short int);

  /* face --> edges connectivity */
  BFT_MALLOC(cm->f2e_idx, cm->n_max_fbyc + 1, short int);
  BFT_MALLOC(cm->f2e_ids, 2*cm->n_max_ebyc, short int);

  /* edge --> face connectivity */
  BFT_MALLOC(cm->e2f_ids, 2*cm->n_max_ebyc, short int);
  BFT_MALLOC(cm->e2f_sgn, 2*cm->n_max_ebyc, short int);

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_mesh_t structure
 *
 * \param[in, out]  p_cm   pointer of pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_free(cs_cell_mesh_t     **p_cm)
{
  cs_cell_mesh_t  *cm = *p_cm;

  if (cm == NULL)
    return;

  BFT_FREE(cm->vtag);
  BFT_FREE(cm->v_ids);
  BFT_FREE(cm->wvc);
  BFT_FREE(cm->xv);

  BFT_FREE(cm->e_ids);
  BFT_FREE(cm->etag);
  BFT_FREE(cm->edge);
  BFT_FREE(cm->dface);

  BFT_FREE(cm->f_ids);
  BFT_FREE(cm->f_sgn);
  BFT_FREE(cm->face);
  BFT_FREE(cm->dedge);

  BFT_FREE(cm->e2v_ids);
  BFT_FREE(cm->e2v_sgn);

  BFT_FREE(cm->f2e_idx);
  BFT_FREE(cm->f2e_ids);

  BFT_FREE(cm->e2f_ids);
  BFT_FREE(cm->e2f_sgn);

  BFT_FREE(cm);
  *p_cm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_cell_mesh_t structure for a given cell id. According
 *         to the requested level, some quantities may not be defined;
 *
 * \param[in]       c_id      cell id
 * \param[in]       flag      indicate which members are really defined
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  cm        pointer to a cs_cell_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_build(cs_lnum_t                    c_id,
                   cs_flag_t                    flag,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_cell_mesh_t              *cm)
{
  if (cm == NULL)
    return;

  cm->flag = flag;
  cm->c_id = c_id;
  cm->xc = quant->cell_centers + 3*c_id;
  cm->vol_c = quant->cell_vol[c_id];

  if (flag == 0)
    return;

  if (flag & CS_CDO_LOCAL_V) {

    const double  invvol = 1/cm->vol_c;
    const cs_connect_index_t  *c2v = connect->c2v;
    const cs_lnum_t  vstart = c2v->idx[c_id];
    const cs_lnum_t  vend = c2v->idx[c_id+1];

    cm->n_vc = vend - vstart;
    for (cs_lnum_t i = vstart, ii = 0; i < vend; i++, ii++) {

      const cs_lnum_t  v_id = c2v->ids[i];
      const cs_real_t  *xv = quant->vtx_coord + 3*v_id;

      for (int k = 0; k < 3; k++)
        cm->xv[3*ii+k] = xv[k];

      cm->vtag[v_id] = ii;
      cm->v_ids[ii] = v_id;
      cm->wvc[ii] = invvol * quant->dcell_vol[i];

    }

  } // vertex information

  if (flag & CS_CDO_LOCAL_E) {

    const cs_connect_index_t  *c2e = connect->c2e;
    const cs_lnum_t  estart = c2e->idx[c_id];
    const cs_lnum_t  eend = c2e->idx[c_id+1];

    cm->n_ec = eend - estart;
    for (cs_lnum_t i = estart, ii = 0; i < eend; i++, ii++) {

      const cs_lnum_t  e_id = c2e->ids[i];
      const cs_quant_t  qe = quant->edge[e_id];
      const cs_dface_t  qdf = quant->dface[i];

      cm->e_ids[ii] = e_id;
      cm->etag[e_id] = ii;

      const double  meas = qdf.sface[0].meas + qdf.sface[1].meas;
      const double inv_meas = 1/meas;

      cm->dface[ii].meas = meas;
      cm->edge[ii].meas = qe.meas;
      for (int k = 0; k < 3; k++) {
        cm->edge[ii].center[k] = qe.center[k];
        cm->edge[ii].unitv[k] = qe.unitv[k];
        cm->dface[ii].unitv[k] = inv_meas*qdf.vect[k];
      }

    } // Loop on cell edges

    if (flag & CS_CDO_LOCAL_EV) {

      for (cs_lnum_t i = estart, ii = 0; i < eend; i++, ii++) {

        const cs_lnum_t  e_id = cm->e_ids[ii];
        const cs_lnum_t  eshft = 2*e_id, _eshft = 2*ii;
        const cs_lnum_t  *ids = connect->e2v->col_id + eshft;
        const short int  *sgn = connect->e2v->sgn + eshft;

        cm->e2v_ids[_eshft] = cm->vtag[ids[0]];
        cm->e2v_ids[_eshft+1] = cm->vtag[ids[1]];
        cm->e2v_sgn[_eshft] = sgn[0];
        cm->e2v_sgn[_eshft+1] = sgn[1];

      } // Loop on cell edges

    } // edge --> vertices connectivity

  } // edge information

  if (flag & CS_CDO_LOCAL_F) {

    const cs_sla_matrix_t  *c2f = connect->c2f;

    cs_lnum_t  fstart = c2f->idx[c_id];
    cs_lnum_t  fend = c2f->idx[c_id+1];

    cm->n_fc = fend - fstart;
    for (cs_lnum_t i = fstart, ii = 0; i < fend; i++, ii++) {

      const cs_lnum_t  f_id = c2f->col_id[i];

      cm->f_ids[ii] = f_id;
      cm->f_sgn[ii] = c2f->sgn[i];

      const cs_quant_t  pfq = quant->face[f_id];
      const cs_nvec3_t  deq = quant->dedge[i];

      /* Related quantities */
      cm->face[ii].meas = pfq.meas;
      cm->dedge[ii].meas = deq.meas;
      for (int k = 0; k < 3; k++) {
        cm->face[ii].center[k] = pfq.center[k];
        cm->face[ii].unitv[k] = pfq.unitv[k];
        cm->dedge[ii].unitv[k] = deq.unitv[k];
      }

    } // Loop on cell faces

    if (flag & CS_CDO_LOCAL_FE) {

      const cs_sla_matrix_t  *f2e = connect->f2e;

      assert(flag & CS_CDO_LOCAL_E);
      cm->f2e_idx[0] = 0;
      for (int ii = 0; ii < cm->n_fc; ii++) {

        cs_lnum_t  f_id = cm->f_ids[ii];
        cs_lnum_t  festart = f2e->idx[f_id];
        cs_lnum_t  feend = f2e->idx[f_id+1];
        short int  _shft = cm->f2e_idx[ii];

        cm->f2e_idx[ii+1] = _shft + feend - festart;
        for (int j = 0; j < feend - festart; j++) {
          cs_lnum_t  e_id = f2e->col_id[festart + j];
          cm->f2e_ids[_shft + j] = cm->etag[e_id];
        }

      } // Loop on cell faces

    } // face --> edges connectivity

  } // face information

  if (flag & CS_CDO_LOCAL_EF) {

    assert((flag & CS_CDO_LOCAL_E) && (flag & CS_CDO_LOCAL_F));

    short int  *ecount = NULL;
    BFT_MALLOC(ecount, cm->n_ec, short int);
    for (int ii = 0; ii < cm->n_ec; ii++)
      ecount[ii] = 0;

    const cs_sla_matrix_t  *f2e = connect->f2e;

    for (short int f = 0; f < cm->n_fc; f++) {

      cs_lnum_t  f_id = cm->f_ids[f];

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        cs_lnum_t  e_id = f2e->col_id[j];
        short int  e = cm->etag[e_id];
        short int  shft = 2*e + ecount[e];

        assert(ecount[e] < 2);
        cm->e2f_ids[shft] = f;
        cm->e2f_sgn[shft] = f2e->sgn[j];
        ecount[e] += 1;

      } /* Loop on face edges */

    } // Loop on cell faces

    BFT_FREE(ecount);

  } // edge --> faces connectivity

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_t structure
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_face_mesh_create(const cs_cdo_connect_t     *connect)
{
  cs_face_mesh_t  *fm = NULL;

  if (connect == NULL)
    return fm;

  BFT_MALLOC(fm, 1, cs_face_mesh_t);

  fm->n_max_vbyf = connect->n_max_vbyf;

  fm->c_id = -1;
  fm->xc = NULL;

  /* Face-related quantities */
  fm->f_id = -1;
  fm->f_sgn = 0;

  /* Vertex-related quantities */
  fm->n_vf = 0;
  BFT_MALLOC(fm->v_ids, fm->n_max_vbyf, cs_lnum_t);
  BFT_MALLOC(fm->xv, 3*fm->n_max_vbyf, double);
  BFT_MALLOC(fm->wvf, fm->n_max_vbyf, double);

  /* Edge-related quantities */
  fm->n_ef = 0;
  BFT_MALLOC(fm->e_ids, fm->n_max_vbyf, cs_lnum_t);
  BFT_MALLOC(fm->edge,  fm->n_max_vbyf, cs_quant_t);
  BFT_MALLOC(fm->e2v_ids, 2*fm->n_max_vbyf, short int);
  BFT_MALLOC(fm->tef, fm->n_max_vbyf, double);

  return fm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_face_mesh_t structure
 *
 * \param[in, out]  p_fm   pointer of pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_free(cs_face_mesh_t     **p_fm)
{
  cs_face_mesh_t  *fm = *p_fm;

  if (fm == NULL)
    return;

  BFT_FREE(fm->v_ids);
  BFT_FREE(fm->xv);
  BFT_FREE(fm->wvf);

  BFT_FREE(fm->e_ids);
  BFT_FREE(fm->edge);
  BFT_FREE(fm->e2v_ids);
  BFT_FREE(fm->tef);

  BFT_FREE(fm);
  *p_fm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given face/cell id.
 *
 * \param[in]       c_id      cell id
 * \param[in]       f_id      face id in the mesh structure
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build(cs_lnum_t                    c_id,
                   cs_lnum_t                    f_id,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_face_mesh_t              *fm)
{
  if (fm == NULL)
    return;

  /* Sanity checks */
  assert(c_id > -1);
  assert(f_id > -1);

  fm->c_id = c_id;
  fm->xc = quant->cell_centers + 3*c_id;

  /* Face-related quantities */
  const cs_quant_t  pfq = quant->face[f_id];

  fm->f_id = f_id;
  fm->face.meas = pfq.meas;
  for (int k = 0; k < 3; k++) {
    fm->face.center[k] = pfq.center[k];
    fm->face.unitv[k] = pfq.unitv[k];
  }

  const cs_sla_matrix_t  *c2f = connect->c2f;
  bool  found = false;
  for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {
    if (c2f->col_id[i] == f_id) {

      const cs_nvec3_t  deq = quant->dedge[i];

      fm->dedge.meas = deq.meas;
      for (int k = 0; k < 3; k++)
        fm->dedge.unitv[k] = deq.unitv[k];
      fm->f_sgn = c2f->sgn[i];
      found = true;
      break;
    }
  }

  if (!found) // Sanity check
    bft_error(__FILE__, __LINE__, 0,
              _(" Face %d not found.\n Stop build a face mesh."), f_id);

  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_lnum_t  festart = f2e->idx[f_id];
  const cs_lnum_t  feend = f2e->idx[f_id+1];

  fm->n_vf = fm->n_ef = feend - festart;
  short int nv = 0;
  for (int i = 0; i < fm->n_vf; i++)
    fm->v_ids[i] = -1;

  for (cs_lnum_t j = festart, je = 0; j < feend; j++, je++) {

    const cs_lnum_t  e_id = f2e->col_id[j];
    const cs_quant_t  qe = quant->edge[e_id];

    fm->e_ids[je] = e_id;
    fm->edge[je].meas = qe.meas;
    for (int k = 0; k < 3; k++) {
      fm->edge[je].center[k] = qe.center[k];
      fm->edge[je].unitv[k] = qe.unitv[k];
    }

    const cs_lnum_t  eshft = 2*e_id;
    cs_lnum_t  v1_id = e2v->col_id[eshft];
    cs_lnum_t  v2_id = e2v->col_id[eshft+1];

    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf; v++) {
      if (fm->v_ids[v] == -1) break; // Reached the end of the list of vertices
      else {
        if (fm->v_ids[v] == v1_id) v1 = v;
        else if (fm->v_ids[v] == v2_id) v2 = v;
      }
    }

    /* Add vertices if not already identified */
    if (v1 == -1) // Not found -> Add v1
      fm->v_ids[nv] = v1_id, v1 = nv++;
    if (v2 == -1) // Not found -> Add v2
      fm->v_ids[nv] = v2_id, v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*je;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } // Loop on face edges

  assert(nv == fm->n_vf); // Sanity check

  /* Update vertex coordinates */
  int  shift = 0;
  for (short int v = 0; v < fm->n_vf; v++) {
    const cs_real_t *xv = quant->vtx_coord + 3*fm->v_ids[v];
    for (int k = 0; k < 3; k++)
      fm->xv[shift++] = xv[k];
  }

  /* Define wvf and tef */
  for (int i = 0; i < fm->n_vf; i++)
    fm->wvf[i] = 0;

  for (short int e = 0; e < fm->n_ef; e++) {

    double  xef_len;
    cs_real_3_t  xef_un, un;

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const cs_quant_t  peq = fm->edge[e];

    cs_math_3_length_unitv(peq.center, pfq.center, &xef_len, xef_un);
    cs_math_3_cross_product(xef_un, peq.unitv, un);

    /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
    const double  tef = 0.5 * xef_len * peq.meas * cs_math_3_norm(un);

    fm->wvf[v1] += tef;
    fm->wvf[v2] += tef;
    fm->tef[e] = tef;

  }

  const double  invf = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++)
    fm->wvf[v] *= invf;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given cell from a
 *         cs_cell_mesh_t structure
 *
 * \param[in]       cm        pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f_id      face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build_from_cell_mesh(const cs_cell_mesh_t    *cm,
                                  short int                f,
                                  cs_face_mesh_t          *fm)
{
  if (fm == NULL || cm == NULL)
    return;

  /* Sanity checks */
  assert(f > -1 && f < cm->n_fc);

  fm->c_id = cm->c_id;
  fm->xc = cm->xc;

  /* Face-related quantities */
  const cs_quant_t  pfq = cm->face[f];

  fm->f_id = f;
  fm->f_sgn = cm->f_sgn[f];
  fm->face.meas = pfq.meas;
  for (int k = 0; k < 3; k++) {
    fm->face.center[k] = pfq.center[k];
    fm->face.unitv[k] = pfq.unitv[k];
  }

  const cs_nvec3_t  deq = cm->dedge[f];

  fm->dedge.meas = deq.meas;
  for (int k = 0; k < 3; k++)
    fm->dedge.unitv[k] = deq.unitv[k];

  const int  festart = cm->f2e_idx[f];
  const int  feend = cm->f2e_idx[f+1];

  fm->n_vf = fm->n_ef = feend - festart;
  short int nv = 0;
  for (int i = 0; i < fm->n_vf; i++)
    fm->v_ids[i] = -1;

  for (int j = festart, je = 0; j < feend; j++, je++) {

    const short int  e = cm->f2e_ids[j];
    const cs_quant_t  peq = cm->edge[e];

    fm->e_ids[je] = e;
    fm->edge[je].meas = peq.meas;
    for (int k = 0; k < 3; k++) {
      fm->edge[je].center[k] = peq.center[k];
      fm->edge[je].unitv[k] = peq.unitv[k];
    }

    const cs_lnum_t  eshft = 2*e;
    cs_lnum_t  v1c_id = cm->e2v_ids[eshft];
    cs_lnum_t  v2c_id = cm->e2v_ids[eshft+1];

    /* Compact vertex numbering to this face */
    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf; v++) {
      if (fm->v_ids[v] == -1) break; // Reached the end of the list of vertices
      else {
        if (fm->v_ids[v] == v1c_id) v1 = v;
        else if (fm->v_ids[v] == v2c_id) v2 = v;
      }
    }

    /* Add vertices if not already identified */
    if (v1 == -1) // Not found -> Add v1
      fm->v_ids[nv] = v1c_id, v1 = nv++;
    if (v2 == -1) // Not found -> Add v2
      fm->v_ids[nv] = v2c_id, v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*je;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } // Loop on face edges

  assert(nv == fm->n_vf); // Sanity check

  /* Update vertex coordinates */
  int  shift = 0;
  for (short int v = 0; v < fm->n_vf; v++) {
    const cs_real_t *xv = cm->xv + 3*fm->v_ids[v];
    for (int k = 0; k < 3; k++)
      fm->xv[shift++] = xv[k];
  }

  /* Define wvf and tef */
  for (int i = 0; i < fm->n_vf; i++)
    fm->wvf[i] = 0;

  for (short int e = 0; e < fm->n_ef; e++) {

    double  xef_len;
    cs_real_3_t  xef_un, un;

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const cs_quant_t  peq = fm->edge[e];

    cs_math_3_length_unitv(peq.center, pfq.center, &xef_len, xef_un);
    cs_math_3_cross_product(xef_un, peq.unitv, un);

    /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
    const double  tef = 0.5 * xef_len * peq.meas * cs_math_3_norm(un);

    fm->wvf[v1] += tef;
    fm->wvf[v2] += tef;
    fm->tef[e] = tef;

  }

  const double  invf = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++)
    fm->wvf[v] *= invf;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
