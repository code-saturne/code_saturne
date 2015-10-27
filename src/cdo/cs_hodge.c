/*============================================================================
 * Build discrete Hodge operators
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
#include <math.h>
#include <float.h>
#include <limits.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_sort.h"
#include "cs_evaluate.h"
#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hodge.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_hodge.c

  \brief Build discrete Hodge operators

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_HODGE_DBG 0

/* Main structure used to define a discrete Hodge operator */
struct _hodge_builder_t {

  cs_lnum_t   n_ent;             /* Number of entities */
  int         n_maxloc_ent;      /* Max local number of entities by primal
                                    cells (use for allocation) */

  cs_param_hodge_t  h_info;   /* Pointer the set of parameters related
                                 to the discrete Hodge op. to build.
                                 Not owned by this structure. */

  cs_real_t         t_cur;   /* Current physical time */
  bool              uniform; /* True if material property is uniform */
  cs_real_33_t      matval;  /* Tensor related to the material property */


  cs_locmat_t      *hloc;    /* Local dense matrix related to a local
                                discrete Hodge op. */

  void             *algoq;   /* Quantities used during the definition of
                                the local discrete Hodge op.
                                This structure is attached to each type
                                of algorithm */

};

/* Geometric quantities related to the construction of local discrete
   Hodge op. when the algo. is either COST (included DGA, SUSHI, Generalized
   Crouzeix-Raviart) and the type is  between edges and faces (primal or dual).
*/

struct _cost_quant_t {

  double     *invsvol;    /* size = n_ent*/
  double     *qmq;        /* symmetric dense matrix of size n_ent */
  double     *T;          /* dense matrix of size n_ent (not symmetric) */

  cs_qvect_t   *pq;       /* primal geometric quantity (size: n_ent) */
  cs_qvect_t   *dq;       /* dual geometric quantity (size: n_ent) */

};

/* Quantities related to the construction of a local discrete  Hodge op.
   when the Whitney Barycentric Subdivision algo. is employed.
   Used only for Vp --> Cd hodge (up to now) */

struct _wbs_quant_t {

  cs_connect_index_t  *f2v; /* face (border or interior) to vertices connect. */

  short int           *tag; /* size = n_vertices; default value = -1 (not used)
                               otherwise local id in [0, n_ent] */

  double   *wf;    /* weights related to each vertex of a given face */
  double   *wc;    /* weights related to each vertex of a given cell */
  double   *cumul; /* sum of temporary contribution (same size of wf and wc) */

};

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

static const double  invdim = 1./3.;

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize by default the matrix related to a discrete
 *          Hodge op. based on vertices
 *          Note: values are filled in a second step
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_hodge_vertex(const cs_cdo_connect_t     *connect,
                   const cs_cdo_quantities_t  *quant)
{
  int  i, j, shift;

  cs_connect_index_t  *v2v = NULL, *v2c = NULL;

  const cs_connect_index_t  *c2v = connect->c2v;
  const int  n_vertices = quant->n_vertices;

  /* Allocate and initialize the matrix */
  cs_sla_matrix_t  *h_mat = cs_sla_matrix_create(n_vertices,
                                                 n_vertices,
                                                 1,
                                                 CS_SLA_MAT_MSR,
                                                 false);

  /* Initialize index (v2v connectivity) */
  v2c = cs_index_transpose(n_vertices, c2v);
  v2v = cs_index_compose(n_vertices, v2c, c2v);
  cs_index_free(&v2c);

  cs_index_sort(v2v);
  h_mat->flag |= CS_SLA_MATRIX_SORTED;

  /* Update index */
  h_mat->idx[0] = 0;
  for (i = 0; i < n_vertices; i++)
    h_mat->idx[i+1] = h_mat->idx[i] + v2v->idx[i+1]-v2v->idx[i]-1;

  /* Fill column num */
  BFT_MALLOC(h_mat->col_id, h_mat->idx[n_vertices], cs_lnum_t);
  shift = 0;
  for (i = 0; i < n_vertices; i++)
    for (j = v2v->idx[i]; j < v2v->idx[i+1]; j++)
      if (v2v->ids[j] != i)
        h_mat->col_id[shift++] = v2v->ids[j];

  /* Sanity check */
  assert(shift == h_mat->idx[n_vertices]);

  /* Free temporary memory */
  cs_index_free(&v2v);

  /* Allocate and initialize value array */
  BFT_MALLOC(h_mat->val, h_mat->idx[n_vertices], double);
  for (i = 0; i < h_mat->idx[n_vertices]; i++)
    h_mat->val[i] = 0.0;

  for (i = 0; i < n_vertices; i++)
    h_mat->diag[i] = 0.0;

  return h_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize by default the matrix related to a discrete
 *          Hodge op. based on edges
 *          Note: values are filled in a second step
 *
 * \param[in]    connect  pointer to a cs_cdo_connect_t structure
 * \param[in]    quant    pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_hodge_edge(const cs_cdo_connect_t     *connect,
                 const cs_cdo_quantities_t  *quant)
{
  int  c_id, eid, eid2, i, j, count, shift;

  int  *etag = NULL;
  cs_connect_index_t  *e2f =  NULL, *f2c = NULL, *e2c = NULL, *e2e = NULL;

  const int  n_edges = quant->n_edges;
  const int  n_cells = quant->n_cells;
  const cs_connect_index_t  *c2e = connect->c2e;

  /* Allocate and initialize the matrix */
  cs_sla_matrix_t  *h_mat = cs_sla_matrix_create(n_edges,
                                                 n_edges,
                                                 1,
                                                 CS_SLA_MAT_MSR,
                                                 false);

  /* Build a edge -> cell connectivity */
  e2f = cs_index_map(connect->e2f->n_rows,
                     connect->e2f->idx,
                     connect->e2f->col_id);
  f2c = cs_index_map(connect->f2c->n_rows,
                     connect->f2c->idx,
                     connect->f2c->col_id);

  e2c = cs_index_compose(n_cells, e2f, f2c);

  /* Count nnz in H */
  BFT_MALLOC(etag, n_edges, int);
  for (eid = 0; eid < n_edges; eid++)
    etag[eid] = -1;

  for (eid = 0; eid < n_edges; eid++) {

    count = 0;
    for (i = e2c->idx[eid]; i < e2c->idx[eid+1]; i++) {
      c_id = e2c->ids[i];
      for (j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {
        eid2 = c2e->ids[j];
        if (eid != eid2 && etag[eid2] != eid) {
          etag[eid2] = eid;
          count++;
        }

      } /* Loop on edges sharing this cell */

    } /* Loop on cells sharing this edge (eid) */

    h_mat->idx[eid+1] = count;

  } /* End of loop on edges */

  /* Update index */
  for (i = 0; i < n_edges; i++)
    h_mat->idx[i+1] = h_mat->idx[i+1] + h_mat->idx[i];

  /* Fill column num */
  BFT_MALLOC(h_mat->col_id, h_mat->idx[n_edges], cs_lnum_t);
  for (i = 0; i < h_mat->idx[n_edges]; i++)
    h_mat->col_id[i] = -1;

  for (eid = 0; eid < n_edges; eid++)
    etag[eid] = -1;

  for (eid = 0; eid < n_edges; eid++) {

    shift = h_mat->idx[eid];
    for (i = e2c->idx[eid]; i < e2c->idx[eid+1]; i++) {

      c_id = e2c->ids[i];
      for (j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {

        eid2 = c2e->ids[j];
        if (eid != eid2 && etag[eid2] != eid) {
          etag[eid2] = eid;
          h_mat->col_id[shift++] = eid2;
        }

      } /* Loop on edges sharing this cell */

    } /* Loop on cells sharing this edge (eid) */

  } /* End of loop on edges */

  /* Order column entries in increasing order */
  e2e = cs_index_map(h_mat->n_rows, h_mat->idx, h_mat->col_id);
  cs_index_sort(e2e);
  h_mat->flag |= CS_SLA_MATRIX_SORTED;

  /* Partial buffer free */
  BFT_FREE(etag);
  cs_index_free(&e2e); /* Not owner. Only delete the link with Hodge index */
  cs_index_free(&f2c);
  cs_index_free(&e2f);
  cs_index_free(&e2c);

  /* Allocate and initialize value array */
  for (i = 0; i < n_edges; i++)
    h_mat->diag[i] = 0.0;

  BFT_MALLOC(h_mat->val, h_mat->idx[h_mat->n_rows], double);
  for (i = 0; i < h_mat->idx[n_edges]; i++)
    h_mat->val[i] = 0.0;

  return h_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize by default the matrix related to a discrete
 *          Hodge op. based on faces
 *          Note: values are filled in a second step
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_hodge_face(const cs_cdo_connect_t     *connect,
                 const cs_cdo_quantities_t  *quant)
{
  int  f_id, j, shift;

  cs_connect_index_t  *f2f = NULL, *c2f = NULL, *f2c = NULL;

  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_sla_matrix_t *mc2f = connect->c2f;
  const cs_sla_matrix_t *mf2c = connect->f2c;

  /* Allocate and initialize the matrix */
  cs_sla_matrix_t  *h_mat = cs_sla_matrix_create(n_faces,
                                                 n_faces,
                                                 1,
                                                 CS_SLA_MAT_MSR,
                                                 false);

  /* Build a face -> face connectivity */
  f2c = cs_index_map(mf2c->n_rows, mf2c->idx, mf2c->col_id);
  c2f = cs_index_map(mc2f->n_rows, mc2f->idx, mc2f->col_id);
  f2f = cs_index_compose(n_faces, f2c, c2f);
  cs_index_sort(f2f);
  h_mat->flag |= CS_SLA_MATRIX_SORTED;

  /* Update index: f2f has the diagonal entry. Remove it for the Hodge index */
  h_mat->idx[0] = 0;
  for (f_id = 0; f_id < n_faces; f_id++)
    h_mat->idx[f_id+1] = h_mat->idx[f_id] + f2f->idx[f_id+1]-f2f->idx[f_id]-1;

  /* Fill column num */
  BFT_MALLOC(h_mat->col_id, h_mat->idx[n_faces], cs_lnum_t);
  shift = 0;
  for (f_id = 0; f_id < n_faces; f_id++)
    for (j = f2f->idx[f_id]; j < f2f->idx[f_id+1]; j++)
      if (f2f->ids[j] != f_id)
        h_mat->col_id[shift++] = f2f->ids[j];

  /* Sanity check */
  assert(shift == h_mat->idx[n_faces]);

  /* Free temporary memory */
  cs_index_free(&f2f);
  cs_index_free(&f2c);
  cs_index_free(&c2f);

  /* Allocate and initialize value array */
  for (f_id = 0; f_id < n_faces; f_id++)
    h_mat->diag[f_id] = 0.0;

  BFT_MALLOC(h_mat->val, h_mat->idx[n_faces], double);
  for (f_id = 0; f_id < h_mat->idx[n_faces]; f_id++)
    h_mat->val[f_id] = 0.0;

  return h_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a _cost_quant_t structure
 *
 * \param[in]  n_max_ent    max number of entities by primal cell
 *
 * \return  a pointer to a new allocated _cost_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static struct _cost_quant_t *
_init_cost_quant(int    n_max_ent)
{
   struct _cost_quant_t *hq = NULL;

  BFT_MALLOC(hq, 1, struct _cost_quant_t);

  hq->invsvol = NULL;
  hq->qmq = NULL;
  hq->T = NULL;
  hq->pq = NULL;
  hq->dq = NULL;

  if (n_max_ent > 0) {

    int  msize = n_max_ent*n_max_ent;
    int  tot_size = n_max_ent + 2*msize;

    /* Allocate invsvol with the total requested size and then reference
       other pointers from this one */
    BFT_MALLOC(hq->invsvol, tot_size, double);
    for (int i = 0; i < tot_size; i++)
      hq->invsvol[i] = 0;

    hq->qmq = hq->invsvol + n_max_ent;
    hq->T = hq->invsvol + n_max_ent + msize;

    BFT_MALLOC(hq->pq, n_max_ent, cs_qvect_t);
    BFT_MALLOC(hq->dq, n_max_ent, cs_qvect_t);

  }

  return hq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_hodge_costq_t structure
 *
 * \param[in]  hb    pointer to the cs_hodge_costq_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static struct _cost_quant_t *
_free_cost_quant(struct _cost_quant_t  *hq)
{
  if (hq == NULL)
    return hq;

  BFT_FREE(hq->invsvol); /* Free in the same time qmq and T */
  BFT_FREE(hq->pq);
  BFT_FREE(hq->dq);

  BFT_FREE(hq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo. and edge/face quantities
 *
 * \param[in]      n_loc_ent  number of local entities
 * \param[in]      matval     values of the tensor related to the material pty
 * \param[in]      pq         pointer to the first set of quantities
 * \param[in]      dq         pointer to the second set of quantities
 * \param[in, out] hq         pointer to a _cost_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cost_quant(int                     n_loc_ent,
                    const cs_real_33_t      matval,
                    const cs_qvect_t       *pq,
                    const cs_qvect_t       *dq,
                    struct _cost_quant_t   *hq)
{
  int  i, j, ii, ij, ji, jj;
  double  dpq, tmp_val;
  cs_real_3_t  mdq_i;

  /* Compute T and qmq matrices */
  for (i = 0; i < n_loc_ent; i++) {

    /* Compute invsvol related to each entity */
    dpq = pq[i].meas*dq[i].meas * _dp3(dq[i].unitv, pq[i].unitv);
    hq->invsvol[i] = 3./dpq; /* 1/subvol where subvol = 1/d * dpq */

    /* Compute diagonal entries */
    _mv3(matval, dq[i].unitv, &mdq_i);
    ii = i*n_loc_ent+i;
    hq->qmq[ii] = dq[i].meas * dq[i].meas * _dp3(dq[i].unitv, mdq_i);
    hq->T[ii] = dpq;

    for (j = i+1; j < n_loc_ent; j++) {

      ij = i*n_loc_ent+j, ji = j*n_loc_ent+i, jj = j*n_loc_ent+j;

      /* Compute qmq (symmetric) */
      tmp_val = dq[j].meas * dq[i].meas * _dp3(dq[j].unitv, mdq_i);
      hq->qmq[ji] = hq->qmq[ij] = tmp_val;

      /* Compute T (not symmetric) */
      hq->T[ij] = dq[i].meas*pq[j].meas * _dp3(dq[i].unitv, pq[j].unitv);
      hq->T[ji] = dq[j].meas*pq[i].meas * _dp3(dq[j].unitv, pq[i].unitv);

    } /* Loop on entities (J) */

  } /* Loop on entities (I) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic COST algo.
 *
 * \param[in]      cid        cell id
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 * \param[in, out] hq         pointer to a _cost_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_using_cost(int                         cid,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  cs_hodge_builder_t         *hb,
                  struct _cost_quant_t       *hq)
{
  int  i, j, k;
  double  invsurf;

  int  n_ent = 0;
  cs_locmat_t  *hloc = hb->hloc;

  const cs_param_hodge_t  h_info = hb->h_info;
  const cs_real_33_t ptyval = {
    {hb->matval[0][0], hb->matval[0][1], hb->matval[0][2]},
    {hb->matval[1][0], hb->matval[1][1], hb->matval[1][2]},
    {hb->matval[2][0], hb->matval[2][1], hb->matval[2][2]}};

  /* Set numbering and geometrical quantities Hodge builder */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    {
      const cs_connect_index_t  *c2e = connect->c2e;

      for (i = c2e->idx[cid]; i < c2e->idx[cid+1]; i++) {

        const cs_lnum_t  e_id = c2e->ids[i];
        const cs_dface_t  fd = quant->dface[i];   /* Dual face quantities */
        const cs_quant_t  ep = quant->edge[e_id]; /* Edge quantities */

        hloc->ids[n_ent] = e_id;

        /* Primal and dual vector quantities are split into
           a measure and a unit vector in order to achieve a better accuracy */
        hq->pq[n_ent].meas = ep.meas;
        hq->dq[n_ent].meas = fd.meas[0] + fd.meas[1];
        invsurf = 1/hq->dq[n_ent].meas;
        for (k = 0; k < 3; k++) {
          hq->dq[n_ent].unitv[k] = invsurf * fd.vect[k];
          hq->pq[n_ent].unitv[k] = ep.unitv[k];
        }
        n_ent++;

      } /* Loop on cell edges */

    }
    break;

  case CS_PARAM_HODGE_TYPE_FPED:
    {
      const cs_sla_matrix_t  *c2f = connect->c2f;

      for (i = c2f->idx[cid]; i < c2f->idx[cid+1]; i++) {

        const cs_lnum_t  f_id = c2f->col_id[i];
        const double  *ed = &(quant->dedge[4*i]); /* Dual edge quantities */
        const cs_quant_t  fp = quant->face[f_id];  /* Face quantities */

        hloc->ids[n_ent] = f_id;

        /* Primal and dual vector quantities are split into
           a measure and a unit vector in order to achieve a better accuracy */
        hq->dq[n_ent].meas = ed[0];
        hq->pq[n_ent].meas = fp.meas;
        for (k = 0; k < 3; k++) {
          hq->pq[n_ent].unitv[k] = fp.unitv[k];
          hq->dq[n_ent].unitv[k] = ed[1+k];
        }
        n_ent++;

      } /* Loop on cell faces */

    }
    break;

  case CS_PARAM_HODGE_TYPE_EDFP:
    {
      const cs_sla_matrix_t  *c2f = connect->c2f;

      for (i = c2f->idx[cid]; i < c2f->idx[cid+1]; i++) {

        const cs_lnum_t  f_id = c2f->col_id[i];
        const short int  sgn = c2f->sgn[i];
        const double  *ed = &(quant->dedge[4*i]);  /* Dual edge quantities */
        const cs_quant_t  fp = quant->face[f_id];  /* Face quantities */

        hloc->ids[n_ent] = f_id;

        /* Primal and dual vector quantities are split into
           a measure and a unit vector in order to achieve a better accuracy */
        hq->dq[n_ent].meas = ed[0];
        hq->pq[n_ent].meas = fp.meas;
        for (k = 0; k < 3; k++) {
          hq->pq[n_ent].unitv[k] = sgn*fp.unitv[k];
          hq->dq[n_ent].unitv[k] = sgn*ed[1+k];
        }
        n_ent++;

      } /* Loop on cell faces */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" This type of discrete Hodge operator is not covered.\n"));

  } /* End of switch */

  assert(n_ent < hloc->n_max_ent + 1);
  hloc->n_ent = n_ent;

  /* Compute additional geometrical quantities: invsvol, qmq and T
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  /* PRIMAL --> DUAL */
  if (h_info.type == CS_PARAM_HODGE_TYPE_FPED ||
      h_info.type == CS_PARAM_HODGE_TYPE_EPFD)
    _compute_cost_quant(n_ent, ptyval, hq->pq, hq->dq, hq);

  /* DUAL --> PRIMAL */
  else if (h_info.type == CS_PARAM_HODGE_TYPE_EDFP)
    _compute_cost_quant(n_ent, ptyval, hq->dq, hq->pq, hq);

  /* Initialize local Hodge matrix */
  for (i = 0; i < n_ent*n_ent; i++)
    hloc->mat[i] = 0.0;

  /* Coefficients related to the value of beta */
  const double  beta = h_info.coef;
  const double  beta2 = beta*beta;
  const double  invcvol = 1 / quant->cell_vol[cid];
  const double  invcvol2 = invcvol*invcvol;
  const double  coef1 = beta*invcvol2;
  const double  coef2 = (1+ 2*beta)*invcvol;
  const double  coef3 = coef2 - 6*beta2*invcvol;
  const double  coef4 = beta2*invcvol;

  /* Add contribution from each sub-volume related to each edge */
  for (k = 0; k < n_ent; k++) { /* Loop over sub-volumes */

    int  kk =  k*n_ent+k;
    double  val_kk = beta * hq->qmq[kk] * hq->invsvol[k];

    for (i = 0; i < n_ent; i++) { /* Loop on cell entities I */

      int  ik = i*n_ent+k;
      double  Tik = hq->T[ik];

      hloc->mat[i*n_ent+i] += Tik*(Tik*val_kk - 2*hq->qmq[ik]);

      for (j = i + 1; j < n_ent; j++) { /* Loop on cell entities J */

        int  jk = j*n_ent+k, ij = i*n_ent+j;
        double  Tjk = hq->T[jk];

        hloc->mat[ij] += Tik*Tjk*val_kk - Tik*hq->qmq[jk] - Tjk*hq->qmq[ik];

      } /* End of loop on J entities */

    } /* End of loop on I entities */

  } /* End of loop on P sub-regions */

  /* Add contribution independent of the sub-region */
  for (i = 0; i < n_ent; i++) {/* Loop on cell entities I */

    int  ii = i*n_ent+i;
    double  miis = hq->qmq[ii]*hq->invsvol[i];

    hloc->mat[ii] = coef1*hloc->mat[ii] + coef3*hq->qmq[ii] + beta2*miis;

    for (j = i + 1; j < n_ent; j++) { /* Loop on cell entities J */

      int  jj = j*n_ent+j, ij = i*n_ent+j, ji = j*n_ent+i;
      double  mjjs = hq->qmq[jj]*hq->invsvol[j];
      double  contrib =  hq->T[ij]*mjjs + hq->T[ji]*miis;

      hloc->mat[ij] = coef1*hloc->mat[ij] + coef2*hq->qmq[ij] - coef4*contrib;
      hloc->mat[ji] = hloc->mat[ij];

    } /* End of loop on J entities */

  } /* End of loop on I entities */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a structure used to compute a discrete Hodge op. when using
 *          WBS algo.
 *
 * \param[in]   n_ent_max    max number of local entities
 * \param[in]   n_entities   number of entities associated to the mesh
 * \param[in]   connect      pointer to a cs_cdo_connect_t struct.
 *
 * \return a pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static struct _wbs_quant_t *
_init_wbs_quant(int                       n_ent_max,
                cs_lnum_t                 n_entities,
                const cs_cdo_connect_t   *connect)
{
  cs_lnum_t  i;

  cs_connect_index_t  *f2e = NULL, *e2v = NULL;
  struct _wbs_quant_t  *hq = NULL;

  const cs_sla_matrix_t  *mf2e = connect->f2e;
  const cs_sla_matrix_t  *me2v = connect->e2v;

  /* Allocate structure */
  BFT_MALLOC(hq, 1, struct _wbs_quant_t);

  /* Weights and cumulators */
  BFT_MALLOC(hq->wf, 3*n_ent_max, double);
  for (i = 0; i < 3*n_ent_max; i++)
    hq->wf[i] = 0;
  hq->wc = hq->wf + n_ent_max;
  hq->cumul = hq->wf + 2*n_ent_max;

  /* Tag */
  BFT_MALLOC(hq->tag, n_entities, short int);
  for (i = 0; i < n_entities; i++)
    hq->tag[i] = -1;

  /* Build f2v connect (both interior and border faces) */
  f2e = cs_index_map(mf2e->n_rows, mf2e->idx, mf2e->col_id);
  e2v = cs_index_map(me2v->n_rows, me2v->idx, me2v->col_id);
  hq->f2v = cs_index_compose(n_entities, f2e, e2v);

  /* Free mapping */
  cs_index_free(&f2e);
  cs_index_free(&e2v);

  return  hq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a structure used to compute a discrete Hodge op. when using
 *          WBS algo.
 *
 * \param[in]   hq     pointer to a _wbs_quant_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static struct _wbs_quant_t *
_free_wbs_quant(struct _wbs_quant_t  *hq)
{
  if (hq == NULL)
    return hq;

  BFT_FREE(hq->wf); /* Deallocate in the same time hvq->wc and hvq->cumul */
  BFT_FREE(hq->tag);
  cs_index_free(&(hq->f2v));

  BFT_FREE(hq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each face a weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *
 * \param[in]      f_id      id of the face
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantites_t structure
 * \param[in, out] hq        pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_get_face_weights(cs_lnum_t                   f_id,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  struct _wbs_quant_t        *hq)
{
  cs_lnum_t  i;
  double  contrib, len;
  cs_real_3_t  un, cp;

  const short int  *loc_ids = hq->tag;
  const cs_quant_t  fq = quant->face[f_id];
  const double  invsurf = 1/fq.meas;

  /* Reset weights */
  for (i = 0; i < connect->n_max_vbyc; i++) hq->wf[i] = 0;

  /* Compute a weight for each vertex of the current face */
  for (i = connect->f2e->idx[f_id]; i < connect->f2e->idx[f_id+1]; i++) {

    cs_lnum_t  e_id = connect->f2e->col_id[i];
    cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
    cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];
    cs_quant_t  eq = quant->edge[e_id];

    _lenunit3(eq.center, fq.center, &len, &un);
    _cp3(un, eq.unitv, &cp);

    contrib = 0.25 * eq.meas * len * _n3(cp) * invsurf;

    hq->wf[loc_ids[v1_id]] += contrib;
    hq->wf[loc_ids[v2_id]] += contrib;

  } /* End of loop on face edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using a conforming algo.
 *          based on the barycentric subdivision of a polyhedron.
 *          This construction is cellwise.
 *          Note: the local matrix is stored inside hb->hloc
 *
 * \param[in]      cid        cell id
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 * \param[in, out] hq         pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_using_wbs(int                         cid,
                 const cs_cdo_connect_t     *connect,
                 const cs_cdo_quantities_t  *quant,
                 cs_hodge_builder_t         *hb,
                 struct _wbs_quant_t        *hq)
{
  int  i, j, l;
  double  contrib, len;
  cs_real_3_t  un;
  cs_lnum_t  ii, jj, v_id;

  cs_locmat_t  *hl = hb->hloc;

  const double  volc = quant->cell_vol[cid];
  const double  ovcell = 1/volc;
  const cs_connect_index_t  *c2v = connect->c2v, *c2e = connect->c2e;
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_lnum_t  vshift = c2v->idx[cid];

  /* Local initializations */
  for (i = 0, j = vshift; j < c2v->idx[cid+1]; i++, j++) {
    v_id = c2v->ids[j];
    hl->ids[i] = v_id;
    hq->tag[v_id] = i;
    hq->wc[i] = ovcell*quant->dcell_vol[j];
  }
  hl->n_ent = i;
  assert(hl->n_ent <= hl->n_max_ent);

  /* Reset local matrix */
  for (i = 0; i < hl->n_ent*hl->n_ent; i++)
    hl->mat[i] = 0.;

  /* Contributions on extra-diag entries resulting from edges.
     The diagonal contribution can be computed more efficiently when one
     considers the contribution of each vertex */
  for (l = c2e->idx[cid]; l < c2e->idx[cid+1]; l++) {

    cs_lnum_t  e_id = c2e->ids[l];
    cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
    cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];
    int  _v1 = hq->tag[v1_id], _v2 = hq->tag[v2_id];

    const cs_quant_t  peq = quant->edge[e_id];
    const cs_dface_t  dfq = quant->dface[l];

    /* Sanity checks */
    assert(_v1 > -1 && _v1 < hl->n_max_ent);
    assert(_v2 > -1 && _v2 < hl->n_max_ent);

    /* Compute pvol_{e,c} */
    contrib  = dfq.meas[0] * _dp3(peq.unitv, &(dfq.unitv[0]));
    contrib += dfq.meas[1] * _dp3(peq.unitv, &(dfq.unitv[3]));

    double  pvol_ec = invdim * peq.meas * contrib;

    i = _v1, j = _v2;
    if (_v1 > _v2)  i = _v2, j = _v1;

    /* Extra-diag contribution:
       e-e: 1/10*1/4 sur |p_ec|
       v-e: 1/20*1/2 sur 1/2|p_ec| for v and v' i.e (x2) => same contrib. as e-e
       e-e + e-v => 1/10 * 1/2 * |p_ec|
     */
    hl->mat[i*hl->n_ent+j] += 0.5*pvol_ec;

  } /* End of loop on cell edges */

  /* Contributions resulting from faces */
  for (l = c2f->idx[cid]; l < c2f->idx[cid+1]; l++) {

    cs_lnum_t  f_id = c2f->col_id[l];
    cs_quant_t  pfq = quant->face[f_id];

    /* Compute pvol_{f,c} */
    _lenunit3(&(quant->cell_centers[3*cid]), pfq.center, &len, &un);

    double  pvol_fc = invdim * len * pfq.meas * _dp3(pfq.unitv, un);

    /* Compute a weight for each vertex of the current face */
    _get_face_weights(f_id, connect, quant, hq);

    /* Compute contributions */
    for (ii = hq->f2v->idx[f_id]; ii < hq->f2v->idx[f_id+1]; ii++) {

      cs_lnum_t  v1_id = hq->f2v->ids[ii];
      int  _v1 = hq->tag[v1_id];
      double  w1f = hq->wf[_v1];

      hq->cumul[_v1] += pvol_fc*w1f;

      /* Diagonal contributions:
         v-f: 1/10*w_v1f^2 sur |pvol_fc|
         e-f: 1/10*w_v1f^2 sur |pvol_fc|
         f-f: 1/10*w_v1f^2 sur |pvol_fc|
         v-f + f-f + e-f => x3
      */
      hl->mat[_v1*hl->n_ent+_v1] += 3*pvol_fc*w1f*w1f;

      for (jj = ii+1; jj < hq->f2v->idx[f_id+1]; jj++) {

        cs_lnum_t  v2_id = hq->f2v->ids[jj];
        int  _v2 = hq->tag[v2_id];
        double  w2f = hq->wf[_v2];

        i = _v1, j = _v2;
        if (_v1 > _v2) i = _v2, j = _v1;

        /* Extra-diag. contributions
           v-f: 1/10*(w_v1f*w_v2f) sur |pvol_fc|
           e-f: 1/10*(w_v1f*w_v2f) sur |pvol_fc|
           f-f: 1/10*(w_v1f*w_v2f) sur |pvol_fc|
         */
        hl->mat[i*hl->n_ent+j] += 3*w1f*w2f*pvol_fc;

      } /* End of loop on face vertices ii != jj */

    } /* End of loop on face vertices */

  } /* End of loop on cell faces */

  /* Other contributions (cell bubble, vertices, edges for diagonal term...) */
  for (i = 0; i < hl->n_ent; i++) {

    /* Diagonal:
       v-v: 1/10*w_vc*|c|
       e-e: 1/10*1/2*w_vc*|c|
       e-v+v-e: 1/20*w_vc*|c|
       v-v + e-e + e-v => 1/10*2*w_vc*|c|

       v-c: 1/10*w_vc^2*|c|
       e-c: 1/10*w_vc^2*|c|
       c-c: 1/10*w_vc^2*|c|
       e-c + v-c + c-c => 1/10*3*w_vc*w_vc*|c|

       f-c: 1/10*Cumul*wc

    */
    contrib = hq->cumul[i]*hq->wc[i];
    hl->mat[i*hl->n_ent+i] += volc*hq->wc[i]*(2 + 3*hq->wc[i]) + contrib;

    /* Extra-diag. contribution
       c-c: 1/10*w1c*w2c sur |c|
       v-c: 1/20*w1c*w2c sur |c| (x2)
       e-c: 1/20*w1c*w2c sur |c| (x2)

       f-c: 1/10* 1/2*(w1c*Cumul_2 + w2c*Cumul_1)
     */
    for (j = i+1; j < hl->n_ent; j++) {

      contrib = 0.5*(hq->cumul[i]*hq->wc[j] + hq->cumul[j]*hq->wc[i]);
      hl->mat[i*hl->n_ent+j] += 3*volc*hq->wc[i]*hq->wc[j] + contrib;

    } /* End of loop vtx_j */

  } /* End of loop vtx_i */

  /* Matrix is symmetric by construction */
  for (i = 0; i < hl->n_ent; i++) {
    hl->mat[i*hl->n_ent+i] *= 0.1; /* => int_tetra Theta_LAG = 1/10 |tet| */
    for (j = i+1; j < hl->n_ent; j++) {
      hl->mat[i*hl->n_ent+j] *= 0.1;
      hl->mat[j*hl->n_ent+i] = hl->mat[i*hl->n_ent+j];
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge op. using the Voronoi algo.
 *
 * \param[in]       cid       cell id
 * \param[in]       connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]       quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in, out]  hb        pointer to a cs_hodge_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_using_voronoi(cs_lnum_t                    c_id,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     cs_hodge_builder_t          *hb)
{
  int  ii, k;
  cs_lnum_t  i;
  double  contrib;
  cs_real_3_t  un, mv;

  cs_locmat_t  *hl = hb->hloc;

  const cs_param_hodge_t  h_info = hb->h_info;
  const cs_real_33_t ptymat = {
    {hb->matval[0][0], hb->matval[0][1], hb->matval[0][2]},
    {hb->matval[1][0], hb->matval[1][1], hb->matval[1][2]},
    {hb->matval[2][0], hb->matval[2][1], hb->matval[2][2]}};

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    {
      const cs_connect_index_t  *c2e = connect->c2e;
      const cs_lnum_t  start = c2e->idx[c_id];
      const cs_lnum_t  end = c2e->idx[c_id+1];

      hl->n_ent = end - start;

      /* Loop on cell edges */
      for (i = start, ii = 0; i < end; i++, ii++) {

        cs_dface_t  dfq = quant->dface[i];
        cs_lnum_t  e_id = c2e->ids[i];
        cs_real_t  len = quant->edge[e_id].meas;

        hl->ids[ii] = e_id;

        /* First sub-triangle contribution */
        _mv3(ptymat, dfq.unitv, &mv);
        contrib = dfq.meas[0] * _dp3(mv, dfq.unitv);
        /* Second sub-triangle contribution */
        _mv3(ptymat, dfq.unitv + 3, &mv);
        contrib += dfq.meas[1] * _dp3(mv, dfq.unitv + 3);

        /* Only a diagonal term */
        hl->mat[ii*hl->n_ent+ii] = contrib/len;

      } /* End of loop on cell edges */

    } /* EpFd */
  case CS_PARAM_HODGE_TYPE_FPED:
    {
      const cs_sla_matrix_t *c2f = connect->c2f;
      const cs_lnum_t  start = c2f->idx[c_id];
      const cs_lnum_t  end = c2f->idx[c_id+1];

      hl->n_ent = end - start;

      for (i = start, ii = 0; i < end; i++, ii++) {

        cs_real_t  len = quant->dedge[4*i];
        cs_lnum_t  f_id = c2f->col_id[i];
        cs_real_t  surf = quant->face[f_id].meas;

        hl->ids[ii] = f_id;

        for (k = 0; k < 3; k++)
          un[k] = quant->dedge[1+4*i+k];
        _mv3(ptymat, un, &mv);

        /* Only a diagonal term */
        hl->mat[ii*hl->n_ent+ii] = len * _dp3(mv, un) / surf;

      } /* End of loop on cell faces */

    } /* FpEd */

  case CS_PARAM_HODGE_TYPE_VPCD:
    {
      const cs_real_t  ptyval = ptymat[0][0]; // Must be isotropic
      const cs_connect_index_t  *c2v = connect->c2v;
      const cs_lnum_t  start = c2v->idx[c_id];
      const cs_lnum_t  end = c2v->idx[c_id+1];

      hl->n_ent = end - start;

      for (i = start, ii = 0; i < end; i++, ii++) {

        cs_lnum_t  v_id = c2v->ids[i];

        hl->ids[ii] = v_id;

        /* Only a diagonal term */
        hl->mat[ii*hl->n_ent+ii] = ptyval * quant->dcell_vol[i];

      } // Loop on cell vertices

    } /* VpCd */
  default:
    break;

  } // End of switch

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_hodge_builder_t structure
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  h_info        algorithm used to build the discrete Hodge op.
 *
 * \return  a new allocated cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_init(const cs_cdo_connect_t   *connect,
                      const cs_time_step_t     *time_step,
                      cs_param_hodge_t          h_info)
{
  cs_hodge_builder_t  *hb = NULL;

  BFT_MALLOC(hb, 1, cs_hodge_builder_t);

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    assert(cs_param_pty_get_type(h_info.pty_id) == CS_PARAM_PTY_ISO);
    hb->n_maxloc_ent = connect->n_max_vbyc;
    hb->n_ent = connect->v_info->n_ent;
    break;
  case CS_PARAM_HODGE_TYPE_EPFD:
    hb->n_maxloc_ent = connect->n_max_ebyc;
    hb->n_ent = connect->e_info->n_ent;
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    hb->n_maxloc_ent = connect->n_max_fbyc;
    hb->n_ent = connect->f_info->n_ent;
    break;
  default:
    hb->n_ent = 0;
    hb->n_maxloc_ent = 0;
    break;

  }

  /* Allocate the local dense matrix storing the coefficient of the local
     discrete Hodge op. associated to a cell */
  hb->hloc = cs_locmat_create(hb->n_maxloc_ent);

  /* Allocate the structure used to stored quantities used during the build
     of the local discrete Hodge op. */
  switch (h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    hb->algoq = _init_cost_quant(hb->n_maxloc_ent);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);
    hb->algoq = _init_wbs_quant(hb->n_maxloc_ent, hb->n_ent, connect);
    break;

  default:
    hb->algoq = NULL;
    break;

  }

  hb->t_cur = time_step->t_cur;
  hb->uniform = cs_param_pty_is_uniform(h_info.pty_id);

  if (hb->uniform) { /* Material property is uniform */

    cs_real_3_t  xyz = {0., 0., 0.};

    cs_evaluate_pty(h_info.pty_id,
                    hb->t_cur,       // When ?
                    xyz,             // Anywhere since uniform
                    h_info.inv_pty,  // Need to inverse the tensor ?
                    &(hb->matval));
  }

  /* Copy h_info */
  hb->h_info.pty_id = h_info.pty_id;
  hb->h_info.inv_pty = h_info.inv_pty;
  hb->h_info.type = h_info.type;
  hb->h_info.algo = h_info.algo;
  hb->h_info.coef = h_info.coef;

  return hb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_hodge_builder_t structure
 *
 * \param[in]  hb    pointer to the cs_hodge_builder_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_free(cs_hodge_builder_t  *hb)
{
  if (hb == NULL)
    return hb;

  hb->hloc = cs_locmat_free(hb->hloc);

  switch (hb->h_info.algo) {
  case CS_PARAM_HODGE_ALGO_COST:
    hb->algoq = _free_cost_quant((struct _cost_quant_t *)hb->algoq);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    hb->algoq = _free_wbs_quant((struct _wbs_quant_t *)hb->algoq);
    break;

  default:
    hb->algoq = NULL;
    break;
  }

  BFT_FREE(hb);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge
 *
 * \param[in]      c_id       cell id
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 *
 * \return a pointer to a cs_locmat_t struct. (local dense matrix)
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_hodge_build_local(int                         c_id,
                     const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     cs_hodge_builder_t         *hb)
{
  /* Sanity check */
  assert(hb != NULL);

  if (!hb->uniform) { /* Material property is not uniform */

    cs_real_3_t  xc;
    for (int k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    cs_evaluate_pty(hb->h_info.pty_id,
                    hb->t_cur,          // When ?
                    xc,                 // Where ?
                    hb->h_info.inv_pty,
                    &(hb->matval));

  }

  switch (hb->h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    _build_using_cost(c_id, connect, quant, hb,
                      (struct _cost_quant_t *)hb->algoq);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    _build_using_wbs(c_id, connect, quant, hb,
                     (struct _wbs_quant_t *)hb->algoq);
    break;

  case CS_PARAM_HODGE_ALGO_VORONOI:
    _build_using_voronoi(c_id, connect, quant, hb);
    break;

  default:
    break;

  }

#if CS_HODGE_DBG > 2
  cs_locmat_dump(c_id, hb->hloc);
#endif

  return hb->hloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a discrete Hodge operator
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]  time_step  pointer to a time step structure
 * \param[in]  h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_compute(const cs_cdo_connect_t      *connect,
                 const cs_cdo_quantities_t   *quant,
                 const cs_time_step_t        *time_step,
                 const cs_param_hodge_t       h_info)
{
  bool  only_diag = (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI) ? true : false;
  cs_sla_matrix_t  *h_mat = NULL;

  /* Allocate and initialize a cs_hodge_builder_t structure */
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, time_step, h_info);

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      h_mat = cs_sla_matrix_create(quant->n_vertices, quant->n_vertices, 1,
                                   CS_SLA_MAT_MSR, false);
    else
      h_mat = _init_hodge_vertex(connect, quant);
    break;

  case CS_PARAM_HODGE_TYPE_EPFD:
    if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      h_mat = cs_sla_matrix_create(quant->n_edges, quant->n_edges, 1,
                                   CS_SLA_MAT_MSR, false);
    else
      h_mat = _init_hodge_edge(connect, quant);
    break;

  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      h_mat = cs_sla_matrix_create(quant->n_faces, quant->n_faces, 1,
                                   CS_SLA_MAT_MSR, false);
    else
      h_mat = _init_hodge_face(connect, quant);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of Hodge operator.\n"));


  } /* End switch */

  /* Fill the matrix related to the discrete Hodge operator.
     Proceed cellwise and then perform an assembly. */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* The local (dense) matrix is stored inside hloc
       n_ent = number of entities related to the current cell */
    cs_hodge_build_local(c_id, connect, quant, hb);

    /* Assemble the local matrix into the system matrix */
    cs_sla_assemble_msr_sym(hb->hloc, h_mat, only_diag);

  } /* End of loop on cells */

  /* Free temporary memory */
  hb = cs_hodge_builder_free(hb);

  return h_mat;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
