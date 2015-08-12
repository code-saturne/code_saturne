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

#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_sort.h"

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

#define HODGE_DBG 2

/* Macros used to build the flag related to each discrete Hodge operator */
#define  CS_HODGE_INV_PTY  (1 << 0)
#define  CS_HODGE_VPCD     (1 << 1)
#define  CS_HODGE_EPFD     (1 << 2)
#define  CS_HODGE_FPED     (1 << 3)
#define  CS_HODGE_CPVD     (1 << 4)
#define  CS_HODGE_VORONOI  (1 << 5)
#define  CS_HODGE_WHITNEY  (1 << 6)
#define  CS_HODGE_WBS      (1 << 7)
#define  CS_HODGE_COST     (1 << 8)

/* Geometric quantities related to the construction of local discrete
   Hodge op. when the algo. is either COST (included DGA, SUSHI, Generalized
   Crouzeix-Raviart) and the type is  between edges and faces (primal or dual).
*/

struct _hodge_builder_t {

  int         n_max_ent;  /* max number of entities by primal cells
                             (use for allocation) */
  int         n_ent;

  double     *invsvol;    /* size = n_ent*/
  double     *qmq;        /* symmetric dense matrix of size n_ent */
  double     *T;          /* dense matrix of size n_ent (not symmetric) */

  cs_qvect_t   *pq;       /* primal geometric quantity (size: n_ent) */
  cs_qvect_t   *dq;       /* dual geometric quantity (size: n_ent) */

};

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

static const double  invdim = 1./3.;
static const double  _hepfd_tet_coef = 0.05/36;

/*! \endcond (end ignore by Doxygen) */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo. and edge/face quantities
 *
 * \param[in]    matval   values of the tensor related to the material pty
 * \param[in]    pq       pointer to the first set of quantities
 * \param[in]    dq       pointer to the second set of quantities
 * \param[inout] hb       pointer to a _hodge_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_builder(const cs_real_33_t    matval,
                 const cs_qvect_t     *pq,
                 const cs_qvect_t     *dq,
                 cs_hodge_builder_t   *hb)
{
  int  i, j, ii, ij, ji, jj;
  double  dpq, tmp_val;
  cs_real_3_t  mdq_i;

  const int  n_ent = hb->n_ent;

  /* Compute T and qmq matrices */
  for (i = 0; i < n_ent; i++) {

    /* Compute invsvol related to each entity */
    dpq = pq[i].meas*dq[i].meas * _dp3(dq[i].unitv, pq[i].unitv);
    hb->invsvol[i] = 3./dpq; /* 1/subvol where subvol = 1/d * dpq */

    /* Compute diagonal entries */
    _mv3(matval, dq[i].unitv, &mdq_i);
    ii = i*n_ent+i;
    hb->qmq[ii] = dq[i].meas * dq[i].meas * _dp3(dq[i].unitv, mdq_i);
    hb->T[ii] = dpq;

    for (j = i+1; j < n_ent; j++) {

      ij = i*n_ent+j, ji = j*n_ent+i, jj = j*n_ent+j;

      /* Compute qmq (symmetric) */
      tmp_val = dq[j].meas * dq[i].meas * _dp3(dq[j].unitv, mdq_i);
      hb->qmq[ji] = hb->qmq[ij] = tmp_val;

      /* Compute T (not symmetric) */
      hb->T[ij] = dq[i].meas*pq[j].meas * _dp3(dq[i].unitv, pq[j].unitv);
      hb->T[ji] = dq[j].meas*pq[i].meas * _dp3(dq[j].unitv, pq[i].unitv);

    } /* Loop on entities (J) */

  } /* Loop on entities (I) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Prepare the construction of the local hodge related to cell cid
 *          in the case of the COST algo.
 *
 * \param[in]    cid        cell id
 * \param[in]    connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]    quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]    h_info     pointer to a cs_param_hodge_t struct.
 * \param[inout] hloc       pointer to a cs_toolbox_locmat_t struct. to fill
 * \param[inout] hb         pointer to a cs_hodge_builder_t struct. to fill
 */
/*----------------------------------------------------------------------------*/

static void
_prepare_local_hodge_cost(int                          cid,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_param_hodge_t       h_info,
                          cs_toolbox_locmat_t         *hloc,
                          cs_hodge_builder_t          *hb)
{
  int  i, k;
  double  invsurf;
  cs_real_3_t  xyz;

  int  n_ent = 0;

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
        hb->pq[n_ent].meas = ep.meas;
        hb->dq[n_ent].meas = fd.meas[0] + fd.meas[1];
        invsurf = 1/hb->dq[n_ent].meas;
        for (k = 0; k < 3; k++) {
          hb->dq[n_ent].unitv[k] = invsurf * fd.vect[k];
          hb->pq[n_ent].unitv[k] = ep.unitv[k];
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
        hb->dq[n_ent].meas = ed[0];
        hb->pq[n_ent].meas = fp.meas;
        for (k = 0; k < 3; k++) {
          hb->pq[n_ent].unitv[k] = fp.unitv[k];
          hb->dq[n_ent].unitv[k] = ed[1+k];
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
        hb->dq[n_ent].meas = ed[0];
        hb->pq[n_ent].meas = fp.meas;
        for (k = 0; k < 3; k++) {
          hb->pq[n_ent].unitv[k] = sgn*fp.unitv[k];
          hb->dq[n_ent].unitv[k] = sgn*ed[1+k];
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
  hb->n_ent = n_ent;

  /* Compute additional geometrical quantities: invsvol, qmq and T
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  /* Get the material property for this cell */
  for (k = 0; k < 3; k++)
    xyz[k] = quant->cell_centers[3*cid+k];

  cs_real_33_t  matval;
  cs_param_pty_get_val(h_info.pty_id,
                       -1.0, // when ?
                       xyz,  // where ?
                       h_info.inv_pty,
                       &matval);

  const cs_real_33_t ptyval = {{matval[0][0], matval[0][1], matval[0][2]},
                               {matval[1][0], matval[1][1], matval[1][2]},
                               {matval[2][0], matval[2][1], matval[2][2]}};
  /* PRIMAL --> DUAL */
  if (h_info.type == CS_PARAM_HODGE_TYPE_FPED ||
      h_info.type == CS_PARAM_HODGE_TYPE_EPFD)
    _compute_builder(ptyval, hb->pq, hb->dq, hb);

  /* DUAL --> PRIMAL */
  else if (h_info.type == CS_PARAM_HODGE_TYPE_EDFP)
    _compute_builder(ptyval, hb->dq, hb->pq, hb);

  /* Initialize local Hodge matrix */
  for (i = 0; i < n_ent*n_ent; i++)
    hloc->mat[i] = 0.0;

}

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
 * \brief   Assemble a local discrete Hodge operator
 *
 * \param[in]     hloc       pointer to a cs_toolbox_locmat_t struct.
 * \param[inout]  h_mat    pointer to a cs_sla_matrix_t
 */
/*----------------------------------------------------------------------------*/

static void
_assemble_local_hodge(const cs_toolbox_locmat_t  *hloc,
                      cs_sla_matrix_t            *h_mat)
{
  int  i, j, k;

  const int  n_ent = hloc->n_ent;

  /* Assemble local H contribution into H */
  for (i = 0; i < n_ent; i++) {

    int  iid = hloc->ids[i];
    cs_lnum_t  start = h_mat->idx[iid];

    /* Add diagonal term : H(i,i) */
    h_mat->diag[iid] += hloc->mat[i*n_ent+i];

    for (j = i + 1; j < n_ent; j++) {

      int  ij = i*n_ent+j;

      if (fabs(hloc->mat[ij]) > cs_base_zthreshold) { /* Not zero */

        int  jid = hloc->ids[j];

        /* First add: H(i,j) */
        for (k = start;
             k < h_mat->idx[iid+1] && h_mat->col_id[k] != jid; k++);
        assert(k < h_mat->idx[iid+1]);
        h_mat->val[k] += hloc->mat[ij];
        start = k;

        /* Second add: H(j,i) */
        for (k = h_mat->idx[jid];
             k < h_mat->idx[jid+1] && h_mat->col_id[k] != iid; k++);
        assert(h_mat->col_id[k] == iid);
        h_mat->val[k] += hloc->mat[ij];

      } /* Hloc[ij] != 0.0 */

    } /* Loop on j (Add extra-diag terms) */

  } /* Loop on i (Add diagonal term) */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize by default the matrix related to a discrete
 *          Hodge operator
 *
 * \param[in]    connect  pointer to a cs_cdo_connect_t structure
 * \param[in]    quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]    type     type of the discrete Hodge op. to initiliaze
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_init_matrix(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     const cs_param_hodge_type_t  type)
{
  cs_sla_matrix_t  *h_mat = NULL;

  switch (type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    h_mat = _init_hodge_vertex(connect, quant);
    break;

  case CS_PARAM_HODGE_TYPE_EPFD:
    h_mat = _init_hodge_edge(connect, quant);
    break;

  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    h_mat = _init_hodge_face(connect, quant);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of discrete Hodge operator\n"));

  } /* Hodge type */

  return h_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_hodge_builder_t structure
 *
 * \param[in]  n_max_ent    max number of entities by primal cell
 *
 * \return  a new allocated cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_init(int   n_max_ent)
{
  cs_hodge_builder_t  *hb = NULL;

  BFT_MALLOC(hb, 1, cs_hodge_builder_t);

  hb->n_max_ent = n_max_ent;
  hb->n_ent = 0;
  hb->invsvol = NULL;
  hb->qmq = NULL;
  hb->T = NULL;
  hb->pq = NULL;
  hb->dq = NULL;

  if (n_max_ent > 0) {

    int  i;
    int  msize = n_max_ent*n_max_ent;
    int  tot_size = n_max_ent + 2*msize;

    /* Allocate invsvol with the total requested size and then reference
       other pointers from this one */
    BFT_MALLOC(hb->invsvol, tot_size, double);
    for (i = 0; i < tot_size; i++)
      hb->invsvol[i] = 0;

    hb->qmq = hb->invsvol + n_max_ent;
    hb->T = hb->invsvol + n_max_ent + msize;

    BFT_MALLOC(hb->pq, n_max_ent, cs_qvect_t);
    BFT_MALLOC(hb->dq, n_max_ent, cs_qvect_t);

  }

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

  if (hb->n_max_ent > 0) {
    BFT_FREE(hb->invsvol); /* Free in the same time qmq and T */
    BFT_FREE(hb->pq);
    BFT_FREE(hb->dq);
  }

  BFT_FREE(hb);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic COST algo.
 *
 * \param[in]     cid        cell id
 * \param[in]     connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]     quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]     h_info     pointer to a cs_param_hodge_t struct.
 * \param[in,out] hl         pointer to a cs_toolbox_locmat_t struct.
 * \param[in,out] hb         pointer to a cs_hodge_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_cost_build_local(int                         cid,
                          const cs_cdo_connect_t     *connect,
                          const cs_cdo_quantities_t  *quant,
                          const cs_param_hodge_t      h_info,
                          cs_toolbox_locmat_t        *hl,
                          cs_hodge_builder_t         *hb)
{
  int  i, j, k;
  double  val_kk, contrib;

  const double  invcvol = 1 / quant->cell_vol[cid];
  const double  invcvol2 = invcvol*invcvol;

  /* Coefficients related to the value of beta */
  const double  beta = h_info.coef;
  const double  beta2 = beta*beta;
  const double  coef1 = beta*invcvol2;
  const double  coef2 = (1+ 2*beta)*invcvol;
  const double  coef3 = coef2 - 6*beta2*invcvol;
  const double  coef4 = beta2*invcvol;

  /* Sanity checks */
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  /* Fill Hodge builder and initialize local Hodge matrix */
  _prepare_local_hodge_cost(cid, connect, quant, h_info, hl, hb);

  /* Add contribution from each sub-volume related to each edge */
  for (k = 0; k < hl->n_ent; k++) { /* Loop over sub-volumes */

    int  kk =  k*hl->n_ent+k;

    for (i = 0; i < hl->n_ent; i++) { /* Loop on cell entities I */

      int  ik = i*hl->n_ent+k;
      double  Tik = hb->T[ik];

      val_kk = beta * hb->qmq[kk] * hb->invsvol[k];
      contrib = Tik*val_kk - 2*hb->qmq[ik];
      hl->mat[i*hl->n_ent+i] += Tik*contrib;

      for (j = i + 1; j < hl->n_ent; j++) { /* Loop on cell entities J */

        int  jk = j*hl->n_ent+k, ij = i*hl->n_ent+j;
        double  Tjk = hb->T[jk];

        hl->mat[ij] += Tik*Tjk*val_kk - Tik*hb->qmq[jk] - Tjk*hb->qmq[ik];

      } /* End of loop on J entities */

    } /* End of loop on I entities */

  } /* End of loop on P sub-regions */

  /* Add contribution independent of the sub-region */
  for (i = 0; i < hl->n_ent; i++) {/* Loop on cell entities I */

    int  ii = i*hl->n_ent+i;
    double  mii = hb->qmq[ii];
    double  miis = mii*hb->invsvol[i];

    hl->mat[ii] = coef1*hl->mat[ii] + coef3*mii + beta2*miis;

    for (j = i + 1; j < hl->n_ent; j++) { /* Loop on cell entities J */

      int  jj = j*hl->n_ent+j, ij = i*hl->n_ent+j, ji = j*hl->n_ent+i;
      double mjjs = hb->qmq[jj]*hb->invsvol[j];

      contrib = hb->T[ij]*mjjs + hb->T[ji]*miis;
      hl->mat[ij] = coef1*hl->mat[ij] + coef2*hb->qmq[ij] - coef4*contrib;
      hl->mat[ji] = hl->mat[ij];

    } /* End of loop on J entities */

  } /* End of loop on I entities */

#if HODGE_DBG > 2
  cs_toolbox_locmat_dump(cid, hl);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a discrete Hodge operator using a generic reconstruction
 *          algorithm: Reconstruction with Consistent and Coercive part
 *          stabilized by Beta (COST)
 *          A call to cs_hodge_builder_init() should have been done before.
 *
 * \param[in]    connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]    quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]    h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_cost_build(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant,
                    const cs_param_hodge_t       h_info)
{
  cs_lnum_t  c_id;

  cs_toolbox_locmat_t  *hloc = NULL;

  /* Sanity checks */
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  /* Build a local hodge matrix */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    hloc = cs_toolbox_locmat_create(connect->n_max_ebyc);
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    hloc = cs_toolbox_locmat_create(connect->n_max_fbyc);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of Hodge operator.\n"));

  } /* End switch */

  /* Initialize Hodge matrix structure */
  cs_sla_matrix_t  *h_mat = cs_hodge_init_matrix(connect, quant, h_info.type);

  /* Initialize temporary structures to build the Hodge operator */
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(hloc->n_max_ent);

  /* Fill H matrix */
  for (c_id = 0; c_id < quant->n_cells; c_id++) { /* Loop on cells */

    /* The local (dense) matrix is stored inside hloc
       n_ent = number of entities related to the current cell */
    cs_hodge_cost_build_local(c_id, connect, quant, h_info, hloc, hb);
    _assemble_local_hodge(hloc, h_mat);

  } /* End of loop on cells */

  /* Free temporary memory */
  hb = cs_hodge_builder_free(hb);
  hloc = cs_toolbox_locmat_free(hloc);

  return h_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   H is an operator from primal edges to dual faces. It mimics a Hodge
 *          operator from primal mesh to dual mesh.
 *          Use voronoi algorithm.
 *
 * \param[in]   connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]   quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]    h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_voronoi_build(const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_param_hodge_t       h_info)
{
  cs_lnum_t  c_id;
  int  entid, i, k;
  double  meas, contrib;
  cs_real_3_t  dv, xyz, mv, xc;
  cs_real_33_t  matval;

  cs_sla_matrix_t  *h_mat = NULL;

  /* Sanity checks */
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);

  /* Build a local hodge matrix */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    h_mat = cs_sla_matrix_create(quant->n_edges, quant->n_edges, 1,
                                 CS_SLA_MAT_MSR, false);
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
    h_mat = cs_sla_matrix_create(quant->n_faces, quant->n_faces, 1,
                                 CS_SLA_MAT_MSR, false);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of Hodge operator.\n"));

  } /* End switch */

  if (cs_param_pty_is_uniform(h_info.pty_id))
    cs_param_pty_get_val(h_info.pty_id,
                         -1.0, // when ?
                         xyz,  // where ?
                         h_info.inv_pty,
                         &matval);

  /* Fill H */
  for (c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Get material at the center of each cell in case of non-uniformity */
    if (!cs_param_pty_is_uniform(h_info.pty_id)) {
      for (k = 0; k < 3; k++)
        xc[k] = quant->cell_centers[3*c_id+k];
      cs_param_pty_get_val(h_info.pty_id,
                           -1.0, // when ?
                           xyz,  // where ?
                           h_info.inv_pty,
                           &matval);
    }

    const cs_real_33_t ptyval = {{matval[0][0], matval[0][1], matval[0][2]},
                                 {matval[1][0], matval[1][1], matval[1][2]},
                                 {matval[2][0], matval[2][1], matval[2][2]}};

    if (h_info.type == CS_PARAM_HODGE_TYPE_EPFD) { /* Loop on cell edges */

      const cs_connect_index_t  *c2e = connect->c2e;

      for (i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; i++) {

        cs_dface_t  dfq = quant->dface[i];

        entid = c2e->ids[i];  /* edge id in full system */

        /* First sub-triangle contribution */
        _mv3(ptyval, dfq.unitv, &mv);
        contrib = dfq.meas[0] * _dp3(mv, dfq.unitv);
        /* Second sub-triangle contribution */
        _mv3(ptyval, dfq.unitv + 3, &mv);
        contrib += dfq.meas[1] * _dp3(mv, dfq.unitv + 3);

        h_mat->diag[entid] += contrib;

      } /* End of loop on cell edges */

    }
    else { /* FPED */

      const cs_lnum_t  *idx = connect->c2f->idx;

      for (i = idx[c_id]; i < idx[c_id+1]; i++) {

        meas = quant->dedge[3*i];
        for (k = 0; k < 3; k++)
          dv[k] = quant->dedge[1+3*i+k];
        _mv3(ptyval, dv, &mv);

        entid = connect->c2f->col_id[i];
        h_mat->diag[entid] += meas * _dp3(mv, dv);

      } /* End of loop on cell faces */

    } /* FPED */

  } /* End of loop on cells */

  if (h_info.type == CS_PARAM_HODGE_TYPE_EPFD) { /* Loop on cell edges */
    for (entid = 0; entid < quant->n_edges; entid++)
      h_mat->diag[entid] /= quant->edge[entid].meas;
  }
  else { /* Loop on cell faces */
    for (entid = 0; entid < quant->n_faces; entid++)
      h_mat->diag[entid] /= quant->face[entid].meas;
  }

  return h_mat;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
