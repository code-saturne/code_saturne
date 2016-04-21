/*============================================================================
 * Build discrete Hodge operators
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
#include <math.h>
#include <float.h>
#include <limits.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_sort.h"
#include "cs_timer_stats.h"

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
 * Local Macro definitions
 *============================================================================*/

#define CS_HODGE_DBG 0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/* Main structure used to define a discrete Hodge operator */
struct _hodge_builder_t {

  int         n_maxent_byc;      /* Max local number of entities by primal
                                    cells (use for allocation) */

  cs_param_hodge_t  h_info;   /* Set of parameters related to the discrete
                                 Hodge operator to build. */

  cs_real_33_t      ptymat;   /* Tensor related to the material property.
                                 for EpFd, FpEd, EdFp Hodge op.
                                 Set by default to identity */

  cs_real_t         ptyval;   /* Value related to the material property
                                 for VpCd or CpVd hodge op.
                                 Set by default to unity */

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

  double    *kappa;     /* symmetric dense matrix of size n_ent */
  double    *alpha;     /* dense matrix of size n_ent (not symmetric) */

  cs_real_3_t  *pq;     /* primal geometric quantity (size: n_ent) */
  cs_real_3_t  *dq;     /* dual geometric quantity (size: n_ent) */

};

/* Quantities related to the construction of a local discrete  Hodge op.
   when the Whitney Barycentric Subdivision algo. is employed.
   Used only for Vp --> Cd hodge (up to now) */

struct _wbs_quant_t {

  /* Buffers of size n_max_vbyc */
  double   *w_vol;  /* 1/20 * |p_{e,f}| for each face */
  double   *wvf;    /* weights related to each vertex for a face */

};

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

/* Id related to cs_timer_stats structure used for monitoring */
static int  hodge_ts_id = -1;
static int  hodge_cost_ts_id = -1;
static int  hodge_wbs_ts_id = -1;
static int  hodge_vor_ts_id = -1;

/*! \endcond (end ignore by Doxygen) */

static const double  cs_hodge_wbs_const = 1/60.; // 1/20 * 1/3

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
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t iv = 0; iv < n_vertices; iv++)
    h_mat->diag[iv] = 0.0;

  size_t  nnz = h_mat->idx[n_vertices];
  BFT_MALLOC(h_mat->val, nnz, double);
# pragma omp parallel for if (nnz > CS_THR_MIN)
  for (size_t v = 0; v < nnz; v++)
    h_mat->val[v] = 0.0;

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
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t ie = 0; ie < n_edges; ie++)
    etag[ie] = -1;

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
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t ie = 0; ie < n_edges; ie++)
    h_mat->diag[ie] = 0.0;

  size_t  nnz = h_mat->idx[h_mat->n_rows];
  BFT_MALLOC(h_mat->val, nnz, double);
# pragma omp parallel for if (nnz > CS_THR_MIN)
  for (size_t e = 0; e < nnz; e++)
    h_mat->val[e] = 0.0;

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
#pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t jf = 0; jf < n_faces; jf++)
    h_mat->diag[jf] = 0.0;

  size_t  nnz = h_mat->idx[n_faces];
  BFT_MALLOC(h_mat->val, nnz, double);
# pragma omp parallel for if (nnz > CS_THR_MIN)
  for (size_t f = 0; f < nnz; f++)
    h_mat->val[f] = 0.0;

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

  hq->kappa = NULL;
  hq->alpha = NULL;
  hq->pq = NULL;
  hq->dq = NULL;

  if (n_max_ent > 0) {

    const int  tot_size = n_max_ent*(1 + n_max_ent);

    BFT_MALLOC(hq->kappa, tot_size, double);
    for (int i = 0; i < tot_size; i++)
      hq->kappa[i] = 0;

    hq->alpha = hq->kappa + n_max_ent;

    BFT_MALLOC(hq->pq, n_max_ent, cs_real_3_t);
    BFT_MALLOC(hq->dq, n_max_ent, cs_real_3_t);

  }

  return hq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_hodge_costq_t structure
 *
 * \param[in]  hq    pointer to the cs_hodge_costq_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

static struct _cost_quant_t *
_free_cost_quant(struct _cost_quant_t  *hq)
{
  if (hq == NULL)
    return hq;

  BFT_FREE(hq->kappa); /* Free in the same time alpha */
  BFT_FREE(hq->pq);
  BFT_FREE(hq->dq);

  BFT_FREE(hq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo.
 *          Initialize the local discrete Hodge op. with the consistency part
 *
 * \param[in]      n_ent      number of local entities
 * \param[in]      invcvol    1/|c|
 * \param[in]      ptymat     values of the tensor related to the material pty
 * \param[in]      pq         pointer to the first set of quantities
 * \param[in]      dq         pointer to the second set of quantities
 * \param[in, out] hloc       pointer to a cs_locmat_t struct.
 * \param[in, out] hq         pointer to a _cost_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cost_quant(const int               n_ent,
                    const double            invcvol,
                    const cs_real_33_t      ptymat,
                    const cs_real_3_t      *pq,
                    const cs_real_3_t      *dq,
                    struct _cost_quant_t   *hq,
                    cs_locmat_t            *hloc)
{
  /* Compute several useful quantities
     alpha_ij = delta_ij - pq_j.Consist_i where Consist_i = 1/|c| dq_i
     qmq_ii = dq_i.mat.dq_i
     kappa_i = qmq_ii / |subvol_i|
  */

  cs_real_3_t  mdq_i;

  for (int i = 0; i < n_ent; i++) {

    const double  dsvol_i = _dp3(dq[i], pq[i]);

    double  *alpha_i = hq->alpha + i*n_ent;
    double  *mi = hloc->val + i*n_ent;

    alpha_i[i] = 1 - invcvol * dsvol_i;
    cs_math_33_3_product(ptymat, dq[i], mdq_i);

    const double  qmq_ii = _dp3(dq[i], mdq_i);

    mi[i] = invcvol * qmq_ii;
    hq->kappa[i] = 3. * qmq_ii / dsvol_i;

    for (int j = i+1; j < n_ent; j++) {

      /* Initialize the lower left part of hloc with the consistency part */
      mi[j] = invcvol * _dp3(dq[j], mdq_i);

      /* Compute T (not symmetric) */
      alpha_i[j] = -invcvol * _dp3(pq[j], dq[i]);
      hq->alpha[j*n_ent + i] = -invcvol * _dp3(pq[i], dq[j]);

    } /* Loop on entities (J) */

  } /* Loop on entities (I) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic COST algo.
 *          and cellwise view of the mesh
 *          COST means COnsistency + STabilization
 *
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hb         pointer to a cs_hodge_builder_t structure
 * \param[in, out] hq         pointer to a _cost_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_cellwise_build_with_cost(const cs_cdo_locmesh_t      *lm,
                          cs_hodge_builder_t          *hb,
                          struct _cost_quant_t        *hq)
{
  if (hodge_cost_ts_id > -1)
    cs_timer_stats_start(hodge_cost_ts_id);

  cs_locmat_t  *hloc = hb->hloc;

  const int  n_ent = hloc->n_ent;
  const cs_param_hodge_t  h_info = hb->h_info;

  /* Set numbering and geometrical quantities Hodge builder */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:

    for (int _e = 0; _e < lm->n_ec; _e++) {

      const cs_nvec3_t  dfq = lm->dface[_e];
      const cs_quant_t  peq = lm->edge[_e];

      hloc->ids[_e] = lm->e_ids[_e];
      for (int k = 0; k < 3; k++) {
        hq->dq[_e][k] = dfq.meas * dfq.unitv[k];
        hq->pq[_e][k] = peq.meas * peq.unitv[k];
      }

    } /* Loop on cell edges */
    break;

  case CS_PARAM_HODGE_TYPE_FPED:

    for (int _f = 0; _f < lm->n_fc; _f++) {

      const cs_nvec3_t  deq = lm->dedge[_f];
      const cs_quant_t  pfq = lm->face[_f];

      hloc->ids[_f] = lm->f_ids[_f];
      for (int k = 0; k < 3; k++) {
        hq->pq[_f][k] = pfq.meas * pfq.unitv[k];
        hq->dq[_f][k] = deq.meas * deq.unitv[k];
      }

    } /* Loop on cell faces */
    break;

  case CS_PARAM_HODGE_TYPE_EDFP:

    for (int _f = 0; _f < lm->n_fc; _f++) {

      const short int  sgn = lm->f_sgn[_f];
      const cs_nvec3_t  deq = lm->dedge[_f];
      const cs_quant_t  pfq = lm->face[_f];

      hloc->ids[_f] = lm->f_ids[_f];
      for (int k = 0; k < 3; k++) {
        hq->pq[_f][k] = sgn * deq.meas * deq.unitv[k];
        hq->dq[_f][k] = sgn * pfq.meas * pfq.unitv[k];
      }

    } /* Loop on cell faces */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" This type of discrete Hodge operator is not covered.\n"));

  } /* End of switch */

  /* Sanity checks */
  assert(n_ent < hloc->n_max_ent + 1);

  /* Compute additional geometrical quantities: qmq and T
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  const double  invcvol = 1 / lm->vol_c;
  const double  beta2 = h_info.coef * h_info.coef;

  /* PRIMAL --> DUAL */
  if (h_info.type == CS_PARAM_HODGE_TYPE_FPED ||
      h_info.type == CS_PARAM_HODGE_TYPE_EPFD)
    _compute_cost_quant(n_ent, invcvol,
                        (const cs_real_3_t *)hb->ptymat,
                        (const cs_real_t (*)[3])hq->pq,
                        (const cs_real_t (*)[3])hq->dq,
                        hq, hloc);

  /* DUAL --> PRIMAL */
  else if (h_info.type == CS_PARAM_HODGE_TYPE_EDFP)
    _compute_cost_quant(n_ent, invcvol,
                        (const cs_real_3_t *)hb->ptymat,
                        (const cs_real_t (*)[3])hq->dq,
                        (const cs_real_t (*)[3])hq->pq,
                        hq, hloc);

  double  stab_part;

  for (int i = 0; i < n_ent; i++) { /* Loop on cell entities I */

    const int  shift_i = i*n_ent;
    const double  *alpha_i = hq->alpha + shift_i;
    double  *mi = hloc->val + shift_i;

    /* Add contribution from the stabilization part for
       each sub-volume related to a primal entity */
    stab_part = 0;
    for (int k = 0; k < n_ent; k++) /* Loop over sub-volumes */
      stab_part += hq->kappa[k] * alpha_i[k] * alpha_i[k];

    mi[i] += beta2 * stab_part; // Consistency part has already been computed

    /* Compute extra-diag entries */
    for (int j = i + 1; j < n_ent; j++) { /* Loop on cell entities J */

      const int  shift_j = j*n_ent;
      const double  *alpha_j = hq->alpha + shift_j;
      double  *mj = hloc->val + shift_j;

      /* Add contribution from the stabilization part for
         each sub-volume related to a primal entity */
      stab_part = 0;
      for (int k = 0; k < n_ent; k++) /* Loop over sub-volumes */
        stab_part += hq->kappa[k] * alpha_i[k] * alpha_j[k];

      mi[j] += beta2 * stab_part;
      mj[i] = mi[j]; // Symmetric by construction

    } /* End of loop on J entities */

  } /* End of loop on I entities */

  if (hodge_cost_ts_id > -1)
    cs_timer_stats_stop(hodge_cost_ts_id);
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

  if (hodge_cost_ts_id > -1)
    cs_timer_stats_start(hodge_cost_ts_id);

  int  n_ent = 0;
  cs_locmat_t  *hloc = hb->hloc;

  const cs_param_hodge_t  h_info = hb->h_info;

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
        for (k = 0; k < 3; k++) {
          hq->dq[n_ent][k] = fd.vect[k];
          hq->pq[n_ent][k] = ep.meas * ep.unitv[k];
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
        const cs_nvec3_t  ed = quant->dedge[i]; /* Dual edge quantities */
        const cs_quant_t  fp = quant->face[f_id];  /* Face quantities */

        hloc->ids[n_ent] = f_id;
        for (k = 0; k < 3; k++) {
          hq->pq[n_ent][k] = fp.meas * fp.unitv[k];
          hq->dq[n_ent][k] = ed.meas * ed.unitv[k];
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
        const cs_nvec3_t  ed = quant->dedge[i];    /* Dual edge quantities */
        const cs_quant_t  fp = quant->face[f_id];  /* Face quantities */

        hloc->ids[n_ent] = f_id;
        for (k = 0; k < 3; k++) {
          hq->pq[n_ent][k] = sgn * fp.meas * fp.unitv[k];
          hq->dq[n_ent][k] = sgn * ed.meas * ed.unitv[k];
        }
        n_ent++;

      } /* Loop on cell faces */

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" This type of discrete Hodge operator is not covered.\n"));

  } /* End of switch */

  /* Sanity checks */
  assert(n_ent < hloc->n_max_ent + 1);
  assert(n_ent == hloc->n_ent);

  /* Compute additional geometrical quantities: qmq and T
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  const double  invcvol = 1 / quant->cell_vol[cid];
  const double  beta2 = h_info.coef * h_info.coef;

  /* PRIMAL --> DUAL */
  if (h_info.type == CS_PARAM_HODGE_TYPE_FPED ||
      h_info.type == CS_PARAM_HODGE_TYPE_EPFD)
    _compute_cost_quant(n_ent, invcvol,
                        (const cs_real_3_t *)hb->ptymat,
                        (const cs_real_t (*)[3])hq->pq,
                        (const cs_real_t (*)[3])hq->dq,
                        hq, hloc);

  /* DUAL --> PRIMAL */
  else if (h_info.type == CS_PARAM_HODGE_TYPE_EDFP)
    _compute_cost_quant(n_ent, invcvol,
                        (const cs_real_3_t *)hb->ptymat,
                        (const cs_real_t (*)[3])hq->dq,
                        (const cs_real_t (*)[3])hq->pq,
                        hq, hloc);

  double  stab_part;

  for (i = 0; i < n_ent; i++) { /* Loop on cell entities I */

    const int  shift_i = i*n_ent;
    const double  *alpha_i = hq->alpha + shift_i;
    double  *mi = hloc->val + shift_i;

    /* Add contribution from the stabilization part for
       each sub-volume related to a primal entity */
    stab_part = 0;
    for (k = 0; k < n_ent; k++) /* Loop over sub-volumes */
      stab_part += hq->kappa[k] * alpha_i[k] * alpha_i[k];

    mi[i] += beta2 * stab_part; // Consistency part has already been computed

    /* Compute extra-diag entries */
    for (j = i + 1; j < n_ent; j++) { /* Loop on cell entities J */

      const int  shift_j = j*n_ent;
      const double  *alpha_j = hq->alpha + shift_j;
      double  *mj = hloc->val + shift_j;

      /* Add contribution from the stabilization part for
         each sub-volume related to a primal entity */
      stab_part = 0;
      for (k = 0; k < n_ent; k++) /* Loop over sub-volumes */
        stab_part += hq->kappa[k] * alpha_i[k] * alpha_j[k];

      mi[j] += beta2 * stab_part;
      mj[i] = mi[j]; // Symmetric by construction

    } /* End of loop on J entities */

  } /* End of loop on I entities */

  if (hodge_cost_ts_id > -1)
    cs_timer_stats_stop(hodge_cost_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 * \param[in, out] sloc       pointer to a local matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_stiffness_using_cost(const cs_cdo_locmesh_t     *lm,
                            cs_hodge_builder_t         *hb,
                            cs_locmat_t                *sloc)
{
  if (hodge_cost_ts_id > -1)
    cs_timer_stats_start(hodge_cost_ts_id);

  int  n_ent = 0;
  cs_locmat_t  *hloc = hb->hloc;
  struct _cost_quant_t  *hq = (struct _cost_quant_t *)hb->algoq;

  const cs_param_hodge_t  h_info = hb->h_info;

  /* Set numbering and geometrical quantities Hodge builder */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    {
      n_ent = lm->n_ec;
      hloc->n_ent = lm->n_ec;
      for (int ii = 0; ii < n_ent; ii++) {

        cs_nvec3_t  dfq = lm->dface[ii];
        cs_quant_t  peq = lm->edge[ii];

        for (int k = 0; k < 3; k++) {
          hq->dq[ii][k] = dfq.meas * dfq.unitv[k];
          hq->pq[ii][k] = peq.meas * peq.unitv[k];
        }

      } /* Loop on cell edges */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" This type of discrete Hodge operator is not covered.\n"));

  } /* End of switch */

  /* Compute additional geometrical quantities.
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  const double  invcvol = 1 / lm->vol_c;
  const double  beta2 = h_info.coef * h_info.coef;

  _compute_cost_quant(n_ent, invcvol,
                      (const cs_real_3_t *)hb->ptymat,
                      (const cs_real_t (*)[3])hq->pq,
                      (const cs_real_t (*)[3])hq->dq,
                      hq, hloc);

  double  stab_part;

  for (int ei = 0; ei < n_ent; ei++) { /* Loop on cell entities I */

    const int  shift_i = ei*n_ent;
    const double  *alpha_i = hq->alpha + shift_i;
    const short int  _2ei = 2*ei;
    const short int  i1ei = lm->e2v_sgn[_2ei];
    const short int  i2ei = lm->e2v_sgn[_2ei+1];
    const short int  i1 = lm->e2v_ids[_2ei];
    const short int  i2 = lm->e2v_ids[_2ei+1];
    const double  *hi = hloc->val + shift_i;

    double  *si1 = sloc->val + i1*sloc->n_ent;
    double  *si2 = sloc->val + i2*sloc->n_ent;

    /* Add contribution from the stabilization part for
       each sub-volume related to a primal entity */
    stab_part = 0;
    for (int ek = 0; ek < n_ent; ek++) /* Loop over sub-volumes */
      stab_part += hq->kappa[ek] * alpha_i[ek] * alpha_i[ek];

    /* Diagonal value: consistency part has already been computed */
    const double  dval = hi[ei] + beta2 * stab_part;

    si1[i1] += dval;
    si2[i2] += dval;
    if (i1 < i2)
      si1[i2] -= dval;
    else
      si2[i1] -= dval;

    /* Compute extra-diag entries */
    for (int ej = ei + 1; ej < n_ent; ej++) { /* Loop on cell entities J */

      const int  shift_j = ej*n_ent;
      const double  *alpha_j = hq->alpha + shift_j;
      const short int  _2ej = 2*ej;
      const short int  j1ej = lm->e2v_sgn[_2ej];
      const short int  j2ej = lm->e2v_sgn[_2ej+1];
      const short int  j1 = lm->e2v_ids[_2ej];
      const short int  j2 = lm->e2v_ids[_2ej+1];

      double  *sj1 = sloc->val + j1*sloc->n_ent;
      double  *sj2 = sloc->val + j2*sloc->n_ent;

      /* Add contribution from the stabilization part for
         each sub-volume related to a primal entity */
      stab_part = 0;
      for (int ek = 0; ek < n_ent; ek++) /* Loop over sub-volumes */
        stab_part += hq->kappa[ek] * alpha_i[ek] * alpha_j[ek];

      /* Extra-diagonal value */
      const double xval = hi[ej] + beta2 * stab_part;

      /* Vertex i1 */
      const double  val1 = xval * i1ei;
      if (i1 == j1)
        si1[j1] += 2*val1 * j1ej;
      else if (i1 < j1)
        si1[j1] += val1 * j1ej;
      else
        sj1[i1] += val1 * j1ej;

      if (i1 == j2)
        si1[j2] += 2*val1 * j2ej;
      else if (i1 < j2)
        si1[j2] += val1 * j2ej;
      else
        sj2[i1] += val1 * j2ej;

      /* Vertex i2 */
      const double  val2 = xval * i2ei;
      if (i2 == j1)
        si2[j1] += 2*val2 * j1ej;
      else if (i2 < j1)
        si2[j1] += val2 * j1ej;
      else
        sj1[i2] += val2 * j1ej;

      if (i2 == j2)
        si2[j2] += 2*val2 * j2ej;
      else if (i2 < j2)
        si2[j2] += val2 * j2ej;
      else
        sj2[i2] += val2 * j2ej;

    } /* End of loop on J entities */

  } /* End of loop on I entities */

  /* Stiffness matrix is symmetric by construction */
  for (int ei = 0; ei < sloc->n_ent; ei++) {
    double *si = sloc->val + ei*sloc->n_ent;
    for (int ej = 0; ej < ei; ej++)
      si[ej] = sloc->val[ej*sloc->n_ent + ei];
  }

  if (hodge_cost_ts_id > -1)
    cs_timer_stats_stop(hodge_cost_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a structure used to compute a discrete Hodge op. when using
 *          WBS algo.
 *
 * \param[in]   n_max_vbyc    max number of vertices in a cell
 *
 * \return a pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static struct _wbs_quant_t *
_init_wbs_quant(int    n_max_vbyc)
{
  struct _wbs_quant_t  *hq = NULL;

  /* Allocate structure */
  BFT_MALLOC(hq, 1, struct _wbs_quant_t);

  BFT_MALLOC(hq->w_vol, n_max_vbyc, double);
  BFT_MALLOC(hq->wvf, n_max_vbyc, double);

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

  BFT_FREE(hq->w_vol);
  BFT_FREE(hq->wvf);

  BFT_FREE(hq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each face a weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *         Set also the vertices coordinates
 *
 * \param[in]      f         face id in the cellwise numbering
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      lm        pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hq        pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_wbs_quant(short int                   f,
               const cs_cdo_quantities_t  *quant,
               const cs_cdo_locmesh_t     *lm,
               struct _wbs_quant_t        *hq)
{
  double  len;
  cs_real_3_t  un, cp;

  const cs_real_t  *xc = quant->cell_centers + 3*lm->c_id;
  const cs_quant_t  pfq = lm->face[f];
  const double  ovf = 0.5/pfq.meas;

  /* Reset weights */
  for (short int v = 0; v < lm->n_vc; v++) hq->wvf[v] = 0;

  /* Compute the height of the pyramid of base f */
  cs_math_3_length_unitv(pfq.center, xc, &len, un);

  /* One wants to compute 1/20 * |p_{e,f}| */
  const double  hf = cs_hodge_wbs_const * fabs(len * _dp3(un, pfq.unitv));

  /* Compute a weight for each vertex of the current face */
  for (int i = lm->f2e_idx[f], ii = 0; i < lm->f2e_idx[f+1]; i++, ii++) {

    const short int  e = lm->f2e_ids[i];
    const cs_quant_t  peq = lm->edge[e];
    const short int  v1 = lm->e2v_ids[2*e];
    const short int  v2 = lm->e2v_ids[2*e+1];

    cs_math_3_length_unitv(peq.center, pfq.center, &len, un);
    cs_math_3_cross_product(un, peq.unitv, cp);

    const double  tef = 0.5*peq.meas*len * cs_math_3_norm(cp);

    hq->wvf[v1] += tef;
    hq->wvf[v2] += tef;
    hq->w_vol[ii] = hf * tef;

  } /* End of loop on face edges */

  /* wvf = |dual_cell(v) cap f| / |f| */
  for (short int v = 0; v < lm->n_vc; v++) hq->wvf[v] *= ovf;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic WBS algo.
 *          and cellwise view of the mesh
 *          WBS means Whitney Barycentric Subdivision
 *
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hb         pointer to a cs_hodge_builder_t structure
 * \param[in, out] hq         pointer to a _wbs_quant_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_cellwise_build_with_wbs(const cs_cdo_quantities_t   *quant,
                         const cs_cdo_locmesh_t      *lm,
                         cs_hodge_builder_t          *hb,
                         struct _wbs_quant_t         *hq)
{
  /* Sanity check */
  assert(hb->h_info.algo == CS_PARAM_HODGE_ALGO_WBS);
  assert(hb->h_info.type == CS_PARAM_HODGE_TYPE_VPCD);

  if (hodge_cost_ts_id > -1)
    cs_timer_stats_start(hodge_wbs_ts_id);

  cs_locmat_t  *m = hb->hloc;

  const double  c_coef = 0.1*lm->vol_c;

  /* Initialize the upper part of the local Hodge matrix */
  for (short int vi = 0; vi < lm->n_vc; vi++) {

    double  *mi = m->val + vi*lm->n_vc;
    const double  vi_coef = c_coef * lm->wvc[vi];

    m->ids[vi] = lm->v_ids[vi];
    mi[vi] = vi_coef * lm->wvc[vi];   // Diagonal entry
    for (short int vj = vi+1; vj < lm->n_vc; vj++)
      mi[vj] = vi_coef * lm->wvc[vj]; // Extra-diagonal entries

  } // Loop on cell vertices

  /* Loop on each pef and add the contribution */
  for (short int f = 0; f < lm->n_fc; f++) {

    /* Define useful quantities for WBS algo. */
    _set_wbs_quant(f, quant, lm, hq);

    for (int i = lm->f2e_idx[f], ii = 0; i < lm->f2e_idx[f+1]; i++, ii++) {

      const short int  e = lm->f2e_ids[i];
      const short int  v1 = lm->e2v_ids[2*e];
      const short int  v2 = lm->e2v_ids[2*e+1];

      /* Sanity check */
      assert(v1 > -1 && v2 > -1);

      /* Add local contribution */
      for (short int vi = 0; vi < lm->n_vc; vi++) {

        double  *mi = m->val + vi*lm->n_vc;

        const bool  is_vi = (vi == v1 || vi == v2) ? true : false;

        const double  wic = lm->wvc[vi];
        const double  wif = hq->wvf[vi];

        double  dval = 2*wif*(wif + wic);
        if (is_vi)
          dval += 2*(1 + wic + wif);
        /* 1/20 * |tet| (cf. Rapport HI-A7/7561 in 1991) */
        dval *= hq->w_vol[ii];
        mi[vi] += dval;  /* Add diagonal entry */

        /* Extra-diagonal entries */
        for (short int vj = vi+1; vj < lm->n_vc; vj++) {

          const bool  is_vj = (vj == v1 || vj == v2) ? true : false;

          const double  wjc = lm->wvc[vj];
          const double  wjf = hq->wvf[vj];

          double xval = 2*wif*wjf + wif*wjc + wic*wjf;
          if (is_vi)
            xval += wjf + wjc;
          if (is_vj)
            xval += wif + wic;
          if (is_vi && is_vj)
            xval += 1;
          xval *= hq->w_vol[ii];

          mi[vj] += xval; /* Add extra-diag. entry */

        } // Extra-diag entries

      } // Loop on cell vertices

    } // Loop on face edges

  } // Loop on cell faces

  /* Take into account the value of the associated property */
  if (fabs(hb->ptyval - 1.0) > cs_math_get_machine_epsilon()) {
    for (short int vi = 0; vi < lm->n_vc; vi++) {
      double  *mi = m->val + vi*lm->n_vc;
      for (short int vj = vi; vj < lm->n_vc; vj++)
        mi[vj] *= hb->ptyval;
    }
  }

  /* Local matrix is symmetric by construction. Set the lower part. */
  for (short int vj = 0; vj < lm->n_vc; vj++) {
    double  *mj = m->val + vj*lm->n_vc;
    for (short int vi = vj+1; vi < lm->n_vc; vi++)
      m->val[vi*lm->n_vc + vj] = mj[vi];
  }

  if (hodge_wbs_ts_id > -1)
    cs_timer_stats_stop(hodge_wbs_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 * \param[in, out] sloc       pointer to a local matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_stiffness_using_voronoi(const cs_cdo_locmesh_t    *lm,
                               cs_hodge_builder_t        *hb,
                               cs_locmat_t               *sloc)
{
  cs_real_3_t  mv;

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_start(hodge_vor_ts_id);

  const cs_param_hodge_t  h_info = hb->h_info;

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    {
      /* Loop on cell edges */
      for (int ii = 0; ii < lm->n_ec; ii++) {

        cs_nvec3_t  dfq = lm->dface[ii];
        cs_quant_t  peq = lm->edge[ii];

        cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, dfq.unitv, mv);

        /* Only a diagonal term */
        const double  dval = _dp3(mv, dfq.unitv) * dfq.meas/peq.meas;
        const short int  _2ii = 2*ii;
        const short int  vi = lm->e2v_ids[_2ii];
        const short int  vj = lm->e2v_ids[_2ii+1];
        const short int  sgn_i = lm->e2v_sgn[_2ii];
        const short int  sgn_j = lm->e2v_sgn[_2ii+1];

        double  *si = sloc->val + vi*sloc->n_ent;
        double  *sj = sloc->val + vj*sloc->n_ent;

        si[vi] += dval;
        sj[vj] += dval;
        si[vj] = sj[vi] = dval * sgn_j * sgn_i;

      } /* End of loop on cell edges */

    } /* EpFd */
    break;

  default:
    break;

  } // End of switch

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_stop(hodge_vor_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge op. using the Voronoi algo.
 *
 * \param[in]       lm     pointer to a cs_cdo_locmesh_t structure
 * \param[in, out]  hb     pointer to a cs_hodge_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_cellwise_build_with_voronoi(const cs_cdo_locmesh_t    *lm,
                             cs_hodge_builder_t        *hb)
{
  cs_real_3_t  mv;
  cs_locmat_t  *hmat = hb->hloc;

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_start(hodge_vor_ts_id);

  /* Only a diagonal term by entity to compute */
  switch (hb->h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:

    for (int e = 0; e < lm->n_ec; e++) { /* Loop on cell edges */

      cs_nvec3_t  dfq = lm->dface[e];
      cs_quant_t  peq = lm->edge[e];

      hmat->ids[e] = lm->e_ids[e];
      cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, dfq.unitv, mv);
      hmat->val[e*lm->n_ec+e] = _dp3(mv, dfq.unitv) * dfq.meas/peq.meas;

    }
    break;

  case CS_PARAM_HODGE_TYPE_FPED:

    for (int f = 0; f < lm->n_fc; f++) { /* Loop on cell faces */

      cs_nvec3_t  deq = lm->dedge[f];
      cs_quant_t  pfq = lm->face[f];

      hmat->ids[f] = lm->f_ids[f];
      cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, deq.unitv, mv);
      hmat->val[f*lm->n_fc+f] = _dp3(mv, deq.unitv) * deq.meas/pfq.meas;

    }
    break;

  case CS_PARAM_HODGE_TYPE_VPCD:

    for (int v = 0; v < lm->n_vc; v++) { /* Loop on cell vertices */

      hmat->ids[v] = lm->v_ids[v];
      hmat->val[v*lm->n_vc+v] = hb->ptyval * lm->wvc[v] * lm->vol_c;

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of discrete Hodge operator for the Voronoi algo.");
    break;

  } // End of switch

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_stop(hodge_vor_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge op. using the Voronoi algo.
 *
 * \param[in]       c_id      cell id
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
  int  ii;
  cs_lnum_t  i;
  double  contrib;
  cs_real_3_t  mv;

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_start(hodge_vor_ts_id);

  cs_locmat_t  *hl = hb->hloc;

  const cs_param_hodge_t  h_info = hb->h_info;

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    {
      const cs_connect_index_t  *c2e = connect->c2e;

      /* Loop on cell edges */
      for (i = c2e->idx[c_id], ii = 0; i < c2e->idx[c_id+1]; i++, ii++) {

        cs_dface_t  dfq = quant->dface[i];
        cs_nvec3_t  df0q = dfq.sface[0], df1q = dfq.sface[1];
        cs_lnum_t  e_id = c2e->ids[i];
        cs_real_t  len = quant->edge[e_id].meas;

        hl->ids[ii] = e_id;

        /* First sub-triangle contribution */
        cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, df0q.unitv, mv);
        contrib = df0q.meas * _dp3(mv, df0q.unitv);
        /* Second sub-triangle contribution */
        cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, df1q.unitv, mv);
        contrib += df1q.meas * _dp3(mv, df1q.unitv);

        /* Only a diagonal term */
        hl->val[ii*hl->n_ent+ii] = contrib/len;

      } /* End of loop on cell edges */

    } /* EpFd */
    break;

  case CS_PARAM_HODGE_TYPE_FPED:
    {
      const cs_sla_matrix_t *c2f = connect->c2f;

      for (i = c2f->idx[c_id], ii = 0; i < c2f->idx[c_id+1]; i++, ii++) {

        cs_nvec3_t  deq = quant->dedge[i];
        cs_lnum_t  f_id = c2f->col_id[i];
        cs_real_t  surf = quant->face[f_id].meas;

        hl->ids[ii] = f_id;
        /* Only a diagonal term */
        cs_math_33_3_product((const cs_real_3_t *)hb->ptymat, deq.unitv, mv);
        hl->val[ii*hl->n_ent+ii] = deq.meas * _dp3(mv, deq.unitv) / surf;

      } /* End of loop on cell faces */

    } /* FpEd */
    break;

  case CS_PARAM_HODGE_TYPE_VPCD:
    {
      const cs_connect_index_t  *c2v = connect->c2v;

      for (i = c2v->idx[c_id], ii = 0; i < c2v->idx[c_id+1]; i++, ii++) {

        hl->ids[ii] = c2v->ids[i];
        /* Only a diagonal term */
        hl->val[ii*hl->n_ent+ii] = hb->ptyval * quant->dcell_vol[i];

      } // Loop on cell vertices

    } /* VpCd */
    break;

  default:
    break;

  } // End of switch

  if (hodge_vor_ts_id > -1)
    cs_timer_stats_stop(hodge_vor_ts_id);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize cs_timer_stats_t structure for monitoring purpose
 *
 * \param[in]  level      level of details requested
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_set_timer_stats(int   level)
{
  if (level < 1)
    return;

  /* Timer statistics */
  hodge_ts_id = cs_timer_stats_create("operations", "hodge", "hodge");

  if (level > 1) {
    hodge_cost_ts_id = cs_timer_stats_create("hodge", "hodgeC", "hodgeC");
    hodge_wbs_ts_id = cs_timer_stats_create("hodge", "hodgeW", "hodgeW");
    hodge_vor_ts_id = cs_timer_stats_create("hodge", "hodgeV", "hodgeV");
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_hodge_builder_t structure
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t struct.
 * \param[in]  h_info        algorithm used to build the discrete Hodge op.
 *
 * \return  a new allocated cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_init(const cs_cdo_connect_t   *connect,
                      cs_param_hodge_t          h_info)
{
  cs_hodge_builder_t  *hb = NULL;

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  BFT_MALLOC(hb, 1, cs_hodge_builder_t);

  /* Metadata */
  hb->h_info.inv_pty = h_info.inv_pty;
  hb->h_info.type    = h_info.type;
  hb->h_info.algo    = h_info.algo;
  hb->h_info.coef    = h_info.coef;

  /* Initialize by default the property values */
  hb->ptyval = 1.0; /* for VpCd, CpVd */

  hb->ptymat[0][0] = 1., hb->ptymat[0][1] = hb->ptymat[0][2] = 0.;
  hb->ptymat[1][1] = 1., hb->ptymat[1][0] = hb->ptymat[1][2] = 0.;
  hb->ptymat[2][2] = 1., hb->ptymat[2][1] = hb->ptymat[2][0] = 0.;

  /* Allocate the local dense matrix storing the coefficient of the local
     discrete Hodge op. associated to a cell */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    hb->n_maxent_byc = connect->n_max_vbyc;
    hb->hloc = cs_locmat_create(hb->n_maxent_byc);
    hb->hloc->n_ent = connect->v_info->n_elts;
    break;
  case CS_PARAM_HODGE_TYPE_EPFD:
    hb->n_maxent_byc = connect->n_max_ebyc;
    hb->hloc = cs_locmat_create(hb->n_maxent_byc);
    hb->hloc->n_ent = connect->e_info->n_elts;
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    hb->n_maxent_byc = connect->n_max_fbyc;
    hb->hloc = cs_locmat_create(hb->n_maxent_byc);
    hb->hloc->n_ent = connect->f_info->n_elts;
    break;
  default:
    hb->n_maxent_byc = 0;
    hb->hloc = NULL;
    break;

  }

  /* Allocate the structure used to stored quantities used during the build
     of the local discrete Hodge op. */
  switch (h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    hb->algoq = _init_cost_quant(hb->n_maxent_byc);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);
    hb->algoq = _init_wbs_quant(hb->n_maxent_byc);
    break;

  default:
    hb->algoq = NULL;
    break;

  }

  if (hodge_ts_id > -1)
    cs_timer_stats_stop(hodge_ts_id);

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

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

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

  if (hodge_ts_id > -1)
    cs_timer_stats_stop(hodge_ts_id);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the property attached to a hodge builder
 *
 * \param[in, out]  hb       pointer to a cs_hodge_builder_t structure
 * \param[in]       ptyval   value of the property
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_builder_set_val(cs_hodge_builder_t    *hb,
                         cs_real_t              ptyval)
{
  if (hb == NULL)
    return;

  hb->ptyval = ptyval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the property attached to a hodge builder
 *
 * \param[in, out]  hb       pointer to a cs_hodge_builder_t structure
 * \param[in]       ptymat   values of the tensor related to a property
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_builder_set_tensor(cs_hodge_builder_t     *hb,
                            const cs_real_33_t      ptymat)
{
  if (hb == NULL)
    return;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      hb->ptymat[i][j] = ptymat[i][j];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix from a local discrete Hodge H and
 *          the local discrete gradient and divergence
 *          S = Gloc^t * H * Gloc
 *
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t struct.
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 * \param[in, out] sloc       pointer to a local stiffness matrix struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_build_local_stiffness(const cs_cdo_locmesh_t     *lm,
                               cs_hodge_builder_t         *hb,
                               cs_locmat_t                *sloc)
{
  int  n_ent = 0;

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  /* Sanity check */
  assert(hb != NULL);

  const cs_param_hodge_t  h_info = hb->h_info;

  /* Set n_ent and reset local hodge matrix */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_EPFD:
    n_ent = lm->n_ec;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of discrete Hodge operator for stiffness matrix.");
  }

  hb->hloc->n_ent = n_ent;
  for (int i = 0; i < n_ent*n_ent; i++)
    hb->hloc->val[i] = 0;

  /* Switch according to the requested type of algorithm to use */
  switch (h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    _build_stiffness_using_cost(lm, hb, sloc);
    break;

  case CS_PARAM_HODGE_ALGO_VORONOI:
    _build_stiffness_using_voronoi(lm, hb, sloc);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of algorithm to build a discrete Hodge"
                " operator\n"));
    break;

  }


  if (hodge_ts_id > -1)
    cs_timer_stats_stop(hodge_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge using a cell-wise view of the mesh
 *
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      lm         pointer to a cs_cdo_locmesh_t structure
 * \param[in, out] hb         pointer to a cs_hodge_builder_t structure
 *
 * \return a pointer to a cs_locmat_t struct. (local dense matrix)
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_hodge_build_cellwise(const cs_cdo_quantities_t   *quant,
                        const cs_cdo_locmesh_t      *lm,
                        cs_hodge_builder_t          *hb)
{
  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  /* Sanity check */
  assert(hb != NULL);

  /* Set n_ent and reset local hodge matrix */
  int  n_ent = 0;
  switch (hb->h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    n_ent = lm->n_vc;
    break;
  case CS_PARAM_HODGE_TYPE_EPFD:
    n_ent = lm->n_ec;
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    n_ent = lm->n_fc;
    break;
  case CS_PARAM_HODGE_TYPE_CPVD:
    n_ent = 1;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of discrete Hodge operator.");
  }

  hb->hloc->n_ent = n_ent;
  for (int i = 0; i < n_ent*n_ent; i++)
    hb->hloc->val[i] = 0;

  /* Switch according to the requested type of algorithm to use */
  switch (hb->h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    _cellwise_build_with_cost(lm, hb, (struct _cost_quant_t *)hb->algoq);
    break;

  case CS_PARAM_HODGE_ALGO_VORONOI:
    _cellwise_build_with_voronoi(lm, hb);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    _cellwise_build_with_wbs(quant, lm, hb, (struct _wbs_quant_t *)hb->algoq);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of algorithm to build a discrete Hodge"
                " operator\n"));
    break;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  cs_locmat_dump(lm->c_id, hb->hloc);
#endif

  if (hodge_ts_id > -1)
    cs_timer_stats_stop(hodge_ts_id);

  return hb->hloc;
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
  int  n_ent = 0;

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  /* Sanity check */
  assert(hb != NULL);

  const cs_param_hodge_t  h_info = hb->h_info;

  /* Set n_ent and reset local hodge matrix */
  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    n_ent = connect->c2v->idx[c_id+1] - connect->c2v->idx[c_id];
    break;
  case CS_PARAM_HODGE_TYPE_EPFD:
    n_ent = connect->c2e->idx[c_id+1] - connect->c2e->idx[c_id];
    break;
  case CS_PARAM_HODGE_TYPE_FPED:
  case CS_PARAM_HODGE_TYPE_EDFP:
    n_ent = connect->c2f->idx[c_id+1] - connect->c2f->idx[c_id];
    break;
  case CS_PARAM_HODGE_TYPE_CPVD:
    n_ent = 1;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of discrete Hodge operator.");
  }

  hb->hloc->n_ent = n_ent;
  for (int i = 0; i < n_ent*n_ent; i++)
    hb->hloc->val[i] = 0;

  /* Switch according to the requested type of algorithm to use */
  switch (h_info.algo) {

  case CS_PARAM_HODGE_ALGO_COST:
    _build_using_cost(c_id, connect, quant, hb,
                      (struct _cost_quant_t *)hb->algoq);
    break;

  case CS_PARAM_HODGE_ALGO_VORONOI:
    _build_using_voronoi(c_id, connect, quant, hb);
    break;

  case CS_PARAM_HODGE_ALGO_WBS:
    bft_error(__FILE__, __LINE__, 0,
              _(" Please change your function call with the following one:\n"
                " cs_hodge_build_cellwise(quant, cell_mesh, hodge_builder)"));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of algorithm to build a discrete Hodge"
                " operator\n"));
    break;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  cs_locmat_dump(c_id, hb->hloc);
#endif

  if (hodge_ts_id > -1)
    cs_timer_stats_stop(hodge_ts_id);

  return hb->hloc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the full matrix related to a discrete Hodge operator
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]  pty        pointer to a cs_property_t struct.
 * \param[in]  h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_compute(const cs_cdo_connect_t      *connect,
                 const cs_cdo_quantities_t   *quant,
                 const cs_property_t         *pty,
                 const cs_param_hodge_t       h_info)
{
  bool  only_diag = (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI) ? true : false;
  cs_sla_matrix_t  *h_mat = NULL;
  cs_cdo_locmesh_t  *c_mesh = NULL;
  cs_flag_t  lm_flag = 0;

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  /* Allocate and initialize a cs_hodge_builder_t structure */
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, h_info);

  bool  update_pty = true;
  bool  pty_is_uniform = true;

  if (pty == NULL)
    update_pty = false;
  else
    pty_is_uniform = cs_property_is_uniform(pty);

  switch (h_info.type) {

  case CS_PARAM_HODGE_TYPE_VPCD:
    if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      h_mat = cs_sla_matrix_create(quant->n_vertices, quant->n_vertices, 1,
                                   CS_SLA_MAT_MSR, false);
    else
      h_mat = _init_hodge_vertex(connect, quant);

    if (h_info.algo == CS_PARAM_HODGE_ALGO_WBS) {
      lm_flag = CS_CDO_LOCAL_V | CS_CDO_LOCAL_E | CS_CDO_LOCAL_EV;
      lm_flag |= CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;
      c_mesh = cs_cdo_locmesh_create(connect);
    }
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

    /* Update the value of the property if needed */
    if (update_pty) {
      if (h_info.type == CS_PARAM_HODGE_TYPE_VPCD)
        hb->ptyval = cs_property_get_cell_value(c_id, pty);
      else
        cs_property_get_cell_tensor(c_id, pty, h_info.inv_pty, hb->ptymat);

      if (pty_is_uniform)
        update_pty = false;
    }

    /* The local (dense) matrix is stored inside hloc
       n_ent = number of entities related to the current cell */
    if (h_info.algo == CS_PARAM_HODGE_ALGO_WBS) {

      /* Set the local mesh structure for the current cell */
      cs_cdo_locmesh_build(c_id, lm_flag, connect, quant, c_mesh);
      cs_hodge_build_cellwise(quant, c_mesh, hb);

    }
    else
      cs_hodge_build_local(c_id, connect, quant, hb);

    /* Assemble the local matrix into the system matrix */
    cs_sla_assemble_msr_sym(hb->hloc, h_mat, only_diag);

  } /* End of loop on cells */

  /* Free temporary memory */
  hb = cs_hodge_builder_free(hb);
  if (h_info.algo == CS_PARAM_HODGE_ALGO_WBS)
    cs_cdo_locmesh_free(&c_mesh);

  if (hodge_ts_id > -1)
    cs_timer_stats_start(hodge_ts_id);

  return h_mat;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
