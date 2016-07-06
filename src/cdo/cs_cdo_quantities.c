/*============================================================================
 * Manage geometrical quantities needed in CDO schemes
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

#include <limits.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  CDO_QUANTITIES_DBG 0  /* Switch off/on debug information */

/* Redefined names of function from cs_math to get shorter names */
#define _n3  cs_math_3_norm
#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/* Temporary structures to build mesh quantities */

typedef struct {

  int         XYZ[3]; /* Direct permutation of the ref. axis such that nZ
                         is maximal */
  cs_nvec3_t  q;      /* face surface and its unit normal */
  double      omega;  /* P = Point belonging to the face omega = - < n, P> */

} _cdo_fspec_t;

typedef struct { /* Face sub-quantities */

  double  F1;
  double  Fa;
  double  Fb;
  double  Fc;
  double  Fa2;
  double  Fb2;
  double  Fc2;

} _cdo_fsubq_t;

typedef struct { /* These quantities are the integral of q on the plane
                    (alpha, beta) where the face is projected */

  double  p1;     /* q = 1 */
  double  pa;     /* q = alpha */
  double  pb;     /* q = beta */
  double  pc;     /* q = gamma */
  double  pab;    /* q = alpha * beta */
  double  pa2;    /* q = alpha^2 */
  double  pb2;    /* q = beta^2 */

} _cdo_projq_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static const double one_12 = 1/12.;
static const double one_24 = 1/24.;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Several private functions for volume and centroid computation
 * ---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Define an unitary normal to the current face
 * Compute omega = - <n, P> where P belongs to the face
 * Choose projection axis in order to maximize the projected area
 * Define a direct basis (alpha, beta, gamma) with this choice
 * ---------------------------------------------------------------------------*/

static _cdo_fspec_t
_get_fspec(cs_lnum_t                    f_id,
           const cs_mesh_t             *m,
           const cs_mesh_quantities_t  *mq)
{
  cs_lnum_t  f, j, k, v, e, s;
  double  inv_n, nx, ny, nz;
  double  P[3]; /* Point belonging to the current face */
  _cdo_fspec_t  fspec;

  const int X = 0, Y = 1, Z = 2;

  /* Treatment according to the kind of face (interior or border) */

  if (f_id < m->n_i_faces) { /* Interior face */

    /* Choose a vertex belonging to this face */
    f = f_id;
    s = m->i_face_vtx_idx[f];
    e = m->i_face_vtx_idx[f+1];
    inv_n = 1.0 / (e - s);

    for (k = 0; k < 3; k++)
      P[k] = 0.0;

    for (j = s; j < e; j++) {
      v = m->i_face_vtx_lst[j];
      for (k = 0; k < 3; k++)
        P[k] += m->vtx_coord[3*v+k];
    }

    for (k = 0; k < 3; k++)
      P[k] *= inv_n;

    cs_nvec3(&(mq->i_face_normal[3*f]), &(fspec.q));

  }
  else { /* Border face */

    /* Choose a vertex belonging to this face */
    f = f_id - m->n_i_faces;
    s = m->b_face_vtx_idx[f];
    e = m->b_face_vtx_idx[f+1];
    inv_n = 1.0 / (e - s);

    for (k = 0; k < 3; k++)
      P[k] = 0.0;

    for (j = s; j < e; j++) {
      v = m->b_face_vtx_lst[j];
      for (k = 0; k < 3; k++)
        P[k] += m->vtx_coord[3*v+k];
    }

    for (k = 0; k < 3; k++)
      P[k] *= inv_n;

    cs_nvec3(&(mq->b_face_normal[3*f]), &(fspec.q));

  }

  /* Define omega = -<n,P>*/
  fspec.omega = - _dp3(P, fspec.q.unitv);

  /* Define a direct basis such that n[Z] maximal */
  nx = fabs(fspec.q.unitv[X]);
  ny = fabs(fspec.q.unitv[Y]);
  nz = fabs(fspec.q.unitv[Z]);

  if (nx > ny && nx > nz)
    fspec.XYZ[Z] = X;
  else
    fspec.XYZ[Z] = (ny > nz) ? Y : Z;

  fspec.XYZ[X] = (fspec.XYZ[Z] + 1) % 3;
  fspec.XYZ[Y] = (fspec.XYZ[X] + 1) % 3;

#if CDO_QUANTITIES_DBG > 1
  printf("\n F: (%d, %d) >> surf: %e; omega: %e; XYZ: %d%d%d; [%e, %e, %e]\n",
         f_id, f, fspec.q.meas, fspec.omega,
         fspec.XYZ[0], fspec.XYZ[1], fspec.XYZ[2],
         fspec.q.unitv[0], fspec.q.unitv[1], fspec.q.unitv[2]);
#endif

  return fspec;
}

/* ---------------------------------------------------------------------------*
 * Compute projected integrals and quantities
 * ---------------------------------------------------------------------------*/

static _cdo_projq_t
_get_proj_quantities(cs_lnum_t                f_id,
                     const cs_cdo_connect_t  *connect,
                     const cs_real_t          coords[],
                     const int                axis[])
{
  cs_lnum_t  v_id[2];
  cs_real_t  a[2], b[2], a2[2], b2[2];

  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  /* Initialize structure */
  _cdo_projq_t  projq;

  /* These quantities are the integral of q on the plane
     (alpha, beta) where the face is projected */

  projq.p1 = projq.pa = projq.pb =  projq.pc = 0.0;
  projq.pab = projq.pa2 = projq.pb2 = 0.0;

  /* Scan edges which belong to the current face */
  for (cs_lnum_t  i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    short int  e_sgn = f2e->sgn[i];
    cs_lnum_t  e_id = f2e->col_id[i];
    cs_lnum_t  s = e2v->idx[e_id];

    if (e_sgn > 0)
      v_id[0] = e2v->col_id[s], v_id[1] = e2v->col_id[s+1];
    else
      v_id[0] = e2v->col_id[s+1], v_id[1] = e2v->col_id[s];

    /* Vertices in the plane (alpha, beta) */
    for (int k = 0; k < 2; k++) {
      a[k] = coords[3*v_id[k] + axis[0]];
      a2[k] = a[k]*a[k];
      b[k] = coords[3*v_id[k] + axis[1]];
      b2[k] = b[k]*b[k];
    }

    /* Related variables */
    cs_real_t  a0_3 = a2[0] * a[0];
    cs_real_t  b0_3 = b2[0] * b[0];
    cs_real_t  da  = a[1] - a[0], db  = b[1] - b[0];
    cs_real_t  C1  = a[0] + a[1];
    cs_real_t  Ca  = C1 * a[1] + a2[0];
    cs_real_t  Cb  = b2[1] + b[1]*b[0] + b2[0];
    cs_real_t  Ca2 = a[1] * Ca + a0_3;
    cs_real_t  Cb2 = b[1] * Cb + b0_3;
    cs_real_t  Cab = 3*a2[1] + 2*a[1]*a[0] + a2[0];
    cs_real_t  Kab = a2[1] + 2*a[1]*a[0] + 3*a2[0];

    projq.p1  += db * C1;
    projq.pa  += db * Ca;
    projq.pb  += da * Cb;
    projq.pa2 += db * Ca2;
    projq.pb2 += da * Cb2;
    projq.pab += db * (b[1] * Cab + b[0] * Kab);

  } /* Loop on face edges */

  projq.p1  *=  0.5;
  projq.pa  *=  cs_math_onesix;
  projq.pb  *= -cs_math_onesix;
  projq.pab *=  one_24;
  projq.pa2 *=  one_12;
  projq.pb2 *= -one_12;

  return  projq;
}

/* ---------------------------------------------------------------------------*/

static _cdo_fsubq_t
_get_fsub_quantities(cs_lnum_t                 f_id,
                     const cs_cdo_connect_t   *connect,
                     const cs_real_t          *coord,
                     _cdo_fspec_t              fspec)
{
  _cdo_fsubq_t  fsubq;

  double  na = fspec.q.unitv[fspec.XYZ[0]];
  double  nb = fspec.q.unitv[fspec.XYZ[1]];
  double  nc = fspec.q.unitv[fspec.XYZ[2]];
  double  k1 = 1./nc;
  double  k2 = k1 * k1;
  double  k3 = k2 * k1;

  /* Compute projected quantities */
  _cdo_projq_t  projq = _get_proj_quantities(f_id, connect, coord, fspec.XYZ);

#if CDO_QUANTITIES_DBG > 1
  printf(" F: %d >> p1: %.4e, pa: %.4e, pb: %.4e, pc: %.4e,"
         " pab: %.4e, pa2: %.4e, pb2: %.4e\n",
         f_id, projq.p1, projq.pa, projq.pb, projq.pc,
         projq.pab, projq.pa2, projq.pb2);
#endif

  /* Compute face sub-quantities */
  fsubq.F1 = k1*projq.p1;
  fsubq.Fa = k1 * projq.pa;
  fsubq.Fb = k1 * projq.pb;
  fsubq.Fc = -k2 * (projq.pa * na + projq.pb * nb + fspec.omega * projq.p1);

  fsubq.Fa2 = k1 * projq.pa2;
  fsubq.Fb2 = k1 * projq.pb2;
  fsubq.Fc2 = k3 * (na*na * projq.pa2 + 2*na*nb * projq.pab +
                    nb*nb * projq.pb2 +
                    fspec.omega * (2*na * projq.pa +
                                   2*nb * projq.pb +
                                   fspec.omega * projq.p1));

  return fsubq;
}

/*----------------------------------------------------------------------------
 * Build edge centers and edge vectors
 * ---------------------------------------------------------------------------*/

static void
_compute_edge_quantities(const cs_mesh_t         *mesh,
                         const cs_cdo_connect_t  *topo,
                         cs_cdo_quantities_t     *iq)  /* In/out */
{
  int  i, j, k;
  int  v_id[2];
  double  xaxb[3];
  double  xva, xvb, len, invlen;
  cs_quant_t  eq;

  const int n_edges = iq->n_edges;

  /* Sanity check */
  assert(topo->e2v != NULL);

  /* Build edge centers and edge vectors */
  BFT_MALLOC(iq->edge, n_edges, cs_quant_t);

  for (i = 0; i < n_edges; i++) {

    /* Get the two vertex ids related to the current edge */
    for (j = topo->e2v->idx[i], k = 0; j < topo->e2v->idx[i+1]; j++, k++)
      v_id[k] = topo->e2v->col_id[j];
    assert(k == 2);

    for (k = 0; k < 3; k++) {
      xva = mesh->vtx_coord[3*v_id[0]+k];
      xvb = mesh->vtx_coord[3*v_id[1]+k];
      xaxb[k] = xvb - xva;
      eq.center[k] = 0.5 * ( xva + xvb );
    }
    len = _n3(xaxb);
    assert(len > 0);
    invlen = 1/len;
    eq.meas = len;

    if (v_id[1] > v_id[0]) /* vb > va */
      for (k = 0; k < 3; k++)
        eq.unitv[k] = invlen * xaxb[k];
    else  /* Change orientation */
      for (k = 0; k < 3; k++)
        eq.unitv[k] = -invlen * xaxb[k];

    iq->edge[i] = eq;

  } /* End of loop on edges */

}

/*----------------------------------------------------------------------------
 * Compute dual cell volumes related to primal vertices
 * Storage based on c2v connectivity
 * ---------------------------------------------------------------------------*/

static void
_compute_dcell_quantities(const cs_cdo_connect_t  *topo,
                          cs_cdo_quantities_t     *quant)    /* In/out */
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_connect_index_t  *c2v = topo->c2v;
  const cs_sla_matrix_t  *c2f = topo->c2f;
  const cs_sla_matrix_t  *f2e = topo->f2e;

  /* Compute part of dual volume related to each primal cell */
  const cs_lnum_t  c2v_idx_size = c2v->idx[quant->n_cells];
  BFT_MALLOC(quant->dcell_vol, c2v_idx_size, double);
# pragma omp parallel for if (c2v_idx_size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < c2v_idx_size; i++)
    quant->dcell_vol[i] = 0.0;

  /* Link between the mesh numbering and the cellwise numbering */
  short int  *vtag = NULL;
  BFT_MALLOC(vtag, quant->n_vertices, short int);
# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
    vtag[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  *xc = quant->cell_centers + 3*c_id;
    double  *v_vol = quant->dcell_vol + c2v->idx[c_id];

    /* Define vtag */
    for (cs_lnum_t i = c2v->idx[c_id], ii = 0; i < c2v->idx[c_id+1]; i++, ii++)
      vtag[c2v->ids[i]] = ii;

    for (cs_lnum_t jf = c2f->idx[c_id]; jf < c2f->idx[c_id+1]; jf++) {

      const cs_lnum_t  f_id = topo->c2f->col_id[jf];
      const cs_quant_t  pfq = quant->face[f_id];

      for (cs_lnum_t je = f2e->idx[f_id]; je < f2e->idx[f_id+1]; je++) {

        const cs_lnum_t  e_id = f2e->col_id[je];
        const cs_lnum_t  eshft = 2*e_id;
        const cs_lnum_t  v1_id = topo->e2v->col_id[eshft];
        const cs_lnum_t  v2_id = topo->e2v->col_id[eshft+1];

        const double pvol = 0.5 * cs_math_voltet(m->vtx_coord + 3*v1_id,
                                                 m->vtx_coord + 3*v2_id,
                                                 pfq.center,
                                                 xc);

        v_vol[vtag[v1_id]] += pvol;
        v_vol[vtag[v2_id]] += pvol;

      } // Loop on face edges

    } // Loop on cell faces

  }  // Loop on cells

  /* Free buffer */
  BFT_FREE(vtag);
}

/*----------------------------------------------------------------------------
 * Compute dual face normals (face crossed by primal edges).
 * Given a cell and an edge, there are two faces attached to the
 * couple (cell, edge)
 * The triplet (edge, face, cell) induces an elementary triangle s(e,f,c)
 * The dual face is the union of these two triangles.
 * Storage based on c2e connectivity
 * ---------------------------------------------------------------------------*/

static void
_compute_dface_quantities(const cs_cdo_connect_t  *topo,
                          cs_cdo_quantities_t     *iq)  /* In/out */
{
  cs_lnum_t  c_id, i, j, k, size, shift, parent;
  cs_nvec3_t  nvec;
  cs_real_3_t  trinorm, xexf, xexc;

  cs_lnum_t  *tag_shift = NULL;

  /* Sanity check */
  assert(topo->e2f != NULL);
  assert(topo->f2c != NULL);
  assert(topo->c2e != NULL);

  const cs_connect_index_t  *c2e = topo->c2e;

  /* Allocate and initialize arrays */
  size = c2e->idx[iq->n_cells];
  BFT_MALLOC(iq->dface, size, cs_dface_t);

  BFT_MALLOC(tag_shift, iq->n_edges, cs_lnum_t);
  for (i = 0; i < iq->n_edges; i++)
    tag_shift[i] = 0;

  for (c_id = 0; c_id < iq->n_cells; c_id++) {

    /* Tag cell edges */
    for (i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; i++)
      tag_shift[c2e->ids[i]] = i+1;

    /* Get cell center */
    const cs_real_t  *xc = iq->cell_centers + 3*c_id;

    for (i = topo->c2f->idx[c_id]; i < topo->c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = topo->c2f->col_id[i];
      const cs_quant_t  f_q = iq->face[f_id]; /* Face quantities */

      for (j = topo->f2e->idx[f_id]; j < topo->f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  e_id = topo->f2e->col_id[j];
        const cs_quant_t  e_q = iq->edge[e_id]; /* Edge quantities */

        /* Compute the vectorial area for the triangle : xc, xf, xe */
        for (k = 0; k < 3; k++) {
          xexf[k] = f_q.center[k] - e_q.center[k];
          xexc[k] = xc[k] - e_q.center[k];
        }
        cs_math_3_cross_product(xexf, xexc, trinorm);
        cs_nvec3(trinorm, &nvec);

        /* One should have (trinorm, te) > 0 */
        const double  orient = _dp3(nvec.unitv, e_q.unitv);
        assert(fabs(orient) > 0);

        if (tag_shift[e_id] > 0) /* First time */
          shift = tag_shift[e_id]-1, tag_shift[e_id] *= -1, parent = 0;
        else /* Second time (<0) */
          tag_shift[e_id] *= -1, shift = tag_shift[e_id]-1, parent = 1;

        /* Store the computed data */
        iq->dface[shift].parent_id[parent] = f_id;
        iq->dface[shift].sface[parent].meas = 0.5*nvec.meas;
        if (orient < 0)
          for (k = 0; k < 3; k++)
            iq->dface[shift].sface[parent].unitv[k] = -nvec.unitv[k];
        else
          for (k = 0; k < 3; k++)
            iq->dface[shift].sface[parent].unitv[k] =  nvec.unitv[k];

      } /* Loop on face edges */

    } /* Loop on cell faces */

  } /* Loop on cells */

  BFT_FREE(tag_shift);

  /* Compute the dual face normal from the two elementary contributions */
  for (c_id = 0; c_id < iq->n_cells; c_id++) {
    for (i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; i++) {

      cs_nvec3_t  t1 = iq->dface[i].sface[0];
      cs_nvec3_t  t2 = iq->dface[i].sface[1];

      for (k = 0; k < 3; k++)
        iq->dface[i].vect[k] = t1.meas*t1.unitv[k] + t2.meas*t2.unitv[k];

    }
  } /* End of loop on cells */

}

/*----------------------------------------------------------------------------
 * Define the cs_quant_info_t structures related to cells, faces and edges
 * ---------------------------------------------------------------------------*/

static void
_compute_quant_info(cs_cdo_quantities_t     *quant)    /* In/out */
{
  assert(quant != NULL); // Sanity check

  /* Cell info (set default values) */
  quant->cell_info.min_id = quant->cell_info.max_id = -1;
  quant->cell_info.h_min = quant->cell_info.meas_min = DBL_MAX;
  quant->cell_info.h_max = quant->cell_info.meas_max = -DBL_MAX;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    const double  meas = quant->cell_vol[c_id];

    if (meas > quant->cell_info.meas_max) {
      quant->cell_info.meas_max = meas;
      quant->cell_info.h_max = pow(meas, cs_math_onethird);
      quant->cell_info.max_id = c_id;
    }
    if (meas < quant->cell_info.meas_min) {
      quant->cell_info.meas_min = meas;
      quant->cell_info.h_min = pow(meas, cs_math_onethird);
      quant->cell_info.min_id = c_id;
    }

  } // Loop on cells

  /* Face info (set default values) */
  quant->face_info.min_id = quant->face_info.max_id = -1;
  quant->face_info.h_min = quant->face_info.meas_min = DBL_MAX;
  quant->face_info.h_max = quant->face_info.meas_max = -DBL_MAX;

  for (cs_lnum_t  f_id = 0; f_id < quant->n_faces; f_id++) {

    const double  meas = quant->face[f_id].meas;

    if (meas > quant->face_info.meas_max) {
      quant->face_info.meas_max = meas;
      quant->face_info.h_max = sqrt(meas);
      quant->face_info.max_id = f_id;
    }
    if (meas < quant->face_info.meas_min) {
      quant->face_info.meas_min = meas;
      quant->face_info.h_min = sqrt(meas);
      quant->face_info.min_id = f_id;
    }

  } // Loop on faces

  /* Edge info (set default values) */
  quant->edge_info.min_id = quant->edge_info.max_id = -1;
  quant->edge_info.h_min = quant->edge_info.meas_min = DBL_MAX;
  quant->edge_info.h_max = quant->edge_info.meas_max = -DBL_MAX;

  for (cs_lnum_t  e_id = 0; e_id < quant->n_edges; e_id++) {

    const double  meas = quant->edge[e_id].meas;

    if (meas > quant->edge_info.meas_max) {
      quant->edge_info.meas_max = meas;
      quant->edge_info.h_max = meas;
      quant->edge_info.max_id = e_id;
    }
    if (meas < quant->edge_info.meas_min) {
      quant->edge_info.meas_min = meas;
      quant->edge_info.h_min = meas;
      quant->edge_info.min_id = e_id;
    }

  } // Loop on edges

}
/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : copy data from Saturne structure
  ----------------------------------------------------------------------------*/

static void
_saturn_algorithm(const cs_mesh_t             *mesh,
                  const cs_mesh_quantities_t  *mq,
                  cs_cdo_quantities_t         *cdoq) /* In/out */
{
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_faces = n_i_faces + n_b_faces;

  assert(mq != NULL);

  BFT_MALLOC(cdoq->face, n_faces, cs_quant_t);
  BFT_MALLOC(cdoq->cell_centers, 3*n_cells, cs_real_t);
  BFT_MALLOC(cdoq->cell_vol, n_cells, cs_real_t);

  /* Copy cell volumes and compute vol_tot */
  memcpy(cdoq->cell_centers, mq->cell_cen, 3*n_cells*sizeof(cs_real_t));
  memcpy(cdoq->cell_vol, mq->cell_vol, n_cells*sizeof(cs_real_t));

  cdoq->vol_tot = 0.0;
  for (cs_lnum_t i = 0; i < n_cells; i++)
    cdoq->vol_tot += mq->cell_vol[i];

  /* Concatenate face centers for interior and border faces */
# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    for (int k = 0; k < 3; k++)
      cdoq->face[f_id].center[k] = mq->i_face_cog[3*f_id+k];

  for (cs_lnum_t j = 0, f_id = n_i_faces; j < n_b_faces; j++, f_id++) {
    for (int k = 0; k < 3; k++)
      cdoq->face[f_id].center[k] = mq->b_face_cog[3*j+k];
  }
}

/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : cell centers are computed as the
  vertex average over cell vertices. Other quantities are copied from those
  computed by Code_Saturne (default algorithm)
  ----------------------------------------------------------------------------*/

static void
_vtx_algorithm(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mq,
               const cs_cdo_connect_t      *connect,
               cs_cdo_quantities_t         *quant) /* In/out */
{
  const cs_lnum_t  n_cells = mesh->n_cells;
  const cs_lnum_t  n_i_faces = mesh->n_i_faces;
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_lnum_t  n_faces = n_i_faces + n_b_faces;
  const cs_connect_index_t  *c2v = connect->c2v;

  assert(mq != NULL);

  BFT_MALLOC(quant->face, n_faces, cs_quant_t);
  BFT_MALLOC(quant->cell_centers, 3*n_cells, cs_real_t);
  BFT_MALLOC(quant->cell_vol, n_cells, cs_real_t);

  /* Copy cell volumes and compute vol_tot */
  memcpy(quant->cell_vol, mq->cell_vol, n_cells*sizeof(cs_real_t));
  quant->vol_tot = 0.0;
  for (cs_lnum_t i = 0; i < n_cells; i++)
    quant->vol_tot += mq->cell_vol[i];

  /* Concatenate face centers for interior and border faces */
# pragma omp parallel for if (n_i_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++)
    for (int k = 0; k < 3; k++)
      quant->face[f_id].center[k] = mq->i_face_cog[3*f_id+k];

  for (cs_lnum_t j = 0, f_id = n_i_faces; j < n_b_faces; j++, f_id++)
    for (int k = 0; k < 3; k++)
      quant->face[f_id].center[k] = mq->b_face_cog[3*j+k];

  /* Compute cell centers */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_lnum_t  vs = c2v->idx[c_id];
    const cs_lnum_t  ve = c2v->idx[c_id+1];

    /* Sanity checks */
    assert(ve - vs > 0);
    assert(ve - vs < SHRT_MAX);

    const double  coef = 1./(ve-vs);
    double  *xc = quant->cell_centers + 3*c_id;

    xc[0] = xc[1] = xc[2] = 0;
    for (cs_lnum_t jv = vs; jv < ve; jv++) {

      const cs_real_t  *xv = mesh->vtx_coord + 3*c2v->ids[jv];

      xc[0] += coef * xv[0];
      xc[1] += coef * xv[1];
      xc[2] += coef * xv[2];

    } // Loop on cell vertices

  } // Loop on cells

}

/*----------------------------------------------------------------------------
 * Algorithm for computing cell barycenters inspired from the article
 * "Fast and accurate computation of polyhedral mass properties"
 * Journal of Graphics, 1997 by Brian Mirtich
 *
 * Compute also : face centers, cell volumes.
 * ---------------------------------------------------------------------------*/

static void
_mirtich_algorithm(const cs_mesh_t             *mesh,
                   const cs_mesh_quantities_t  *mq,
                   const cs_cdo_connect_t      *connect,
                   cs_cdo_quantities_t         *quant) /* In/out */
{
  cs_lnum_t  i, k, c_id, f_id, A, B, C, sgn;
  double  Fvol, inv_surf;
  _cdo_fspec_t  fspec;
  _cdo_fsubq_t  fsubq;

  const int X = 0, Y = 1, Z = 2;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_faces = quant->n_faces;

  /* Sanity check */

  assert(connect->f2c != NULL);
  assert(connect->c2f != NULL);

  /* Allocate and initialize cell quantities */

  BFT_MALLOC(quant->face, n_faces, cs_quant_t);
  BFT_MALLOC(quant->cell_centers, 3*n_cells, cs_real_t);
  BFT_MALLOC(quant->cell_vol, n_cells, cs_real_t);

  for (i = 0; i < n_cells; i++) {
    quant->cell_vol[i] = 0.0;
    for (k = 0; k < 3; k++)
      quant->cell_centers[3*i+k] = 0.0;
  }

  for (f_id = 0; f_id < n_faces; f_id++) {   /* Loop on faces */

    /* Choose gamma to maximize normal according gamma (x, y, or z)
       Define a direct basis (alpha, beta, gamma) with this choice
       Compute omega = - <n, P> where P belongs to the face */

    fspec = _get_fspec(f_id, mesh, mq);

    A = fspec.XYZ[X];
    B = fspec.XYZ[Y];
    C = fspec.XYZ[Z];

    fsubq = _get_fsub_quantities(f_id, connect, mesh->vtx_coord, fspec);

    inv_surf = 1.0/fsubq.F1;
    quant->face[f_id].center[A] = inv_surf * fsubq.Fa;
    quant->face[f_id].center[B] = inv_surf * fsubq.Fb;
    quant->face[f_id].center[C] = inv_surf * fsubq.Fc;

    /* Update cell quantities */
    for (i = connect->f2c->idx[f_id]; i < connect->f2c->idx[f_id+1]; i++) {

      c_id = connect->f2c->col_id[i];
      sgn = connect->f2c->sgn[i];

      Fvol = ( (fspec.XYZ[X] == 0) ? fsubq.Fa :
               ( (fspec.XYZ[Y] == 0) ? fsubq.Fb : fsubq.Fc) );

      quant->cell_vol[c_id] += sgn * fspec.q.unitv[0] * Fvol;
      quant->cell_centers[3*c_id + A] += sgn * fspec.q.unitv[A] * fsubq.Fa2;
      quant->cell_centers[3*c_id + B] += sgn * fspec.q.unitv[B] * fsubq.Fb2;
      quant->cell_centers[3*c_id + C] += sgn * fspec.q.unitv[C] * fsubq.Fc2;

    } /* End of loop on cell faces */

  } /* End of loop on faces */

  /* Compute cell center of gravity and total volume */
  quant->vol_tot = 0.0;

  for (i = 0; i < n_cells; i++) {

    double  inv_vol2 = 0.5 / quant->cell_vol[i];

    quant->vol_tot += quant->cell_vol[i];

    for (k = 0; k < 3; k++)
      quant->cell_centers[3*i+k] *= inv_vol2;

#if CDO_QUANTITIES_DBG > 1
    printf("\n (%4d) volINNOV: % -12.5e | volSAT % -12.5e | %-12.5e\n",
           i+1, quant->cell_vol[i], mq->cell_vol[i],
           fabs(quant->cell_vol[i] - mq->cell_vol[i]));
    printf(" INNOVccog (% -.4e, % -.4e, % -.4e)",
           quant->cell_centers[3*i], quant->cell_centers[3*i+1],
           quant->cell_centers[3*i+2]);
    printf(" SATccog (% -.4e, % -.4e, % -.4e) Delta (% -.4e, % -.4e, % -.4e)\n",
           mq->cell_cen[3*i],mq->cell_cen[3*i+1],mq->cell_cen[3*i+2],
           fabs(quant->cell_centers[3*i]-mq->cell_cen[3*i]),
           fabs(quant->cell_centers[3*i+1]-mq->cell_cen[3*i+1]),
           fabs(quant->cell_centers[3*i+2]-mq->cell_cen[3*i+2]));
#endif

  }

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a cs_cdo_quantities_t structure
 *
 * \param[in]  m           pointer to a cs_mesh_t structure
 * \param[in]  mq          pointer to a cs_mesh_quantities_t structure
 * \param[in]  topo        pointer to a cs_cdo_connect_t structure
 *
 * \return  a new allocated pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_quantities_t *
cs_cdo_quantities_build(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *mq,
                        const cs_cdo_connect_t      *topo)
{
  cs_cdo_quantities_t  *cdoq = NULL;
  cs_cdo_cc_algo_t  cc_algo = CS_CDO_CC_SATURNE; // default value

  /* Sanity check */
  assert(topo != NULL);
  assert(topo->c2f != NULL);

  /* Build cs_cdo_quantities_t structure */
  BFT_MALLOC(cdoq, 1, cs_cdo_quantities_t);

  cdoq->n_cells = m->n_cells;
  cdoq->n_i_faces = m->n_i_faces;
  cdoq->n_b_faces = m->n_b_faces;
  cdoq->n_faces = m->n_i_faces + m->n_b_faces;
  cdoq->n_vertices = m->n_vertices;
  cdoq->vtx_coord = m->vtx_coord;

  if (topo->e2v != NULL)
    cdoq->n_edges = topo->e2v->n_rows;
  else /* Not used by the numerical scheme */
    cdoq->n_edges = -1;

  /* Initialize quantities not used by all schemes */
  cdoq->edge = NULL;
  cdoq->dcell_vol = NULL;
  cdoq->dface = NULL;
  cdoq->dedge = NULL;

  /* User can modify the default value */
  cc_algo = cs_user_cdo_geometric_settings();

  switch (cc_algo) { /* Compute face/cell centers and cell volumes */

  case CS_CDO_CC_BARYC:
    bft_printf(" -cdo- Cell.Center.Algo >> Mirtich\n");
    /* Compute (real) barycentric centers, face centers and cell volumes */
    _mirtich_algorithm(m, mq, topo, cdoq);
    break;

  case CS_CDO_CC_MEANV:
    bft_printf(" -cdo- Cell.Center.Algo >> Vertices.MeanValue\n");
    /* Compute cell centers, face centers and cell volumes and total volume */
    _vtx_algorithm(m, mq, topo, cdoq);
    break;

  case CS_CDO_CC_SATURNE:
    bft_printf(" -cdo- Cell.Center.Algo >> Original\n");
    _saturn_algorithm(m, mq, cdoq);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Unkwown algorithm for cell center computation\n"));

  } /* switch cc_algo */

  /* Finalize definition of cs_quant_t struct. for faces.
     Define face normal with unitary norm and face area */
# pragma omp parallel for if (m->n_i_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < m->n_i_faces; f_id++) {
    cs_nvec3_t  v;
    cs_nvec3((mq->i_face_normal + 3*f_id), &v);
    cdoq->face[f_id].meas = v.meas;
    for (int k = 0; k < 3; k++)
      cdoq->face[f_id].unitv[k] = v.unitv[k];
  }

  for (cs_lnum_t i = 0, f_id = m->n_i_faces; i < m->n_b_faces; i++, f_id++) {
    cs_nvec3_t  v;
    cs_nvec3((mq->b_face_normal + 3*i), &v);
    cdoq->face[f_id].meas = v.meas;
    for (int k = 0; k < 3; k++)
      cdoq->face[f_id].unitv[k] = v.unitv[k];
  }

  /* Compute dual edge quantities */
  const cs_lnum_t  idx_size = topo->c2f->idx[cdoq->n_cells];
  const cs_sla_matrix_t  *c2f = topo->c2f;

  BFT_MALLOC(cdoq->dedge, idx_size, cs_nvec3_t);

# pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    cs_real_t  *xc = cdoq->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      cs_lnum_t  f_id = c2f->col_id[i];
      short int  sgn = c2f->sgn[i];

      cs_math_3_length_unitv(xc, cdoq->face[f_id].center,
                             &(cdoq->dedge[i].meas),
                             cdoq->dedge[i].unitv);

      for (int k = 0; k < 3; k++)
        cdoq->dedge[i].unitv[k] *= sgn;

    } // Loop on cell faces

  } /* End of loop on cells */

  /* Compute edge quantities if needed */
  if (cdoq->n_edges > 0) {
    _compute_edge_quantities(m, topo, cdoq);
    _compute_dface_quantities(topo, cdoq);
  }

  _compute_dcell_quantities(topo, cdoq);

  /* Define cs_quant_info_t structure */
  _compute_quant_info(cdoq);

  return cdoq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdo_quantities_t structure
 *
 * \param[in]  q        pointer to the cs_cdo_quantities_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_quantities_t *
cs_cdo_quantities_free(cs_cdo_quantities_t   *q)
{
  if (q == NULL)
    return q;

  BFT_FREE(q->cell_centers);
  BFT_FREE(q->cell_vol);

  BFT_FREE(q->face);
  BFT_FREE(q->dedge);

  BFT_FREE(q->edge);
  BFT_FREE(q->dface);

  BFT_FREE(q->dcell_vol);

  /* vtx_coord is free when the structure cs_mesh_t is destroyed */

  BFT_FREE(q);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Summarize generic information about the cdo mesh quantities
 *
 * \param[in]  cdoq     pointer to cs_cdo_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_summary(const cs_cdo_quantities_t  *quant)
{
  /* Output */
  bft_printf("\n CDO mesh quantities information:\n");
  bft_printf(" --cdo-- h_cell  %6.4e %6.4e (min/max)\n",
             quant->cell_info.h_min, quant->cell_info.h_max);
  bft_printf(" --cdo-- h_face  %6.4e %6.4e (min/max)\n",
             quant->face_info.h_min, quant->face_info.h_max);
  bft_printf(" --cdo-- h_edge  %6.4e %6.4e (min/max)\n",
             quant->edge_info.h_min, quant->edge_info.h_max);
  bft_printf("\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_cdo_quantities_t structure
 *
 * \param[in]  cdoq     pointer to cs_cdo_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_dump(const cs_cdo_quantities_t  *cdoq)
{
  cs_lnum_t  i, p;

  FILE  *fdump = NULL;

  fdump = fopen("cdo_quantities_dump.dat", "w");

  if (cdoq == NULL) {
    fprintf(fdump, "Empty structure.\n");
    fclose(fdump);
    return;
  }

  fprintf(fdump, "\n Quantities structure: %p\n\n", (const void *)cdoq);

  fprintf(fdump, " -cdoq- n_cells =    %d\n", cdoq->n_cells);
  fprintf(fdump, " -cdoq- n_faces =    %d\n", cdoq->n_faces);
  fprintf(fdump, " -cdoq- n_edges =    %d\n", cdoq->n_edges);
  fprintf(fdump, " -cdoq- n_vertices = %d\n", cdoq->n_vertices);
  fprintf(fdump, " -cdoq- Total volume = %.6e\n\n", cdoq->vol_tot);

  fprintf(fdump, "\n *** Cell Quantities ***\n");
  fprintf(fdump, "-msg- num.; volume ; center (3)\n");
  for (i = 0; i < cdoq->n_cells; i++) {
    p = 3*i;
    fprintf(fdump, " [%6d] | %12.8e | % -12.8e | % -12.8e |% -12.8e\n",
            i+1, cdoq->cell_vol[i], cdoq->cell_centers[p],
            cdoq->cell_centers[p+1], cdoq->cell_centers[p+2]);
  }

  fprintf(fdump, "\n\n *** Face Quantities ***\n");
  fprintf(fdump, "-msg- num. ; measure ; unitary vector (3) ; center (3)\n");
  for (i = 0; i < cdoq->n_faces; i++)
    cs_quant_dump(fdump, i+1, cdoq->face[i]);

  fprintf(fdump, "\n\n *** Edge Quantities ***\n");
  fprintf(fdump, "-msg- num. ; measure ; unitary vector (3) ; center (3)\n");
  for (i = 0; i < cdoq->n_edges; i++)
    cs_quant_dump(fdump, i+1, cdoq->edge[i]);

  fclose(fdump);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_quant_t structure
 *
 * \param[in]  f         FILE struct (stdout if NULL)
 * \param[in]  num       entity number related to this quantity struct.
 * \param[in]  q         cs_quant_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_quant_dump(FILE             *f,
              cs_lnum_t         num,
              const cs_quant_t  q)
{
  FILE  *_f = f;

  if (_f == NULL) _f = stdout;

  fprintf(_f, " -cdoq-  [%8d] | % -10.6e | % -10.6e | % -10.6e | % -10.6e |"
          " % -10.6e | % -10.6e | % -10.6e\n", num, q.meas,
          q.unitv[0], q.unitv[1], q.unitv[2], q.center[0], q.center[1],
          q.center[2]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex the dual cell volume which is also
 *
 *                sum    |celld(v) cap c| = pvol_v
 *              c in C_v
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantites_t structure
 * \param[in, out] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_vtx(const cs_cdo_connect_t     *connect,
                    const cs_cdo_quantities_t  *quant,
                    double                     *p_pvol[])
{
  cs_lnum_t  i, j;

  double  *pvol = *p_pvol;

  const cs_connect_index_t  *c2v = connect->c2v;

  /* Allocate if needed and initialize */
  if (pvol == NULL)
    BFT_MALLOC(pvol, quant->n_vertices, double);
  for (i = 0; i < quant->n_vertices; i++)
    pvol[i] = 0;

  for (i = 0; i < quant->n_cells; i++)
    for (j = c2v->idx[i]; j < c2v->idx[i+1]; j++)
      pvol[c2v->ids[j]] += quant->dcell_vol[j];

  /* Return pointer */

  *p_pvol = pvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each edge a related volume pvol_e which constitutes
 *         a partition of unity
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantites_t structure
 * \param[in, out] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_edge(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     double                      *p_pvol[])
{
  cs_lnum_t  i, j;
  double  dvol;

  double  *pvol = *p_pvol;

  /* Allocate if needed and initialize */
  if (pvol == NULL)
    BFT_MALLOC(pvol, quant->n_edges, double);
  for (i = 0; i < quant->n_edges; i++)
    pvol[i] = 0;

  for (i = 0; i < quant->n_cells; i++) {
    for (j = connect->c2e->idx[i]; j < connect->c2e->idx[i+1]; j++) {

      const cs_lnum_t  e_id = connect->c2e->ids[j];
      const cs_quant_t  peq = quant->edge[e_id];
      const cs_nvec3_t  df0q = quant->dface[j].sface[0];
      const cs_nvec3_t  df1q = quant->dface[j].sface[1];

      dvol  = df0q.meas * _dp3(peq.unitv, df0q.unitv);
      dvol += df1q.meas * _dp3(peq.unitv, df1q.unitv);
      pvol[e_id] += dvol * cs_math_onethird * peq.meas;

    }
  }

  /* Return pointer */
  *p_pvol = pvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each face a related volume pvol_f which constitutes
 *         a partition of unity
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantites_t structure
 * \param[in, out] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_face(const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     double                     *p_pvol[])
{
  double  *pvol = *p_pvol;

  const cs_sla_matrix_t  *c2f = connect->c2f;

  /* Allocate if needed */
  if (pvol == NULL)
    BFT_MALLOC(pvol, quant->n_faces, double);
  /* Initialize */
  for (cs_lnum_t i = 0; i < quant->n_faces; i++)
    pvol[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->col_id[j];
      const cs_quant_t  pfq = quant->face[f_id];
      const cs_nvec3_t  deq = quant->dedge[j];

      /* Compute volume of the pyramid p_fc */
      pvol[f_id] +=
        cs_math_onethird * pfq.meas * deq.meas * _dp3(pfq.unitv, deq.unitv);

    }
  }

  /* Return pointer */
  *p_pvol = pvol;
}

/*----------------------------------------------------------------------------*/

#undef _dp3
#undef _n3

END_C_DECLS
