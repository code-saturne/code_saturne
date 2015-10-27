/*============================================================================
 * Manage geometrical quantities needed in CDO schemes
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

#include <string.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_cdo_toolbox.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

#define  CDO_QUANTITIES_DBG 0  /* Switch off/on debug information */

/* Temporary structures to build mesh quantities */

typedef struct {

  int         XYZ[3]; /* Direct permutation of the ref. axis such that nZ
                         is maximal */
  cs_qvect_t  q;      /* face surface and its unit normal */
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

static const double one_3 = 1/3.;
static const double one_6 = 1/6.;
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

    cs_qvect(&(mq->i_face_normal[3*f+k]), &(fspec.q));

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

    cs_qvect(&(mq->b_face_normal[3*f+k]), &(fspec.q));

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
  projq.pa  *=  one_6;
  projq.pb  *= -one_6;
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
                          cs_cdo_quantities_t     *iq)    /* In/out */
{
  int  e_id, f_id, c_id, v_id, i, j, k, l, shift;
  double  vol, xc[3], xv[3];
  cs_quant_t  eq, fq;

  cs_lnum_t  *cell_shift = NULL;

  const cs_mesh_t  *m = cs_glob_mesh;

  /* Compute part of dual volume related to each primal cell */
  BFT_MALLOC(iq->dcell_vol, topo->c2v->idx[iq->n_cells], double);
  for (i = 0; i < topo->c2v->idx[iq->n_cells]; i++)
    iq->dcell_vol[i] = 0.0;

  BFT_MALLOC(cell_shift, iq->n_cells, cs_lnum_t);
  for (i = 0; i < iq->n_cells; i++)
    cell_shift[i] = topo->c2v->idx[i];

  /* Fill array : one needs c2v->lst to be ordered. It should have been done. */
  for (v_id = 0; v_id < iq->n_vertices; v_id++) {

    for (k = 0; k < 3; k++)
      xv[k] = m->vtx_coord[3*v_id+k];

    for (i = topo->v2e->idx[v_id]; i < topo->v2e->idx[v_id+1]; i++) {

      e_id = topo->v2e->col_id[i];
      eq = iq->edge[e_id]; /* Get quantities related to this edge */

      for (j = topo->e2f->idx[e_id]; j < topo->e2f->idx[e_id+1]; j++) {

        f_id = topo->e2f->col_id[j];
        fq = iq->face[f_id]; /* Get quantities related to this face */

        /* Get related cell(s) */
        for (l = topo->f2c->idx[f_id]; l < topo->f2c->idx[f_id+1]; l++) {

          /* Get cell center */
          c_id = topo->f2c->col_id[l];
          for (k = 0; k < 3; k++)
            xc[k] = iq->cell_centers[3*c_id+k];

          vol = cs_voltet(xv, eq.center, fq.center, xc);

          shift = cell_shift[c_id];
          if (topo->c2v->ids[shift] != v_id)
            cell_shift[c_id] += 1, shift++;
          iq->dcell_vol[shift] += vol;

        } /* End of loop on cells related to this face */

      } /* End of loop on faces related to this edge */

    } /* End of loop on edge related to this vertex */

  } /* End of loop on vertices */

  /* Free buffer */
  BFT_FREE(cell_shift);
}

/*----------------------------------------------------------------------------
 * Compute dual face normals (face crossed by primal edges).
 * Given a cell, a face and an edge, there are two dual face normals
 * Storage based on c2e connectivity
 * ---------------------------------------------------------------------------*/

static void
_compute_dface_quantities(const cs_cdo_connect_t  *topo,
                          cs_cdo_quantities_t     *iq)  /* In/out */
{
  cs_lnum_t  c_id, i, j, k, size, shift, parent;
  double orient, area, inv;
  cs_real_3_t  trinorm, xexf, xexc, xc;

  cs_lnum_t  *tag_shift = NULL;

  /* Sanity check */
  assert(topo->e2f != NULL);
  assert(topo->f2c != NULL);
  assert(topo->c2e != NULL);

  /* Allocate and initialize arrays */
  size = topo->c2e->idx[iq->n_cells];
  BFT_MALLOC(iq->dface, size, cs_dface_t);

  BFT_MALLOC(tag_shift, iq->n_edges, cs_lnum_t);
  for (i = 0; i < iq->n_edges; i++)
    tag_shift[i] = 0;

  for (c_id = 0; c_id < iq->n_cells; c_id++) {

    /* Tag cell edges */
    for (i = topo->c2e->idx[c_id]; i < topo->c2e->idx[c_id+1]; i++)
      tag_shift[topo->c2e->ids[i]] = i+1;

    /* Get cell center */
    for (k = 0; k < 3; k++)
      xc[k] = iq->cell_centers[3*c_id+k];

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
        _cp3(xexf, xexc, &trinorm);

        /* One should have (trinorm, te) > 0 */
        orient = _dp3(trinorm, e_q.unitv);
        assert(fabs(orient) > 0);

        if (tag_shift[e_id] > 0) /* First time */
          shift = tag_shift[e_id]-1, tag_shift[e_id] *= -1, parent = 0;
        else /* Second time (<0) */
          tag_shift[e_id] *= -1, shift = tag_shift[e_id]-1, parent = 1;

        iq->dface[shift].parent_id[parent] = f_id;
        if (orient < 0) {
          for (k = 0; k < 3; k++)
            iq->dface[shift].unitv[3*parent+k] = -0.5 * trinorm[k];
        }
        else {
          for (k = 0; k < 3; k++)
            iq->dface[shift].unitv[3*parent+k] =  0.5 * trinorm[k];
        }

      } /* Loop on face edges */

    } /* Loop on cell faces */

  } /* Loop on cells */

  /* Normalize dual face normal and fill dual face area  */
  for (c_id = 0; c_id < iq->n_cells; c_id++) {
    for (i = topo->c2e->idx[c_id]; i < topo->c2e->idx[c_id+1]; i++) {

      for (k = 0; k < 3; k++)
        iq->dface[i].vect[k] = iq->dface[i].unitv[k] + iq->dface[i].unitv[3+k];

      for (parent = 0; parent < 2; parent++) {
        area = _n3(&(iq->dface[i].unitv[3*parent]));
        iq->dface[i].meas[parent] = area;
        inv = 1 / area;
        for (k = 0; k < 3; k++)
          iq->dface[i].unitv[3*parent+k] *= inv;
      }

    }
  } /* End of loop on cells */

  BFT_FREE(tag_shift);
}

/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : copy data from Saturne structure
  ----------------------------------------------------------------------------*/

static void
_saturn_algorithm(const cs_mesh_t             *mesh,
                  const cs_mesh_quantities_t  *mq,
                  cs_cdo_quantities_t         *cdoq) /* In/out */
{
  cs_lnum_t  i, j, k;

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
  for (i = 0; i < n_cells; i++)
    cdoq->vol_tot += mq->cell_vol[i];

  /* Concatenate face centers for interior and border faces */
  for (i = 0; i < n_i_faces; i++)
    for (k = 0; k < 3; k++)
      cdoq->face[i].center[k] = mq->i_face_cog[3*i+k];
  for (j = 0; j < n_b_faces; j++, i++)
    for (k = 0; k < 3; k++)
      cdoq->face[i].center[k] = mq->b_face_cog[3*j+k];

}

/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : cell centers are computed as the
  vertex average over cell vertices. Same thing for face centers.
  ----------------------------------------------------------------------------*/

static void
_vtx_algorithm(const cs_mesh_t             *mesh,
               const cs_mesh_quantities_t  *mq,
               const cs_cdo_connect_t      *connect,
               cs_cdo_quantities_t         *iq) /* In/out */
{
  int  i, j, k, c_id, f_id, v_id, s, e, n_face_vertices, count;
  double  xf[3], xc[3], inv;

  int  *vtag = NULL;

  const double  *xv = mesh->vtx_coord;
  const int  n_cells = mesh->n_cells;
  const int  n_vertices = mesh->n_vertices;
  const int  n_i_faces = mesh->n_i_faces;
  const int  n_b_faces = mesh->n_b_faces;
  const int  n_faces = n_i_faces + n_b_faces;

  /* Sanity check */
  assert(connect->f2c != NULL);
  assert(connect->c2f != NULL);

  /* Build cell centers and face centers */
  BFT_MALLOC(iq->face, n_faces, cs_quant_t);
  BFT_MALLOC(iq->cell_centers, 3*n_cells, double);
  BFT_MALLOC(vtag, n_vertices, int);
  for (i = 0; i < n_vertices; i++)
    vtag[i] = -1;

  for (c_id = 0; c_id < n_cells; c_id++) {   /* Loop on cells */

    count = 0; /* number of cell vertices */

    for (k = 0; k < 3; k++)
      xc[k] = 0.0;

    for (i = connect->c2f->idx[c_id]; i < connect->c2f->idx[c_id+1]; i++) {

      f_id = connect->c2f->col_id[i];
      for (k = 0; k < 3; k++)
        xf[k] = 0.0;

      if (f_id < n_i_faces) { /* Interior faces */

        s = mesh->i_face_vtx_idx[f_id], e = mesh->i_face_vtx_idx[f_id+1];
        n_face_vertices = e - s;

        for (j = s; j < e; j++) {

          v_id = mesh->i_face_vtx_lst[j];
          for (k = 0; k < 3; k++)
            xf[k] += xv[3*v_id+k];

          if (vtag[v_id] != c_id) { /* Not already computed */
            count++;
            vtag[v_id] = c_id;
            for (k = 0; k < 3; k++)
              xc[k] += xv[3*v_id+k];
          }

        } /* End of loop on face vertices */

      }
      else { /* Border faces */

        int _f_id = f_id - n_i_faces;
        s = mesh->b_face_vtx_idx[_f_id], e = mesh->b_face_vtx_idx[_f_id+1];
        n_face_vertices = e - s;

        for (j = s; j < e; j++) {

          v_id = mesh->b_face_vtx_lst[j];
          for (k = 0; k < 3; k++)
            xf[k] += xv[3*v_id+k];

          if (vtag[v_id] != c_id) { /* Not already computed */
            count++;
            vtag[v_id] = c_id;
            for (k = 0; k < 3; k++)
              xc[k] += xv[3*v_id+k];
          }

        } /* End of loop on face vertices */

      } /* End of border faces */

      inv = 1/n_face_vertices;
      for (k = 0; k < 3; k++)
        iq->face[f_id].center[k] = inv*xf[k];

    } /* End of loop on cell faces */

    inv = 1/count;
    for (k = 0; k < 3; k++)
      iq->cell_centers[3*c_id+k] = inv*xc[k];

  } /* End of loop on cells */

  /* Free memory */
  BFT_FREE(vtag);

  /* Copy cell volumes and compute vol_tot */
  BFT_MALLOC(iq->cell_vol, n_cells, double);
  iq->vol_tot = 0.0;
  for (i = 0; i < n_cells; i++) {
    iq->cell_vol[i] = mq->cell_vol[i];
    iq->vol_tot += mq->cell_vol[i];
  }

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
                   cs_cdo_quantities_t         *iq) /* In/out */
{
  cs_lnum_t  i, k, c_id, f_id, A, B, C, sgn;
  double  Fvol, inv_surf;
  _cdo_fspec_t  fspec;
  _cdo_fsubq_t  fsubq;

  const int X = 0, Y = 1, Z = 2;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_faces = iq->n_faces;

  /* Sanity check */

  assert(connect->f2c != NULL);
  assert(connect->c2f != NULL);

  /* Allocate and initialize cell quantities */

  BFT_MALLOC(iq->face, n_faces, cs_quant_t);
  BFT_MALLOC(iq->cell_centers, 3*n_cells, cs_real_t);
  BFT_MALLOC(iq->cell_vol, n_cells, cs_real_t);

  for (i = 0; i < n_cells; i++) {
    iq->cell_vol[i] = 0.0;
    for (k = 0; k < 3; k++)
      iq->cell_centers[3*i+k] = 0.0;
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
    iq->face[f_id].center[A] = inv_surf * fsubq.Fa;
    iq->face[f_id].center[B] = inv_surf * fsubq.Fb;
    iq->face[f_id].center[C] = inv_surf * fsubq.Fc;

    /* Update cell quantities */
    for (i = connect->f2c->idx[f_id]; i < connect->f2c->idx[f_id+1]; i++) {

      c_id = connect->f2c->col_id[i];
      sgn = connect->f2c->sgn[i];

      Fvol = ( (fspec.XYZ[X] == 0) ? fsubq.Fa :
               ( (fspec.XYZ[Y] == 0) ? fsubq.Fb : fsubq.Fc) );

      iq->cell_vol[c_id] += sgn * fspec.q.unitv[0] * Fvol;
      iq->cell_centers[3*c_id + A] += sgn * fspec.q.unitv[A] * fsubq.Fa2;
      iq->cell_centers[3*c_id + B] += sgn * fspec.q.unitv[B] * fsubq.Fb2;
      iq->cell_centers[3*c_id + C] += sgn * fspec.q.unitv[C] * fsubq.Fc2;

    } /* End of loop on cell faces */

  } /* End of loop on faces */

  /* Compute cell center of gravity and total volume */
  iq->vol_tot = 0.0;

  for (i = 0; i < n_cells; i++) {

    double  inv_vol2 = 0.5 / iq->cell_vol[i];

    iq->vol_tot += iq->cell_vol[i];

    for (k = 0; k < 3; k++)
      iq->cell_centers[3*i+k] *= inv_vol2;

#if CDO_QUANTITIES_DBG > 1
    printf("\n (%4d) volINNOV: % -12.5e | volSAT % -12.5e | %-12.5e\n",
           i+1, iq->cell_vol[i], mq->cell_vol[i],
           fabs(iq->cell_vol[i] - mq->cell_vol[i]));
    printf(" INNOVccog (% -.4e, % -.4e, % -.4e)",
           iq->cell_centers[3*i], iq->cell_centers[3*i+1],
           iq->cell_centers[3*i+2]);
    printf(" SATccog (% -.4e, % -.4e, % -.4e) Delta (% -.4e, % -.4e, % -.4e)\n",
           mq->cell_cen[3*i],mq->cell_cen[3*i+1],mq->cell_cen[3*i+2],
           fabs(iq->cell_centers[3*i]-mq->cell_cen[3*i]),
           fabs(iq->cell_centers[3*i+1]-mq->cell_cen[3*i+1]),
           fabs(iq->cell_centers[3*i+2]-mq->cell_cen[3*i+2]));
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
  short int  sgn;
  int  i, j, k, ii, nnz, c_id, f_id;
  double  mes, inv, vec[3], xc[3];

  cs_cdo_quantities_t  *cdoq = NULL;
  cs_cdo_cc_algo_t cc_algo = CS_CDO_CC_SATURNE; // default value

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
  for (i = 0; i < m->n_i_faces; i++) {

    for (k = 0; k < 3; k++)
      vec[k] = mq->i_face_normal[3*i+k];
    mes = _n3(vec);
    inv = 1.0/mes;

    cdoq->face[i].meas = mes;
    for (k = 0; k < 3; k++)
      cdoq->face[i].unitv[k] = inv*vec[k];

  }

  for (i = 0, j = m->n_i_faces; i < m->n_b_faces; i++, j++) {

    for (k = 0; k < 3; k++)
      vec[k] = mq->b_face_normal[3*i+k];
    mes = _n3(vec);
    inv = 1.0/mes;

    cdoq->face[j].meas = mes;
    for (k = 0; k < 3; k++)
      cdoq->face[j].unitv[k] = inv*vec[k];

  }

  /* Compute dual edge quantities */
  nnz = topo->c2f->idx[cdoq->n_cells];
  BFT_MALLOC(cdoq->dedge, 4*nnz, double);

  for (c_id = 0; c_id < cdoq->n_cells; c_id++) {

    for (k = 0; k < 3; k++)
      xc[k] = cdoq->cell_centers[3*c_id+k];

    for (i = topo->c2f->idx[c_id]; i < topo->c2f->idx[c_id+1]; i++) {

      f_id = topo->c2f->col_id[i];
      sgn = topo->c2f->sgn[i];
      for (k = 0; k < 3; k++)
        vec[k] = cdoq->face[f_id].center[k] - xc[k];
      mes = _n3(vec);

      /* Fill dedge */
      cdoq->dedge[4*i] = mes;
      ii = 4*i+1, inv = 1.0/mes;
      for (k = 0; k < 3; k++)
        cdoq->dedge[ii+k] = sgn*inv*vec[k];

    }

  } /* End of loop on cells */

  /* Compute edge quantities if needed */
  if (cdoq->n_edges > 0)
    _compute_edge_quantities(m, topo, cdoq);

  _compute_dface_quantities(topo, cdoq);
  _compute_dcell_quantities(topo, cdoq);

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
      const cs_quant_t  edgeq = quant->edge[e_id];
      const cs_dface_t  dface = quant->dface[j];

      dvol  = dface.meas[0] * _dp3(edgeq.unitv, &(dface.unitv[0]));
      dvol += dface.meas[1] * _dp3(edgeq.unitv, &(dface.unitv[3]));
      pvol[e_id] += dvol * one_3 * edgeq.meas;

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
  cs_lnum_t  i, k, pos;
  double  height;
  cs_real_3_t  xfc, xc;

  double  *pvol = *p_pvol;

  const cs_sla_matrix_t  *c2f = connect->c2f;

  /* Allocate if needed and initialize */
  if (pvol == NULL)
    BFT_MALLOC(pvol, quant->n_faces, double);
  for (i = 0; i < quant->n_faces; i++)
    pvol[i] = 0;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    for (k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    for (pos = c2f->idx[c_id]; pos < c2f->idx[c_id+1]; pos++) {

      cs_lnum_t  f_id = c2f->col_id[pos];
      const cs_quant_t  fq = quant->face[f_id];

      /* Step 1: Compute h_(f,c) the height of the pyramid of base face f */
      for (k = 0; k < 3; k++)
        xfc[k] = xc[k] - fq.center[k];
      height = fabs(_dp3(fq.unitv, xfc));

      /* Step 2: Compute volume of the pyramid */
      pvol[f_id] += one_3 * fq.meas * height;

    }
  }

  /* Return pointer */
  *p_pvol = pvol;
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
 * \param[in]      loc_ids   indirection to a local numbering
 * \param[in, out] wf        already allocated to n_max_vbyc (reset)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_face_weights(cs_lnum_t                   f_id,
                        const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *quant,
                        const short int             loc_ids[],
                        double                      wf[])
{
  cs_lnum_t  ii;
  double  contrib, len;
  cs_real_3_t  un, cp;

  const cs_quant_t  pfq = quant->face[f_id];
  const double  invsurf = 1/pfq.meas;

  /* Reset weighting */
  for (ii = 0; ii < connect->n_max_vbyc; ii++) wf[ii] = 0;

  /* Compute a weight for each vertex of the current face */
  for (ii = connect->f2e->idx[f_id]; ii < connect->f2e->idx[f_id+1]; ii++) {

    cs_lnum_t  e_id = connect->f2e->col_id[ii];
    cs_quant_t  peq = quant->edge[e_id];
    cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
    cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];

    _lenunit3(peq.center, pfq.center, &len, &un);
    _cp3(un, peq.unitv, &cp);

    contrib = 0.25 * peq.meas * len * _n3(cp) * invsurf;

    wf[loc_ids[v1_id]] += contrib;
    wf[loc_ids[v2_id]] += contrib;

  } /* End of loop on face edges */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
