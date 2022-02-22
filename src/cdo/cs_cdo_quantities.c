/*============================================================================
 * Manage geometrical quantities needed in CDO schemes
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_array_reduce.h"
#include "cs_log.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_param_cdo.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define  CS_CDO_QUANTITIES_DBG  0  /* Switch off/on debug information */

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

static long long  cs_cdo_quantities_time = 0;

/* Store in a flag which quantities have to be computed */

cs_flag_t  cs_cdo_quantities_flag = 0;

/* Algorithm used to compute the cell center */

cs_cdo_quantities_cell_center_algo_t
cs_cdo_cell_center_algo = CS_CDO_QUANTITIES_BARYC_CENTER;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and initialize a \ref cs_mesh_quantities_t structure
 *
 * \return  a pointer to the newly allocated structure
 */
/*----------------------------------------------------------------------------*/

static cs_cdo_quantities_t *
_create_cdo_quantities(void)
{
  cs_cdo_quantities_t  *cdoq = NULL;

  /* Build cs_cdo_quantities_t structure */

  BFT_MALLOC(cdoq, 1, cs_cdo_quantities_t);

  cdoq->remove_boundary_faces = false;

  cdoq->vol_tot = 0.;

  /* Cell-based quantities */

  cdoq->cell_info.h_min = cdoq->cell_info.meas_min =  DBL_MAX;
  cdoq->cell_info.h_max = cdoq->cell_info.meas_max = -DBL_MAX;

  cdoq->n_cells = 0;
  cdoq->n_g_cells = 0;
  cdoq->cell_centers = NULL;
  cdoq->cell_vol = NULL;

  /* Face-based quantities */

  cdoq->face_info.h_min = cdoq->face_info.meas_min = DBL_MAX;
  cdoq->face_info.h_max = cdoq->face_info.meas_max = -DBL_MAX;

  cdoq->n_faces = cdoq->n_i_faces = cdoq->n_b_faces = 0;
  cdoq->n_g_faces = 0;
  cdoq->dedge_vector = NULL;
  cdoq->pvol_fc = NULL;

  /* Edge-based quantities */

  cdoq->edge_info.h_min = cdoq->edge_info.meas_min = DBL_MAX;
  cdoq->edge_info.h_max = cdoq->edge_info.meas_max = -DBL_MAX;

  cdoq->n_edges = 0;
  cdoq->n_g_edges = 0;
  cdoq->edge_vector = NULL;
  cdoq->pvol_ec = NULL;
  cdoq->dface_normal = NULL;

  /* Vertex-based quantities */

  cdoq->n_vertices = 0;
  cdoq->n_g_vertices = 0;
  cdoq->dcell_vol = NULL;

  /* Shared pointers are not initialized at this stage */

  return cdoq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function related to the Mirtich algorithm
 *         - Define an unitary normal to the current face
 *         - Compute omega = - <n, P> where P belongs to the face
 *         - Choose projection axis in order to maximize the projected area
 *         - Define a direct basis (alpha, beta, gamma) with this choice
 *
 * \param[in]  f_id      id of the face to treat
 * \param[in]  m         pointer to a cs_mesh_t structure
 * \param[in]  mq        pointer to a cs_mesh_quantities_t structure
 *
 * \return a _cdo_fspec_t structure storing the computed quantities
 */
/*----------------------------------------------------------------------------*/

static _cdo_fspec_t
_get_fspec(cs_lnum_t                    f_id,
           const cs_mesh_t             *m,
           const cs_mesh_quantities_t  *mq)
{
  const int  X = 0, Y = 1, Z = 2;

  cs_lnum_t  f, j, k, v, e, s;
  double  inv_n, nx, ny, nz;
  double  P[3]; /* Point belonging to the current face */

  _cdo_fspec_t  fspec =
    { .XYZ = {X, Y, Z},
      .omega = 0.,
      .q.meas = 0.,
      .q.unitv[0] = 0, .q.unitv[1] = 0, .q.unitv[2] = 0,
    };

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

#if CS_CDO_QUANTITIES_DBG > 1 && defined(DEBUG) && !defined(NDEBUG)
  printf("\n F: (%d, %d) >> surf: %e; omega: %e; XYZ: %d%d%d; [%e, %e, %e]\n",
         f_id, f, fspec.q.meas, fspec.omega,
         fspec.XYZ[0], fspec.XYZ[1], fspec.XYZ[2],
         fspec.q.unitv[0], fspec.q.unitv[1], fspec.q.unitv[2]);
#endif

  return fspec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function related to the Mirtich algorithm
 *         Compute projected integrals and quantities
 *
 * \param[in]  f_id      id of the face to treat
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  coords    coordinates of the mesh vertices
 * \param[in]  axis      local basis to the given face
 *
 * \return a _cdo_projq_t structure storing the computed quantities
 */
/*----------------------------------------------------------------------------*/

static _cdo_projq_t
_get_proj_quantities(cs_lnum_t                f_id,
                     const cs_cdo_connect_t  *connect,
                     const cs_real_t          coords[],
                     const int                axis[])
{
  cs_lnum_t  v_id[2];
  cs_real_t  a[2], b[2], a2[2], b2[2];

  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  /* Initialize structure */

  _cdo_projq_t  projq;

  /* These quantities are the integral of q on the plane
     (alpha, beta) where the face is projected */

  projq.p1 = projq.pa = projq.pb =  projq.pc = 0.0;
  projq.pab = projq.pa2 = projq.pb2 = 0.0;

  /* Scan edges which belong to the current face */

  for (cs_lnum_t  i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    short int  e_sgn = f2e->sgn[i];
    cs_lnum_t  e_id = f2e->ids[i];
    cs_lnum_t  s = 2*e_id;

    if (e_sgn > 0)
      v_id[0] = e2v->ids[s], v_id[1] = e2v->ids[s+1];
    else
      v_id[0] = e2v->ids[s+1], v_id[1] = e2v->ids[s];

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
  projq.pa  *=  cs_math_1ov6;
  projq.pb  *= -cs_math_1ov6;
  projq.pab *=  cs_math_1ov24;
  projq.pa2 *=  cs_math_1ov12;
  projq.pb2 *= -cs_math_1ov12;

  return  projq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Function related to the Mirtich algorithm
 *         Compute quantities related to the sub-faces of a given face
 *
 * \param[in]  f_id      id of the face to treat
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  coord     coordinates of the mesh vertices
 *
 *
 * \return a _cdo_fsubq_t structure storing the computed quantities
 */
/*----------------------------------------------------------------------------*/

static _cdo_fsubq_t
_get_fsub_quantities(cs_lnum_t                 f_id,
                     const cs_cdo_connect_t   *connect,
                     const cs_real_t          *coord,
                     const _cdo_fspec_t        fspec)
{
  const double  na = fspec.q.unitv[fspec.XYZ[0]];
  const double  nb = fspec.q.unitv[fspec.XYZ[1]];
  const double  nc = fspec.q.unitv[fspec.XYZ[2]];
  const double  k1 = 1./nc, k2 = k1 * k1, k3 = k2 * k1;

  /* Compute projected quantities */

  const _cdo_projq_t  projq = _get_proj_quantities(f_id, connect, coord, fspec.XYZ);

#if CS_CDO_QUANTITIES_DBG > 1 && defined(DEBUG) && !defined(NDEBUG)
  printf(" F: %d >> p1: %.4e, pa: %.4e, pb: %.4e, pc: %.4e,"
         " pab: %.4e, pa2: %.4e, pb2: %.4e\n",
         f_id, projq.p1, projq.pa, projq.pb, projq.pc,
         projq.pab, projq.pa2, projq.pb2);
#endif

  _cdo_fsubq_t  fsubq;

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute additional quantities related to faces like dual edges
 *         (segment between x_f and x_c and scanned with c2f adjacency)
 *
 * \param[in]      topo              pointer to a cs_cdo_connect_t structure
 * \param[in, out] cdoq              pointer to cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_face_based_quantities(const cs_cdo_connect_t  *topo,
                               cs_cdo_quantities_t     *cdoq)
{
  /* Compute dual edge quantities */

  const cs_lnum_t  n_cells = cdoq->n_cells;
  const cs_adjacency_t  *c2f = topo->c2f;

  BFT_MALLOC(cdoq->dedge_vector, 3*c2f->idx[n_cells], cs_real_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t  *xc = cdoq->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->ids[i];
      const short int  sgn = c2f->sgn[i];
      const cs_real_t  *xf = (f_id < cdoq->n_i_faces) ?
        cdoq->i_face_center + 3*f_id :
        cdoq->b_face_center + 3*(f_id - cdoq->n_i_faces);

      for (int k = 0; k < 3; k++)
        cdoq->dedge_vector[3*i+k] = sgn * (xf[k] - xc[k]);

    } /* Loop on cell faces */

  } /* End of loop on cells */

  /* Compute the volume of the pyramid */

  cs_flag_t  masks[3] = { CS_CDO_QUANTITIES_FB_SCHEME,
                          CS_CDO_QUANTITIES_HHO_SCHEME,
                          CS_CDO_QUANTITIES_VCB_SCHEME };

  if (cs_flag_at_least(cs_cdo_quantities_flag, 3, masks))
    cs_cdo_quantities_compute_pvol_fc(cdoq, topo->c2f, &(cdoq->pvol_fc));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute dual face normals (face crossed by primal edges).
 *         Given a cell c and an edge e, there are two faces attached to the
 *         couple (c, e)
 *         For a face f in c, the triplet (e, f, c) induces a triangle s(e,f,c)
 *         The dual face is the union of these two triangles.
 *         Storage is based on the c2e connectivity
 *
 * \param[in]      topo             pointer to a cs_cdo_connect_t structure
 * \param[in, out] quant            pointer to cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_edge_based_quantities(const cs_cdo_connect_t  *topo,
                               cs_cdo_quantities_t     *quant)
{
  assert(topo->e2v != NULL);
  assert(topo->f2e != NULL && topo->f2c != NULL && topo->c2e != NULL);

  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_edges = quant->n_edges;

  cs_real_t  *edge_center = NULL;

  /* Build edge centers and edge vectors */

  BFT_MALLOC(edge_center, 3*n_edges, cs_real_t);
  BFT_MALLOC(quant->edge_vector, 3*n_edges, cs_real_t);

# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

    /* Get the two vertex ids related to the current edge */

    const cs_lnum_t  *v_ids = topo->e2v->ids + 2*e_id;
    const cs_real_t  *xa = quant->vtx_coord + 3*v_ids[0];
    const cs_real_t  *xb = quant->vtx_coord + 3*v_ids[1];

    cs_real_t  *xe = edge_center + 3*e_id;
    cs_real_t  *vect = quant->edge_vector + 3*e_id;
    if (v_ids[1] > v_ids[0]) { /* vb > va */

      for (int k = 0; k < 3; k++) {
        vect[k] = xb[k] - xa[k];
        xe[k] = 0.5* (xb[k] + xa[k]) ;
      }

    }
    else {

      for (int k = 0; k < 3; k++) {
        vect[k] = xa[k] - xb[k];
        xe[k] = 0.5* (xb[k] + xa[k]) ;
      }

    }

  } /* End of loop on edges */

  cs_flag_t  masks[3] = { CS_CDO_QUANTITIES_EB_SCHEME,
                          CS_CDO_QUANTITIES_VB_SCHEME,
                          CS_CDO_QUANTITIES_VCB_SCHEME };
  if (!cs_flag_at_least(cs_cdo_quantities_flag, 3, masks)) {
    BFT_FREE(edge_center);
    return;
  }

  /* Allocate and initialize array
   * a) Compute the two vector areas composing each dual face
   * b) Compute the volume associated to each edge in a cell */

  BFT_MALLOC(quant->pvol_ec, topo->c2e->idx[n_cells], cs_real_t);
  BFT_MALLOC(quant->dface_normal, 3*topo->c2e->idx[n_cells], cs_real_t);

# pragma omp parallel shared(quant, topo, edge_center)
  { /* OMP Block */

    const cs_adjacency_t  *c2f = topo->c2f, *f2e = topo->f2e;

#     pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      const cs_lnum_t  *c2e_idx = topo->c2e->idx + c_id;
      const cs_lnum_t  *c2e_ids = topo->c2e->ids + c2e_idx[0];
      const short int  n_ec = c2e_idx[1] - c2e_idx[0];

      /* Initialize cell_dface */

      cs_real_t  *cell_dface = quant->dface_normal + 3*c2e_idx[0];
      memset(cell_dface, 0, 3*n_ec*sizeof(cs_real_t));

      /* Get the cell center */

      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];

        /* Compute xf -> xc */

        const cs_real_t  *xf = (f_id < quant->n_i_faces) ?
          quant->i_face_center + 3* f_id :
          quant->b_face_center + 3*(f_id - quant->n_i_faces);
        const cs_real_3_t  xfxc = { xf[0] - xc[0],
                                    xf[1] - xc[1],
                                    xf[2] - xc[2] };

        for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

          const cs_lnum_t  e_id = topo->f2e->ids[j];
          const cs_real_t  *xe = edge_center + 3*e_id;

          /* Compute the vectorial area for the triangle : xc, xf, xe */

          cs_real_3_t  tria_vect, xexc;
          for (int k = 0; k < 3; k++)
            xexc[k] = xc[k] - xe[k];
          cs_math_3_cross_product(xfxc, xexc, tria_vect);

          /* Find the corresponding local id for this cell edge */

          short int e = n_ec;
          for (short int _e = 0; _e < n_ec; _e++) {
            if (c2e_ids[_e] == e_id) {
              e = _e;
              break;
            }
          }
          CS_CDO_OMP_ASSERT(e < n_ec);

          /* One should have (tria.unitv, edge.unitv) > 0 */

          cs_nvec3_t  edge = cs_quant_set_edge_nvec(e_id, quant);
          cs_nvec3_t  tria;
          cs_nvec3(tria_vect, &tria);
          const double  orient = _dp3(tria.unitv, edge.unitv);
          CS_CDO_OMP_ASSERT(fabs(orient) > 0);

          /* Store this portion of dual face area at the right place */

          cs_real_t  *_dface = cell_dface + 3*e;
          if (orient < 0)
            for (int k = 0; k < 3; k++) _dface[k] -= 0.5 * tria_vect[k];
          else
            for (int k = 0; k < 3; k++) _dface[k] += 0.5 * tria_vect[k];

        } /* Loop on face edges */

      } /* Loop on cell faces */

      /* Compute pvol_ec */

      cs_real_t  *_pvol = quant->pvol_ec + c2e_idx[0];
      for (short int e = 0; e < n_ec; e++) {
        _pvol[e] = cs_math_1ov3 * _dp3(cell_dface + 3*e,
                                       quant->edge_vector + 3*c2e_ids[e]);
        CS_CDO_OMP_ASSERT(_pvol[e] > 0);
      }

    } /* End of loop on cells */

  } /* End of OpenMP block */

  BFT_FREE(edge_center);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute dual cell volumes related to primal vertices
 *         Storage is based on the c2v connectivity
 *
 * \param[in]      topo     pointer to a cs_cdo_connect_t structure
 * \param[in, out] quant    pointer to cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dcell_quantities(const cs_cdo_connect_t  *topo,
                          cs_cdo_quantities_t     *quant)
{
  cs_flag_t  masks[2] = { CS_CDO_QUANTITIES_VB_SCHEME |
                          CS_CDO_QUANTITIES_VCB_SCHEME };
  if (!cs_flag_at_least(cs_cdo_quantities_flag, 2, masks))
    return;

  assert(topo->f2e != NULL && topo->f2c != NULL && topo->c2v != NULL);

  const cs_adjacency_t  *c2f = topo->c2f, *f2e = topo->f2e;

  /* Allocate and initialize arrays */

  BFT_MALLOC(quant->dcell_vol, topo->c2v->idx[quant->n_cells], double);

# pragma omp parallel for shared(quant, topo, c2f, f2e) \
  CS_CDO_OMP_SCHEDULE
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Compute part of dual volume related to each primal cell */

    const cs_lnum_t  *c2v_idx = topo->c2v->idx + c_id;
    const cs_lnum_t  *c2v_ids = topo->c2v->ids + c2v_idx[0];
    const short int  n_vc = c2v_idx[1] - c2v_idx[0];
    const cs_real_t  *xc = quant->cell_centers + 3*c_id;

    /* Initialize */

    double  *vol_vc = quant->dcell_vol + c2v_idx[0];
    for (short int v = 0; v < n_vc; v++) vol_vc[v] = 0.0;

    for (cs_lnum_t jf = c2f->idx[c_id]; jf < c2f->idx[c_id+1]; jf++) {

      const cs_lnum_t  f_id = topo->c2f->ids[jf];
      const cs_lnum_t  bf_id = f_id - quant->n_i_faces;

      const cs_real_t  *xf;
      if (bf_id > -1)
        xf = quant->b_face_center + 3*bf_id;
      else
        xf = quant->i_face_center + 3*f_id;

      for (cs_lnum_t je = f2e->idx[f_id]; je < f2e->idx[f_id+1]; je++) {

        const cs_lnum_t  e_id = f2e->ids[je];
        const cs_lnum_t *v_ids = topo->e2v->ids + 2*e_id;
        const cs_lnum_t  v1_id = v_ids[0], v2_id = v_ids[1];
        const double pvol = 0.5 * cs_math_voltet(quant->vtx_coord + 3*v1_id,
                                                 quant->vtx_coord + 3*v2_id,
                                                 xf,
                                                 xc);

        /* Find the corresponding local cell edge */

        short int _v1 = n_vc, _v2 = n_vc;
        for (short int _v = 0; _v < n_vc; _v++) {
          if (c2v_ids[_v] == v1_id) _v1 = _v;
          if (c2v_ids[_v] == v2_id) _v2 = _v;
        }
        CS_CDO_OMP_ASSERT(_v1 < n_vc && _v2 < n_vc);
        vol_vc[_v1] += pvol;
        vol_vc[_v2] += pvol;

      } /* Loop on face edges */

    } /* Loop on cell faces */

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the cs_quant_info_t structures related to cells, faces and
 *        edges. For edges, it depends on scheme flags
 *
 * \param[in, out] quant      pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_quant_info(cs_cdo_quantities_t     *quant)
{
  if (quant == NULL)
    return;

  /* Cell info */

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const double  meas = quant->cell_vol[c_id];

    if (meas > quant->cell_info.meas_max) {
      quant->cell_info.meas_max = meas;
      quant->cell_info.h_max = cbrt(meas);
    }
    if (meas < quant->cell_info.meas_min) {
      quant->cell_info.meas_min = meas;
      quant->cell_info.h_min = cbrt(meas);
    }

  } /* Loop on cells */

  /* Face info */

  for (cs_lnum_t  f_id = 0; f_id < quant->n_i_faces; f_id++) {

    const cs_real_t  meas = quant->i_face_surf[f_id];

    if (meas > quant->face_info.meas_max) {
      quant->face_info.meas_max = meas;
      quant->face_info.h_max = sqrt(meas);
    }
    if (meas < quant->face_info.meas_min) {
      quant->face_info.meas_min = meas;
      quant->face_info.h_min = sqrt(meas);
    }

  } /* Loop on interior faces */

  for (cs_lnum_t  f_id = 0; f_id < quant->n_b_faces; f_id++) {

    const cs_real_t  meas = quant->b_face_surf[f_id];

    if (meas > quant->face_info.meas_max) {
      quant->face_info.meas_max = meas;
      quant->face_info.h_max = sqrt(meas);
    }
    if (meas < quant->face_info.meas_min) {
      quant->face_info.meas_min = meas;
      quant->face_info.h_min = sqrt(meas);
    }

  } /* Loop on border faces */

  /* Edge info */

  if (quant->edge_vector != NULL) {

    for (cs_lnum_t  e_id = 0; e_id < quant->n_edges; e_id++) {

      cs_nvec3_t  edge = cs_quant_set_edge_nvec(e_id, quant);

      if (edge.meas > quant->edge_info.meas_max) {
        quant->edge_info.meas_max = edge.meas;
        quant->edge_info.h_max = edge.meas;
      }
      if (edge.meas < quant->edge_info.meas_min) {
        quant->edge_info.meas_min = edge.meas;
        quant->edge_info.h_min = edge.meas;
      }

    } /* Loop on edges */

  }

  if (cs_glob_n_ranks > 1) { /* Synchronization across ranks */

    double  buf[12] =
      {  quant->cell_info.meas_max,  quant->cell_info.h_max,
        -quant->cell_info.meas_min, -quant->cell_info.h_min,
         quant->face_info.meas_max,  quant->face_info.h_max,
        -quant->face_info.meas_min, -quant->face_info.h_min,
         quant->edge_info.meas_max,  quant->edge_info.h_max,
        -quant->edge_info.meas_min, -quant->edge_info.h_min };

    cs_parall_max(12, CS_DOUBLE, buf);

    quant->cell_info.meas_max =  buf[0];
    quant->cell_info.h_max    =  buf[1];
    quant->cell_info.meas_min = -buf[2];
    quant->cell_info.h_min    = -buf[3];
    quant->face_info.meas_max =  buf[4];
    quant->face_info.h_max    =  buf[5];
    quant->face_info.meas_min = -buf[6];
    quant->face_info.h_min    = -buf[7];
    quant->edge_info.meas_max =  buf[8];
    quant->edge_info.h_max    =  buf[9];
    quant->edge_info.meas_min = -buf[10];
    quant->edge_info.h_min    = -buf[11];

  } /* Synchronization */
}

/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : cell centers are computed as the
  vertex average over cell vertices. Other quantities are copied from those
  computed by code_saturne (default algorithm)
  ----------------------------------------------------------------------------*/

static void
_vtx_algorithm(const cs_cdo_connect_t      *connect,
               cs_cdo_quantities_t         *quant) /* In/out */
{
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_adjacency_t  *c2v = connect->c2v;

  /* Compute cell centers */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_lnum_t  vs = c2v->idx[c_id];
    const cs_lnum_t  ve = c2v->idx[c_id+1];

    assert(ve - vs > 0);
    assert(ve - vs < SHRT_MAX);

    double  *xc = quant->cell_centers + 3*c_id;

    xc[0] = xc[1] = xc[2] = 0;
    for (cs_lnum_t jv = vs; jv < ve; jv++) {

      const cs_real_t  *xv = quant->vtx_coord + 3*c2v->ids[jv];

      xc[0] += xv[0];
      xc[1] += xv[1];
      xc[2] += xv[2];

    } /* Loop on cell vertices */

    const double  coef = 1./(ve-vs);
    for (int k = 0; k < 3; k++) xc[k] *= coef;

  } /* Loop on cells */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Algorithm for computing the real cell barycenters inspired from the
 *        article "Fast and accurate computation of polyhedral mass properties"
 *        Journal of Graphics, 1997 by Brian Mirtich
 *
 * \param[in]     m         pointer to a cs_mesh_t structure
 * \param[in]     mq        pointer to a cs_mesh_quantities_t structure
 * \param[in]     connect   pointer to a cs_cdo_connect_t structure
 * \param[in,out] quant     pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_mirtich_algorithm(const cs_mesh_t             *mesh,
                   const cs_mesh_quantities_t  *mq,
                   const cs_cdo_connect_t      *connect,
                   cs_cdo_quantities_t         *quant)
{
  const int  X = 0, Y = 1, Z = 2;
  const cs_lnum_t n_cells = quant->n_cells;
  const cs_lnum_t n_faces = quant->n_faces;

  assert(connect->f2c != NULL);
  assert(connect->c2f != NULL);

# pragma omp parallel for if (3*n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_cells; i++)
    quant->cell_centers[i] = 0.0;

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {   /* Loop on faces */

    /* Choose gamma to maximize normal according gamma (x, y, or z)
       Define a direct basis (alpha, beta, gamma) with this choice
       Compute omega = - <n, P> where P belongs to the face */

    _cdo_fspec_t  fspec = _get_fspec(f_id, mesh, mq);

    cs_lnum_t  A = fspec.XYZ[X];
    cs_lnum_t  B = fspec.XYZ[Y];
    cs_lnum_t  C = fspec.XYZ[Z];

    _cdo_fsubq_t  fsubq = _get_fsub_quantities(f_id,
                                               connect, quant->vtx_coord,
                                               fspec);

    /* Update cell quantities */

    const cs_adjacency_t  *f2c = connect->f2c;
    for (cs_lnum_t i = f2c->idx[f_id]; i < f2c->idx[f_id+1]; i++) {

      const cs_lnum_t  c_id = f2c->ids[i];
      const short int  sgn = f2c->sgn[i];

      quant->cell_centers[3*c_id + A] += sgn * fspec.q.unitv[A] * fsubq.Fa2;
      quant->cell_centers[3*c_id + B] += sgn * fspec.q.unitv[B] * fsubq.Fb2;
      quant->cell_centers[3*c_id + C] += sgn * fspec.q.unitv[C] * fsubq.Fc2;

    } /* End of loop on cell faces */

  } /* End of loop on faces */

  /* Compute cell center of gravity and total volume */

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    const double  inv_2vol = 0.5 / quant->cell_vol[i];
    for (int k = 0; k < 3; k++)
      quant->cell_centers[3*i+k] *= inv_2vol;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the time elapsed to build the cs_cdo_quantities_t structure
 *
 * \return the value of the time elapsed in ns
 */
/*----------------------------------------------------------------------------*/

long long
cs_cdo_quantities_get_time_perfo(void)
{
  return cs_cdo_quantities_time;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set which quantities have to be computed. Additionnal quantities
 *         are added to cs_cdo_quantities_flag (static variable)
 *
 * \param[in]  option_flag     flag to set geometrical quantities to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_set(cs_flag_t   option_flag)
{
  cs_cdo_quantities_flag |= option_flag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of algorithm to use for computing the cell center
 *
 * \param[in]  algo     type of algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_set_algo_ccenter(cs_cdo_quantities_cell_center_algo_t   algo)
{
  cs_cdo_cell_center_algo = algo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a cs_cdo_quantities_t structure. Some quantities are shared
 *        with the \ref cs_mesh_quantities_t structure and other are not
 *        built according to the given flags in cs_cdo_quantities_flag.
 *
 * \param[in]  m                 pointer to a cs_mesh_t structure
 * \param[in]  mq                pointer to a cs_mesh_quantities_t structure
 * \param[in]  topo              pointer to a cs_cdo_connect_t structure
 *
 * \return  a new allocated pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_quantities_t *
cs_cdo_quantities_build(const cs_mesh_t             *m,
                        const cs_mesh_quantities_t  *mq,
                        const cs_cdo_connect_t      *topo)
{
  cs_timer_t t0 = cs_timer_time();

  /* Sanity checks */

  assert(topo != NULL && topo->c2f != NULL);

  cs_cdo_quantities_t  *cdoq = _create_cdo_quantities();

  /* Compute the volume of the whole domain */

  cdoq->vol_tot = mq->tot_vol;

  /* Is there a removal of boundary faces to speed-up 2D computations */

  if (m->n_g_b_faces_all > m->n_g_b_faces) {

    cdoq->remove_boundary_faces = true;

    /* The computation of the cell barycenter with the Mirtich algorithm is
       false when boundary faces are removed. Switch to a correct one. */

    if (cs_cdo_cell_center_algo == CS_CDO_QUANTITIES_BARYC_CENTER)
      cs_cdo_cell_center_algo = CS_CDO_QUANTITIES_MEANV_CENTER;

  }

  /* 1) Initialize shared quantities */
  /*    ============================ */

  /* Face-related quantities */
  /* ----------------------- */

  /* Shared quantities related to faces (interior and border) */

  cdoq->n_i_faces = m->n_i_faces;
  cdoq->i_face_normal = mq->i_face_normal;
  cdoq->i_face_center = mq->i_face_cog;
  cdoq->i_face_surf = mq->i_face_surf;

  cdoq->n_b_faces = m->n_b_faces;
  cdoq->b_face_normal = mq->b_face_normal;
  cdoq->b_face_center = mq->b_face_cog;
  cdoq->b_face_surf = mq->b_face_surf;

  cdoq->n_faces = m->n_i_faces + m->n_b_faces;
  cdoq->n_g_faces = m->n_g_i_faces + m->n_g_b_faces;

  /* Vertex-related quantities */
  /* ------------------------- */

  cdoq->n_vertices = m->n_vertices;
  cdoq->n_g_vertices = m->n_g_vertices;
  cdoq->vtx_coord = m->vtx_coord;

  /* Edge-related quantities */
  /* ----------------------- */

  cdoq->n_edges = topo->n_edges;
  cdoq->n_g_edges = topo->n_g_edges;

  /* Cell-related quantities */
  /* ----------------------- */

  const cs_lnum_t  n_cells = m->n_cells;

  cdoq->n_cells = n_cells;
  cdoq->n_g_cells = m->n_g_cells;
  cdoq->cell_vol = mq->cell_vol;

  /* Compute the cell centers */

  switch (cs_cdo_cell_center_algo) {

  case CS_CDO_QUANTITIES_SATURNE_CENTER:
    cdoq->cell_centers = mq->cell_cen; /* shared */
    break;

  case CS_CDO_QUANTITIES_BARYC_CENTER:
    BFT_MALLOC(cdoq->cell_centers, 3*n_cells, cs_real_t);
    _mirtich_algorithm(m, mq, topo, cdoq);
    break;

  case CS_CDO_QUANTITIES_MEANV_CENTER:
    BFT_MALLOC(cdoq->cell_centers, 3*n_cells, cs_real_t);
    _vtx_algorithm(topo, cdoq);
    break;

  } /* Cell center algorithm */

  /* 2) Define quantities available for CDO schemes */
  /*    =========================================== */

  /* Face-related quantities */
  /* ----------------------- */

  _compute_face_based_quantities(topo, cdoq);

  /* Vertex-related quantities */
  /* ----------------------- */

  /* Compute dual cell volume attached to each vertex in a cell */
  _compute_dcell_quantities(topo, cdoq);

  /* Edge-related quantities */
  /* ----------------------- */

  _compute_edge_based_quantities(topo, cdoq);

  /* 3) Define metadata */
  /*    =============== */

  /* Define cs_quant_info_t structure */

  _compute_quant_info(cdoq);

  /* Monitoring */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);

  cs_cdo_quantities_time += time_count.nsec;

  return cdoq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdo_quantities_t structure
 *
 * \param[in]  cdoq      pointer to the cs_cdo_quantities_t struct. to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_quantities_t *
cs_cdo_quantities_free(cs_cdo_quantities_t   *cdoq)
{
  if (cdoq == NULL)
    return cdoq;

  /* Cell-related quantities */

  if (cs_cdo_cell_center_algo != CS_CDO_QUANTITIES_SATURNE_CENTER)
    BFT_FREE(cdoq->cell_centers);

  /* Face-related quantities */

  BFT_FREE(cdoq->dedge_vector);
  BFT_FREE(cdoq->pvol_fc);

  /* Edge-related quantities */

  BFT_FREE(cdoq->edge_vector);
  BFT_FREE(cdoq->dface_normal);
  BFT_FREE(cdoq->pvol_ec);

  /* Vertex-related quantities
   * vtx_coord is free when the structure cs_mesh_t is destroyed */

  BFT_FREE(cdoq->dcell_vol);

  BFT_FREE(cdoq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Summarize generic information about the cdo mesh quantities
 *
 * \param[in]  quant     pointer to cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_summary(const cs_cdo_quantities_t  *quant)
{
  cs_log_printf(CS_LOG_SETUP, "\n## CDO quantities settings\n");

  switch (cs_cdo_cell_center_algo) {

  case CS_CDO_QUANTITIES_SATURNE_CENTER:
    cs_log_printf(CS_LOG_SETUP, " * Cell.Center.Algo: Shared with FV\n");
    break;

  case CS_CDO_QUANTITIES_BARYC_CENTER:
    cs_log_printf(CS_LOG_SETUP, " * Cell.Center.Algo: Barycenter (Mirtich)\n");
    break;

  case CS_CDO_QUANTITIES_MEANV_CENTER:
    cs_log_printf(CS_LOG_SETUP, " * Cell.Center.Algo: Mean-value of vertices\n");
    break;

  } /* Cell center algorithm */

  /* Output */

  cs_log_printf(CS_LOG_DEFAULT, "\n CDO mesh quantities information:\n");
  cs_log_printf(CS_LOG_DEFAULT,
                " --cdo-- h_cell  %6.4e %6.4e (min/max)\n"
                " --cdo-- h_face  %6.4e %6.4e (min/max)\n",
                quant->cell_info.h_min, quant->cell_info.h_max,
                quant->face_info.h_min, quant->face_info.h_max);

  if (quant->edge_vector != NULL)
    cs_log_printf(CS_LOG_DEFAULT,
                  " --cdo-- h_edge  %6.4e %6.4e (min/max)\n",
                  quant->edge_info.h_min, quant->edge_info.h_max);
  else
    cs_log_printf(CS_LOG_DEFAULT, "\n");

#if CS_CDO_QUANTITIES_DBG > 0 && defined(DEBUG) && !defined(NDEBUG)
  cs_cdo_quantities_dump(quant);
#endif
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
  int  lname = strlen("DumpQuantities.dat") + 1;

  /* Define the name of the dump file */

  char *fname = NULL;
  if (cs_glob_n_ranks > 1) {
    lname += 6;
    BFT_MALLOC(fname, lname, char);
    sprintf(fname, "DumpQuantities.%05d.dat", cs_glob_rank_id);
  }
  else {
    BFT_MALLOC(fname, lname, char);
    sprintf(fname, "DumpQuantities.dat");
  }
  FILE  *fdump = fopen(fname, "w");

  if (cdoq == NULL) {
    fprintf(fdump, "Empty structure.\n");
    fclose(fdump);
    return;
  }

  fprintf(fdump, "\n Quantities structure: %p\n\n", (const void *)cdoq);

  fprintf(fdump, " -cdoq- n_cells =    %ld\n", (long)cdoq->n_cells);
  fprintf(fdump, " -cdoq- n_faces =    %ld\n", (long)cdoq->n_faces);
  fprintf(fdump, " -cdoq- n_edges =    %ld\n", (long)cdoq->n_edges);
  fprintf(fdump, " -cdoq- n_vertices = %ld\n", (long)cdoq->n_vertices);
  fprintf(fdump, " -cdoq- Total volume = %.6e\n\n", cdoq->vol_tot);

  fprintf(fdump, "\n *** Cell Quantities ***\n");
  fprintf(fdump, "-msg- num.; volume ; center (3)\n");
  for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
    cs_lnum_t  p = 3*i;
    fprintf(fdump, " [%6ld] | %12.8e | % -12.8e | % -12.8e |% -12.8e\n",
            (long)i+1, cdoq->cell_vol[i], cdoq->cell_centers[p],
            cdoq->cell_centers[p+1], cdoq->cell_centers[p+2]);
  }

  fprintf(fdump, "\n\n *** Interior Face Quantities ***\n");
  fprintf(fdump, "-msg- num. ; measure ; unitary vector (3) ; center (3)\n");
  for (cs_lnum_t f_id = 0; f_id < cdoq->n_i_faces; f_id++) {
    cs_quant_t  q = cs_quant_set_face(f_id, cdoq);
    cs_quant_dump(fdump, f_id+1, q);
  }

  fprintf(fdump, "\n\n *** Border   Face Quantities ***\n");
  fprintf(fdump, "-msg- num. ; measure ; unitary vector (3) ; center (3)\n");
  for (cs_lnum_t f_id = cdoq->n_i_faces; f_id < cdoq->n_faces; f_id++) {
    cs_quant_t  q = cs_quant_set_face(f_id, cdoq);
    cs_quant_dump(fdump, f_id-cdoq->n_i_faces+1, q);
  }

  fprintf(fdump, "\n\n *** Edge Quantities ***\n");
  fprintf(fdump, "-msg- num. ; measure ; unitary vector (3) ; center (3)\n");
  for (cs_lnum_t i = 0; i < cdoq->n_edges; i++) {
    const cs_nvec3_t  e_vect = cs_quant_set_edge_nvec(i, cdoq);
    fprintf(fdump, " -cdoq-  [%8ld] | % -10.6e | % -10.6e | % -10.6e |"
            " % -10.6e |\n", (long)i+1, e_vect.meas,
            e_vect.unitv[0], e_vect.unitv[1], e_vect.unitv[2]);
  }

  fclose(fdump);
  BFT_FREE(fname);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the portion of volume surrounding each face of a cell.
 *        This volume corresponds to a pyramid with base f and apex x_f
 *        The computed quantity is scanned with the c2f adjacency
 *
 * \param[in]      cdoq        pointer to cs_cdo_quantities_t structure
 * \param[in]      c2f         pointer to the cell --> edges connectivity
 * \param[in, out] p_pvol_fc   double pointer to the face volume in each cell
 *                             If not allocated before calling this function,
 *                             one allocates the array storing the volumes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_pvol_fc(const cs_cdo_quantities_t    *cdoq,
                                  const cs_adjacency_t         *c2f,
                                  cs_real_t                   **p_pvol_fc)
{
  if (cdoq == NULL || c2f == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: A mandatory structure is not allocated.\n", __func__);

  const cs_lnum_t  n_cells = cdoq->n_cells;

  cs_real_t  *pvol_fc = *p_pvol_fc;

  /* Initialize array */

  if (pvol_fc == NULL)
    BFT_MALLOC(pvol_fc, c2f->idx[n_cells], cs_real_t);

#if defined(DEBUG) && !defined(NDEBUG)
  memset(pvol_fc, 0, c2f->idx[n_cells]*sizeof(cs_real_t));
#endif

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->ids[j];
      const cs_nvec3_t  fp_nvec = cs_quant_set_face_nvec(f_id, cdoq);
      const cs_nvec3_t  ed_nvec = cs_quant_set_dedge_nvec(j, cdoq);

      cs_real_t  p_fc = _dp3(fp_nvec.unitv, ed_nvec.unitv);
      p_fc *= cs_math_1ov3 * fp_nvec.meas * ed_nvec.meas;

      pvol_fc[j] = p_fc;

    } /* Loop on cell faces */
  } /* Loop on cells */

  /* Return pointer */

  *p_pvol_fc = pvol_fc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the portion of volume surrounding each edge of a cell
 *        The computed quantity is scanned with the c2e adjacency
 *
 * \param[in]      cdoq        pointer to cs_cdo_quantities_t structure
 * \param[in]      c2e         pointer to the cell --> edges connectivity
 * \param[in, out] p_pvol_ec   double pointer to the edge volume in each cell
 *                             If not allocated before calling this function,
 *                             one allocates the array storing the volumes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_pvol_ec(const cs_cdo_quantities_t   *cdoq,
                                  const cs_adjacency_t        *c2e,
                                  cs_real_t                  **p_pvol_ec)
{
  if (cdoq == NULL || c2e == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: A mandatory structure is not allocated.\n", __func__);

  const cs_lnum_t  n_cells = cdoq->n_cells;

  cs_real_t  *pvol_ec = *p_pvol_ec;

  /* Initialize array */

  if (pvol_ec == NULL)
    BFT_MALLOC(pvol_ec, c2e->idx[n_cells], cs_real_t);

  if (cdoq->pvol_ec != NULL)
    memcpy(pvol_ec, cdoq->pvol_ec, c2e->idx[n_cells]*sizeof(cs_real_t));

  else {

#if defined(DEBUG) && !defined(NDEBUG)
    memset(pvol_ec, 0, c2e->idx[n_cells]*sizeof(cs_real_t));
#endif

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

      const cs_lnum_t  start = c2e->idx[c_id];
      const cs_real_t  *sface = cdoq->dface_normal + 3*start;
      const cs_lnum_t  *c2e_ids = c2e->ids + start;

      cs_real_t  *_pvol = pvol_ec + start;
      for (cs_lnum_t j = 0; j < c2e->idx[c_id+1]-start; j++)
        _pvol[j] = cs_math_1ov3 * _dp3(sface + 3*j,
                                       cdoq->edge_vector + 3*c2e_ids[j]);

    } /* Loop on cells */

  } /* Not already existing */

  /* Return pointer */

  *p_pvol_ec = pvol_ec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the dual volume surrounding each vertex
 *
 * \param[in]      cdoq       pointer to cs_cdo_quantities_t structure
 * \param[in]      c2v        pointer to the cell-->vertices connectivity
 * \param[in, out] dual_vol   dual volumes related to each vertex
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_dual_volumes(const cs_cdo_quantities_t   *cdoq,
                                       const cs_adjacency_t        *c2v,
                                       cs_real_t                   *dual_vol)
{
  if (dual_vol == NULL)
    return;

  assert(cdoq != NULL && c2v != NULL);

  memset(dual_vol, 0, cdoq->n_vertices*sizeof(cs_real_t));

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++)
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      dual_vol[c2v->ids[j]] += cdoq->dcell_vol[j];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the triangles with basis each edge of the face
 *         and apex the face center.
 *         Case of interior faces.
 *         Storage in agreement with the if2v adjacency structure
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]       f_id      interior face id
 * \param[in, out]  tef       quantities to compute (pre-allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_i_tef(const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *cdoq,
                                cs_lnum_t                     f_id,
                                cs_real_t                     tef[])
{
  if (tef == NULL) return;

  const cs_real_t  *xf = cdoq->i_face_center + 3*f_id;
  const cs_lnum_t  *idx = connect->if2v->idx + f_id;
  const cs_lnum_t  *ids = connect->if2v->ids + idx[0];
  const int  n_ef = idx[1] - idx[0]; /* n_ef = n_vf */

  cs_lnum_t  _v0, _v1;
  for (int  e = 0; e < n_ef; e++) { /* Compute */

    if (e < n_ef - 1)
      _v0 = ids[e], _v1 = ids[e+1];
    else
      _v0 = ids[n_ef-1], _v1 = ids[0];

    tef[e] = cs_math_surftri(cdoq->vtx_coord + 3*_v0,
                             cdoq->vtx_coord + 3*_v1,
                             xf);

  } /* Loop on face edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the triangles with basis each edge of the face
 *         and apex the face center.
 *         Case of boundary faces.
 *         Storage in agreement with the bf2v adjacency structure
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]       bf_id     border face id
 * \param[in, out]  tef       quantities to compute (pre-allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_b_tef(const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *cdoq,
                                cs_lnum_t                     bf_id,
                                cs_real_t                     tef[])
{
  if (tef == NULL) return;

  const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;
  const cs_lnum_t  *idx = connect->bf2v->idx + bf_id;
  const cs_lnum_t  *ids = connect->bf2v->ids + idx[0];
  const int  n_ef = idx[1] - idx[0]; /* n_ef = n_vf */

  cs_lnum_t  _v0, _v1;
  for (int  e = 0; e < n_ef; e++) { /* Compute */

    if (e < n_ef - 1)
      _v0 = ids[e], _v1 = ids[e+1];
    else
      _v0 = ids[n_ef-1], _v1 = ids[0];

    tef[e] = cs_math_surftri(cdoq->vtx_coord + 3*_v0,
                             cdoq->vtx_coord + 3*_v1,
                             xf);

  } /* Loop on face edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weight related to each vertex of a face. This weight
 *         ensures a 2nd order approximation if the face center is the face
 *         barycenter.
 *         Case of interior faces.
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]       f_id      interior face id
 * \param[in, out]  wvf       quantities to compute (pre-allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_i_wvf(const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *cdoq,
                                cs_lnum_t                     f_id,
                                cs_real_t                     wvf[])
{
  if (wvf == NULL) return;

  const cs_real_t  *xf = cdoq->i_face_center + 3*f_id;
  const cs_lnum_t  *idx = connect->if2v->idx + f_id;
  const cs_lnum_t  *ids = connect->if2v->ids + idx[0];
  const int  n_vf = idx[1] - idx[0];

  for (cs_lnum_t  v = 0; v < n_vf; v++) wvf[v] = 0.; /* Init */

  int  _v0, _v1;
  for (cs_lnum_t  v = 0; v < n_vf; v++) { /* Compute */

    if (v < n_vf - 1)
      _v0 = v, _v1 = v+1;
    else
      _v0 = n_vf-1, _v1 = 0;

    const double  tef = cs_math_surftri(cdoq->vtx_coord + 3*ids[_v0],
                                        cdoq->vtx_coord + 3*ids[_v1],
                                        xf);
    wvf[_v0] += tef;
    wvf[_v1] += tef;

  } /* Loop on face vertices */

  const cs_real_t  coef = 0.5/cdoq->i_face_surf[f_id];
  for (cs_lnum_t  v = 0; v < n_vf; v++) wvf[v] *= coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weight related to each vertex of a face. This weight
 *         ensures a 2nd order approximation if the face center is the face
 *         barycenter.
 *         Case of boundary faces.
 *
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq      pointer to a cs_cdo_quantities_t structure
 * \param[in]       bf_id     border face id
 * \param[in, out]  wvf       quantities to compute (pre-allocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_compute_b_wvf(const cs_cdo_connect_t       *connect,
                                const cs_cdo_quantities_t    *cdoq,
                                cs_lnum_t                     bf_id,
                                cs_real_t                     wvf[])
{
  if (wvf == NULL) return;

  const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;
  const cs_lnum_t  *idx = connect->bf2v->idx + bf_id;
  const cs_lnum_t  *ids = connect->bf2v->ids + idx[0];
  const int  n_vf = idx[1] - idx[0];

  for (cs_lnum_t  v = 0; v < n_vf; v++) wvf[v] = 0.; /* Init */

  int  _v0, _v1;
  for (cs_lnum_t  v = 0; v < n_vf; v++) { /* Compute */

    if (v < n_vf - 1)
      _v0 = v, _v1 = v+1;
    else
      _v0 = n_vf-1, _v1 = 0;

    const double  tef = cs_math_surftri(cdoq->vtx_coord + 3*ids[_v0],
                                        cdoq->vtx_coord + 3*ids[_v1],
                                        xf);
    wvf[_v0] += tef;
    wvf[_v1] += tef;

  } /* Loop on face vertices */

  const cs_real_t  coef = 0.5/cdoq->b_face_surf[bf_id];
  for (cs_lnum_t  v = 0; v < n_vf; v++) wvf[v] *= coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a cs_quant_t structure for a primal face (interior or border)
 *
 * \param[in]  f_id     id related to the face (f_id > n_i_face -> border face)
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a initialize structure
 */
/*----------------------------------------------------------------------------*/

cs_quant_t
cs_quant_set_face(cs_lnum_t                    f_id,
                  const cs_cdo_quantities_t   *cdoq)
{
  cs_quant_t  q = {.meas = 0.,
                   .unitv[0] = 0., .unitv[1] = 0., .unitv[2] = 0.,
                   .center[0] = 0., .center[1] = 0., .center[2] = 0.};
  cs_nvec3_t  nv;

  if (f_id < cdoq->n_i_faces) { /* Interior face */

    q.meas = cdoq->i_face_surf[f_id];

    cs_nvec3(cdoq->i_face_normal + 3*f_id, &nv);
    for (int k = 0; k < 3; k++)
      q.unitv[k] = nv.unitv[k];

    const cs_real_t  *xf = cdoq->i_face_center + 3*f_id;
    for (int k = 0; k < 3; k++)
      q.center[k] = xf[k];

  }
  else { /* Border face */

    const cs_lnum_t  bf_id = f_id - cdoq->n_i_faces;

    q.meas = cdoq->b_face_surf[bf_id];
    cs_nvec3(cdoq->b_face_normal + 3*bf_id, &nv);
    for (int k = 0; k < 3; k++)
      q.unitv[k] = nv.unitv[k];

    const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;
    for (int k = 0; k < 3; k++)
      q.center[k] = xf[k];

  }

  return q;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the face surface and its unit normal vector for a primal
 *        face (interior or border)
 *
 * \param[in]  f_id     id related to the face (f_id > n_i_face -> border face)
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to the face normalized vector
 */
/*----------------------------------------------------------------------------*/

cs_nvec3_t
cs_quant_set_face_nvec(cs_lnum_t                    f_id,
                       const cs_cdo_quantities_t   *cdoq)
{
  cs_nvec3_t  nv;

  if (f_id < cdoq->n_i_faces)  /* Interior face */
    cs_nvec3(cdoq->i_face_normal + 3*f_id, &nv);
  else  /* Border face */
    cs_nvec3(cdoq->b_face_normal + 3*(f_id - cdoq->n_i_faces), &nv);

  return nv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the normalized vector associated to a primal edge
 *
 * \param[in]  e_id     id related to an edge
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return  a pointer to the edge normalized vector
 */
/*----------------------------------------------------------------------------*/

cs_nvec3_t
cs_quant_set_edge_nvec(cs_lnum_t                    e_id,
                       const cs_cdo_quantities_t   *cdoq)
{
  cs_nvec3_t  nv;
  cs_nvec3(cdoq->edge_vector + 3*e_id, &nv);

  return nv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the two normalized vector associated to a dual edge
 *
 * \param[in]  shift    position in c2f_idx
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return  a pointer to the dual edge normalized vector
 */
/*----------------------------------------------------------------------------*/

cs_nvec3_t
cs_quant_set_dedge_nvec(cs_lnum_t                     shift,
                        const cs_cdo_quantities_t    *cdoq)
{
  cs_nvec3_t  nv;
  cs_nvec3(cdoq->dedge_vector + 3*shift, &nv);

  return nv;
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

  fprintf(_f, " -cdoq-  [%8ld] | % -10.6e | % -10.6e | % -10.6e | % -10.6e |"
          " % -10.6e | % -10.6e | % -10.6e\n", (long)num, q.meas,
          q.unitv[0], q.unitv[1], q.unitv[2], q.center[0], q.center[1],
          q.center[2]);
}

/*----------------------------------------------------------------------------*/

#undef _dp3
#undef _n3

END_C_DECLS
