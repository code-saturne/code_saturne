/*============================================================================
 * Manage geometrical quantities needed in CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_cdo_quantities_algo_ccenter_t  cs_cdo_quantities_cc_algo =
  CS_CDO_QUANTITIES_SATURNE_CENTER;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

#if CS_CDO_QUANTITIES_DBG > 1 && defined(DEBUG) && !defined(NDEBUG)
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
  projq.pa  *=  cs_math_onesix;
  projq.pb  *= -cs_math_onesix;
  projq.pab *=  cs_math_one24;
  projq.pa2 *=  cs_math_onetwelve;
  projq.pb2 *= -cs_math_onetwelve;

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

#if CS_CDO_QUANTITIES_DBG > 1 && defined(DEBUG) && !defined(NDEBUG)
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute dual face normals (face crossed by primal edges).
 *         Given a cell c and an edge e, there are two faces attached to the
 *         couple (c, e)
 *         For a face f in c, the triplet (e, f, c) induces a triangle s(e,f,c)
 *         The dual face is the union of these two triangles.
 *         Storage is based on the c2e connectivity
 *
 * \param[in]      topo     pointer to a cs_cdo_connect_t structure
 * \param[in, out] quant    pointer to cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_edge_based_quantities(const cs_cdo_connect_t  *topo,
                               cs_cdo_quantities_t     *quant)
{
  const int n_edges = quant->n_edges;

  /* Sanity check */
  assert(topo->e2v != NULL);
  assert(topo->f2e != NULL && topo->f2c != NULL && topo->c2e != NULL);

  cs_real_t  *edge_center = NULL;

  /* Build edge centers and edge vectors */
  BFT_MALLOC(quant->edge_vector, 3*n_edges, cs_real_t);
  BFT_MALLOC(edge_center, 3*n_edges, cs_real_t);

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
        xe[k] = 0.5* ( xb[k] + xa[k]) ;
      }

    }
    else {

      for (int k = 0; k < 3; k++) {
        vect[k] = xa[k] - xb[k];
        xe[k] = 0.5* ( xb[k] + xa[k]) ;
      }

    }

  } /* End of loop on edges */

  /* Allocate and initialize array */
  BFT_MALLOC(quant->sface_normal, 6*topo->c2e->idx[quant->n_cells], cs_real_t);

  /* Manage openMP cache */
  short int  **parent_thread_array = NULL;
  BFT_MALLOC(parent_thread_array, cs_glob_n_threads, short int *);
  for (int i = 0; i < cs_glob_n_threads; i++)
    parent_thread_array[i] = NULL;

# pragma omp parallel default(none) \
  shared(quant, topo, parent_thread_array, edge_center, cs_glob_n_threads)
  { // OMP Block
    const cs_adjacency_t  *c2f = topo->c2f, *f2e = topo->f2e;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int t_id = omp_get_thread_num();
#else
    int t_id = 0;
#endif

    short int  *parent = parent_thread_array[t_id];
    BFT_MALLOC(parent, topo->n_max_ebyc, short int);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_lnum_t  *c2e_idx = topo->c2e->idx + c_id;
      const cs_lnum_t  *c2e_ids = topo->c2e->ids + c2e_idx[0];
      const short int  n_ec = c2e_idx[1] - c2e_idx[0];

      /* Get cell center */
      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      /* Initialize parent array */
      for (short int e = 0; e < n_ec; e++) parent[e] = 0;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];

        /* Compute xf -> xc */
        cs_real_3_t  xfxc;

        if (f_id < quant->n_i_faces) {
          for (int k = 0; k < 3; k++)
            xfxc[k] = quant->i_face_center[3*f_id+k] - xc[k];
        }
        else {
          const cs_lnum_t  bf_id = f_id - quant->n_i_faces;
          assert(bf_id > -1);
          for (int k = 0; k < 3; k++)
            xfxc[k] = quant->b_face_center[3*bf_id+k] - xc[k];
        }

        for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

          const cs_lnum_t  e_id = topo->f2e->ids[j];
          const cs_real_t  *xe = edge_center + 3*e_id;

          /* Compute the vectorial area for the triangle : xc, xf, xe */
          cs_real_3_t  tria_vect, xexc;
          for (int k = 0; k < 3; k++)
            xexc[k] = xc[k] - xe[k];
          cs_math_3_cross_product(xfxc, xexc, tria_vect);

          cs_nvec3_t  tria;
          cs_nvec3(tria_vect, &tria);
          cs_nvec3_t  edge = cs_quant_set_edge_nvec(e_id, quant);

          /* Find the corresponding local cell edge */
          short int e = n_ec;
          for (short int _e = 0; _e < n_ec; _e++) {
            if (c2e_ids[_e] == e_id) {
              e = _e;
              break;
            }
          }
          assert(e < n_ec);

          /* Portion of dual faces to consider */
          cs_real_t  *sface = quant->sface_normal + 6*(c2e_idx[0]+e);

          /* One should have (normal_tria, tangent_e) > 0 */
          const double  orient = _dp3(tria.unitv, edge.unitv);
          assert(fabs(orient) > 0);

          /* Store the computed data */
          if (orient < 0) {
            for (int k = 0; k < 3; k++)
              sface[3*parent[e] + k] = -0.5 * tria_vect[k];
          }
          else {
            for (int k = 0; k < 3; k++)
              sface[3*parent[e] + k] =  0.5 * tria_vect[k];
          }

#if defined(DEBUG) && !defined(NDEBUG)
          cs_nvec3_t  df_nvec;
          cs_nvec3(sface + 3*parent[e], &df_nvec);
          assert(fabs(_dp3(df_nvec.unitv, edge.unitv)) > 0);
#endif
          parent[e] += 1;

        } /* Loop on face edges */

      } /* Loop on cell faces */

#if defined(DEBUG) && !defined(NDEBUG)
      for (short int e = 0; e < n_ec; e++)
        if (parent[e] != 2)
          bft_error(__FILE__, __LINE__, 0,
                    " Connectivity error detected while building dual face"
                    " quantity for cell %d\n"
                    " Each edge should have 2 adjacent faces in a cell.\n"
                    " Here, there is (are) %d face(s).", c_id, parent[e]);
#endif
    } /* End of loop on cells */

    BFT_FREE(parent);

  } // End of OpenMP block

  BFT_FREE(parent_thread_array);
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
  /* Sanity checks */
  assert(topo->f2e != NULL && topo->f2c != NULL && topo->c2v != NULL);

  const cs_adjacency_t  *c2f = topo->c2f, *f2e = topo->f2e;

  /* Allocate and initialize arrays */
  BFT_MALLOC(quant->dcell_vol, topo->c2v->idx[quant->n_cells], double);

# pragma omp parallel for default(none) shared(quant, topo, c2f, f2e)   \
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
        assert(_v1 < n_vc && _v2 < n_vc);
        vol_vc[_v1] += pvol;
        vol_vc[_v2] += pvol;

      } // Loop on face edges

    } // Loop on cell faces

  }  // Loop on cells

}

/*----------------------------------------------------------------------------
 * Define the cs_quant_info_t structures related to cells, faces and edges
 * ---------------------------------------------------------------------------*/

static void
_compute_quant_info(cs_cdo_quantities_t     *quant)    /* In/out */
{
  assert(quant != NULL); // Sanity check

  /* Cell info (set default values) */
  quant->cell_info.h_min = quant->cell_info.meas_min = DBL_MAX;
  quant->cell_info.h_max = quant->cell_info.meas_max = -DBL_MAX;

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    const double  meas = quant->cell_vol[c_id];

    if (meas > quant->cell_info.meas_max) {
      quant->cell_info.meas_max = meas;
      quant->cell_info.h_max = cbrt(meas);
    }
    if (meas < quant->cell_info.meas_min) {
      quant->cell_info.meas_min = meas;
      quant->cell_info.h_min = cbrt(meas);
    }

  } // Loop on cells

  /* Face info (set default values) */
  quant->face_info.h_min = quant->face_info.meas_min = DBL_MAX;
  quant->face_info.h_max = quant->face_info.meas_max = -DBL_MAX;

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

  } // Loop on faces

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

  } // Loop on faces

  /* Edge info (set default values) */
  quant->edge_info.h_min = quant->edge_info.meas_min = DBL_MAX;
  quant->edge_info.h_max = quant->edge_info.meas_max = -DBL_MAX;

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

  } // Loop on edges

}

/*----------------------------------------------------------------------------
  Algorithm for computing mesh quantities : cell centers are computed as the
  vertex average over cell vertices. Other quantities are copied from those
  computed by Code_Saturne (default algorithm)
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

    /* Sanity checks */
    assert(ve - vs > 0);
    assert(ve - vs < SHRT_MAX);

    double  *xc = quant->cell_centers + 3*c_id;

    xc[0] = xc[1] = xc[2] = 0;
    for (cs_lnum_t jv = vs; jv < ve; jv++) {

      const cs_real_t  *xv = quant->vtx_coord + 3*c2v->ids[jv];

      xc[0] += xv[0];
      xc[1] += xv[1];
      xc[2] += xv[2];

    } // Loop on cell vertices

    const double  coef = 1./(ve-vs);
    for (int k = 0; k < 3; k++) xc[k] *= coef;

  } // Loop on cells

}

/*----------------------------------------------------------------------------
 * Algorithm for computing cell barycenters inspired from the article
 * "Fast and accurate computation of polyhedral mass properties"
 * Journal of Graphics, 1997 by Brian Mirtich
 *
 * Compute also : cell volumes.
 * ---------------------------------------------------------------------------*/

static void
_mirtich_algorithm(const cs_mesh_t             *mesh,
                   const cs_mesh_quantities_t  *mq,
                   const cs_cdo_connect_t      *connect,
                   cs_cdo_quantities_t         *quant) /* In/out */
{
  const int X = 0, Y = 1, Z = 2;
  const cs_lnum_t n_cells = quant->n_cells;
  const cs_lnum_t n_faces = quant->n_faces;

  /* Sanity checks */
  assert(connect->f2c != NULL);
  assert(connect->c2f != NULL);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_cells; i++) {
    quant->cell_vol[i] = 0.0;
    for (int k = 0; k < 3; k++)
      quant->cell_centers[3*i+k] = 0.0;
  }

  for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++) {   /* Loop on faces */

    /* Choose gamma to maximize normal according gamma (x, y, or z)
       Define a direct basis (alpha, beta, gamma) with this choice
       Compute omega = - <n, P> where P belongs to the face */

    _cdo_fspec_t  fspec = _get_fspec(f_id, mesh, mq);

    cs_lnum_t  A = fspec.XYZ[X];
    cs_lnum_t  B = fspec.XYZ[Y];
    cs_lnum_t  C = fspec.XYZ[Z];

    _cdo_fsubq_t  fsubq = _get_fsub_quantities(f_id,
                                               connect,
                                               quant->vtx_coord,
                                               fspec);

    /* Update cell quantities */
    const cs_adjacency_t  *f2c = connect->f2c;
    for (cs_lnum_t i = f2c->idx[f_id]; i < f2c->idx[f_id+1]; i++) {

      const cs_lnum_t  c_id = f2c->ids[i];
      const short int  sgn = f2c->sgn[i];
      const double Fvol = ( (fspec.XYZ[X] == 0) ? fsubq.Fa :
                            ( (fspec.XYZ[Y] == 0) ? fsubq.Fb : fsubq.Fc) );

      quant->cell_vol[c_id] += sgn * fspec.q.unitv[0] * Fvol;
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

/*----------------------------------------------------------------------------
 * Compute dual face normals (face crossed by primal edges).
 * Given a cell and an edge, there are two faces attached to the
 * couple (cell, edge)
 * The triplet (edge, face, cell) induces an elementary triangle s(e,f,c)
 * The dual face is the union of these two triangles.
 * Storage based on c2e connectivity
 * ---------------------------------------------------------------------------*/

static void
_define_cell_flag(const cs_cdo_connect_t  *topo,
                  cs_cdo_quantities_t     *cdoq)
{
  assert(topo->cell_type != NULL); /* Sanity check */

  const double  ortho_threshold = 1e-10;

  /* Allocate flag related to each cell */
  BFT_MALLOC(cdoq->cell_flag, cdoq->n_cells, cs_flag_t);
  for (cs_lnum_t i = 0; i < cdoq->n_cells; i++)
    cdoq->cell_flag[i] = 0;

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    if (topo->cell_type[c_id] == FVM_CELL_HEXA) { // Check if orthogonal

      double  ortho_crit = 0.;
      for (cs_lnum_t i = topo->c2f->idx[c_id]; i < topo->c2f->idx[c_id+1];
           i++) {

        const cs_lnum_t  f_id = topo->c2f->ids[i];
        const cs_nvec3_t  de_nv = cs_quant_set_dedge_nvec(i, cdoq);
        const cs_nvec3_t  pf_nv = cs_quant_set_face_nvec(f_id, cdoq);

        ortho_crit += fabs(1 - _dp3(de_nv.unitv, pf_nv.unitv));

      } // Loop on cell faces

      if (ortho_crit < ortho_threshold)
        cdoq->cell_flag[c_id] |= CS_CDO_ORTHO;

    } // Hexahedron

  } // Loop on cells

}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the triangle of base given by q (related to a
 *         segment) with apex located at xa
 *
 * \param[in]  qa   pointer to a cs_quant_t structure related to a segment
 * \param[in]  xb   coordinates of the apex to consider
 *
 * \return the value the area of the triangle
 */
/*----------------------------------------------------------------------------*/

double
cs_compute_area_from_quant(const cs_quant_t   qa,
                           const cs_real_t   *xb)
{
  double  xab[3], xab_un[3], cp[3];
  xab[0] = xb[0] - qa.center[0];
  xab[1] = xb[1] - qa.center[1];
  xab[2] = xb[2] - qa.center[2];

  const double  xab_len = cs_math_3_norm(xab);
  const double  inv_len = 1/xab_len;

  xab_un[0] = inv_len * xab[0];
  xab_un[1] = inv_len * xab[1];
  xab_un[2] = inv_len * xab[2];

  cs_math_3_cross_product(xab_un, qa.unitv, cp);

  /* tab = ||(qb.center -xa) x qb||/2 */
  return 0.5 * xab_len * qa.meas * cs_math_3_norm(cp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of algorithm to use for computing the cell center
 *
 * \param[in]  algo     type of algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_set_algo_ccenter(cs_cdo_quantities_algo_ccenter_t   algo)
{
  cs_cdo_quantities_cc_algo = algo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a cs_cdo_quantities_t structure
 *
 * \param[in]  m           pointer to a cs_mesh_t structure
 * \param[in]  mq          pointer to a cs_mesh_quantities_t structure
 * \param[in]  topo        pointer to a cs_cdo_connect_t structure
 *
 * \return  a new allocated pointer to a cs_cdo_quantities_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_quantities_t *
cs_cdo_quantities_build(const cs_mesh_t              *m,
                        const cs_mesh_quantities_t   *mq,
                        const cs_cdo_connect_t       *topo)
{
  cs_timer_t t0 = cs_timer_time();

  cs_cdo_quantities_t  *cdoq = NULL;

  /* Sanity checks */
  assert(topo != NULL);
  assert(topo->c2f != NULL);

  /* Build cs_cdo_quantities_t structure */
  BFT_MALLOC(cdoq, 1, cs_cdo_quantities_t);

  /* Dimension of each type of entities */
  const cs_lnum_t  n_cells = m->n_cells;
  cdoq->n_cells = n_cells;
  cdoq->n_g_cells = m->n_g_cells;

  /* Compute the volume of the whole domain */
  cdoq->vol_tot = mq->tot_vol;

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

  cdoq->n_vertices = m->n_vertices;
  cdoq->n_g_vertices = m->n_g_vertices;
  cdoq->vtx_coord = m->vtx_coord;

  /* Compute cell centers and cell volumes */
  BFT_MALLOC(cdoq->cell_centers, 3*n_cells, cs_real_t);
  BFT_MALLOC(cdoq->cell_vol, n_cells, cs_real_t);

  switch (cs_cdo_quantities_cc_algo) {

  case CS_CDO_QUANTITIES_MEANV_CENTER:
    cs_log_printf(CS_LOG_SETUP,
                  " -cdo- Cell.Center.Algo >> Vertices.MeanValue\n");

    /* Copy cell volumes */
    memcpy(cdoq->cell_vol, mq->cell_vol, n_cells*sizeof(cs_real_t));

    /* Compute the cell center */
    _vtx_algorithm(topo, cdoq);
    break;

  case CS_CDO_QUANTITIES_BARYC_CENTER:
    cs_log_printf(CS_LOG_SETUP,
                  " -cdo- Cell.Center.Algo >> Mirtich\n");
    /* Compute (real) the barycentric centers and cell volumes */
    _mirtich_algorithm(m, mq, topo, cdoq);
    break;


  case CS_CDO_QUANTITIES_SATURNE_CENTER:
    cs_log_printf(CS_LOG_SETUP,
                  " -cdo- Cell.Center.Algo >> Original\n");

    /* Copy cell volumes and cell centers */
    memcpy(cdoq->cell_centers, mq->cell_cen, 3*n_cells*sizeof(cs_real_t));
    memcpy(cdoq->cell_vol, mq->cell_vol, n_cells*sizeof(cs_real_t));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _("Unkwown algorithm for cell center computation\n"));

  } /* switch according to cs_cdo_quantities_cc_algo */

  if (topo->e2v != NULL) {
    cdoq->n_edges = topo->n_edges;

    if (cs_glob_n_ranks == 1)
      cdoq->n_g_edges = cdoq->n_edges;

    else { /* Compute the global number of edges */

      const cs_lnum_t  *e2v_ids = topo->e2v->ids;

      cs_gnum_t  *e2v_gnum = NULL;
      BFT_MALLOC(e2v_gnum, 2*cdoq->n_edges, cs_gnum_t);

#     pragma omp parallel for if (cdoq->n_edges > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < cdoq->n_edges; e++) {

        cs_gnum_t  v1 = m->global_vtx_num[e2v_ids[2*e]];
        cs_gnum_t  v2 = m->global_vtx_num[e2v_ids[2*e+1]];
        if (v1 < v2)
          e2v_gnum[2*e] = v1, e2v_gnum[2*e+1] = v2;
        else
          e2v_gnum[2*e+1] = v2, e2v_gnum[2*e+1] = v1;

      }

      cs_lnum_t  *order = NULL;
      BFT_MALLOC(order, cdoq->n_edges, cs_lnum_t);
      cs_order_gnum_allocated_s(NULL, e2v_gnum, 2, order, cdoq->n_edges);

      cs_gnum_t  *order_couples = NULL;
      BFT_MALLOC(order_couples, 2*cdoq->n_edges, cs_gnum_t);
#     pragma omp parallel for if (cdoq->n_edges > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < cdoq->n_edges; e++) {
        const cs_lnum_t  o_id = 2*order[e];
        order_couples[2*e] = e2v_gnum[o_id];
        order_couples[2*e+1] = e2v_gnum[o_id+1];
      }

      fvm_io_num_t *edge_io_num
        = fvm_io_num_create_from_adj_s(NULL, order_couples, cdoq->n_edges, 2);

      cdoq->n_g_edges = fvm_io_num_get_global_count(edge_io_num);

      /* Free memory */
      BFT_FREE(order);
      BFT_FREE(e2v_gnum);
      BFT_FREE(order_couples);
      fvm_io_num_destroy(edge_io_num);

    } // parallel run

  }
  else /* Not used by the numerical scheme */
    cdoq->n_edges = cdoq->n_g_edges = -1;

  /* Vertex-based schemes */
  cdoq->dedge_vector = NULL;

  /* Compute dual edge quantities */
  /* ---------------------------- */

  const cs_lnum_t  idx_size = topo->c2f->idx[n_cells];
  const cs_adjacency_t  *c2f = topo->c2f;

  BFT_MALLOC(cdoq->dedge_vector, 3*idx_size, cs_real_t);

# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t  *xc = cdoq->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->ids[i];
      const short int  sgn = c2f->sgn[i];
      const cs_lnum_t  bf_id = f_id - cdoq->n_i_faces;
      if (bf_id > -1) {

        const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;
        for (int k = 0; k < 3; k++)
          cdoq->dedge_vector[3*i+k] = sgn * (xf[k] - xc[k]);

      }
      else {

        const cs_real_t  *xf = cdoq->i_face_center + 3*f_id;
        for (int k = 0; k < 3; k++)
          cdoq->dedge_vector[3*i+k] = sgn * (xf[k] - xc[k]);

      }

    } // Loop on cell faces

  } /* End of loop on cells */

  /* Compute edge-related quantities if needed */
  cdoq->edge_vector = NULL;
  cdoq->sface_normal = NULL;

  if (cdoq->n_edges > 0)
    _compute_edge_based_quantities(topo, cdoq);

  /* Compute dual cell volume attached to each vertex in a cell */
  cdoq->dcell_vol = NULL;
  _compute_dcell_quantities(topo, cdoq);

  /* Define cs_quant_info_t structure */
  _compute_quant_info(cdoq);

  /* Define a flag for each cell */
  _define_cell_flag(topo, cdoq);

  /* Monitoring */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_t  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f s\n",
                "<CDO/Quantities> Runtime", time_count.wall_nsec*1e-9);

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

  /* Cell-related quantities */
  BFT_FREE(q->cell_centers);
  BFT_FREE(q->cell_vol);
  BFT_FREE(q->cell_flag);

  /* Face-related quantities */
  BFT_FREE(q->dedge_vector);

  /* Edge-related quantities */
  BFT_FREE(q->edge_vector);
  BFT_FREE(q->sface_normal);

  /* Vertex-related quantities */
  BFT_FREE(q->dcell_vol);

  /* vtx_coord is free when the structure cs_mesh_t is destroyed */

  BFT_FREE(q);

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
  cs_log_printf(CS_LOG_DEFAULT, "\n CDO mesh quantities information:\n");

  /* Information about activated flags */
  cs_gnum_t  n_ortho_cells = 0;
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    if (quant->cell_flag[c_id] & CS_CDO_ORTHO)
      n_ortho_cells += 1;
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_GNUM_TYPE, &n_ortho_cells);

  cs_log_printf(CS_LOG_DEFAULT,
                " --cdo-- n_ortho_cells  %9lu\n", n_ortho_cells);

  /* Output */
  double  h_min[3] = {quant->cell_info.h_min,
                      quant->face_info.h_min,
                      quant->edge_info.h_min};
  double  h_max[3] = {quant->cell_info.h_max,
                      quant->face_info.h_max,
                      quant->edge_info.h_max};

  if (cs_glob_n_ranks > 1) {
    cs_parall_min(3, CS_DOUBLE, h_min);
    cs_parall_max(3, CS_DOUBLE, h_max);
  }

  cs_log_printf(CS_LOG_DEFAULT,
                " --cdo-- h_cell  %6.4e %6.4e (min/max)\n"
                " --cdo-- h_face  %6.4e %6.4e (min/max)\n"
                " --cdo-- h_edge  %6.4e %6.4e (min/max)\n\n",
                h_min[0], h_max[0], h_min[1], h_max[1], h_min[2], h_max[2]);

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

  fprintf(fdump, " -cdoq- n_cells =    %d\n", cdoq->n_cells);
  fprintf(fdump, " -cdoq- n_faces =    %d\n", cdoq->n_faces);
  fprintf(fdump, " -cdoq- n_edges =    %d\n", cdoq->n_edges);
  fprintf(fdump, " -cdoq- n_vertices = %d\n", cdoq->n_vertices);
  fprintf(fdump, " -cdoq- Total volume = %.6e\n\n", cdoq->vol_tot);

  fprintf(fdump, "\n *** Cell Quantities ***\n");
  fprintf(fdump, "-msg- num.; volume ; center (3)\n");
  for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
    cs_lnum_t  p = 3*i;
    fprintf(fdump, " [%6d] | %12.8e | % -12.8e | % -12.8e |% -12.8e\n",
            i+1, cdoq->cell_vol[i], cdoq->cell_centers[p],
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
    fprintf(fdump, " -cdoq-  [%8d] | % -10.6e | % -10.6e | % -10.6e |"
            " % -10.6e |\n", i+1, e_vect.meas,
            e_vect.unitv[0], e_vect.unitv[1], e_vect.unitv[2]);
  }

  fclose(fdump);
  BFT_FREE(fname);
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

  const cs_lnum_t  bf_id = f_id - cdoq->n_i_faces;
  if (bf_id > -1) { // Border face

    q.meas = cdoq->b_face_surf[bf_id];

    cs_nvec3(cdoq->b_face_normal + 3*bf_id, &nv);
    for (int k = 0; k < 3; k++)
      q.unitv[k] = nv.unitv[k];

    const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;
    for (int k = 0; k < 3; k++)
      q.center[k] = xf[k];

  }
  else { // Interior face

    q.meas = cdoq->i_face_surf[f_id];

    cs_nvec3(cdoq->i_face_normal + 3*f_id, &nv);
    for (int k = 0; k < 3; k++)
      q.unitv[k] = nv.unitv[k];

    const cs_real_t  *xf = cdoq->i_face_center + 3*f_id;
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
  const cs_lnum_t  bf_id = f_id - cdoq->n_i_faces;
  if (bf_id > -1)  // Border face
    cs_nvec3(cdoq->b_face_normal + 3*bf_id, &nv);
  else             // Interior face
    cs_nvec3(cdoq->i_face_normal + 3*f_id, &nv);

  return nv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the two normalized vector associated to a primal edge
 *
 * \param[in]  e_id     id related to the face (f_id > n_i_face -> border face)
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
 * \param[in]  f_shift    position in c2f_idx
 * \param[in]  cdoq       pointer to a cs_cdo_quantities_t structure
 *
 * \return  a pointer to the dual edge normalized vector
 */
/*----------------------------------------------------------------------------*/

cs_nvec3_t
cs_quant_set_dedge_nvec(cs_lnum_t                     f_shift,
                        const cs_cdo_quantities_t    *cdoq)
{
  cs_nvec3_t  nv;
  cs_nvec3(cdoq->dedge_vector + 3*f_shift, &nv);

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

  fprintf(_f, " -cdoq-  [%8d] | % -10.6e | % -10.6e | % -10.6e | % -10.6e |"
          " % -10.6e | % -10.6e | % -10.6e\n", num, q.meas,
          q.unitv[0], q.unitv[1], q.unitv[2], q.center[0], q.center[1],
          q.center[2]);
}

/*----------------------------------------------------------------------------*/

#undef _dp3
#undef _n3

END_C_DECLS
