#ifndef __CS_CDO_QUANTITIES_H__
#define __CS_CDO_QUANTITIES_H__

/*============================================================================
 * Manage geometrical quantities needed in CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_flag.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Information useful to get simpler algo. */
#define CS_CDO_ORTHO      (1 << 0)  // Orthogonality condition is checked

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of algorithm used to compute the cell center */
typedef enum {

  /* Center is computed as the mean of cell vertices */
  CS_CDO_QUANTITIES_MEANV_CENTER,

  /* Center is computed as the real cell barycenter */
  CS_CDO_QUANTITIES_BARYC_CENTER,

  /* Use the cell center computed in cs_mesh_quantities.c (Default behavior) */
  CS_CDO_QUANTITIES_SATURNE_CENTER,

  CS_CDO_QUANTITIES_N_CENTER_ALGOS

} cs_cdo_quantities_algo_ccenter_t;

/* Structure storing information about variation of entities across the
   mesh for a given type of entity (cell, face and edge) */
typedef struct {

  /* Measure is either a volume for cells, a surface for faces or a length
     for edges */
  double   meas_min;  /* Min. value of the entity measure */
  double   meas_max;  /* Max. value of the entity measure */
  double   h_min;     /* Estimation of the min. value of the diameter */
  double   h_max;     /* Estimation of the max. value of the diameter */

} cs_quant_info_t;

/* For primal vector quantities (edge or face) */
typedef struct {

  double  meas;       /* length or area */
  double  unitv[3];   /* unitary vector: tangent or normal to the element */
  double  center[3];

} cs_quant_t;

typedef struct { /* Specific mesh quantities */

  /* Global mesh quantities */
  double           vol_tot;

  /* Cell-based quantities */
  /* ===================== */

  cs_lnum_t         n_cells;        /* Local number of cells */
  cs_gnum_t         n_g_cells;      /* Global number of cells */
  cs_real_t        *cell_centers;
  cs_real_t        *cell_vol;
  cs_flag_t        *cell_flag;      /* Flag attached to cell to associate
                                       metadata like boundary cell or
                                       orthogonality */

  cs_quant_info_t   cell_info;

  /* Face-based quantities */
  /* ===================== */

  cs_lnum_t         n_faces;        /* n_i_faces + n_b_faces */
  cs_gnum_t         n_g_faces;      /* Global number of faces */

  /* cs_quant_t structure attached to a face (interior or border) can be built
     on-the-fly: cs_quant_get(flag, f_id, quant). In order to reduce the memory
     consumption one shares face quantities with the ones defined in the legacy
     part and stored in the cs_mesh_quantities_t structure */

  cs_lnum_t         n_i_faces;      /* Local number of interior faces */
  const cs_real_t  *i_face_normal;  /* Shared with cs_mesh_quantities_t */
  const cs_real_t  *i_face_center;  /* Shared with cs_mesh_quantities_t */
  const cs_real_t  *i_face_surf;    /* Shared with cs_mesh_quantities_t */

  cs_lnum_t         n_b_faces;      /* Local number of border faces */
  const cs_real_t  *b_face_normal;  /* Shared with cs_mesh_quantities_t */
  const cs_real_t  *b_face_center;  /* Shared with cs_mesh_quantities_t */
  const cs_real_t  *b_face_surf;    /* Shared with cs_mesh_quantities_t */

  /* cs_nvec3_t structure attached to a dual edge can be built on-the-fly to
     access to its length and its unit tangential vector. One
     recalls that a dual edge is associated to a primal face and is shared with
     two cells if the face lies inside the domain and shared with one cell if
     the face lies on the domain boundary.
     Scan this quantity with the c2f connectivity
  */

  cs_real_t        *dedge_vector;   /* Allocation to 3*c2f->idx[n_faces] */

  cs_quant_info_t   face_info;

  /* Edge-based quantities */
  /* ===================== */

  cs_lnum_t         n_edges;        /* Local number of edges */
  cs_gnum_t         n_g_edges;      /* Global number of edges */

  cs_real_t        *edge_vector;    /* Allocation to 3*n_edges
                                       Norm of the vector is equal to the
                                       distance between two vertices.
                                       Unit vector is the tangential direction
                                       attached to the edge */

  /* For each edge e belonging to a cell c, two contributions coming from 2
     triangles  s(x_c, x_f, x_e) for a face f in Face_e are considered.
     Scan this quantity with the c2e connectivity */

  cs_real_t        *sface_normal;   /* 2 triangle-face normals by edge in a
                                       cell */

  cs_quant_info_t   edge_info;

  /* Vertex-based quantities */
  /* ======================= */

  cs_lnum_t         n_vertices;      /* Local number of vertices */
  cs_gnum_t         n_g_vertices;    /* Global number of vertices */

  cs_real_t        *dcell_vol;       /* Dual volume related to each vertex.
                                        Scan with the c2v connectivity */
  const cs_real_t  *vtx_coord;       /* Shared with the cs_mesh_t structure */

} cs_cdo_quantities_t;

/*============================================================================
 * Global variables
 *============================================================================*/

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

static inline double
cs_compute_area_from_quant(const cs_quant_t    qa,
                           const cs_real_t    *xb)
{
  const double  xab[3] = {xb[0] - qa.center[0],
                          xb[1] - qa.center[1],
                          xb[2] - qa.center[2]};
  const double  cp[3] = {qa.unitv[1]*xab[2] - qa.unitv[2]*xab[1],
                         qa.unitv[2]*xab[0] - qa.unitv[0]*xab[2],
                         qa.unitv[0]*xab[1] - qa.unitv[1]*xab[0]};

  return 0.5 * qa.meas * cs_math_3_norm(cp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of algorithm to use for computing the cell center
 *
 * \param[in]  algo     type of algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_set_algo_ccenter(cs_cdo_quantities_algo_ccenter_t   algo);

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
                        const cs_cdo_connect_t      *topo);

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
cs_cdo_quantities_free(cs_cdo_quantities_t   *q);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Summarize generic information about the cdo mesh quantities
 *
 * \param[in]  quant   pointer to cs_cdo_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_summary(const cs_cdo_quantities_t  *quant);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_cdo_quantities_t structure
 *
 * \param[in]  cdoq   pointer to cs_cdo_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_dump(const cs_cdo_quantities_t  *cdoq);

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
                                       cs_real_t                   *dual_vol);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the area of the triangles with basis each edge of the face
 *         and apex the face center.
 *         Case of interior faces.
 *         Storage in agreement with the bf2v adjacency structure
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
                                cs_real_t                     tef[]);

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
                                cs_real_t                     tef[]);

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
                                cs_real_t                     wvf[]);

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
                                cs_real_t                     wvf[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the face vector which the face_area * face_normal for a
 *        primal face (interior or border)
 *
 * \param[in]  f_id     id related to the face (f_id > n_i_face -> border face)
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to the face vector
 */
/*----------------------------------------------------------------------------*/

inline static const cs_real_t *
cs_quant_get_face_vector_area(cs_lnum_t                    f_id,
                              const cs_cdo_quantities_t   *cdoq)
{
  if (f_id < cdoq->n_i_faces)   /* Interior face */
    return cdoq->i_face_normal + 3*f_id;
  else                          /* Border face */
    return cdoq->b_face_normal + 3*(f_id - cdoq->n_i_faces);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the face center for a primal face (interior or border)
 *
 * \param[in]  f_id     id related to the face (f_id > n_i_face -> border face)
 * \param[in]  cdoq     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to the face center coordinates
 */
/*----------------------------------------------------------------------------*/

inline static const cs_real_t *
cs_quant_get_face_center(cs_lnum_t                    f_id,
                         const cs_cdo_quantities_t   *cdoq)
{
  if (f_id < cdoq->n_i_faces)   /* Interior face */
    return cdoq->i_face_center + 3*f_id;
  else                          /* Border face */
    return cdoq->b_face_center + 3*(f_id - cdoq->n_i_faces);
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
                  const cs_cdo_quantities_t   *cdoq);

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
                       const cs_cdo_quantities_t   *cdoq);

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
                       const cs_cdo_quantities_t   *cdoq);

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
                        const cs_cdo_quantities_t    *cdoq);

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
              const cs_quant_t  q);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_QUANTITIES_H__ */
