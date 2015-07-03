#ifndef __CS_CDO_QUANTITIES_H__
#define __CS_CDO_QUANTITIES_H__

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_cdo.h"
#include "cs_cdo_connect.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Type of algorithm to compute geometrical quantities */
typedef enum {

  CS_CDO_CC_NONE,   /* Default */
  CS_CDO_CC_MEANV,  /* Cell center is computed as the mean of cell vertices */
  CS_CDO_CC_BARYC,  /* Cell center is computed as the real cell barycenter */
  CS_CDO_CC_SATUR,  /* Cell center is given by Code_Saturne */
  CS_CDO_CC_ORTHO,  /* Cell center is optimized to enforce orthogonality
                         between cell-face edge and face plane */
  CS_CDO_N_CC_ALGOS

} cs_cdo_cc_algo_t;

typedef struct {

  double  meas;
  double  unitv[3];

} cs_qvect_t;

/* For primal vector quantities (edge or face) */
typedef struct {

  double  meas;       /* length or area */
  double  unitv[3];   /* unitary vector: tangent or normal to the element */
  double  center[3];

} cs_quant_t;

/* For dual face quantities. Add also a link to entity id related to this dual
   quantities */

typedef struct { /* TODO: remove what is the less necessary in order to
                    save memory comsumption */

  cs_lnum_t  parent_id[2];  /* parent entity id of (primal) faces */
  double     meas[2];       /* length or area related to each parent */
  double     unitv[6];      /* unitary vector related to each parent */
  double     vect[3];       /* dual face vector (TO REMOVE) */

} cs_dface_t;

typedef struct { /* Specific mesh quantities */

  /* Global mesh quantities */
  double      vol_tot;

  /* Vertex-based quantities */
  cs_lnum_t   n_vertices;
  double     *dcell_vol;    /* dual volume related to each vertex.
                               Scan with the c2v connectivity */

  /* Cell-based quantities */
  cs_lnum_t   n_cells;
  cs_real_t  *cell_centers;
  cs_real_t  *cell_vol;

  /* Face-based quantities */
  cs_lnum_t     n_i_faces;
  cs_lnum_t     n_b_faces;
  cs_lnum_t     n_faces;    /* n_i_faces + n_b_faces */
  cs_quant_t   *face;       /* face quantities */
  double       *dedge;      /* dual edge quantities. 4 values by entry
                               length and unitv[3].
                               Scan with the c2f connectivity */

  /* Edge-based quantities */
  cs_lnum_t     n_edges;
  cs_quant_t   *edge;       /* edge quantities */
  cs_dface_t   *dface;      /* for each edge belonging to a cell, two
                               contributions coming from 2 triangles
                               s(x_cell, x_face, x_edge) for face\in\Face_edge
                               are considered.
                               Scan with the c2e connectivity */

} cs_cdo_quantities_t;

/*============================================================================
 * Global variables
 *============================================================================*/

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
 * \brief Dump a cs_cdo_quantities_t structure
 *
 * \param[in]  iq     pointer to cs_cdo_quantities_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_quantities_dump(const cs_cdo_quantities_t  *iq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_quant_t structure
 *
 * \param[in]  f         FILE struct (stdout if NULL)
 * \param[in]  num       entity number related to this quantity struct.
 * \param[in]  quant     cs_quant_t structure to dump
 */
/*----------------------------------------------------------------------------*/

void
cs_quant_dump(FILE             *f,
              cs_lnum_t         num,
              const cs_quant_t  q);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex the dual cell volume which is also
 *
 *                sum    |celld(v) cap c| = pvol_v
 *              c\in\C_v
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantites_t structure
 * \param[inout] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_vtx(const cs_cdo_connect_t     *connect,
                    const cs_cdo_quantities_t  *quant,
                    double                     *p_pvol[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each edge a related volume pvol_e which constitutes
 *         a partition of unity
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantites_t structure
 * \param[inout] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_edge(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     double                      *p_pvol[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each face a related volume pvol_f which constitutes
 *         a partition of unity
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantites_t structure
 * \param[inout] p_pvol    pvol (if NULL, allocated in this routine)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_pvol_face(const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     double                     *p_pvol[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each face a weight related to each vertex w_{v,f}
 *         This weight is equal to |dc(v) cap f|/|f| so that the sum of the
 *         weights is equal to 1.
 *
 * \param[in]    f_id      id of the face
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantites_t structure
 * \param[in]    loc_ids   indirection to a local numbering
 * \param[inout] wf        already allocated to n_max_vbyc (reset)
 */
/*----------------------------------------------------------------------------*/

void
cs_compute_face_weights(cs_lnum_t                   f_id,
                        const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *quant,
                        const short int             loc_ids[],
                        double                      wf[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_QUANTITIES_H__ */
