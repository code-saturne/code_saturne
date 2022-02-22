#ifndef __CS_CDO_CONNECT_H__
#define __CS_CDO_CONNECT_H__

/*============================================================================
 * Manage connectivity (Topological features of the mesh)
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

#include "cs_base.h"
#include "cs_flag.h"
#include "cs_param_types.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Additional macros */

#define CS_TRIANGLE_CASE          3 /* Number of vertices in a triangle */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Connectivity structure */

typedef struct {

  /* Vertex-related members */
  /* ---------------------- */

  cs_lnum_t          n_vertices;

  /* range set and interface set structures for scalar-valued vertex DoFs
     (these structures may be not allocated according to the settings) */

  cs_range_set_t      *vtx_rset;
  cs_interface_set_t  *vtx_ifs;


  /* Edge-related members */
  /* -------------------- */

  cs_lnum_t          n_edges;
  cs_gnum_t          n_g_edges;

  cs_adjacency_t      *e2v;             /* edge --> vertices connectivity */

  /* range set and interface set structures for scalar-valued edge DoFs
     (these structures may be not allocated according to the settings) */

  cs_range_set_t      *edge_rset;
  cs_interface_set_t  *edge_ifs;

  /* Face-related members */
  /* -------------------- */

  cs_lnum_t          n_faces[3];  /* 0: all, 1: border, 2: interior */

  cs_adjacency_t    *f2c;         /* face --> cells connectivity */
  cs_adjacency_t    *f2e;         /* face --> edges connectivity */
  cs_adjacency_t    *bf2v;        /* border face --> vertices connectivity
                                     (map from cs_mesh_t) */
  cs_adjacency_t    *if2v;        /* interior face --> vertices connectivity
                                     (map from cs_mesh_t) */

  /* range set and interface set structures for scalar-valued edge DoFs
     (these structures may be not allocated according to the settings) */

  cs_range_set_t      *face_rset;
  cs_interface_set_t  *face_ifs;

  /* Cell-related members */
  /* -------------------- */

  cs_lnum_t          n_cells;

  fvm_element_t     *cell_type;   /* type of cell */
  cs_flag_t         *cell_flag;   /* Flag (Border/Solid) */
  cs_adjacency_t    *c2f;         /* cell --> faces connectivity */
  cs_adjacency_t    *c2e;         /* cell --> edges connectivity */
  cs_adjacency_t    *c2v;         /* cell --> vertices connectivity */

  /* Delta of ids between the min./max. values of entities related to a cell
     Useful to store compactly the link between mesh ids and cell mesh ids
     needed during the cell mesh definition */

  cs_lnum_t  e_max_cell_range;
  cs_lnum_t  v_max_cell_range;

  /* Max. connectitivy size for cells */

  int  n_max_vbyc;   /* max. number of vertices in a cell */
  int  n_max_ebyc;   /* max. number of edges in a cell */
  int  n_max_fbyc;   /* max. number of faces in a cell */
  int  n_max_vbyf;   /* max. number of vertices in a face */
  int  n_max_v2fc;   /* max. number of faces connected to a vertex in a cell */
  int  n_max_v2ec;   /* max. number of edges connected to a vertex in a cell */

  /* Adjacency related to linear systems (allocated only if needed) */

  cs_adjacency_t       *v2v;    /* vertex to vertices through cells */
  cs_adjacency_t       *f2f;    /* face to faces through cells */
  cs_adjacency_t       *e2e;    /* edge to edges through cells */

} cs_cdo_connect_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Static Inline Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the next three vertices in a row from a face to edge connectivity
 *         and a edge to vertex connectivity
 *
 * \param[in]       f2e_ids     face-edge connectivity
 * \param[in]       e2v_ids     edge-vertex connectivity
 * \param[in]       start_idx   index from which the current cell infos start
 * \param[in, out]  v0          id of the first vertex
 * \param[in, out]  v1          id of the second vertex
 * \param[in, out]  v2          id of the third vertex
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_connect_get_next_3_vertices(const cs_lnum_t   *f2e_ids,
                               const cs_lnum_t   *e2v_ids,
                               const cs_lnum_t    start_idx,
                               cs_lnum_t         *v0,
                               cs_lnum_t         *v1,
                               cs_lnum_t         *v2)
{
  const cs_lnum_t _2e0  = 2*f2e_ids[start_idx],
                  _2e1  = 2*f2e_ids[start_idx+1];
  const cs_lnum_t tmp = e2v_ids[_2e1];

  *v0 = e2v_ids[_2e0];
  *v1 = e2v_ids[_2e0+1];
  *v2 = ((tmp != *v0) && (tmp != *v1)) ? tmp : e2v_ids[_2e1+1];
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the time elapsed to build the cs_cdo_connect_t structure
 *
 * \return the value of the time elapsed in ns
 */
/*----------------------------------------------------------------------------*/

long long
cs_cdo_connect_get_time_perfo(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a \ref cs_range_set_t structure and a
 *        \ref cs_interface_set_t structure for schemes with DoFs at vertices
 *
 * \param[in]       mesh          pointer to a cs_mesh_t structure
 * \param[in]       n_vtx_dofs    number of DoFs per vertex
 * \param[in]       interlaced    false means a block viewpoint
 * \param[in, out]  p_ifs         pointer of  pointer to a cs_interface_set_t
 * \param[in, out]  p_rs          pointer of  pointer to a cs_range_set_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_assign_vtx_ifs_rs(const cs_mesh_t       *mesh,
                                 int                    n_vtx_dofs,
                                 bool                   interlaced,
                                 cs_interface_set_t   **p_ifs,
                                 cs_range_set_t       **p_rs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a \ref cs_range_set_t structure and a
 *        \ref cs_interface_set_t structure for schemes with DoFs at faces.
 *
 * \param[in]       mesh          pointer to a cs_mesh_t structure
 * \param[in]       n_faces       number of faces (interior + border)
 * \param[in]       n_face_dofs   number of DoFs per face
 * \param[in, out]  p_ifs         pointer of  pointer to a cs_interface_set_t
 * \param[in, out]  p_rs          pointer of  pointer to a cs_range_set_t
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_assign_face_ifs_rs(const cs_mesh_t       *mesh,
                                  cs_lnum_t              n_faces,
                                  int                    n_face_dofs,
                                  cs_interface_set_t   **p_ifs,
                                  cs_range_set_t       **p_rs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and define a new cs_interface_set_t structure on faces
 *
 * \param[in]  mesh          pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_interface_set_t structure
 */
/*----------------------------------------------------------------------------*/

cs_interface_set_t *
cs_cdo_connect_define_face_interface(const cs_mesh_t       *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a new cs_cdo_connect_t structure
 *        Range sets and interface sets are allocated and defined according to
 *        the value of the different scheme flags.
 *        cs_range_set_t structure related to vertices is shared the cs_mesh_t
 *        structure (the global one)
 *
 * \param[in, out]  mesh              pointer to a cs_mesh_t structure
 * \param[in]       eb_scheme_flag    metadata for Edge-based schemes
 * \param[in]       fb_scheme_flag    metadata for Face-based schemes
 * \param[in]       vb_scheme_flag    metadata for Vertex-based schemes
 * \param[in]       vcb_scheme_flag   metadata for Vertex+Cell-based schemes
 * \param[in]       hho_scheme_flag   metadata for HHO schemes
 *
 * \return  a pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_init(cs_mesh_t      *mesh,
                    cs_flag_t       eb_scheme_flag,
                    cs_flag_t       fb_scheme_flag,
                    cs_flag_t       vb_scheme_flag,
                    cs_flag_t       vcb_scheme_flag,
                    cs_flag_t       hho_scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a cs_cdo_connect_t structure
 *
 * \param[in]  connect     pointer to the cs_cdo_connect_t struct. to destroy
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_free(cs_cdo_connect_t   *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the discrete curl operator across each primal faces.
 *        From an edge-based array (seen as circulations) compute a face-based
 *        array (seen as fluxes)
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]      edge_values  array of values at edges
 * \param[in, out] curl_values  array storing the curl across faces (allocated
 *                              if necessary)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_discrete_curl(const cs_cdo_connect_t    *connect,
                             const cs_real_t           *edge_values,
                             cs_real_t                **p_curl_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the discrete curl operator across each primal faces.
 *        From an edge-based array (seen as circulations) compute a face-based
 *        array (seen as fluxes)
 *
 * \param[in]      connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]      edge_values  array of values at edges
 * \param[in, out] curl_values  array storing the curl across faces (allocated
 *                              if necessary)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_discrete_curl(const cs_cdo_connect_t    *connect,
                             const cs_real_t           *edge_values,
                             cs_real_t                **p_curl_values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of connectivity information
 *
 * \param[in]  connect           pointer to cs_cdo_connect_t structure
 * \param[in]  eb_scheme_flag    metadata for Edge-based schemes
 * \param[in]  vb_scheme_flag    metadata for Vertex-based schemes
 * \param[in]  vcb_scheme_flag   metadata for Vertex+Cell-based schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_summary(const cs_cdo_connect_t  *connect,
                       cs_flag_t                eb_scheme_flag,
                       cs_flag_t                vb_scheme_flag,
                       cs_flag_t                vcb_scheme_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_cdo_connect_t structure
 *
 * \param[in]  connect     pointer to cs_cdo_connect_t structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_dump(const cs_cdo_connect_t  *connect);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_CONNECT_H__ */
