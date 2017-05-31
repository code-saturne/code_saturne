#ifndef __CS_CDO_CONNECT_H__
#define __CS_CDO_CONNECT_H__

/*============================================================================
 * Manage connectivity (Topological features of the mesh)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"
#include "cs_mesh.h"
#include "cs_range_set.h"
#include "cs_sla.h"

#include "fvm_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Connectivity structure */
typedef struct {

  /* Vertex-related members */
  cs_lnum_t          n_vertices;
  cs_sla_matrix_t   *v2e;         // vertex --> edges connectivity

  /* Edge-related members */
  cs_lnum_t          n_edges;
  cs_sla_matrix_t   *e2f;         // edge --> faces connectivity
  cs_sla_matrix_t   *e2v;         // edge --> vertices connectivity

  /* Face-related members */
  cs_lnum_t          n_faces[3];  // 0: all, 1: border, 2: interior
  cs_sla_matrix_t   *f2c;         // face --> cells connectivity
  cs_sla_matrix_t   *f2e;         // face --> edges connectivity

  /* Cell-related members */
  cs_lnum_t            n_cells;
  fvm_element_t       *cell_type; // type of cell
  cs_flag_t           *cell_flag; // Flag (Border)
  cs_sla_matrix_t     *c2f;       // cell --> faces connectivity
  cs_connect_index_t  *c2e;       // cell -> edges connectivity
  cs_connect_index_t  *c2v;       // cell -> vertices connectivity

  /* Delta of ids between the min./max. values of entities related to a cell
     Useful to store compactly the link between mesh ids and cell mesh ids
     needed during the cell mesh definition */
  cs_lnum_t  e_max_cell_range;
  cs_lnum_t  v_max_cell_range;

  /* Max. connectitivy size for cells */
  int  n_max_vbyc;   // max. number of vertices in a cell
  int  n_max_ebyc;   // max. number of edges in a cell
  int  n_max_fbyc;   // max. number of faces in a cell
  int  n_max_vbyf;   // max. number of vertices in a face
  int  n_max_v2fc;   // max. number of faces connected to a vertex in a cell
  int  n_max_v2ec;   // max. number of edges connected to a vertex in a cell

  /* Structures to handle parallelism */
  const cs_range_set_t   *v_rs;  // range set related to vertices
  cs_range_set_t         *f_rs;  // range set related to faces

} cs_cdo_connect_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define a new cs_cdo_connect_t structure
 *        Range sets related to vertices and faces are computed inside and
 *        set as members of the cs_mesh_t structure
 *
 * \param[in, out]  mesh          pointer to a cs_mesh_t structure
 * \param[in]       scheme_flag   flag storing requested space schemes
 *
 * \return  a pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_init(cs_mesh_t      *mesh,
                    cs_flag_t       scheme_flag);

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
 * \brief  Summary of connectivity information
 *
 * \param[in]  connect     pointer to cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_summary(const cs_cdo_connect_t  *connect);

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
