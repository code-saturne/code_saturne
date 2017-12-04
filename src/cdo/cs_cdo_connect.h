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

#include "fvm_defs.h"

#include "cs_base.h"
#include "cs_flag.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_range_set.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Main categories to consider for high-level structures
   Remark: HHO-P1 and CDO-Fb vector-valued shares the same structures
*/
#define CS_CDO_CONNECT_VTX_SCA   0 /* Vb or VCb scalar-valued eq. */
#define CS_CDO_CONNECT_FACE_SP0  1 /* Fb or HHO-P0 scalar-valued eq. */
#define CS_CDO_CONNECT_FACE_SP1  2 /* HHO-P1 scalar-valued */
#define CS_CDO_CONNECT_FACE_VP0  2 /* Fb vector-valued eq. */
#define CS_CDO_CONNECT_FACE_SP2  3 /* HHO-P2 scalar-valued eq. */
#define CS_CDO_CONNECT_N_CASES   4

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Connectivity structure */
typedef struct {

  cs_lnum_t          n_vertices;
  cs_lnum_t          n_edges;
  cs_lnum_t          n_faces[3];  // 0: all, 1: border, 2: interior
  cs_lnum_t          n_cells;

  /* Edge-related members */
  cs_adjacency_t    *e2v;         // edge --> vertices connectivity

  /* Face-related members */
  cs_adjacency_t    *f2c;         // face --> cells connectivity
  cs_adjacency_t    *f2e;         // face --> edges connectivity

  /* Cell-related members */
  fvm_element_t     *cell_type;   // type of cell
  cs_flag_t         *cell_flag;   // Flag (Border)
  cs_adjacency_t    *c2f;         // cell --> faces connectivity
  cs_adjacency_t    *c2e;         // cell --> edges connectivity
  cs_adjacency_t    *c2v;         // cell --> vertices connectivity

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

  /* Structures to handle parallelism/assembler */
  cs_range_set_t      *range_sets[CS_CDO_CONNECT_N_CASES];
  cs_interface_set_t   *interfaces[CS_CDO_CONNECT_N_CASES];

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
 *        Range sets and interface sets are allocated and defined according to
 *        the value of the different scheme flags.
 *        cs_range_set_t structure related to vertices is shared the cs_mesh_t
 *        structure (the global one)
 *
 * \param[in, out]  mesh             pointer to a cs_mesh_t structure
 * \param[in]       vb_scheme_flag   metadata for Vb schemes
 * \param[in]       vcb_scheme_flag  metadata for V+C schemes
 * \param[in]       fb_scheme_flag   metadata for Fb schemes
 * \param[in]       hho_scheme_flag  metadata for HHO schemes
 *
 * \return  a pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_init(cs_mesh_t      *mesh,
                    cs_flag_t       vb_scheme_flag,
                    cs_flag_t       vcb_scheme_flag,
                    cs_flag_t       fb_scheme_flag,
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
