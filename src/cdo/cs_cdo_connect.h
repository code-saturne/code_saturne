#ifndef __CS_CDO_CONNECT_H__
#define __CS_CDO_CONNECT_H__

/*============================================================================
 * Manage connectivity (Topological features of the mesh)
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

#include "cs_cdo.h"
#include "cs_cdo_toolbox.h"
#include "cs_sla.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* First level of information */
#define CS_CDO_CONNECT_IN    (1 << 0)  /* Interior entity */
#define CS_CDO_CONNECT_BD    (1 << 1)  /* Border entity */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  cs_flag_t  *flag;   // CS_CDO_CONNECT_IN/BD

  cs_lnum_t   n_elts;    // number of entities
  cs_lnum_t   n_i_elts;  // number of interior entities
  cs_lnum_t   n_b_elts;  // number of border entities

} cs_connect_info_t;

typedef struct { /* Connectivity structure */

  /* Upward oriented connectivity */
  cs_sla_matrix_t   *v2e;  // vertex --> edges connectivity
  cs_sla_matrix_t   *e2f;  // edge --> faces connectivity
  cs_sla_matrix_t   *f2c;  // face --> cells connectivity

  /* Downward oriented connectivity */
  cs_sla_matrix_t   *e2v;  // edge --> vertices connectivity
  cs_sla_matrix_t   *f2e;  // face --> edges connectivity
  cs_sla_matrix_t   *c2f;  // cell --> faces connectivity

  /* Specific CDO connect : not oriented (same spirit as Code_Saturne
     historical connectivity).
     Use this connectivity to scan dual quantities */
  cs_connect_index_t  *c2e;  // cell -> edges connectivity
  cs_connect_index_t  *c2v;  // cell -> vertices connectivity
  cs_connect_index_t  *v2v;  // vertex --> vertices through cell connectivity

  /* Max. connectitivy size for cells */
  cs_lnum_t  n_max_vbyc;    // max. number of vertices in a cell
  cs_lnum_t  n_max_ebyc;    // max. number of edges in a cell
  cs_lnum_t  n_max_fbyc;    // max. number of faces in a cell
  cs_lnum_t  n_max_vbyf;    // max. number of vertices in a face

  /* Status internal or border entity */
  cs_connect_info_t  *v_info;  // count of interior/border vertices
  cs_connect_info_t  *e_info;  // count of interior/border edges
  cs_connect_info_t  *f_info;  // count of interior/border faces
  cs_connect_info_t  *c_info;  /* count of interior/border cells
                                  a border cell has at least one border face */

} cs_cdo_connect_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  String related to flag in cs_connect_info_t
 *
 * \param[in]  flag     retrieve name for this flag
 */
/*----------------------------------------------------------------------------*/

const char *
cs_cdo_connect_flagname(short int  flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate and define by default a cs_cdo_connect_t structure
 *
 * \param[in]  m            pointer to a cs_mesh_t structure
 *
 * \return  a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_connect_t *
cs_cdo_connect_init(const cs_mesh_t      *m);

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
 * \brief Update a cs_cdo_connect_t structure with respect to the requested
 *        space scheme
 *
 * \param[in, out]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]       scheme_flag  type of numerical schemes requested
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_connect_update(cs_cdo_connect_t       *connect,
                      cs_flag_t               scheme_flag);

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
