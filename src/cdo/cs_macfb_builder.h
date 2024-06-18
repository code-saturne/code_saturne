#ifndef __CS_MACFB_BUILDER_H__
#define __CS_MACFB_BUILDER_H__

/*============================================================================
 * Low-level functions and structures used to build the algebraic system with
 * a cellwise process when MAC-fb schemes are set for the space
 * discretization
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_property.h"
#include "cs_sdm.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Cellwise builder for MAC-fb discretization */
typedef struct {

  /* Cell informations */

  cs_lnum_t c_id; /*!< cell id  */

  short int n_dofs; /*!< number of dofs in a cell */

  short int n_max_dofs; /*!< maximum number of dofs in a cell */

  cs_lnum_t dof_ids[30]; /*!< face DoF ids */

  /* Face informations */

  short int n_fc; /*!< number of faces in a cell (inner and outer) > 6 */

  cs_lnum_t f_ids[30]; /*!< face ids on this rank */

  cs_flag_cartesian_axis_t f_axis[30]; /*!< axis of local faces */
  short int f_sgn_axis[6]; /*!< axis incidence number between f and c */

  cs_real_t f_vol_cv[6];  /*!< volume of the control volume for a face */
  cs_real_t f_h_cv[6];    /*!< normal length of the control volume for a face */
  short int f_opp_idx[6]; /*!< opposite face index to the
                               face inside the cell */

  short int f2f_idx[24]; /*!< cellwise face->face index with same
                            direction and sharing an edge with the current face
                            Negative index means to skip entry
                            4 entries by face */

  cs_lnum_t f2f_ids[24]; /*!< cellwise face->face ids with same
                             direction and sharing an edge with the current face
                             Negative index means to that is not a face
                             4 entries by face  */

  short int f2fo_idx[48]; /*!< cellwise face->face index with orthornal
                          direction and sharing an edge with the current face
                          Negative index means to skip entry
                          8 = (4*2) entries by face */

  cs_real_t f2f_h[24]; /*!< distance beetwen barycenter of each face->face
                            4 entries by face */
  cs_real_t f2f_surf_cv_c[24]; /*!< surface of each face of the controle
                                volume face->face (restricted to the cell) */

  cs_lnum_t f2e_ids[24]; /*!< cellwise face->edges ids, same order as f2f */

  cs_real_t dir_values[24]; /*!< dirichlet values at edges */

} cs_macfb_builder_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Pointer of pointers to global structures */

extern cs_macfb_builder_t **cs_mac_builders;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate global structures used for MAC builder
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to \ref cs_macfb_builder_t
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_macfb_builder_t structure corresponding to id
 *
 * \param[in]   id   id in the array of pointer to cs_macfb_builder_t struct.
 *
 * \return a pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_builder_t *cs_macfb_get_builder(int id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate  a cs_macfb_builder_t structure
 *
 * \return a pointer to a new allocated cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_macfb_builder_t *cs_macfb_builder_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_macfb_builder_t structure
 *
 * \param[in, out] p_builder  pointer of pointer on a cs_macfb_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_free(cs_macfb_builder_t **p_builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize to invalid values a cs_macfb_builder_t structure
 *
 * \param[in, out]  macb       pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_reset(cs_macfb_builder_t *macb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set-up face informations needed to build operators
 *
 * \param[in]       cm         pointer to a cs_cell_mesh_t structure
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       quant      pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  macb       pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_cellwise_setup(const cs_cell_mesh_t      *cm,
                                     const cs_cdo_connect_t    *connect,
                                     const cs_cdo_quantities_t *quant,
                                     cs_macfb_builder_t        *macb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_macfb_builder_t structure
 *
 * \param[in]    macb    pointer to a cs_macfb_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void cs_macfb_builder_dump(const cs_macfb_builder_t *macb);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MACFB_BUILDER_H__ */
