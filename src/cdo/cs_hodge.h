#ifndef __CS_HODGE_H__
#define __CS_HODGE_H__

/*============================================================================
 * Build discrete Hodge operators
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
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

#include "cs_cdo_toolbox.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_sla.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro and type definitions
 *============================================================================*/

typedef struct _hodge_builder_t cs_hodge_builder_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add or associate a discrete Hodge operator to a given setting
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 * \param[in]  quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]  param     set of parameters related a discrete Hodge operator
 *
 * \return  an id to the cs_hodge_t structure related to the given setting
 */
/*----------------------------------------------------------------------------*/

int
cs_hodge_add(const cs_cdo_connect_t     *connect,
             const cs_cdo_quantities_t  *quant,
             cs_param_hodge_t            param);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize by default the matrix related to a discrete
 *          Hodge operator
 *
 * \param[in]    connect  pointer to a cs_cdo_connect_t structure
 * \param[in]    quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]    type     type of the discrete Hodge op. to initiliaze
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_init_matrix(const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     const cs_param_hodge_type_t  type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_hodge_builder_t structure
 *
 * \param[in]  n_max_ent    max number of entities by primal cell
 *
 * \return  a new allocated cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_init(int   n_max_ent);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_hodge_builder_t structure
 *
 * \param[in]  hcq    pointer to the cs_hodge_builder_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_free(cs_hodge_builder_t  *hcq);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic COST algo.
 *
 * \param[in]     cid        cell id
 * \param[in]     connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]     quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]     h_info     pointer to a cs_param_hodge_t struct.
 * \param[in,out] hl         pointer to a cs_toolbox_locmat_t struct.
 * \param[in,out] hb         pointer to a cs_hodge_builder_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_cost_build_local(int                         cid,
                          const cs_cdo_connect_t     *connect,
                          const cs_cdo_quantities_t  *quant,
                          const cs_param_hodge_t      h_info,
                          cs_toolbox_locmat_t        *hl,
                          cs_hodge_builder_t         *hb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a discrete Hodge operator using a generic reconstruction
 *          algorithm: Reconstruction with Consistent and Coercive part
 *          stabilized by Beta (COST)
 *          A call to cs_hodge_builder_init() should have been done before.
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]  h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_cost_build(const cs_cdo_connect_t      *connect,
                    const cs_cdo_quantities_t   *quant,
                    const cs_param_hodge_t       h_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   H is an operator from primal edges to dual faces. It mimics a Hodge
 *          operator from primal mesh to dual mesh.
 *          Use voronoi algorithm.
 *
 * \param[in]   connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]   quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]    h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_voronoi_build(const cs_cdo_connect_t      *connect,
                       const cs_cdo_quantities_t   *quant,
                       const cs_param_hodge_t       h_info);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HODGE_H__ */

