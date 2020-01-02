#ifndef __CS_STATIC_CONDENSATION_H__
#define __CS_STATIC_CONDENSATION_H__

/*============================================================================
 * Routines to handle manipulations related to the static condensation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and store
 *          information inside the rc_tilda and acx_tilda to be able to compute
 *          the values at cell centers
 *          rc_tilda = Acc^-1 * cell_rhs
 *          Case of scalar-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in, out] rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in, out] acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_scalar_eq(const cs_adjacency_t    *c2x,
                                 cs_real_t               *rc_tilda,
                                 cs_real_t               *acx_tilda,
                                 cs_cell_builder_t       *cb,
                                 cs_cell_sys_t           *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Opposite process of the static condensation.
 *          Define the field at cells given the field at x locations and arrays
 *          storing the static condensation.
 *          Case of scalar-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in]      rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in]      acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in]      px          values of the fields at x locations
 * \param[in, out] pc          values of the field at cells
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_recover_scalar(const cs_adjacency_t    *c2x,
                                      const cs_real_t         *rc_tilda,
                                      const cs_real_t         *acx_tilda,
                                      const cs_real_t         *px,
                                      cs_real_t               *pc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and store
 *          information inside the rc_tilda and acx_tilda to be able to compute
 *          the values at cell centers
 *          rc_tilda = Acc^-1 * cell_rhs
 *          Case of vector-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in, out] rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in, out] acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in, out] cb          pointer to a cs_cell_builder_t structure
 * \param[in, out] csys        pointer to a cs_cell_sys_t structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_vector_eq(const cs_adjacency_t    *c2x,
                                 cs_real_t               *rc_tilda,
                                 cs_real_t               *acx_tilda,
                                 cs_cell_builder_t       *cb,
                                 cs_cell_sys_t           *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Opposite process of the static condensation.
 *          Define the field at cells given the field at x locations and arrays
 *          storing the static condensation.
 *          Case of vector-valued CDO equations
 *
 * \param[in]      c2x         pointer to a cs_adjacency_t structure
 * \param[in]      rc_tilda    pointer to the rhs related to cell DoFs (Acc-1
 * \param[in]      acx_tilda   pointer to an unrolled matrix Acc^-1 * Acx
 * \param[in]      px          values of the fields at x locations
 * \param[in, out] pc          values of the field at cells
 */
/*----------------------------------------------------------------------------*/

void
cs_static_condensation_recover_vector(const cs_adjacency_t    *c2x,
                                      const cs_real_t         *rc_tilda,
                                      const cs_real_t         *acx_tilda,
                                      const cs_real_t         *px,
                                      cs_real_t               *pc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_STATIC_CONDENSATION_H__ */
