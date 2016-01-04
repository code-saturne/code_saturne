#ifndef __CS_HODGE_H__
#define __CS_HODGE_H__

/*============================================================================
 * Build discrete Hodge operators
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
#include "cs_time_step.h"
#include "cs_cdo_toolbox.h"
#include "cs_sla.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_param.h"
#include "cs_property.h"

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
 * \brief   Allocate and initialize a cs_hodge_builder_t structure
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  h_info        algorithm used to build the discrete Hodge op.
 *
 * \return  a new allocated cs_hodge_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_init(const cs_cdo_connect_t   *connect,
                      cs_param_hodge_t          h_info);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_hodge_builder_t structure
 *
 * \param[in]  hb    pointer to the cs_hodge_builder_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_hodge_builder_t *
cs_hodge_builder_free(cs_hodge_builder_t  *hb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the property attached to a hodge builder
 *
 * \param[in, out]  hb       pointer to a cs_hodge_builder_t structure
 * \param[in]       ptyval   value of the property
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_builder_set_val(cs_hodge_builder_t    *hb,
                         cs_real_t              ptyval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the value of the property attached to a hodge builder
 *
 * \param[in, out]  hb       pointer to a cs_hodge_builder_t structure
 * \param[in]       ptymat   values of the tensor related to a property
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_builder_set_tensor(cs_hodge_builder_t     *hb,
                            const cs_real_33_t      ptymat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge
 *
 * \param[in]      c_id       cell id
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] hb         pointer to a cs_hodge_builder_t struct.
 *
 * \return a pointer to a cs_locmat_t struct. (local dense matrix)
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_hodge_build_local(int                         c_id,
                     const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     cs_hodge_builder_t         *hb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a discrete Hodge operator
 *
 * \param[in]  connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]  pty        pointer to a cs_property_t struct.
 * \param[in]  h_info     pointer to a cs_param_hodge_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

cs_sla_matrix_t *
cs_hodge_compute(const cs_cdo_connect_t      *connect,
                 const cs_cdo_quantities_t   *quant,
                 const cs_property_t         *pty,
                 const cs_param_hodge_t       h_info);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_HODGE_H__ */

