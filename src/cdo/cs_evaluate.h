#ifndef __CS_EVALUATE_H__
#define __CS_EVALUATE_H__

/*============================================================================
 * Functions and structures to deal with source term computation
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

#include "cs_cdo.h"
#include "cs_cdo_quantities.h"
#include "cs_quadrature.h"
#include "cs_param.h"

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
 * \brief  Compute the contribution related to a source term or a
 *         boundary condition for instance
 *
 * \param[in]    m          pointer to a cs_mesh_t struct.
 * \param[in]    quant      additional mesh quantities struct.
 * \param[in]    connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]    tcur      current physical time of the simulation
 * \param[in]    dof_flag   indicate where the evaluation has to be done
 * \param[in]    loc_id     id related to a cs_mesh_location_t struct.
 * \param[in]    def_type   type of definition
 * \param[in]    quad_type  type of quadrature (not always used)
 * \param[in]    def        access to the definition of the values
 * \param[inout] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate(const cs_mesh_t              *m,
            const cs_cdo_quantities_t    *quant,
            const cs_cdo_connect_t       *connect,
            double                        tcur,
            cs_flag_t                     dof_flag,
            int                           loc_id,
            cs_param_def_type_t           def_type,
            cs_quadra_type_t              quad_type,
            cs_def_t                      def,
            double                       *p_values[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EVALUATE_H__ */
