#ifndef __CS_EVALUATE_H__
#define __CS_EVALUATE_H__

/*============================================================================
 * Functions and structures to deal with source term computation
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
#include "cs_time_step.h"

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
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      dof_flag   indicate where the evaluation has to be done
 * \param[in]      loc_id     id related to a cs_mesh_location_t struct.
 * \param[in]      def_type   type of definition
 * \param[in]      quad_type  type of quadrature (not always used)
 * \param[in]      use_subdiv consider or not the subdivision into tetrahedra
 * \param[in]      def        access to the definition of the values
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate(const cs_cdo_quantities_t    *quant,
            const cs_cdo_connect_t       *connect,
            const cs_time_step_t         *time_step,
            cs_flag_t                     dof_flag,
            int                           loc_id,
            cs_param_def_type_t           def_type,
            cs_quadra_type_t              quad_type,
            bool                          use_subdiv,
            cs_def_t                      def,
            double                       *p_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the 3x3 matrix related to a general material property.
 *         This value is computed at location (x,y,z) and time tcur.
 *
 * \param[in]     pty_id    id related to the material property to deal with
 * \param[in]     tcur      time at which we evaluate the material property
 * \param[in]     xyz       location at which  we evaluate the material property
 * \param[in]     invers    true or false
 * \param[in,out] matval    pointer to the 3x3 matrix to return
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_pty(int                pty_id,
                cs_real_t          tcur,
                const cs_real_3_t  xyz,
                bool               invers,
                cs_real_33_t      *matval);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the vector related to an advection field.
 *         This value is computed at location (x,y,z) and time tcur.
 *
 * \param[in]     adv_id    id related to the advection field to deal with
 * \param[in]     tcur      time at which we evaluate the material property
 * \param[in]     xyz       location at which  we evaluate the material property
 * \param[in,out] matval    pointer to the 3x3 matrix to return
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_adv_field(int                 adv_id,
                      cs_real_t           tcur,
                      const cs_real_3_t   xyz,
                      cs_real_3_t        *vect);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EVALUATE_H__ */
