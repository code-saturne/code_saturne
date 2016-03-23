#ifndef __CS_EVALUATE_H__
#define __CS_EVALUATE_H__

/*============================================================================
 * Functions and structures to deal with evaluation of quantities
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
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect,
                                const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         The value defined by the analytic function is by unity of volume
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      ml_id       id related to a cs_mesh_location_t structure
 * \param[in]      ana         accessor to values thanks to a function pointer
 * \param[in]      quad_type   type of quadrature (not always used)
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_from_analytic(cs_flag_t              dof_flag,
                                  int                    ml_id,
                                  cs_analytic_func_t    *ana,
                                  cs_quadra_type_t       quad_type,
                                  double                 retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         Accessor to the value is by unit of volume
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value to set
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_from_value(cs_flag_t       dof_flag,
                               int             ml_id,
                               cs_get_t        get,
                               double          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a quantity defined by analytic
 *         function for all the degrees of freedom
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      ml_id       id related to a cs_mesh_location_t structure
 * \param[in]      ana         accessor to values thanks to a function pointer
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_analytic(cs_flag_t              dof_flag,
                                    int                    ml_id,
                                    cs_analytic_func_t    *ana,
                                    double                 retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF in the case of a potential field in order
 *         to put a given quantity inside the volume associated to ml_id
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value related to the
 *                           quantity to put in the volume spanned by ml_id
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_qov(cs_flag_t       dof_flag,
                               int             ml_id,
                               cs_get_t        get,
                               double          retval[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store the value related to each DoF in the case of a potential field
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value to set
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_value(cs_flag_t       dof_flag,
                                 int             ml_id,
                                 cs_get_t        get,
                                 double          retval[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EVALUATE_H__ */
