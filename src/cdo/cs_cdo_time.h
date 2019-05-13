#ifndef __CS_CDO_TIME_H__
#define __CS_CDO_TIME_H__

/*============================================================================
 * Routines to handle common features related to the time scheme when using
 * CDO schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include "cs_cdo_local.h"
#include "cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the time discretization to a local system
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_cdo_time_scheme_t)(const cs_equation_param_t  *eqp,
                       const double                tpty_val,
                       const cs_sdm_t             *mass_mat,
                       const cs_flag_t             system_flag,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys);

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_structure_t according
 *         to the space scheme
 *
 * \param[in]  sys_flag       metadata about how is set the algebraic system
 * \param[in]  eqp            pointer to a cs_equation_param_t
 *
 * \return  a pointer to the function handling the time discretization
 */
/*----------------------------------------------------------------------------*/

cs_cdo_time_scheme_t *
cs_cdo_time_get_scheme_function(const cs_flag_t             sys_flag,
                                const cs_equation_param_t  *eqp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the RHS with the previously computed array of values
 *         Do not use OpenMP inside this function since it may be called from
 *         an OpenMP block
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t
 * \param[in]      stride      number of entries for each DoF
 * \param[in]      n_dofs      number of DoF to deal with
 * \param[in]      dof_ids     list of DoF ids or NULL if no indirection
 * \param[in]      values      array of values
 * \param[in, out] rhs         right-hand side to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_update_rhs(const cs_equation_param_t   *eqp,
                       int                          stride,
                       cs_lnum_t                    n_dofs,
                       const cs_lnum_t             *dof_ids,
                       const cs_real_t             *values,
                       cs_real_t                   *rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an implicit time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_imp(const cs_equation_param_t  *eqp,
                     const double                tpty_val,
                     const cs_sdm_t             *mass_mat,
                     const cs_flag_t             system_flag,
                     cs_cell_builder_t          *cb,
                     cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an implicit time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_imp(const cs_equation_param_t  *eqp,
                const double                tpty_val,
                const cs_sdm_t             *mass_mat,
                const cs_flag_t             system_flag,
                cs_cell_builder_t          *cb,
                cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an explicit time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_exp(const cs_equation_param_t  *eqp,
                     const double                tpty_val,
                     const cs_sdm_t             *mass_mat,
                     const cs_flag_t             system_flag,
                     cs_cell_builder_t          *cb,
                     cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system an explicit time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_exp(const cs_equation_param_t  *eqp,
                const double                tpty_val,
                const cs_sdm_t             *mass_mat,
                const cs_flag_t             system_flag,
                cs_cell_builder_t          *cb,
                cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system a "theta" time discretization when
 *          a CDO scheme is used and the mass matrix related to the time
 *          discretization is diagonal (lumping or Voronoi Hodge)
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_diag_theta(const cs_equation_param_t  *eqp,
                       const double                tpty_val,
                       const cs_sdm_t             *mass_mat,
                       const cs_flag_t             system_flag,
                       cs_cell_builder_t          *cb,
                       cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply to the local system a "theta" time discretization when
 *          a CDO scheme is used
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      tpty_val     current value of the time property
 * \param[in]      system_flag  indicate what is needed to build the system
 * \param[in]      mass_mat     pointer to a discrete Hodge op.
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 * \param[in, out] csys         pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_time_theta(const cs_equation_param_t  *eqp,
                  const double                tpty_val,
                  const cs_sdm_t             *mass_mat,
                  const cs_flag_t             system_flag,
                  cs_cell_builder_t          *cb,
                  cs_cell_sys_t              *csys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_TIME_H__ */
