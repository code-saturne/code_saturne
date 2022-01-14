#ifndef __CS_CDOVB_SCALSYS_H__
#define __CS_CDOVB_SCALSYS_H__

/*============================================================================
 * Build an algebraic CDO vertex-based system of equations. These equations
 * corresponds to scalar-valued unsteady convection diffusion reaction
 * equations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_cdo_quantities.h"
#include "cs_equation.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_equation_system_param.h"
#include "cs_mesh.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to the main shared structures
 *
 * \param[in]  mesh        basic mesh structure
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_init_common(const cs_mesh_t              *mesh,
                             const cs_cdo_connect_t       *connect,
                             const cs_cdo_quantities_t    *quant,
                             const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize equation builders and scheme context for each
 *        equation which are in the extra-diagonal blocks related to a system
 *        of equations. Structures associated to diagonal blocks should be
 *        already initialized during the treatment of the classical full
 *        equations.
 *
 *        Case of scalar-valued CDO-Vb scheme in each block
 *
 * \param[in]      n_eqs            number of equations
 * \param[in, out] block_factories  array of the core members for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_init_structures(int                        n_eqs,
                                 cs_equation_core_t       **block_factories);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free an array of structures (equation parameters, equation builders
 *        or scheme context) for each equation which are in the extra-diagonal
 *        blocks related to a system of equations. Structures associated to
 *        diagonal blocks are freed during the treatment of the classical full
 *        equations.
 *
 *        Case of scalar-valued CDO-Vb scheme in each block
 *
 * \param[in]      n_eqs    number of equations
 * \param[in, out] blocks   array of the core structures for an equation
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_free_structures(int                        n_eqs,
                                 cs_equation_core_t       **blocks);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build and solve the linear system of equations. The number of rows in
 *        the system is equal to the number of equations. Thus there are
 *        n_eqs*n_eqs blocks in the system. Each block corresponds potentially
 *        to a scalar-valued unsteady convection/diffusion/reaction equation
 *        with a CDO-Vb scheme using an implicit time scheme.
 *
 * \param[in]      cur2prev  true="current to previous" operation is performed
 * \param[in]      n_eqs     number of equations
 * \param[in]      sysp      set of paremeters for the system of equations
 * \param[in, out] blocks    array of the core members for an equation
 * \param[in, out] p_ms      double pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scalsys_solve_implicit(bool                           cur2prev,
                                int                            n_equations,
                                cs_equation_system_param_t    *sysp,
                                cs_equation_core_t           **blocks,
                                cs_matrix_structure_t        **p_ms);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDOVB_SCALSYS_H__ */
