#ifndef __CS_SADDLE_SOLVER_SETUP_H__
#define __CS_SADDLE_SOLVER_SETUP_H__

/*============================================================================
 * Routines to handle the SLES settings for a saddle-point problem
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_saddle_solver.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_param_saddle_setup.h

  \brief Routines to handle the setup of SLES relying on a \ref
         cs_param_saddle_t structure
*/

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new solver for solving a saddle-point problem.
 *
 * \param[in] n1_elts         number of elements associated to the (1,1)-block
 * \param[in] n1_dofs_by_elt  number of DoFs by elements in the (1,1)-block
 * \param[in] n2_elts         number of elements associated to the (2,2)-block
 * \param[in] n2_dofs_by_elt  number of DoFs by elements in the (2,2)-block
 * \param[in] saddlep         set of parameters for the saddle-point solver
 * \param[in] sh              pointer to a system helper structure
 * \param[in] main_sles       pointer to the main SLES structure related to
 *                            this saddle-point problem
 *
 * \return a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_saddle_solver_t *
cs_saddle_solver_add(cs_lnum_t                 n1_elts,
                     int                       n1_dofs_by_elt,
                     cs_lnum_t                 n2_elts,
                     int                       n2_dofs_by_elt,
                     const cs_param_saddle_t  *saddlep,
                     cs_cdo_system_helper_t   *sh,
                     cs_sles_t                *main_sles);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all remaining structures related to saddle-point solvers
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the SLES structures (potentially several) related to a
 *        saddle-point system
 */
/*----------------------------------------------------------------------------*/

void
cs_saddle_solver_setup_sles(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SADDLE_SOLVER_SETUP_H__ */
