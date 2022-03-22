#ifndef __CS_EQUATION_SYSTEM_PARAM_H__
#define __CS_EQUATION_SYSTEM_PARAM_H__

/*============================================================================
 * Functions to handle a set of coupled equations hinging on the cs_equation_t
 * structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

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

#include "cs_equation_param.h"
#include "cs_iter_algo.h"
#include "cs_param_types.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \enum cs_equation_system_sles_t
 *  \brief High-level information about the way of solving the system of
 *         equations
 *
 * \var CS_EQUATION_SYSTEM_SLES_MUMPS
 *      Associated keyword: "mumps"
 *      Direct solver to solve the full system
 *
 */

typedef enum {

  CS_EQUATION_SYSTEM_SLES_MUMPS,

  CS_EQUATION_SYSTEM_N_SLES_TYPES

} cs_equation_system_sles_strategy_t;

/*! \struct cs_equation_system_param_t
 *  \brief Main structure storing the parameter settings
 */

typedef struct {

  /*!
   * @name Generic metadata
   * @{
   *
   * \var name
   *      Name of the system of equations
   *
   * \var verbosity
   *      Level of verbosity (quantity of output information)
   *
   * \var space_scheme
   *      Associated space discretization. One assumes that all blocks share
   *      the same space discretization.
   *
   * \var block_var_dim
   *      Dimension of the variable in each block
   */

  char *restrict            name;

  int                       verbosity;

  cs_param_space_scheme_t   space_scheme;

  int                       block_var_dim;

  /*!
   * @}
   * @name Linear algebra (SLES)
   * @{
   *
   * \var sles_setup_done
   *      Indicator to know if the setup stage has be done
   *
   * \var sles_strategy
   *      Type of strategy used to solve the resulting system
   *
   * \var linear_solver
   *      Set of parameters to control the (main) linear solver
   */

  bool                                   sles_setup_done;

  cs_equation_system_sles_strategy_t     sles_strategy;

  cs_iter_algo_param_t                   linear_solver;

  /*!
   * @}
   */

} cs_equation_system_param_t;


/*! \enum cs_equation_system_key_t
 *  \brief List of available keys for setting the parameters of a system
 *         of equations
 *
 * \var CS_SYSKEY_LINEAR_SOLVER_ATOL
 *      Absolute tolerance for which the (main) linear solver stops iterating
 *
 * \var CS_SYSKEY_LINEAR_SOLVER_DTOL
 *      Divergence tolerance for which the (main) linear solver stops iterating
 *
 * \var CS_SYSKEY_LINEAR_SOLVER_RTOL
 *      Relative tolerance for which the (main) linear solver stops iterating
 *
 * \var CS_SYSKEY_LINEAR_SOLVER_VERBOSITY
 *      Level of verbosity for the (main) linear solver
 *
 * \var CS_SYSKEY_LINEAR_SOLVER_MAX_ITER
 *      Maximal number of iterations for the (main) linear solver
 *
 * \var CS_SYSKEY_SLES_STRATEGY
 *      Strategy for solving the linear system arising from the discretization
 *      of the system of equations
 *
 * \var CS_SYSKEY_VERBOSITY
 *      Level of verbosity related to the system of equations (a specific
 *      verbosity is available for the linear algebra;
 *      cf. CS_SYSKEY_LINEAR_SOLVER_VERBOSITY)
 */

typedef enum {

  CS_SYSKEY_LINEAR_SOLVER_ATOL,
  CS_SYSKEY_LINEAR_SOLVER_DTOL,
  CS_SYSKEY_LINEAR_SOLVER_RTOL,
  CS_SYSKEY_LINEAR_SOLVER_VERBOSITY,
  CS_SYSKEY_LINEAR_SOLVER_MAX_ITER,
  CS_SYSKEY_SLES_STRATEGY,
  CS_SYSKEY_VERBOSITY,

  CS_SYSKEY_N_KEYS

} cs_equation_system_key_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new cs_equation_system_param_t structure
 *
 * \param[in]  name            name of system of equations
 * \param[in]  block_var_dim   dimension of the variable in each block
 *
 * \return a pointer to a newly initialized cs_equation_system_param_t
 */
/*----------------------------------------------------------------------------*/

cs_equation_system_param_t *
cs_equation_system_param_create(const char       *name,
                                int               block_var_dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_system_param_t structure
 *
 * \param[in, out]  sysp     pointer to the structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_system_param_t *
cs_equation_system_param_free(cs_equation_system_param_t    *sysp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log the setup gathered in the structure cs_equation_system_param_t
 *
 * \param[in] sysp     pointer to a parameter structure to log
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_param_log(const cs_equation_system_param_t    *sysp);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter related to a keyname in a cs_equation_system_param_t
 *         structure
 *
 * \param[in, out] sysp     pointer to a parameter structure to set
 * \param[in]      key      key related to the member of eq to set
 * \param[in]      keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_param_set(cs_equation_system_param_t    *sysp,
                             cs_equation_system_key_t       key,
                             const char                    *keyval);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_SYSTEM_PARAM_H__ */
