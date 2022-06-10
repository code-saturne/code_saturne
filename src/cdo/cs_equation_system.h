#ifndef __CS_EQUATION_SYSTEM_H__
#define __CS_EQUATION_SYSTEM_H__

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

#include "cs_equation.h"
#include "cs_equation_priv.h"
#include "cs_equation_system_param.h"
#include "cs_param_types.h"

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
 * \brief Create and initialize equation builders and scheme context for each
 *        equation which are in the extra-diagonal blocks related to a system
 *        of equations. Structures associated to diagonal blocks should be
 *        already initialized during the treatment of the classical equations.
 *
 *        Generic prototype to define the function pointer.
 *
 * \param[in]      n_eqs       number of equations
 * \param[in]      sysp        set of parameters to specify a system of eqs
 * \param[in, out] core_array  array of the core structures for an equation
 * \param[out]     p_sh        double pointer to a system helper to define
 *
 * \return a pointer to a new allocated system context structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_system_define_t)(int                                 n_eqs,
                              const cs_equation_system_param_t   *sysp,
                              cs_equation_core_t                **core_array,
                              cs_cdo_system_helper_t            **p_sh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free an array of structures (equation parameters, equation builders
 *        or scheme context) for each equation which are in the extra-diagonal
 *        blocks related to a system of equations. Structures associated to
 *        diagonal blocks are freed during the treatment of the classical full
 *        equations.
 *
 *        Generic prototype to define the function pointer.
 *
 * \param[in]      n_eqs        number of equations
 * \param[in, out] core_array   array of the core structures for an equation
 * \param[in, out] sys_context  pointer to a context structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_system_free_t)(int                      n_eqs,
                            cs_equation_core_t     **core_array,
                            void                    *sys_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve a linear system within the CDO framework
 *
 * \param[in]      c2p         do a "current to previous" operation performed ?
 * \param[in]      n_eqs       number of equations in the system to solve
 * \param[in]      sysp        set of paremeters for the system of equations
 * \param[in, out] blocks      array of the core members for an equation
 * \param[in, out] sys_context pointer to a sructure cast on-the-fly ?
 * \param[in, out] sh          pointer to a system helper structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_system_solve_t)(bool                           c2p,
                             int                            n_eqs,
                             cs_equation_system_param_t    *sysp,
                             cs_equation_core_t           **blocks,
                             void                          *sys_context,
                             cs_cdo_system_helper_t        *sh);


/*! \struct cs_equation_system_t
 *  \brief Main structure to handle a set of coupled equations
 */

typedef struct {

  /*!
   * @name Metadata
   * @{
   *
   * \var param
   *      Set of parameters to specify the settings of the system of equations
   *
   */

  cs_equation_system_param_t   *param;

  cs_timer_counter_t            timer;     /*!< Cumulated elapsed time */
  int                           timer_id;  /*!< Id of the timer statistics */

  /*!
   * @name Structure of the system of equations
   * @{
   *
   * \var system_helper
   *      Set of structures to define the system of equations (rhs, matrix or
   *      matrices, range sets, interface sets, etc.)
   *
   * \var context
   *      pointer to a structure cast on-the-fly which depends on the numerical
   *      scheme
   */

  cs_cdo_system_helper_t       *system_helper;

  void                         *context;

  /*!
   * @name Diagonal block (equations)
   * @{
   *
   * \var n_equations
   *      Number of coupled equations (> 1) composing the system
   *
   * \var equations
   *      Array of pointer to the equations constituting the coupled
   *      system. These equations correspond to the each row and the
   *      cs_equation_param_t associated to an equation corresponds to the
   *      setting of the diagonal block.
   */

  int                           n_equations;
  cs_equation_t               **equations;

  /*!
   * @}
   * @name Crossed terms
   * @{
   *
   * The setting of each block relies on the cs_equation_param_t structure.
   * The cs_equation_param_t structures related to the diagonal blocks are
   * shared with the cs_equation_t structures in the "equations" member and
   * thus not owned by the current structure. The extra-diagonal blocks
   * dealing with the crossed terms (i.e. the coupling between variables) are
   * owned by this structure.
   *
   * By default, there is no crossed term (i.e. params[1] = NULL)
   *
   * The same rationale applies to builder structures and scheme context
   * structures. All these structures are contained in the structure \ref
   * cs_equation_core_t to avoid manipulating void ** structures
   *
   * \var block_factories
   *      Matrix of cs_equation_core_t structures. The size of the matrix is
   *      n_equations (stored as an array of size n_equations^2). These
   *      structures enable to build and solve the system of equations.
   */

  cs_equation_core_t          **block_factories;

  /*!
   * @}
   * @name Pointer to functions
   * @{
   *
   * \var define
   *      Initialize builder and scheme context structures. Pointer of function
   *      given by the prototype cs_equation_system_init_structures_t
   *
   * \var free
   *      Free builder and scheme context structures. Pointer of function given
   *      by the prototype cs_equation_system_free_context_t
   *
   * \var solve_system
   *      Solve the system of equations (unsteady case). Pointer of function
   *      given by the generic prototype cs_equation_system_solve_t
   *
   * \var solve_steady_state_system
   *      Solve the system of equations (steady-state case). Pointer of
   *      function given by the generic prototype cs_equation_system_solve_t
   */

  cs_equation_system_define_t     *define;
  cs_equation_system_free_t       *free;

  cs_equation_system_solve_t      *solve_system;
  cs_equation_system_solve_t      *solve_steady_state_system;

  /*!
   * @}
   */

} cs_equation_system_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of systems of equations
 *
 * \return the number of systems
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_system_get_n_systems(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new structure to handle system of coupled equations
 *
 * \param[in] sysname         name of the system of equations
 * \param[in] n_eqs           number of coupled equations composing the system
 * \param[in] block_var_dim   dimension of the variable in each block
 *
 * \return  a pointer to the new allocated cs_equation_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_system_t *
cs_equation_system_add(const char             *sysname,
                       int                     n_eqs,
                       int                     block_var_dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_equation_system_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_destroy_all(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup for all structures managing systems of equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log_monitoring(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of shared pointer to the main structures
 *
 * \param[in]  mesh        basic mesh structure
 * \param[in]  connect     additional connectivity data
 * \param[in]  quant       additional mesh quantities
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_init_sharing(const cs_mesh_t             *mesh,
                                const cs_cdo_connect_t      *connect,
                                const cs_cdo_quantities_t   *quant,
                                const cs_time_step_t        *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing all the systems of
 *         equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_functions(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the SLES associated to each system of equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_sles(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the builder and scheme context structures associated to all
 *         the systems of equations which have been added.
 *         For the diagonal blocks, one relies on the builder and context of
 *         the related equations. For extra-diagonal blocks, one defines new
 *         builder and context structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_define(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve of a system of coupled equations. Unsteady case.
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in, out] eqsys      pointer to the structure to solve
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_solve(bool                     cur2prev,
                         cs_equation_system_t    *eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign the given equation to the diagonal block located at
 *         position (row_id, row_id) in the matrix of blocks
 *
 * \param[in]      row_id  position in the block matrix
 * \param[in]      eq      pointer to the equation to add
 * \param[in, out] eqsys   pointer to a cs_equation_system_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_assign_equation(int                       row_id,
                                   cs_equation_t            *eq,
                                   cs_equation_system_t     *eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign the given equation parameters to the block with ids
 *         (row_id, col_id) in the block matrix
 *
 * \param[in]      row_id   row position id
 * \param[in]      col_id   column position id
 * \param[in]      eqp      pointer to the equation parameter to add
 * \param[in, out] eqsys    pointer to a cs_equation_system_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_assign_param(int                       row_id,
                                int                       col_id,
                                cs_equation_param_t      *eqp,
                                cs_equation_system_t     *eqsys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_SYSTEM_H__ */
