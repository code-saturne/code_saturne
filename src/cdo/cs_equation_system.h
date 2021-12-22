#ifndef __CS_EQUATION_SYSTEM_H__
#define __CS_EQUATION_SYSTEM_H__

/*============================================================================
 * Functions to handle a set of coupled equations hinging on the cs_equation_t
 * structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_equation_system_t
 *  \brief Main structure to handle a set of coupled equations
 */

typedef struct {

  /*!
   * @name Metadata
   * @{
   *
   * \var name
   * name of the system of equations
   *
   * \var n_equations
   * Number of coupled equations (> 1) composing the system
   *
   * \var equations
   * Array of pointer to the equations constituting the coupled system. These
   * equations correspond to the each row and the cs_equation_param_t
   * associated to an equation corresponds to the setting of the diagonal
   * block.
   */

  char *restrict    name;
  int               n_equations;
  cs_equation_t   **equations;

  int               timer_id;   /*!< Id of the timer statistics */

  /*!
   * @}
   * @name Cross-terms
   * @{
   *
   * The setting of each block relies on the cs_equation_param_t structure.
   * The cs_equation_param_t structures related to the diagonal blocks are
   * shared with the cs_equation_t structures in the "equations" member and
   * thus not owned by the current structure.  The extra-diagonal blocks
   * dealing with the cross terms (i.e. the coupling between variables) are
   * owned by this structure.
   *
   * By default, there is no cross-term (i.e. params[1] = NULL)
   *
   * \varam params
   * Matrix of size n_coupled_equations (stored as an array of size
   * n_coupled_equations^2)
   */

  cs_equation_param_t   **params;

  /*!
   * @}
   */

} cs_equation_system_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new structure to handle system of coupled equations
 *
 * \param[in] sysname     name of the system of equations
 * \param[in] n_eqs       number of coupled equations composing the system
 *
 * \return  a pointer to the new allocated cs_equation_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_system_t *
cs_equation_system_create(const char   *sysname,
                          int           n_eqs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a structure used to couple equations
 *
 * \param[in, out] p_eqsys    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_free(cs_equation_system_t  **p_eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log a structure used to couple equations
 *
 * \param[in] eqsys    pointer to the structure to log
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log(cs_equation_system_t  *eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve of a system of coupled equations. Unsteady case.
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in, out] eqsys      pointer to the structure to log
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_solve(bool                     cur2prev,
                         const cs_mesh_t         *mesh,
                         cs_equation_system_t    *eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the given equation and associate it to the row block with id
 *         equal to row_id. The equation parameter is also set.
 *
 * \param[in]      row_id  position in the block rows
 * \param[in]      eq      pointer to the equation to add
 * \param[in, out] eqsys   pointer to a cs_equation_system_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_equation(int                       row_id,
                                cs_equation_t            *eq,
                                cs_equation_system_t     *eqsys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the given equation parameters and associate it to the matrix of
 *         equation parameters at (row_id, col_id)
 *
 * \param[in]      row_id   row position id
 * \param[in]      col_id   column position id
 * \param[in]      eqp      pointer to the equation parameter to add
 * \param[in, out] eqsys    pointer to a cs_equation_system_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_param(int                       row_id,
                             int                       col_id,
                             cs_equation_param_t      *eqp,
                             cs_equation_system_t     *eqsys);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_SYSTEM_H__ */
