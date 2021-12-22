/*============================================================================
 * Functions to handle the cs_equation_system_t structure
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_equation_param.h"
#include "cs_equation_priv.h"
#include "cs_timer_stats.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_EQUATION_SYSTEM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
                          int           n_eqs)
{
  cs_equation_system_t  *eqsys = NULL;

  if (n_eqs < 2)
    return NULL;

  BFT_MALLOC(eqsys, 1, cs_equation_system_t);

  /* Store sysname */

  size_t  len = strlen(sysname);
  BFT_MALLOC(eqsys->name, len + 1, char);
  strncpy(eqsys->name, sysname, len);
  eqsys->name[len] = '\0';

  /* Set timer statistic structure to a default value */

  eqsys->timer_id = cs_timer_stats_id_by_name(sysname);
  if (eqsys->timer_id < 0)
    eqsys->timer_id = cs_timer_stats_create(NULL, /* new root */
                                            sysname,
                                            sysname);

  eqsys->n_equations = n_eqs;

  BFT_MALLOC(eqsys->equations, n_eqs, cs_equation_t *);
  for (int i = 0; i < n_eqs; i++)
    eqsys->equations[i] = NULL;

  BFT_MALLOC(eqsys->params, n_eqs*n_eqs, cs_equation_param_t *);
  for (int i = 0; i < n_eqs*n_eqs; i++)
    eqsys->params[i] = NULL;

  return eqsys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a structure used to couple equations
 *
 * \param[in, out] p_eqsys    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_free(cs_equation_system_t  **p_eqsys)
{
  if (p_eqsys == NULL)
    return;

  cs_equation_system_t  *eqsys = *p_eqsys;

  if (eqsys == NULL)
    return;

  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  int  n_eqs = eqsys->n_equations;

  BFT_FREE(eqsys->name);
  BFT_FREE(eqsys->equations);

  /* Free the extra-diagonal cs_equation_param_t */

  for (int i = 0; i < n_eqs; i++) {

    cs_equation_param_t  **eqp_array = eqsys->params + i*n_eqs;

    for (int j = 0; j < n_eqs; j++)
      if (i != j)
        eqp_array[j] = cs_equation_param_free(eqp_array[j]);

  }

  BFT_FREE(eqsys->params);

  BFT_FREE(eqsys);
  *p_eqsys = NULL;

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log a structure used to couple equations
 *
 * \param[in] eqsys    pointer to the structure to log
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log(cs_equation_system_t  *eqsys)
{
  if (eqsys == NULL)
    return;

  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);
  cs_log_printf(CS_LOG_SETUP,
                "## Summary of settings for the system of equations: %s\n",
                eqsys->name);

  cs_log_printf(CS_LOG_SETUP, "  * %s | Number of equations: %d\n",
                eqsys->name, eqsys->n_equations);
  cs_log_printf(CS_LOG_SETUP, "  * %s | Equations: ", eqsys->name);
  for (int i = 0; i < eqsys->n_equations; i++) {

    cs_equation_t  *eq = eqsys->equations[i];
    if (eq != NULL)
      cs_log_printf(CS_LOG_SETUP, " %s (id=%d);", cs_equation_get_name(eq), i);

  }
  cs_log_printf(CS_LOG_SETUP, "\n");

  /* Log the setting of the extra-diagonal blocks */

  cs_log_printf(CS_LOG_SETUP, "  * %s | Settings for extra-diagonal blocks\n",
                eqsys->name);
  cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h2);

  for (int i = 0; i < eqsys->n_equations; i++) {
    for (int j = 0; j < eqsys->n_equations; j++) {
      if (i != j)
        cs_equation_param_log(eqsys->params[i*eqsys->n_equations + j]);
    }
  }

  cs_log_printf(CS_LOG_SETUP, "%s\n\n", cs_sep_h1);

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

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
                         cs_equation_system_t    *eqsys)
{
  if (eqsys == NULL)
    return;

  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  /* TODO */

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

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
                                cs_equation_system_t     *eqsys)
{
  if (eqsys == NULL)
    return;

  if (row_id >= eqsys->n_equations)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid row id %d (max. possible is %d)\n",
              __func__, row_id, eqsys->n_equations-1);

  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  eqsys->equations[row_id] = eq;

  if (eq == NULL)
    eqsys->params[row_id*eqsys->n_equations + row_id] = NULL;
  else
    eqsys->params[row_id*eqsys->n_equations + row_id] = eq->param;

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

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
                             cs_equation_system_t     *eqsys)
{
  if (eqsys == NULL)
    return;

  if (row_id >= eqsys->n_equations)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid row_id %d (max. possible is %d)\n",
              __func__, row_id, eqsys->n_equations-1);
  if (col_id >= eqsys->n_equations)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid col_id %d (max. possible is %d)\n",
              __func__, col_id, eqsys->n_equations-1);

  eqsys->params[row_id*eqsys->n_equations + col_id] = eqp;

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
