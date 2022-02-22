/*============================================================================
 * Functions to handle the cs_equation_system_param_t structure
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_system_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define CS_EQUATION_SYSTEM_PARAM_DBG  0

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
                                int               block_var_dim)
{
  cs_equation_system_param_t  *sysp = NULL;

  BFT_MALLOC(sysp, 1, cs_equation_system_param_t);

  sysp->verbosity = 1;

  /* Store name */

  size_t  len = strlen(name);
  BFT_MALLOC(sysp->name, len + 1, char);
  strncpy(sysp->name, name, len);
  sysp->name[len] = '\0';

  /* Dimension of the variable for each block */

  assert(block_var_dim > 0);
  sysp->block_var_dim = block_var_dim;

  /* Space discretization */

  sysp->space_scheme = CS_SPACE_SCHEME_CDOVB;

  /* Linear algebra settings by default */

  sysp->sles_setup_done = false;
  sysp->sles_strategy = CS_EQUATION_SYSTEM_SLES_MUMPS;

  sysp->linear_solver.n_max_algo_iter = 100;
  sysp->linear_solver.rtol = 1e-06;
  sysp->linear_solver.atol = 1e-08;
  sysp->linear_solver.dtol = 1e3;
  sysp->linear_solver.verbosity = 1;

  return sysp;
}

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
cs_equation_system_param_free(cs_equation_system_param_t    *sysp)
{
  if (sysp == NULL)
    return sysp;

  BFT_FREE(sysp->name);

  BFT_FREE(sysp);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Log the setup gathered in the structure cs_equation_system_param_t
 *
 * \param[in] sysp     pointer to a parameter structure to log
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_param_log(const cs_equation_system_param_t    *sysp)
{
  if (sysp == NULL)
    return;

  char desc[128];
  sprintf(desc, "  * %s |", sysp->name);

  cs_log_printf(CS_LOG_SETUP, "%s Verbosity: %d\n", desc, sysp->verbosity);
  cs_log_printf(CS_LOG_SETUP, "%s Common space scheme: %s\n",
                desc, cs_param_get_space_scheme_name(sysp->space_scheme));
  cs_log_printf(CS_LOG_SETUP, "%s Common variable dimension: %d\n",
                desc, sysp->block_var_dim);

  cs_log_printf(CS_LOG_SETUP, "%s Linear algebra setup\n", desc);
  switch (sysp->sles_strategy) {

  case CS_EQUATION_SYSTEM_SLES_MUMPS:
    cs_log_printf(CS_LOG_SETUP, "%s Strategy: Full MUMPS\n", desc);
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Strategy: Unknown\n", desc);
    break;

  } /* Switch on strategy */

  cs_log_printf(CS_LOG_SETUP, "%s Tolerances of the linear solver:"
                " rtol: %5.3e; atol: %5.3e; dtol: %5.3e; verbosity: %d\n",
                desc, sysp->linear_solver.rtol, sysp->linear_solver.atol,
                sysp->linear_solver.dtol, sysp->linear_solver.verbosity);
  cs_log_printf(CS_LOG_SETUP, "%s Max number of iterations: %d\n",
                desc, sysp->linear_solver.n_max_algo_iter);
}

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
                             const char                    *keyval)
{
  if (sysp == NULL)
    return;
  if (keyval == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty key value.\n", __func__);

  /* Conversion of the string to lower case */

  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_SYSKEY_LINEAR_SOLVER_ATOL:
    sysp->linear_solver.atol = atof(val);
    if (sysp->linear_solver.atol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the absolute tolerance"
                " of the linear solver\n", __func__);
    break;

  case CS_SYSKEY_LINEAR_SOLVER_DTOL:
    sysp->linear_solver.dtol = atof(val);
    if (sysp->linear_solver.dtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the divergence tolerance"
                " of the linear solver\n", __func__);
    break;

  case CS_SYSKEY_LINEAR_SOLVER_RTOL:
    sysp->linear_solver.rtol = atof(val);
    if (sysp->linear_solver.rtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the divergence tolerance"
                " of the linear solver\n", __func__);
    break;

  case CS_SYSKEY_LINEAR_SOLVER_VERBOSITY:
    sysp->linear_solver.verbosity = atoi(val);
    break;

  case CS_SYSKEY_LINEAR_SOLVER_MAX_ITER:
    sysp->linear_solver.n_max_algo_iter = atoi(val);
    break;

  case CS_SYSKEY_SLES_STRATEGY:
    if (strcmp(val, "mumps") == 0) {
#if defined(HAVE_MUMPS)
      sysp->sles_strategy = CS_EQUATION_SYSTEM_SLES_MUMPS;
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
      sysp->sles_strategy = CS_EQUATION_SYSTEM_SLES_MUMPS;
#else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting \"%s\" key\n"
                " MUMPS is not available with your installation.\n"
                " Please check your installation settings.\n",
                __func__, "CS_SYSKEY_SLES_STRATEGY");
#endif  /* PETSC_HAVE_MUMPS */
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid val %s related to key CS_SYSKEY_SLES_STRATEGY\n"
                " Choice between: mumps\n",
                __func__, _val);
    }
    break;

  case CS_SYSKEY_VERBOSITY:
    sysp->verbosity = atoi(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid key for setting the equation system \"%s\".",
              __func__, sysp->name);

  } /* End of switch */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
