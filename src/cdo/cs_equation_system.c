/*============================================================================
 * Functions to handle the cs_equation_system_t structure
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

static int  _n_equation_systems = 0;
static cs_equation_system_t  **_equation_systems = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set metadata associated to a cs_equatino_system_t structure. This
 *         has to be done when cs_equation_parameter_t structure have been all
 *         set.
 *
 * \param[in] eqsys     pointer to the system of equations to set
 */
/*----------------------------------------------------------------------------*/

static void
_set_common_metadata(cs_equation_system_t  *eqsys)
{
  if (eqsys == NULL)
    return;

  cs_param_space_scheme_t  common_space_scheme = CS_SPACE_N_SCHEMES;
  int  common_var_dim = -1;

  for (int i = 0; i < eqsys->n_equations; i++) {
    for  (int j = 0; j < eqsys->n_equations; j++) {

      cs_equation_param_t  *eqp = eqsys->params[i*eqsys->n_equations + j];

      if (common_var_dim == -1)
        common_var_dim = eqp->dim;
      else {

        if (common_var_dim != eqp->dim)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Incompatible var. dim. (current: %d; previous: %d)\n",
                    __func__, eqp->dim, common_var_dim);

      }

      if (common_space_scheme == CS_SPACE_N_SCHEMES)
        common_space_scheme = eqp->space_scheme;
      else {

        if (common_space_scheme != eqp->space_scheme)
          bft_error(__FILE__, __LINE__, 0,
                    "%s: Incompatible space scheme (current: %s; previous: %s)"
                    "\n", __func__,
                    cs_param_get_space_scheme_name(common_space_scheme),
                    cs_param_get_space_scheme_name(eqp->space_scheme));

      }

    }
  }

  eqsys->space_scheme = common_space_scheme;
  eqsys->block_var_dim = common_var_dim;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to handle system of coupled equations
 *
 * \param[in] sysname       name of the system of equations
 * \param[in] n_eqs         number of coupled equations composing the system
 *
 * \return  a pointer to the new allocated cs_equation_system_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_equation_system_t *
_equation_system_create(const char                *sysname,
                        int                        n_eqs)
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
  eqsys->space_scheme = CS_SPACE_N_SCHEMES; /* not set */
  eqsys->block_var_dim = -1;                /* not set */

  BFT_MALLOC(eqsys->equations, n_eqs, cs_equation_t *);
  for (int i = 0; i < n_eqs; i++)
    eqsys->equations[i] = NULL;

  BFT_MALLOC(eqsys->params, n_eqs*n_eqs, cs_equation_param_t *);
  for (int i = 0; i < n_eqs*n_eqs; i++)
    eqsys->params[i] = NULL;

  BFT_MALLOC(eqsys->builders, n_eqs*n_eqs, cs_equation_builder_t *);
  for (int i = 0; i < n_eqs*n_eqs; i++)
    eqsys->builders[i] = NULL;

  BFT_MALLOC(eqsys->context_structures, n_eqs*n_eqs, void *);
  for (int i = 0; i < n_eqs*n_eqs; i++)
    eqsys->context_structures[i] = NULL;

  eqsys->init_context = NULL;
  eqsys->free_context = NULL;
  eqsys->solve_system = NULL;
  eqsys->solve_steady_state_system = NULL;

  return eqsys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a structure used to couple equations
 *
 * \param[in, out] p_eqsys    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

static void
_equation_system_free(cs_equation_system_t  **p_eqsys)
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

  if (n_eqs > 0) {

    /* Free the extra-diagonal cs_equation_param_t structures, builders and
       scheme context structures */

    for (int i = 0; i < n_eqs; i++) {

      cs_equation_param_t  **eqp_array = eqsys->params + i*n_eqs;
      cs_equation_builder_t  **builders = eqsys->builders + i*n_eqs;
      void  **context_structures = eqsys->context_structures + i*n_eqs;

      for (int j = 0; j < n_eqs; j++) {
        if (i != j) {

          eqp_array[j] = cs_equation_param_free(eqp_array[j]);
          cs_equation_builder_free(builders + j);
          eqsys->free_context(context_structures + j);

        }
      } /* Loop on equations (j) */

    } /* Loop on equations (i) */

    BFT_FREE(eqsys->params);
    BFT_FREE(eqsys->builders);
    BFT_FREE(eqsys->context_structures);

    BFT_FREE(eqsys->equations);

  } /* n_eqs > 0 */

  BFT_FREE(eqsys);
  *p_eqsys = NULL;

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new structure to handle system of coupled equations
 *
 * \param[in] sysname       name of the system of equations
 * \param[in] n_eqs         number of coupled equations composing the system
 *
 * \return  a pointer to the new allocated cs_equation_system_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_system_t *
cs_equation_system_add(const char                *sysname,
                       int                        n_eqs)
{
  int  sys_id = _n_equation_systems;

  cs_equation_system_t  *eqsys = _equation_system_create(sysname, n_eqs);

  _n_equation_systems++;
  BFT_REALLOC(_equation_systems, _n_equation_systems, cs_equation_system_t *);

  _equation_systems[sys_id] = eqsys;

  return eqsys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all structures used to couple equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_destroy_all(void)
{
  for (int i = 0; i < _n_equation_systems; i++)
    _equation_system_free(_equation_systems + i);

  BFT_FREE(_equation_systems);
  _equation_systems = NULL;
  _n_equation_systems = 0;
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
  cs_log_printf(CS_LOG_SETUP, "  * %s | Common space scheme: %s\n",
                eqsys->name,
                cs_param_get_space_scheme_name(eqsys->space_scheme));
  cs_log_printf(CS_LOG_SETUP, "  * %s | Common variable dimension: %d\n",
                eqsys->name, eqsys->block_var_dim);
  cs_log_printf(CS_LOG_SETUP, "  * %s | Equations (diagonal blocks): ",
                eqsys->name);

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
 * \brief  Assign a set of pointer functions for managing the
 *         cs_equation_system_t structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_functions(void)
{
  for (int i = 0; i < _n_equation_systems; i++) {

    cs_equation_system_t  *eqsys = _equation_systems[i];

    if (eqsys == NULL)
      bft_error(__FILE__, __LINE__, 0, "%s: System not allocated.", __func__);

    if (eqsys->n_equations < 1)
      return;

    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    /* Set the space scheme and block_var_dim for the system of equations */

    _set_common_metadata(eqsys);

    /* Now set the function pointers */

    switch (eqsys->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (eqsys->block_var_dim == 1) { /* Each block is scalar-valued  */

        /* Set the solve functions */

        eqsys->solve_steady_state_system = NULL; /* Not used up to now */
        eqsys->solve_system = NULL;              /* To be set */

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid block_var_dim (=%d) for system \"%s\".\n"
                  "%s: Only scalar-valued (=1) blocks are handled.\n",
                  __func__, eqsys->block_var_dim, eqsys->name, __func__);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid space scheme (%s) for system \"%s\"\n",
                __func__, cs_param_get_space_scheme_name(eqsys->space_scheme),
                eqsys->name);
      break;

    } /* Switch on space scheme */

    /* The first equation should have the same init_context and free_context
       functions as all other equations since they share the space
       discretization and the same dimension of variable */

    eqsys->init_context = eqsys->equations[0]->init_context;
    eqsys->free_context = eqsys->equations[0]->free_context;

    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize builder and scheme context structures associated to all
 *         the systems of equations which have been added
 *
 * \param[in]       mesh      pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_initialize(const cs_mesh_t             *mesh)
{
  for (int sys_id = 0; sys_id < _n_equation_systems; sys_id++) {

    cs_equation_system_t  *eqsys = _equation_systems[sys_id];

    if (eqsys == NULL)
      bft_error(__FILE__, __LINE__, 0, "%s: System not allocated.", __func__);

    const int n_eqs = eqsys->n_equations;

    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    for (int i = 0; i < n_eqs; i++) {
      for (int j = 0; j < n_eqs; j++) {

        if (i == j) { /* Diagonal block => associated to an equation */

          cs_equation_t  *eq = eqsys->equations[i];
          assert(eq != NULL);

          if (eq->builder == NULL)  /* This should not be the case */
            eq->builder = cs_equation_builder_init(eq->param, mesh);

          eqsys->builders[i*n_eqs+i] = eq->builder;

          if (eq->scheme_context == NULL) /* This should not be the case */
            eq->scheme_context = eq->init_context(eq->param,
                                                  eq->field_id,
                                                  eq->boundary_flux_id,
                                                  eq->builder);

          eqsys->context_structures[i*n_eqs+i] = eq->scheme_context;

        }
        else { /* Extra-diagonal block */

          int  ij = i*n_eqs+j;
          cs_equation_param_t  *eqp = eqsys->params[ij];

          assert(eqp != NULL);
          assert(eqsys->builders[ij] == NULL);
          assert(eqsys->context_structures[ij] == NULL);

          cs_equation_builder_t  *eqb = cs_equation_builder_init(eqp, mesh);

          eqsys->builders[ij] = eqb;
          eqsys->context_structures[ij] = eqsys->init_context(eqp,
                                                              -1, /* No field */
                                                              -1, /* No field */
                                                              eqb);

        }

      } /* column j */
    } /* row i */

    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */
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

  if (eqsys->solve_system == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: No solve function set for system \"%s\"\n",
              __func__, eqsys->name);

  eqsys->solve_system(cur2prev,
                      eqsys->n_equations,
                      mesh,
                      eqsys->params,
                      eqsys->builders,
                      eqsys->context_structures);

  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the given equation and associate it to the row block with id
 *         equal to row_id. The equation parameter is also set.
 *
 * \param[in]      row_id  position in the block matrix
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

  else {

    eqsys->params[row_id*eqsys->n_equations + row_id] = eq->param;

  }

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

  if (eqp != NULL) {

    if (eqsys->space_scheme == CS_SPACE_N_SCHEMES)
      eqsys->space_scheme = eqp->space_scheme;
    else
      if (eqsys->space_scheme != eqp->space_scheme)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Stop adding an equation to the system \"%s\".\n"
                  "%s: Space discretization differs.",
                  __func__, eqsys->name, __func__);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
