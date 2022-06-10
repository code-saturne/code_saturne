/*============================================================================
 * Functions to handle the cs_equation_system_t structure
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdovb_scalsys.h"
#include "cs_equation_system_sles.h"
#include "cs_equation_param.h"
#include "cs_timer_stats.h"

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
 * \brief Check the coherency of the settings between a cs_equation_system_t
 *        structure and each parameter structure associated to a block.
 *
 * \param[in] eqsys     pointer to the system of equations to set
 */
/*----------------------------------------------------------------------------*/

static void
_check_common_metadata(cs_equation_system_t  *eqsys)
{
  if (eqsys == NULL)
    return;

  cs_param_space_scheme_t  common_space_scheme = CS_SPACE_N_SCHEMES;
  int  common_var_dim = -1;

  int n_eqs = eqsys->n_equations;

  for (int i = 0; i < n_eqs; i++) {
    for  (int j = 0; j < n_eqs; j++) {

      cs_equation_core_t  *block = eqsys->block_factories[i*n_eqs+j];
      assert(block != NULL);

      const cs_equation_param_t  *eqp = block->param;

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

  if (eqsys->param->space_scheme != common_space_scheme)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible space scheme (system: %s; equations: %s)\n",
              __func__,
              cs_param_get_space_scheme_name(eqsys->param->space_scheme),
              cs_param_get_space_scheme_name(common_space_scheme));

  if (eqsys->param->block_var_dim != common_var_dim)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Incompatible var. dim. (system: %d; equations: %d)\n",
              __func__, eqsys->param->block_var_dim, common_var_dim);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to handle system of coupled equations
 *
 * \param[in] n_eqs      number of coupled equations composing the system
 * \param[in] sysname    name of the system of equations
 *
 * \return  a pointer to the new allocated cs_equation_system_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_equation_system_t *
_equation_system_create(int            n_eqs,
                        const char    *sysname)
{
  cs_equation_system_t  *eqsys = NULL;

  if (n_eqs < 2)
    return NULL;

  BFT_MALLOC(eqsys, 1, cs_equation_system_t);

  eqsys->n_equations = n_eqs;

  /* Monitoring:
   * Set timer statistic structure to a default value */

  CS_TIMER_COUNTER_INIT(eqsys->timer);

  eqsys->timer_id = cs_timer_stats_id_by_name(sysname);
  if (eqsys->timer_id < 0)
    eqsys->timer_id = cs_timer_stats_create(NULL, /* new root */
                                            sysname,
                                            sysname);

  /* Metadata */

  eqsys->param = NULL;

  /* Structures */

  eqsys->system_helper = NULL;

  BFT_MALLOC(eqsys->equations, n_eqs, cs_equation_t *);
  for (int i = 0; i < n_eqs; i++)
    eqsys->equations[i] = NULL;

  BFT_MALLOC(eqsys->block_factories, n_eqs*n_eqs, cs_equation_core_t *);
  for (int i = 0; i < n_eqs*n_eqs; i++)
    eqsys->block_factories[i] = NULL;

  /* Function pointers */

  eqsys->define = NULL;
  eqsys->free = NULL;
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

  int  n_eqs = eqsys->n_equations;

  eqsys->param = cs_equation_system_param_free(eqsys->param);

  /* Free all structures inside array of structures */

  eqsys->context = eqsys->free(n_eqs, eqsys->block_factories, eqsys->context);

  cs_cdo_system_helper_free(&(eqsys->system_helper));

  BFT_FREE(eqsys->block_factories);
  BFT_FREE(eqsys->equations);

  BFT_FREE(eqsys);
  *p_eqsys = NULL;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
cs_equation_system_get_n_systems(void)
{
  return _n_equation_systems;
}

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
                       int                     block_var_dim)
{
  if (n_eqs < 2)
    return NULL;

  /* Create a new structure */

  cs_equation_system_t  *eqsys = _equation_system_create(n_eqs, sysname);

  /* Add a set of parameters by default */

  eqsys->param = cs_equation_system_param_create(sysname, block_var_dim);

  int  sys_id = _n_equation_systems;
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
 * \brief Log the setup for all structures managing systems of equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log_setup(void)
{
  if (_n_equation_systems < 1)
    return;

  cs_log_printf(CS_LOG_SETUP, "\nSettings for systems of equations\n");
  cs_log_printf(CS_LOG_SETUP, "%s\n", cs_sep_h1);

  for (int sys_id = 0; sys_id < _n_equation_systems; sys_id++) {

    cs_equation_system_t  *eqsys = _equation_systems[sys_id];

    if (eqsys == NULL)
      continue;

    cs_timer_t  t1 = cs_timer_time();
    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    const char  *sysname = eqsys->param->name;
    const int  n_eqs = eqsys->n_equations;

    cs_log_printf(CS_LOG_SETUP,
                  "\nSummary of settings for the system of equations: %s\n",
                  sysname);
    cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h2);

    cs_log_printf(CS_LOG_SETUP, "  * %s | Number of equations: %d\n",
                  sysname, eqsys->n_equations);

    cs_equation_system_param_log(eqsys->param);

    cs_log_printf(CS_LOG_SETUP, "  * %s | Equations (diagonal blocks):\n",
                  sysname);

    for (int i = 0; i < n_eqs; i++) {

      cs_equation_t  *eq = eqsys->equations[i];
      if (eq != NULL)
        cs_log_printf(CS_LOG_SETUP, "\t%s (block_row_id=%d)\n",
                      cs_equation_get_name(eq), i);

    }

    /* Log the setting of the extra-diagonal blocks */

    cs_log_printf(CS_LOG_SETUP,
                  "\nSystem \"%s\": Settings for extra-diagonal blocks\n",
                  sysname);
    cs_log_printf(CS_LOG_SETUP, "%s", cs_sep_h2);

    for (int i = 0; i < n_eqs; i++) {
      for (int j = 0; j < n_eqs; j++) {
        if (i != j) {

          cs_equation_core_t  *block = eqsys->block_factories[i*n_eqs+j];

          cs_equation_param_log(block->param);

        }
      }
    }

    cs_timer_t  t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a synthesis of the monitoring information in the performance
 *         file
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_log_monitoring(void)
{
  cs_log_printf(CS_LOG_PERFORMANCE, "\n%-43s %9s\n", " ", "Total");

  for (int i = 0; i < _n_equation_systems; i++) {

    cs_equation_system_t  *eqsys = _equation_systems[i];
    assert(eqsys != NULL);
    cs_equation_system_param_t  *sysp = eqsys->param;

    /* Display high-level timer counter related to the current equation
       before deleting the structure */

    cs_log_printf(CS_LOG_PERFORMANCE,
                  " <CDO system/%20s> Runtime  %9.3f seconds\n",
                  sysp->name, eqsys->timer.nsec*1e-9);

  } /* Loop on equations */
}

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
                                const cs_time_step_t        *time_step)
{
  for (int i = 0; i < _n_equation_systems; i++) {

    cs_equation_system_t  *eqsys = _equation_systems[i];

    if (eqsys == NULL)
      bft_error(__FILE__, __LINE__, 0, "%s: System not allocated.", __func__);

    if (eqsys->n_equations < 1)
      return;

    cs_timer_t  t1 = cs_timer_time();
    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    /* Check if there is no issue with the settings between the system of
       equation and the settings related to each block. Check the coherency for
       the space scheme and block_var_dim */

    _check_common_metadata(eqsys);

    /* Now set the function pointers */

    cs_equation_system_param_t  *sysp = eqsys->param;
    assert(sysp != NULL);

    switch (sysp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (sysp->block_var_dim == 1)
        cs_cdovb_scalsys_init_sharing(mesh, connect, quant, time_step);
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid block_var_dim (=%d) for system \"%s\".\n"
                  "%s: Only scalar-valued (=1) blocks are handled.\n",
                  __func__, sysp->block_var_dim, sysp->name, __func__);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid space scheme (%s) for system \"%s\"\n",
                __func__, cs_param_get_space_scheme_name(sysp->space_scheme),
                sysp->name);
      break;

    } /* Switch on space scheme */

    cs_timer_t  t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing all the systems of
 *         equations
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

    cs_timer_t  t1 = cs_timer_time();
    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    /* Set associated function pointers */

    cs_equation_system_param_t  *sysp = eqsys->param;
    assert(sysp != NULL);

    switch (sysp->space_scheme) {

    case CS_SPACE_SCHEME_CDOVB:
      if (sysp->block_var_dim == 1) { /* Each block is scalar-valued  */

        eqsys->define = cs_cdovb_scalsys_define;
        eqsys->free = cs_cdovb_scalsys_free;

        eqsys->solve_steady_state_system = NULL; /* Not used up to now */
        eqsys->solve_system = cs_cdovb_scalsys_solve_implicit;

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid block_var_dim (=%d) for system \"%s\".\n"
                  "%s: Only scalar-valued (=1) blocks are handled.\n",
                  __func__, sysp->block_var_dim, sysp->name, __func__);
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid space scheme (%s) for system \"%s\"\n",
                __func__, cs_param_get_space_scheme_name(sysp->space_scheme),
                sysp->name);
      break;

    } /* Switch on space scheme */

    cs_timer_t  t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the SLES associated to each system of equations
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_system_set_sles(void)
{
  for (int sys_id = 0; sys_id < _n_equation_systems; sys_id++) {

    cs_equation_system_t  *eqsys = _equation_systems[sys_id];

    if (eqsys == NULL)
      bft_error(__FILE__, __LINE__, 0, "%s: System not allocated.", __func__);

    cs_equation_system_param_t  *sysp = eqsys->param;
    assert(sysp != NULL);

    cs_timer_t  t1 = cs_timer_time();
    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    cs_equation_system_sles_init(eqsys->n_equations,
                                 sysp,
                                 eqsys->block_factories);

    cs_timer_t  t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems */
}

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
cs_equation_system_define(void)
{
  for (int sys_id = 0; sys_id < _n_equation_systems; sys_id++) {

    cs_equation_system_t  *eqsys = _equation_systems[sys_id];

    if (eqsys == NULL)
      bft_error(__FILE__, __LINE__, 0, "%s: System not allocated.", __func__);

    const int  n_eqs = eqsys->n_equations;
    const cs_equation_system_param_t  *sysp = eqsys->param;
    assert(sysp != NULL);

    cs_timer_t  t1 = cs_timer_time();
    if (eqsys->timer_id > -1)
      cs_timer_stats_start(eqsys->timer_id);

    for (int i = 0; i < n_eqs; i++) {

      int ii = i*n_eqs+i;

      const cs_equation_t  *eq = eqsys->equations[i];
      assert(eq != NULL);

      cs_equation_core_t  *block_ii = eqsys->block_factories[ii];
      assert(block_ii != NULL);

      cs_equation_define_core(eq, &block_ii);

    } /* Loop on equations (Diagonal blocks) */

    eqsys->context = eqsys->define(n_eqs, sysp, eqsys->block_factories,
                                   &eqsys->system_helper);

    cs_timer_t  t2 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
    if (eqsys->timer_id > -1)
      cs_timer_stats_stop(eqsys->timer_id);

  } /* Loop on systems of equations */

  /* Allocate or update the low-level assemble structures */

  cs_cdo_system_allocate_assembly();
}

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
                         cs_equation_system_t    *eqsys)
{
  if (eqsys == NULL)
    return;

  cs_timer_t  t1 = cs_timer_time();
  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  if (eqsys->solve_system == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: No solve function set for system \"%s\"\n",
              __func__, (eqsys->param == NULL) ? NULL : eqsys->param->name);

  /* One assumes that by default the matrix structure is not stored. So one
     has to build this structure before each solving step */

  eqsys->solve_system(cur2prev,
                      eqsys->n_equations,
                      eqsys->param,
                      eqsys->block_factories,
                      eqsys->context,
                      eqsys->system_helper);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

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
                                   cs_equation_system_t     *eqsys)
{
  if (eqsys == NULL)
    return;

  cs_timer_t  t1 = cs_timer_time();
  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  int  n_eqs = eqsys->n_equations;

  if (row_id >= n_eqs)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid row id %d (max. possible is %d)\n",
              __func__, row_id, n_eqs-1);

  eqsys->equations[row_id] = eq;

  /* Set what is already available as structure pointers */

  cs_equation_core_t  *block_ii = NULL;
  cs_equation_define_core(eq, &block_ii);
  eqsys->block_factories[row_id*n_eqs + row_id] = block_ii;

  block_ii->param->flag |= CS_EQUATION_INSIDE_SYSTEM;

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

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
                                cs_equation_system_t     *eqsys)
{
  if (eqsys == NULL)
    return;
  if (eqp == NULL)
    return;

  cs_timer_t  t1 = cs_timer_time();
  if (eqsys->timer_id > -1)
    cs_timer_stats_start(eqsys->timer_id);

  int  n_eqs = eqsys->n_equations;

  if (row_id >= n_eqs)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid row_id %d (max. possible is %d)\n",
              __func__, row_id, n_eqs-1);
  if (col_id >= n_eqs)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid col_id %d (max. possible is %d)\n",
              __func__, col_id, n_eqs-1);

  const char  *sysname = (eqsys->param == NULL) ? NULL : eqsys->param->name;

  if (eqsys->block_factories[row_id*n_eqs + col_id] != NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The block (%d, %d) has already been assigned in"
              " system \"%s\"\n", __func__, row_id, col_id, sysname);

  cs_equation_core_t  *block_ij = NULL;

  BFT_MALLOC(block_ij, 1, cs_equation_core_t);

  eqp->flag |= CS_EQUATION_INSIDE_SYSTEM;

  block_ij->param = eqp;
  block_ij->builder = NULL;
  block_ij->scheme_context = NULL;

  eqsys->block_factories[row_id*n_eqs + col_id] = block_ij;

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqsys->timer), &t1, &t2);
  if (eqsys->timer_id > -1)
    cs_timer_stats_stop(eqsys->timer_id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
