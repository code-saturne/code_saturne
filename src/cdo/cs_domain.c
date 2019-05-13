/*============================================================================
 * Manage a computational domain within the CDO framework
 *  - Properties and advection fields attached to this domain
 *  - Equations to solve on this domain
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"
#include "cs_quadrature.h"
#include "cs_restart.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_domain.c

  \brief  Manage a computational domain
    - Settings, fields, connectivities and geometrical quantities
    - Properties and advection fields attached to this domain
    - Equations to solve on this domain
*/

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_domain_t  *cs_glob_domain = NULL; /* Pointer to the main computational
                                        domain */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

static const char _err_empty_domain[] =
  " Stop setting an empty cs_domain_t structure.\n"
  " Please check your settings.\n";

static double  cs_domain_kahan_time_compensation = 0.0;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the context for CDO/HHO schemes
 *
 * \param[in]  cdo_mode         type of activation for the CDO/HHO module
 *
 * \return a pointer to a new allocated cs_domain_cdo_context_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_cdo_context_t *
_create_cdo_context(int     cdo_mode)
{
  cs_domain_cdo_context_t  *cc = NULL;

  BFT_MALLOC(cc, 1, cs_domain_cdo_context_t);

  cc->mode = cdo_mode;

  cc->force_advfield_update = false;

  /* Metadata related to each family of schemes */
  cc->vb_scheme_flag = 0;
  cc->vcb_scheme_flag = 0;
  cc->fb_scheme_flag = 0;
  cc->hho_scheme_flag = 0;

  return cc;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize by default a cs_domain_t structure
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_create(void)
{
  cs_domain_t  *domain = NULL;

  BFT_MALLOC(domain, 1, cs_domain_t);

  domain->mesh = NULL;
  domain->mesh_quantities = NULL;
  domain->connect = NULL;
  domain->cdo_quantities = NULL;
  domain->cdo_context = NULL;

  /* Default initialization of the time step */
  domain->only_steady = true;
  domain->is_last_iter = false;
  domain->time_step_def = NULL;

  /* Global structure for time step management */
  domain->time_step = cs_get_glob_time_step();

  /* Other options */
  domain->restart_nt = 0;
  domain->output_nt = -1;
  domain->verbosity = 1;

  /* By default: CDO-HHO schemes are not activated */
  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_OFF);

  /* By default a wall is defined for the whole boundary of the domain */
  cs_glob_boundaries = cs_boundary_create(CS_BOUNDARY_WALL);
  domain->boundaries = cs_glob_boundaries;
  domain->ale_boundaries = cs_boundary_create(CS_BOUNDARY_ALE_FIXED);

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(domain->tcp); /* domain post */
  CS_TIMER_COUNTER_INIT(domain->tcs); /* domain setup */

  /* Initialization of several modules */
  cs_math_set_machine_epsilon(); /* Compute and set machine epsilon */
  cs_quadrature_setup();         /* Compute constant used in quadrature rules */

  return domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_domain_t structure
 *
 * \param[in, out]   p_domain    pointer of pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_free(cs_domain_t   **p_domain)
{
  if (p_domain == NULL)
    return;

  cs_domain_t  *domain = *p_domain;

  /* cs_mesh_t and cs_mesh_quantities_t structure are not freed since they
     are only shared */
  domain->mesh = NULL;
  domain->mesh_quantities = NULL;

  domain->time_step_def = cs_xdef_free(domain->time_step_def);
  domain->time_step = NULL;

  if (domain->cdo_context != NULL)
    BFT_FREE(domain->cdo_context);

  /* Free arrays related to the domain boundary */
  cs_boundary_free(&(domain->boundaries));
  cs_boundary_free(&(domain->ale_boundaries));

  /* Free CDO structures related to geometric quantities and connectivity */
  domain->cdo_quantities = cs_cdo_quantities_free(domain->cdo_quantities);
  domain->connect = cs_cdo_connect_free(domain->connect);

  BFT_FREE(domain);
  *p_domain = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the global variable storing the mode of activation to apply
 *          to CDO/HHO schemes
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 * \param[in]        mode      type of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_cdo_mode(cs_domain_t    *domain,
                       int             mode)
{
  if (domain == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: domain is not allocated.",
              __func__);

  if (domain->cdo_context == NULL)
    domain->cdo_context = _create_cdo_context(mode);
  else
    domain->cdo_context->mode = mode;

  CS_PROCF(set_cdo_mode, SET_CDO_MODE)(&mode);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the mode of activation for the CDO/HHO schemes
 *
 * \param[in]   domain       pointer to a cs_domain_t structure
 *
 * \return the mode of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

int
cs_domain_get_cdo_mode(const cs_domain_t   *domain)
{
  if (domain == NULL)
    return CS_DOMAIN_CDO_MODE_OFF;
  if (domain->cdo_context == NULL)
    return CS_DOMAIN_CDO_MODE_OFF;

  return domain->cdo_context->mode;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one needs to continue iterations in time
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_iteration(cs_domain_t  *domain)
{
  bool  one_more_iter = true;

  cs_time_step_t  *ts = domain->time_step;

  if (ts->nt_max > 0) // nt_max has been set
    if (ts->nt_cur > ts->nt_max)
      one_more_iter = false;

  if (ts->t_max > 0) // t_max has been set
    if (ts->t_cur >= ts->t_max)
      one_more_iter = false;

  if (domain->only_steady)
    one_more_iter = false;

  if (!domain->only_steady && ts->nt_max <= 0 && ts->t_max <= 0)
    one_more_iter = false;

  return one_more_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if an ouput is requested according to the domain setting
 *
 * \param[in]   domain    pointer to a cs_domain_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_log(const cs_domain_t      *domain)
{
  const cs_time_step_t  *ts = domain->time_step;

  if (domain->verbosity < 0)
    return false;

  if (domain->only_steady)
    return true;

  if (domain->output_nt > 0)
    if (ts->nt_cur % domain->output_nt == 0)
      return true;

  if (domain->is_last_iter)
    return true;

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the current time step for this new time iteration
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_define_current_time_step(cs_domain_t   *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  if (domain->only_steady)
    return;

  cs_time_step_t  *ts = domain->time_step;
  cs_xdef_t  *ts_def = domain->time_step_def;

  if (ts_def == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Please check your settings: Unsteady computation but no"
              " current time step defined.\n", __func__);

  const double  t_cur = ts->t_cur;
  const int  nt_cur = ts->nt_cur;

  if (ts_def->type != CS_XDEF_BY_VALUE) { /* dt_cur may change */

    ts->dt[2] = ts->dt[1];
    ts->dt[1] = ts->dt[0];

    if (ts_def->type == CS_XDEF_BY_TIME_FUNCTION) {

      /* Compute current time step */
      cs_xdef_time_func_input_t  *param
        = (cs_xdef_time_func_input_t *)ts_def->input;
      param->func(nt_cur, t_cur, param->input, &(ts->dt[0]));

      /* Update time_options */
      double  dtmin = CS_MIN(domain->time_options.dtmin, ts->dt[0]);
      double  dtmax = CS_MAX(domain->time_options.dtmax, ts->dt[0]);

      domain->time_options.dtmin = dtmin;
      domain->time_options.dtmax = dtmax;

      /* TODO: Check how the following value is set in FORTRAN
       * domain->time_options.dtref = 0.5*(dtmin + dtmax); */
      if (domain->time_step->dt_ref < 0) /* Should be the initial val. */
        domain->time_step->dt_ref = ts->dt[0];

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid way of defining the current time step.\n"
                " Please modify your settings.", __func__);

  }

  /* Check if this is the last iteration */
  if (ts->t_max > 0) /* t_max has been set */
    if (t_cur + ts->dt[0] > ts->t_max)
      domain->is_last_iter = true;
  if (ts->nt_max > 0) /* nt_max has been set */
    if (nt_cur + 1 > ts->nt_max)
      domain->is_last_iter = true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the simulated time after one temporal iteration
 *
 * \param[in, out]  domain     pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_increment_time(cs_domain_t  *domain)
{
  cs_time_step_t  *ts = domain->time_step;

  /* Use Kahan's trick to limit the truncation error */
  double  z = ts->dt[0] - cs_domain_kahan_time_compensation;
  double  t = ts->t_cur + z;

  cs_domain_kahan_time_compensation = (t - ts->t_cur) - z;
  ts->t_cur = t;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a welcome message indicating which mode of CDO is activated
 *
 * \param[in]  domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_cdo_log(const cs_domain_t   *domain)
{
  if (domain == NULL) bft_error(__FILE__, __LINE__, 0, _err_empty_domain);

  int  cdo_mode = CS_DOMAIN_CDO_MODE_OFF;
  if (domain->cdo_context != NULL)
    cdo_mode = domain->cdo_context->mode;

  switch (cdo_mode) {

  case CS_DOMAIN_CDO_MODE_ONLY:
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n -msg- CDO/HHO module is activated *** Experimental ***"
                  "\n -msg- CDO/HHO module is in a stand-alone mode\n");
    break;

  case CS_DOMAIN_CDO_MODE_WITH_FV:
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n -msg- CDO/HHO module is activated *** Experimental ***"
                  "\n -msg- CDO/HHO module with FV schemes mode\n");
    break;

  default:
  case CS_DOMAIN_CDO_MODE_OFF:
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n -msg- CDO/HHO module is not activated\n");
    break;

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
