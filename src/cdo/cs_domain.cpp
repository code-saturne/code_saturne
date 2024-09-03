/*============================================================================
 * Manage a computational domain
 *  - Mesh quantities and connectivities
 *  - Properties and advection fields attached to this domain
 *  - Equations to solve on this domain
 *============================================================================*/

/* VERS */

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
#include "cs_coupling.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_property.h"
#include "cs_prototypes.h"
#include "cs_quadrature.h"
#include "cs_restart.h"
#include "cs_time_step.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_domain.h"

/*----------------------------------------------------------------------------*/

/*!
  \file cs_domain.c

  \brief  Manage a computational domain
    - Settings, fields, connectivities and geometrical quantities
    - Properties and advection fields attached to this domain
    - Equations to solve on this domain
*/

BEGIN_C_DECLS

/*============================================================================
 * Static global variables
 *============================================================================*/

cs_domain_t *cs_glob_domain = nullptr; /* Pointer to the main computational
                                       domain */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

static double  cs_domain_kahan_time_compensation = 0.0;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function defintitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a new string with the name equal to "Unknown"
 *
 * \return a newly allocated string
 */
/*----------------------------------------------------------------------------*/

static char *
_set_to_unknown(void)
{
  char *name = nullptr;
  int len = strlen("Unknown");

  BFT_MALLOC(name, len + 1, char);
  sprintf(name, "%s", "Unknown");

  return name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the context for CDO/HHO schemes
 *
 * \return a pointer to a new allocated cs_domain_cdo_context_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_domain_cdo_context_t *
_create_cdo_context(void)
{
  cs_domain_cdo_context_t *cc = nullptr;

  BFT_MALLOC(cc, 1, cs_domain_cdo_context_t);

  /* Metadata related to each family of schemes. By default, not used */

  cc->vb_scheme_flag  = 0;
  cc->vcb_scheme_flag = 0;
  cc->eb_scheme_flag  = 0;
  cc->fb_scheme_flag  = 0;
  cc->cb_scheme_flag  = 0;

  cc->hho_scheme_flag = 0;
  cc->mac_scheme_flag = 0;

  return cc;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create and initialize by default a cs_domain_t structure
 *
 * \return a pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

cs_domain_t *
cs_domain_create(void)
{
  cs_domain_t *domain = nullptr;

  /* Initialization of several modules */

  cs_quadrature_setup();         /* Compute constant used in quadrature rules */

  /* Add two predefined properties which can be called from everywhere:
   *  1. the unity property
   *  2. the time_step property
   *
   * Simply call cs_property_by_name("pty_name"); to retrieve the pointer to
   * the related property structure
   */

  cs_property_t  *unity = cs_property_add("unity", CS_PROPERTY_ISO);
  cs_property_def_constant_value(unity, 1.0);

  cs_property_t  *dt_pty = cs_property_add("time_step", CS_PROPERTY_ISO);
  cs_property_set_reference_value(dt_pty, -1); /* Default=-1 => steady-state */

  /* Create the domain structure and proceed to a default initialization */

  BFT_MALLOC(domain, 1, cs_domain_t);

  /* Working directory names */

  domain->run_id     = nullptr;
  domain->case_name  = nullptr;
  domain->study_name = nullptr;

  cs_base_get_run_identity(&(domain->run_id),
                           &(domain->case_name),
                           &(domain->study_name));

  if (domain->run_id == nullptr)
    domain->run_id = _set_to_unknown();
  if (domain->case_name == nullptr)
    domain->case_name = _set_to_unknown();
  if (domain->study_name == nullptr)
    domain->study_name = _set_to_unknown();

  /* Quantitities and connectivities associated to the mesh and a numerical
     scheme */

  domain->mesh            = nullptr;
  domain->mesh_quantities = nullptr;
  domain->connect         = nullptr;
  domain->cdo_quantities  = nullptr;

  /* By default a wall is defined for the whole boundary of the domain */

  cs_glob_boundaries = cs_boundary_create(CS_BOUNDARY_CATEGORY_FLOW,
                                          CS_BOUNDARY_WALL);
  domain->boundaries = cs_glob_boundaries;
  domain->ale_boundaries = cs_boundary_create(CS_BOUNDARY_CATEGORY_ALE,
                                              CS_BOUNDARY_ALE_FIXED);

  /* Default initialization of the time step */

  domain->only_steady = true;
  domain->is_last_iter = false;
  domain->stage = CS_DOMAIN_STAGE_BEFORE_STEADY_COMPUTATION;

  /* Global structure for time step management */

  domain->time_step = cs_get_glob_time_step();

  domain->time_options.iptlro = 0;
  domain->time_options.idtvar
    = CS_TIME_STEP_CONSTANT; /* constant time step by default */
  domain->time_options.coumax = 1.;
  domain->time_options.cflmmx = 0.99;
  domain->time_options.foumax = 10.;
  domain->time_options.varrdt = 0.1;
  domain->time_options.dtmin  = -1.e13;
  domain->time_options.dtmax  = -1.e13;
  domain->time_options.relxst = 0.7; /* Not used in CDO schemes */

  /* Other options */

  domain->restart_nt = CS_RESTART_INTERVAL_DEFAULT;
  domain->output_nt = -10;  /* automatic with 10 first iterations */
  domain->verbosity = 1;    /* default verbosity */

  /* By default: CDO/HHO schemes are not activated (see cs_param_cdo.c) */

  domain->cdo_context = _create_cdo_context();

  /* Monitoring */

  CS_TIMER_COUNTER_INIT(domain->tcp); /* domain post */
  CS_TIMER_COUNTER_INIT(domain->tca); /* domain all op. */

  return domain;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a cs_domain_t structure
 *
 * \param[in, out] p_domain    pointer of pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_free(cs_domain_t   **p_domain)
{
  if (p_domain == nullptr)
    return;

  cs_domain_t  *domain = *p_domain;

  /* Working directory names */

  BFT_FREE(domain->run_id);
  BFT_FREE(domain->case_name);
  BFT_FREE(domain->study_name);

  /* cs_mesh_t and cs_mesh_quantities_t structure are not freed since they
     are only shared */

  domain->mesh            = nullptr;
  domain->mesh_quantities = nullptr;

  domain->time_step = nullptr;

  if (domain->cdo_context != nullptr)
    BFT_FREE(domain->cdo_context);

  /* Free arrays related to the domain boundary */

  cs_boundary_free(&(domain->boundaries));
  cs_boundary_free(&(domain->ale_boundaries));

  BFT_FREE(domain);
  *p_domain = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the global variable storing the mode of activation to apply to
 *        CDO/HHO schemes. Deprecated way to set the CDO mode.
 *
 * \param[in, out] domain    pointer to a cs_domain_t structure
 * \param[in]      mode      type of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_cdo_mode(cs_domain_t *domain, cs_param_cdo_mode_t mode)
{
  CS_NO_WARN_IF_UNUSED(domain);
  cs_param_cdo_mode_set(mode);  /* New way */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the mode of activation for the CDO/HHO schemes. Deprecated way
 *        to retrieve the CDO mode.
 *
 * \param[in] domain       pointer to a cs_domain_t structure
 *
 * \return the mode of activation for the CDO/HHO module
 */
/*----------------------------------------------------------------------------*/

int
cs_domain_get_cdo_mode(const cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);
  return cs_param_cdo_mode_get();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters related to the way output/logging is done
 *
 * \param[in, out] domain      pointer to a cs_domain_t structure
 * \param[in]      restart_nt  frequency for the restart process
 * \param[in]      log_nt      output frequency into the log (> 0: constant with
 *                             the given frequency, < 0: automatic with log_nt
 *                             first iterations detailed, 0: do nothing)
 * \param[in]      verbosity   level of information displayed
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_output_param(cs_domain_t  *domain,
                           int           restart_nt,
                           int           log_nt,
                           int           verbosity)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The domain structure is not allocated.", __func__);

  domain->restart_nt = restart_nt;
  domain->output_nt = log_nt;
  domain->verbosity = verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the computation stage in the domain structure
 *
 * \param[in, out] domain    pointer to a cs_domain_t structure
 * \param[in]      stage     stage in the computation run
 */
/*----------------------------------------------------------------------------*/

void
cs_domain_set_stage(cs_domain_t         *domain,
                    cs_domain_stage_t    stage)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The domain structure is not allocated.", __func__);

  domain->stage = stage;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the computation stage from the domain structure
 *
 * \param[in] domain    pointer to a cs_domain_t structure
 *
 * \return the current stage in the computation run
 */
/*----------------------------------------------------------------------------*/

cs_domain_stage_t
cs_domain_get_stage(const cs_domain_t    *domain)
{
  if (domain == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The domain structure is not allocated.", __func__);

  return domain->stage;
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

  cs_coupling_sync_apps(0,      /* flags */
                        ts->nt_cur,
                        &(ts->nt_max),
                        &(ts->dt_ref));

  if (ts->nt_max > 0) /* nt_max has been set */
    if (ts->nt_cur >= ts->nt_max)
      one_more_iter = false;

  if (ts->t_max > 0) /* t_max has been set */
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
 * \brief Check if an output is requested according to the domain setting
 *
 * \param[in] domain  pointer to a cs_domain_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_domain_needs_log(const cs_domain_t  *domain)
{
  if (domain->verbosity < 0)
    return false;

  if (domain->only_steady)
    return true;

  return cs_log_default_is_active();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the simulated time after one temporal iteration
 *
 * \param[in, out] domain     pointer to a cs_domain_t structure
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

END_C_DECLS
