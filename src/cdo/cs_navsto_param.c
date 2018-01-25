/*============================================================================
 * Routines to handle cs_navsto_param_t structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_PARAM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char
cs_navsto_param_model_name[CS_NAVSTO_N_MODELS][CS_BASE_STRING_LEN] =
  { N_("Stokes velocity-pressure system"),
    N_("Incompressible Navier-Stokes velocity-pressure system")
  };

static const char
cs_navsto_param_time_state_name[CS_NAVSTO_N_TIME_STATES][CS_BASE_STRING_LEN] =
  { N_("Fully unsteady"),
    N_("Fully steady"),
    N_("Steaty-state as the limit of an unsteady computation")
  };

static const char
cs_navsto_param_coupling_name[CS_NAVSTO_N_COUPLINGS][CS_BASE_STRING_LEN] =
  { N_("Uzawa-Augmented Lagrangian coupling"),
    N_("Artificial compressibility algorithm"),
    N_("Incremental projection algorithm")
  };

static const char _err_empty_nsp[] =
  N_(" Stop setting an empty cs_navsto_param_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in]  model          model related to the NS system to solve
 * \param[in]  time_state     state of the time for the NS equations
 * \param[in]  algo_coupling  algorithm used for solving the NS system
*
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(cs_navsto_param_model_t        model,
                       cs_navsto_param_time_state_t   time_state,
                       cs_navsto_param_coupling_t     algo_coupling)
{
  cs_navsto_param_t  *param = NULL;
  BFT_MALLOC(param, 1, cs_navsto_param_t);

  param->verbosity = 1;

  /* Which equations are solved and which terms are needed */
  param->model = model;
  param->has_gravity = false;
  param->gravity[0] = param->gravity[1] = param->gravity[2] = 0.;
  param->time_state = time_state;
  param->coupling = algo_coupling;
  param->ac_zeta_coef = 1.0;    /* Default value if not set by the user */

  param->space_scheme = CS_SPACE_SCHEME_CDOFB;

  return param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param)
{
  if (param == NULL)
    return param;

  BFT_FREE(param);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a cs_navsto_param_t
 *         structure
 *
 * \param[in, out] nsp      pointer to a cs_navsto_param_t structure to set
 * \param[in]      key      key related to the member of eq to set
 * \param[in]      keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set(cs_navsto_param_t    *nsp,
                    cs_navsto_key_t       key,
                    const char           *keyval)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_nsp);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_NSKEY_AC_ZETA_COEF:
    nsp->ac_zeta_coef = atof(val);
    break;

  case CS_NSKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_fb") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_CDOFB;
    }
    else if (strcmp(val, "hho_p0") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P0;
    }
    else if (strcmp(val, "hho_p1") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P1;
    }
    else if (strcmp(val, "hho_p2") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P2;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_SPACE_SCHEME\n"
                  " Choice between hho_{p0, p1, p2} or cdo_fb"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_VERBOSITY:
    nsp->verbosity = atoi(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid key for setting the Navier-Stokes system."),
              __func__);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    return;

  /* Sanity checks */
  if (nsp->model == CS_NAVSTO_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);
  if (nsp->coupling == CS_NAVSTO_N_COUPLINGS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid way of coupling the Navier-Stokes equations.\n",
              __func__);

  cs_log_printf(CS_LOG_SETUP, " <NavSto/Verbosity> %d\n", nsp->verbosity);

  const char *space_scheme = cs_param_get_space_scheme_name(nsp->space_scheme);
  if (nsp->space_scheme != CS_SPACE_N_SCHEMES)
    cs_log_printf(CS_LOG_SETUP, " <NavSto/Space scheme> %s\n", space_scheme);
  else
    bft_error(__FILE__, __LINE__, 0,
              " Undefined space scheme for the Navier-Stokes system");

  cs_log_printf(CS_LOG_SETUP, " <NavSto/Model> %s\n",
                cs_navsto_param_model_name[nsp->model]);
  cs_log_printf(CS_LOG_SETUP, " <NavSto/Time status> %s\n",
                cs_navsto_param_time_state_name[nsp->time_state]);
  cs_log_printf(CS_LOG_SETUP, " <NavSto/Coupling> %s\n",
                cs_navsto_param_coupling_name[nsp->coupling]);
  cs_log_printf(CS_LOG_SETUP, " <NavSto/Gravity effect> %s",
                cs_base_strtf(nsp->has_gravity));
  if (nsp->has_gravity)
    cs_log_printf(CS_LOG_SETUP, " vector: [% 5.3e; % 5.3e; % 5.3e]\n",
                  nsp->gravity[0], nsp->gravity[1], nsp->gravity[2]);
  else
    cs_log_printf(CS_LOG_SETUP, "\n");
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
