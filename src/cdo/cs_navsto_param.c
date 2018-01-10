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
 * Local variables
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char
cs_navsto_param_model_name[CS_NAVSTO_PARAM_N_MODELS][CS_BASE_STRING_LEN] =
  { N_("Stokes velocity-pressure system"),
    N_("Navier-Stokes velocity-pressure system")
  };

static const char
cs_navsto_param_coupling_name[CS_NAVSTO_PARAM_N_COUPLINGS][CS_BASE_STRING_LEN] =
  { N_("Full system coupling"),
    N_("Artificial compressibility algorithm"),
    N_("Incremental projection algorithm")
  };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes system
 *
 * \param[in]  model     model related to the system to solve
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(cs_navsto_param_model_t        model)
{
  cs_navsto_param_t  *param = NULL;
  BFT_MALLOC(param, 1, cs_navsto_param_t);

  param->model = model;
  param->coupling = CS_NAVSTO_PARAM_COUPLING_PROJECTION;

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
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  param    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *param)
{
  if (param == NULL)
    return;

  /* Sanity checks */
  if (param->model == CS_NAVSTO_PARAM_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);
  if (param->coupling == CS_NAVSTO_PARAM_N_COUPLINGS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid way of coupling the Navier-Stokes equations.\n",
              __func__);

  cs_log_printf(CS_LOG_SETUP, " <NavSto/Model> %s\n",
                cs_navsto_param_model_name[param->model]);
  cs_log_printf(CS_LOG_SETUP, " <NavSto/Coupling> %s\n",
                cs_navsto_param_coupling_name[param->coupling]);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
