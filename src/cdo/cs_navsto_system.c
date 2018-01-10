/*============================================================================
 * Routines to handle cs_navsto_system_t structure
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

#include "cs_equation.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_system.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_SYSTEM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

static cs_navsto_system_t  *cs_navsto_system = NULL;

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

struct _cs_navsto_system_t {

  /* Set of parameters to handle the Navier-Stokes system */
  cs_navsto_param_t  *param;

  /* Set of equations to solve (rely on the model choice, i.e. all equations are
     not necessary solved) */
  cs_equation_t      *momentum; /* Momentum balance equation (vector-valued) */
  cs_equation_t      *mass;     /* Mass balance equation (scalar-valued) */
  cs_equation_t      *energy;   /* Energy balance equation (scalar-valued) */

  /* Set of fields (resolved variables): fields are created according to the
     choice of model for Navier-Stokes */
  cs_field_t         *velocity;
  cs_field_t         *pressure;
  cs_field_t         *temperature;

  /* Set of properties: properties and their related fields are allocated
     according to the choice of model for Navier-Stokes */
  cs_property_t      *density;
  cs_property_t      *viscosity;

  /* Additional structures allocated if needed i.e. according to the choice
     of model for the Navier-Stokes system */
  void               *context;

  /* Pointer of functions related to the main stages */

  /* TODO: init_setup, finalize_setup, build, solve */

};

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the resolution of the Navier-Stokes system has been
 *        activated
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_navsto_system_is_activated(void)
{
  if (cs_navsto_system == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize the Navier-Stokes system
 *
 * \param[in]  model     type of model related to the Navier-Stokes system
 *
 * \return a pointer to a new allocated groundwater flow structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_system_t *
cs_navsto_system_activate(cs_navsto_param_model_t        model)
{
  cs_navsto_system_t  *navsto = NULL;

  /* Sanity checks */
  if (model == CS_NAVSTO_PARAM_N_MODELS)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);

  BFT_MALLOC(navsto, 1, cs_navsto_system_t);

  navsto->param = cs_navsto_param_create(model);

  /* Main set of equations */
  navsto->momentum = NULL;
  navsto->mass = NULL;
  navsto->energy = NULL; /* Up to now this equation is not handled */

  /* Main set of variables */
  navsto->velocity = NULL;
  navsto->pressure = NULL;
  navsto->temperature = NULL;

  /* Main set of properties */
  navsto->density = NULL;
  navsto->viscosity = NULL;

  /* Additional data fitting the choice of model */
  navsto->context = NULL;

  /* Set the static variable */
  cs_navsto_system = navsto;

  return navsto;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the main structure related to the Navier-Stokes system
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_destroy(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return;

  /* Free the context according to the model choice */
  assert(navsto->context == NULL);

  /* Set of numerical parameters */
  navsto->param = cs_navsto_param_free(navsto->param);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_system_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_system_log_setup(void)
{
  cs_navsto_system_t  *navsto = cs_navsto_system;

  if (navsto == NULL)
    return;

  /* Main set of numerical parameters */
  cs_navsto_param_log(navsto->param);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
