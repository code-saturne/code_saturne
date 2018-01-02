/*============================================================================
 * Caller interaction with Lagrangian module.
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

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_parall.h"

#include "cs_lagr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_query.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local types and structures
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return Lagranian module status.
 *
 * \return 0 if module is not active, > 0 if active
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_model_type(void)
{
  int retval = 0;

  if (cs_glob_lagr_time_scheme != NULL)
    retval = cs_glob_lagr_time_scheme->iilagr;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return Lagranian particle restart status.
 *
 * \return 1 if particles restart is available, 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_lagr_particle_restart(void)
{
  int retval = 0;

  if (cs_glob_lagr_time_scheme != NULL)
    retval = cs_glob_lagr_time_scheme->isuila;

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
