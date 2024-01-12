/*============================================================================
 * User functions for input of calculation parameters.
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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-ctwr.c
 *
 * \brief Cooling towers parameters example
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Activate cooling tower model.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* Note: packing zones are defined in cs_user_zones.c. */

  /*! [ctwr_user_model_1] */
  {
    /* Activate cooling tower model */
    cs_glob_physical_model_flag[CS_COOLING_TOWERS] = 1;

    /* Evaporation model */
    cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
    ct_opt->evap_model = CS_CTWR_POPPE;
  }
  /*! [ctwr_user_model_1] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* Activate compressibility */
  {
    cs_velocity_pressure_model_t *vp_model
      = cs_get_glob_velocity_pressure_model();
    vp_model->idilat = 2;
  }

  /* Authorize variable density */
  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
    fp->irovar = 1;
  }

  /* Define humid air properties */
  {
    cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
    // Used to compute the humid air density as a function of (P,T,humidity)
    fp->ro0 = 1.2; //1.293

    // Humid air viscosity
    fp->viscl0 = 1.765e-05;

    cs_air_fluid_props_t *air_prop = cs_glob_air_props;
    // Dry air and water vapor properties
    air_prop->cp_a = 1006.0;
    air_prop->cp_v = 1831.0;

    // Initial absolute humidity
    air_prop->humidity0 = 5.626e-03;//34.5% relative humidity

    // Humid air conductivity - considered constant in the modelling
    air_prop->lambda_h = 2.493;

    // Liquid water properties
    air_prop->rho_l = 997.85615;
    air_prop->cp_l = 4179.0;
    air_prop->lambda_l = 0.02493;

    // Phase change properties
    air_prop->hv0 = 2501600.0;

    air_prop->droplet_diam = 0.005;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
