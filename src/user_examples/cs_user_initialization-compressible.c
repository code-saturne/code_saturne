/*============================================================================
 * User initialization prior to solving time steps.
 *============================================================================*/

/* VERS */

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
#include <math.h>

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
 * \file cs_user_initialization.c
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain)
{
  /*! [init_compressible] */
  /* If this is restarted computation, do not reinitialize values */
  if (domain->time_step->nt_prev > 0)
    return;

  const cs_lnum_t n_cells = domain->mesh->n_cells;

  cs_real_t *cvar_pr = CS_F_(p)->val;
  cs_real_t *cpro_rho = CS_F_(rho)->val;
  cs_real_t *cvar_energ = CS_F_(e_tot)->val;
  cs_real_t *cvar_tempk = CS_F_(t_kelvin)->val;
  cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

  const cs_real_t p0 = cs_glob_fluid_properties->p0;
  const cs_real_t t0 = cs_glob_fluid_properties->t0;
  const cs_real_t ro0 = cs_glob_fluid_properties->ro0;
  const cs_real_t cv0 = cs_glob_fluid_properties->cv0;

  /* Velocity all components */
  cs_array_set_value_real(n_cells, 3, 0., (cs_real_t*)cvar_vel);

  /* User defined scalars */
  cs_real_t * cvar_scal = cs_field_by_name("user_name")->val;

  /* Initialize each cell value */
  cs_array_set_value_real(n_cells, 1, 0., cvar_scal);

  /* Pressure, Density, Temperature, Total Energy

   * Only 2 out of these 4 variables are independent: one may choose to
   * initialize any pair of variables picked out of these 4, except
   * (Temperature-Energy). The remaining 2 variables will be deduced
   * automatically.

   * Initialize 2 and only 2 variables

   *   To do so, set iutile=1 for each of the 2 selected variables
   *             and iutile=0 for each of the 2 others

   *   In the example provided below, Pressure and Temperature are
   *   initialized.


   * ithvar indicates which variables have been set:
   *   it is completed automatically for each variable and
   *   it must not be modified. */

  cs_cf_model_t *cf_model = cs_get_glob_cf_model();
  const int ithvar = cf_model->ithvar;

  /* 1. Pressure (Pa) */
  if (true) {
    cf_model->ithvar = 2*ithvar;
    cs_array_set_value_real(n_cells, 1, p0, cvar_pr);
  }

  /*  2. Density (kg.m-3) */
  if (false) {
    cf_model->ithvar = 3*ithvar;
    cs_array_set_value_real(n_cells, 1, ro0, cpro_rho);
  }

  /* 3. Temperature (K -- Warning: Kelvin) */
  if (true) {
    cf_model->ithvar = 5*ithvar;
    cs_array_set_value_real(n_cells, 1, t0, cvar_tempk);
  }

  /* 4. Total Energy (J/kg) */
  if (false) {
    cf_model->ithvar = 7*ithvar;
    cs_array_set_value_real(n_cells, 1, cv0*t0, cvar_energ);
  }
  /*! [init_compressible] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
