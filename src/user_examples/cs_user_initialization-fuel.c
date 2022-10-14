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
 * \file cs_user_initialization-fuel.c
 *
 * \brief Initialization prior to solving time steps.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
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
  /*![init]*/
  /* If this is restarted computation, do not reinitialize values */
  if (domain->time_step->nt_prev > 0)
    return;

  const cs_real_t n_cells = domain->mesh->n_cells;

  const cs_combustion_model_t  *cm = cs_glob_combustion_model;
  const int io2 = cm->io2 - 1;
  const int in2 = cm->in2 - 1;
  const int ico2 = cm->ico2 - 1;
  const int ih2o = cm->ih2o - 1;

  const cs_real_t *wmole = cm->wmole;
  const cs_real_t *oxyo2 = cm->oxyo2;
  const cs_real_t *oxyn2 = cm->oxyn2;
  const cs_real_t *oxyh2o = cm->oxyh2o;
  const cs_real_t *oxyco2 = cm->oxyco2;

  /* The domain is filled with air at t1init;

     Computation of h1init */

  cs_real_t t1init = 1000.0;
  cs_real_t t2init = 1000.0;

  /* Transported variables for droplets */

  const cs_real_t trefth = 25.0 + cs_physical_constants_celsius_to_kelvin;

  cs_real_t h2init = cm->fuel.h02fol + cm->fuel.cp2fol*(t2init-trefth);

  /* Transported variables for the mix (droplets and carrying gases) */

  cs_real_t coefe[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  for (int ige = 0; ige < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; ige++)
    coefe[ige] = 0.;

  /* Consider oxydant 1 */

  coefe[io2] =    wmole[io2]*oxyo2[1]
               / (  wmole[io2]*oxyo2[1]  + wmole[in2]*oxyn2[1]
                  + wmole[ih2o]*oxyh2o[1] + wmole[ico2]*oxyco2[1]);
  coefe[ih2o] =   wmole[ih2o]*oxyh2o[1]
                / (  wmole[io2]*oxyo2[1] + wmole[in2] * oxyn2[1]
                   + wmole[ih2o]*oxyh2o[1] + wmole[ico2]*oxyco2[1]);
  coefe[ico2] =   wmole[ico2]*oxyco2[1]
                / (  wmole[io2]*oxyo2[1] +wmole[in2]*oxyn2[1]
                   + wmole[ih2o]*oxyh2o[1]+wmole[ico2]*oxyco2[1]);
  coefe[in2] = 1. - coefe[io2] - coefe[ih2o] - coefe[ico2];

  cs_real_t h1init = cs_fuel_t2h_gas(coefe, t1init);

  cs_array_set_value_real(n_cells, 1, h1init, CS_F_(h)->val);

  /* Transported variables for the mix (passive scalars, variance)
   * Variables not present here are initialized to 0. */

  if (CS_F_(yco2) != NULL) {

    /* Consider oxydant 1 */
    const int ioxy = 0;
    const cs_real_t wmo2 = wmole[io2];
    const cs_real_t wmco2 = wmole[ico2];
    const cs_real_t wmh2o = wmole[ih2o];
    const cs_real_t wmn2 = wmole[in2];
    const cs_real_t dmas = ( oxyo2 [ioxy]*wmo2  + oxyn2[ioxy]*wmn2
                           + oxyh2o[ioxy]*wmh2o + oxyco2[ioxy]*wmco2);

    const cs_real_t co2 = oxyco2[ioxy]*wmco2/dmas;
    cs_array_set_value_real(n_cells, 1, co2, CS_F_(yco2)->val);

  }
  /*![init]*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
