/*============================================================================
 * User initialization prior to solving time steps.
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
 * \file cs_user_initialization-pulverized_coal.c
 *
 * \brief Initialization prior to solving time steps.
 *        Pulverized coal example.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
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
  /*! [loc_var_dec] */
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_fluid_properties_t *fprops = cs_glob_fluid_properties;

  cs_coal_model_t *cm = cs_glob_coal_model;

  cs_real_t coefe[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
  /*! [loc_var_dec] */

  /*! [init] */
  /* If this is restarted computation, do not reinitialize values */

  if (domain->time_step->nt_prev > 0)
    return;

  /* Control Print */

  bft_printf("%s: settings for pulverized coal\n", __func__);

  /* All the domain is filled with the first oxidizer at TinitK
     ---------------------------------------------------------- */

  /* Transported variables for the mix (solid+carrying gas)^2 */

  const int ico2 = cm->ico2 - 1;
  const int ih2o = cm->ih2o - 1;
  const int in2  = cm->in2 - 1;
  const int io2  = cm->io2 - 1;

  for (int ige = 0; ige < CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS; ige++)
    coefe[ige] = 0.;

  /* Oxidizers are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
     the composition of the fisrt oxidiser is taken in account. */

  const int ioxy = 0;

  cs_real_t dmas =   cm->wmole[io2]  * cm->oxyo2[ioxy]
                   + cm->wmole[in2]  * cm->oxyn2[ioxy]
                   + cm->wmole[ih2o] * cm->oxyh2o[ioxy]
                   + cm->wmole[ico2] * cm->oxyco2[ioxy];

  coefe[io2]  = cm->wmole[io2]  * cm->oxyo2[ioxy ]/dmas;
  coefe[ih2o] = cm->wmole[ih2o] * cm->oxyh2o[ioxy]/dmas;
  coefe[ico2] = cm->wmole[ico2] * cm->oxyco2[ioxy]/dmas;
  coefe[in2]  = cm->wmole[in2]  * cm->oxyn2[ioxy ]/dmas;

  /* Computation of h1init and t1init */

  cs_real_t t1init = fprops->t0;
  cs_real_t h1init
    = cs_coal_ht_convert_h_to_t_gas_by_yi(t1init, coefe);

  cs_real_t *cvar_h = CS_F_(h)->val;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
    cvar_h[cell_id] = h1init;

  /* Transported variables for the mix (passive scalars, variance).
     Variables not present here are initialized to 0. */

  if (cm->ieqco2 >= 1) {
    cs_real_t *cvar_yco2 = cs_field_by_name("x_c_co2")->val;

    cs_real_t xco2 = cm->oxyco2[ioxy] * cm->wmole[ico2] / dmas;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      cvar_yco2[cell_id] = xco2;
  }

  if (cm->ieqnox == 1) {
    cs_real_t *cvar_nox = cs_field_by_name("x_c_h_ox")->val;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++)
      cvar_nox[cell_id] = h1init;
  }
  /*! [init] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
