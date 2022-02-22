/*============================================================================
 * Gas combustion model.
 *============================================================================*/

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_combustion_model.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_combustion_gas.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_combustion_gas.c
        Gas combustion model.
*/

/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Compute molar and mass fractions of
 *         elementary species Ye, Xe (fuel, O2, CO2, H2O, N2) from
 *         global species Yg (fuel, oxidant, products)
 *
 * \param[in]     yg            global mass fraction
 * \param[out]    ye            elementary mass fraction
 * \param[out]    xe            elementary molar fraction
 */
/*----------------------------------------------------------------------------*/

void
cs_combustion_gas_yg2xye(const cs_real_t  yg[],
                         cs_real_t        ye[],
                         cs_real_t        xe[])
{
  const cs_combustion_model_t *cm = cs_glob_combustion_model;
  const int n_gas_e = cm->n_gas_el_comp;
  const int n_gas_g = cm->n_gas_species;

  /* Yg -> Ye conversion */

  for (int i = 0; i < n_gas_e; i++) {
    ye[i] = 0;
    for (int j = 0; j < n_gas_g; j++)
      ye[i] += cm->gas.coefeg[j][i] * yg[j];
  }

  /* Verification */

  cs_real_t ytot = 0;
  for (int i = 0; i < n_gas_e; i++)
    ytot += ye[i];

  if (ytot < 0 || (1-ytot) < -cs_math_epzero)
    bft_printf(_(" Warning:\n"
                 " --------\n"
                 "   %s; mass fraction sum outside [0, 1] bounds\n"
                 "   sum_1=1,%d Yi = %g\n\n"),
               __func__, n_gas_e, ytot);

  /* Molar mass mixture */

  cs_real_t mm = 0;
  for (int i = 0; i < n_gas_e; i++) {
    mm += ye[i] / cm->wmole[i];
  }
  mm = 1 / mm;

  /* Ye -> Xe conversion */
  for (int i = 0; i < n_gas_e; i++) {
    xe[i] = ye[i] * mm / cm->wmole[i];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
