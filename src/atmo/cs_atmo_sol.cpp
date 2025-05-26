/*============================================================================
 * Atmospheric soil module - Build constants and variables to describe ground model
 *============================================================================*/

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_physical_constants.h"
#include "atmo/cs_atmo.h"
#include "base/cs_field.h"
#include "mesh/cs_mesh.h"
#include "atmo/cs_air_props.h"
#include "base/cs_math.h"
#include "alge/cs_divergence.h"
#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "atmo/cs_atmo_solcat.h"
#include "atmo/cs_atmo_soliva.h"
#include "atmo/cs_atmo_solmoy.h"
#include "base/cs_parall.h"


/*----------------------------------------------------------------------------
*  Header for the current file
*----------------------------------------------------------------------------*/

#include "atmo/cs_atmo_sol.h"

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build constants and variables to describe ground model
 */
/*----------------------------------------------------------------------------*/


void
cs_atmsol(void)
{
  /* Local variables */
  int error;
  cs_lnum_t n_elts;
  int n_soil_cat;
  const cs_lnum_t *elt_ids;
  cs_f_atmo_get_soil_zone(&n_elts, &n_soil_cat, &elt_ids);

  /* Get the number of elements in the soil zone */
  const cs_atmo_option_t *at_opt = cs_glob_atmo_option;

  cs_parall_sum(1, CS_INT_TYPE,&n_elts);

  /* There are some soil faces on some ranks */
  /* Note: we can use soil categories without soil model */
  /* (which solve temperature and humidity) */
  if (n_elts > 0) {
    /* Second pass, print and check soil categories parameters */
    cs_atmo_soil_cat(2);
    cs_solmoy(&error);
    if (error != 0) {
      bft_printf("Allocation error of atmodsol::solmoy\n");
      cs_exit(1);
    }

    /* Initialization of soil variables */
    /* Only if soil is activated */
    if (at_opt->soil_model >= 0) {
      cs_soliva();
    }
  } /* End of second call */
}

/*----------------------------------------------------------------------------*/
