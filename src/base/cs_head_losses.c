/*============================================================================
 * Head losses computation..
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field_pointer.h"
#include "cs_gui.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parall.h"
#include "cs_mesh_location.h"
#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_head_losses.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_head_losses.c
        Head losses computation.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Global variables
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
/*!
 * \brief Compute head loss coefficients.
 *
 * \param[out]  cku  head loss coefficients for all zones
 */
/*----------------------------------------------------------------------------*/

void
cs_head_losses_compute(cs_real_6_t cku[])
{
  const int n_zones = cs_volume_zone_n_zones();

  if (n_zones == 0)
    return;

  /* We first use an interleaved definition, then switch to the
     non-interleaved variant still used in Fortran */

  cs_lnum_t n_loc_cells = 0;
  cs_lnum_t n_hl_cells = 0;

  for (int i = 0; i < n_zones; i++) {
    const cs_zone_t  *z = cs_volume_zone_by_id(i);
    if (z->type & CS_VOLUME_ZONE_HEAD_LOSS) {
      n_hl_cells += z->n_elts;
      if (z->n_elts > n_loc_cells)
        n_loc_cells = z->n_elts;
    }
  }

  const cs_real_3_t *cvara_vel = (const cs_real_3_t *)(CS_F_(vel)->val_pre);

  /* Loop on head loss zones */

  cs_lnum_t n_p_cells = 0;

  for (int i = 0; i < n_zones; i++) {

    const cs_zone_t  *z = cs_volume_zone_by_id(i);
    if (z->type & CS_VOLUME_ZONE_HEAD_LOSS) {

      const cs_lnum_t n_z_cells = z->n_elts;
      cs_real_6_t *_cku = cku + n_p_cells;

      /* Initialize */

      for (cs_lnum_t j = 0; j < n_z_cells; j++) {
        for (cs_lnum_t k = 0; k < 6; k++)
          _cku[j][k] = 0.;
      }

      /* GUI definitions go first, then user function definitions */

      cs_gui_head_losses(z, cvara_vel, _cku);
      cs_user_head_losses(z, _cku);

      /* update previous cells accumulator */

      n_p_cells += n_z_cells;

    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
