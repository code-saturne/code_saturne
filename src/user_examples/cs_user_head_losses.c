/*============================================================================
 * User head loss definitions.
 *============================================================================*/

/* VERS */

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_parameters.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_head_losses.c
 *
 * \brief User head loss definitions.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define head losses for a given volume zone.
 *
 * Head loss tensor coefficients for each cell are organized as follows:
 * ck11, ck22, ck33, ck12, ck23, ck13.
 *
 * Coefficients are set to zero (then computed based on definitions provided
 * through the GUI if this is the case) before calling this function, so
 * setting values to zero is usually not necessary, unless we want to fully
 * overwrite a GUI-based definition.
 *
 * Diagonal coefficients must be positive; the calculation may diverge
 * if this is not the case.
 *
 * \param[in]       zone  pointer to zone structure
 * \param[in, out]  cku   head loss coefficients
 */
/*----------------------------------------------------------------------------*/

void
cs_user_head_losses(const  cs_volume_zone_t  *zone,
                    cs_real_t                 cku[][6])
{
  /*! [map_field_arrays] */
  const cs_real_3_t *cvara_vel = (const cs_real_3_t *)(CS_F_(u)->val_pre);
  /*! [map_field_arrays] */

  /* Note that in the following examples, we check the zone name, so we
     know which zone we are dealing with using in case of multiple zones.
     The zone must have been defined either in the GUI or in
     \ref cs_user_zones.c. */

  /* Example: diagonal tensor for head losses in direction x */

  /*! [head_loss_1] */
  {
    if (strcmp(zone->name, "head_loss_1") == 0) {
      for (cs_lnum_t i = 0; i < zone->n_cells; i++) {
        cs_lnum_t c_id = zone->cell_ids[i];
        cs_real_t v = cs_math_3_norm(cvara_vel[c_id]);
        cku[i][0] = 10.0 * v;
        cku[i][1] = 0.0;
        cku[i][2] = 0.0;
      }
    }
  }
  /*! [head_loss_1] */

  /* Example: 3x3 tensor
   *
   * Example of head losses at alpha = 45 degres x,y
   * direction x resists by ck0 and y by ck1
   * ck1 = 0 represents vanes as follows: /////
   * in coordinate system x y

   *   Y|  /y
   *    | /
   *    |/
   *    \--------------- X
   *     \ / ALPHA
   *      \
   *       \ x
  */

  /*! [head_loss_2] */
  {
    if (strcmp(zone->name, "head_loss_1") == 0) {

      /* define rotation matrix outside of loop on cells */

      cs_real_t alpha  = cs_math_pi/4.0;
      cs_real_t cosa = cos(alpha);
      cs_real_t sina = sin(alpha);
      cs_real_t ck0 = 10.0;
      cs_real_t ck1 =  0.0;

      cs_real_t a11 = cs_math_sq(cosa)*ck0 + cs_math_sq(sina)*ck1;
      cs_real_t a22 = cs_math_sq(sina)*ck0 + cs_math_sq(cosa)*ck1;
      cs_real_t a12 = cosa * sina * (ck0 - ck1);

      /* compute local coefficients */

      for (cs_lnum_t i = 0; i < zone->n_cells; i++) {

        cs_lnum_t c_id = zone->cell_ids[i];
        cs_real_t v = cs_math_3_norm(cvara_vel[c_id]);

        cku[i][0] = a11 * v;
        cku[i][1] = a22 * v;
        cku[i][2] = 0.;
        cku[i][3] = a12 * v;
        cku[i][4] = 0.;
        cku[i][5] = 0.;

      }
    }
  }
  /*! [head_loss_2] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
