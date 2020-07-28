/*============================================================================
 * Lagrangian model options.
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

#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_lagr.h"
#include "cs_lagr_post.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_particle.h"
#include "cs_lagr_prototypes.h"
#include "cs_prototypes.h"

/*---------------------------------------------------------------------------*/
/*
 * \brief User function of the Lagrangian particle-tracking module
 *
 *  User input of physical, numerical and post-processing options.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_lagr_model(void)
{

  /* ==========================================================================
   * 4. Calculation features for the dispersed phases
   * ========================================================================== */

  /* 4.4 Volume statistics
     --------------------- */

  /* Threshold for the management of volume statistics
     -------------------------------------------------
     * the value of the seuil variable is a statistical weight.
     * each cell of the mesh contains a statistical weight
     (sum of the statistical weights of all the particles
     located in the cell); seuil is the minimal value from
     which the contribution in statistical weight of a particle
     is not taken into account anymore in the full model
     of turbulent dispersion, in the resolution of the
     Poisson equation of correction of the mean velocities, and
     in the writing of the listing and post-processing. */

  cs_glob_lagr_stat_options->threshold = 0.0;

  cs_glob_lagr_stat_options->nstist = 120000;

  /* 4.4.2 Volume statistical variables  */
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* Activation of the calculation of the particle volume fraction     */
  /* Name of the mean : Part_vol_frac    */

  cs_lagr_stat_activate(CS_LAGR_STAT_VOLUME_FRACTION);

  /* Activation of the calculation of the particle velocity */

  cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY);

  cs_lagr_stat_activate_attr(CS_LAGR_VELOCITY_SEEN);

  /* Activation of the calculation of the particle residence time */

  cs_lagr_stat_activate_attr(CS_LAGR_RESIDENCE_TIME);

  /* Activation of the calculation of the weight */

  cs_lagr_stat_activate_attr(CS_LAGR_STAT_WEIGHT);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
