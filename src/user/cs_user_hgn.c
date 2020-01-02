/*============================================================================
 * User properties for homogeneous two-phase compressible model
 *============================================================================*/

/* VERS */

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_hgn.c
 *
 * \brief Define user properties for two-phase homogeneous compressible model.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computation of the relaxation time-scale.
 *
 * This function computes the value of the relaxation time-scale (for the
 * return to equilibrium).
 *
 * \param[in]  mesh       pointer to mesh
 * \param[in]  alpha_eq   equilibrium volume fraction
 * \param[in]  y_eq       equilibrium mass fraction
 * \param[in]  z_eq       equilibrium energy fraction
 * \param[in]  ei         specific internal energy
 * \param[in]  v          specific volume
 * \param[in]  relax_tau  relaxation time scale towards equilibrium
 */
/*----------------------------------------------------------------------------*/

void
cs_user_hgn_thermo_relax_time(const cs_mesh_t *mesh,
                              const cs_real_t *alpha_eq,
                              const cs_real_t *y_eq,
                              const cs_real_t *z_eq,
                              const cs_real_t *ei,
                              const cs_real_t *v,
                              cs_real_t       *relax_tau)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
