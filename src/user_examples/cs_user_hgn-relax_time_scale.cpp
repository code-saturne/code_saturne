/*============================================================================
 * User properties for homogeneous two-phase compressible model
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Compute the relaxation time-scale to equilibrium in the
 * frame of the homogeneous two-phase model..
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
cs_user_hgn_thermo_relax_time([[maybe_unused]] const cs_mesh_t  *mesh,
                              [[maybe_unused]] const cs_real_t  *alpha_eq,
                              [[maybe_unused]] const cs_real_t  *y_eq,
                              [[maybe_unused]] const cs_real_t  *z_eq,
                              [[maybe_unused]] const cs_real_t  *ei,
                              [[maybe_unused]] const cs_real_t  *v,
                              [[maybe_unused]] cs_real_t        *relax_tau)
{
  //-----------------
  // pour run 1
  //-----------------

  /*alpha_c=0.35;
    lbd0=4.0e-4;
    delta_sat=cs::abs(alpha-alpha_eq);
    epsilon_a=alpha/alpha_c;
    tps_relax=lbd0   *   exp(-delta_sat*delta_sat/1.0e-3)   *   exp(-epsilon_a*epsilon_a*epsilon_a);*/

  //-----------------
  // pour run 2
  //-----------------

  /*alpha_c=0.42;
    lbd0=3.0e-5;
    delta_sat=cs::abs(alpha-alpha_eq);
    epsilon_a=alpha/alpha_c;
    tps_relax=lbd0   *   exp(-delta_sat*delta_sat/5.e-4)  *   exp(-epsilon_a*epsilon_a);*/
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
