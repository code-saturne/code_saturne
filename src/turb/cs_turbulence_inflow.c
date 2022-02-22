/*============================================================================
 * Turbulent inflow generation
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_equation_param.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_turbulence_inflow.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mass injection for turbulent quantities based
 *        on k and epsilon values.
 *
 * \param[in]  zone_name  name of zone to which injection should be added
 * \param[in]  k          turbulent kinetic energy
 * \param[in]  eps        turbulent dissipation
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_inflow_volume_mass_injection_k_eps(const char  *zone_name,
                                                 double       k,
                                                 double       eps)
{
  cs_turb_model_type_t  iturb = cs_glob_turb_model->iturb;
  int                   itytur = cs_glob_turb_model->itytur;

  if (itytur == 2) {

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(k)), zone_name, &k);

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(eps)), zone_name, &eps);

  }
  else if (itytur == 3) {
    cs_real_t val[6] = {2./3.*k, 2./3.*k, 2./3.*k, 0, 0, 0};

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(rij)), zone_name, val);

  }
  else if (iturb == CS_TURB_V2F_PHI) {

    double twothirds = 2./3.;

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(k)), zone_name, &k);

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(eps)), zone_name, &eps);

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(phi)), zone_name, &twothirds);

    /* There is no mass source term in the equation for f_bar */

  }
  else if (iturb == CS_TURB_K_OMEGA) {

    double omega_in = eps / cs_turb_cmu / k;

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(k)), zone_name, &k);

    cs_equation_add_volume_mass_injection_by_value
      (cs_field_get_equation_param(CS_F_(omg)), zone_name, &omega_in);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define mass injection for turbulent quantities based
 *        on a hydraulic diameter and reference velocity.
 *
 * \param[in]  zone_name  name of zone to which injection should be added
 * \param[in]  uref2      square of the reference flow velocity
 * \param[in]  dh         hydraulic diameter \f$ D_H \f$
 * \param[in]  rho        mass density \f$ \rho \f$
 * \param[in]  mu         dynamic viscosity \f$ \nu \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_turbulence_inflow_volume_mass_injection_ke_hyd_diam(const char  *zone_name,
                                                       double       uref2,
                                                       double       dh,
                                                       double       rho,
                                                       double       mu)
{
  cs_real_t ustar2 = 0, k = cs_math_epzero, eps = cs_math_epzero;

  /* Turbulence values */

  cs_turbulence_bc_ke_hyd_diam(uref2,
                               dh,
                               rho,
                               mu,
                               &ustar2,
                               &k,
                               &eps);

  cs_turbulence_inflow_volume_mass_injection_k_eps(zone_name,
                                                   k,
                                                   eps);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
