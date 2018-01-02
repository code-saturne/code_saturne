/*============================================================================
 * Physical properties management for groundwater flow module.
 *============================================================================*/

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


/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_gwf_parameters.h"
#include "cs_math.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"


/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gwf_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update physical properties of the Ground Water Flow module.
 *
 * Species transport is delayed by retention in solid phase.
 * This delay is computed as follows:
 * \f[
 * R = 1 + \rho \dfrac{K_d}{\theta} ;
 * \f]
 * where \f$ R \f$ is the delay factor, \f$ \rho \f$ the soil (bulk) density,
 * \f$ K_d \f$ the contaminant distribution coefficient and \f$ \theta \f$
 * the moisture content.
 *
 * Please refer to the
 * <a href="../../theory.pdf#gwf_sp_transp"><b>dedicated section</b></a>
 * of the theory guide for more informations.
 *
 */
 /*---------------------------------------------------------------------------*/

void cs_gwf_physical_properties(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_field_t *sca, *delay, *sat, *rosoil, *kd;
  cs_gwf_soilwater_partition_t sorption_scal;

  sat = cs_field_by_name("saturation");
  rosoil = cs_field_by_name("soil_density");

  int key_sorbed_c_id = cs_field_key_id("gwf_sorbed_concentration_id");
  int key_kinet_soilwater_partition_id =
    cs_field_key_id("gwf_soilwater_partition");

  /* Loop on user scalars */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
     sca = cs_field_by_id(f_id);
     if ((sca->type & CS_FIELD_VARIABLE) && (sca->type & CS_FIELD_USER)) {
      /* Sorption properties are needed to update delay */
      cs_field_get_key_struct(sca,
                              key_kinet_soilwater_partition_id,
                              &sorption_scal);

      /* Update of delay */
      kd = cs_field_by_id(sorption_scal.ikd);
      delay = cs_field_by_id(sorption_scal.idel);
      for (int iel = 0; iel < n_cells; iel++) {
        delay->val[iel] = 1. + rosoil->val[iel] * kd->val[iel] / sat->val[iel];
      }

      /* Update of sorbed concentration
         only if kinetic model is enabled for the scalar */
      if (sorption_scal.kinetic == 1) {
        int sorbed_c_id = cs_field_get_key_int(sca, key_sorbed_c_id);
        cs_real_t *sorb = cs_field_by_id(sorbed_c_id)->val;
        cs_real_t *kp = cs_field_by_id(sorption_scal.ikp)->val;
        cs_real_t *km = cs_field_by_id(sorption_scal.ikm)->val;

        cs_gwf_sorbed_concentration_update(sca->val, kp, km, sorb);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update sorbed concentration for scalars with kinetic sorption.
 *
 * It is estimated by the following analytical expression :
 * \f[
 * S^{n+1} = S^n \exp(- k^{-} \Delta t) - C^n \dfrac{k^{+}}{k^{-}}
 * \left(\exp(- k^{-} \Delta t) - 1 \right)
 * \f]
 *
 * \param[in]   c_scal      concentration field
 * \param[in]   kp          kplus field
 * \param[in]   km          kminus field
 * \param[in]   sorb        sorbed concentration field
 */
 /*----------------------------------------------------------------------------*/

void cs_gwf_sorbed_concentration_update(cs_real_t *c_scal,
                                        cs_real_t *kp,
                                        cs_real_t *km,
                                        cs_real_t *sorb)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t *dt = CS_F_(dt)->val;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    if (km[c_id] > cs_math_epzero) {
      cs_real_t expkdt = exp(-km[c_id] * dt[c_id]);
      cs_real_t kpskm = kp[c_id] / km[c_id];
      sorb[c_id] =  expkdt * sorb[c_id]
                  - (expkdt-1.) * kpskm * c_scal[c_id];
    }
    else
      sorb[c_id] =  sorb[c_id] + dt[c_id]*kp[c_id]*c_scal[c_id];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
