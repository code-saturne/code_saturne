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
#include "cs_mesh_quantities.h"

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
 * \brief  Update delay of the Groundwater Flow module.
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
void cs_gwf_delay_update(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_field_t *sca, *delay, *sat, *rosoil, *kd;
  cs_gwf_soilwater_partition_t sorption_scal;
  int key_part = cs_field_key_id("gwf_soilwater_partition");

  /* Get soil property fields */
  sat = cs_field_by_name("saturation");
  rosoil = cs_field_by_name("soil_density");

  /* Loop on user scalars */
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    sca = cs_field_by_id(f_id);
    if ((sca->type & CS_FIELD_VARIABLE) && (sca->type & CS_FIELD_USER)) {
      /* Sorption properties are needed to update delay */
      cs_field_get_key_struct(sca, key_part, &sorption_scal);

      /* Update of delay */
      kd = cs_field_by_id(sorption_scal.ikd);
      delay = cs_field_by_id(sorption_scal.idel);
      for (int c_id = 0; c_id < n_cells; c_id++) {
        delay->val[c_id] = 1. + rosoil->val[c_id] * kd->val[c_id] / sat->val[c_id];
      }
    }
  }
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add first-order decay in left-hand member (implicit)
 *
 * dc/dt = -decay_rate * c
 */
 /*---------------------------------------------------------------------------*/
void cs_gwf_decay_rate(int        f_id,
                       cs_real_t *ts_imp)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  const cs_real_t *vol = mesh_quantities->cell_vol;
  cs_real_t decay_rate = 0.;
  cs_field_t *sca;

  sca = cs_field_by_id(f_id);
  decay_rate = cs_field_get_key_double(sca, cs_field_key_id("fo_decay_rate"));
  /* First-order decay_rate in implict term */
  for (int c_id = 0; c_id < n_cells; c_id++)
    ts_imp[c_id] -= decay_rate * vol[c_id];
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Update sorbed concentration for scalars with kinetic sorption.
 *
 * It is estimated by the following analytical expression :
 * \f[
 * S^{n+1} = S^n \exp(- (k^{-} + decay_rate) * \Delta t) - C^n \dfrac{k^{+}}{k^{-}}
 * \left(\exp(- (k^{-} + decay_rate) * \Delta t) - 1 \right)
 * \f]
 *
 * \param[in]   f_id      scalar index
 */
 /*---------------------------------------------------------------------------*/
void cs_gwf_sorbed_concentration_update(int f_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t *dt = CS_F_(dt)->val;
  cs_real_t decay_rate = 0.;
  cs_field_t *sca, *kd, *kp, *km, *sorb;
  cs_gwf_soilwater_partition_t sorption_scal;
  int key_part = cs_field_key_id("gwf_soilwater_partition");
  int key_sorb = cs_field_key_id("gwf_sorbed_concentration_id");

  /* Get scalar and sorbed concentration fields, and first-order decay rate */
  sca = cs_field_by_id(f_id);
  sorb = cs_field_by_id(cs_field_get_key_int(sca, key_sorb));
  decay_rate = cs_field_get_key_double(sca, cs_field_key_id("fo_decay_rate"));

  /* Soil-water partition structure */
  cs_field_get_key_struct(sca, key_part, &sorption_scal);
  kp = cs_field_by_id(sorption_scal.ikp);
  km = cs_field_by_id(sorption_scal.ikm);

  /* Update sorbed concentration */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    /* Case of reversible sorption or decay rate term */
    if (km->val[c_id] + decay_rate > cs_math_epzero) {
      cs_real_t expkdt = exp(-(km->val[c_id] + decay_rate) * dt[c_id]);
      cs_real_t kpskm = kp->val[c_id] / (km->val[c_id] + decay_rate);
      sorb->val[c_id] =  expkdt * sorb->val[c_id]
        - (expkdt-1.) * kpskm * sca->val[c_id];
    }
    else /* Irreversible sorption and no decay rate */
      sorb->val[c_id] =  sorb->val[c_id]
        + dt[c_id]*kp->val[c_id]*sca->val[c_id];
}
/*---------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Update liquid concentration according to precipitation phenomenon.
 *
 * If liquid concentration exceeds solubility index, then it is clipped and
 * transferred in precipitated concentration
 *
 * \param[in]   f_id      scalar index
 */
 /*---------------------------------------------------------------------------*/
void cs_gwf_precipitation(int f_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t ctot_tmp = 0.;
  cs_field_t *sca, *precip, *solub;
  cs_gwf_soilwater_partition_t sorption_scal;
  int key_part = cs_field_key_id("gwf_soilwater_partition");
  int key_pre = cs_field_key_id("gwf_precip_concentration_id");

  /* Get scalar, precipitation and solubility index fields */
  sca = cs_field_by_id(f_id);
  cs_field_get_key_struct(sca, key_part, &sorption_scal);
  precip = cs_field_by_id(cs_field_get_key_int(sca, key_pre));
  solub = cs_field_by_id(sorption_scal.imxsol);

  /* Clipping of solute concentration, update of precipitation concentration */
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    {
      ctot_tmp = sca->val[c_id] + precip->val[c_id];
      sca->val[c_id] = CS_MIN(ctot_tmp, solub->val[c_id]);
      precip->val[c_id] = CS_MAX(0., ctot_tmp - solub->val[c_id]);
    }
}
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief Impact of kinetic chemical reaction on evolution equation
 *        of total concentration.
 *
 *
 * S^{n+1} = S^n \exp(- (k^{-} + decay_rate) * \Delta t) - C^n \dfrac{k^{+}}{k^{-}}
 * \left(\exp(- (k^{-} + decay_rate) * \Delta t) - 1 \right)
 * \f]
 *
 * \param[in]   f_id      scalar index
 */
 /*---------------------------------------------------------------------------*/
void cs_gwf_kinetic_reaction(int f_id, cs_real_t *ts_imp, cs_real_t *ts_exp)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  const cs_real_t *vol = mesh_quantities->cell_vol;
  cs_real_t *dt = CS_F_(dt)->val;
  cs_real_t decay_rate = 0.;
  cs_real_t expkdt, kplskm, rokpl;
  cs_field_t *sca, *kd, *kp, *km, *sorb, *rosoil;
  cs_gwf_soilwater_partition_t sorption_scal;
  int key_part = cs_field_key_id("gwf_soilwater_partition");
  int key_sorb = cs_field_key_id("gwf_sorbed_concentration_id");

  /* Get soil density */
  rosoil = cs_field_by_name("soil_density");

  /* Get scalar and sorbed concentration fields, and first-order decay rate */
  sca = cs_field_by_id(f_id);
  sorb = cs_field_by_id(cs_field_get_key_int(sca, key_sorb));
  decay_rate = cs_field_get_key_double(sca, cs_field_key_id("fo_decay_rate"));

  /* Soil-water partition structure */
  cs_field_get_key_struct(sca, key_part, &sorption_scal);
  kp = cs_field_by_id(sorption_scal.ikp);
  km = cs_field_by_id(sorption_scal.ikm);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
  /* General case (reversible sorption) */
    if (km->val[c_id] > cs_math_epzero) {
      cs_real_t expkdt = exp(-(km->val[c_id] + decay_rate) * dt[c_id]);
      cs_real_t kpskm = kp->val[c_id] / (km->val[c_id] + decay_rate);
      ts_exp[c_id] += - vol[c_id] *
        (decay_rate * rosoil->val[c_id] * sorb->val[c_id]
          + rosoil->val[c_id]/dt[c_id] * (1-expkdt)
         *(kpskm*sca->val[c_id] - sorb->val[c_id]));
      ts_imp[c_id] += vol[c_id] * rosoil->val[c_id] / dt[c_id]
         * (1-expkdt)*kpskm;
    }
    else { /* Case of irreversible sorption without decay rate */
      rokpl = rosoil->val[c_id] * kp->val[c_id];
      ts_exp[c_id] += - vol[c_id] * rokpl * sca->val[c_id];
      ts_imp[c_id] += + vol[c_id] * rokpl;
    }
}

END_C_DECLS
