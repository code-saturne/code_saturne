/*============================================================================
 * Physical properties management for groundwater flow module.
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
 * \brief  Update delay of all transported species (user scalars)
 *
 * Species transport is delayed by retention in solid phase.
 * This delay is computed as follows:
 * \f[
 *    R = 1 + \rho \dfrac{K_d}{\theta} ;
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
        delay->val[c_id] = 1. +  rosoil->val[c_id]
                               * kd->val[c_id] / sat->val[c_id];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add first-order decay to implicit part of source term array
 *
 * Corresponding EDO for decay phenomenon is:
 * \f[
 *    \dfrac{dc}{dt} = - decay_rate  c
 * \f]
 *
 * \param[in]     f_id      field index of scalar on which decay is set
 * \param[in,out] ts_imp    source term implicit part for scalar of index f_id
 */
 /*---------------------------------------------------------------------------*/

void cs_gwf_decay_rate(const int        f_id,
                       cs_real_t       *ts_imp)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  const cs_real_t *vol = mesh_quantities->cell_vol;

  cs_field_t *sca = cs_field_by_id(f_id);
  const int dr_id = cs_field_key_id("fo_decay_rate");
  cs_real_t decay_rate = cs_field_get_key_double(sca, dr_id);

  /* First-order decay_rate in implict term */
  for (int c_id = 0; c_id < n_cells; c_id++)
    ts_imp[c_id] -= decay_rate * vol[c_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update sorbed concentration for scalars with kinetic sorption.
 *
 * It is estimated by the following analytical expression :
 * \f[
 * S^{n+1} = S^n \exp(- (k^{-} + decay_rate) * \Delta t)
 *           - C^n \dfrac{k^{+}}{k^{-}}
 *             \left(\exp(- (k^{-} + decay_rate) * \Delta t) - 1 \right)
 * \f]
 *
 * \param[in]   f_id   field index of scalar which properties are updated
 */
 /*---------------------------------------------------------------------------*/

void cs_gwf_sorbed_concentration_update(const int f_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_real_t *dt = CS_F_(dt)->val;

  /* Get scalar, sorbed concentration fields */
  cs_field_t *sca = cs_field_by_id(f_id);
  int key_sorb = cs_field_key_id("gwf_sorbed_concentration_id");
  cs_field_t *sorb = cs_field_by_id(cs_field_get_key_int(sca, key_sorb));

  /* Get first-order decay rate */
  const int dr_id = cs_field_key_id("fo_decay_rate");
  cs_real_t decay_rate = cs_field_get_key_double(sca, dr_id);

  /* Get soil-water partition structure */
  cs_gwf_soilwater_partition_t sorption_scal;
  const int key_part = cs_field_key_id("gwf_soilwater_partition");
  cs_field_get_key_struct(sca, key_part, &sorption_scal);

  /* Get k+ and k- fields */
  cs_field_t *kp = cs_field_by_id(sorption_scal.ikp);
  cs_field_t *km = cs_field_by_id(sorption_scal.ikm);

  /* Update sorbed concentration */

  /* First choice : analytical resolution */
  if (sorption_scal.anai){
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      /* Case of reversible sorption or decay rate term */
      if (km->val[c_id] + decay_rate > cs_math_epzero) {
        cs_real_t expkdt = exp(-(km->val[c_id] + decay_rate) * dt[c_id]);
        cs_real_t kpskm = kp->val[c_id] / (km->val[c_id] + decay_rate);
        sorb->val[c_id] =  expkdt * sorb->val[c_id]
                         - (expkdt-1.) * kpskm * sca->val[c_id];
      }
      else { /* Irreversible sorption and no decay rate */
        sorb->val[c_id] =  sorb->val[c_id]
                         + dt[c_id]*kp->val[c_id]*sca->val[c_id];
      }
  }
  else { /* Second choice : explicit */
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      sorb->val[c_id] +=  dt[c_id] * (kp->val[c_id] * sca->val[c_id]
                        - (km->val[c_id] + decay_rate) * sorb->val[c_id]);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update liquid concentration according to precipitation phenomenon.
 *
 * If liquid concentration exceeds solubility index, then it is clipped and
 * transferred in precipitated concentration.
 *
 * Note that the decay rate is applied before the clipping.
 *
 * \param[in]   f_id    field index of scalar which properties are updated
 */
 /*---------------------------------------------------------------------------*/

void cs_gwf_precipitation(const int f_id)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_real_t *dt = CS_F_(dt)->val;
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

  /* Get first-order decay rate */
  const int dr_id =  cs_field_key_id("fo_decay_rate");
  cs_real_t decay_rate = cs_field_get_key_double(sca, dr_id);

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    /* Apply the radioactive decay rate to precipitation concentration before
       clipping. To be consistent with decay treatment of liquid concentration,
       an implicit scheme is applied. */
    precip->val[c_id] *= 1. / (1. + decay_rate * dt[c_id]);

    /* Clipping of solute concentration, update of precip. concentration */
    ctot_tmp          = sca->val[c_id] + precip->val[c_id];
    sca->val[c_id]    = CS_MIN(ctot_tmp, solub->val[c_id]);
    precip->val[c_id] = CS_MAX(0., ctot_tmp - solub->val[c_id]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take into account kinetic chemical reaction in evolution equation
 *        of total liquid concentration.
 *
 * \f[
 * S^{n+1} = S^n \exp(- (k^{-} + decay_rate) * \Delta t)
 *           - C^n \dfrac{k^{+}}{k^{-}}
 *             \left(\exp(- (k^{-} + decay_rate) * \Delta t) - 1 \right)
 * \f]
 *
 * \param[in]   f_id   field index of scalar which properties are updated
 */
 /*---------------------------------------------------------------------------*/

void cs_gwf_kinetic_reaction(const int  f_id,
                             cs_real_t *ts_imp,
                             cs_real_t *ts_exp)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  const cs_mesh_quantities_t *mesh_quantities = cs_glob_mesh_quantities;
  const cs_real_t *vol = mesh_quantities->cell_vol;

  cs_real_t *dt = CS_F_(dt)->val;

  /* Get soil density */
  cs_field_t *rosoil = cs_field_by_name("soil_density");

  /* Get scalar and sorbed concentration fields */
  cs_field_t *sca = cs_field_by_id(f_id);
  int key_sorb = cs_field_key_id("gwf_sorbed_concentration_id");
  cs_field_t *sorb = cs_field_by_id(cs_field_get_key_int(sca, key_sorb));

  /* Get first-order decay rate */
  const int dr_id =  cs_field_key_id("fo_decay_rate");
  cs_real_t decay_rate = cs_field_get_key_double(sca, dr_id);

  /* Soil-water partition structure */
  int key_part = cs_field_key_id("gwf_soilwater_partition");
  cs_gwf_soilwater_partition_t sorption_scal;
  cs_field_get_key_struct(sca, key_part, &sorption_scal);

  /* Get k+ and k- fields */
  cs_field_t *kp = cs_field_by_id(sorption_scal.ikp);
  cs_field_t *km = cs_field_by_id(sorption_scal.ikm);

  /* First choice : analytical resolution */
  if (sorption_scal.anai) {
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      /* General case (reversible sorption or presence of decay rate) */
      if (km->val[c_id] + decay_rate > cs_math_epzero) {
        cs_real_t expkdt = exp(-(km->val[c_id]+decay_rate) * dt[c_id]);
        cs_real_t kpskm = kp->val[c_id] / (km->val[c_id] + decay_rate);
        ts_exp[c_id] += - vol[c_id]
                         *(decay_rate * rosoil->val[c_id] * sorb->val[c_id]
                           + rosoil->val[c_id]/dt[c_id] * (1-expkdt)
                            *(kpskm*sca->val[c_id] - sorb->val[c_id]));
        ts_imp[c_id] +=  vol[c_id] * rosoil->val[c_id] / dt[c_id]
                        *(1-expkdt)*kpskm;
      }
      else { /* case of irreversible sorption without decay rate */
        cs_real_t rokpl = rosoil->val[c_id] * kp->val[c_id];
        ts_exp[c_id] += - vol[c_id] * rokpl * sca->val[c_id];
        ts_imp[c_id] += + vol[c_id] * rokpl;
      }
  }
  else { /* Second choice : explicit */
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      ts_exp[c_id] += vol[c_id] * rosoil->val[c_id]
                     *(  km->val[c_id] * sorb->val[c_id]
                       - kp->val[c_id] * sca->val[c_id]);
      ts_imp[c_id] += vol[c_id] * rosoil->val[c_id] * kp->val[c_id];
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
