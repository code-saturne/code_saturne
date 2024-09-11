/*============================================================================
 * Cooling towers related functions
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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_atmo.h"
#include "cs_atmo_profile_std.h"
#include "cs_base.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_zone.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_intprf.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_thermal_model.h"
#include "cs_volume_zone.h"

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute cell reference pressure
 *
 * \param[in]     cell_id       Cell index
 * \param[in]     p0            Fluid properties reference pressure
 * \param[in]     ref_ressure   Atmospheric reference pressure
 * \return        pphy          Reference pressure
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ctwr_compute_reference_pressure(cs_lnum_t  cell_id,
                                   cs_real_t  p0,
                                   cs_field_t *ref_pressure)
{
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_OFF) {
    return p0;
  }

  else {
    cs_real_t pphy = 0;
    cs_real_t dum = 0;
    if (cs_glob_atmo_option->meteo_profile == 0) {
      cs_atmo_profile_std(cell_cen[cell_id][2], &pphy, &dum, &dum);
    }
    else if (cs_glob_atmo_option->meteo_profile == 1) {
      int nbmett = cs_glob_atmo_option->met_1d_nlevels_t;
      int nbmetm = cs_glob_atmo_option->met_1d_ntimes;
      pphy = cs_intprf(nbmett,
          nbmetm,
          cs_glob_atmo_option->z_temp_met,
          cs_glob_atmo_option->time_met,
          cs_glob_atmo_option->hyd_p_met,
          cell_cen[cell_id][2],
          cs_glob_time_step->t_cur);
    }
    else {
      pphy = ref_pressure->val[cell_id];
    }
    return pphy;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reset the field variables based on the restart values
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     humidity0   Reference humidity
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_restart_field_vars(cs_real_t  rho0,
                           cs_real_t  t0,
                           cs_real_t  p0,
                           cs_real_t  humidity0,
                           cs_real_t  molmassrat)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_halo_t *halo = m->halo;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;

  /* Initialize the fields - based on map */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;     /* Humid air (bulk) Cp */
  cs_real_t *t_h = NULL;
  cs_real_t *t_h_a = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    t_h = cs_field_by_name("real_temperature")->val; /* Humid air temp */
    t_h_a = cs_field_by_name("real_temperature")->val_pre; /* Humid air temp */
  }
  else {
    t_h = cs_field_by_name("temperature")->val; /* Humid air temp */
    t_h_a = cs_field_by_name("temperature")->val_pre; /* Humid air temp */
  }

  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* Humid air enthalpy */
  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val;    /* Water mass fraction in
                                                        humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;     /* Saturated humidity */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;     /* Absolute humidity in
                                                        humid air (bulk) */

  /* Packing liquid quantities */
  cs_real_t *t_l_p = (cs_real_t *)CS_F_(t_l_pack)->val;    /* Liquid temperature */
  cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;  /* Liquid enthalpy */
  cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per
                                                         unit cell volume */

  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid vertical
                                                            velocity */

  /* Rain variables */
  cs_field_t *cfld_yp = cs_field_by_name_try("ym_l_r"); /* Rain mass fraction */
  cs_field_t *cfld_taup = cs_field_by_name_try("ym_l_r_drift_tau");
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("ym_l_r_drift_vel");

  cs_real_t *cpro_taup = NULL;
  if (cfld_taup != NULL)
    cpro_taup = cfld_taup->val;
  else
    BFT_MALLOC(cpro_taup, n_cells_with_ghosts, cs_real_t);

  /* Get ct zones information */

  int *_n_ct_zones = cs_get_glob_ctwr_n_zones();
  cs_ctwr_zone_t **_ct_zone = cs_get_glob_ctwr_zone();

  /* Check if there are any leaking packing zones, if yes, there is rain */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  for (int ict = 0; ict < *_n_ct_zones && !(ct_opt->has_rain); ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];
    if (ct->xleak_fac > 0.0)
      ct_opt->has_rain = true;
  }

  const cs_air_fluid_props_t  *air_prop = cs_glob_air_props;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = air_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  /* Recompute the initial values which were used in the initialization of
   * the calculation which is being restarted */
  cs_real_t ym_w_ini = humidity0 / (1.0 + humidity0); // From 'ctiniv'
  if (ym_w_ini < 0.0)
    ym_w_ini = 0;

  if (ym_w_ini >= 1.0)
    ym_w_ini = 1. - cs_math_epzero;

  cs_real_t x_ini = ym_w_ini/(1.0-ym_w_ini);

  cs_real_t t_h_ini = t0 - cs_physical_constants_celsius_to_kelvin;

  cs_real_t rho_h_ini = cs_air_rho_humidair(x_ini,
                                            rho0,
                                            p0,
                                            t0,
                                            molmassrat,
                                            t_h_ini);

  /* Clip counters for humidity / water variables */
  int nclip_yw_min = 0;
  int nclip_yw_max = 0;

  /* Initialize the cooling towers variables */

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field */

    /* Clippings of water mass fraction */
    if (ym_w[cell_id] < 0.0) {
      ym_w[cell_id] = 0;
      nclip_yw_min += 1;
    }

    if (ym_w[cell_id] >= 1.0) {
      ym_w[cell_id] = 1. - cs_math_epzero;
      nclip_yw_max += 1;
    }
    x[cell_id] = ym_w[cell_id]/(1.0 - ym_w[cell_id]);

    /* Bulk humid air temperature at the reference temperature
       This is only calculated once at the beginning so same as
       'cs_ctwr_init_field_vars'
       No, this would be the value at the previous time step -
       At present, it is not stored in the restart file, so for lack
       of information initialize it with the present value of the temperature */
    t_h_a[cell_id] = t_h[cell_id];

    /* Update the humid air enthalpy based on the solved value of T_h */
    //FIXME Need to use the method of 'cs_ctwr_phyvar_update'

    x_s[cell_id] = cs_air_x_sat(t_h[cell_id],p0);

    cp_h[cell_id] = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h[cell_id],
                                     x[cell_id],
                                     x_s[cell_id],
                                     t_h[cell_id]);

    /* Update the liquid temperature based on the solved liquid enthalpy
     * NB: May not be required as it is also done in 'cs_ctwr_phyvar_update'?
     * No, it must be done here because here we sweep over the entire
     * computational domain whereas 'cs_ctwr_phyvar_update' updates
     * T_l only over the packing zones */
    t_l_p[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
    if (y_l_p[cell_id] > 0.) {
      cs_real_t h_liq = yh_l_p[cell_id] / y_l_p[cell_id];
      t_l_p[cell_id] = cs_liq_h_to_t(h_liq);
    }

    /* Initialize the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coefficient:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim =   pow(droplet_diam, 2.) * rho_l / (18. * visc)
                      * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;

    /* Use the same humid air density which was used at the beginning of the
     * calculation being restarted, otherwise since rho_h changes during the
     * calculation, reynolds, v_lim and cpro_taup will end up being different
     * from the initial values used in the calculation being restarted */

    /* Droplet Reynolds number */
    //    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
    cs_real_t reynolds = rho_h_ini * v_lim * droplet_diam / visc;

    for (int sweep = 0;
         sweep < 100 && fabs(reynolds - reynolds_old) > 0.001;
         sweep++) {
      reynolds_old = reynolds;
      v_lim =   pow(droplet_diam, 2.) * rho_l
              / (18. * visc * (1 + 0.15 * pow(reynolds, 0.687)))
              * cs_math_3_norm(gravity);
      //      reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
      reynolds = rho_h_ini * v_lim * droplet_diam / visc;
    }

    cpro_taup[cell_id] = v_lim / cs_math_3_norm(gravity);

    /* Initialize rain variable */
    if (ct_opt->has_rain) {//FIXME useless
      cs_real_3_t *drift_vel = (cs_real_3_t *restrict)(cfld_drift_vel->val);
      drift_vel[cell_id][0] = cpro_taup[cell_id] * gravity[0];
      drift_vel[cell_id][1] = cpro_taup[cell_id] * gravity[1];
      drift_vel[cell_id][2] = cpro_taup[cell_id] * gravity[2];
    }

  }

  cs_gnum_t n_g_clip_yw_min = nclip_yw_min;
  cs_gnum_t n_g_clip_yw_max = nclip_yw_max;

  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_min);
  cs_parall_sum(1, CS_GNUM_TYPE, &n_g_clip_yw_max);

  /* Printing clips in listing */
  if (n_g_clip_yw_min >= 1 || n_g_clip_yw_max >= 1) {
    bft_printf("WARNING : clipping on water mass fraction at restart in"
               "cs_ctwr_restart_field_vars : min_clip = %lu, max_clip = %lu\n",
                n_g_clip_yw_min, n_g_clip_yw_max);
  }

  /* Loop over exchange zones */
  for (int ict = 0; ict < *_n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Initialize the liquid vertical velocity component
       * this is correct for droplet and extended for other packing zones */
      vel_l[cell_id] = ct->v_liq_pack;
    }
  }

  //Check enthalpies
  cs_real_t h_max = -1.e12;
  cs_real_t h_min = 1.e12;
  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    h_min = CS_MIN(h_min,h_h[cell_id]);
    h_max = CS_MAX(h_max,h_h[cell_id]);
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_taup);
    if (cfld_yp != NULL)
      cs_halo_sync_var(halo, CS_HALO_STANDARD, cfld_yp->val);
    if (cfld_drift_vel != NULL) {
      cs_halo_sync_var_strided(halo, CS_HALO_STANDARD, cfld_drift_vel->val, 3);
      if (m->n_init_perio > 0)
        cs_halo_perio_sync_var_vect(halo, CS_HALO_STANDARD,
                                    cfld_drift_vel->val, 3);
    }
  }

  /* Free memory */
  if (cfld_taup != NULL)
    BFT_FREE(cpro_taup);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update the thermo physical properties fields for the humid air and
 *        the liquid.
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_phyvar_update(cs_real_t  rho0,
                      cs_real_t  t0,
                      cs_real_t  p0)
{
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);
  const cs_lnum_t *b_face_cells
    = (const cs_lnum_t *)(cs_glob_mesh->b_face_cells);
  const cs_halo_t *halo = cs_glob_mesh->halo;

  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  /* Fields necessary for humid atmosphere model */
  cs_field_t *meteo_pressure = cs_field_by_name_try("meteo_pressure");
  cs_field_t *yw_liq = cs_field_by_name_try("liquid_water");
  cs_field_t *real_temp = cs_field_by_name_try("real_temperature");
  cs_field_t *beta_h = cs_field_by_name_try("thermal_expansion");

  /* Water / air molar mass ratio */
  const cs_real_t molmassrat = air_prop->molmass_rat;

  cs_real_t *rho_m = (cs_real_t *)CS_F_(rho)->val;    /* Humid air + rain
                                                         (bulk) density */
  cs_real_t *rho_h = cs_field_by_name("rho_humid_air")->val; /* Humid air
                                                                density */
  cs_real_t rho_l = air_prop->rho_l; /* Liquid density */
  cs_real_t *cp_h = (cs_real_t *)CS_F_(cp)->val;      /* Humid air (bulk) Cp */

  /* Fields based on maps */

  cs_real_t *t_h = NULL;
  cs_real_t *theta_liq = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    t_h = cs_field_by_name("real_temperature")->val; /* Humid air temp */
    theta_liq = cs_field_by_name("temperature")->val; /* Liq. pot. temp. */
  }
  else
    t_h = cs_field_by_name("temperature")->val; /* Humid air temp */

  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;       /* Humid air enthalpy */
  cs_real_t *therm_diff_h = cs_field_by_name("thermal_conductivity")->val;
  cs_real_t *cpro_x1 = cs_field_by_name("x_c")->val;
  cs_real_t *bpro_x1 = cs_field_by_name("b_x_c")->val;
  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val;     /* Water mass fraction
                                                         in humid air */
  cs_real_t *vol_f_c = cs_field_by_name("vol_f_c")->val; /* Vol frac.
                                                              cont. phase */
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val;      /* Absolute humidity
                                                         in bulk air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;      /* Saturated humidity */
  cs_real_t *x_rel = cs_field_by_name("x_rel")->val;  /* Relative humidity */

/* Packing zone variables */
  cs_real_t *t_l_p = (cs_real_t *)CS_F_(t_l_pack)->val;     /* Liquid temp */
  cs_real_t *h_l_p = cs_field_by_name("h_l_packing")->val;  /* Liquid enthalpy */
  cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;   /* Ylp.hlp */
  cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per unit
                                                         cell volume*/
  cs_real_t *ym_l_p = cs_field_by_name("ym_l_packing")->val; /* Liquid mass
                                                                * fraction in
                                                                * packing */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid vertical
                                                            velocity */
  cs_real_t *mf_l = cs_field_by_name("mass_flux_l")->val; /* Liquid mass flux */

  cs_real_t *liq_mass_flow
    = cs_field_by_name("inner_mass_flux_y_l_packing")->val; //FIXME

  /* Variable and properties for rain zones */
  cs_field_t *cfld_yp = cs_field_by_name_try("ym_l_r");   /* Rain mass fraction */
  cs_real_t *vol_f_r = cs_field_by_name("vol_f_r")->val; /* Vol frac. rain */

  cs_real_t *ymh_l_r = cs_field_by_name("ymh_l_r")->val; /* Ylr time hlr */
  cs_real_t *h_l_r = cs_field_by_name("h_l_r")->val; /* Rain enthalpy */
  cs_real_t *t_l_r = cs_field_by_name("temp_l_r")->val; /* Rain temperature */

  cs_real_t *ym_l_r = NULL;
  if (cfld_yp != NULL)
    ym_l_r = cfld_yp->val;

  cs_lnum_t n_cells = cs_glob_mesh->n_cells;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  cs_real_t lambda_h = cs_glob_air_props->lambda_h;

  /* Clipping counters for water / humidity variables */
  int nclip_yw_min = 0;
  int nclip_yw_max = 0;
  int nclip_yr_min = 0;
  int nclip_yr_max = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Clippings of water mass fraction */
    if (ym_w[cell_id] < 0.0) {
      ym_w[cell_id] = 0;
      nclip_yw_min += 1;
    }
    if (ym_w[cell_id] >= 1.0) {
      ym_w[cell_id] = 1. - cs_math_epzero;
      nclip_yw_max += 1;
    }

    /* Clippings of rain mass fraction */
    if (ym_l_r != NULL) {
      if (ym_l_r[cell_id] < 0.0) {
        ym_l_r[cell_id] = 0;
        nclip_yr_min += 1;
      }
      if ((ym_l_r[cell_id] + ym_w[cell_id]) >= 1.0) {
        ym_l_r[cell_id] = 1. - ym_w[cell_id] - cs_math_epzero;
        nclip_yr_max += 1;
      }

      /* Recompute rain enthalpy and temperature from Ylr.hlr */
      if (ym_l_r[cell_id] > 1.e-4) {
        h_l_r[cell_id] = ymh_l_r[cell_id]/ym_l_r[cell_id];
        t_l_r[cell_id] = cs_liq_h_to_t(h_l_r[cell_id]);
      }
      else {
        h_l_r[cell_id] = h_h[cell_id];
        t_l_r[cell_id] = t_h[cell_id];
      }
    }

    /* Continuous phase mass fraction */
    cpro_x1[cell_id] = 1. - ym_l_r[cell_id];

    //TODO not one for rain zones - Why not?
    //If it represents the humid air, then it should be one?  If it represents
    //the dry air, then it should account for both ym_l_r and ym_w
    //

    /* Compute cell reference pressure */
    cs_real_t pphy = cs_ctwr_compute_reference_pressure(cell_id, p0, meteo_pressure);
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
      cs_rho_humidair(ym_w[cell_id],
                      theta_liq[cell_id],
                      pphy,
                      &(yw_liq->val[cell_id]),
                      &(real_temp->val[cell_id]),
                      &(rho_h[cell_id]),
                      &(beta_h->val[cell_id]));
    }
    else {
      rho_h[cell_id] = cs_air_rho_humidair(x[cell_id],
                                           rho0,
                                           pphy,
                                           t0,
                                           molmassrat,
                                           t_h[cell_id]);
    }

    /* Homogeneous mixture density */
    rho_m[cell_id] = 1. / ((1. - ym_l_r[cell_id])/rho_h[cell_id]
                           + ym_l_r[cell_id]/rho_l);

    /* Update volume fractions for each phase */
    vol_f_c[cell_id] = cpro_x1[cell_id] * rho_m[cell_id] / rho_h[cell_id];

    vol_f_r[cell_id] = ym_l_r[cell_id] * rho_m[cell_id] / rho_l;

    /* Update humidity field */
    x[cell_id] = ym_w[cell_id]/(1.0 - ym_w[cell_id] - ym_l_r[cell_id]);

    // FIXME for drops - This should be the proportion of 'gaseous' water
    // (dissolved and condensate) in the humid air:
    //   Y(dry air)+ Y(gasesous water) + Y(drops) = 1 in all computational cells
    //   If we do that, then the density needs to be revised as well and the
    //   temperatures of both the bulk (dry air + gaseous water +drops) and the
    //   humid air must be solved for.
    // Here, the approximation is that Y(drops) is negligible
    // This is NOT generally true : on MISTRAL, we reach Y(drops) > 0.5

    /* Saturated humidity */
    x_s[cell_id] = cs_air_x_sat(t_h[cell_id], p0);

    /*Relative humidity */
    x_rel[cell_id] = x[cell_id] / x_s[cell_id];


    /* Update the humid air temperature using new enthalpy but old
     * Specific heat */

    cp_h[cell_id] = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h[cell_id],
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    // Update the humid air enthalpy diffusivity lambda_h if solve for T_h?
    // Need to update since a_0 is variable as a function of T and humidity
    therm_diff_h[cell_id] = lambda_h;
  }

  cs_gnum_t n_g_clip_count[4] = {(cs_gnum_t)nclip_yw_min,
                                 (cs_gnum_t)nclip_yw_max,
                                 (cs_gnum_t)nclip_yr_min,
                                 (cs_gnum_t)nclip_yr_max};

  cs_parall_sum(4, CS_GNUM_TYPE, n_g_clip_count);

  /* Printing clips in listing */
  if (n_g_clip_count[0] >= 1 || n_g_clip_count[1] >= 1) {
    bft_printf("WARNING : clipping on rain mass fraction"
               "in cs_ctwr_phyvar_update : min_clip = %lu, max_clip = %lu\n",
               n_g_clip_count[0], n_g_clip_count[1]);
  }
  if (n_g_clip_count[2] >= 1 || n_g_clip_count[2] >= 1) {
    bft_printf("WARNING : clipping on water mass fraction"
                "in cs_ctwr_phyvar_update : min_clip = %lu, max_clip = %lu\n",
                n_g_clip_count[2], n_g_clip_count[3]);
  }

  /* If solving rain velocity */
  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();
  if (ct_opt->solve_rain_velocity) {
    int class_id = 1;
    char vg_lim_name[80];
    sprintf(vg_lim_name, "vg_lim_p_%02d", class_id);

    /* Drops terminal velocity fields */
    cs_field_t *vg_lim_p = cs_field_by_name(vg_lim_name);
    cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
      cs_glob_physical_constants->gravity[1],
      cs_glob_physical_constants->gravity[2]};
    cs_field_t *cfld_taup = cs_field_by_name_try("ym_l_r_drift_tau");
    cs_real_t *cpro_taup = NULL;
    if (cfld_taup != NULL)
      cpro_taup = cfld_taup->val;

    /* Continuous phase drift velocity */
    cs_field_t *vd_c = cs_field_by_name("vd_c");

    /* Continuous phase velocity */
    cs_field_t *v_c = cs_field_by_name("v_c");

    /* Rain drift velocity variables */
    char f_name[80];
    sprintf(f_name, "vd_p_%02d", class_id);
    cs_field_t *vd_p = cs_field_by_name(f_name);
    sprintf(f_name, "v_p_%02d", class_id);
    cs_field_t *vp = cs_field_by_name(f_name);
    cs_real_3_t *vel = (cs_real_3_t *)CS_F_(vel)->val;

    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (cs_lnum_t i = 0; i < 3; i++) {

      /* Drops terminal velocity */
      vg_lim_p->val[cell_id * 3 + i] = cpro_taup[cell_id] * gravity[i];

      /* Continuous phase drift velocity */
      vd_c->val[cell_id * 3 + i] = 0.;

      /* Rain drops drift velocity calculation */
      if (ym_l_r[cell_id] > 1.e-9) {
        //vp->val[cell_id*3 + i] /= CS_MAX(ym_l_r[cell_id], 1.0e-9);
        vd_p->val[cell_id*3 + i] = vp->val[cell_id*3 + i] - vel[cell_id][i];
      }
    else {
        vd_p->val[cell_id * 3 + i] = 0.;
      }
      vd_c->val[cell_id * 3 + i] -= ym_l_r[cell_id]
        * vd_p->val[cell_id * 3 + i];
      vd_c->val[cell_id * 3 + i] /= cpro_x1[cell_id];
      v_c->val[cell_id * 3 + i] = vel[cell_id][i] + vd_c->val[cell_id * 3 + i];
      }
    }
  }

  /* Cooling tower zones */
  cs_ctwr_zone_t **_ct_zone = cs_get_glob_ctwr_zone();
  const int *_n_ct_zones = cs_get_glob_ctwr_n_zones();

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < *_n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;

    /* Packing zone */
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Update the injected liquid temperature
       * NB: (y_l_p.h_l) is transported and not (h_l) */
      if (y_l_p[cell_id] > 0.) {
        h_l_p[cell_id] = yh_l_p[cell_id] / y_l_p[cell_id];
        t_l_p[cell_id] = cs_liq_h_to_t(h_l_p[cell_id]);
        mf_l[cell_id] = y_l_p[cell_id] * rho_m[cell_id] * vel_l[cell_id];
        ym_l_p[cell_id] = y_l_p[cell_id] / (1 + y_l_p[cell_id]);
      }
    }

    /* Update Inlet packing zone temperature if imposed */
    if (ct->delta_t > 0) {
      /* Recompute outgoing temperature */
      ct->t_l_out = 0.0;

      /* Compute liquid water quantities
       * And humid air quantities at liquid outlet */
      for (cs_lnum_t i = 0; i < ct->n_outlet_faces; i++) {

        cs_lnum_t face_id = ct->outlet_faces_ids[i];
        cs_lnum_t cell_id_l;

        /* Convention: outlet is positive mass flux
         * Then upwind cell for liquid is i_face_cells[][0] */
        int sign = 1;
        if (liq_mass_flow[face_id] < 0) {
          sign = -1;
          cell_id_l = i_face_cells[face_id][1];
        }
        else {
          cell_id_l = i_face_cells[face_id][0];
        }

        /* h_l is in fact (y_l. h_l),
         * and the transport field is (y_l*liq_mass_flow) */
        ct->t_l_out += sign * t_l_p[cell_id_l]
          * y_l_p[cell_id_l] * liq_mass_flow[face_id];
        ct->q_l_out += sign * y_l_p[cell_id_l] * liq_mass_flow[face_id];
      }

      cs_parall_sum(1, CS_REAL_TYPE, &(ct->t_l_out));
      cs_parall_sum(1, CS_REAL_TYPE, &(ct->q_l_out));

      ct->t_l_out /= ct->q_l_out;

      /* Relaxation of ct->t_l_bc */
      ct->t_l_bc = (1. - ct->relax) * ct->t_l_bc
                 + ct->relax * (ct->t_l_out + ct->delta_t);

      /* Clipping between 0 and 100 */
      ct->t_l_bc = CS_MAX(CS_MIN(ct->t_l_bc, 100.), 0.);
    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, x);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, x_s);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_x1);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cp_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, h_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rho_m);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_l_p);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, yh_l_p);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, h_l_p);
  }

  for (cs_lnum_t face_id = 0; face_id < n_b_faces; face_id++) {
    bpro_x1[face_id] = cpro_x1[b_face_cells[face_id]];
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
