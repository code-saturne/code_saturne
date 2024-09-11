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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_air_props.h"
#include "cs_array.h"
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
#include "cs_restart.h"
#include "cs_thermal_model.h"
#include "cs_velocity_pressure.h"
#include "cs_volume_zone.h"

#include "cs_ctwr.h"
#include "cs_ctwr_physical_properties.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr_initialize.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Initialize cooling towers fields, stage 0
 *----------------------------------------------------------------------------*/

void
cs_ctwr_fields_init0(void)
{
  int has_restart = cs_restart_present();
  cs_halo_t *halo = cs_glob_mesh->halo;

  /* Fluid properties and physical variables */
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_air_fluid_props_t *air_prop = cs_glob_air_props;

  cs_field_t *t_h = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
    t_h = cs_field_by_name("real_temperature"); /* Humid air temp */
  else
    t_h = cs_field_by_name("temperature"); /* Humid air temp */

  cs_field_t *y_l_p = cs_field_by_name("y_l_packing");
  cs_field_t *ym_w  = cs_field_by_name("ym_water");
  cs_field_t *ym_l_r  = cs_field_by_name("ym_l_r");
  cs_field_t *yh_l_p = cs_field_by_name("yh_l_packing");
  cs_field_t *t_l_p = CS_F_(t_l_pack);

  cs_real_t tkelvin = cs_physical_constants_celsius_to_kelvin;
  const cs_real_t xhum = air_prop->humidity0;

  /* Only if the simulation is not a restart from another one */
  if (has_restart == 0) {
    for (cs_lnum_t cell_id = 0; cell_id < cs_glob_mesh->n_cells; cell_id++) {
      /* Humid air */

      if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_OFF)
        t_h->val[cell_id] = fp->t0 - tkelvin;

      ym_w->val[cell_id] = (1 - ym_l_r->val[cell_id]) * xhum / (1. + xhum);

      /* Liquid in packing */
      t_l_p->val[cell_id] = t_h->val[cell_id];
      y_l_p->val[cell_id] = 0.;

    }
    if (halo != NULL) {
      cs_halo_sync_var(halo, CS_HALO_STANDARD, t_h->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, ym_w->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, t_l_p->val);
      cs_halo_sync_var(halo, CS_HALO_STANDARD, y_l_p->val);
    }

    /* Diffusivities of the dry air and the injected liquid
     * TODO : check if overwrites what users have specified */
    const int kvisl0 = cs_field_key_id("diffusivity_ref");

    cs_field_set_key_double(ym_w, kvisl0, 1.e-12);
    cs_field_set_key_double(y_l_p, kvisl0, 1.e-12);

    /* Initialize :
     * - the enthalpies, which are the solution variables
     * - the humidity, which users might have modified if they changed the
     *   mass fraction of the dry air in the humid air */

    cs_ctwr_init_field_vars(fp->ro0, fp->t0, fp->p0, air_prop->molmass_rat);

    if (air_prop->cp_l <= 0 || air_prop->lambda_l <= 0)
      bft_error(__FILE__,__LINE__, 0, _("Negative lambda or cp for liquid"));

    else
      cs_field_set_key_double(yh_l_p, kvisl0,
                              air_prop->lambda_l / air_prop->cp_l);

  }

  else {
    /* TODO (from old ctiniv0 subroutine) Add restarts */
    const int kvisl0 = cs_field_key_id("diffusivity_ref");

    /* Diffusivities of the dry air and the injected liquid */
    cs_field_set_key_double(ym_w, kvisl0, 1.e-12);
    cs_field_set_key_double(y_l_p, kvisl0, 1.e-12);

    /* Restarts - recompute the required properties based on the saved solution
     * variables. For example : the humidity, liquid vertical velocity, etc. */
    cs_ctwr_restart_field_vars(fp->ro0, fp->t0, fp->p0, air_prop->humidity0,
                               air_prop->molmass_rat);
  }
}

/*----------------------------------------------------------------------------
 * Initialize cooling towers fields, stage 1
 *----------------------------------------------------------------------------*/

void
cs_ctwr_fields_init1(void)
{
  cs_halo_t *halo = cs_glob_mesh->halo;

  cs_field_t *t_h = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID)
    t_h = cs_field_by_name("real_temperature"); /* Humid air temp */
  else
    t_h = cs_field_by_name("temperature"); /* Humid air temp */

  cs_field_t *y_l_p = cs_field_by_name("y_l_packing");
  cs_field_t *ym_w  = cs_field_by_name("ym_water");
  cs_field_t *t_l_p = CS_F_(t_l_pack);

  /* Liquid inner mass flux */
  cs_lnum_t iflmas =
    cs_field_get_key_int(y_l_p, cs_field_key_id("inner_mass_flux_id"));
  cs_real_t *i_mass_flux = cs_field_by_id(iflmas)->val;

  /* Liquid boundary mass flux */
  cs_lnum_t iflmab =
    cs_field_get_key_int(y_l_p, cs_field_key_id("boundary_mass_flux_id"));
  cs_real_t *b_mass_flux = cs_field_by_id(iflmab)->val;

  cs_ctwr_init_flow_vars(i_mass_flux);

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_h->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, ym_w->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, t_l_p->val);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, y_l_p->val);
  }

  for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_b_faces; face_id++) {
    b_mass_flux[face_id] = 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the field variables
 *
 * \param[in]     rho0        Reference density of humid air
 * \param[in]     t0          Reference temperature of humid air
 * \param[in]     p0          Reference pressure
 * \param[in]     molmassrat  Dry air to water vapor molecular mass ratio
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_field_vars(cs_real_t  rho0,
                        cs_real_t  t0,
                        cs_real_t  p0,
                        cs_real_t  molmassrat)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_halo_t *halo = m->halo;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_with_ghosts = m->n_cells_with_ghosts;
  const cs_mesh_quantities_t *fvq   = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)fvq->cell_cen;

  /* Fields necessary for humid atmosphere model */
  cs_field_t *meteo_pressure = cs_field_by_name_try("meteo_pressure");
  cs_field_t *yw_liq = cs_field_by_name_try("liquid_water");
  cs_field_t *real_temp = cs_field_by_name_try("real_temperature");
  cs_field_t *beta_h = cs_field_by_name_try("thermal_expansion");

  /* Initialize the fields - based on map */
  cs_real_t *rho_m = (cs_real_t *)CS_F_(rho)->val; /* Air + rain (mixture)
                                                      density */
  cs_real_t *rho_h = cs_field_by_name("rho_humid_air")->val;
  cs_real_t *t_h = NULL;
  cs_real_t *t_h_a = NULL;
  cs_real_t *theta_liq = NULL;
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {
    t_h = cs_field_by_name("real_temperature")->val; /* Humid air temp */
    t_h_a = cs_field_by_name("real_temperature")->val_pre; /* Humid air temp */
    theta_liq = cs_field_by_name("temperature")->val; /* Liq. pot. temp. */
  }
  else {
    t_h = cs_field_by_name("temperature")->val; /* Humid air temp */
    t_h_a = cs_field_by_name("temperature")->val_pre; /* Humid air temp */
  }
  cs_real_t *h_h = (cs_real_t *)CS_F_(h)->val;     /* Humid air enthalpy */
  cs_real_t *ym_w = (cs_real_t *)CS_F_(ym_w)->val;  /* Water mass fraction in
                                                      humid air */
  cs_real_t *x_s = cs_field_by_name("x_s")->val;
  cs_real_t *x = (cs_real_t *)CS_F_(humid)->val; /* Humidity in air (bulk) */

  cs_real_t *t_l_p = (cs_real_t *)CS_F_(t_l_pack)->val;  /* Liquid temperature */
  cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;  /* Liquid enthalpy */
  cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass per unit */

  /* Packing zone liquid vertical velocity component */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val;

  /* Rain drops variables */
  cs_field_t *cfld_yp = cs_field_by_name_try("ym_l_r"); /* Rain mass fraction */
  cs_real_t *ym_l_r = cfld_yp->val; /* Rain mass fraction */
  cs_field_t *cfld_taup = cs_field_by_name_try("ym_l_r_drift_tau");
  cs_field_t *cfld_drift_vel = cs_field_by_name_try("ym_l_r_drift_vel");

  const cs_real_3_t *vel = (const cs_real_3_t *)(CS_F_(vel)->val);
  cs_field_t *cfld_vc = cs_field_by_name_try("v_c");
  cs_real_3_t *v_c;
  if (cfld_vc != NULL)
    v_c = (cs_real_3_t *)cfld_vc->val;

  /* Phases volume fractions */
  cs_real_t *vol_f_c = cs_field_by_name("vol_f_c")->val;
  cs_real_t *vol_f_r = cs_field_by_name("vol_f_r")->val;

  cs_ctwr_option_t *ct_opt = cs_get_glob_ctwr_option();

  cs_real_t *cpro_taup = NULL;
  if (cfld_taup != NULL)
    cpro_taup = cfld_taup->val;
  else
    BFT_MALLOC(cpro_taup, n_cells_with_ghosts, cs_real_t);

  const cs_air_fluid_props_t  *air_prop = cs_glob_air_props;
  cs_real_t rho_l = air_prop->rho_l;
  cs_real_t visc = cs_glob_fluid_properties->viscl0;
  cs_real_t droplet_diam = air_prop->droplet_diam;

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  /* Count clippings for rain / humidity variables */
  cs_gnum_t nclip_ym_w_min = 0;
  cs_gnum_t nclip_ym_w_max = 0;

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {

    /* Update humidity field in case users have updated the initial
       dry air mass fraction.
       Note: this is a bit dubious as users could also have chosen
       to reset the humidity ? */

    /* Clippings of water mass fraction */
    if (ym_w[cell_id] < 0.0) {
      ym_w[cell_id] = 0;
      nclip_ym_w_min += 1; //TODO : print it
    }

    if (ym_w[cell_id] >= 1.0) {
      ym_w[cell_id] = 1. - cs_math_epzero;
      nclip_ym_w_max += 1; //TODO : print it
    }
    x[cell_id] = ym_w[cell_id]/(1.0 - ym_w[cell_id] - ym_l_r[cell_id]);

    /* Update the humid air density */
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] == CS_ATMO_HUMID) {

      /* Compute local reference pressure for cell_id*/
      cs_real_t pphy = 0;
      cs_real_t dum = 0;
      if (cs_glob_atmo_option->meteo_profile == 0)
        cs_atmo_profile_std(cell_cen[cell_id][2], &pphy, &dum, &dum);

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
        pphy = meteo_pressure->val[cell_id];
      }

      cs_rho_humidair(ym_w[cell_id],
                      theta_liq[cell_id],
                      pphy,
                      &(yw_liq->val[cell_id]),
                      &(real_temp->val[cell_id]),
                      &(rho_h[cell_id]),
                      &(beta_h->val[cell_id]));
    }
    else {

      /* Bulk humid air temperature */
      t_h[cell_id] = t0 - cs_physical_constants_celsius_to_kelvin;
      t_h_a[cell_id] = t_h[cell_id];

      rho_h[cell_id] = cs_air_rho_humidair(x[cell_id],
                       rho0,
                       p0,
                       t0,
                       molmassrat,
                       t_h[cell_id]);
    }

    rho_m[cell_id] = 1. / ((1. - ym_l_r[cell_id])/rho_h[cell_id]
                           + ym_l_r[cell_id]/rho_l);

    /* Update volume fractions */
    vol_f_c[cell_id] = (1 - ym_l_r[cell_id]) * rho_m[cell_id] / rho_h[cell_id];
    vol_f_r[cell_id] = ym_l_r[cell_id] * rho_m[cell_id] / rho_l;

    /* Update the humid air enthalpy */
    x_s[cell_id] = cs_air_x_sat(t_h[cell_id],p0);

    cs_real_t cp_h = cs_air_cp_humidair(x[cell_id], x_s[cell_id]);

    h_h[cell_id] = cs_air_h_humidair(cp_h,
                                      x[cell_id],
                                      x_s[cell_id],
                                      t_h[cell_id]);

    /* Initialize the liquid vertical velocity component
     * this is correct for droplet and extended for other packing zones
     * NB: this value is derived from the drag coefficient:
     * C_D = 24 / Re * (1 + 0.15 * Re^0.687)
     * See ZOPLU HT-31-08-06 */

    cs_real_t v_lim =   cs_math_pow2(droplet_diam) * rho_l / (18. * visc)
                      * cs_math_3_norm(gravity);

    cs_real_t reynolds_old = 0.;
    cs_real_t reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;

    // FIXME make it global for the zone as restart...
    for (int sweep = 0;
         sweep < 100 && CS_ABS(reynolds - reynolds_old) > 0.001;
         sweep++) {
      reynolds_old = reynolds;
      v_lim =   pow(droplet_diam, 2.) * rho_l
              / (18. * visc * (1. + 0.15 * pow(reynolds, 0.687)))
              * cs_math_3_norm(gravity);
      reynolds = rho_h[cell_id] * v_lim * droplet_diam / visc;
    }

    cpro_taup[cell_id] = v_lim / cs_math_3_norm(gravity);

    /* Initialize rain variables (note that Yp is already set to 0) */
    if (ct_opt->has_rain) {
      cs_real_3_t *drift_vel = (cs_real_3_t *restrict)(cfld_drift_vel->val);
      drift_vel[cell_id][0] = cpro_taup[cell_id] * gravity[0];
      drift_vel[cell_id][1] = cpro_taup[cell_id] * gravity[1];
      drift_vel[cell_id][2] = cpro_taup[cell_id] * gravity[2];
    }

    /* Initialize rain velocity variables */
    if (cfld_vc != NULL) {
      v_c[cell_id][0] = vel[cell_id][0];
      v_c[cell_id][1] = vel[cell_id][1];
      v_c[cell_id][2] = vel[cell_id][2];
    }
  }
  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, rho_h);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, cpro_taup);
  }

  /* Cooling tower zones */
  cs_ctwr_zone_t **_ct_zone = cs_get_glob_ctwr_zone();
  int *_n_ct_zones = cs_get_glob_ctwr_n_zones();

  /* Loop over exchange zones */
  for (int ict = 0; ict < *_n_ct_zones; ict++) {

    cs_ctwr_zone_t *ct = _ct_zone[ict];

    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (cs_lnum_t i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];

      /* Initialize with the injection water temperature */
      t_l_p[cell_id] = ct->t_l_bc;

      /* Update the injected liquid enthalpy */
      yh_l_p[cell_id] = cs_liq_t_to_h(t_l_p[cell_id]);

      /* Initialize the liquid vertical velocity component
       * this is correct for droplet and extended for other packing zones */
      vel_l[cell_id] = ct->v_liq_pack;

      /* Note that rho_h * Y_l * vel_l * Stot = q_l_bc */
      cs_real_t y_l_bc =   ct->q_l_bc
                         / (rho_h[cell_id] * vel_l[cell_id] * ct->surface);

      /* Initialize the liquid transported variables:
         liquid mass and enthalpy corrected by the density ratio */
      y_l_p[cell_id] = y_l_bc;

      /* The transported value is (y_l.h_l) and not (h_l) */
      yh_l_p[cell_id] *= y_l_p[cell_id];
    }
  }

  /* Parallel synchronization */
  if (halo != NULL) {
    cs_halo_sync_var(halo, CS_HALO_STANDARD, vel_l);
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
  if (cfld_taup == NULL)
    BFT_FREE(cpro_taup);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the flow variables relevant to the cooling tower scalars
 * inside the packing zones
 *
 * \param[in,out] liq_mass_flow Liquid mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_ctwr_init_flow_vars(cs_real_t  liq_mass_flow[])
{
  cs_real_t *y_l_p = (cs_real_t *)CS_F_(y_l_pack)->val; /* Liquid mass fraction
                                                         in packing */
  cs_real_t *yh_l_p = (cs_real_t *)CS_F_(yh_l_pack)->val;      /* Liquid enthalpy
                                                         in packing */
  cs_real_t *t_l_p = (cs_real_t *)CS_F_(t_l_pack)->val;      /* Liquid temperature
                                                         in packing */

  cs_real_t *rho_h = (cs_real_t *)CS_F_(rho)->val; /* Humid air
                                                      (bulk) density */
  cs_real_t *vel_l = cs_field_by_name("vertvel_l")->val; /* Liquid velocity
                                                            in packing */

  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->i_face_normal;
  const cs_lnum_2_t *i_face_cells =
    (const cs_lnum_2_t *)(cs_glob_mesh->i_face_cells);

  const cs_lnum_t n_cells_with_ghosts = cs_glob_mesh->n_cells_with_ghosts;
  const cs_lnum_t n_i_faces = cs_glob_mesh->n_i_faces;

  const cs_halo_t *halo = cs_glob_mesh->halo;
  cs_lnum_t *packing_cell;

  /* Normalized gravity vector */

  cs_real_t gravity[] = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  cs_real_t g_dir[3];
  cs_math_3_normalize(gravity, g_dir);

  /* Initialize the liquid mass flux to null */
  cs_array_real_fill_zero(n_i_faces, liq_mass_flow);

  /* Tag and initialize the ct values in the packing zone cells */

  BFT_MALLOC(packing_cell, n_cells_with_ghosts, int);

  cs_array_int_set_value(n_cells_with_ghosts, -1, packing_cell);

  /* Cooling tower zones */
  cs_ctwr_zone_t **_ct_zone = cs_get_glob_ctwr_zone();
  const int *_n_ct_zones = cs_get_glob_ctwr_n_zones();

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < *_n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_MALLOC(ct->inlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_faces_ids, n_i_faces, cs_lnum_t);
    BFT_MALLOC(ct->outlet_cells_ids, n_i_faces, cs_lnum_t);
    const cs_lnum_t *ze_cell_ids = cs_volume_zone_by_name(ct->name)->elt_ids;
    for (int i = 0; i < ct->n_cells; i++) {
      cs_lnum_t cell_id = ze_cell_ids[i];
      if (ct->type != CS_CTWR_INJECTION)
        packing_cell[cell_id] = ict;
    }
  }

  /* Parallel synchronization */
  if (halo != NULL)
    cs_halo_sync_untyped(halo, CS_HALO_STANDARD, sizeof(int), packing_cell);

  /* Initialize the liquid mass flux at packing zone faces
   * and the ghost cells for the liquid mass and enthalpy
   * Initialize the couples (inlet faces, upwind cells) and
   * (outlet faces, upwind cells) arrays */

  for (cs_lnum_t face_id = 0; face_id < n_i_faces; face_id++) {

    cs_lnum_t cell_id_1 = i_face_cells[face_id][0];
    cs_lnum_t cell_id_2 = i_face_cells[face_id][1];

    /* one of neigh. cells is in packing */
    if (packing_cell[cell_id_1] != -1 || packing_cell[cell_id_2] != -1) {

      int ct_id = CS_MAX(packing_cell[cell_id_1], packing_cell[cell_id_2]);
      cs_ctwr_zone_t *ct = _ct_zone[ct_id];

      /* Vertical (align with gravity) component of the surface vector */
      cs_real_t liq_surf = cs_math_3_dot_product(g_dir,
                                                 i_face_normal[face_id]);

      /* Face mass flux of the liquid */
      cs_lnum_t cell_id;
      if (liq_surf > 0.) { /* cell_id_1 is upwind cell for liq. flow */
        if (packing_cell[cell_id_1] != -1) /* cell_id_1 in the packing */
          cell_id = cell_id_1;
        else /* cell_id_1 in HALO of the packing and outside of it */
          cell_id = cell_id_2;
      }
      else { /* cell_id_2 is upwind cell for liq. flow */
        if (packing_cell[cell_id_2] != -1) /* cell_id_2 in the packing */
          cell_id = cell_id_2;
        else /* cell_id_2 in HALO of the packing and outside of it */
          cell_id = cell_id_1;
      }

      cs_real_t y_l_bc = ct->q_l_bc / (  rho_h[cell_id] * vel_l[cell_id]
                                       * ct->surface);
      liq_mass_flow[face_id] = rho_h[cell_id] * vel_l[cell_id] * liq_surf;

      /* Initialize a band of ghost cells on the top side of the
         packing zone in order to impose boundary values
         Take the upwind value for initialization */

      /* cell_id_1 in packing and not cell_id_2 */
      if (packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] == -1) {

        /* cell_id_2 is an inlet halo */
        if (liq_mass_flow[face_id] < 0.0) {

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
          y_l_p[cell_id_2] = y_l_bc;
          t_l_p[cell_id_2] = ct->t_l_bc;
          yh_l_p[cell_id_2] = cs_liq_t_to_h(ct->t_l_bc);
          /* The transported value is (y_l.h_l) and not (h_l) */
          yh_l_p[cell_id_2] *= y_l_p[cell_id_2];
        }
        /* face_id is an outlet */
        else {

          /* cell_id_2 is an outlet halo */
          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_2;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;

          ct->surface_out += liq_surf;
        }
      }
      else if (packing_cell[cell_id_1] == -1 && packing_cell[cell_id_2] >= 0) {

        /* cell_id_1 is an inlet halo */
        if (liq_mass_flow[face_id] > 0.0) {

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
          y_l_p[cell_id_1] = y_l_bc;
          t_l_p[cell_id_1] = ct->t_l_bc;
          yh_l_p[cell_id_1] = cs_liq_t_to_h(ct->t_l_bc);
          /* The transported value is (y_l.h_l) and not (h_l) */
          yh_l_p[cell_id_1] *= y_l_p[cell_id_1];
        }
        /* cell_id_1 is an outlet */
        else {
          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_1;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;

          ct->surface_out += liq_surf;
        }

        /* Neighbouring zones, inlet for one, outlet fot the other */
      }
      else if (   packing_cell[cell_id_1] >= 0 && packing_cell[cell_id_2] >= 0
               && packing_cell[cell_id_1] != packing_cell[cell_id_2]) {

        /* cell_id_1 is an inlet for CT2, an outlet for CT1 */
        if (liq_mass_flow[face_id] > 0.0) {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_1;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;
          ct->surface_out += liq_surf;

        }
        /* cell_id_2 is an inlet for CT1, an outlet for CT2 */
        else {
          /* CT2 */
          ct = _ct_zone[packing_cell[cell_id_2]];

          ct->outlet_faces_ids[ct->n_outlet_faces] = face_id;
          ct->outlet_cells_ids[ct->n_outlet_cells] = cell_id_2;

          ct->n_outlet_faces ++;
          ct->n_outlet_cells ++;
          ct->surface_out += liq_surf;

          /* CT1 */
          ct = _ct_zone[packing_cell[cell_id_1]];

          ct->inlet_faces_ids[ct->n_inlet_faces] = face_id;

          ct->n_inlet_faces ++;
          ct->surface_in += liq_surf;
        }

      }
    }
    else {
      liq_mass_flow[face_id] = 0.0;
    }
  }

  /* Loop over Cooling tower zones */
  for (int ict = 0; ict < *_n_ct_zones; ict++) {
    cs_ctwr_zone_t *ct = _ct_zone[ict];

    BFT_REALLOC(ct->inlet_faces_ids, ct->n_inlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_faces_ids, ct->n_outlet_faces, cs_lnum_t);
    BFT_REALLOC(ct->outlet_cells_ids, ct->n_outlet_cells, cs_lnum_t);

    cs_parall_sum(1, CS_REAL_TYPE, &(ct->surface_in));
    cs_parall_sum(1, CS_REAL_TYPE, &(ct->surface_out));
  }

  BFT_FREE(packing_cell);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
