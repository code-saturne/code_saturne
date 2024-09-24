/*============================================================================
 * User initialization prior to solving time steps.
 *============================================================================*/

/* VERS */

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization-atmospheric.c
 *
 * \brief Initialization example for atmospheric model.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Initialize variables.
 *
 * This function is called at beginning of the computation
 * (restart or not) before the time step loop.
 *
 * This is intended to initialize or modify (when restarted)
 * variable and time step values.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_initialization(cs_domain_t     *domain)
{
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_3_t *cell_cen
    = (cs_real_3_t *)domain->mesh_quantities->cell_cen;

  /* Example to initialize percentage of type of soil for the soil module */

  /* Initialization of the atmo_soil_percentages field
   * It is of dimension "soil_number +1"
   * default = 0
   *  For 7 soil categories:
   *   water  = 1
   *   forest = 2
   *   divers = 3
   *   rocks  = 4
   *   diffuse building = 5
   *   mix building     = 6
   *   dense building   = 7
   *  For 5 soil categories
   *   water  = 1
   *   forest = 2
   *   divers = 3
   *   mineral = 4
   *   building = 5
   */

  /*! [atmo_soil_init] */
  if (cs_glob_atmo_option->soil_zone_id > -1) {
    const cs_zone_t *z
      = cs_boundary_zone_by_id(cs_glob_atmo_option->soil_zone_id);

    cs_field_t *f = cs_field_by_name("atmo_soil_percentages");

    for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {
      for (cs_lnum_t soil_id = 0; soil_id < f->dim; soil_id++)
        f->val[soil_id + f->dim * elt_id] = 0.;

      /* 100% of mineral */
      f->val[4 + f->dim * elt_id] = 100.;
    }

  }
  /*! [atmo_soil_init] */

  /* Example to initialize variables using 1-D meteo data */

  cs_real_t *cvar_k = nullptr;
  cs_real_t *cvar_ep = nullptr;
  cs_real_t *cvar_fb = nullptr;
  cs_real_t *cvar_omg = nullptr;
  cs_real_t *cvar_phi = nullptr;
  cs_real_t *cvar_nusa = nullptr;
  cs_real_6_t *cvar_rij = nullptr;

  cs_real_3_t *cvar_vel = (cs_real_3_t *)CS_F_(vel)->val;

  cs_field_t *th_f = cs_thermal_model_field();

  /*  Initialize variables using an input meteo profile
   *   (only if we are not doing a restart */

  if (cs_glob_turb_model->itytur == 2) {
    cvar_k = CS_F_(k)->val;
    cvar_ep = CS_F_(eps)->val;
  }
  else if (cs_glob_turb_model->itytur == 3) {
    cvar_ep = CS_F_(eps)->val;
    cvar_rij = (cs_real_6_t *)(CS_F_(rij)->val);
  }
  else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
    cvar_k = CS_F_(k)->val;
    cvar_ep = CS_F_(eps)->val;
    cvar_phi = CS_F_(phi)->val;
    cvar_fb = CS_F_(f_bar)->val;
  }
  else if (cs_glob_turb_model->iturb ==  CS_TURB_K_OMEGA) {
    cvar_k = CS_F_(k)->val;
    cvar_omg = CS_F_(omg)->val;
  }
  else if (cs_glob_turb_model->iturb ==  CS_TURB_SPALART_ALLMARAS) {
    cvar_nusa = CS_F_(nusa)->val;
  }

  const cs_real_t d2s3 = 2.0/3.0;

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {

    const cs_real_t xuent = cs_intprf(cs_glob_atmo_option->met_1d_nlevels_d,
                                      cs_glob_atmo_option->met_1d_ntimes,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->u_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xvent = cs_intprf(cs_glob_atmo_option->met_1d_nlevels_d,
                                      cs_glob_atmo_option->met_1d_ntimes,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->v_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xkent = cs_intprf(cs_glob_atmo_option->met_1d_nlevels_d,
                                      cs_glob_atmo_option->met_1d_ntimes,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->ek_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xeent = cs_intprf(cs_glob_atmo_option->met_1d_nlevels_d,
                                      cs_glob_atmo_option->met_1d_ntimes,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->ep_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    cvar_vel[c_id][0] = xuent;
    cvar_vel[c_id][1] = xvent;
    cvar_vel[c_id][2] = 0.;

    if (cs_glob_turb_model->itytur == 2) {
      cvar_k[c_id] = xkent;
      cvar_ep[c_id] = xeent;
    }
    else if (cs_glob_turb_model->itytur == 3) {
      cvar_ep[c_id] = xeent;
      for (int ii = 0; ii < 3; ii++)
        cvar_rij[c_id][ii] = d2s3*xkent;
      for (int ii = 3; ii < 6; ii++)
        cvar_rij[c_id][ii] = 0;
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_V2F_PHI) {
      cvar_fb[c_id] = 0;
      cvar_k[c_id] = xkent;
      cvar_ep[c_id] = xeent;
      cvar_phi[c_id] = d2s3;
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_K_OMEGA) {
      cvar_k[c_id] = xkent;
      cvar_omg[c_id] = xeent/cs_turb_cmu/xkent;
    }
    else if (cs_glob_turb_model->iturb == CS_TURB_SPALART_ALLMARAS) {
      cvar_nusa[c_id] = cs_turb_cmu*cs_math_pow2(xkent)/xeent;
    }

    if (th_f != nullptr) {

      /* Assume the scalar is a potential temperature */
      const cs_real_t tpent = cs_intprf(cs_glob_atmo_option->met_1d_nlevels_max_t,
                                        cs_glob_atmo_option->met_1d_ntimes,
                                        cs_glob_atmo_option->z_temp_met,
                                        cs_glob_atmo_option->time_met,
                                        cs_glob_atmo_option->pot_t_met,
                                        cell_cen[c_id][2],
                                        cs_glob_time_step->t_cur);
      th_f->val[c_id] = tpent;

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_initialization.c
 *
 * \brief Fill in vertical profiles of atmospheric properties prior to solve
 *        1D radiative transfers.
 *
 * \param[in, out] preray        pressure vertical profile
 * \param[in, out] temray        real temperature vertical profile
 * \param[in, out] romray        density vertical profile
 * \param[in, out] qvray         water vapor content vertical profile
 * \param[in, out] qlray         water liquid content vertical profile
 * \param[in, out] ncray         droplets density vertical profile
 * \param[in, out] aeroso        aerosol concentration vertical profile
 */
/*----------------------------------------------------------------------------*/

void
cs_user_atmo_1d_rad_prf(cs_real_t   preray[],
                        cs_real_t   temray[],
                        cs_real_t   romray[],
                        cs_real_t   qvray[],
                        cs_real_t   qlray[],
                        cs_real_t   ncray[],
                        cs_real_t   aeroso[])
{

  /*! [humid_aerosols_atmo] */

  cs_atmo_option_t *at_opt = cs_glob_atmo_option;
  const cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();

  cs_real_t *zvert = at_opt->rad_1d_z;
  cs_real_t *zray  = at_opt->rad_1d_zray;

  const int kvert = at_opt->rad_1d_nlevels;
  const int kmx   = at_opt->rad_1d_nlevels_max;

  const cs_real_t rvsra = phys_pro->rvsra;
  const cs_real_t rair = phys_pro->r_pg_cnst;
  const cs_real_t gz = cs_glob_physical_constants->gravity[2];
  const cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  aeroso[0] = 10.0;

  for (int k = 1; k < kvert; k++) {
    zray[k] = zvert[k];

    const cs_real_t rhum = rair * (1.0 + (rvsra-1.0)*qvray[k]);
    const cs_real_t tmean = 0.50 * (temray[k-1] + temray[k]) + tkelvi;
    const cs_real_t rap = -fabs(gz) * (zray[k]-zray[k-1]) / rhum / tmean;
    preray[k] = preray[k-1] * exp(rap);

    if (zray[k] < 50.0) {
      aeroso[k] = aeroso[0];
    }
    else {
      aeroso[k] = aeroso[0]*exp(-(zray[k]-50.0) / 1.25e3);
      if (aeroso[k] < 5.0)
        aeroso[k] = 5.0;
    }
  }

  /* Filling the additional levels above meshed domain
     (at these levels, pressure, temperature, density profiles have been
     initialized with standard atmosphere profiles) */

  for (int k = kvert; k < kmx; k++) {
    zray[k] = zvert[k];

    /* read meteo data for temperature, water wapor and liquid content in
       upper layers for example to fill temray, qvray, qlray */

    const cs_real_t rhum = rair*(1.0+(rvsra-1.0)*qvray[k]);
    const cs_real_t tmean = 0.50*(temray[k-1]+temray[k]) + tkelvi;
    const cs_real_t rap = -fabs(gz)*(zray[k]-zray[k-1]) / rhum / tmean;
    preray[k] = preray[k-1]*exp(rap);
    romray[k] = preray[k] / (temray[k]+tkelvi) / rhum;

    /* nc not known above the meshed domain
       droplets radius is assumed of mean volume radius = 5 microns */
    ncray[k]
      = 1.e-6*(3.0*romray[k]*qlray[k])/(4.0*cs_math_pi*1.e3*pow(5.e-6, 3.0));

    // similarly, aerosol concentration not known
    aeroso[k] = aeroso[0]*exp(-(zray[k]-50.0) / 1.25e3);
    if (aeroso[k] < 5.0)
      aeroso[k] = 5.0;
  }

  /*! [humid_aerosols_atmo] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
