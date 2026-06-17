/*============================================================================
 * User initialization prior to solving time steps.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * Initialize variables.
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
cs_user_initialization([[maybe_unused]] cs_domain_t  *domain)
{
  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_3_t *cell_cen
    = (cs_real_3_t *)domain->mesh_quantities->cell_cen;

  /* Example to initialize percentage of type of ground for the ground module */

  /* Initialization of the atmo_ground_percentages field
   * It is of dimension "ground_number +1"
   * default = 0
   *  For 7 ground categories:
   *   water  = 1
   *   forest = 2
   *   divers = 3
   *   rocks  = 4
   *   diffuse building = 5
   *   mix building     = 6
   *   dense building   = 7
   *  For 5 ground categories
   *   water  = 1
   *   forest = 2
   *   divers = 3
   *   mineral = 4
   *   building = 5
   */

  /*! [atmo_ground_init] */
  if (cs_glob_atmo_option->ground_zone_id > -1) {
    const cs_zone_t *z
      = cs_boundary_zone_by_id(cs_glob_atmo_option->ground_zone_id);

    cs_field_t *f = cs_field("atmo_ground_percentages");

    for (cs_lnum_t elt_id = 0; elt_id < z->n_elts; elt_id++) {
      for (cs_lnum_t ground_id = 0; ground_id < f->dim; ground_id++)
        f->val[ground_id + f->dim * elt_id] = 0.;

      /* 100% of mineral */
      f->val[4 + f->dim * elt_id] = 100.;
    }

  }
  /*! [atmo_ground_init] */

  /* Example to initialize variables using 1-D meteo data */

  cs_span<cs_real_t> cvar_k;
  cs_span<cs_real_t> cvar_ep;
  cs_span<cs_real_t> cvar_fb;
  cs_span<cs_real_t> cvar_omg;
  cs_span<cs_real_t> cvar_phi;
  cs_span<cs_real_t> cvar_nusa;
  cs_span_2d<cs_real_t> cvar_rij;

  auto cvar_vel = CS_F_(vel)->get_val_v();

  cs_field_t *th_f = cs_thermal_model_field();

  /*  Initialize variables using an input meteo profile
   *   (only if we are not doing a restart */

  if (cs_glob_turb_model->itytur == 2) {
    cvar_k = CS_F_(k)->get_val_s();
    cvar_ep = CS_F_(eps)->get_val_s();
  }
  else if (cs_glob_turb_model->order == CS_TURB_SECOND_ORDER) {
    cvar_ep = CS_F_(eps)->get_val_s();
    cvar_rij = CS_F_(rij)->get_val_t();
  }
  else if (cs_glob_turb_model->model == CS_TURB_V2F_PHI) {
    cvar_k = CS_F_(k)->get_val_s();
    cvar_ep = CS_F_(eps)->get_val_s();
    cvar_phi = CS_F_(phi)->get_val_s();
    cvar_fb = CS_F_(f_bar)->get_val_s();
  }
  else if (cs_glob_turb_model->model ==  CS_TURB_K_OMEGA) {
    cvar_k = CS_F_(k)->get_val_s();
    cvar_omg = CS_F_(omg)->get_val_s();
  }
  else if (cs_glob_turb_model->model ==  CS_TURB_SPALART_ALLMARAS) {
    cvar_nusa = CS_F_(nusa)->get_val_s();
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

    cvar_vel(c_id, 0) = xuent;
    cvar_vel(c_id, 1) = xvent;
    cvar_vel(c_id, 2) = 0.;

    if (cs_glob_turb_model->itytur == 2) {
      cvar_k[c_id] = xkent;
      cvar_ep[c_id] = xeent;
    }
    else if (cs_glob_turb_model->order == CS_TURB_SECOND_ORDER) {
      cvar_ep[c_id] = xeent;
      for (int ii = 0; ii < 3; ii++)
        cvar_rij(c_id, ii) = d2s3*xkent;
      for (int ii = 3; ii < 6; ii++)
        cvar_rij(c_id, ii) = 0;
    }
    else if (cs_glob_turb_model->model == CS_TURB_V2F_PHI) {
      cvar_fb[c_id] = 0;
      cvar_k[c_id] = xkent;
      cvar_ep[c_id] = xeent;
      cvar_phi[c_id] = d2s3;
    }
    else if (cs_glob_turb_model->model == CS_TURB_K_OMEGA) {
      cvar_k[c_id] = xkent;
      cvar_omg[c_id] = xeent/cs_turb_cmu/xkent;
    }
    else if (cs_glob_turb_model->model == CS_TURB_SPALART_ALLMARAS) {
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
