/*============================================================================
 * User initialization prior to solving time steps.
 *============================================================================*/

/* VERS */

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
  if (cs_glob_atmo_option->soil_model == 1) {
    cs_zone_t *z = cs_boundary_zone_by_id(cs_glob_atmo_option->soil_zone_id);

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

  cs_real_t *cvar_k = NULL;
  cs_real_t *cvar_ep = NULL;
  cs_real_t *cvar_fb = NULL;
  cs_real_t *cvar_omg = NULL;
  cs_real_t *cvar_phi = NULL;
  cs_real_t *cvar_nusa = NULL;
  cs_real_6_t *cvar_rij = NULL;

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

    const cs_real_t xuent = cs_intprf(cs_glob_atmo_option->nbmetd,
                                      cs_glob_atmo_option->nbmetm,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->u_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xvent = cs_intprf(cs_glob_atmo_option->nbmetd,
                                      cs_glob_atmo_option->nbmetm,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->v_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xkent = cs_intprf(cs_glob_atmo_option->nbmetd,
                                      cs_glob_atmo_option->nbmetm,
                                      cs_glob_atmo_option->z_temp_met,
                                      cs_glob_atmo_option->time_met,
                                      cs_glob_atmo_option->ek_met,
                                      cell_cen[c_id][2],
                                      cs_glob_time_step->t_cur);

    const cs_real_t xeent = cs_intprf(cs_glob_atmo_option->nbmetd,
                                      cs_glob_atmo_option->nbmetm,
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

    if (th_f != NULL) {

      /* Assume the scalar is a potential temperature */
      const cs_real_t tpent = cs_intprf(cs_glob_atmo_option->nbmaxt,
                                        cs_glob_atmo_option->nbmetm,
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

END_C_DECLS
