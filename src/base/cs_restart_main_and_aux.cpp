/*============================================================================
 * write/read main and auxiliary checkpoint files
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


/* These two macros only work when "r" and "dummy_XXX" exist.
 * Both have the merit of simplifying the calls -> Only the section name
 * is needed
 */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_ale.h"
#include "cs_array.h"
#include "cs_atmo.h"
#include "cs_boundary_conditions.h"
#include "cs_coal.h"
#include "cs_combustion_gas.h"
#include "cs_elec_model.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_mobile_structures.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbulence_model.h"
#include "cs_turbomachinery.h"
#include "cs_velocity_pressure.h"
#include "cs_vof.h"
#include "cs_wall_condensation.h"
#include "cs_wall_condensation_1d_thermal.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_restart_main_and_aux.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_restart_checkpoint_io.c
        Read and write functions for main and auxiliary checkpoint files.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _WRITE_INT_VAL(_sec) \
  cs_restart_write_section(r, _sec, 0, 1, \
                           CS_TYPE_int, &dummy_int);

#define _WRITE_REAL_VAL(_sec) \
  cs_restart_write_section(r, _sec, 0, 1, \
                           CS_TYPE_cs_real_t, &dummy_real);

#define _READ_INT_VAL(_sec) \
  cs_restart_read_section(r, _sec, 0, 1, \
                          CS_TYPE_int, &dummy_int);

#define _READ_REAL_VAL(_sec) \
  cs_restart_read_section(r, _sec, 0, 1, \
                          CS_TYPE_cs_real_t, &dummy_real);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write main checkpoint file
 *----------------------------------------------------------------------------*/

static void
_write_main_checkpoint(void)
{
  cs_log_printf(CS_LOG_DEFAULT, "** Writing the main restart file\n");
  cs_log_printf(CS_LOG_DEFAULT, "   -----------------------------\n");

  /* Open restart file */
  cs_restart_t *r = cs_restart_create("main.csc", nullptr, CS_RESTART_MODE_WRITE);

  cs_log_printf(CS_LOG_DEFAULT, " Start writing\n");

  /* Write restart version (for version x.y.z, xxyyzz) */
  int dummy_int = 400000;
  _WRITE_INT_VAL("code_saturne:checkpoint:main:version");

  /* Main fields metadata */
  cs_restart_write_field_info(r);

  /* Number of time steps, previous time */
  dummy_int = cs_glob_time_step->nt_cur;
  _WRITE_INT_VAL("nbre_pas_de_temps");

  cs_real_t dummy_real = cs_glob_time_step->t_cur;
  _WRITE_REAL_VAL("instant_precedent");

  /* Turbulence model */
  dummy_int = cs_glob_turb_model->iturb;
  _WRITE_INT_VAL("turbulence_model");

  /* ALE */
  dummy_int = (int)cs_glob_ale;
  _WRITE_INT_VAL("methode_ALE");

  /* VoF */
  dummy_int = (int)cs_glob_vof_parameters->vof_model;
  _WRITE_INT_VAL("vof");

  /* Turbumachinery */
  cs_turbomachinery_restart_write(r);

  /* Atmo chemistry */
  if (cs_glob_atmo_chemistry->model > 0 ||
      cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF) {
    dummy_int = cs_atmo_chemistry_need_initialization();
    _WRITE_INT_VAL("atmospheric_chem");
  }

  cs_log_printf(CS_LOG_DEFAULT, " End writing the options\n");

  /* Main variables */

  cs_restart_write_variables(r, 0, nullptr);
  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
      if (eqp->ibdtso > 1) {
        for (int i = 1; i < eqp->ibdtso; i++)
          cs_restart_write_field_vals(r, f_id, i);
      }
    }
  }

  cs_restart_write_fields(r, CS_RESTART_MAIN);

  /* Notebook variables */
  cs_restart_write_notebook_variables(r);

  /* Close file */
  cs_restart_destroy(&r);

  cs_log_printf(CS_LOG_DEFAULT, " End writing\n");
}

/*----------------------------------------------------------------------------
 * Read main checkpoint file
 *----------------------------------------------------------------------------*/

static cs_map_name_to_id_t *
_read_main_checkpoint(void)
{
  cs_log_printf(CS_LOG_DEFAULT, "   READING THE MAIN RESTART FILE\n");

  /* Open file */
  cs_restart_t *r = cs_restart_create("main.csc", nullptr, CS_RESTART_MODE_READ);

  cs_log_printf(CS_LOG_DEFAULT, " Start reading\n");

  cs_real_t dummy_real = -99999.;

  /* Check restart file */

  int retval
    = cs_restart_read_section_compat(r,
                                     "code_saturne:checkpoint:main:version",
                                     "version_fichier_suite_principal",
                                     0, 1,
                                     CS_TYPE_int,
                                     &dummy_real);

  if (retval != CS_RESTART_SUCCESS) {
    retval = cs_restart_check_if_restart_from_ncfd(r);
    if (retval == 0) {
      cs_log_warning(_("STOP WHILE READING THE MAIN RESTART FILE\n\n"
                       "  WRONG FILE TYPE\n\n"
                       "The file \"main.csc\" does not look like a proper\n"
                       "  main restart file.\n\n"
                       "The calculation cannot be executed.\n\n"
                       "Please make sure the file used as a restart file\n"
                       "  actually is a correct main restart file.\n"));
      cs_exit(EXIT_FAILURE);
    }
  }

  /* Check base locations */
  bool match_cell, match_i_face, match_b_face, match_vertex;
  cs_restart_check_base_location(r,
                                 &match_cell,
                                 &match_i_face,
                                 &match_b_face,
                                 &match_vertex);

  /* For the moment only cell location is tested */
  if (match_cell == false) {
    cs_log_warning(_("STOP READING THE MAIN RESTART FILE\n\n"
                     "  INCONSISTANT RESTART AND CHECKPOINT DATA\n\n"
                     "The number of cells has changed\n\n"
                     "The calculation cannot be executed.\n\n"
                     "Please make sure the file used as restart file does\n"
                     "    correspond to your case\n"));
    cs_exit(EXIT_FAILURE);
  }

  /* Read field info */
  cs_map_name_to_id_t *old_field_map = nullptr;
  cs_restart_read_field_info(r, &old_field_map);

  cs_log_printf(CS_LOG_DEFAULT, " Reading dimensions complete\n");

  /* ------------------------------- */
  /* Read options and stop if needed */
  /* ------------------------------- */
  int nb_errors = 0;

  int dummy_int = 0;
  retval = _READ_INT_VAL("nbre_pas_de_temps");

  // If section doesnt exist, check if it is a restart from neptune
  if (retval != CS_RESTART_SUCCESS) {
    retval = _READ_INT_VAL("ntcabs");
  }

  nb_errors += retval;

  dummy_real = 0.;
  retval = _READ_REAL_VAL("instant_precedent");

  // If section doesnt exist, check if it is a restart from neptune
  if (retval != CS_RESTART_SUCCESS) {
    retval = _READ_REAL_VAL("ttcabs");
  }

  nb_errors += retval;

  /* Set both values still stored in the two dummy variables */
  cs_time_step_define_prev(dummy_int, dummy_real);

  /* Exit if error */
  if (nb_errors != 0) {
    cs_log_warning(_("STOP WHILE READING THE MAIN RESTART FILE\n\n"
                     "  ERROR READING THE TEMPORAL INFORMATION\n\n"
                     "The computation cannot be executed.\n\n"
                     "Please check the integrity of the file used as\n"
                     "    restart file\n"));
    cs_exit(EXIT_FAILURE);
  }

  /* Stop if requested time/iterations number is too small */
  if (cs_glob_time_step->t_max >= 0.) {
    if (cs_glob_time_step->t_prev > cs_glob_time_step->t_max) {
      cs_log_warning(_("STOP READING THE MAIN RESTART FILE\n\n"
                       "  PREVIOUS TIME t_prev = %12.4e\n"
                       "  TIME WANTED   t_max  = %12.4e\n\n"
                       "The requested time, t_max, has to be greater than\n"
                       "  the already simulated time t_prev.\n\n"
                       "The calculation cannot be executed.\n\n"
                       "Please check (increase) t_max.\n"
                       "Please make sure the file used as restart file does\n"
                       "      correspond to your case.\n"),
                      cs_glob_time_step->t_prev,
                      cs_glob_time_step->t_max);
      cs_exit(EXIT_FAILURE);
    }
  }
  else if (cs_glob_time_step->nt_prev > cs_glob_time_step->nt_max) {
    cs_log_warning(_("STOP READING THE MAIN RESTART FILE\n\n"
                     "  NUMBER OF THE PREVIOUS TIME STEP nt_prev = %d\n"
                     "  NUMBER OF TIME STEPS WANTED      nt_max  = %d\n\n"
                     "The requested number of time steps (absolute),\n"
                     "  nt_max, has to be greater than\n"
                     "  the number of time steps already run, t_prev.\n\n"
                     "The calculation cannot be executed.\n\n"
                     "Please check (increase) nt_max.\n"
                     "Please make sure the file used as restart file does\n"
                     "      correspond to your case.\n"),
                    cs_glob_time_step->nt_prev,
                    cs_glob_time_step->nt_max);
    cs_exit(EXIT_FAILURE);
  }

  cs_log_printf(CS_LOG_DEFAULT,
                _("  Reading the previous time step number \n"
                  "  (restarting computation)  nt_prev = %d\n"),
                cs_glob_time_step->nt_prev);

  cs_log_printf(CS_LOG_DEFAULT,
                _("  Reading the previous time step value \n"
                  "  (restarting computation)   t_prev = %12.4e\n"),
                cs_glob_time_step->t_prev);

  /* ALE */
  dummy_int = 0;
  retval = _READ_INT_VAL("methode_ALE");

  cs_ale_type_t ale_ =  CS_ALE_NONE;

  if (retval != CS_RESTART_SUCCESS) {
    if (cs_glob_ale != CS_ALE_NONE) {
      cs_log_warning(_("STOP READING THE MAIN RESTART FILE\n\n"
                       "  ERROR AT READING THE INDICATOR OF ALE METHOD\n\n"
                       "The read restart file might come from a previous\n"
                       "  version of Code Saturne, without ALE.\n"
                       "The calculation will be executed but\n"
                       "  ALE data will be reset.\n"
                       "Please check the integrity of the file used as\n"
                       "    restart file, however.\n"));
    }
  }
  else
    ale_ = (cs_ale_type_t)dummy_int;

  /* Auxiliary file needs to be read if previous computation was already
   * using ALE.
   */
  if (cs_glob_ale != CS_ALE_NONE &&
      ale_ != CS_ALE_NONE &&
      cs_glob_restart_auxiliary->read_auxiliary != 1) {
      cs_log_warning(_("STOP READING THE MAIN RESTART FILE\n\n"
                     "  IALE INDICATOR OF THE PREVIOUS CALCULATION = %d\n"
                     "  IALE INDICATOR OF THE CURRECT CALCULATION  = %d\n\n"
                     "The coordinates of the mesh nodes need to be read again.\n"
                     "  They are stored in the auxiliary restart file.\n"
                     "Therefore the ILEAUX indicator needs to be equal to 1.\n"
                     "Its current value is ILEAUX = %d\n\n"
                     "The calculation cannot be executed.\n\n"
                     "Please check the value of ILEAUX.\n"),
                   ale_,
                   cs_glob_ale,
                   cs_glob_restart_auxiliary->read_auxiliary);
    cs_exit(EXIT_FAILURE);
  }

  /* VoF */
  dummy_int = 0;
  retval = _READ_INT_VAL("vof");

  if (retval != CS_RESTART_SUCCESS) {
    if (cs_glob_vof_parameters->vof_model & CS_VOF_ENABLED)
      cs_log_warning(_("ERROR WHILE READING THE MAIN RESTART FILE\n\n"
                       "  ERROR READING THE INDICATOR OF THE\n"
                       "  VOLUME OF FLUID METHOD\n\n"
                       "The read restart file might come from a previous\n"
                       "  version of code_saturne, without VOF.\n"
                       "The calculation will be executed but\n"
                       "  Volume of Fluid method data will be reset.\n"
                       "Please check the integrity of the file used as\n"
                       "    restart file, however.\n"));
  }

  /* Instant de maillage mobile precedent (rotor/stator) */
  if (cs_turbomachinery_get_model() != CS_TURBOMACHINERY_NONE)
    cs_turbomachinery_restart_read(r);

  cs_log_printf(CS_LOG_DEFAULT, " Reading options complete\n");

  /* -------------- */
  /* Read variables */
  /* -------------- */

  cs_restart_read_variables(r, old_field_map, 0, nullptr);

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp->ibdtso > 1) {
        int ierr = 0;
        for (int i = 1; i < eqp->ibdtso; i++) {
          retval = cs_restart_read_field_vals(r, f->id, i);
          if (retval != CS_RESTART_SUCCESS)
            ierr += 1;
        }
        if (ierr != 0)
          eqp->ibdtso = -eqp->ibdtso;
      }
    }
  }

  cs_restart_read_fields(r, CS_RESTART_MAIN);

  /* ------------------------------- */
  /* Read atmospheric chemistry data */
  /* ------------------------------- */
  if (cs_glob_atmo_chemistry->model > 0 ||
      cs_glob_atmo_chemistry->aerosol_model != CS_ATMO_AEROSOL_OFF) {
    dummy_int = 1;
    retval = _READ_INT_VAL("atmospheric_chem");

    if (retval == CS_RESTART_SUCCESS && dummy_int > 0)
      cs_atmo_chemistry_initialization_deactivate();
  }

  /* ---------- */
  /* Close file */
  /* ---------- */

  cs_restart_destroy(&r);

  return old_field_map;
}

/*----------------------------------------------------------------------------
 * Write auxiliary checkpoint file
 *----------------------------------------------------------------------------*/

static void
_write_auxiliary_checkpoint(void)
{
  cs_log_printf(CS_LOG_DEFAULT, "** Writing the auxiliary restart file\n");
  cs_log_printf(CS_LOG_DEFAULT, "   ----------------------------------\n");

  cs_restart_t *r = cs_restart_create("auxiliary.csc",
                                      nullptr,
                                      CS_RESTART_MODE_WRITE);

  cs_log_printf(CS_LOG_DEFAULT, " Start writing\n");

  /* Restart version (for version x.y.z, xxyyzz) */
  int dummy_int = 400000;
  _WRITE_INT_VAL("code_saturne:checkpoint:auxiliary:version");

  /* ---------- */
  /* Dimensions */
  /* ---------- */

  /* Variable time step indicator */
  dummy_int = (int)cs_glob_time_step_options->idtvar;
  _WRITE_INT_VAL("indic_dt_variable");

  dummy_int = (int)cs_glob_ale;
  _WRITE_INT_VAL("methode_ALE");

  dummy_int = (int)cs_glob_vof_parameters->vof_model;
  _WRITE_INT_VAL("vof");

  cs_log_printf(CS_LOG_DEFAULT, " End writing the dimensions and options\n");

  /* ----------------- */
  /* Writing variables */
  /* ----------------- */

  /* Point de reference pour la pression totale
   * On n'ecrit que si XYZP0 a ete specifie par l'utilisateur ou
   * calcule a partir de faces de sorties ou de Dirichlet
   */
  cs_fluid_properties_t *cgfp = cs_get_glob_fluid_properties();
  dummy_int = cgfp->ixyzp0;
  if (cgfp->ixyzp0== 1) {
    cs_restart_write_section(r,"ref_presstot01", 0 , 3,
                             CS_TYPE_cs_real_t,
                             cgfp->xyzp0);
  }

  /* The physical variables here below are required for the low-Mach algorithm */
  if (cs_glob_velocity_pressure_model->idilat == 3 ||
      cgfp->ipthrm == 1) {
    cs_real_t dummy_real = cgfp->ro0;
    _WRITE_REAL_VAL("ro001");

    dummy_real = cgfp->pther;
    _WRITE_REAL_VAL("pther01");
  }

  cs_restart_write_linked_fields(r, "diffusivity_id", nullptr);

  cs_log_printf(CS_LOG_DEFAULT, " End writing the physical properties\n");

  /* Time step */
  if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) {
    cs_real_t dummy_real = CS_F_(dt)->val[0];
    _WRITE_REAL_VAL("dt_variable_temps");
  }

  cs_log_printf(CS_LOG_DEFAULT, " End writing the time step\n");

  /* Mass fluxes */
  cs_restart_write_linked_fields(r, "inner_mass_flux_id", nullptr);
  cs_restart_write_linked_fields(r, "boundary_mass_flux_id", nullptr);

  /* Boundary condition coefficients */
  cs_restart_write_bc_coeffs(r);

  /* Source terms when extrapolated */
  int n_written_fields =
    cs_restart_write_linked_fields(r, "source_term_prev_id", nullptr);

  if (n_written_fields > 0)
    cs_log_printf(CS_LOG_DEFAULT, " End writing the source terms\n");

  /* Time moments */
  cs_time_moment_restart_write(r);

  /* Wall distance */
  // Currently not done

  /* Wall temperature associated to the condensation model
   * with or without the 1D thermal model tag1D
   */
  cs_wall_condensation_t *wco = cs_get_glob_wall_condensation();
  if (wco->icondb == 0) {
    cs_real_t *tmp =  nullptr;
    if (wco->nztag1d == 1) {
      cs_wall_cond_1d_thermal_t *wco1d = cs_get_glob_wall_cond_1d_thermal();

      BFT_MALLOC(tmp, wco1d->znmurx * cs_glob_mesh->n_b_faces, cs_real_t);
      cs_array_real_fill_zero(wco1d->znmurx * cs_glob_mesh->n_b_faces,
                              tmp);

      for (cs_lnum_t e_id = 0; e_id < wco->nfbpcd; e_id++) {
        cs_lnum_t f_id = wco->ifbpcd[e_id];
        cs_lnum_t z_id = wco->izzftcd[e_id];
        for (int i = 0; i < wco1d->znmur[z_id]; i++)
          tmp[f_id * wco1d->znmurx + i] =
            wco1d->ztmur[e_id * wco1d->znmurx + i];
      }

      cs_restart_write_section(r,
                               "tmur_bf_prev",
                               3,
                               wco1d->znmurx,
                               CS_TYPE_cs_real_t,
                               tmp);
    }
    else {
      BFT_MALLOC(tmp, cs_glob_mesh->n_b_faces, cs_real_t);
      cs_array_real_fill_zero(cs_glob_mesh->n_b_faces, tmp);

      for (cs_lnum_t e_id = 0; e_id < wco->nfbpcd; e_id++) {
        cs_lnum_t f_id = wco->ifbpcd[e_id];
        cs_lnum_t z_id = wco->izzftcd[e_id];
        tmp[f_id] = wco->ztpar[z_id];
      }

      cs_restart_write_section(r,
                               "tpar_bf_prev",
                               3,
                               1,
                               CS_TYPE_cs_real_t,
                               tmp);
    }

    BFT_FREE(tmp);
  }

  /* ALE */
  if (cs_glob_ale != CS_ALE_NONE) {
    cs_ale_restart_write(r);
    cs_mobile_structures_restart_write(r);

    cs_log_printf(CS_LOG_DEFAULT, " End writing the ALE data\n");
  }

  /* ---------------------------------------- */
  /* Combustion related fields and structures */
  /* ---------------------------------------- */

  /* 3 points model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0) {
    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;
    cs_real_t dummy_real = cm->hinfue;
    _WRITE_REAL_VAL("hinfue_cod3p");

    dummy_real = cm->hinoxy;
    _WRITE_REAL_VAL("hinoxy_cod3p");

    dummy_real = cm->tinfue;
    _WRITE_REAL_VAL("tinfue_cod3p");

    dummy_real = cm->tinoxy;
    _WRITE_REAL_VAL("tinoxy_cod3p");

    // numéro des zones
    cs_restart_write_section(r,
                             "num_zone_fb_cod3p",
                             3,
                             1,
                             CS_TYPE_cs_real_t,
                             cs_glob_bc_pm_info->izfppp);

    // Entree Fuel (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientfu_zone_bord_cod3p",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientfu + 1);

    // Entree oxydant (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientox_zone_bord_cod3p",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientox + 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  " End writing combustion information (COD3P)\n");
  }

  /* SLFM model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0) {
    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    cs_real_t dummy_real = cm->hinfue;
    _WRITE_REAL_VAL("hinfue_slfm");

    dummy_real = cm->hinoxy;
    _WRITE_REAL_VAL("hinoxy_slfm");

    dummy_real = cm->tinfue;
    _WRITE_REAL_VAL("tinfue_slfm");

    dummy_real = cm->tinoxy;
    _WRITE_REAL_VAL("tinoxy_slfm");

    // numéro des zones
    cs_restart_write_section(r,
                             "num_zone_fb_slfm",
                             3,
                             1,
                             CS_TYPE_cs_real_t,
                             cs_glob_bc_pm_info->izfppp);

    // Entree Fuel (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientfu_zone_bord_slfm",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientfu + 1);

    // Entree oxydant (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientox_zone_bord_slfm",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientox + 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  " End writing combustion information (SLFM)\n");
  }

  /* EBU model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0) {
    cs_real_t dummy_real = cs_glob_bc_pm_info->tgf;
    _WRITE_REAL_VAL("temperature_gaz_frais_ebu");

    dummy_real = cs_glob_bc_pm_info->frmel;
    _WRITE_REAL_VAL("frmel_ebu");

    // numéro des zones
    cs_restart_write_section(r,
                             "num_zone_fb_ebu",
                             3,
                             1,
                             CS_TYPE_cs_real_t,
                             cs_glob_bc_pm_info->izfppp);

    // Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientgb_zone_bord_ebu",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientgb + 1);

    // Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientgf_zone_bord_ebu",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientgf + 1);

    // FMENT (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "fment_zone_bord_ebu",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->fment + 1);

    // TKENT (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "tkent_zone_bord_ebu",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->tkent + 1);

    cs_log_printf(CS_LOG_DEFAULT,
                  " End writing the combustion information (EBU)\n");
  }

  /* LWC model */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_LW] >= 0) {
    const cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    cs_real_t dummy_real = cm->fmin_lwc;
    _WRITE_REAL_VAL("fmin_lwc");

    dummy_real = cm->fmax_lwc;
    _WRITE_REAL_VAL("fmax_lwc");

    dummy_real = cm->hmin_lwc;
    _WRITE_REAL_VAL("hmin_lwc");

    dummy_real = cm->hmax_lwc;
    _WRITE_REAL_VAL("hmax_lwc");

    // numéro des zones
    cs_restart_write_section(r,
                             "num_zone_fb_lwc",
                             3,
                             1,
                             CS_TYPE_cs_real_t,
                             cs_glob_bc_pm_info->izfppp);

    // Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientgb_zone_bord_lwc",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientgb + 1);

    // Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "ientgf_zone_bord_lwc",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->ientgf + 1);

    // FMENT (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "fment_zone_bord_lwc",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->fment + 1);

    // TKENT (si ce n'est pas NOZPPM, erreur)
    cs_restart_write_section(r,
                             "tkent_zone_bord_lwc",
                             0,
                             CS_MAX_BC_PM_ZONE_NUM,
                             CS_TYPE_int,
                             cs_glob_bc_pm_info->tkent + 1);

    cs_log_printf(CS_LOG_DEFAULT, " End writing combustion information (LWC)\n");
  }

  /* Pulverized coal combustion model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0) {
    cs_coal_model_t  *cm = cs_glob_coal_model;
    const char *_prefix = "masse_volumique_charbon";
    int _len = strlen(_prefix) + 3;

    for (int i = 0; i < cm->n_coals; i++) {
      char _rub[64] = "";

      if (i < 100) // Hard coded limit, is it still needed ?
        snprintf(_rub, _len, "%s%02d", _prefix, i);
      else
        snprintf(_rub, _len, "%sXX", _prefix);

      cs_real_t dummy_real = cm->rhock[i];
      _WRITE_REAL_VAL(_rub);
    }

    cs_log_printf(CS_LOG_DEFAULT, " End writing combustion information (CP)\n");
  }

  /* ------------------------------------------------------- */
  /* Grandeurs complementaires pour les versions electriques */
  /* ------------------------------------------------------- */

  if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0 ||
      cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0) {

    cs_elec_option_t *ce = cs_get_glob_elec_option();

    if (ce->ielcor == 1) {
      cs_real_t dummy_real = ce->pot_diff;
      _WRITE_REAL_VAL("ddpot_recalage_arc_elec");

      dummy_real = ce->elcou;
      _WRITE_REAL_VAL("elcou_recalage_arc_elec");

      if (cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0) {
        dummy_real = ce->coejou;
        _WRITE_REAL_VAL("coeff_recalage_joule");
      }

      cs_log_printf(CS_LOG_DEFAULT, " End writing the electric information\n");
    }
  }

  /* Write fields */
  cs_restart_write_fields(r, CS_RESTART_AUXILIARY);

  /* Close file */
  cs_restart_destroy(&r);

  cs_log_printf(CS_LOG_DEFAULT, " End writing\n");
}

/*----------------------------------------------------------------------------
 * Read auxiliary checkpoint file
 *----------------------------------------------------------------------------*/

static void
_read_auxiliary_checkpoint(cs_map_name_to_id_t *old_field_map)
{
  cs_log_printf(CS_LOG_DEFAULT, "      READING THE AUXILIARY RESTART FILE\n");

  cs_restart_t *r = cs_restart_create("auxiliary.csc",
                                      nullptr,
                                      CS_RESTART_MODE_READ);

  cs_log_printf(CS_LOG_DEFAULT, " Start reading\n");

  /* Check restart file */
  cs_real_t dummy_real = -99999.;

  int retval
    = cs_restart_read_section_compat(r,
                                     "code_saturne:checkpoint:auxiliary:version",
                                     "version_fichier_suite_auxiliaire",
                                     0, 1,
                                     CS_TYPE_int,
                                     &dummy_real);

  if (retval != CS_RESTART_SUCCESS) {
    cs_log_warning(_("STOP READING THE AUXILIARY RESTART FILE\n\n"
                     "  WRONG FILE TYPE\n\n"
                     "The file \"auxiliary.csc\" does not look like a proper\n"
                     "  auxiliary restart file.\n\n"
                     "The calculation cannot be executed.\n\n"
                     "Please make sure the file used as a restart file\n"
                     "  actually is a correct auxiliary restart file.\n"
                     "If necessary, it is possible to deactivate the reading\n"
                     "  of the auxiliary restart file by ILEAUX = 0.\n"));
    cs_exit(EXIT_FAILURE);
  }

  /* Check base locations */
  bool match_cell, match_i_face, match_b_face, match_vertex;
  cs_restart_check_base_location(r,
                                 &match_cell,
                                 &match_i_face,
                                 &match_b_face,
                                 &match_vertex);

  if (!match_cell) {
    cs_log_warning(_("STOP READING THE AUXILIARY RESTART FILE\n\n"
                     "  INCOHERENT PREVIOUS NAD ACTUAL DATA\n\n"
                     "The number of cells was modified\n\n"
                     "The run can not be executed.\n\n"
                     "Verify that the restart file used corresponds to\n"
                     "    the present case.\n"
                     "If necessary, it is possible to deactivate the reading\n"
                     "  of the auxiliary restart file by ILEAUX = 0.\n"));
    cs_exit(EXIT_FAILURE);
  }

  bool face_states_[2] = {match_i_face, match_b_face};
  const char *face_names[2] = {"internal", "boundary"};

  for (int i = 0; i < 2; i++) {
    if (!face_states_[i]) {
      cs_log_warning
        (_("STOP READING THE AUXILIARY RESTART FILE\n\n"
           "  PREVIOUS and PRESENT INPUT DATA ARE DIFFERENT\n\n"
           "The number of %s faces has been modified\n\n"
           "The run can continue but the data on the\n"
           "  %s faces will not be reread in the suite file.\n"
           "They will be initialized by the default values.\n\n"
           " This situation can occur when the restart file\n"
           "  originates from a run using different options\n"
           "  to join the grids or when the periodicity boundary\n"
           "  conditions have been modified.\n"
           " This situation can also be generated when the\n"
           "  run is conducted on a different machine\n"
           "  in which case the precision of the machine modifies\n"
           "  the number of faces generated when joinning the grids.\n\n"
           " Finally, this situation can be due to the fact that\n"
           "  the auxiliary restart file does not correspond to\n"
           "  the present case.\n\n"
           "Verify that the auxiliary restart file being used\n"
           "  corresponds to the present case.\n\n"
           " The run will continue...\n"),
         face_names[i], face_names[i]);
    }
  }

  /* ALE */
  int dummy_int = 1;
  retval = _READ_INT_VAL("methode_ALE");

  int ale_aux_id = (retval == CS_RESTART_SUCCESS) ? dummy_int : 0;

  if (retval != CS_RESTART_SUCCESS) {
    if (cs_glob_ale != CS_ALE_NONE)
      cs_log_warning
        (_("ERROR WHILE READING THE AUXILIARY RESTART FILE\n\n"
           "  ERROR WHEN READING THE INDICATOR OF THE ALE METHOD\n\n"
           "It is possible that the file read corresponds to an old\n"
           "  version of code_saturne, without the ALE method.\n"
           "The run will be executed with reinitialising all ALE data.\n\n"
           "Verify that the restart file used has not been damaged.\n"));
  }

  if (cs_glob_ale_need_init == -999) {
    if (cs_glob_ale != CS_ALE_NONE && ale_aux_id > 0)
      cs_glob_ale_need_init = 0;
    else if (cs_glob_ale != CS_ALE_NONE)
      cs_glob_ale_need_init = 1;
    else
      cs_glob_ale_need_init = 0;
  }

  /* VoF */
  dummy_int = 0;
  retval = _READ_INT_VAL("vof");

  int vof_aux_id = (retval == CS_RESTART_SUCCESS) ? dummy_int : 0;

  if (retval != CS_RESTART_SUCCESS) {
    if (cs_glob_vof_parameters->vof_model & CS_VOF_ENABLED)
      cs_log_warning
        (_("ERROR WHILE READING THE AUXILIARY RESTART FILE\n\n"
           "  ERROR WHEN READING THE INDICATOR OF THE VOF MODEL\n\n"
           "It is possible that the file read corresponds to an old\n"
           "  version of code_saturne, without the VOF model.\n"
           "The run will be executed with reinitialising all\n"
           "  VOF model data.\n\n"
           "Verify that the restart file used has not been damaged.\n"));
  }

  cs_log_printf(CS_LOG_DEFAULT, " Finished reading options.\n");

  /* ------------------- */
  /* Physical properties */
  /* ------------------- */

  cs_fluid_properties_t *cgfp = cs_get_glob_fluid_properties();

  /* Pressure reference point */
  if (cgfp->ixyzp0 == -1) {
    retval = cs_restart_read_section(r, "ref_presstot01", 0, 3,
                                     CS_TYPE_cs_real_t,
                                     cgfp->xyzp0);
    if (retval == CS_RESTART_SUCCESS) {
      cgfp->ixyzp0 = 1;
      cs_log_printf(CS_LOG_DEFAULT,
                    _("   Apdatation of the reference point for the"
                      " total pressure\n"
                      "       by reading the restart file\n"
                      "    XYZP0 = %14.5e, %14.5e, %14.5e \n"),
                    cgfp->xyzp0[0],
                    cgfp->xyzp0[1],
                    cgfp->xyzp0[2]);
    }
  }

  /* Here the physical variables below are required for the low-Mach algorithm */

  if (cs_glob_velocity_pressure_model->idilat == 3 ||
      cgfp->ipthrm == 1) {
    dummy_real = 0.;
    retval = _READ_REAL_VAL("ro001");
    if (retval == CS_RESTART_SUCCESS)
      cgfp->ro0 = dummy_real;

    retval = _READ_REAL_VAL("pther01");
    if (retval == CS_RESTART_SUCCESS)
      cgfp->pther = dummy_real;
  }

  /* Density */
  if (cgfp->irovar == 1 ||
      (cs_glob_vof_parameters->vof_model & CS_VOF_ENABLED && vof_aux_id > 0)) {
    int read_rho_ok = 1;
    if (cs_restart_get_field_read_status(CS_F_(rho)->id) == 0)
      read_rho_ok = 0;

    if (match_b_face) {
      if (cs_restart_get_field_read_status(CS_F_(rho_b)->id) == 0)
        read_rho_ok = 0;
    }

    if (match_b_face && read_rho_ok == 1)
      cs_parameters_set_init_state_on(1);
  }
  else
    cs_parameters_set_init_state_on(1); // 1 is density


  /* Viscosity */
  int key_t_ext_id = cs_field_key_id_try("time_extrapolated");

  if (cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id) > 0 ||
      (cs_glob_vof_parameters->vof_model & CS_VOF_ENABLED && vof_aux_id >= 0)) {
    int read_mu_t_ok = 1;

    if (cs_restart_get_field_read_status(CS_F_(mu)->id) == 0)
      read_mu_t_ok = 0;

    if (cs_restart_get_field_read_status(CS_F_(mu_t)->id) == 0)
      read_mu_t_ok = 0;

    if (read_mu_t_ok == 1)
      cs_parameters_set_init_state_on(0); // 0 is viscosity
  }

  /* Cp */
  if (cgfp->icp > -1) {
    if (cs_field_get_key_int(CS_F_(cp), key_t_ext_id) > 0 ||
        cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 1) {

      if (cs_restart_get_field_read_status(CS_F_(cp)->id) == 1)
        cs_parameters_set_init_state_on(2); // 2 is cp id
    }
  }

  /* Read diffusivities if needed */
  cs_restart_read_linked_fields(r, old_field_map, "diffusivity_id", nullptr);

  cs_log_printf(CS_LOG_DEFAULT, " Finished reading physical properties.\n");

  /* -------------------- */
  /* Time step quantities */
  /* -------------------- */

  dummy_int = 0;
  retval = _READ_INT_VAL("indic_dt_variable");

  if (retval != CS_RESTART_SUCCESS) {
    cs_log_warning("Error while reading th time stepping mode\n");
    cs_exit(EXIT_FAILURE);
  }

  if (cs_glob_time_step_options->idtvar != dummy_int) {
    cs_log_warning(_("Warning: computation was restarted with time stepping\n"
                     "option idtvar = %d while the previous run used\n"
                     "option idtvar = %d.\n"),
                   cs_glob_time_step_options->idtvar,
                   dummy_int);
  }
  else if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_ADAPTIVE) {
    retval = _READ_REAL_VAL("dt_variable_temps");
    if (retval == CS_RESTART_SUCCESS)
      cs_array_real_set_scalar(cs_glob_mesh->n_cells,
                               dummy_real,
                               CS_F_(dt)->val);
  }
  else if (cs_glob_time_step_options->idtvar == CS_TIME_STEP_LOCAL) {
    retval = cs_restart_read_field_vals(r, CS_F_(dt)->id, 0);
    if (retval != CS_RESTART_SUCCESS)
      cs_log_warning(_("Reading time step field values failed.\n"
                       "Continuing with default values.\n"));
  }

  /* ----------- */
  /* Mass fluxes */
  /* ----------- */

  if (match_i_face || match_b_face) {
    /* Read fluxes */
    cs_restart_read_linked_fields(r, old_field_map,
                                  "inner_mass_flux_id", nullptr);
    cs_restart_read_linked_fields(r, old_field_map,
                                  "boundary_mass_flux_id", nullptr);

    /* Initiliaze void fraction if needed */
    if (cs_glob_vof_parameters->vof_model & CS_VOF_ENABLED && vof_aux_id < 0) {
      const int kimasf = cs_field_key_id("inner_mass_flux_id");

      cs_field_t *mflux =
        cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kimasf));
      cs_field_t *vof_flux =
        cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kimasf));

      const cs_real_t oo_rho1 = 1. / cs_glob_vof_parameters->rho1;

      /* Use wscalar option, with weight being mflux */
      cs_array_real_set_wscalar(cs_glob_mesh->n_i_faces, oo_rho1,
                                mflux->val, vof_flux->val);

      if (mflux->n_time_vals > 1)
        cs_array_real_set_wscalar(cs_glob_mesh->n_i_faces, oo_rho1,
                                  mflux->val_pre, vof_flux->val_pre);

      /* Boundary values */
      const int kbmasf = cs_field_key_id("boundary_mass_flux_id");

      cs_field_t *b_mflux
        = cs_field_by_id(cs_field_get_key_int(CS_F_(vel), kbmasf));
      cs_field_t *b_vof_flux
        = cs_field_by_id(cs_field_get_key_int(CS_F_(void_f), kbmasf));

      cs_array_real_set_wscalar(cs_glob_mesh->n_b_faces, oo_rho1,
                                b_mflux->val, b_vof_flux->val);

      if (mflux->n_time_vals > 1)
        cs_array_real_set_wscalar(cs_glob_mesh->n_b_faces, oo_rho1,
                                  b_mflux->val_pre, b_vof_flux->val_pre);
    }
  }

  /* ------------------- */
  /* Boundary conditions */
  /* ------------------- */

  if (match_b_face) {
    cs_restart_read_bc_coeffs(r);
    cs_log_printf(CS_LOG_DEFAULT, " Finished reading boundary conditions.\n");
  }

  /* ------------ */
  /* Source terms */
  /* ------------ */

  cs_restart_read_linked_fields(r, old_field_map, "source_term_prev_id", nullptr);
  cs_log_printf(CS_LOG_DEFAULT, " Finished reading source terms.\n");

  /* ------------ */
  /* Time moments */
  /* ------------ */

  cs_time_moment_restart_read(r);

  /* -------------------------------------------------------------------- */
  /* Wall temperature associated to the condensation model with or without
   * the 1D thermal model tag1D */
  /* -------------------------------------------------------------------- */

  cs_wall_condensation_t *wco = cs_get_glob_wall_condensation();
  if (wco->icondb == 0 ) {
    cs_real_t *tmp =  nullptr;
    if (wco->nztag1d == 1) {
      cs_wall_cond_1d_thermal_t *wco1d = cs_get_glob_wall_cond_1d_thermal();

      BFT_MALLOC(tmp, wco1d->znmurx * cs_glob_mesh->n_b_faces, cs_real_t);
      retval = cs_restart_read_section(r,
                                       "tmur_bf_prev",
                                       3,
                                       wco1d->znmurx,
                                       CS_TYPE_cs_real_t,
                                       tmp);
      if (retval == CS_RESTART_SUCCESS) {
        for (cs_lnum_t e_id = 0; e_id < wco->nfbpcd; e_id++) {
          cs_lnum_t f_id = wco->ifbpcd[e_id];
          cs_lnum_t z_id = wco->izzftcd[e_id];
          for (int i = 0; i < wco1d->znmur[z_id]; i++)
            wco1d->ztmur[e_id * wco1d->znmurx + i]
              = tmp[f_id * wco1d->znmurx + i];
        }
      }

    }
    else {
      BFT_MALLOC(tmp, cs_glob_mesh->n_b_faces, cs_real_t);

      retval = cs_restart_read_section(r,
                                       "tpar_bf_prev",
                                       3,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       tmp);

      if (retval == CS_RESTART_SUCCESS) {
        for (cs_lnum_t e_id = 0; e_id < wco->nfbpcd; e_id++) {
          cs_lnum_t f_id = wco->ifbpcd[e_id];
          cs_lnum_t z_id = wco->izzftcd[e_id];

          wco->ztpar[z_id] = tmp[f_id];
        }
      }
    }

    BFT_FREE(tmp);
  }

  /* ------------------- */
  /* ALE vertex movement */
  /* ------------------- */

  if (cs_glob_ale != CS_ALE_NONE && ale_aux_id > 0) {
    cs_ale_restart_read(r);
    cs_mobile_structures_restart_read(r);

    cs_log_printf(CS_LOG_DEFAULT, " Finished reading ALE information.\n");
  }

  /* ----------------------- */
  /* Combustion related data */
  /* ----------------------- */

  /* 3 Points model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0) {
    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    retval = _READ_REAL_VAL("hinfue_cod3p");
    if (retval == CS_RESTART_SUCCESS)
      cm->hinfue = dummy_real;

    retval = _READ_REAL_VAL("hinoxy_cod3p");
    if (retval == CS_RESTART_SUCCESS)
      cm->hinoxy = dummy_real;

    retval = _READ_REAL_VAL("tinfue_cod3p");
    if (retval == CS_RESTART_SUCCESS)
      cm->tinfue = dummy_real;

    retval = _READ_REAL_VAL("tinoxy_cod3p");
    if (retval == CS_RESTART_SUCCESS)
      cm->tinoxy = dummy_real;

    /* boundary faces data is only read if the number did not change */
    if (match_b_face) {

      // numéro des zones
      retval = cs_restart_read_section(r,
                                       "num_zone_fb_cod3p",
                                       3,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       cs_glob_bc_pm_info->izfppp);

      // Entree Fuel (si ce n'est pas NOZPPM, erreur)
      retval = cs_restart_read_section(r,
                                       "ientfu_zone_bord_cod3p",
                                       0,
                                       CS_MAX_BC_PM_ZONE_NUM,
                                       CS_TYPE_int,
                                       cs_glob_bc_pm_info->ientfu + 1);

      // Entree oxydant (si ce n'est pas NOZPPM, erreur)
      retval += cs_restart_read_section(r,
                                        "ientox_zone_bord_cod3p",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->ientox + 1);

      /* We check that reading call worked for both ientfu and ientox.
       * if not, izfppp is reinitialized.
       */
      if (retval != CS_RESTART_SUCCESS) {
        for (cs_lnum_t f_id = 0; f_id < cs_glob_mesh->n_b_faces; f_id++)
          cs_glob_bc_pm_info->izfppp[f_id] = 0;
      }
    }
  }

  /* SLFM model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0) {
    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    retval = _READ_REAL_VAL("hinfue_slfm");
    if (retval == CS_RESTART_SUCCESS)
      cm->hinfue = dummy_real;

    retval = _READ_REAL_VAL("hinoxy_slfm");
    if (retval == CS_RESTART_SUCCESS)
      cm->hinoxy = dummy_real;

    retval = _READ_REAL_VAL("tinfue_slfm");
    if (retval == CS_RESTART_SUCCESS)
      cm->tinfue = dummy_real;

    retval = _READ_REAL_VAL("tinoxy_slfm");
    if (retval == CS_RESTART_SUCCESS)
      cm->tinoxy = dummy_real;

    if (match_b_face) {
      // numéro des zones
      retval = cs_restart_read_section(r,
                                       "num_zone_fb_slfm",
                                       3,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       cs_glob_bc_pm_info->izfppp);

      // Entree Fuel (si ce n'est pas NOZPPM, erreur)
      retval = cs_restart_read_section(r,
                                       "ientfu_zone_bord_cod3p",
                                       0,
                                       CS_MAX_BC_PM_ZONE_NUM,
                                       CS_TYPE_int,
                                       cs_glob_bc_pm_info->ientfu + 1);

      // Entree oxydant (si ce n'est pas NOZPPM, erreur)
      retval += cs_restart_read_section(r,
                                        "ientox_zone_bord_cod3p",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->ientox + 1);

      /* We check that reading call worked for both ientfu and ientox.
       * if not, izfppp is reinitialized.
       */
      if (retval != CS_RESTART_SUCCESS) {
        for (cs_lnum_t f_id = 0; f_id < cs_glob_mesh->n_b_faces; f_id++)
          cs_glob_bc_pm_info->izfppp[f_id] = 0;
      }
    }
  }

  /* EBU Model */

  if (cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0) {
    retval = _READ_REAL_VAL("temperature_gaz_frais_ebu");
    if (retval == CS_RESTART_SUCCESS)
      cs_glob_bc_pm_info->tgf = dummy_real;

    retval = _READ_REAL_VAL("frmel_ebu");
    if (retval == CS_RESTART_SUCCESS)
      cs_glob_bc_pm_info->frmel = dummy_real;

    if (match_b_face) {
      // numéro des zones
      retval = cs_restart_read_section(r,
                                       "num_zone_fb_ebu",
                                       3,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       cs_glob_bc_pm_info->izfppp);

      // Entree gas brule
      retval = cs_restart_read_section(r,
                                       "ientgb_zone_bord_ebu",
                                       0,
                                       CS_MAX_BC_PM_ZONE_NUM,
                                       CS_TYPE_int,
                                       cs_glob_bc_pm_info->ientgb + 1);

      // Entree gas frais
      retval += cs_restart_read_section(r,
                                        "ientgf_zone_bord_ebu",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->ientgf + 1);

      // fment
      retval += cs_restart_read_section(r,
                                        "fment_zone_bord_ebu",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->fment + 1);
      // tkent
      retval += cs_restart_read_section(r,
                                        "tkent_zone_bord_ebu",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->tkent + 1);

      if (retval != CS_RESTART_SUCCESS) {
        for (cs_lnum_t f_id = 0; f_id < cs_glob_mesh->n_b_faces; f_id++)
          cs_glob_bc_pm_info->izfppp[f_id] = 0;
      }
    }
  }

  /* LWC */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_LW] >= 0) {
    cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

    retval = _READ_REAL_VAL("fmin_lwc");
    if (retval == CS_RESTART_SUCCESS)
      cm->fmin_lwc = dummy_real;

    retval = _READ_REAL_VAL("fmax_lwc");
    if (retval == CS_RESTART_SUCCESS)
      cm->fmax_lwc = dummy_real;

    retval = _READ_REAL_VAL("hmin_lwc");
    if (retval == CS_RESTART_SUCCESS)
      cm->hmin_lwc = dummy_real;

    retval = _READ_REAL_VAL("hmax_lwc");
    if (retval == CS_RESTART_SUCCESS)
      cm->hmax_lwc = dummy_real;

    if (match_b_face) {
      // numéro des zones
      retval = cs_restart_read_section(r,
                                       "num_zone_fb_lwc",
                                       3,
                                       1,
                                       CS_TYPE_cs_real_t,
                                       cs_glob_bc_pm_info->izfppp);

      // Entree gas brule
      retval = cs_restart_read_section(r,
                                       "ientgb_zone_bord_lwc",
                                       0,
                                       CS_MAX_BC_PM_ZONE_NUM,
                                       CS_TYPE_int,
                                       cs_glob_bc_pm_info->ientgb + 1);

      // Entree gas frais
      retval += cs_restart_read_section(r,
                                        "ientgf_zone_bord_lwc",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->ientgf + 1);

      // fment
      retval += cs_restart_read_section(r,
                                        "fment_zone_bord_lwc",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->fment + 1);
      // tkent
      retval += cs_restart_read_section(r,
                                        "tkent_zone_bord_lwc",
                                        0,
                                        CS_MAX_BC_PM_ZONE_NUM,
                                        CS_TYPE_int,
                                        cs_glob_bc_pm_info->tkent + 1);

      if (retval != CS_RESTART_SUCCESS) {
        for (cs_lnum_t f_id = 0; f_id < cs_glob_mesh->n_b_faces; f_id++)
          cs_glob_bc_pm_info->izfppp[f_id] = 0;
      }
    }
  }

  /* Pulverized coal */
  if (cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0) {
    cs_coal_model_t  *cm = cs_glob_coal_model;
    const char *_prefix = "masse_volumique_charbon";
    int _len = strlen(_prefix) + 3;

    for (int i = 0; i < cm->n_coals; i++) {
      char _rub[64] = "";

      if (i < 100) // Hard coded limit, is it still needed ?
        snprintf(_rub, _len, "%s%02d", _prefix, i);
      else
        snprintf(_rub, _len, "%sYY", _prefix);

      retval = _READ_REAL_VAL(_rub);
      if (retval == CS_RESTART_SUCCESS)
        cm->rhock[i] = dummy_real;
    }
  }

  /* Electric arc model */

  if (cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] > 0 ||
      cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0) {

    cs_elec_option_t *ce = cs_get_glob_elec_option();

    if (ce->ielcor == 1) {
      retval = _READ_REAL_VAL("ddpot_recalage_arc_elec");
      if (retval == CS_RESTART_SUCCESS)
        ce->pot_diff = dummy_real;

      retval = _READ_REAL_VAL("elcou_recalage_arc_elec");
      if (retval == CS_RESTART_SUCCESS)
        ce->elcou = dummy_real;

      if (cs_glob_physical_model_flag[CS_JOULE_EFFECT] > 0) {
        retval = _READ_REAL_VAL("coeff_recalage_joule");
        if (retval == CS_RESTART_SUCCESS)
          ce->coejou = dummy_real;
      }
    }

    cs_log_printf(CS_LOG_DEFAULT, " Finished reading electric information.\n");
  }


  /* -------------------------------- */
  /* Read fields based on restart key */
  /* -------------------------------- */

  cs_restart_read_fields(r, CS_RESTART_AUXILIARY);

  /* ---------- */
  /* Close file */
  /* ---------- */

  cs_restart_destroy(&r);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write main and auxiliary checkpoint files
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_main_and_aux_write(void)
{
  /* Write main checkpoint file */
  _write_main_checkpoint();

  /* Auxiliary file */
  if (cs_glob_restart_auxiliary->write_auxiliary == 1)
    _write_auxiliary_checkpoint();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read main and auxiliary checkpoint files
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_main_and_aux_read(void)
{
  /* Read main checkpoint file */
  cs_map_name_to_id_t *old_field_map = _read_main_checkpoint();

  /* Read auxiliary file */
  if (cs_glob_restart_auxiliary->read_auxiliary == 1)
    _read_auxiliary_checkpoint(old_field_map);

  /* Cleanup */
  cs_map_name_to_id_destroy(&old_field_map);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
