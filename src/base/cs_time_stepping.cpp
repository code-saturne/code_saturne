/*============================================================================
 * Main time loop.
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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_1d_wall_thermal.h"
#include "cs_1d_wall_thermal_check.h"
#include "cs_ale.h"
#include "cs_at_data_assim.h"
#include "cs_atmo.h"
#include "cs_boundary_conditions.h"
#include "cs_boundary_conditions_set_coeffs.h"
#include "cs_cdo_main.h"
#include "cs_cf_boundary_conditions.h"
#include "cs_control.h"
#include "cs_coupling.h"
#include "cs_domain_op.h"
#include "cs_domain_setup.h"
#include "cs_field_pointer.h"
#include "cs_gas_mix.h"
#include "cs_gui.h"
#include "cs_turbulence_htles.h"
#include "cs_ibm.h"
#include "cs_initialize_fields.h"
#include "cs_lagr.h"
#include "cs_lagr_lec.h"
#include "cs_les_balance.h"
#include "cs_les_inflow.h"
#include "cs_log_iteration.h"
#include "cs_mesh.h"
#include "cs_mesh_save.h"
#include "cs_mobile_structures.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_physical_properties_default.h"
#include "cs_porosity_from_scan.h"
#include "cs_porous_model.h"
#include "cs_post.h"
#include "cs_post_default.h"
#include "cs_prototypes.h"
#include "cs_rad_transfer.h"
#include "cs_rad_transfer_restart.h"
#include "cs_resource.h"
#include "cs_restart.h"
#include "cs_restart_default.h"
#include "cs_restart_main_and_aux.h"
#include "cs_restart_map.h"
#include "cs_sat_coupling.h"
#include "cs_solve_all.h"
#include "cs_runaway_check.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_timer_stats.h"
#include "cs_turbomachinery.h"
#include "cs_turbulence_bc.h"
#include "cs_turbulence_model.h"
#include "cs_volume_mass_injection.h"
#include "cs_vof.h"
#include "cs_wall_condensation.h"
#include "cs_wall_condensation_1d_thermal.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_time_stepping.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * External function prototypes
 *============================================================================*/

/* Bindings to Fortran routines */

void
cs_f_boundary_conditions_init(void);

void
cs_f_init_chemistry_reacnum(void);

void
cs_f_atmo_models_boundary_conditions_map(void);

void
cs_f_combustion_models_boundary_conditions_map(void);

void
cs_f_initialization_variables(void);

void
cs_f_atmsol(void);

void
cs_f_user_extra_operations_wrapper(cs_real_t dt[]);

void
cs_f_finalize_meteo(void);

void
cs_f_finalize_imbrication(void);

void
cs_f_finalize_chemistry(void);

void
cs_f_finalize_steady_laminar_flamelet_library(void);

void
cs_f_pp_models_bc_map(void);

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_time_stepping.c
        Main time loop.
*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Main time loop.
 */
/*----------------------------------------------------------------------------*/

void
cs_time_stepping(void)
{
  cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;
  cs_time_step_t *ts = cs_get_glob_time_step();
  const cs_turb_model_t *turb_model = cs_get_glob_turb_model();

  int idtvar = cs_glob_time_step_options->idtvar;
  cs_turbomachinery_model_t iturbo = cs_turbomachinery_get_model();

  /* Initialization
     -------------- */

  /* Define timer stats based on options */

  int lagr_stats_id = -1;

  if (cs_glob_lagr_time_scheme->iilagr > 0) {
    lagr_stats_id = cs_timer_stats_create("stages",
                                          "lagrangian_stage",
                                          "Lagrangian Module");
    cs_timer_stats_create("lagrangian_stage",
                          "particle_displacement_stage",
                          "particle displacement");
  }

  /* End of modules initialization
     ----------------------------- */

  cs_turb_init_ref_quantities();

  /* Zone definition for head-loss, mass sources term,
     condensation sources term and 1D-wall module
     ------------------------------------------------- */

  /* Mass source terms
     ----------------- */

  /* Total number of cells with mass source term */

  cs_lnum_t ncetsm = 0;

  cs_volume_mass_injection_get_arrays(nullptr,
                                      &ncetsm,
                                      nullptr,
                                      nullptr,
                                      nullptr,
                                      nullptr);

  cs_lnum_t n_cells_mst_tot = ncetsm;

  if (cs_glob_rank_id > -1)
    cs_parall_sum(1, CS_LNUM_TYPE, &n_cells_mst_tot);

  if (n_cells_mst_tot > 0) {
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n\nMass source terms process activated on a total"
         " of %llu cells.\n"), (unsigned long long)n_cells_mst_tot);

    cs_log_separator(CS_LOG_DEFAULT);
  }

  /* Condensation mass source terms
     ------------------------------ */

  cs_user_wall_condensation(1);

  /* Total number of cells with condensation source term */
  const cs_lnum_t nftcdt = cs_glob_wall_condensation->nfbpcd;
  cs_gnum_t n_cells_cst_tot = nftcdt;

  if (cs_glob_rank_id > -1)
    cs_parall_sum(1, CS_LNUM_TYPE, &n_cells_cst_tot);

  if (n_cells_cst_tot > 0) {
    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n\nCondensation source terms process activated on a total"
         " of %llu cells.\n"), (unsigned long long)n_cells_cst_tot);

    cs_log_separator(CS_LOG_DEFAULT);
  }

  /* 1D-wall module
     -------------- */

  int isuit1;
  if (cs_get_glob_1d_wall_thermal()->use_restart && cs_restart_present())
    isuit1 = 1;
  else if (cs_restart_present())
    isuit1 = cs_restart_present();
  else
    isuit1 = 0;

  cs_1d_wall_thermal_create();
  cs_user_1d_wall_thermal(1);

  cs_get_glob_1d_wall_thermal()->nfpt1t = cs_glob_1d_wall_thermal->nfpt1d;
  if (cs_glob_rank_id > -1)
    cs_parall_counter(&cs_get_glob_1d_wall_thermal()->nfpt1t, 1);

  if (cs_get_glob_1d_wall_thermal()->nfpt1t > 0) {

    cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n\n1D-wall thermal module activated on a total"
         " of %llu boundary faces (%llu local boundary faces)\n"),
       (unsigned long long)cs_get_glob_1d_wall_thermal()->nfpt1t,
       (unsigned long long)cs_glob_1d_wall_thermal->nfpt1d);

    cs_log_separator(CS_LOG_DEFAULT);

  }

  cs_1d_wall_thermal_check(1);

  /* Memory management
     ----------------- */

  cs_f_boundary_conditions_init();

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {
    if (cs_glob_atmo_chemistry->model > 0) {
      cs_f_init_chemistry_reacnum();
    }
  }

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0)
    cs_cf_boundary_conditions_init();

  if (   cs_glob_wall_condensation->icondb == 0
      || cs_glob_wall_condensation->icondv == 0) {
    cs_wall_condensation_create();
  }

  if (cs_get_glob_1d_wall_thermal()->nfpt1t > 0)
    cs_1d_wall_thermal_local_models_create();

  /* Map arrays from Lagrangian module */
  if (cs_glob_lagr_time_scheme->iilagr > 0)
    cs_lagr_init_arrays();

  if (cs_glob_les_balance->i_les_balance > 0)
    cs_les_balance_create();

  /* Default initializations
     ----------------------- */

  const int *bc_type = cs_glob_bc_type;

  cs_field_map_and_init_bcs();

  cs_field_allocate_or_map_all();

  cs_initialize_fields_stage_0();

  if (cs_glob_porous_model >= 1) {
    /* Make fluid surfaces of mesh quantity point to the created fields */
    cs_porous_model_set_has_disable_flag(1);
    cs_porous_model_init_fluid_quantities();
  }

  /* Possible restart
     ---------------- */

  /* Timer statistics */

  int restart_stats_id = cs_timer_stats_id_by_name("checkpoint_restart_stage");
  int post_stats_id = cs_timer_stats_id_by_name("postprocessing_stage");

  if (cs_restart_present() == 1) {

    cs_restart_initialize_fields_read_status();

    cs_timer_stats_start(restart_stats_id);

    cs_restart_map_build();

    /* In the case of points cloud, porosity is a variable
       and readed here */
    cs_restart_main_and_aux_read();

    /* Radiative module restart */
    if (cs_glob_rad_transfer_params->type > 0)
      cs_rad_transfer_read();

    /* Lagrangian module restart (particles) */
    if (cs_glob_lagr_time_scheme->iilagr > 0)
      cs_restart_lagrangian_checkpoint_read();

    cs_les_synthetic_eddy_restart_read();

    cs_porous_model_read();

    /* TODO
       cs_restart_map_free may not be called yet, because
       cs_lagr_solve_initialize and the first call of cs_lagr_solve_time_step
       may also need restart data for particles and statistics respectively.
       This should be solved by moving the corresponding stages at least to
       cs_lagr_solve_initialize sor as to free mapping data before the time loop.
    */

    if (cs_glob_lagr_time_scheme->iilagr < 1)
      cs_restart_map_free();

    cs_timer_stats_stop(restart_stats_id);

  }

  /* Test presence of control_file to modify nt_max if required */

  int ntmsav = ts->nt_max;

  cs_control_check_file();

  if (idtvar == 1 && ntmsav > ts->nt_max && ts->nt_max == ts->nt_cur) {
    if (cs_coupling_is_sync_active())
      ts->nt_max++;
  }

  /* Compute the porosity if needed */
  if (cs_glob_porous_model >= 1) {
    /* Compute porosity from scan */
    if (cs_glob_porosity_from_scan_opt->compute_porosity_from_scan) {

      if (!(cs_glob_porosity_from_scan_opt->use_restart)) {
        cs_log_printf(CS_LOG_DEFAULT,
                      _(" Compute porosity field from scan\n"
                        " WARNING: user porosity will be ignored"
                        " (GUI, cs_user_porosity.c)"));

        cs_compute_porosity_from_scan();
      }
      /* Save pre-process for restart */
      cs_porous_model_write();
      cs_porous_model_fluid_surfaces_preprocessing();
    }
    /* Note using porosity from scan: give the hand to the user */
    else if (cs_glob_porosity_ibm_opt->porosity_mode > 0) {

      cs_log_printf(CS_LOG_DEFAULT,
                    _(" Compute porosity field from immersed boundaries\n"));

      cs_immersed_boundaries(m, mq);
      cs_porous_model_fluid_surfaces_preprocessing();
    }
    else {

      cs_gui_porosity();
      cs_user_porosity(cs_glob_domain);
      cs_porous_model_clip();

    }

    if (cs_glob_porous_model == 3) {
      /* Compute solid quantities and update fluid volume and porosity */
      if (!(cs_glob_porosity_from_scan_opt->use_staircase)) {
        cs_porous_model_mesh_quantities_update();
      }
    }

    cs_mesh_quantities_fluid_vol_reductions(m, mq);
  }

  /* Initialize wall condensation model */
  if (   cs_glob_wall_condensation->icondb == 0
      || cs_glob_wall_condensation->icondv == 0)
    cs_wall_condensation_initialize();

  /* Initialization for the Synthetic turbulence Inlets
     -------------------------------------------------- */

  cs_les_inflow_initialize();

  /* Initializations (user and additional)
     dt rom romb viscl visct viscls (tpucou with periodicity)
     -------------------------------------------------------- */

  /* BC mappings for specific physical models (deprecated) */
  cs_f_pp_models_bc_map();

  if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0
      || cs_glob_physical_model_flag[CS_COMBUSTION_LW] >= 0)
    cs_f_combustion_models_boundary_conditions_map();

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0)
    cs_f_atmo_models_boundary_conditions_map();

  cs_f_initialization_variables();

  if (cs_glob_param_cdo_mode >= CS_PARAM_CDO_MODE_OFF) {  // CDO mode
    assert(cs_glob_domain != nullptr);
    cs_domain_initialize_systems(cs_glob_domain);
  }

  int iterns = -1;
  cs_physical_properties_update(iterns);

  /* Initialization for the atmospheric soil model
     --------------------------------------------- */

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0)
    cs_f_atmsol();

  /* Initialization for the Hybrid Temporal LES model (HTLES)
     -------------------------------------------------------- */

  if (turb_model->hybrid_turb == CS_HYBRID_HTLES)
    cs_htles_initialization();

  /* Initializations for the 1D thermal wall module
     ---------------------------------------------- */

  /* It is assumed that all phases see the same wall temperature.
     Since information is written into a continuation file, a portion
     of memory is needed even after the time loop
     -> IFPT1D and TPPT1D).

     We call cs_user_1d_wall_thermal when there are at least some
     boundary faces with a 1D thermal module. */

  if (cs_get_glob_1d_wall_thermal()->nfpt1t > 0) {

    /* Second call: filling in the geometry definition and
       initialization arrays."(IFPT1D,NPPT1D,EPPT1D,RGPT1D,TPPT1D)
    */
    cs_user_1d_wall_thermal(2);

    cs_1d_wall_thermal_check(2);

    if (isuit1 == 1)
      cs_1d_wall_thermal_read();
    else
      /* Create mesh, initialize temperature */
      cs_1d_wall_thermal_mesh_create();
  }

  /* First pass for the BCs:
     - initilalize bc_type, reference pressure point...
     -------------------------------------------------- */

  /* First pass for initialization BC types
     -- Couplage code_saturne/code_saturne */
  cs_sat_coupling_initialize();

  cs_boundary_conditions_set_coeffs_init();

  /* Arrays for time block, to discard afterwards
     -------------------------------------------- */

  /* Build volume mass injection cell lists when present on at least one rank.
     This is a collective call for consistency, in case the user requires it. */

  if (cs_volume_zone_n_type_zones(CS_VOLUME_ZONE_MASS_SOURCE_TERM) > 0)
    cs_volume_mass_injection_build_lists();

  /* ALE mobile structures */

  if (cs_glob_ale > CS_ALE_NONE)
    cs_mobile_structures_initialize();

  /* Lagrangian initialization */

  if (cs_glob_lagr_time_scheme->iilagr > 0) {

    cs_timer_stats_start(lagr_stats_id);

    cs_lagr_solve_initialize(CS_F_(dt)->val);

    cs_timer_stats_stop(lagr_stats_id);

  }

  /* Solve CDO module(s) or user-defined equations using CDO schemes
     --------------------------------------------------------------- */

  if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_WITH_FV) {
    /* FV and CDO activated */
    cs_cdo_solve_steady_state_domain();
  }

  /* Logging of initial values */

  cs_log_iteration();

  /* Start of time loop
     ------------------ */

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n\n"
       "===============================================================\n\n\n\n"
       "                       MAIN CALCULATION\n"
       "                       ================\n\n\n"
       "===============================================================\n\n\n"));

  ts->nt_cur = ts->nt_prev;
  ts->t_cur = ts->t_prev;

  cs_log_separator(CS_LOG_DEFAULT);

  /* Number of ALE iterations (related to the current execution)
     If CS_GLOB_ALE_NEED_INIT=1, we perform an initialization iteration.
     (if CS_GLOB_ALE_NEED_INIT=-999, it means we have restarted computation
     without rereading "lecamx" -> as if we were performing a restart
     of a computation without ALE)
  */

  if (cs_glob_ale_need_init == -999)
    cs_glob_ale_need_init = 1;

  int itrale = 1;
  if (cs_glob_ale_need_init == 1) {
    itrale = 0;
    cs_log_printf(CS_LOG_DEFAULT,
                  _("\n INSTANT %18.9f    ALE INITIALIZATION\n"
                    "============================================="
                    "================="), ts->t_cur);
  }

  /* In case of code coupling, sync status with other codes. */

  if (itrale > 0) {

    /* Synchronization in dttvar if idtvar = 1
       (i.e. keep coupled codes waiting until time step is computed
       only when needed).
       In case the coupling modifies the reference time step, make sure
       the matching field is updated. Do not do this after initialization.
       except for the adaptive time step (idtvar = 1), handled in dttvar.
    */

    if (idtvar != 1) {
      cs_real_t *dt = CS_F_(dt)->val;

      cs_coupling_sync_apps(0,      /* flags */
                            ts->nt_cur,
                            &(ts->nt_max),
                            &(ts->dt_ref));

      for (cs_lnum_t c_id = 0; c_id < n_cells_ext; c_id++)
        dt[c_id] = ts->dt_ref;
    }

    if (ts->nt_max == ts->nt_cur && ts->nt_max > ts->nt_prev)
      bft_error(__FILE__, __LINE__, 0,
                _("nt_max == nt_cur && nt_max > nt_prev"));

  }

 /* Possible postprocessing of initialization values */

  cs_timer_stats_start(post_stats_id);

  cs_post_activate_by_time_step(ts);

  cs_post_default_write_variables();

  cs_timer_stats_stop(post_stats_id);

  cs_real_t dt_cpl;

  /* Start time loop */

  do {

    if (ts->t_max > 0 && ts->t_max > ts->t_cur) {
      ts->nt_max = ts->nt_cur + (int)((ts->t_max-ts->t_cur)/ts->dt_ref);

      if (ts->nt_max <= ts->nt_cur)
        ts->nt_max = ts->nt_cur + 1;
    }

    if (itrale > 0 && ts->nt_max > ts->nt_prev) {
      cs_timer_stats_increment_time_step();
      /* Time step computed in dttvar if idtvar = 1. */
      if (idtvar != 1)
        cs_time_step_increment(ts->dt_ref);
      else
        cs_time_step_increment(ts->dt[0]);
    }

    /* Step forward in time
       -------------------- */

    /* Test presence of control_file to modify nt_max if required */
    cs_control_check_file();

    if ((idtvar == 0 || idtvar == 1) && (ts->t_max > 0)) {
      if (ts->t_cur >= ts->t_max)
        ts->nt_max = ts->nt_cur;
      else if (ts->nt_max < 0)   /* Changed by control_file */
        ts->nt_max = ts->nt_cur + 1;
    }

    /* Check for runaway (diverging) computation */
    cs_runaway_check();

    /* Set default logging */
    cs_log_iteration_set_active();

    if (idtvar != 1 && ts->nt_max > ts->nt_prev && itrale > 0) {
      if (cs_log_default_is_active())
        cs_log_printf
          (CS_LOG_DEFAULT,
           _("\n INSTANT %18.9f    TIME STEP NUMBER %15d\n"
             " ============================================="
             "=================\n\n\n"),
           ts->t_cur, ts->nt_cur);
    }

    bool mesh_modified = false;
    cs_volume_zone_build_all(mesh_modified);
    cs_boundary_zone_build_all(mesh_modified);

    cs_real_t titer1 = cs_timer_wtime();

    cs_log_iteration_prepare();

    cs_solve_all(itrale);

    cs_1d_wall_thermal_log();

    if (ts->nt_max > ts->nt_prev && itrale > 0) {

      /* Solve CDO module(s) or user-defined equations using CDO schemes
         --------------------------------------------------------------- */

      if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_WITH_FV)
        /* FV and CDO activated */
        cs_cdo_solve_unsteady_state_domain();

      /* Lagrangian module
         ----------------- */

      if (cs_glob_lagr_time_scheme->iilagr > 0) {

        cs_timer_stats_start(lagr_stats_id);

        cs_lagr_solve_time_step(bc_type, CS_F_(dt)->val);

        cs_timer_stats_stop(lagr_stats_id);

      }

      /* Update gradients needed in LES balance computation
         -------------------------------------------------- */

      if (cs_glob_les_balance->i_les_balance > 0)
        cs_les_balance_update_gradients();

      /* Compute temporal means (accumulation)
         ------------------------------------- */

      cs_time_moment_update_all();

    }

    /* Update mesh (ALE)
       ----------------- */

    if (cs_glob_ale > CS_ALE_NONE && ts->nt_max > ts->nt_prev) {

      if (itrale == 0 || itrale > cs_glob_ale_n_ini_f)
        cs_ale_update_mesh(itrale);

    }

    /* Optional processing by user
       --------------------------- */

    if (itrale > 0) {

      cs_timer_stats_start(post_stats_id);

      /* 1D profiles postprocessing output */

      cs_gui_balance_by_zone();
      cs_gui_pressure_drop_by_zone();

      cs_f_user_extra_operations_wrapper(CS_F_(dt)->val);

      cs_user_extra_operations(cs_glob_domain);

      if (cs_glob_les_balance->i_les_balance > 0)
        cs_les_balance_compute();

      cs_timer_stats_stop(post_stats_id);

    }

    /* Stop tests
       ---------- */

    /* Test for lack of remaining time */

    cs_resource_get_max_timestep(ts->nt_cur, &ts->nt_max);

    /* Stop test for couplings */

    if (idtvar != 1) {  /* synchronization in dttvar if idtvar = 1 */
      dt_cpl = ts->dt_ref;

      cs_coupling_sync_apps(0,      /* flags */
                            ts->nt_cur,
                            &(ts->nt_max),
                            &dt_cpl);

    }

    /* Possible output of checkpoint files
       ----------------------------------- */

    bool restart_checkpoint_required = cs_restart_checkpoint_required(ts);

    if (ts->nt_cur < ts->nt_max && itrale == 0)
      restart_checkpoint_required = false;

    if (restart_checkpoint_required) {

      cs_timer_stats_start(restart_stats_id);

      const char *info = (ts->nt_cur == ts->nt_max) ? "final" : "intermediate";

      cs_log_printf
        (CS_LOG_DEFAULT,
         _("\n\n Write %s restart files\n"
           "   Checkpoint at iteration %d, physical time %15.5f\n\n"),
         info, ts->nt_cur, ts->t_cur);

      bool checkpoint_mesh = false;
      if (   iturbo == CS_TURBOMACHINERY_TRANSIENT
          && cs_glob_restart_auxiliary->write_auxiliary == 1)
        checkpoint_mesh = true;

      cs_time_stepping_write_checkpoint(checkpoint_mesh);

      /* Indicate that a chechpoint has been done */
      cs_restart_checkpoint_done(ts);

      /* Remove all unnecessary previous dumps of checkpoint files */
      cs_restart_clean_multiwriters_history();

      cs_timer_stats_stop(restart_stats_id);

    } // End test on restart_checkpoint_required

    /* Test to determine if a visualization output is generated
       -------------------------------------------------------- */

    cs_timer_stats_start(post_stats_id);

    cs_post_activate_by_time_step(ts);

    /* When geometry has not been output yet, deactivate all writers. */
    if (itrale == 0)
      cs_post_activate_writer(0, false);

    /* Standard visualization output
       ----------------------------- */

    cs_post_default_write_variables();

    /* CDO module (user-defined equations)
       ----------------------------------- */
    if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_WITH_FV)
      /* FV and CDO activated */
      cs_domain_post(cs_glob_domain);

    /* Write to "run_solver.log" periodically
       -------------------------------------- */

    if (cs_log_default_is_active()) {

      cs_log_equation_convergence_info_write();

      cs_log_iteration();

      cs_log_iteration_l2residual();

    }

    cs_timer_stats_stop(post_stats_id);

    cs_real_t titer2 = cs_timer_wtime();

    if (itrale <= 0) {

      cs_log_printf(CS_LOG_DEFAULT,
                    _("\n Time for ALE initialization:        %14.5f s.\n"),
                    titer2 - titer1);

      cs_log_separator(CS_LOG_DEFAULT);
    }

    /* End of time loop
       ---------------- */

    itrale = itrale + 1;

  } while (ts->nt_cur < ts->nt_max);

  /* Final synchronization for time step.
     This is done after exiting the main time loop, hence telling other codes
     that code_saturne is finished. */

  cs_coupling_sync_apps(0,      /* flags */
                        ts->nt_cur,
                        &(ts->nt_max),
                        &dt_cpl);

  /* Free intermediate arrays */

  if (cs_restart_present() && cs_glob_lagr_time_scheme->iilagr > 0) {
    cs_timer_stats_start(restart_stats_id);
    cs_restart_map_free();
    cs_timer_stats_stop(restart_stats_id);
  }

  /* Finalize probes
     --------------- */

  cs_log_printf
      (CS_LOG_DEFAULT,
       _("\n\n"
         "=============================================================\n\n\n\n"
         "                 FINAL STAGE OF THE CALCULATION\n"
         "                 ==============================\n\n\n"
         " ===========================================================\n\n\n"));

  if (cs_glob_1d_wall_thermal->nfpt1d > 0)
    cs_1d_wall_thermal_free();

  /* Free main arrays */

  cs_restart_finalize_fields_read_status();

  cs_rad_transfer_finalize();

  cs_turbulence_bc_free_pointers();
  cs_boundary_conditions_free();

  cs_f_finalize_meteo();

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 0) {

    if (cs_glob_atmo_imbrication->imbrication_flag)
      cs_f_finalize_imbrication();

    cs_at_data_assim_finalize();

    if (cs_glob_atmo_chemistry->model > 0)
      cs_f_finalize_chemistry();

  }

  if (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0)
    cs_f_finalize_steady_laminar_flamelet_library();

  if (cs_glob_physical_model_flag[CS_GAS_MIX] >= 0)
    cs_gas_mix_finalize();

  if (cs_glob_ale >= 1)
    cs_mobile_structures_finalize();

  if (   cs_glob_1d_wall_thermal->nfpt1d > 0
      || cs_get_glob_1d_wall_thermal()->nfpt1t == 0)
    cs_1d_wall_thermal_finalize();

  if (cs_glob_les_balance->i_les_balance > 0)
    cs_les_balance_finalize();

  cs_log_printf
    (CS_LOG_DEFAULT,
     _("\n\n"
       " ===========================================================\n\n\n\n"
       "                      END OF CALCULATION\n"
       "                      ==================\n\n\n\n"
       "=============================================================\n"));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Output a checkpoint.
 *
 * If needed, the mesh is also output in the checkpoint directory,
 * exect if this function is called for checkpoint serialized in memory
 * (which is a special case for FMI exchange).
 *
 * \param[in]  checkpoint_mesh  also save mesh in checkpoint directory
 */
/*----------------------------------------------------------------------------*/

void
cs_time_stepping_write_checkpoint(bool  checkpoint_mesh)
{
  cs_restart_main_and_aux_write();

  if (checkpoint_mesh)
    cs_mesh_save(cs_glob_mesh, nullptr, "checkpoint", "mesh.csm");

  if (cs_get_glob_1d_wall_thermal()->nfpt1t > 0)
    cs_1d_wall_thermal_write();

  cs_les_synthetic_eddy_restart_write();

  if (cs_glob_lagr_time_scheme->iilagr > 0)
    cs_restart_lagrangian_checkpoint_write();

  if (cs_glob_rad_transfer_params->type > 0)
    cs_rad_transfer_write();

  if (cs_glob_les_balance->i_les_balance > 0)
    cs_les_balance_write_restart();
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
