/*============================================================================
 * Main program
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_atmo.h"
#include "cs_all_to_all.h"
#include "cs_ast_coupling.h"
#include "cs_base.h"

#include "cs_base_fortran.h"
#include "cs_benchmark.h"
#include "cs_boundary_zone.h"
#include "cs_boundary_conditions.h"
#include "cs_calcium.h"
#include "cs_cdo_main.h"
#include "cs_cell_to_vertex.h"
#include "cs_control.h"
#include "cs_coupling.h"
#include "cs_ctwr.h"
#include "cs_domain_setup.h"
#include "cs_ext_library_info.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_fp_exception.h"
#include "cs_function.h"
#include "cs_gradient.h"
#include "cs_gui.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_conjugate_heat_transfer.h"
#include "cs_gui_mesh.h"
#include "cs_gui_mobile_mesh.h"
#include "cs_gui_output.h"
#include "cs_gui_particles.h"
#include "cs_gui_radiative_transfer.h"
#include "cs_gui_util.h"
#include "cs_io.h"
#include "cs_ibm.h"
#include "cs_join.h"
#include "cs_les_inflow.h"
#include "cs_log.h"
#include "cs_log_setup.h"
#include "cs_log_iteration.h"
#include "cs_matrix_default.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_coherency.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_mesh_smoother.h"
#include "cs_notebook.h"
#include "cs_opts.h"
#include "cs_parall.h"
#include "cs_param_cdo.h"
#include "cs_paramedmem_coupling.h"
#include "cs_parameters.h"
#include "cs_physical_properties.h"
#include "cs_post.h"
#include "cs_post_default.h"
#include "cs_preprocess.h"
#include "cs_preprocessor_data.h"
#include "cs_probe.h"
#include "cs_property.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_restart.h"
#include "cs_restart_map.h"
#include "cs_runaway_check.h"
#include "cs_sles.h"
#include "cs_sles_default.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"
#include "cs_sys_coupling.h"
#include "cs_system_info.h"
#include "cs_time_moment.h"
#include "cs_time_table.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_tree.h"
#include "cs_turbomachinery.h"
#include "cs_utilities.h"
#include "cs_vertex_to_cell.h"
#include "cs_volume_mass_injection.h"
#include "cs_volume_zone.h"

#if defined(HAVE_CUDA)
#include "cs_base_cuda.h"
#include "cs_blas_cuda.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_opts_t  opts;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handle first SIGINT by setting the max number of iterations to
 * the current value.
 *
 * \param[in]  signum  signal number.
 */
/*----------------------------------------------------------------------------*/

static void
_sigint_handler(int signum)
{
  cs_time_step_define_nt_max(cs_glob_time_step->nt_cur);

  cs_log_printf(CS_LOG_DEFAULT,
                _("Signal %d received.\n"
                  "--> computation interrupted by environment.\n"
                  "\n"
                  "    maximum time step number set to: %d\n"),
                signum,
                cs_glob_time_step->nt_max);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Main function for code_saturne run.
 *
 * This function is called by main().
 *----------------------------------------------------------------------------*/

static void
_run(void)
{
  int  check_mask = 0;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  cs_base_sigint_handler_set(_sigint_handler);

  /* System information */

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif
  cs_ext_library_info();

#if defined(HAVE_CUDA)
  cs_base_cuda_select_default_device();
#endif

  if (cs_get_device_id() < 0)
    cs_halo_set_buffer_alloc_mode(CS_ALLOC_HOST);

  cs_timer_stats_initialize();
  cs_timer_stats_define_defaults();

  if (cs_glob_tree == NULL)
    cs_glob_tree = cs_tree_node_create(NULL);

  cs_gui_parallel_io();
  if (cs_glob_n_ranks > 1)
    cs_user_parallel_io();
  cs_file_defaults_info();

  cs_gui_mpi_algorithms();

  /* Load user time tables */
  cs_gui_time_tables();
  cs_user_time_table();

  bft_printf("\n");

  cs_base_update_status("initializing\n");

  /* Initialize random number generator
     (used in only some cases, but safe to do, and inexpensive) */

  cs_random_seed(cs_glob_rank_id + 1);

  /* Initialize global structures for main mesh */

  cs_mesh_location_initialize();
  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();
  cs_boundary_zone_initialize();
  cs_volume_zone_initialize();

  cs_preprocess_mesh_define();

  cs_turbomachinery_define();

  /* Call main calculation initialization function or help */

  cs_io_log_initialize();

  cs_field_define_keys_base();
  cs_parameters_define_field_keys();

  cs_sles_initialize();
  cs_sles_set_default_verbosity(cs_sles_default_get_verbosity);

  cs_preprocessor_data_read_headers(cs_glob_mesh,
                                    cs_glob_mesh_builder,
                                    false);

  cs_gui_zones();
  cs_user_zones();

  cs_ctwr_define_zones();

  /* Create a new structure for the computational domain */

  cs_glob_domain = cs_domain_create();

  /* Define MPI-based Couplings if applicable */

  cs_gui_syrthes_coupling();
  cs_user_syrthes_coupling();
  cs_user_saturne_coupling();

  /* Initialize Fortran API and calculation setup */

  if ((opts.preprocess | opts.verif) == false && opts.benchmark <= 0) {

    int _rank_id = cs_glob_rank_id, _n_ranks = cs_glob_n_ranks;

    cs_base_fortran_bft_printf_to_f();

    const char default_restart_mesh[] = "restart_mesh_input";
    if (cs_file_isreg(default_restart_mesh))
      cs_restart_map_set_mesh_input(default_restart_mesh);

    CS_PROCF(csinit, CSINIT)(&_rank_id, &_n_ranks);

    CS_PROCF(initi1, INITI1)();

    if (cs_parameters_need_extended_neighborhood())
      halo_type = CS_HALO_EXTENDED;

    if (cs_glob_ale > CS_ALE_NONE) {
      cs_restart_map_set_locations(true, true);
      cs_gui_mobile_mesh_get_boundaries(cs_glob_domain);
      if (cs_glob_mesh->time_dep < CS_MESH_TRANSIENT_COORDS)
        cs_glob_mesh->time_dep = CS_MESH_TRANSIENT_COORDS;
    }

    /* Initialization step for the setup of CDO schemes
     * - Perform the initialization for each activated module
     * - Create predefined properties and advection fields
     * - Create fields associated to equations or modules
     * - Create fields associated to user-defined equations
     */

    cs_cdo_initialize_setup(cs_glob_domain);

    /* Setup linear solvers */

    cs_gui_linear_solvers();
    cs_user_linear_solvers();

    cs_base_fortran_bft_printf_to_c();

    cs_timer_stats_set_start_time(cs_glob_time_step->nt_cur);

  }
  else if (opts.verif)
    halo_type = CS_HALO_EXTENDED;

  /* Discover applications visible through MPI (requires communication);
     this is done after main calculation initialization so that the user
     may have the option of assigning a name to this instance. */

#if defined(HAVE_MPI)
  cs_coupling_discover_mpi_apps(opts.app_name, NULL);
#endif

  if (opts.app_name != NULL)
    BFT_FREE(opts.app_name);

  /* Initialize couplings and communication if necessary */

  cs_syr_coupling_all_init();
  cs_sat_coupling_all_init();
  cs_sys_coupling_all_init();

  cs_paramedmem_coupling_all_init();

  /* Initialize main post-processing */

  cs_gui_postprocess_writers();
  cs_user_postprocess_writers();
  cs_post_init_writers();

  cs_gui_postprocess_meshes();
  cs_user_postprocess_meshes();
  cs_user_postprocess_probes();

  /* Print info on fields and associated keys and other setup options */

  if (opts.verif == false && opts.preprocess == false && opts.benchmark <= 0)
    cs_log_setup();

  /* Preprocess mesh */

  cs_preprocess_mesh(halo_type);
  cs_mesh_adjacencies_initialize();

  /* Initialization for turbomachinery computations */

  cs_turbomachinery_initialize();

  /* Initialization of internal coupling */

  cs_internal_coupling_initialize();

  cs_internal_coupling_dump();

  /* Initialize meshes for the main post-processing */

  check_mask = ((opts.preprocess | opts.verif) == true) ? 2 + 1 : 0;

  cs_post_init_meshes(check_mask);

  /* Compute iterations or quality criteria depending on verification options */

  if (opts.verif == true) {
    bft_printf(_("\n Computing quality criteria\n"));
    cs_mesh_quality(cs_glob_mesh, cs_glob_mesh_quantities);
    cs_mesh_coherency_check();
    cs_mesh_bad_cells_postprocess(cs_glob_mesh, cs_glob_mesh_quantities);
  }
  else if (opts.preprocess == true)
    cs_mesh_coherency_check();

  if (check_mask && cs_syr_coupling_n_couplings())
    bft_error(__FILE__, __LINE__, 0,
              _("Coupling with SYRTHES is not possible in mesh preprocessing\n"
                "or verification mode."));

  if (opts.preprocess == false) {

    cs_mesh_adjacencies_update_mesh();

#if defined(HAVE_ACCEL)
    cs_preprocess_mesh_update_device(cs_alloc_mode);
#endif

  }

  if (opts.benchmark > 0) {
    int mpi_trace_mode = (opts.benchmark == 2) ? 1 : 0;
    cs_benchmark(mpi_trace_mode);
  }

  if (opts.preprocess == false && opts.benchmark <= 0) {

    /* Check that mesh seems valid */

    cs_mesh_quantities_check_vol(cs_glob_mesh,
                                 cs_glob_mesh_quantities,
                                 (opts.verif ? 1 : 0));

    /* Initialization related to CDO/HHO schemes */

    cs_cdo_initialize_structures(cs_glob_domain,
                                 cs_glob_mesh,
                                 cs_glob_mesh_quantities);

    /* Initialize gradient computation */

    cs_gradient_initialize();

    if (opts.verif == false) {

      /* Initialize sparse linear systems resolution */

      cs_user_matrix_tuning();

      cs_matrix_initialize();

      /* Update Fortran mesh sizes and quantities */

      cs_base_fortran_bft_printf_to_f();
      cs_preprocess_mesh_update_fortran();

      /* Choose between standard and user solver */

      if (cs_user_solver_set() == 0) {

        if (cs_glob_param_cdo_mode == CS_PARAM_CDO_MODE_ONLY) {

          /* Only C language is called within CDO */

          cs_base_fortran_bft_printf_to_c();

          /*----------------------------------------------
           * Call main calculation function (CDO Kernel)
           *----------------------------------------------*/

          cs_cdo_main(cs_glob_domain);

          /* Return to the default behavior */

          cs_base_fortran_bft_printf_to_f();

        }
        else {

          /* Additional initializations required by some models */

          cs_fan_build_all(cs_glob_mesh, cs_glob_mesh_quantities);

          cs_ctwr_build_all();

          cs_volume_mass_injection_flag_zones();

          /* Setup couplings and fixed-mesh postprocessing */

          cs_syr_coupling_init_meshes();

          cs_paramedmem_coupling_define_mesh_fields();

          cs_post_default_write_meshes();

          cs_turbomachinery_restart_mesh();

          /*----------------------------------------------
           * Call main calculation function (code Kernel)
           *----------------------------------------------*/

          /* Maybe should be allocate in caltri.c */
          cs_ale_allocate();

          CS_PROCF(caltri, CALTRI)();

        }

      }
      else {

          /*--------------------------------
           * Call user calculation function
           *--------------------------------*/

          cs_user_solver(cs_glob_mesh,
                         cs_glob_mesh_quantities);

      }

      /* Finalize user extra operations */

      cs_user_extra_operations_finalize(cs_glob_domain);

    }

    /* Finalize gradient computation */

    cs_gradient_finalize();

    /* Finalize synthetic inlet condition generation */

    cs_les_inflow_finalize();

  }

  /* Finalize linear system resolution */

  cs_sles_default_finalize();

  /* Finalize matrix API */

  cs_matrix_finalize();

#if defined(HAVE_CUDA)
  cs_blas_cuda_finalize();
#endif

  /* Switch logging back to C (may be moved depending on Fortran dependencies) */

  cs_base_fortran_bft_printf_to_c();

  bft_printf(_("\n Destroying structures and ending computation\n"));
  bft_printf_flush();

  /* Free the boundary conditions face type and face zone arrays */

  cs_boundary_conditions_free();

  /* Final stage for CDO/HHO schemes */

  cs_cdo_finalize(cs_glob_domain);

  /* Free cs_domain_structure */

  cs_domain_free(&cs_glob_domain);

  /* Free coupling-related data */

#if defined(HAVE_MPI)
  cs_ast_coupling_finalize();
  cs_syr_coupling_all_finalize();
  cs_sat_coupling_all_finalize();
  cs_sys_coupling_all_finalize();
  cs_paramedmem_coupling_all_finalize();
  cs_coupling_finalize();
#endif

  cs_control_finalize();

  /* Free remapping/intersector related structures (stl or medcoupling) */

  cs_utilities_destroy_all_remapping();

  /* Free the checkpoint multiwriter structure */

  cs_restart_multiwriters_destroy_all();

  /* Print some mesh statistics */

  cs_gui_usage_log();
  cs_mesh_selector_stats(cs_glob_mesh);

  /* Finalizations related to some models */

  cs_atmo_finalize();
  cs_ctwr_all_destroy();
  cs_fan_destroy_all();

  /* Free internal coupling */

  cs_internal_coupling_finalize();

  /* Free memory related to properties */

  cs_property_destroy_all();
  cs_thermal_table_finalize();

  /* Free immersed boundaries related structures */

  cs_ibm_finalize();

  /* Free turbomachinery related structures */

  cs_turbomachinery_finalize();
  cs_join_finalize();

  /* Free post processing or logging related structures */

  cs_probe_finalize();
  cs_post_finalize();
  cs_log_iteration_destroy_all();

  cs_function_destroy_all();

  /* Free moments info */

  cs_time_moment_destroy_all();

  /* Free field info */

  cs_gui_radiative_transfers_finalize();
  cs_gui_finalize();

  cs_notebook_destroy_all();
  cs_time_table_destroy_all();

  cs_field_pointer_destroy_all();
  cs_field_destroy_all();
  cs_field_destroy_all_keys();

  /* Free Physical model related structures
     TODO: extend cs_base_atexit_set mechanism to allow registering
     model-specific cleanup functions at their activation point,
     avoiding the need for modification of this function, which
     should be more "generic". */

  cs_base_finalize_sequence();

  /* Free main mesh after printing some statistics */

  cs_cell_to_vertex_free();
  cs_vertex_to_cell_free();
  cs_mesh_adjacencies_finalize();

  cs_boundary_zone_finalize();
  cs_volume_zone_finalize();
  cs_mesh_location_finalize();
  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* Free parameters tree info */

  cs_tree_node_free(&cs_glob_tree);

  /* CPU times and memory management finalization */

  cs_all_to_all_log_finalize();
  cs_io_log_finalize();

  cs_timer_stats_finalize();

  cs_file_free_defaults();

  cs_base_time_summary();

#if defined(HAVE_ACCEL)
  int n_alloc_hd_remain = cs_get_n_allocations_hd();
  if (n_alloc_hd_remain > 0)
    bft_printf(_("Warning: %d remaining host-device allocations\n"
                 "         (possible memory leak)\n"), n_alloc_hd_remain);
#endif

  cs_base_mem_finalize();

  cs_log_printf_flush(CS_LOG_N_TYPES);

  cs_runaway_check_finalize();
}

/*============================================================================
 * Main program
 *============================================================================*/

int
main(int    argc,
     char  *argv[])
{
  /* Initialize wall clock timer */

  (void)cs_timer_wtime();

  /* First analysis of the command line to determine if MPI is required,
     and MPI initialization if it is. */

#if defined(HAVE_MPI)
  cs_base_mpi_init(&argc, &argv);
#endif

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  {
    int t_id;
#pragma omp parallel private(t_id)
    {
      t_id = omp_get_thread_num();
      if (t_id == 0)
        cs_glob_n_threads = omp_get_max_threads();
    }
  }

  if (cs_glob_n_threads > 1)
    cs_glob_e2n_sum_type = CS_E2N_SUM_SCATTER;
#endif

  /* Default initialization */

#if defined(_CS_ARCH_Linux)

  if (getenv("LANG") != NULL)
    setlocale(LC_ALL,"");
  else
    setlocale(LC_ALL, "C");
  setlocale(LC_NUMERIC, "C");

#endif

  /* Trap floating-point exceptions on most systems */

#if defined(DEBUG)
  cs_fp_exception_enable_trap();
#endif

  /* Initialize memory management */

  cs_base_mem_init();

  /* Initialize internationalization */

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, cs_base_get_localedir());
  textdomain(PACKAGE);
#endif

  /* Parse command line */

  cs_opts_define(argc, argv, &opts);

  /* Initialize error handling */

  cs_base_error_init(opts.sig_defaults);

  /* Open 'run_solver.log' (log) files */

  cs_base_trace_set(opts.trace);
  cs_base_fortran_bft_printf_set("run_solver", opts.logrp);

  /* Log-file header and command line arguments recap */

  cs_base_logfile_head(argc, argv);

  /* Load setup parameters if present */

  const char s_param[] = "setup.xml";
  if (cs_file_isreg(s_param)) {
    cs_gui_load_file(s_param);
    cs_notebook_load_from_file();
  }

  /* Call main run() method */

  _run();

  /* Return */

  cs_exit(EXIT_SUCCESS);

  /* Never called, but avoids compiler warning */
  return 0;
}

/*----------------------------------------------------------------------------*/
