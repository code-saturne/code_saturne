/*============================================================================
 * Main program
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

#include "fvm_selector.h"

#include "cs_all_to_all.h"
#include "cs_base.h"
#include "cs_base_fortran.h"
#include "cs_benchmark.h"
#include "cs_boundary_zone.h"
#include "cs_calcium.h"
#include "cs_cdo_main.h"
#include "cs_control.h"
#include "cs_coupling.h"
#include "cs_ctwr.h"
#include "cs_fan.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_file.h"
#include "cs_fp_exception.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_gui.h"
#include "cs_gui_output.h"
#include "cs_gui_particles.h"
#include "cs_gui_radiative_transfer.h"
#include "cs_io.h"
#include "cs_join.h"
#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
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
#include "cs_opts.h"
#include "cs_parameters.h"
#include "cs_partition.h"
#include "cs_physical_properties.h"
#include "cs_post.h"
#include "cs_post_default.h"
#include "cs_preprocess.h"
#include "cs_preprocessor_data.h"
#include "cs_probe.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_sles_default.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"
#include "cs_system_info.h"
#include "cs_time_moment.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_les_inflow.h"
#include "cs_turbomachinery.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Main function for Code_Saturne run.
 *
 * This function is called either by main() in the standard case, or by
 * a SALOME standalone module's yacsinit() function. As the latter cannot
 * return, almost all of the code's execution steps (except command-line
 * initialization, required to know if we are running in standard mode or
 * pluging a SALOME module) are done here.
 *
 * As yacsinit() can take no arguments, the command-line options must be
 * defined as a static global variable by main() so as to be usable
 * by run().
 *----------------------------------------------------------------------------*/

void
cs_run(void);

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_opts_t  opts;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Main function for Code_Saturne run.
 *
 * This function is called either by main() in the standard case, or by
 * a SALOME standalone module's yacsinit() function. As the latter cannot
 * return, almost all of the code's execution steps (except command-line
 * initialization, required to know if we are running in standard mode or
 * pluging a SALOME module) are done here.
 *
 * As yacsinit() can take no arguments, the command-line options must be
 * defined as a static global variable by main() so as to be usable
 * by run().
 *----------------------------------------------------------------------------*/

void
cs_run(void)
{
  cs_int_t  ivoset = 0;

  int  check_mask = 0;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  /* System information */

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif

  cs_timer_stats_initialize();
  cs_timer_stats_define_defaults();

  cs_gui_parallel_io();
  cs_user_parallel_io();
  cs_file_defaults_info();

  cs_partition_external_library_info();

  bft_printf("\n");

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

  /* Call main calculation initialization function or help */

  cs_io_log_initialize();

  cs_field_define_keys_base();
  cs_parameters_define_field_keys();

  cs_sles_initialize();
  cs_sles_set_default_verbosity(cs_sles_default_get_verbosity);

  cs_preprocessor_data_read_headers(cs_glob_mesh,
                                    cs_glob_mesh_builder);
  if (!opts.cdo)
    opts.cdo = cs_user_cdo_activated();

  /* Initialize Fortran API and calculation setup */

  if ((opts.preprocess | opts.verif) == false && opts.benchmark <= 0) {

    cs_int_t _rank_id = cs_glob_rank_id, _n_ranks = cs_glob_n_ranks;

    cs_base_fortran_bft_printf_to_f();

    cs_gui_init();

    cs_gui_zones();
    cs_user_zones();

    if (opts.cdo)
      cs_cdo_initialize_setup();

    else {

      CS_PROCF(csinit, CSINIT)(&_rank_id, &_n_ranks);

      CS_PROCF(initi1, INITI1)();

      CS_PROCF (haltyp, HALTYP) (&ivoset);
      if (ivoset)
        halo_type = CS_HALO_EXTENDED;

    }

    cs_base_fortran_bft_printf_to_c();

    cs_ctwr_build_zones();

    cs_timer_stats_set_start_time(cs_glob_time_step->nt_cur);

  }
  else if (opts.verif)
    halo_type = CS_HALO_EXTENDED;

  /* Discover applications visible through MPI (requires communication);
     this is done after main calculation initialization so that the user
     may have the option of assigning a name to this instance. */

#if defined(HAVE_MPI)
  cs_coupling_discover_mpi_apps(opts.app_name);
#endif

  if (opts.app_name != NULL)
    BFT_FREE(opts.app_name);

  /* Initialize couplings and communication if necessary */

  cs_syr_coupling_all_init();
  cs_sat_coupling_all_init();

  /* Initialize main post-processing */

  cs_gui_postprocess_writers();
  cs_user_postprocess_writers();

  cs_gui_postprocess_meshes();
  cs_user_postprocess_meshes();
  cs_user_postprocess_probes();

  cs_post_init_writers();

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

  if (opts.benchmark > 0) {
    int mpi_trace_mode = (opts.benchmark == 2) ? 1 : 0;
    cs_benchmark(mpi_trace_mode);
  }

  if (check_mask && cs_syr_coupling_n_couplings())
    bft_error(__FILE__, __LINE__, errno,
              _("Coupling with SYRTHES is not possible in mesh preprocessing\n"
                "or verification mode."));

  if (opts.preprocess == false && opts.benchmark <= 0) {

    /* Check that mesh seems valid */

    cs_mesh_quantities_check_vol(cs_glob_mesh,
                                 cs_glob_mesh_quantities,
                                 (opts.verif ? 1 : 0));

    cs_mesh_adjacencies_update_mesh();

    /* Initialization related to CDO/HHO schemes */

    cs_cdo_initialize_structures(cs_glob_mesh, cs_glob_mesh_quantities);

    /* Initialize gradient computation */

    cs_gradient_initialize();
    cs_gradient_perio_initialize();

    if (opts.verif == false) {

      /* Initialize sparse linear systems resolution */

      cs_user_matrix_tuning();

      cs_matrix_initialize();

      /* Update Fortran mesh sizes and quantities */

      cs_base_fortran_bft_printf_to_f();
      cs_preprocess_mesh_update_fortran();

      /* Choose between standard and user solver */

      if (cs_user_solver_set() == 0) {

        if (opts.cdo) {

          /* Only C language is called within CDO */

          cs_base_fortran_bft_printf_to_c();

          /*----------------------------------------------
           * Call main calculation function (CDO Kernel)
           *----------------------------------------------*/

          cs_cdo_main();

          /* Return to the default behavior */

          cs_base_fortran_bft_printf_to_f();

        }
        else {

          /* Additional initializations required by some models */

          cs_fan_build_all(cs_glob_mesh, cs_glob_mesh_quantities);

          cs_ctwr_build_all();

          /* Setup couplings and fixed-mesh postprocessing */

          cs_syr_coupling_init_meshes();

          cs_post_default_write_meshes();

          cs_turbomachinery_restart_mesh();

          /*----------------------------------------------
           * Call main calculation function (code Kernel)
           *----------------------------------------------*/

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

      /* Finalize sparse linear systems resolution */

      cs_sles_default_finalize();

      cs_matrix_finalize();

    }

    /* Finalize gradient computation */

    cs_gradient_perio_finalize();
    cs_gradient_finalize();

    /* Finalize synthetic inlet condition generation */

    cs_inflow_finalize();

  }

  /* Switch logging back to C (may be moved depending on Fortran dependencies) */

  cs_base_fortran_bft_printf_to_c();

  bft_printf(_("\n Destroying structures and ending computation\n"));
  bft_printf_flush();

  /* Final stage for CDO/HHO schemes */

  cs_cdo_finalize();

  /* Free coupling-related data */

  cs_syr_coupling_all_finalize();
#if defined(HAVE_MPI)
  cs_sat_coupling_all_finalize();
  cs_coupling_finalize();
#endif

  cs_control_finalize();

  /* Print some mesh statistics */

  cs_gui_usage_log();
  cs_mesh_selector_stats(cs_glob_mesh);

  /* Finalizations related to some models */

  cs_ctwr_all_destroy();
  cs_fan_destroy_all();

  /* Free internal coupling */

  cs_internal_coupling_finalize();

  /* Free thermal physical properties */

  cs_thermal_table_finalize();

  /* Free turbomachinery related structures */

  cs_turbomachinery_finalize();
  cs_join_finalize();

  /* Free post processing or logging related structures */

  cs_probe_finalize();
  cs_post_finalize();
  cs_log_iteration_destroy_all();

  /* Free moments info */

  cs_time_moment_destroy_all();

  /* Free field info */

  cs_gui_radiative_transfers_finalize();
  cs_gui_finalize();

  cs_field_pointer_destroy_all();
  cs_field_destroy_all();
  cs_field_destroy_all_keys();

  /* Free Lagrangian related structures */

  cs_lagr_finalize();

  /* Free main mesh after printing some statistics */

  cs_mesh_adjacencies_finalize();

  cs_boundary_zone_finalize();
  cs_volume_zone_finalize();
  cs_mesh_location_finalize();
  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* CPU times and memory management finalization */

  cs_all_to_all_log_finalize();
  cs_io_log_finalize();

  cs_timer_stats_finalize();

  cs_file_free_defaults();

  cs_base_time_summary();
  cs_base_mem_finalize();
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

  /* Open 'listing' (log) files */

  cs_base_fortran_bft_printf_set("listing", opts.ilisr0, opts.ilisrp);

  /* Log-file header and command line arguments recap */

  cs_base_logfile_head(argc, argv);

  /* Running as a standalone SALOME component, load YACS component
     library and run yacsinit() component initialization and event loop,
     which should itself include the standard run routine */

  if (opts.yacs_module != NULL) {
    cs_calcium_load_yacs(opts.yacs_module);
    BFT_FREE(opts.yacs_module);
    cs_calcium_start_yacs(); /* Event-loop does not return as of this version */
    cs_calcium_unload_yacs();
  }

  /* In standard case, simply call regular run() method */

  else
    cs_run();

  /* Return */

  cs_exit(EXIT_SUCCESS);

  /* Never called, but avoids compiler warning */
  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
