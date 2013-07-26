/*============================================================================
 * Main program
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

/* On glibc-based systems, define _GNU_SOURCE so as to enable floating-point
   error exceptions; on Itanium, optimized code may raise such exceptions
   due to speculative execution, so we only enable raising of such exceptions
   for code compiled in debug mode, where reduced optimization should not lead
   to such exceptions, and locating the "true" origin of floating-point
   exceptions is helpful.
   _GNU_SOURCE must be defined before including any headers, to ensure
   the correct feature macros are defined first. */

#if defined(__linux__) || defined(__linux) || defined(linux)
#if    (!defined(__ia64__) && !defined(__blrts__) && !defined(__bg__)) \
    || defined(DEBUG)
#define CS_FPE_TRAP
#define _GNU_SOURCE
#endif
#endif

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(CS_FPE_TRAP)
#include <fenv.h>
#endif

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
#include "cs_calcium.h"
#include "cs_coupling.h"
#include "cs_ctwr.h"
#include "cs_field.h"
#include "cs_file.h"
#include "cs_gradient.h"
#include "cs_gradient_perio.h"
#include "cs_gui.h"
#include "cs_gui_mesh.h"
#include "cs_gui_output.h"
#include "cs_gui_particles.h"
#include "cs_io.h"
#include "cs_join.h"
#include "cs_lagr_tracking.h"
#include "cs_log.h"
#include "cs_log_iteration.h"
#include "cs_mesh.h"
#include "cs_mesh_coherency.h"
#include "cs_mesh_from_builder.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_bad_cells.h"
#include "cs_mesh_save.h"
#include "cs_mesh_smoother.h"
#include "cs_mesh_to_builder.h"
#include "cs_mesh_warping.h"
#include "cs_multigrid.h"
#include "cs_opts.h"
#include "cs_parameters.h"
#include "cs_partition.h"
#include "cs_post.h"
#include "cs_preprocessor_data.h"
#include "cs_prototypes.h"
#include "cs_renumber.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"
#include "cs_system_info.h"
#include "cs_timer.h"
#include "cs_les_inflow.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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

#if defined(CS_FPE_TRAP)
static int _fenv_set = 0;    /* Indicates if behavior modified */
static fenv_t _fenv_old;     /* Old exception mask */
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
  cs_int_t  ivoset;
  double  t1, t2;

  bool partition_preprocess = false;
  int  check_mask = 0;
  int  cwf_post = 0;
  double  cwf_threshold = -1.0;
  cs_halo_type_t halo_type = CS_HALO_STANDARD;

  /* System information */

#if defined(HAVE_MPI)
  cs_system_info(cs_glob_mpi_comm);
#else
  cs_system_info();
#endif

  cs_gui_parallel_io();
  cs_user_parallel_io();
  cs_file_defaults_info();

  cs_partition_external_library_info();

  bft_printf("\n");

  /* Initialize global structures for main mesh */

  cs_mesh_location_initialize();
  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();

  /* Define meshes to read */

  cs_user_mesh_input();

  /* Define joining and periodicity parameters if requested
     Must be done before initi1 for the sake of verification */

  cs_gui_mesh_define_joinings();
  cs_user_join();

  cs_gui_mesh_define_periodicities();
  cs_user_periodicity();

  cs_gui_mesh_warping();
  cs_user_mesh_warping();

  /* Call main calculation initialization function or help */

  cs_io_log_initialize();

  cs_field_define_keys_base();
  cs_parameters_define_field_keys();

  cs_preprocessor_data_read_headers(cs_glob_mesh,
                                    cs_glob_mesh_builder);

  /* Initialize Fortran API and calculation setup */

  if ((opts.preprocess | opts.verif) == false && opts.benchmark <= 0) {

    cs_int_t _n_threads = cs_glob_n_threads;
    cs_int_t _rank_id = cs_glob_rank_id, _n_ranks = cs_glob_n_ranks;

    cs_base_fortran_bft_printf_to_f();

    CS_PROCF(csinit, CSINIT)(&_rank_id, &_n_ranks, &_n_threads);

    CS_PROCF(initi1, INITI1)();

    CS_PROCF (haltyp, HALTYP) (&ivoset);
    if (ivoset)
      halo_type = CS_HALO_EXTENDED;

    cs_base_fortran_bft_printf_to_c();

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

  /* Initialize SYRTHES couplings and communication if necessary */

  cs_syr_coupling_all_init();

  /* Initialize Code_Saturne couplings and communication if necessary */

  cs_sat_coupling_all_init();

  /* Set partitioning options */

  {
    int j_id;
    bool join = false;
    bool join_periodic = false;

    for (j_id = 0; j_id < cs_glob_n_joinings; j_id++) {
      if ((cs_glob_join_array[j_id])->param.perio_type == FVM_PERIODICITY_NULL)
        join = true;
      else
        join_periodic = true;
    }

    cs_partition_set_preprocess_hints(join, join_periodic);
    cs_gui_partition();
    cs_user_partition();
  }

  /* Read Preprocessor output */

  cs_preprocessor_data_read_mesh(cs_glob_mesh,
                                 cs_glob_mesh_builder);

  /* Initialize main post-processing */

  cs_gui_postprocess_writers();
  cs_user_postprocess_writers();
  cs_post_init_writers();

  /* Print info on fields and associated keys */

  cs_field_log_defs();
  cs_field_log_key_defs();
  cs_field_log_all_key_vals(true);
  cs_log_printf_flush(CS_LOG_SETUP);

  /* Join meshes / build periodicity links if necessary */

  cs_join_all();

  /* Insert thin walls if necessary */

  cs_user_mesh_thinwall(cs_glob_mesh);

  /* Initialize extended connectivity, ghost cells and other remaining
     parallelism-related structures */

  cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type);
  cs_mesh_update_auxiliary(cs_glob_mesh);

  /* Possible geometry modification */

  cs_user_mesh_modify(cs_glob_mesh);

  /* Discard isolated faces if present */

  cs_post_add_free_faces();
  cs_mesh_discard_free_faces(cs_glob_mesh);

  /* Smoothe mesh if required */

  cs_gui_mesh_smoothe(cs_glob_mesh);
  cs_user_mesh_smoothe(cs_glob_mesh);

  /* Triangulate warped faces if necessary */

  cs_mesh_warping_get_defaults(&cwf_threshold, &cwf_post);

  if (cwf_threshold >= 0.0) {

    t1 = cs_timer_wtime();
    cs_mesh_warping_cut_faces(cs_glob_mesh, cwf_threshold, cwf_post);
    t2 = cs_timer_wtime();

    bft_printf(_("\n Cutting warped faces (%.3g s)\n"), t2-t1);

  }

  /* Now that mesh modification is finished, save mesh if modified */

  cs_user_mesh_save(cs_glob_mesh); /* Disable or force */

  partition_preprocess = cs_partition_get_preprocess();
  if (cs_glob_mesh->modified > 0 || partition_preprocess) {
    if (partition_preprocess) {
      if (cs_glob_mesh->modified > 0)
        cs_mesh_save(cs_glob_mesh, cs_glob_mesh_builder, "mesh_output");
      else
        cs_mesh_to_builder(cs_glob_mesh, cs_glob_mesh_builder, true, NULL);
      cs_partition(cs_glob_mesh, cs_glob_mesh_builder, CS_PARTITION_MAIN);
      cs_mesh_from_builder(cs_glob_mesh, cs_glob_mesh_builder);
      cs_mesh_init_halo(cs_glob_mesh, cs_glob_mesh_builder, halo_type);
      cs_mesh_update_auxiliary(cs_glob_mesh);
    }
    else
      cs_mesh_save(cs_glob_mesh, NULL, "mesh_output");
  }

  /* Destroy the temporary structure used to build the main mesh */

  cs_mesh_builder_destroy(&cs_glob_mesh_builder);

  /* Renumber mesh based on code options */

  cs_user_numbering();

  bft_printf(_("\n Renumbering mesh:\n"));
  bft_printf_flush();
  cs_renumber_mesh(cs_glob_mesh,
                   cs_glob_mesh_quantities);

  /* Initialize group classes */

  cs_mesh_init_group_classes(cs_glob_mesh);

  /* Print info on mesh */

  cs_mesh_print_info(cs_glob_mesh, _("Mesh"));

  /* Compute geometric quantities related to the mesh */

  bft_printf_flush();

  t1 = cs_timer_wtime();
  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_mesh_bad_cells_detect(cs_glob_mesh, cs_glob_mesh_quantities);
  cs_user_mesh_bad_cells_tag(cs_glob_mesh, cs_glob_mesh_quantities);
  t2 = cs_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  /* Initialize selectors and locations for the mesh */
  cs_mesh_init_selectors();
  cs_mesh_location_build(cs_glob_mesh, -1);

  /* Initialize meshes for the main post-processing */

  check_mask = ((opts.preprocess | opts.verif) == true) ? 2 + 1 : 0;

  cs_gui_postprocess_meshes();
  cs_user_postprocess_meshes();
  cs_post_init_meshes(check_mask);

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* For debugging purposes */
  cs_mesh_dump(cs_glob_mesh);
  cs_mesh_quantities_dump(cs_glob_mesh, cs_glob_mesh_quantities);
#endif

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

  if (opts.preprocess == false && opts.benchmark <= 0) {

    /* Check that mesh seems valid */

    cs_mesh_quantities_check_vol(cs_glob_mesh,
                                 cs_glob_mesh_quantities,
                                 (opts.verif ? 1 : 0));

    /* Initialize gradient computation */

    cs_gradient_initialize();
    cs_gradient_perio_initialize();

    if (opts.verif == false) {

      /* Initialize sparse linear systems resolution */

      cs_matrix_initialize();
      cs_sles_initialize();
      cs_multigrid_initialize();

      /* Update Fortran mesh sizes and quantities */

      {
        cs_int_t n_g_cells = cs_glob_mesh->n_g_cells;
        cs_int_t n_g_i_faces = cs_glob_mesh->n_g_i_faces;
        cs_int_t n_g_b_faces = cs_glob_mesh->n_g_b_faces;
        cs_int_t n_g_vertices = cs_glob_mesh->n_g_vertices;
        cs_int_t nthrdi = 1;
        cs_int_t nthrdb = 1;
        cs_int_t ngrpi = 1;
        cs_int_t ngrpb = 1;
        const cs_int_t *idxfi = NULL;
        const cs_int_t *idxfb = NULL;

        if (cs_glob_mesh->i_face_numbering != NULL) {
          const cs_numbering_t *_n = cs_glob_mesh->i_face_numbering;
          nthrdi = _n->n_threads;
          ngrpi = _n->n_groups;
          idxfi = _n->group_index;
        }

        if (cs_glob_mesh->b_face_numbering != NULL) {
          const cs_numbering_t *_n = cs_glob_mesh->b_face_numbering;
          nthrdb = _n->n_threads;
          ngrpb = _n->n_groups;
          idxfb = _n->group_index;
        }

        cs_base_fortran_bft_printf_to_f();

        CS_PROCF (majgeo, MAJGEO)(&(cs_glob_mesh->n_cells),
                                  &(cs_glob_mesh->n_cells_with_ghosts),
                                  &(cs_glob_mesh->n_i_faces),
                                  &(cs_glob_mesh->n_b_faces),
                                  &(cs_glob_mesh->n_vertices),
                                  &(cs_glob_mesh->i_face_vtx_connect_size),
                                  &(cs_glob_mesh->b_face_vtx_connect_size),
                                  &(cs_glob_mesh->n_b_cells),
                                  &n_g_cells,
                                  &n_g_i_faces,
                                  &n_g_b_faces,
                                  &n_g_vertices,
                                  &(cs_glob_mesh->n_families),
                                  &nthrdi,
                                  &nthrdb,
                                  &ngrpi,
                                  &ngrpb,
                                  idxfi,
                                  idxfb,
                                  cs_glob_mesh->i_face_cells,
                                  cs_glob_mesh->b_face_cells,
                                  cs_glob_mesh->b_face_family,
                                  cs_glob_mesh->cell_family,
                                  cs_glob_mesh->i_face_vtx_idx,
                                  cs_glob_mesh->i_face_vtx_lst,
                                  cs_glob_mesh->b_face_vtx_idx,
                                  cs_glob_mesh->b_face_vtx_lst,
                                  cs_glob_mesh->b_cells,
                                  cs_glob_mesh_quantities->b_sym_flag,
                                  &(cs_glob_mesh_quantities->min_vol),
                                  &(cs_glob_mesh_quantities->max_vol),
                                  &(cs_glob_mesh_quantities->tot_vol),
                                  cs_glob_mesh_quantities->cell_cen,
                                  cs_glob_mesh_quantities->i_face_normal,
                                  cs_glob_mesh_quantities->b_face_normal,
                                  cs_glob_mesh_quantities->i_face_cog,
                                  cs_glob_mesh_quantities->b_face_cog,
                                  cs_glob_mesh->vtx_coord,
                                  cs_glob_mesh_quantities->cell_vol,
                                  cs_glob_mesh_quantities->i_face_surf,
                                  cs_glob_mesh_quantities->b_face_surf,
                                  cs_glob_mesh_quantities->i_dist,
                                  cs_glob_mesh_quantities->b_dist,
                                  cs_glob_mesh_quantities->weight,
                                  cs_glob_mesh_quantities->dijpf,
                                  cs_glob_mesh_quantities->diipb,
                                  cs_glob_mesh_quantities->dofij);
      }

      /* Choose between standard and user solver */

      if (cs_user_solver_set() == 0) {

        /*----------------------------------------------
         * Call main calculation function (code Kernel)
         *----------------------------------------------*/

        CS_PROCF(caltri, CALTRI)();

      }
      else {

        /*--------------------------------
         * Call user calculation function
         *--------------------------------*/

        cs_user_solver(cs_glob_mesh,
                       cs_glob_mesh_quantities);

      }

      /* Finalize sparse linear systems resolution */

      cs_multigrid_finalize();
      cs_sles_finalize();
      cs_matrix_finalize();

    }

    /* Finalize gradient computation */

    cs_gradient_perio_finalize();
    cs_gradient_finalize();

    /* Finalize synthetic inlet condition generation */

    cs_inflow_finalize();

  }

  bft_printf(_("\n Destroying structures and ending computation\n"));
  bft_printf_flush();

  CS_PROCF(memfin, MEMFIN)();

  /* Free coupling-related data */

  cs_syr_coupling_all_finalize();
#if defined(HAVE_MPI)
  cs_sat_coupling_all_finalize();
  cs_coupling_finalize();
#endif

  /* Print some mesh statistics */

  cs_gui_usage_log();
  cs_mesh_selector_stats(cs_glob_mesh);

  /* Free cooling towers related structures */

  cs_ctwr_all_destroy();

  /* Free post processing or logging related structures */

  cs_post_finalize();
  cs_log_iteration_destroy_all();

  /* Free GUI-related data */

  cs_gui_particles_free();

  /* Switch logging back to C (may be moved depending on Fortran dependencies) */

  cs_base_fortran_bft_printf_to_c();

  /* Free field info */

  cs_field_destroy_all();
  cs_field_destroy_all_keys();

  /* Free Lagrangian related structures */

  cs_lagr_destroy();

  /* Free main mesh after printing some statistics */

  cs_mesh_location_finalize();
  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* CPU times and memory management finalization */

  cs_all_to_all_log_finalize();
  cs_io_log_finalize();

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

  (void)cs_timer_wtime();

  /* Trap floating-point exceptions */

#if defined(CS_FPE_TRAP)
  if (_fenv_set == 0) {
    if (fegetenv(&_fenv_old) == 0) {
      feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      _fenv_set = 1;
      /* To revert to initial behavior: fesetenv(&_fenv_old); */
    }
  }
#endif

  /* Initialize memory management and signals */

  cs_base_mem_init();
  cs_base_error_init();

  /* Initialize internationalization */

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, cs_base_get_localedir());
  textdomain(PACKAGE);
#endif

  /* Parse command line */

  cs_opts_define(argc, argv, &opts);

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
