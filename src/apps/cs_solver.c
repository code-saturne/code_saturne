/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Main program
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/* On glibc-based systems, define _GNU_SOURCE so as to enable floating-point
   error exceptions; on Itanium, optimized code may raise such exceptions
   due to speculative execution, so we only enable raising of such exceptions
   for code compiled in debug mode, where reduced optimization should not lead
   to such exceptions, and locating the "true" origin of floating-point
   exceptions is helpful. */

#if defined(__linux__) || defined(__linux) || defined(linux)
#if    (!defined(__ia64__) && !defined(__blrts__) && !defined(__bgp__)) \
    || defined(DEBUG)
#define _GNU_SOURCE
#endif
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_GNU_SOURCE)
#include <fenv.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_timer.h>

#if !defined(_GNU_SOURCE)
#include <bft_fp_trap.h>
#endif

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_selector.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_benchmark.h"
#include "cs_base.h"
#include "cs_calcium.h"
#include "cs_coupling.h"
#include "cs_ctwr.h"
#include "cs_gui.h"
#include "cs_io.h"
#include "cs_join.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_solcom.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_warping.h"
#include "cs_mesh_coherency.h"
#include "cs_multigrid.h"
#include "cs_opts.h"
#include "cs_post.h"
#include "cs_preprocessor_data.h"
#include "cs_prototypes.h"
#include "cs_proxy_comm.h"
#include "cs_renumber.h"
#include "cs_restart.h"
#include "cs_sles.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"

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

static int        app_num = -1;
static cs_opts_t  opts;

#if defined(_GNU_SOURCE)
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
  double  t1, t2;

  cs_int_t  iasize, rasize;
  cs_int_t  nituse, nrtuse, nideve, nrdeve;

  int  _verif = -1;

  cs_int_t  *ia = NULL;
  cs_int_t  *ituser = NULL;
  cs_int_t  *idevel = NULL;

  cs_real_t  *ra = NULL;
  cs_real_t  *rtuser = NULL;
  cs_real_t  *rdevel = NULL;

  /* System information */

  cs_base_system_info();

  cs_io_defaults_info();

  /* Initialize global structures for main mesh */

  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();

  /* Call main calculation initialization function or help */

  _verif = (opts.verif == true || opts.benchmark > 0) ? 1 : 0;

  CS_PROCF(initi1, INITI1)(&_verif);

  /* Discover applications visible through MPI (requires communication);
     this is done after main calculation initialization so that the user
     may have the option of assigning a name to this instance. */

#if defined(HAVE_MPI)
  cs_coupling_discover_mpi_apps(app_num, NULL);
#endif

  /* Initialize SYRTHES couplings and communication if necessary */

  cs_syr_coupling_all_init(opts.syr_socket);

  /* Initialize Code_Saturne couplings and communication if necessary */

  cs_sat_coupling_all_init();

  if (opts.ifoenv == 0) {

    /* Read file in obsolete "SolCom" format */

    cs_mesh_solcom_read(cs_glob_mesh,
                        cs_glob_mesh_quantities);

  }
  else {

    /* Read Preprocessor output */

    cs_preprocessor_data_read_mesh(cs_glob_mesh,
                                   cs_glob_mesh_builder);

  }

  /* Initialize main post-processing */

  cs_post_init_main_writer();

  /* Define joining and periodicity parameters if requested */

  CS_PROCF (usjoin, USJOIN)();

  CS_PROCF (usperi, USPERI)();

  /* Join meshes / build periodicity links if necessary */

  cs_join_all();

  /* Initialize ghost cells and other parallelism-related structures */

  cs_mesh_init_halo(cs_glob_mesh);
  cs_mesh_init_parall(cs_glob_mesh);

  /* Possible geometry modification */

  CS_PROCF (usmodg, USMODG)(&(cs_glob_mesh->dim),
                            &(cs_glob_mesh->n_cells_with_ghosts),
                            &(cs_glob_mesh->n_cells),
                            &(cs_glob_mesh->n_i_faces),
                            &(cs_glob_mesh->n_b_faces),
                            &(cs_glob_mesh->n_families),
                            &(cs_glob_mesh->n_max_family_items),
                            &(cs_glob_mesh->n_vertices),
                            &(cs_glob_mesh->i_face_vtx_connect_size),
                            &(cs_glob_mesh->b_face_vtx_connect_size),
                            cs_glob_mesh->i_face_cells,
                            cs_glob_mesh->b_face_cells,
                            cs_glob_mesh->b_face_family,
                            cs_glob_mesh->cell_family,
                            cs_glob_mesh->family_item,
                            cs_glob_mesh->i_face_vtx_idx,
                            cs_glob_mesh->i_face_vtx_lst,
                            cs_glob_mesh->b_face_vtx_idx,
                            cs_glob_mesh->b_face_vtx_lst,
                            cs_glob_mesh->vtx_coord);

  /* Triangulate warped faces if necessary */

  if (opts.cwf == true) {

    t1 = bft_timer_wtime();
    cs_mesh_warping_cut_faces(cs_glob_mesh, opts.cwf_criterion, opts.cwf_post);
    t2 = bft_timer_wtime();

    bft_printf(_("\n Cutting warped faces (%.3g s)\n"), t2-t1);

  }

  /* Renumber mesh based on code options */

  bft_printf(_("\n Renumbering mesh:\n"));
  bft_printf_flush();
  cs_renumber_mesh(cs_glob_mesh,
                   cs_glob_mesh_quantities);

  /* Initialize meshes for the main post-processing */

  cs_post_init_main_meshes();

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

    CS_PROCF (majgeo, MAJGEO)(&(cs_glob_mesh->n_cells),
                              &(cs_glob_mesh->n_cells_with_ghosts),
                              &(cs_glob_mesh->n_i_faces),
                              &(cs_glob_mesh->n_b_faces),
                              &(cs_glob_mesh->n_vertices),
                              &(cs_glob_mesh->i_face_vtx_connect_size),
                              &(cs_glob_mesh->b_face_vtx_connect_size),
                              &n_g_cells,
                              &n_g_i_faces,
                              &n_g_b_faces,
                              &n_g_vertices,
                              &nthrdi,
                              &nthrdb,
                              &ngrpi,
                              &ngrpb,
                              idxfi,
                              idxfb);
  }

  cs_mesh_print_info(cs_glob_mesh, _("Mesh"));

  /* Destroy the temporary structure used to build the main mesh */

  cs_glob_mesh_builder = cs_mesh_builder_destroy(cs_glob_mesh_builder);

  /* Compute geometric quantities related to the mesh */

  bft_printf_flush();

  t1 = bft_timer_wtime();
  cs_mesh_quantities_compute(cs_glob_mesh, cs_glob_mesh_quantities);
  t2 = bft_timer_wtime();

  bft_printf(_("\n Computing geometric quantities (%.3g s)\n"), t2-t1);

  cs_mesh_info(cs_glob_mesh);

  /* Initialize selectors for the mesh */
  cs_mesh_init_selectors();

#if 0
  /* For debugging purposes */
  cs_mesh_dump(cs_glob_mesh);
  cs_mesh_quantities_dump(cs_glob_mesh, cs_glob_mesh_quantities);
#endif

  /* Compute iterations or quality criteria depending on verification options */

  if (opts.verif == true) {
    bft_printf(_("\n Computing quality criteria\n"));
    cs_mesh_quality(cs_glob_mesh, cs_glob_mesh_quantities);
    cs_mesh_coherency_check();
  }

  if (opts.benchmark > 0) {
    int mpi_trace_mode = (opts.benchmark == 2) ? 1 : 0;
    cs_benchmark(mpi_trace_mode);
  }

  if (opts.benchmark <= 0) {

    /* Check that mesh seems valid */

    cs_mesh_quantities_check_vol(cs_glob_mesh,
                                 cs_glob_mesh_quantities);

    /* Allocate Fortran working arrays */

    CS_PROCF(memini, MEMINI)(&iasize, &rasize,
                             &nideve, &nrdeve, &nituse, &nrtuse);

    bft_printf(_("\n"
                 " --- Main Fortran work arrays:\n"
                 "       longia =   %10d (Number of integers)\n"
                 "       longra =   %10d (Number of reals)\n"
                 "       (%d bytes/integer, %d bytes/real)\n"),
               iasize, rasize,
               sizeof(cs_int_t)/sizeof(char),
               sizeof(cs_real_t)/sizeof(char));

    if (nideve > 0 || nrdeve >0)
      bft_printf(_("\n"
                   " --- Developer Fortran work arrays:\n"
                   "       nideve =   %10d (Number of integer)\n"
                   "       nrdeve =   %10d (Number of reals)\n"),
                 nideve, nrdeve);

    bft_printf(_("\n"
                 " --- User Fortran work arrays:\n"
                 "       nituse =   %10d (Number of integers)\n"
                 "       nrtuse =   %10d (Number of reals)\n\n"),
               nituse, nrtuse);

    cs_base_mem_init_work(iasize, rasize, &ia, &ra);

    BFT_MALLOC(ituser, nituse, cs_int_t);
    BFT_MALLOC(rtuser, nrtuse, cs_real_t);

    BFT_MALLOC(idevel, nideve, cs_int_t);
    BFT_MALLOC(rdevel, nrdeve, cs_real_t);

    /* Initialize sparse linear systems resolution */

    cs_sles_initialize();
    cs_multigrid_initialize();

    /*----------------------------------------------
     * Call main calculation function (code Kernel)
     *----------------------------------------------*/

    CS_PROCF(caltri, CALTRI)(&_verif,
                             &nideve, &nrdeve, &nituse, &nrtuse,
                             cs_glob_mesh->i_face_cells,
                             cs_glob_mesh->b_face_cells,
                             cs_glob_mesh->b_face_family,
                             cs_glob_mesh->cell_family,
                             cs_glob_mesh->family_item,
                             cs_glob_mesh->i_face_vtx_idx,
                             cs_glob_mesh->i_face_vtx_lst,
                             cs_glob_mesh->b_face_vtx_idx,
                             cs_glob_mesh->b_face_vtx_lst,
                             idevel, ituser, ia,
                             cs_glob_mesh_quantities->cell_cen,
                             cs_glob_mesh_quantities->i_face_normal,
                             cs_glob_mesh_quantities->b_face_normal,
                             cs_glob_mesh_quantities->i_face_cog,
                             cs_glob_mesh_quantities->b_face_cog,
                             cs_glob_mesh->vtx_coord,
                             cs_glob_mesh_quantities->cell_vol,
                             rdevel, rtuser, ra);

    /* Finalize sparse linear systems resolution */

    cs_multigrid_finalize();
    cs_sles_finalize();

    /* Free working arrays */

    BFT_FREE(ia);
    BFT_FREE(ra);

    BFT_FREE(ituser);
    BFT_FREE(rtuser);

    BFT_FREE(idevel);
    BFT_FREE(rdevel);

  }

  bft_printf(_("\n Destroying structures and ending computation\n"));
  bft_printf_flush();

  /* Free coupling-related data */

  cs_syr_coupling_all_finalize();
#if defined(HAVE_MPI)
  cs_sat_coupling_all_finalize();
  cs_coupling_finalize();
#endif

  /* Free cooling towers related structures */

  cs_ctwr_all_destroy();

  /* Free post processing related structures */

  cs_post_finalize();

  /* Free main mesh */

  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* End of possible communication with a proxy */

  cs_proxy_comm_finalize();

  /* CPU times and memory management finalization */

  cs_restart_print_stats();

  cs_base_bilan_temps();
  cs_base_mem_fin();
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
  app_num = cs_opts_mpi_init(&argc, &argv);
  if (app_num > -1)
    cs_base_mpi_init(app_num);
#endif

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  {
    int t_id;
#pragma omp parallel private(t_id)
    {
      t_id = omp_get_thread_num();
      if (t_id == 0)
        cs_glob_n_threads = omp_get_num_threads();
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

#if defined(ENABLE_NLS)
  bindtextdomain(PACKAGE, LOCALEDIR);
  textdomain(PACKAGE);
#endif

  (void)bft_timer_wtime();

  /* Trap floating-point exceptions */

#if defined(_GNU_SOURCE)
  if (_fenv_set == 0) {
    if (fegetenv(&_fenv_old) == 0) {
      feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      _fenv_set = 1;
      /* To revert to initial behavior: fesetenv(&_fenv_old); */
    }
  }
#else
  bft_fp_trap_set();
#endif

  /* Initialize memory management and signals */

  cs_base_mem_init();
  cs_base_error_init();

  /* Parse command line */

  cs_opts_define(argc, argv, &opts);

  /* Open 'listing' (log) files */

  {
    cs_int_t _n_threads = cs_glob_n_threads;
    cs_int_t _rank_id = cs_glob_rank_id, _n_ranks = cs_glob_n_ranks;

    CS_PROCF(csinit, CSINIT)(&(opts.ifoenv),
                             &_rank_id,
                             &_n_ranks,
                             &_n_threads,
                             &(opts.ilisr0),
                             &(opts.ilisrp));
  }

  cs_base_bft_printf_set();

  /* Log-file header and command line arguments recap */

  cs_opts_logfile_head(argc, argv);

  /* MPI-IO options */

  cs_io_set_defaults(opts.mpi_io_mode);

  /* In case of use with SALOME, optional connection with CFD_Proxy
     launcher or load and start of YACS module */

  /* For verbose output, uncomment the following */
  /*cs_calcium_set_verbosity(3);*/

  /* Using CFD_Proxy, simply initialize connection before standard run */
  if (opts.proxy_socket != NULL) {
    cs_proxy_comm_initialize(opts.proxy_socket,
                             opts.proxy_key,
                             CS_PROXY_COMM_TYPE_SOCKET);
    BFT_FREE(opts.proxy_socket);
    opts.proxy_key = -1;
    cs_calcium_set_comm_proxy();
  }

  /* Running as a standalone SALOME component, load YACS component
     library and run yacsinit() component initialization and event loop,
     which should itself include the standard run routine */

  if (opts.yacs_module != NULL) {
    cs_calcium_load_yacs(opts.yacs_module);
    BFT_FREE(opts.yacs_module);
    cs_calcium_start_yacs(); /* Event-loop does not return as of this version */
    cs_calcium_unload_yacs();
  }

  /* In standard case or with CFD_Proxy, simply call regular run() method */

  else
    cs_run();

  /* Return */

  cs_exit(EXIT_SUCCESS);

  /* Never called, but avoids compiler warning */
  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
