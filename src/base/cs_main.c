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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <locale.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_config.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_fp_trap.h>
#include <bft_timer.h>

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
#include "cs_ecs_messages.h"
#include "cs_gui.h"
#include "cs_io.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_solcom.h"
#include "cs_mesh_quality.h"
#include "cs_mesh_warping.h"
#include "cs_mesh_coherency.h"
#include "cs_multigrid.h"
#include "cs_opts.h"
#include "cs_post.h"
#include "cs_prototypes.h"
#include "cs_proxy_comm.h"
#include "cs_renumber.h"
#include "cs_sles.h"
#include "cs_suite.h"
#include "cs_sat_coupling.h"
#include "cs_syr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Main program
 *============================================================================*/

int
main(int    argc,
     char  *argv[])
{
  double  t1, t2;

  cs_int_t  iasize, rasize;
  cs_int_t  nituse, nrtuse, nideve, nrdeve;

  cs_opts_t  opts;

  int  app_num = -1;
  int  _verif = -1;

  cs_int_t  *ia = NULL;
  cs_int_t  *ituser = NULL;
  cs_int_t  *idevel = NULL;

  cs_real_t  *ra = NULL;
  cs_real_t  *rtuser = NULL;
  cs_real_t  *rdevel = NULL;

  /* First analysis of the command line to determine if MPI is required,
     and MPI initialization if it is. */

#if defined(_CS_HAVE_MPI)
  app_num = cs_opts_mpi_app_num(&argc, &argv);
  if (app_num > -1)
    cs_base_mpi_init(&argc, &argv, app_num);
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

  bft_fp_trap_set();

  /* Initialize memory management and signals */

  cs_base_mem_init();
  cs_base_error_init();

  /* Parse command line */

  cs_opts_define(argc, argv, &opts);

  /* Open 'listing' (log) files */

  CS_PROCF(csinit, CSINIT)(&(opts.ifoenv),
                           &cs_glob_base_rang,
                           &cs_glob_base_nbr,
                           &(opts.ilisr0),
                           &(opts.ilisrp));
  cs_base_bft_printf_set();

  /* Log-file header and command line arguments recap */

  cs_opts_logfile_head(argc, argv);

  /* Connection with CFD_Proxy launcher */

  if (opts.proxy_socket != NULL) {
    cs_proxy_comm_initialize(opts.proxy_socket,
                             opts.proxy_key,
                             CS_PROXY_COMM_TYPE_SOCKET);
    BFT_FREE(opts.proxy_socket);
    opts.proxy_key = -1;
    cs_calcium_set_comm_proxy();
  }

  /* System information */

  cs_base_system_info();

  /* Initialize global structures for main mesh */

  cs_glob_mesh = cs_mesh_create();
  cs_glob_mesh_builder = cs_mesh_builder_create();
  cs_glob_mesh_quantities = cs_mesh_quantities_create();

  /* Initialize reading of Preprocessor output */

  if (opts.ifoenv != 0) {
#if defined(FVM_HAVE_MPI)
    cs_glob_pp_io = cs_io_initialize("preprocessor_output",
                                     "Face-based mesh definition, R0",
                                     CS_IO_MODE_READ,
                                     0,
                                     CS_IO_ECHO_OPEN_CLOSE,
                                     cs_glob_base_mpi_comm);
#else
    cs_glob_pp_io = cs_io_initialize("preprocessor_output",
                                     "Face-based mesh definition, R0",
                                     CS_IO_MODE_READ,
                                     CS_IO_ECHO_OPEN_CLOSE,
                                     -1);
#endif
  }

  /* Call main calculation initialization function or help */

  _verif = opts.iverif;
  if (opts.benchmark > 0 && _verif < 0)
    _verif = 0;

  CS_PROCF(initi1, INITI1)(&_verif);

  /* Discover applications visible through MPI (requires communication);
     this is done after main calculation initialization so that the user
     may have the option of assigning a name to this instance. */

#if defined(_CS_HAVE_MPI)
  cs_coupling_discover_mpi_apps(app_num, NULL);
#endif

  /* Initialize SYRTHES couplings and communication if necessary */

  cs_syr_coupling_all_init();

  if (opts.ifoenv == 0) {

    /* Read file in obsolete "SolCom" format */

    cs_mesh_solcom_read(cs_glob_mesh,
                        cs_glob_mesh_quantities);

  }
  else {

    /* Read Preprocessor output */

    cs_ecs_messages_read_data(cs_glob_mesh,
                              cs_glob_mesh_builder);

  }

  /* Initialize main post-processing */

  cs_post_init_pcp_writer();

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

  cs_post_init_pcp_maillages();

  /* Update some mesh sizes counts and quantities */

  {
    cs_int_t  n_g_cells, n_g_i_faces, n_g_b_faces, n_g_vertices;

    n_g_cells = cs_glob_mesh->n_g_cells;
    n_g_i_faces = cs_glob_mesh->n_g_i_faces;
    n_g_b_faces = cs_glob_mesh->n_g_b_faces;
    n_g_vertices = cs_glob_mesh->n_g_vertices;

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
                              &n_g_vertices);
  }

  cs_mesh_print_info(cs_glob_mesh);

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

  if (opts.iverif == 0) {
    bft_printf(_("\n Computing quality criteria\n"));
    cs_mesh_quality(cs_glob_mesh, cs_glob_mesh_quantities);
  }

  if (opts.iverif >= 0)
    cs_mesh_coherency_check();

  if (opts.benchmark > 0) {
    int mpi_trace_mode = (opts.benchmark == 2) ? 1 : 0;
    cs_benchmark(mpi_trace_mode);
  }

  if (opts.iverif != 0 && opts.benchmark <= 0) {

    /* Allocate Fortran working arrays */

    CS_PROCF(memini, MEMINI)(&iasize, &rasize,
                             &nideve, &nrdeve, &nituse, &nrtuse);

    bft_printf(_("\n"
                 " --- Main Fortran work arrays:\n"
                 "       LONGIA =   %10d (Number of integers)\n"
                 "       LONGRA =   %10d (Number of reals)\n"
                 "       (%d bytes/integer, %d bytes/real)\n"),
               iasize, rasize,
               sizeof(cs_int_t)/sizeof(char),
               sizeof(cs_real_t)/sizeof(char));

    if (nideve > 0 || nrdeve >0)
      bft_printf(_("\n"
                   " --- Developer Fortran work arrays:\n"
                   "       NIDEVE =   %10d (Number of integer)\n"
                   "       NRDEVE =   %10d (Number of reals)\n"),
                 nideve, nrdeve);

    bft_printf(_("\n"
                 " --- User Fortran work arrays:\n"
                 "       NITUSE =   %10d (Number of integers)\n"
                 "       NRTUSE =   %10d (Number of reals)\n\n"),
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

    CS_PROCF(caltri, CALTRI)(&(opts.iverif),
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
#if defined(_CS_HAVE_MPI)
  cs_sat_coupling_all_finalize();
  cs_coupling_finalize();
#endif

  /* Free cooling towers related structures */

  cs_ctwr_all_destroy();

  /* Free post processing related structures */

  cs_post_detruit();

  /* Free main mesh */

  cs_mesh_quantities_destroy(cs_glob_mesh_quantities);
  cs_mesh_destroy(cs_glob_mesh);

  /* End of possible communication with a proxy */

  cs_proxy_comm_finalize();

  /* CPU times and memory management finalization */

  cs_restart_print_stats();

  cs_base_bilan_temps();
  cs_base_mem_fin();

  /* Return */

  cs_exit(EXIT_SUCCESS);

  /* Never called, but avoids compiler warning */
  return 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
