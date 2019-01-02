/*============================================================================
 * Code_Aster coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"
#include "bft_error.h"

#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"

#include "cs_calcium.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"
#include "cs_block_to_part.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ast_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_ast_coupling_t {

  cs_gnum_t    n_g_faces;
  cs_gnum_t    n_g_vertices;

  cs_lnum_t    n_faces;
  cs_lnum_t    n_vertices;

  cs_lnum_t    *s_vtx_num;

#if defined(HAVE_MPI)

  cs_part_to_block_t *face_p2b;
  cs_block_to_part_t *vtx_b2p;

#endif

};

/*============================================================================
 * Global variables
 *============================================================================*/

cs_ast_coupling_t  *cs_glob_ast_coupling = NULL;

static int  comp_id = 0;

static double  cur_time = 0.;
static double  min_time = 0.;
static double  max_time = 0.;

static cs_calcium_timedep_t  time_dep = CALCIUM_iteration;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Send vertices coordinates and structure numbering of coupled mesh.
 *----------------------------------------------------------------------------*/

void CS_PROCF(astgeo, ASTGEO)
(
  cs_int_t   *nbfast,
  cs_int_t   *lstfac,
  cs_int_t   *idfast,
  cs_int_t   *idnast,
  cs_real_t  *almax
)
{
  cs_lnum_t j, n_faces;
  cs_lnum_t n_vertices;

  fvm_nodal_t  *fsi_mesh;

  cs_lnum_t *faces_color = NULL;
  cs_lnum_t *vertices_color = NULL;

  cs_real_t *face_centers = NULL;
  cs_real_t *vtx_coords = NULL;

  cs_real_t *b_face_cog = cs_glob_mesh_quantities->b_face_cog;

  cs_ast_coupling_t  *ast_coupling = NULL;

  n_faces = *(nbfast);

  fsi_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            "MaillageExtraitAster_1",
                                            false,
                                            0,
                                            n_faces,
                                            NULL,
                                            lstfac);

  /* Creation of the information structure for Code_Saturne/Code_Aster
     coupling */

  BFT_MALLOC(ast_coupling, 1, cs_ast_coupling_t);

  n_vertices = fvm_nodal_get_n_entities(fsi_mesh, 0);
  ast_coupling->n_vertices = n_vertices;
  ast_coupling->n_g_vertices = fvm_nodal_get_n_g_vertices(fsi_mesh);

  ast_coupling->n_faces = n_faces;
  ast_coupling->n_g_faces = n_faces;

  BFT_MALLOC(ast_coupling->s_vtx_num, ast_coupling->n_vertices, cs_lnum_t);

  fvm_nodal_get_parent_num(fsi_mesh, 0, ast_coupling->s_vtx_num);

  BFT_MALLOC(faces_color, n_faces, cs_lnum_t);
  BFT_MALLOC(vertices_color, n_vertices, cs_lnum_t);
  BFT_MALLOC(face_centers, 3*n_faces, cs_real_t);
  BFT_MALLOC(vtx_coords, 3*n_vertices, cs_real_t);

  assert(sizeof(cs_coord_t)==sizeof(cs_real_t));

  fvm_nodal_get_vertex_coords(fsi_mesh, CS_INTERLACE,
                              (cs_coord_t *) vtx_coords);

  for (j = 0; j < n_faces; j++) {

    face_centers[3*j]   = b_face_cog[3*(lstfac[j]-1)];
    face_centers[3*j+1] = b_face_cog[3*(lstfac[j]-1)+1];
    face_centers[3*j+2] = b_face_cog[3*(lstfac[j]-1)+2];

    faces_color[j]      = idfast[j];

  }

  for (j = 0; j < n_vertices; j++) {

    vertices_color[j]      = idnast[j];

  }

  /* In parallel, all YACS/Calcium I/O goes through rank 0;
     This is about 1990's level technology/scalability, but to do better
     (for example switch to PLE or MedCoupling), we'll need to wait for
     more users or more parallism in Code_Aster, as this involves
     both codes. */

#if defined(HAVE_MPI)

  ast_coupling->face_p2b = NULL;
  ast_coupling->vtx_b2p = NULL;

  cs_gnum_t *s_vtx_gnum = NULL;
  cs_part_to_block_t *vtx_p2b = NULL;

  if (cs_glob_n_ranks > 1) {

    /* For faces, which are not shared, use local data, as temporary
       mesh access functions are either based on element type sections,
       or ordered by such. */

    cs_parall_counter(&(ast_coupling->n_g_faces), 1);

    fvm_io_num_t *face_gnum
      = fvm_io_num_create(lstfac,
                          cs_glob_mesh->global_b_face_num,
                          n_faces,
                          0);
    cs_gnum_t *s_face_gnum = fvm_io_num_transfer_global_num(face_gnum);
    assert(ast_coupling->n_g_faces == fvm_io_num_get_global_count(face_gnum));

    face_gnum = fvm_io_num_destroy(face_gnum);

    cs_block_dist_info_t  face_bi;

    face_bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                          cs_glob_n_ranks,
                                          cs_glob_n_ranks, /* all on rank 0 */
                                          0,
                                          ast_coupling->n_g_faces);

    ast_coupling->face_p2b
      = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                        face_bi,
                                        ast_coupling->n_faces,
                                        s_face_gnum);
    cs_part_to_block_transfer_gnum(ast_coupling->face_p2b, s_face_gnum);

    /* For vertices, which may be shared, copy global numbering
       associated with temporary mesh. */

    ast_coupling->n_g_vertices = fvm_nodal_get_n_g_vertices(fsi_mesh);

    BFT_MALLOC(s_vtx_gnum, ast_coupling->n_vertices, cs_gnum_t);
    fvm_nodal_get_global_vertex_num(fsi_mesh, s_vtx_gnum);

    cs_block_dist_info_t  vtx_bi;

    vtx_bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                         cs_glob_n_ranks,
                                         cs_glob_n_ranks, /* all on rank 0 */
                                         0,
                                         ast_coupling->n_g_vertices);

    ast_coupling->vtx_b2p
      = cs_block_to_part_create_by_gnum(cs_glob_mpi_comm,
                                        vtx_bi,
                                        ast_coupling->n_vertices,
                                        s_vtx_gnum);

    vtx_p2b = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                              vtx_bi,
                                              ast_coupling->n_vertices,
                                              s_vtx_gnum);

  }

#endif

  fsi_mesh = fvm_nodal_destroy(fsi_mesh);

  if (cs_glob_rank_id <= 0) {

    int sizes[2] = {ast_coupling->n_g_faces, ast_coupling->n_g_vertices};

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "DONGEO", 2, &(sizes[0]));

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "ALMAXI", 1, almax);

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_real_t *g_face_centers = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_face_centers, 3*ast_coupling->n_g_faces, cs_real_t);
    cs_part_to_block_copy_array(ast_coupling->face_p2b,
                                CS_REAL_TYPE,
                                3,
                                face_centers,
                                g_face_centers);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_double(comp_id, time_dep, cur_time, 0, "COOFAC",
                              3*ast_coupling->n_g_faces, g_face_centers);
      BFT_FREE(g_face_centers);
    }

    cs_real_t *g_vtx_coords = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_vtx_coords, 3*ast_coupling->n_g_vertices, cs_real_t);
    cs_part_to_block_copy_array(vtx_p2b,
                                CS_REAL_TYPE,
                                3,
                                vtx_coords,
                                g_vtx_coords);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_double(comp_id, time_dep, cur_time, 0, "COONOD",
                              3*ast_coupling->n_g_vertices, g_vtx_coords);
      BFT_FREE(g_vtx_coords);
    }

    int *g_face_color = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_face_color, ast_coupling->n_g_faces, int);
    cs_part_to_block_copy_array(ast_coupling->face_p2b,
                                CS_INT_TYPE,
                                1,
                                faces_color,
                                g_face_color);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_int(comp_id, time_dep, cur_time, 0, "COLFAC",
                           ast_coupling->n_g_faces, g_face_color);
      BFT_FREE(g_face_color);
    }

    int *g_vtx_color = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_vtx_color, ast_coupling->n_g_vertices, int);
    cs_part_to_block_copy_array(vtx_p2b,
                                CS_INT_TYPE,
                                1,
                                vertices_color,
                                g_vtx_color);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_int(comp_id, time_dep, cur_time, 0, "COLNOD",
                           ast_coupling->n_g_vertices, g_vtx_color);
      BFT_FREE(g_vtx_color);
    }

    cs_part_to_block_destroy(&vtx_p2b);
    BFT_FREE(s_vtx_gnum);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "COOFAC", 3*n_faces, face_centers);

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "COONOD", 3*n_vertices, vtx_coords);

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "COLFAC", n_faces, faces_color);

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "COLNOD", n_vertices, vertices_color);

  }

  cs_glob_ast_coupling = ast_coupling;

  BFT_FREE(faces_color);
  BFT_FREE(vertices_color);
  BFT_FREE(face_centers);
  BFT_FREE(vtx_coords);
}

/*----------------------------------------------------------------------------
 * Send stresses acting on the fluid/structure interface.
 *----------------------------------------------------------------------------*/

void CS_PROCF(astfor, ASTFOR)
(
 cs_lnum_t    *ntcast,
 cs_lnum_t    *nbfast,
 cs_real_t    *forast
)
{
  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

  const cs_lnum_t n_faces = *nbfast;
  const cs_gnum_t  n_g_faces = ast_coupling->n_g_faces;

  cs_real_t  *g_forast = NULL;

  if (cs_glob_rank_id < 1)
    BFT_MALLOC(g_forast, 3*n_g_faces, cs_real_t);

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    cs_part_to_block_copy_array(ast_coupling->face_p2b,
                                CS_REAL_TYPE,
                                3,
                                forast,
                                g_forast);

#endif

  if (cs_glob_n_ranks == 1) {
    assert((cs_gnum_t)n_faces == n_g_faces);
    for (cs_lnum_t i = 0; i < 3*n_faces; i++)
      g_forast[i] = forast[i];
  }

  if (cs_glob_rank_id < 1) {
    cs_calcium_write_double(comp_id, time_dep, cur_time, *ntcast,
                            "FORSAT", 3*n_g_faces, g_forast);
    BFT_FREE(g_forast);
  }
}

/*----------------------------------------------------------------------------
 * Receive predicted or exact displacement of the fluid/structure interface
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcin, ASTCIN)
(
 cs_int_t    *ntcast,
 cs_real_3_t *disale
)
{
  cs_lnum_t  i;

  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

  const int  n_vertices = ast_coupling->n_vertices;
  const int  n_g_vertices = ast_coupling->n_g_vertices;

  cs_real_t  *g_xast = NULL, *xast = NULL;

  BFT_MALLOC(xast, 3*n_vertices, cs_real_t);

  if (cs_glob_rank_id <= 0) {

    int  n_val_read = 0;

    BFT_MALLOC(g_xast, 3*n_g_vertices, cs_real_t);

    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           ntcast,
                           "DEPSAT", 3*n_g_vertices,
                           &n_val_read, g_xast);

    assert(n_val_read == 3*n_g_vertices);

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    cs_block_to_part_copy_array(ast_coupling->vtx_b2p,
                                CS_REAL_TYPE,
                                3,
                                g_xast,
                                xast);

#endif

  if (cs_glob_n_ranks == 1) {
    assert(n_vertices == n_g_vertices);
    for (i = 0; i < 3*n_vertices; i++)
      xast[i] = g_xast[i];
  }

  if (cs_glob_rank_id <= 0)
    BFT_FREE(g_xast);

  /* On remplit dans DEPALE les valeurs de deplacement aux noeuds
     couples a deplacement impose */

  for (i = 0; i < n_vertices; i++) {

    cs_lnum_t parent_vtx_id = ast_coupling->s_vtx_num[i] - 1;

    disale[parent_vtx_id][0] = xast[3*i];
    disale[parent_vtx_id][1] = xast[3*i + 1];
    disale[parent_vtx_id][2] = xast[3*i + 2];

  }

  /* Liberation du(des) tableau(x) */
  BFT_FREE(xast);
}


/*----------------------------------------------------------------------------
 * Receive coupling parameters
 *----------------------------------------------------------------------------*/

void CS_PROCF(astpar, ASTPAR)
(
 cs_int_t  *nbpdt,
 cs_int_t  *nbsspdt,
 cs_real_t *delta,
 cs_real_t *tt,
 cs_real_t *dt
)
{
  /* Initialisation calcium */

  if (cs_glob_rank_id <= 0) {

    int  n_val_read = 0;
    int  iteration = *nbpdt;

    /* Initialization of communication with calcium */

    double ttanc = 0.;

    char instance[200];

    cs_calcium_connect(comp_id, instance);

    /* Reception des donnees */

    iteration = 0;

    /* commande reception nb iterations */
    cs_calcium_read_int(comp_id, time_dep, &min_time, &max_time, &iteration,
                        "NBPDTM", 1, &n_val_read, nbpdt);
    /* commande reception nb sous-it.  */
    cs_calcium_read_int(comp_id, time_dep, &min_time, &max_time, &iteration,
                        "NBSSIT", 1, &n_val_read, nbsspdt);
    /* commande reception de la variable epsilo */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time, &iteration,
                           "EPSILO", 1, &n_val_read, delta);
    /* commande reception de la variable ttinit */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time, &iteration,
                           "TTINIT", 1, &n_val_read, &ttanc);
    /* commande reception de la variable dtref */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time, &iteration,
                           "PDTREF", 1, &n_val_read, dt);

    if (fabs(*tt - ttanc) > 1.e-16)
      bft_error(__FILE__, __LINE__, 0,
                "Arret du calcul: ttinit different de ttpabs \n");

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    MPI_Bcast(nbpdt,   1, CS_MPI_INT,  0, cs_glob_mpi_comm);
    MPI_Bcast(nbsspdt, 1, CS_MPI_INT,  0, cs_glob_mpi_comm);
    MPI_Bcast(delta,   1, CS_MPI_REAL, 0, cs_glob_mpi_comm);
    MPI_Bcast(dt,      1, CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif

  /* Warning */

  bft_printf("@                                                          \n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@                                                          \n"
               "@ @@ ATTENTION : MODIFICATION DES PARAMETRES UTILISATEURS  \n"
               "@    *********                                             \n"
               "@                                                          \n"
               "@    Presence du couplage Code_Saturne/Code_Aster :        \n"
               "@    Les donnees rentrees dans l'outil 'Milieu'            \n"
               "@    ecrasent les donnees rentrees par l'utilisateur       \n"
               "@                                                          \n"
               "@   Nouvelles valeurs:                                     \n"
               "@      NTMABS = %i                                         \n"
               "@      NALIMX = %i                                         \n"
               "@      EPALIM = %f                                         \n"
               "@      DTREF  = %f                                         \n"
               "@                                                          \n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@                                                          \n",
             *nbpdt, *nbsspdt, *delta, *dt);
}

/*----------------------------------------------------------------------------
 * Exchange time-step
 *----------------------------------------------------------------------------*/

void CS_PROCF(astpdt, ASTPDT)
(
  cs_real_t *dttab,
  cs_int_t  *ncelet,
  cs_int_t  *nbpdt
)
{
  int        i;
  cs_real_t  dttmp = 0.;

  if (cs_glob_rank_id <= 0) {

    int  n_val_read = 0;
    double  dtloc = 0.;

    dttmp = dttab[0];

    /* Send time step (calculated in "ddtvar" Subroutine) to 'Milieu'*/
    cs_calcium_write_double(comp_id, time_dep, 0.0, *nbpdt,
                            "DTSAT", 1, &dttmp);

    /* Receive time step sent by 'Milieu' */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time, nbpdt,
                           "DTCALC", 1, &n_val_read, &(dtloc));

    assert(n_val_read == 1);

    dttmp = dtloc;

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    MPI_Bcast(&dttmp, 1, CS_MPI_REAL, 0, cs_glob_mpi_comm);

#endif


  for (i = 0; i < *ncelet; i++)
    dttab[i] = dttmp;


  bft_printf("@                                                          \n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@                                                          \n"
               "@ @@ ATTENTION : MODIFICATION DE LA VALEUR DU PAS DE TEMPS \n"
               "@    *********                                             \n"
               "@                                                          \n"
               "@  Presence du couplage Saturne/Aster:                     \n"
               "@  les options :                                           \n"
               "@  - pdt uniforme et constant (IDTVAR=0)                   \n"
               "@  - pdt uniforme en espace et variable en temps (IDTVAR=1)\n"
               "@  restent activables                                      \n"
               "@                                                          \n"
               "@  l' option :                                             \n"
               "@  - pdt  variable en espace et en temps  (IDTVAR=2)       \n"
               "@  est desactivee                                          \n"
               "@                                                          \n"
               "@  Valeur du pas de temps retenue pour le calcul couple:   \n"
               "@  dt = %f                                                 \n"
               "@                                                          \n"
               "@                                                          \n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@                                                          \n",
               dttmp);

}

/*----------------------------------------------------------------------------
 * Receive convergence value of Code_Saturne/Code_Aster coupling
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcv1, ASTCV1)
(
 cs_int_t  *ntcast,
 cs_int_t  *icv
)
{
  if (cs_glob_rank_id <= 0) {

    int  n_val_read = 0;

    cs_calcium_read_int(comp_id, time_dep, &min_time, &max_time, ntcast,
                        "ICVEXT", 1, &n_val_read, icv);

    assert(n_val_read == 1);

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    MPI_Bcast(icv, 1, CS_MPI_INT, 0, cs_glob_mpi_comm);

#endif
}

/*----------------------------------------------------------------------------
 * Send global convergence value of FSI calculations
 * (internal and external structures)
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcv2, ASTCV2)
(
 cs_int_t  *ntcast,
 cs_int_t  *icv
)
{
  if (cs_glob_rank_id <= 0) {

    cs_calcium_write_int(comp_id, time_dep, cur_time, *ntcast,
                         "ICV", 1, icv);

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Free coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_ast_coupling_finalize(void)
{
  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

#if defined(HAVE_MPI)

  if (ast_coupling->vtx_b2p != NULL)
    cs_block_to_part_destroy(&(ast_coupling->vtx_b2p));
  if (ast_coupling->face_p2b != NULL)
    cs_part_to_block_destroy(&(ast_coupling->face_p2b));

#endif

  BFT_FREE(ast_coupling->s_vtx_num);

  BFT_FREE(ast_coupling);

  cs_glob_ast_coupling = ast_coupling;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
