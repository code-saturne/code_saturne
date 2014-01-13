/*============================================================================
 * Code_Aster coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

#include "fvm_nodal.h"

#include "fvm_nodal_extract.h"
#include "fvm_nodal_project.h"
#include "fvm_nodal_triangulate.h"

#include "cs_calcium.h"
#include "cs_interface.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ast_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_ast_coupling_t {

  cs_gnum_t    n_g_faces;
  cs_gnum_t    n_g_nodes;

  /* #if defined(HAVE_MPI) */

  cs_int_t    *n_faces_by_domain;
  cs_int_t    *n_nodes_by_domain;

  cs_int_t    *face_index;
  cs_int_t    *node_index;

  cs_gnum_t   *NUM_GLOB_PTS;

  /* #endif */

};

/*============================================================================
 *  Global variables
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

/*----------------------------------------------------------------------------
 * Send nodes coordinates and structure numbering of coupled mesh.
 *----------------------------------------------------------------------------*/

void CS_PROCF(astgeo, ASTGEO)
(
  cs_int_t   *nbfast,
  cs_int_t   *nbnast,
  cs_int_t   *lstfac,
  cs_int_t   *idfast,
  cs_int_t   *idnast,
  cs_real_t  *almax
)
{
  int j, n_faces;
  int n_nodes;

  fvm_nodal_t  *ifs_mesh;

  cs_int_t *faces_color = NULL;
  cs_int_t *nodes_color = NULL;

  cs_real_t *faces_coords = NULL;
  cs_real_t *nodes_coords = NULL;

  cs_real_t *b_face_cog = cs_glob_mesh_quantities->b_face_cog;

  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

  n_faces = *(nbfast);

  ifs_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            "MaillageExtraitAster_1",
                                            false,
                                            0,
                                            n_faces,
                                            NULL,
                                            lstfac);

  n_nodes = (cs_int_t) fvm_nodal_get_n_entities(ifs_mesh, 0);


  assert(cs_glob_n_ranks == 1);

  /* Creation of the information structure for Code_Saturne/Code_Aster
     coupling */

  BFT_MALLOC(ast_coupling, 1, cs_ast_coupling_t);

  ast_coupling->n_g_nodes = fvm_nodal_get_n_g_vertices(ifs_mesh);
  ast_coupling->n_g_faces = n_faces;

  BFT_MALLOC(ast_coupling->n_faces_by_domain, cs_glob_n_ranks, cs_int_t);
  BFT_MALLOC(ast_coupling->n_nodes_by_domain, cs_glob_n_ranks, cs_int_t);

  ast_coupling->n_nodes_by_domain[0] = n_nodes;
  ast_coupling->n_faces_by_domain[0] = n_faces;


  BFT_MALLOC(faces_color, n_faces, cs_int_t);
  BFT_MALLOC(nodes_color, n_nodes, cs_int_t);
  BFT_MALLOC(faces_coords, 3*n_faces, cs_real_t);
  BFT_MALLOC(nodes_coords, 3*n_nodes, cs_real_t);

  assert(sizeof(cs_coord_t)==sizeof(cs_real_t));

  fvm_nodal_get_vertex_coords(ifs_mesh, CS_INTERLACE,
                              (cs_coord_t *) nodes_coords);

  for (j = 0; j < n_faces; j++) {

    faces_coords[3*j]   = b_face_cog[3*(lstfac[j]-1)];
    faces_coords[3*j+1] = b_face_cog[3*(lstfac[j]-1)+1];
    faces_coords[3*j+2] = b_face_cog[3*(lstfac[j]-1)+2];

    faces_color[j]      = idfast[j];

  }

  for (j = 0; j < n_nodes; j++) {

    nodes_color[j]      = idnast[j];

  }

  ifs_mesh = fvm_nodal_destroy(ifs_mesh);

  /* TODO: Adapter l'algo parallel en sequentiel */

  if (cs_glob_rank_id <= 0) {

    cs_int_t *geom = NULL;

    BFT_MALLOC(geom,2,cs_int_t);

    geom[0] = n_faces;
    geom[1] = n_nodes;

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "DONGEO", 2, &(geom[0]));

    BFT_FREE(geom);

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "ALMAXI", 1, almax);

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "COOFAC", 3*n_faces,&(faces_coords[0]));

    cs_calcium_write_double(comp_id, time_dep, cur_time, 0,
                            "COONOD", 3*n_nodes,&(nodes_coords[0]));

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "COLFAC", n_faces, &(faces_color[0]));

    cs_calcium_write_int(comp_id, time_dep, cur_time, 0,
                         "COLNOD", n_nodes, &(nodes_color[0]));

  }

  cs_glob_ast_coupling = ast_coupling;

  BFT_FREE(faces_color);
  BFT_FREE(nodes_color);
  BFT_FREE(faces_coords);
  BFT_FREE(nodes_coords);
}

/*----------------------------------------------------------------------------
 * Send stresses acting on the fluid/structure interface.
 *----------------------------------------------------------------------------*/

void CS_PROCF(astfor, ASTFOR)
(
 cs_int_t    *ntcast,
 cs_int_t    *nbfast,
 cs_real_t   *forast
)
{
  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

  const int  n_faces = *nbfast;
  const int  n_g_faces = ast_coupling->n_g_faces;

  cs_real_t  *_forast = NULL;


  if (cs_glob_rank_id <= 0)
    BFT_MALLOC(_forast, 3*n_g_faces, cs_real_t);


  if (cs_glob_n_ranks == 1) {

    int  i;

    assert(n_faces == n_g_faces);

    for (i = 0; i < 3*n_g_faces; i++)
      _forast[i] = forast[i];

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    MPI_Gatherv(forast, 3*n_faces, CS_MPI_REAL,
                _forast, ast_coupling->n_faces_by_domain,
                ast_coupling->face_index, CS_MPI_REAL,
                0, cs_glob_mpi_comm);

  }

#endif

  if (cs_glob_rank_id <= 0) {

    cs_calcium_write_double(comp_id, time_dep, cur_time, *ntcast,
                            "FORSAT", 3*n_g_faces, _forast);

    BFT_FREE(_forast);

  }

}

/*----------------------------------------------------------------------------
 * Receive predicted or exact displacement of the fluid/structure interface
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcin, ASTCIN)
(
 cs_int_t    *ntcast,
 cs_int_t    *nbfast,
 cs_int_t    *lstfac,
 cs_real_3_t *depale
)
{
  cs_int_t  i;
  cs_int_t *parent_num      = NULL;

  const cs_int_t n_vertices = cs_glob_mesh->n_vertices;
  const cs_int_t local_rank = (cs_glob_rank_id == -1) ? 0 : cs_glob_rank_id;

  cs_ast_coupling_t  *ast_coupling = cs_glob_ast_coupling;

  int  n_nodes;
  const int  n_g_nodes = ast_coupling->n_g_nodes;
  const int  n_faces = *nbfast;

  fvm_nodal_t *mesh_ifs = NULL;

  cs_real_t  *_xast = NULL, *xast = NULL;

  n_nodes = ast_coupling->n_nodes_by_domain[local_rank];

  BFT_MALLOC(xast, 3*n_nodes, cs_real_t);

  if (cs_glob_rank_id <= 0) {

    int  n_val_read = 0;

    BFT_MALLOC(_xast, 3*n_g_nodes, cs_real_t);

    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           ntcast,
                           "DEPSAT", 3*n_g_nodes,
                           &n_val_read, _xast);

    assert(n_val_read == 3*n_g_nodes);

  }

  if (cs_glob_n_ranks == 1) {

    assert(n_nodes == n_g_nodes);

    for (i = 0; i < 3*n_nodes; i++)
      xast[i] = _xast[i];

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    MPI_Scatterv(_xast, ast_coupling->n_nodes_by_domain,
                 ast_coupling->node_index, CS_MPI_REAL,
                 xast, n_nodes,
                 CS_MPI_REAL, 0, cs_glob_mpi_comm);

  }

#endif

  if (cs_glob_rank_id <= 0)
    BFT_FREE(_xast);

  /* Creation maillage */

  mesh_ifs = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                            "MaillageExtraitAster_1",
                                            false,
                                            0,
                                            n_faces,
                                            NULL,
                                            lstfac);

  /* Tableau de numerotation */

  BFT_MALLOC(parent_num, n_nodes, cs_int_t);
  fvm_nodal_get_parent_num(mesh_ifs, 0, parent_num);

  mesh_ifs = fvm_nodal_destroy(mesh_ifs);

  /* On remplit dans DEPALE les valeurs de deplacement aux noeuds
     couples a deplacement impose */

  for (i = 0; i < n_nodes; i++) {

    cs_lnum_t id_node_parent = parent_num[i] - 1;

    depale[id_node_parent][0] = xast[3*i];
    depale[id_node_parent][1] = xast[3*i + 1];
    depale[id_node_parent][2] = xast[3*i + 2];

  }

  /* Liberation du(des) tableau(x) */
  BFT_FREE(parent_num);
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

    char *instance = NULL;

    double ttanc = 0.;

    BFT_MALLOC(instance,200,char);

    cs_calcium_connect(comp_id, instance);

    BFT_FREE(instance);


    /* Reception des donnees */

    iteration = 0;

    /* commande reception nb iterations */
    cs_calcium_read_int(comp_id, time_dep, &min_time, &max_time,
                        &iteration,
                        "NBPDTM", 1, &n_val_read, nbpdt);
    /* commande reception nb sous-it.  */
    cs_calcium_read_int(comp_id, time_dep, &min_time, &max_time,
                        &iteration,
                        "NBSSIT",1,&n_val_read, nbsspdt);
    /* commande reception de la variable epsilo */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           &iteration,
                           "EPSILO",1,&n_val_read, delta);
    /* commande reception de la variable ttinit */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           &iteration,
                           "TTINIT",1,&n_val_read, &ttanc);
    /* commande reception de la variable dtref */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           &iteration,
                           "PDTREF",1,&n_val_read, dt);

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

    /* Send time step (calculated in "ddtvar" Subroutine) to'Milieu'*/
    cs_calcium_write_double(comp_id, time_dep, 0.0, *nbpdt,
                            "DTSAT", 1, &dttmp);

    /* Receive time step sent by 'Milieu' */
    cs_calcium_read_double(comp_id, time_dep, &min_time, &max_time,
                           nbpdt,
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
 * Send global convergence value of IFS calculations
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

/*----------------------------------------------------------------------------*/

END_C_DECLS
