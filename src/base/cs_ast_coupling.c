/*============================================================================
 * code_aster coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"

#include "cs_all_to_all.h"
#include "cs_calcium.h"
#include "cs_coupling.h"
#include "cs_interface.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_mesh_connect.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"
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

  int          root_rank;

  cs_gnum_t    n_g_faces;
  cs_gnum_t    n_g_vertices;

  cs_lnum_t    n_faces;
  cs_lnum_t    n_vertices;

  cs_lnum_t    *s_vtx_num;

#if defined(HAVE_MPI)

  cs_part_to_block_t *face_p2b;
  cs_all_to_all_t *vtx_b2p;

#endif

  int  verbosity;  /* verbosity */
  int  iteration;  /* 0 for initialization, < 0 for disconnect,
                      iteration from (re)start otherwise */

  int  nbssit;     /* number of sub-iterations */

  double  dt;
  double  dtref;   /* reference time step */
  double  epsilo;  /* scheme convergence threshold */

  int     icv1;
  int     icv2;
  double  lref;

  int     s_it_id; /* Sub-iteration id */

  double  *xast;
  double  *xvast;
  double  *xvasa;
  double  *xastp;

  double  *foras;
  double  *foaas;
  double  *fopas;

};

/*============================================================================
 * Global variables
 *============================================================================*/

cs_ast_coupling_t  *cs_glob_ast_coupling = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate and initialize dynamic vectors (double) based on the 'nb_dyn'
 * number of points.
 *----------------------------------------------------------------------------*/

static void
_allocate_arrays(cs_ast_coupling_t  *ast_cpl)
{
  const cs_lnum_t  nb_dyn = ast_cpl->n_vertices;
  const cs_lnum_t  nb_for = ast_cpl->n_faces;

  BFT_MALLOC(ast_cpl->xast, 3*nb_dyn, double);
  BFT_MALLOC(ast_cpl->xvast, 3*nb_dyn, double);
  BFT_MALLOC(ast_cpl->xvasa, 3*nb_dyn, double);
  BFT_MALLOC(ast_cpl->xastp, 3*nb_dyn, double);

  for (cs_lnum_t k = 0; k < nb_dyn; k++) {

    ast_cpl->xast[3*k]   = 0.;
    ast_cpl->xast[3*k+1] = 0.;
    ast_cpl->xast[3*k+2] = 0.;

    ast_cpl->xvast[3*k]   = 0.;
    ast_cpl->xvast[3*k+1] = 0.;
    ast_cpl->xvast[3*k+2] = 0.;

    ast_cpl->xvasa[3*k]   = 0.;
    ast_cpl->xvasa[3*k+1] = 0.;
    ast_cpl->xvasa[3*k+2] = 0.;

    ast_cpl->xastp[3*k]   = 0.;
    ast_cpl->xastp[3*k+1] = 0.;
    ast_cpl->xastp[3*k+2] = 0.;
  }

  BFT_MALLOC(ast_cpl->foras, 3*nb_for, double);
  BFT_MALLOC(ast_cpl->foaas, 3*nb_for, double);
  BFT_MALLOC(ast_cpl->fopas, 3*nb_for, double);

  for (cs_lnum_t k = 0; k < nb_for; k++) {

    ast_cpl->foras[3*k]   = 0.;
    ast_cpl->foras[3*k+1] = 0.;
    ast_cpl->foras[3*k+2] = 0.;

    ast_cpl->foaas[3*k]   = 0.;
    ast_cpl->foaas[3*k+1] = 0.;
    ast_cpl->foaas[3*k+2] = 0.;

    ast_cpl->fopas[3*k]   = 0.;
    ast_cpl->fopas[3*k+1] = 0.;
    ast_cpl->fopas[3*k+2] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Receives displacements and velocities from code_aster at current time step
 *----------------------------------------------------------------------------*/

static void
_recv_dyn(cs_ast_coupling_t  *ast_cpl)
{
  int n_val_read = 0;

  double *buffer = NULL;

  if (cs_glob_n_ranks <= 1)
    buffer = ast_cpl->xast;

  else if (cs_glob_rank_id <= 0)
    BFT_MALLOC(buffer, 3*ast_cpl->n_g_vertices, double);

  /* Read displacements */

  if (cs_glob_rank_id <= 0) {
    cs_calcium_read_double(ast_cpl->root_rank, &(ast_cpl->iteration),
                           "DEPAST", 3*ast_cpl->n_g_vertices,
                           &n_val_read, buffer);

    assert((cs_gnum_t)n_val_read == 3*ast_cpl->n_g_vertices);
  }


#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    cs_all_to_all_copy_array(ast_cpl->vtx_b2p,
                             CS_DOUBLE,
                             3,
                             true, /* reverse */
                             buffer,
                             ast_cpl->xast);

#endif

  /* Read velocities */

  if (cs_glob_n_ranks <= 1)
    buffer = ast_cpl->xvast;

  if (cs_glob_rank_id <= 0) {
    cs_calcium_read_double(ast_cpl->root_rank, &(ast_cpl->iteration),
                           "VITAST", 3*ast_cpl->n_g_vertices,
                           &n_val_read, buffer);

    assert((cs_gnum_t)n_val_read == 3*ast_cpl->n_g_vertices);
  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    cs_all_to_all_copy_array(ast_cpl->vtx_b2p,
                             CS_DOUBLE,
                             3,
                             true, /* reverse */
                             buffer,
                             ast_cpl->xvast);

#endif

  if (cs_glob_n_ranks > 1)
    BFT_FREE(buffer);
}

/*----------------------------------------------------------------------------
 * Send convergence indicator to code_aster
 *----------------------------------------------------------------------------*/

static void
_send_icv2(cs_ast_coupling_t  *ast_cpl,
           int                 icv)
{
  if (cs_glob_rank_id > 0)
    return;

  cs_calcium_write_int(ast_cpl->root_rank, ast_cpl->iteration,
                       "ICVAST", 1, &icv);
}

/*----------------------------------------------------------------------------
 * Predict displacement or forces based on values of the current and
 * previous time step(s)
 *
 * valpre = c1 * val1 + c2 * val2 + c3 * val3
 *----------------------------------------------------------------------------*/

static void
_pred(double    *valpre,
      double    *val1,
      double    *val2,
      double    *val3,
      double     c1,
      double     c2,
      double     c3,
      cs_lnum_t  n)
{
  if (n < 1)
    return;

  /* Update prediction array */
  for (cs_lnum_t i = 0; i < n; i++) {
    valpre[3*i]     = c1*val1[3*i]     + c2*val2[3*i]     + c3*val3[3*i];
    valpre[(3*i)+1] = c1*val1[(3*i)+1] + c2*val2[(3*i)+1] + c3*val3[(3*i)+1];
    valpre[(3*i)+2] = c1*val1[(3*i)+2] + c2*val2[(3*i)+2] + c3*val3[(3*i)+2];
  }
}

/*----------------------------------------------------------------------------
 * Compute the L2 norm of the difference between vectors vect1 and vect2
 *
 * dinorm = sqrt(sum on nbpts i
 *                 (sum on component j
 *                    ((vect1[i,j]-vect2[i,j])^2)))
 *----------------------------------------------------------------------------*/

static double
_dinorm(double  *vect1,
        double  *vect2,
        double   nbpts)
{
  /* Compute the norm of the difference */
  double norm = 0.;
  for (cs_lnum_t i = 0; i < nbpts; i++) {
    norm += (vect1[3*i]-vect2[3*i])*(vect1[3*i]-vect2[3*i]);
    norm += (vect1[3*i+1]-vect2[3*i+1])*(vect1[3*i+1]-vect2[3*i+1]);
    norm += (vect1[3*i+2]-vect2[3*i+2])*(vect1[3*i+2]-vect2[3*i+2]);
  }

  /* Note that for vertices, vertices at shared parallel boundaries
     will appear multiple tiles, so have a higher "weight" than
     others, but the effect on the global norm should be minor,
     so we avoid a more complex test here */

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    double _norm = norm;
    int    _nbpts = nbpts;
    MPI_Allreduce(&_norm, &norm, 1, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);
    MPI_Allreduce(&_nbpts, &nbpts, 1, MPI_DOUBLE, MPI_SUM,
                  cs_glob_mpi_comm);
  }
#endif

  norm = sqrt(norm/nbpts);
  return norm;
}

/*----------------------------------------------------------------------------
 * Convergence test for implicit calculation case
 *
 * returns:
 *   0 if not converged
 *   1 if     converged
 *----------------------------------------------------------------------------*/

static int
_conv(cs_ast_coupling_t  *ast_cpl,
      int                *icv)
{
  const cs_lnum_t  nb_dyn = ast_cpl->n_vertices;

  /* Local variables */
  int iret;
  double delast = 0.;

  if (ast_cpl->lref > 0.) {

    delast = _dinorm(ast_cpl->xast, ast_cpl->xastp, nb_dyn) / ast_cpl->lref;

    if (ast_cpl->verbosity > 0)
      bft_printf("--------------------------------\n"
                 "convergence test:\n"
                 "delast = %4.2le\n",
                 delast);

    if (delast <= ast_cpl->epsilo) {
      *icv = 1;

      if (ast_cpl->verbosity > 0)
        bft_printf("icv = %d\n"
                   "convergence of sub iteration\n"
                   "----------------------------\n",
                   *icv);
    }
    else {
      if (ast_cpl->verbosity > 0)
        bft_printf("icv = %i\n"
                   "non convergence of sub iteration\n"
                   "--------------------------------\n",
                   *icv);
    }

    iret = 0;
  }
  else {
    bft_printf("Value of lref is negative or zero\n"
               "calculation is aborted\n"
               "---------------------------------\n");
    iret = -1;
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Overwrites data from sub-iteration k-1 with data from sub-iteration k
 * dynamic data: velocities
 * efforts:      forces
 *----------------------------------------------------------------------------*/

static void
_val_ant(cs_ast_coupling_t  *ast_cpl)
{
  const cs_lnum_t  nb_dyn = ast_cpl->n_vertices;
  const cs_lnum_t  nb_for = ast_cpl->n_faces;

  /* record efforts */
  for (cs_lnum_t i = 0; i< nb_for; i++) {
    ast_cpl->foaas[3*i]   = ast_cpl->foras[3*i];
    ast_cpl->foaas[3*i+1] = ast_cpl->foras[3*i+1];
    ast_cpl->foaas[3*i+2] = ast_cpl->foras[3*i+2];
  }

  /* record dynamic data */
  for (cs_lnum_t i = 0; i< nb_dyn; i++) {
    ast_cpl->xvasa[3*i]   = ast_cpl->xvast[3*i];
    ast_cpl->xvasa[3*i+1] = ast_cpl->xvast[3*i+1];
    ast_cpl->xvasa[3*i+2] = ast_cpl->xvast[3*i+2];
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Send vertex coordinates and structure numbering of coupled mesh.
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

  cs_lnum_t *faces_color = NULL;
  cs_lnum_t *vertices_color = NULL;

  cs_real_t *face_centers = NULL;
  cs_real_t *vtx_coords = NULL;

  cs_real_t *b_face_cog = cs_glob_mesh_quantities->b_face_cog;

  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  n_faces = *(nbfast);

  fvm_nodal_t
    *fsi_mesh = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                               "FSI_mesh_1",
                                               false,
                                               0,
                                               n_faces,
                                               NULL,
                                               lstfac);

  /* Creation of the information structure for code_saturne/code_aster
     coupling */

  n_vertices = fvm_nodal_get_n_entities(fsi_mesh, 0);
  ast_cpl->n_vertices = n_vertices;
  ast_cpl->n_g_vertices = fvm_nodal_get_n_g_vertices(fsi_mesh);

  ast_cpl->n_faces = n_faces;
  ast_cpl->n_g_faces = n_faces;

  BFT_MALLOC(ast_cpl->s_vtx_num, ast_cpl->n_vertices, cs_lnum_t);

  fvm_nodal_get_parent_num(fsi_mesh, 0, ast_cpl->s_vtx_num);

  BFT_MALLOC(faces_color, n_faces, cs_lnum_t);
  BFT_MALLOC(vertices_color, n_vertices, cs_lnum_t);
  BFT_MALLOC(face_centers, 3*n_faces, cs_real_t);
  BFT_MALLOC(vtx_coords, 3*n_vertices, cs_real_t);

  assert(sizeof(cs_coord_t)==sizeof(cs_real_t));

  fvm_nodal_get_vertex_coords(fsi_mesh, CS_INTERLACE,
                              (cs_coord_t *)vtx_coords);

  for (j = 0; j < n_faces; j++) {

    cs_lnum_t f_id = lstfac[j] - 1;
    face_centers[3*j]   = b_face_cog[3*f_id];
    face_centers[3*j+1] = b_face_cog[3*f_id+1];
    face_centers[3*j+2] = b_face_cog[3*f_id+2];

    faces_color[j]      = idfast[j];

  }

  for (j = 0; j < n_vertices; j++) {
    vertices_color[j]   = idnast[j];
  }

  /* In parallel, all YACS/Calcium I/O goes through rank 0;
     This is about 1990's level technology/scalability, so a rewrite
     (for example switch to PLE or MedCoupling) would be useful */

#if defined(HAVE_MPI)

  ast_cpl->face_p2b = NULL;
  ast_cpl->vtx_b2p = NULL;

  cs_gnum_t *s_vtx_gnum = NULL;
  cs_part_to_block_t *vtx_p2b = NULL;

  if (cs_glob_n_ranks > 1) {

    /* For faces, which are not shared, use local data, as temporary
       mesh access functions are either based on element type sections,
       or ordered by such. */

    cs_parall_counter(&(ast_cpl->n_g_faces), 1);

    fvm_io_num_t *face_gnum
      = fvm_io_num_create(lstfac,
                          cs_glob_mesh->global_b_face_num,
                          n_faces,
                          0);
    cs_gnum_t *s_face_gnum = fvm_io_num_transfer_global_num(face_gnum);
    assert(ast_cpl->n_g_faces == fvm_io_num_get_global_count(face_gnum));

    face_gnum = fvm_io_num_destroy(face_gnum);

    cs_block_dist_info_t  face_bi;

    face_bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                          cs_glob_n_ranks,
                                          cs_glob_n_ranks, /* all on rank 0 */
                                          0,
                                          ast_cpl->n_g_faces);

    ast_cpl->face_p2b
      = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                        face_bi,
                                        ast_cpl->n_faces,
                                        s_face_gnum);
    cs_part_to_block_transfer_gnum(ast_cpl->face_p2b, s_face_gnum);

    /* For vertices, which may be shared, copy global numbering
       associated with temporary mesh. */

    BFT_MALLOC(s_vtx_gnum, ast_cpl->n_vertices, cs_gnum_t);
    fvm_nodal_get_global_vertex_num(fsi_mesh, s_vtx_gnum);

    cs_block_dist_info_t  vtx_bi;

    vtx_bi = cs_block_dist_compute_sizes(cs_glob_rank_id,
                                         cs_glob_n_ranks,
                                         cs_glob_n_ranks, /* all on rank 0 */
                                         0,
                                         ast_cpl->n_g_vertices);

    ast_cpl->vtx_b2p
      = cs_all_to_all_create_from_block(ast_cpl->n_vertices,
                                        CS_ALL_TO_ALL_USE_DEST_ID,
                                        s_vtx_gnum,
                                        vtx_bi,
                                        cs_glob_mpi_comm);

    vtx_p2b = cs_part_to_block_create_by_gnum(cs_glob_mpi_comm,
                                              vtx_bi,
                                              ast_cpl->n_vertices,
                                              s_vtx_gnum);

  }

#endif

  fsi_mesh = fvm_nodal_destroy(fsi_mesh);

  _allocate_arrays(ast_cpl);

  if (cs_glob_rank_id <= 0) {

    int sizes[2] = {ast_cpl->n_g_faces, ast_cpl->n_g_vertices};

    ast_cpl->lref = *almax;

    bft_printf("\n"
               "----------------------------------\n"
               " Geometric parameters\n"
               "   number of coupled faces: %llu\n"
               "   number of coupled nodes: %llu\n"
               "   reference length (m): %4.2le\n"
               "----------------------------------\n\n",
               (unsigned long long)(ast_cpl->n_g_faces),
               (unsigned long long)(ast_cpl->n_g_vertices),
               ast_cpl->lref);

    /* Directly for code_aster */

    cs_calcium_write_int(ast_cpl->root_rank, 0,
                         "NB_DYN", 1, &(sizes[1]));

    cs_calcium_write_int(ast_cpl->root_rank, 0,
                         "NB_FOR", 1, &(sizes[0]));

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    cs_real_t *g_face_centers = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_face_centers, 3*ast_cpl->n_g_faces, cs_real_t);
    cs_part_to_block_copy_array(ast_cpl->face_p2b,
                                CS_REAL_TYPE,
                                3,
                                face_centers,
                                g_face_centers);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_double(ast_cpl->root_rank, 0,
                              "COOFAC", 3*ast_cpl->n_g_faces, g_face_centers);
      BFT_FREE(g_face_centers);
    }

    cs_real_t *g_vtx_coords = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_vtx_coords, 3*ast_cpl->n_g_vertices, cs_real_t);
    cs_part_to_block_copy_array(vtx_p2b,
                                CS_REAL_TYPE,
                                3,
                                vtx_coords,
                                g_vtx_coords);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_double(ast_cpl->root_rank, 0,
                              "COONOD", 3*ast_cpl->n_g_vertices, g_vtx_coords);
      BFT_FREE(g_vtx_coords);
    }

    int *g_face_color = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_face_color, ast_cpl->n_g_faces, int);
    cs_part_to_block_copy_array(ast_cpl->face_p2b,
                                CS_INT_TYPE,
                                1,
                                faces_color,
                                g_face_color);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_int(ast_cpl->root_rank, 0,
                           "COLFAC", ast_cpl->n_g_faces, g_face_color);
      BFT_FREE(g_face_color);
    }

    int *g_vtx_color = NULL;
    if (cs_glob_rank_id == 0)
      BFT_MALLOC(g_vtx_color, ast_cpl->n_g_vertices, int);
    cs_part_to_block_copy_array(vtx_p2b,
                                CS_INT_TYPE,
                                1,
                                vertices_color,
                                g_vtx_color);
    if (cs_glob_rank_id == 0) {
      cs_calcium_write_int(ast_cpl->root_rank, 0,
                           "COLNOD", ast_cpl->n_g_vertices, g_vtx_color);
      BFT_FREE(g_vtx_color);
    }

    cs_part_to_block_destroy(&vtx_p2b);
    BFT_FREE(s_vtx_gnum);

  }

#endif

  if (cs_glob_n_ranks == 1) {

    cs_calcium_write_double(ast_cpl->root_rank, 0,
                            "COOFAC", 3*n_faces, face_centers);

    cs_calcium_write_double(ast_cpl->root_rank, 0,
                            "COONOD", 3*n_vertices, vtx_coords);

    cs_calcium_write_int(ast_cpl->root_rank, 0,
                         "COLFAC", n_faces, faces_color);

    cs_calcium_write_int(ast_cpl->root_rank, 0,
                         "COLNOD", n_vertices, vertices_color);

  }

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
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  if (ast_cpl->iteration < 0)
    return;

  const cs_lnum_t n_faces = *nbfast;

  const cs_lnum_t  nb_for = n_faces;

  for (cs_lnum_t i = 0; i < 3*n_faces; i++)
    ast_cpl->foras[i] = forast[i];

  /* Send prediction
     (no difference between explicit and implicit cases for forces) */

  cs_real_t alpha = 2.0;
  cs_real_t c1    = alpha;
  cs_real_t c2    = 1-alpha;
  cs_real_t c3    = 0.;

  _pred(ast_cpl->fopas,
        ast_cpl->foras,
        ast_cpl->foaas,
        ast_cpl->foaas,
        c1,
        c2,
        c3,
        nb_for);

  if (ast_cpl->verbosity > 0)
    bft_printf("--------------------------------------\n"
               "Forces prediction coefficients\n"
               " C1: %4.2le\n"
               " C2: %4.2le\n"
               " C3: %4.2le\n"
               "--------------------------------------\n\n",
               c1, c2, c3);

  /* send forces */

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    double  *fopas;
    BFT_MALLOC(fopas, 3*ast_cpl->n_g_faces, double);

    cs_part_to_block_copy_array(ast_cpl->face_p2b,
                                CS_DOUBLE,
                                3,
                                ast_cpl->fopas,
                                fopas);

    if (cs_glob_rank_id <= 0) {
      cs_calcium_write_double(ast_cpl->root_rank, ast_cpl->iteration,
                              "FORAST", 3*ast_cpl->n_g_faces, fopas);
    }

    BFT_FREE(fopas);

  }

#endif

  if (cs_glob_n_ranks <= 1) {
    cs_calcium_write_double(ast_cpl->root_rank, ast_cpl->iteration,
                            "FORAST", 3*ast_cpl->n_g_faces, ast_cpl->fopas);
  }

  /* Second stage (TODO: place in another, better named function) */
  /* ------------------------------------------------------------ */

  /* explicit case: no need fo a convergence test */

  int icv = 1;

  if (ast_cpl->nbssit <= 1) {

    /* handle convergence even when no test is done */
    ast_cpl->icv1 = icv;
    _send_icv2(ast_cpl, icv);

    /* receive displacements from code_aster */
    _recv_dyn(ast_cpl);

    /* save previous values */
    _val_ant(ast_cpl);

  }

  /* implicit case: requires a convergence test */

  else if (ast_cpl->nbssit > 1) {

    /* compute icv */

    int ierr = _conv(ast_cpl, &icv);
    ast_cpl->icv1 = icv;
    icv = ast_cpl->icv2;
    _send_icv2(ast_cpl, icv);

    if ((ast_cpl->s_it_id +1 >= ast_cpl->nbssit) || (icv == 1)) {
      /* receive displacements  computed by code_aster */
      if (ierr >= 0) _recv_dyn(ast_cpl);

      /* then use with code_saturne ? the question remains open... */
      /* if (ierr >= 0) _send2_dyn(); */

      /* receive displacements from code_aster */
      if (ierr >= 0)  _recv_dyn(ast_cpl);
    }
    else {
      ast_cpl->s_it_id += 1;
    }

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
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  if (ast_cpl->iteration < 0)
    return;

  const cs_lnum_t  nb_dyn = ast_cpl->n_vertices;

  /* Predict displacements */

  cs_real_t c1, c2, c3, alpha, beta;

  /* separate prediction for explicit/implicit cases */
  if (ast_cpl->s_it_id == 0) {
    alpha = 0.5;
    beta  = 0.;
    c1    = 1.;
    c2    = (alpha + beta) * cs_glob_time_step->dt[0];
    c3    = -beta * cs_glob_time_step->dt[1];
    _pred(ast_cpl->xastp,
          ast_cpl->xast,
          ast_cpl->xvast,
          ast_cpl->xvasa,
          c1,
          c2,
          c3,
          nb_dyn);
  }
  else if (ast_cpl->s_it_id > 0) {
    alpha = 0.5;
    c1    = alpha;
    c2    = 1. - alpha;
    c3    = 0.;
    _pred(ast_cpl->xastp,
          ast_cpl->xast,
          ast_cpl->xastp,
          ast_cpl->xast,
          c1,
          c2,
          c3,
          nb_dyn);
  }

  if (ast_cpl->verbosity > 0) {

    bft_printf("*********************************\n"
               "*     sub - iteration %i        *\n"
               "*********************************\n\n",
               ast_cpl->s_it_id);

    bft_printf("--------------------------------------------\n"
               "Displacement prediction coefficients\n"
               " C1: %4.2le\n"
               " C2: %4.2le\n"
               " C3: %4.2le\n"
               "--------------------------------------------\n\n",
               c1, c2, c3);

  }

  /* Now set displacments */

  const cs_lnum_t  n_vertices = ast_cpl->n_vertices;

  /* Set in disale the values of prescribed displacements */

  for (cs_lnum_t i = 0; i < n_vertices; i++) {

    cs_lnum_t parent_vtx_id = ast_cpl->s_vtx_num[i] - 1;

    disale[parent_vtx_id][0] = ast_cpl->xastp[3*i];
    disale[parent_vtx_id][1] = ast_cpl->xastp[3*i + 1];
    disale[parent_vtx_id][2] = ast_cpl->xastp[3*i + 2];

  }
}

/*----------------------------------------------------------------------------
 * Exchange time-step
 *----------------------------------------------------------------------------*/

void
CS_PROCF(astpdt, ASTPDT)
(
  cs_real_t *dttab,
  cs_int_t  *nbpdt
)
{
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  /* Update verbosity */

  if (cs_glob_time_step->nt_cur % cs_glob_log_frequency == 0)
    ast_cpl->verbosity = 1;
  else
    ast_cpl->verbosity = 0;

  if (ast_cpl->iteration < 0)
    return;

  cs_real_t  dttmp = ast_cpl->dtref;
  double  dt_ast = dttmp;

  if (ast_cpl->iteration < 0)
    return;

  int err_code = 0;

  ast_cpl->iteration += 1;

  if (cs_glob_rank_id <= 0) {

    double  dt_sat = dttab[0];
    int  n_val_read = 0;

    /* Receive time step sent by code_aster */

    err_code = cs_calcium_read_double(ast_cpl->root_rank, &(ast_cpl->iteration),
                                      "DTAST", 1, &n_val_read, &dt_ast);

    if (err_code >= 0) {

      assert(n_val_read == 1);

      /* Choose smallest time step */

      if (dt_ast < dttmp)
        dttmp = dt_ast;
      if (dt_sat < dttmp)
        dttmp = dt_sat;

      err_code = cs_calcium_write_double(ast_cpl->root_rank, ast_cpl->iteration,
                                         "DTCALC", 1, &dttmp);

    }
    else {

      /* In case of error (probably disconnect) stop at next iteration */

      const cs_time_step_t *ts = cs_glob_time_step;
      if (ts->nt_cur < ts->nt_max + 1)
        cs_time_step_define_nt_max(ts->nt_cur + 1);

      ast_cpl->iteration = -1;

      bft_printf("----------------------------------\n"
                 "code_aster coupling: disconnected (finished) or error\n"
                 "--> stop at end of next time step\n"
                 "----------------------------------\n\n");

    }

  }

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    MPI_Bcast(&dttmp, 1, CS_MPI_REAL, 0, cs_glob_mpi_comm);

#endif

  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  for (cs_lnum_t i = 0; i < n_cells_ext; i++)
    dttab[i] = dttmp;

  ast_cpl->dt = dttmp;

  if (ast_cpl->verbosity > 0)
    bft_printf("----------------------------------\n"
               "reference time step:     %4.21e\n"
               "code_saturne time step:  %4.2le\n"
               "code_aster time step:    %4.2le\n"
               "selected time step:      %4.2le \n"
               "----------------------------------\n\n",
               ast_cpl->dtref, dttab[0], dt_ast, ast_cpl->dt);

  /* Reset sub-iteration count */
  ast_cpl->s_it_id = 0;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial exchange with code_aster
 *
 * \param[in]  nalimx  maximum number of implicitation iterations of
 *                     the structure displacement
 * \param[in]  epalim  relative precision of implicitation of
 *                     the structure displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_initialize(int        nalimx,
                           cs_real_t  epalim)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  int     nbpdtm = ts->nt_max;
  double  ttinit = ts->t_prev;

  /* Allocate global coupling structure */

  cs_ast_coupling_t *ast_cpl;

  BFT_MALLOC(ast_cpl, 1, cs_ast_coupling_t);

  ast_cpl->root_rank = -1;

  ast_cpl->verbosity = 1;
  ast_cpl->iteration = 0; /* < 0 for disconnect */

  ast_cpl->nbssit = nalimx; /* number of sub-iterations */

  ast_cpl->dt = 0.;
  ast_cpl->dtref = ts->dt_ref;  /* reference time step */

  ast_cpl->epsilo = epalim;     /* scheme convergence threshold */

  ast_cpl->icv1 = 0;
  ast_cpl->icv2 = 0;
  ast_cpl->lref = 0.;

  ast_cpl->s_it_id = 0; /* Sub-iteration id */

  ast_cpl->xast = NULL;
  ast_cpl->xvast = NULL;
  ast_cpl->xvasa = NULL;
  ast_cpl->xastp = NULL;

  ast_cpl->foras = NULL;
  ast_cpl->foaas = NULL;
  ast_cpl->fopas = NULL;

  cs_glob_ast_coupling = ast_cpl;

  /* Set Calcium verbosity based on environment variable */

  const char *calcium_verbosity = getenv("CS_CALCIUM_VERBOSITY");
  if (calcium_verbosity != NULL)
    cs_calcium_set_verbosity(atoi(calcium_verbosity));

  /* Find root rank of coupling */

#if defined(PLE_HAVE_MPI)

  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps != NULL) {

    int n_apps = ple_coupling_mpi_set_n_apps(mpi_apps);
    int n_ast_apps = 0;

    /* First pass to count available code_aster couplings */

    for (int i = 0; i < n_apps; i++) {
      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
      if (strncmp(ai.app_type, "code_aster", 10) == 0)
        n_ast_apps += 1;
    }

    /* In single-coupling mode, no identification necessary */

    if (n_ast_apps == 1) {

      for (int i = 0; i < n_apps; i++) {
        const ple_coupling_mpi_set_info_t
          ai = ple_coupling_mpi_set_get_info(mpi_apps, i);
        if (strncmp(ai.app_type, "code_aster", 10) == 0)
          ast_cpl->root_rank = ai.root_rank;
      }

    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "Detected %d code_aster instances; can handle exactly 1.",
                n_ast_apps);

  }

#else

  bft_error(__FILE__, __LINE__, 0,
            "code_aster now requires MPI");

#endif

  /* Calcium  (communication) initialization */

  if (cs_glob_rank_id <= 0) {

    bft_printf(" Send calculation parameters to code_aster\n");

    /* Send data */

    cs_calcium_write_int(ast_cpl->root_rank, 0, "NBPDTM", 1, &nbpdtm);
    cs_calcium_write_int(ast_cpl->root_rank, 0, "NBSSIT", 1,
                         &(ast_cpl->nbssit));

    cs_calcium_write_double(ast_cpl->root_rank, 0, "EPSILO", 1,
                            &(ast_cpl->epsilo));

    /* Send isyncp and ntchr (false, removed function) */
    int isyncp = 0, ntchr = -1;
    cs_calcium_write_int(ast_cpl->root_rank, 0, "ISYNCP", 1, &(isyncp));
    cs_calcium_write_int(ast_cpl->root_rank, 0, "NTCHRO", 1, &(ntchr));

    cs_calcium_write_double(ast_cpl->root_rank, 0, "TTINIT", 1, &ttinit);
    cs_calcium_write_double(ast_cpl->root_rank, 0, "PDTREF", 1,
                            &(ast_cpl->dtref));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize exchange with code_aster
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_finalize(void)
{
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  BFT_FREE(ast_cpl->xast);
  BFT_FREE(ast_cpl->xvast);
  BFT_FREE(ast_cpl->xvasa);
  BFT_FREE(ast_cpl->xastp);

  BFT_FREE(ast_cpl->foras);
  BFT_FREE(ast_cpl->foaas);
  BFT_FREE(ast_cpl->fopas);

#if defined(HAVE_MPI)

  if (ast_cpl->vtx_b2p != NULL)
    cs_all_to_all_destroy(&(ast_cpl->vtx_b2p));
  if (ast_cpl->face_p2b != NULL)
    cs_part_to_block_destroy(&(ast_cpl->face_p2b));

#endif

  BFT_FREE(ast_cpl->s_vtx_num);

  BFT_FREE(ast_cpl);

  cs_glob_ast_coupling = ast_cpl;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Receive convergence value of code_saturne/code_aster coupling
 *
 * \return  convergence indicator computed by coupling scheme
 *          (1: converged, 0: not converged)
 */
/*----------------------------------------------------------------------------*/

int
cs_ast_coupling_get_ext_cvg(void)
{
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    MPI_Bcast(&(ast_cpl->icv1), 1, CS_MPI_INT, 0, cs_glob_mpi_comm);

#endif

  return ast_cpl->icv1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send global convergence value of FSI calculations
 *
 * \param[in]  icved  convergence indicator (1: converged, 0: not converged)
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_send_cvg(int  icved)
{
  cs_ast_coupling_t  *ast_cpl = cs_glob_ast_coupling;

  ast_cpl->icv2 = icved;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
