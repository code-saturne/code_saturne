/*============================================================================
 * code_aster coupling
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_paramedmem_coupling.h"
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

#if !defined(HAVE_MPI)

/* Fake structure for compilation without MPI (unusable in current form) */

typedef struct {
  int          root_rank; /* Application root rank in MPI_COMM_WORLD */
} ple_coupling_mpi_set_info_t;

#endif

/* Main code_aster coupling structure */

struct _cs_ast_coupling_t {

  ple_coupling_mpi_set_info_t  aci;  /* code_aster coupling info */

  cs_lnum_t    n_faces;       /* Local number of coupled faces */
  cs_lnum_t    n_vertices;    /* Local number of coupled vertices */

  cs_gnum_t    n_g_faces;     /* Global number of coupled faces */
  cs_gnum_t    n_g_vertices;  /* Global number of coupld vertices */

#if defined(HAVE_PARAMEDMEM)

  cs_paramedmem_coupling_t  *mc_faces;
  cs_paramedmem_coupling_t  *mc_vertices;

#endif

  int  verbosity;      /* verbosity level */
  int  visualization;  /* visualization level */

  fvm_nodal_t  *post_mesh;     /* Optional mesh for post-processing output */
  int           post_mesh_id;  /* 0 if post-processing is not active,
                                  or post-processing mesh_id (< 0) */

  int  iteration;      /* 0 for initialization, < 0 for disconnect,
                          iteration from (re)start otherwise */

  int  nbssit;         /* number of sub-iterations */

  double  dt;
  double  dtref;       /* reference time step */
  double  epsilo;      /* scheme convergence threshold */

  int     icv1;        /* Convergence indicator */
  int     icv2;        /* Convergence indicator (final) */

  double  lref;        /* Characteristic macroscopic domain length */

  int     s_it_id;     /* Sub-iteration id */

  double  *xast;       /* Mesh displacement last received (current iteration) */
  double  *xvast;      /* Mesh velocity last received (current iteration) */
  double  *xvasa;      /* Mesh displacement at previous sub-iteration */
  double  *xastp;      /* Mesh velocity at previous sub-iteration */

  double  *foras;      /* Fluid forces at current sub-iteration */
  double  *foaas;      /* Fluid forces at previous sub-iteration */
  double  *fopas;      /* Predicted fluid forces */

};

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _name_f_f[] = "fluid_forces";
static const char _name_m_d[] = "mesh_displacement";
static const char _name_m_v[] = "mesh_velocity";

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
 * Scatter values of type cs_real_3_t (tuples) based on indirection list
 *
 * parameters:
 *   n_elts   <-- number of elements
 *   elt_ids  <-- element ids, or NULL
 *   v_in     <-- input values, on elt_ids location
 *   v_out    <-> output values, on parent location
 *----------------------------------------------------------------------------*/

static void
_scatter_values_r3(cs_lnum_t          n_elts,
                   const cs_lnum_t    elt_ids[],
                   const cs_real_3_t  v_in[],
                   cs_real_3_t        v_out[])
{
  if (elt_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];

      v_out[j][0] = v_in[i][0];
      v_out[j][1] = v_in[i][1];
      v_out[j][2] = v_in[i][2];
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      v_out[i][0] = v_in[i][0];
      v_out[i][1] = v_in[i][1];
      v_out[i][2] = v_in[i][2];
    }
  }
}

/*----------------------------------------------------------------------------
 * Receive displacements and velocities from code_aster at current time step
 *----------------------------------------------------------------------------*/

static void
_recv_dyn(cs_ast_coupling_t  *ast_cpl)
{
#if defined(HAVE_PARAMEDMEM)

  int verbosity = (cs_log_default_is_active()) ? ast_cpl->verbosity : 0;

  cs_paramedmem_attach_field_by_name(ast_cpl->mc_vertices, _name_m_d);
  cs_paramedmem_attach_field_by_name(ast_cpl->mc_vertices, _name_m_v);

  if (verbosity > 1) {
    bft_printf
      (_("code_aster: starting MEDCouping receive of values "
         "at coupled vertices..."));
    bft_printf_flush();
  }

  cs_paramedmem_recv_data(ast_cpl->mc_vertices);

  if (verbosity > 1) {
    bft_printf(_("[ok]\n"));
    bft_printf_flush();
  }

  cs_paramedmem_field_import_l(ast_cpl->mc_vertices, _name_m_d, ast_cpl->xast);
  cs_paramedmem_field_import_l(ast_cpl->mc_vertices, _name_m_v, ast_cpl->xvast);

#endif /* defined(HAVE_PARAMEDMEM) */

  /* For dry run, reset values to zero to avoid uninitialized values */
  if (ast_cpl->aci.root_rank < 0) {
    const cs_lnum_t  nb_dyn = ast_cpl->n_vertices * 3;
    for (cs_lnum_t k = 0; k < nb_dyn; k++) {
      ast_cpl->xast[k]  = 0.;
      ast_cpl->xvast[k] = 0.;
    }
  }
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

  cs_calcium_write_int(ast_cpl->aci.root_rank, ast_cpl->iteration,
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

  int verbosity = (cs_log_default_is_active()) ? ast_cpl->verbosity : 0;

  delast = _dinorm(ast_cpl->xast, ast_cpl->xastp, nb_dyn) / ast_cpl->lref;

  if (verbosity > 0)
    bft_printf("--------------------------------\n"
               "convergence test:\n"
               "delast = %4.2le\n",
               delast);

  if (delast <= ast_cpl->epsilo) {
    *icv = 1;

    if (verbosity > 0)
      bft_printf("icv = %d\n"
                 "convergence of sub iteration\n"
                 "----------------------------\n",
                 *icv);
  }
  else {
    if (verbosity > 0)
      bft_printf("icv = %i\n"
                 "non convergence of sub iteration\n"
                 "--------------------------------\n",
                 *icv);
  }

  iret = 0;

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

/*----------------------------------------------------------------------------
 * Post process variables associated with code_aster couplings
 *
 * parameters:
 *   coupling        <--  Void pointer to code_aster coupling structure
 *   ts              <--  time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_ast_coupling_post_function(void                  *coupling,
                               const cs_time_step_t  *ts)
{
#if defined(HAVE_PARAMEDMEM)

  const cs_ast_coupling_t  *cpl = coupling;

  if (cpl->post_mesh == NULL)
    return;

  /* Note: since numbering in fvm_nodal_t structures (ordered by
     element type) may not align with the selection order, we need to project
     values on parent faces first */

  const cs_lnum_t *face_ids
    = cs_paramedmem_mesh_get_elt_list(cpl->mc_faces);
  const cs_lnum_t *vtx_ids
    = cs_paramedmem_mesh_get_vertex_list(cpl->mc_vertices);

  cs_real_t *values;
  {
    const cs_mesh_t *m = cs_glob_mesh;
    cs_lnum_t n_vals = CS_MAX(m->n_b_faces, m->n_vertices) * 3;
    BFT_MALLOC(values, n_vals, cs_real_t);
    for (cs_lnum_t i = 0; i < n_vals; i++)
      values[i] = 0.;
  }

  _scatter_values_r3(cpl->n_vertices,
                     vtx_ids,
                     (const cs_real_3_t *)cpl->xast,
                     (cs_real_3_t *)values);

  cs_post_write_vertex_var(cpl->post_mesh_id,
                           CS_POST_WRITER_ALL_ASSOCIATED,
                           "FSI mesh displacement",
                           3,
                           true, /* interlaced */
                           true, /* use parent */
                           CS_POST_TYPE_cs_real_t,
                           values,
                           ts);

  _scatter_values_r3(cpl->n_vertices,
                     vtx_ids,
                     (const cs_real_3_t *)cpl->xvast,
                     (cs_real_3_t *)values);

  cs_post_write_vertex_var(cpl->post_mesh_id,
                           CS_POST_WRITER_ALL_ASSOCIATED,
                           "FSI mesh velocity",
                           3,
                           true, /* interlaced */
                           true, /* on parent */
                           CS_POST_TYPE_cs_real_t,
                           values,
                           ts);

  _scatter_values_r3(cpl->n_faces,
                     face_ids,
                     (const cs_real_3_t *)cpl->foras,
                     (cs_real_3_t *)values);

  cs_post_write_var(cpl->post_mesh_id,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    "Stress",
                    3,
                    true,  /* interlaced */
                    true,  /* on parent */
                    CS_POST_TYPE_cs_real_t,
                    NULL,  /* cell values */
                    NULL,  /* interior face values */
                    values,
                    ts);

  BFT_FREE(values);

#endif /* defined(HAVE_PARAMEDMEM) */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial exchange with code_aster
 *
 * \param[in]  verbosity      verbosity level for code_aster coupling
 * \param[in]  visualization  visualization level for code_aster coupling
 * \param[in]  nalimx         maximum number of implicitation iterations of
 *                            the structure displacement
 * \param[in]  epalim         relative precision of implicitation of
 *                            the structure displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_initialize(int        verbosity,
                           int        visualization,
                           int        nalimx,
                           cs_real_t  epalim)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  int     nbpdtm = ts->nt_max;
  double  ttinit = ts->t_prev;

  /* Allocate global coupling structure */

  cs_ast_coupling_t *cpl;

  BFT_MALLOC(cpl, 1, cs_ast_coupling_t);

  memset(&(cpl->aci), 0, sizeof(ple_coupling_mpi_set_info_t));
  cpl->aci.root_rank = -1;

  cpl->n_faces = 0;
  cpl->n_vertices = 0;

  cpl->n_g_faces = 0;
  cpl->n_g_vertices = 0;

#if defined(HAVE_PARAMEDMEM)

  cpl->mc_faces = NULL;
  cpl->mc_vertices = NULL;

#endif

  cpl->verbosity = verbosity;
  cpl->visualization = visualization;

  cpl->post_mesh = NULL;

  cpl->iteration = 0; /* < 0 for disconnect */

  cpl->nbssit = nalimx; /* number of sub-iterations */

  cpl->dt = 0.;
  cpl->dtref = ts->dt_ref;  /* reference time step */
  cpl->epsilo = epalim;     /* scheme convergence threshold */

  cpl->icv1 = 0;
  cpl->icv2 = 0;
  cpl->lref = 0.;

  cpl->s_it_id = 0; /* Sub-iteration id */

  cpl->xast = NULL;
  cpl->xvast = NULL;
  cpl->xvasa = NULL;
  cpl->xastp = NULL;

  cpl->foras = NULL;
  cpl->foaas = NULL;
  cpl->fopas = NULL;

  cs_glob_ast_coupling = cpl;

  /* Set verbosity based on environment variable */

  const char *verbosity_s = getenv("CS_AST_COUPLING_VERBOSITY");
  if (verbosity_s != NULL) {
    int _verbosity = atoi(verbosity_s);
    cpl->verbosity = _verbosity;
  }
  cs_calcium_set_verbosity(cpl->verbosity);

  /* Find root rank of coupling */

#if defined(PLE_HAVE_MPI) && defined(HAVE_PARAMEDMEM)

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
          cpl->aci = ai;
      }

    }
    else if (n_ast_apps == 0) {
      bft_printf("\n"
                 "Warning: no matching code_aster instance detected.\n"
                 "         dry run in coupling simulation mode.\n");
      bft_printf_flush();
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                "Detected %d code_aster instances; can handle exactly 1.",
                n_ast_apps);

  }

#else

  bft_error(__FILE__, __LINE__, 0,
            "code_aster coupling requires MEDCoupling with MPI support.");

#endif

  /* Calcium  (communication) initialization */

  if (cs_glob_rank_id <= 0) {

    bft_printf(" Send calculation parameters to code_aster\n");

    /* Send data */

    cs_calcium_write_int(cpl->aci.root_rank, 0, "NBPDTM", 1, &nbpdtm);
    cs_calcium_write_int(cpl->aci.root_rank, 0, "NBSSIT", 1,
                         &(cpl->nbssit));

    cs_calcium_write_double(cpl->aci.root_rank, 0, "EPSILO", 1,
                            &(cpl->epsilo));

    /* Send isyncp and ntchr (false, removed function) */
    int isyncp = 0, ntchr = -1;
    cs_calcium_write_int(cpl->aci.root_rank, 0, "ISYNCP", 1, &(isyncp));
    cs_calcium_write_int(cpl->aci.root_rank, 0, "NTCHRO", 1, &(ntchr));

    cs_calcium_write_double(cpl->aci.root_rank, 0, "TTINIT", 1, &ttinit);
    cs_calcium_write_double(cpl->aci.root_rank, 0, "PDTREF", 1,
                            &(cpl->dtref));

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
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

  if (cpl == NULL)
    return;

  BFT_FREE(cpl->xast);
  BFT_FREE(cpl->xvast);
  BFT_FREE(cpl->xvasa);
  BFT_FREE(cpl->xastp);

  BFT_FREE(cpl->foras);
  BFT_FREE(cpl->foaas);
  BFT_FREE(cpl->fopas);

  if (cpl->post_mesh != NULL)
    cpl->post_mesh = fvm_nodal_destroy(cpl->post_mesh);

#if defined(HAVE_PARAMEDMEM)

  cs_paramedmem_coupling_destroy(cpl->mc_vertices);
  cs_paramedmem_coupling_destroy(cpl->mc_faces);

  cpl->mc_vertices = NULL;
  cpl->mc_faces = NULL;

#endif /* defined(HAVE_PARAMEDMEM) */

  BFT_FREE(cpl);

  cs_glob_ast_coupling = cpl;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Extract and exchange mesh information for surfaces coupled with
 *        code_aster.
 *
 * \param[in]  n_faces   number of coupled faces.
 * \param[in]  face_ids  ids of coupled faces (ordered by increasing id)
 * \param[in]  almax     characteristic macroscopic domain length
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_geometry(cs_lnum_t         n_faces,
                         const cs_lnum_t  *face_ids,
                         cs_real_t         almax)
{
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

#if defined(HAVE_PARAMEDMEM)

  if (cpl->aci.root_rank > -1) {
    cpl->mc_faces = cs_paramedmem_coupling_create(NULL,
                                                  cpl->aci.app_name,
                                                  "fsi_face_exchange");
    cpl->mc_vertices = cs_paramedmem_coupling_create(NULL,
                                                     cpl->aci.app_name,
                                                     "fsi_vertices_exchange");
  }
  else {
    cpl->mc_faces
      = cs_paramedmem_coupling_create_uncoupled("fsi_face_exchange");
    cpl->mc_vertices
      = cs_paramedmem_coupling_create_uncoupled("fsi_vertices_exchange");
  }

  cs_paramedmem_add_mesh_from_ids(cpl->mc_faces,
                                  n_faces,
                                  face_ids,
                                  2);

  cs_paramedmem_add_mesh_from_ids(cpl->mc_vertices,
                                  n_faces,
                                  face_ids,
                                  2);

  cpl->n_faces = n_faces;
  cpl->n_vertices = cs_paramedmem_mesh_get_n_vertices(cpl->mc_faces);

#endif

  fvm_nodal_t *fsi_mesh
    = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                     "FSI_mesh_1",
                                     true, /* include families */
                                     0,
                                     n_faces,
                                     NULL,
                                     face_ids);

  cpl->n_g_faces = n_faces;
  cs_parall_counter(&(cpl->n_g_faces), 1);
  cpl->n_g_vertices = fvm_nodal_get_n_g_vertices(fsi_mesh);

  /* Creation of the information structure for code_saturne/code_aster
     coupling */

  assert(sizeof(cs_coord_t)==sizeof(cs_real_t));

  if (cpl->visualization > 0)
    cpl->post_mesh = fsi_mesh;

  else
    fsi_mesh = fvm_nodal_destroy(fsi_mesh);

  _allocate_arrays(cpl);

  if (almax <= 0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: almax = %g, where a positive value is expected.",
              __func__, almax);

  cpl->lref = almax;

  if (cs_glob_rank_id <= 0) {

    bft_printf("\n"
               "----------------------------------\n"
               " Geometric parameters\n"
               "   number of coupled faces: %llu\n"
               "   number of coupled vertices: %llu\n"
               "   reference length (m): %4.2le\n"
               "----------------------------------\n\n",
               (unsigned long long)(cpl->n_g_faces),
               (unsigned long long)(cpl->n_g_vertices),
               cpl->lref);

  }

  /* Define coupled fields */

#if defined(HAVE_PARAMEDMEM)

  cs_paramedmem_def_coupled_field(cpl->mc_vertices,
                                  _name_m_d,
                                  3,
                                  CS_MEDCPL_FIELD_INT_MAXIMUM,
                                  CS_MEDCPL_ON_NODES,
                                  CS_MEDCPL_LINEAR_TIME);

  cs_paramedmem_def_coupled_field(cpl->mc_vertices,
                                  _name_m_v,
                                  3,
                                  CS_MEDCPL_FIELD_INT_MAXIMUM,
                                  CS_MEDCPL_ON_NODES,
                                  CS_MEDCPL_LINEAR_TIME);

  cs_paramedmem_def_coupled_field(cpl->mc_faces,
                                  _name_f_f,
                                  3,
                                  CS_MEDCPL_FIELD_INT_CONSERVATION,
                                  CS_MEDCPL_ON_CELLS,
                                  CS_MEDCPL_LINEAR_TIME);

#endif /* defined(HAVE_PARAMEDMEM) */

  /* Post-processing */

  if (cpl->visualization > 0) {

    const int writer_ids[] = {CS_POST_WRITER_DEFAULT};
    cpl->post_mesh_id = cs_post_get_free_mesh_id();

    cs_post_define_existing_mesh(cpl->post_mesh_id,
                                 cpl->post_mesh,
                                 0,      /* dim_shift */
                                 false,  /* transfer */
                                 false,  /* auto variables */
                                 1,      /* number of associated writers */
                                 writer_ids);

    cs_post_add_time_dep_output(_cs_ast_coupling_post_function,
                              (void *)cpl);

  }
  else
    cpl->post_mesh_id = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Exchange time-step information with code_aster.
 *
 * \param[in, out]  c_dt  time step at each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_exchange_time_step(cs_real_t  c_dt[])
{
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  cs_real_t  dttmp = cpl->dtref;
  double  dt_ast = dttmp;

  if (cpl->iteration < 0)
    return;

  int err_code = 0;

  cpl->iteration += 1;

  if (cs_glob_rank_id <= 0) {

    double  dt_sat = c_dt[0];
    int  n_val_read = 0;

    /* Receive time step sent by code_aster */

    err_code = cs_calcium_read_double(cpl->aci.root_rank, &(cpl->iteration),
                                      "DTAST", 1, &n_val_read, &dt_ast);

    if (err_code >= 0) {

      assert(n_val_read == 1);

      /* Choose smallest time step */

      if (dt_ast < dttmp)
        dttmp = dt_ast;
      if (dt_sat < dttmp)
        dttmp = dt_sat;

      err_code = cs_calcium_write_double(cpl->aci.root_rank, cpl->iteration,
                                         "DTCALC", 1, &dttmp);

    }
    else {

      /* In case of error (probably disconnect) stop at next iteration */

      const cs_time_step_t *ts = cs_glob_time_step;
      if (ts->nt_cur < ts->nt_max + 1)
        cs_time_step_define_nt_max(ts->nt_cur + 1);

      cpl->iteration = -1;

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
    c_dt[i] = dttmp;

  cpl->dt = dttmp;

  int verbosity = (cs_log_default_is_active()) ? cpl->verbosity : 0;
  if (verbosity > 0)
    bft_printf("----------------------------------\n"
               "reference time step:     %4.21e\n"
               "code_saturne time step:  %4.2le\n"
               "code_aster time step:    %4.2le\n"
               "selected time step:      %4.2le \n"
               "----------------------------------\n\n",
               cpl->dtref, c_dt[0], dt_ast, cpl->dt);

  /* Reset sub-iteration count */
  cpl->s_it_id = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send stresses acting on the fluid/structure interface
 *        and receive displacements.
 *
 * \param[in]  fluid_forces  forces from fluid at coupled faces
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_exchange_fields(const cs_real_t  fluid_forces[])
{
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  int verbosity = (cs_log_default_is_active()) ? cpl->verbosity : 0;

  const cs_lnum_t n_faces = cpl->n_faces;

  for (cs_lnum_t i = 0; i < 3*n_faces; i++)
    cpl->foras[i] = fluid_forces[i];

  /* Send prediction
     (no difference between explicit and implicit cases for forces) */

  cs_real_t alpha = 2.0;
  cs_real_t c1    = alpha;
  cs_real_t c2    = 1-alpha;
  cs_real_t c3    = 0.;

  _pred(cpl->fopas,
        cpl->foras,
        cpl->foaas,
        cpl->foaas,
        c1,
        c2,
        c3,
        n_faces);

  if (verbosity > 0)
    bft_printf("--------------------------------------\n"
               "Forces prediction coefficients\n"
               " C1: %4.2le\n"
               " C2: %4.2le\n"
               " C3: %4.2le\n"
               "--------------------------------------\n\n",
               c1, c2, c3);

  /* Send forces */

#if defined(HAVE_PARAMEDMEM)

  cs_paramedmem_field_export_l(cpl->mc_faces, _name_f_f, cpl->fopas);
  cs_paramedmem_attach_field_by_name(cpl->mc_faces, _name_f_f);

  if (verbosity > 1) {
    bft_printf
      (_("code_aster: starting MEDCoupling send of values "
         "at coupled faces..."));
    bft_printf_flush();
  }

  cs_paramedmem_send_data(cpl->mc_faces);

  if (verbosity > 1) {
    bft_printf(_("[ok]\n"));
    bft_printf_flush();
  }

#endif /* defined(HAVE_PARAMEDMEM) */

  /* Second stage (TODO: place in another, better named function) */
  /* ------------------------------------------------------------ */

  /* explicit case: no need fo a convergence test */

  int icv = 1;

  if (cpl->nbssit <= 1) {

    /* handle convergence even when no test is done */
    cpl->icv1 = icv;
    _send_icv2(cpl, icv);

    /* receive displacements from code_aster */
    _recv_dyn(cpl);

    /* save previous values */
    _val_ant(cpl);

  }

  /* implicit case: requires a convergence test */

  else if (cpl->nbssit > 1) {

    /* compute icv */

    int ierr = _conv(cpl, &icv);
    cpl->icv1 = icv;
    icv = cpl->icv2;
    _send_icv2(cpl, icv);

    if ((cpl->s_it_id +1 >= cpl->nbssit) || (icv == 1)) {
      /* receive displacements  computed by code_aster */
      if (ierr >= 0) _recv_dyn(cpl);

      /* then use with code_saturne ? the question remains open... */
      /* if (ierr >= 0) _send2_dyn(); */

      /* receive displacements from code_aster */
      if (ierr >= 0)  _recv_dyn(cpl);
    }
    else {
      cpl->s_it_id += 1;
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute predicted or exact displacement of the
 *        fluid/structure interface.
 *
 * \param[out]  disp  prescribed displacement at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_compute_displacement(cs_real_t  disp[][3])
{
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  const cs_lnum_t  nb_dyn = cpl->n_vertices;

  /* Predict displacements */

  cs_real_t c1, c2, c3, alpha, beta;

  /* separate prediction for explicit/implicit cases */
  if (cpl->s_it_id == 0) {
    alpha = 0.5;
    beta  = 0.;
    c1    = 1.;
    c2    = (alpha + beta) * cs_glob_time_step->dt[0];
    c3    = -beta * cs_glob_time_step->dt[1];
    _pred(cpl->xastp,
          cpl->xast,
          cpl->xvast,
          cpl->xvasa,
          c1,
          c2,
          c3,
          nb_dyn);
  }
  else { /* if (cpl->s_it_id > 0) */
    alpha = 0.5;
    c1    = alpha;
    c2    = 1. - alpha;
    c3    = 0.;
    _pred(cpl->xastp,
          cpl->xast,
          cpl->xastp,
          cpl->xast,
          c1,
          c2,
          c3,
          nb_dyn);
  }

  int verbosity = (cs_log_default_is_active()) ? cpl->verbosity : 0;

  if (verbosity > 0) {

    bft_printf("*********************************\n"
               "*     sub - iteration %i        *\n"
               "*********************************\n\n",
               cpl->s_it_id);

    bft_printf("--------------------------------------------\n"
               "Displacement prediction coefficients\n"
               " C1: %4.2le\n"
               " C2: %4.2le\n"
               " C3: %4.2le\n"
               "--------------------------------------------\n\n",
               c1, c2, c3);

  }

  /* Set in disp the values of prescribed displacements */

#if defined(HAVE_PARAMEDMEM)

  const cs_lnum_t *vtx_ids
    = cs_paramedmem_mesh_get_vertex_list(cpl->mc_vertices);

  _scatter_values_r3(cpl->n_vertices,
                     vtx_ids,
                     (const cs_real_3_t *)cpl->xastp,
                     disp);

#endif /* defined(HAVE_PARAMEDMEM) */
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
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1)
    MPI_Bcast(&(cpl->icv1), 1, MPI_INT, 0, cs_glob_mpi_comm);

#endif

  return cpl->icv1;
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
  cs_ast_coupling_t  *cpl = cs_glob_ast_coupling;

  cpl->icv2 = icved;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
