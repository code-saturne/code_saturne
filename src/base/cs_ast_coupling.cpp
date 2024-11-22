/*============================================================================
 * code_aster coupling
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>
#include <ple_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_extract.h"

#include "cs_all_to_all.h"
#include "cs_array.h"
#include "cs_blas.h"
#include "cs_calcium.h"
#include "cs_coupling.h"
#include "cs_interface.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"
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
  int         root_rank; /* Application root rank in MPI_COMM_WORLD */
  const char *app_type;  /* Application type name (may be empty) */
  const char *app_name;  /* Application instance name (may be empty) */

} ple_coupling_mpi_set_info_t;

#endif

/* Main code_aster coupling structure */

struct _cs_ast_coupling_t {
  ple_coupling_mpi_set_info_t aci; /* code_aster coupling info */

  cs_lnum_t n_faces;    /* Local number of coupled faces */
  cs_lnum_t n_vertices; /* Local number of coupled vertices */

  cs_gnum_t n_g_faces;    /* Global number of coupled faces */
  cs_gnum_t n_g_vertices; /* Global number of coupld vertices */

  cs_paramedmem_coupling_t *mc_faces;
  cs_paramedmem_coupling_t *mc_vertices;

  int verbosity;     /* verbosity level */
  int visualization; /* visualization level */

  fvm_nodal_t *post_mesh;    /* Optional mesh for post-processing output */
  int          post_mesh_id; /* 0 if post-processing is not active,
                                or post-processing mesh_id (< 0) */

  int iteration; /* 0 for initialization, < 0 for disconnect,
                    iteration from (re)start otherwise */

  int nbssit; /* number of sub-iterations */

  cs_real_t dt;
  cs_real_t dtref;  /* reference time step */
  cs_real_t epsilo; /* scheme convergence threshold */

  int icv1; /* Convergence indicator */
  int icv2; /* Convergence indicator (final) */

  cs_real_t rcv1; /* Value of the residual */

  cs_real_t lref; /* Characteristic macroscopic domain length */

  int s_it_id; /* Sub-iteration id */

  cs_real_t *xast_curr[2]; /* Mesh displacement at iteration k+1, k, k-1*/
  cs_real_t *xsat_pred[2]; /* Predicted mesh at iteration k+1, k, k-1*/
  cs_real_t *vast_curr;    /* Mesh velocity last received (current iteration) */
  cs_real_t *vast_prev;    /* Mesh velocity at previous time step n-1 */
  cs_real_t *vast_pprev;   /* Mesh velocity at previous time step n-2 */

  cs_real_t *forc_curr; /* Fluid forces at current sub-iteration */
  cs_real_t *forc_prev; /* Fluid forces at previous time step */
  cs_real_t *forc_pred; /* Predicted fluid forces at current sub-iteration */

  cs_real_t aexxst; /*!< coefficient for the predicted displacement */
  cs_real_t bexxst; /*!< coefficient for the predicted displacement */
  cs_real_t rexxst; /*!< coefficient for the relaxation displacement */

  cs_real_t cfopre; /*!< coefficient for the predicted force */

  cs_real_t *tmp[3]; /* Temporary array */
};

/*============================================================================
 * Static global variables
 *============================================================================*/

static const char _name_f_f[] = "fluid_forces";
static const char _name_m_d[] = "mesh_displacement";
static const char _name_m_v[] = "mesh_velocity";

static int _verbosity     = 1;
static int _visualization = 1;

#define DEBUG_LEVEL_CPL 0

/*============================================================================
 * Global variables
 *============================================================================*/

cs_ast_coupling_t *cs_glob_ast_coupling = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Allocate and initialize dynamic vectors (cs_real_t) based on the 'nb_dyn'
 * number of points.
 *----------------------------------------------------------------------------*/

static int
_get_current_verbosity(const cs_ast_coupling_t *ast_cpl)
{
  return (cs_log_default_is_active()) ? ast_cpl->verbosity : 0;
}

static void
_allocate_arrays(cs_ast_coupling_t *ast_cpl)
{
  const cs_lnum_t nb_dyn = ast_cpl->n_vertices;
  const cs_lnum_t nb_for = ast_cpl->n_faces;

  for (int i = 0; i < 2; i++) {
    BFT_MALLOC(ast_cpl->xast_curr[i], 3 * nb_dyn, cs_real_t);
    BFT_MALLOC(ast_cpl->xsat_pred[i], 3 * nb_dyn, cs_real_t);
  }

  BFT_MALLOC(ast_cpl->vast_curr, 3 * nb_dyn, cs_real_t);
  BFT_MALLOC(ast_cpl->vast_prev, 3 * nb_dyn, cs_real_t);
  BFT_MALLOC(ast_cpl->vast_pprev, 3 * nb_dyn, cs_real_t);

  cs_arrays_set_value<cs_real_t, 1>(3 * nb_dyn,
                                    0.,
                                    ast_cpl->xast_curr[0],
                                    ast_cpl->xast_curr[1],
                                    ast_cpl->xsat_pred[0],
                                    ast_cpl->xsat_pred[1],
                                    ast_cpl->vast_curr,
                                    ast_cpl->vast_prev,
                                    ast_cpl->vast_pprev);

  BFT_MALLOC(ast_cpl->forc_curr, 3 * nb_for, cs_real_t);
  BFT_MALLOC(ast_cpl->forc_prev, 3 * nb_for, cs_real_t);
  BFT_MALLOC(ast_cpl->forc_pred, 3 * nb_for, cs_real_t);

  cs_arrays_set_value<cs_real_t, 1>(3 * nb_for,
                                    0.,
                                    ast_cpl->forc_curr,
                                    ast_cpl->forc_prev,
                                    ast_cpl->forc_pred);

  for (int i = 0; i < 3; i++) {
    BFT_MALLOC(ast_cpl->tmp[i], 3 * CS_MAX(nb_dyn, nb_for), cs_real_t);
  }
}

/*----------------------------------------------------------------------------
 * Scatter values of type cs_real_3_t (tuples) based on indirection list
 *
 * parameters:
 *   n_elts   <-- number of elements
 *   elt_ids  <-- element ids, or nullptr
 *   v_in     <-- input values, on elt_ids location
 *   v_out    <-> output values, on parent location
 *----------------------------------------------------------------------------*/

static void
_scatter_values_r3(cs_lnum_t         n_elts,
                   const cs_lnum_t   elt_ids[],
                   const cs_real_3_t v_in[],
                   cs_real_3_t       v_out[])
{
  assert(v_in != nullptr && v_out != nullptr);

  if (elt_ids != nullptr) {
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
 * Predict displacement or forces based on values of the current and
 * previous time step(s)
 *
 * valpre = c1 * val1 + c2 * val2 + c3 * val3
 *----------------------------------------------------------------------------*/

static void
_pred(cs_real_t       *valpre,
      const cs_real_t *val1,
      const cs_real_t *val2,
      const cs_real_t *val3,
      const cs_real_t  c1,
      const cs_real_t  c2,
      const cs_real_t  c3,
      const cs_lnum_t  n)
{
  if (n < 1)
    return;

  /* Update prediction array */
  const cs_lnum_t size = 3 * n;

  if (val3 != nullptr) {
#pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++) {
      valpre[i] = c1 * val1[i] + c2 * val2[i] + c3 * val3[i];
    }
  }
  else {
#pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++) {
      valpre[i] = c1 * val1[i] + c2 * val2[i];
    }
  }
}

/*----------------------------------------------------------------------------
 * Compute the L2 norm of the difference between vectors vect1 and vect2
 *
 * dinorm = sqrt(sum on nbpts i
 *                 (sum on component j
 *                    ((vect1[i,j]-vect2[i,j])^2)))/nbpts
 *----------------------------------------------------------------------------*/

static cs_real_t
_dinorm(cs_real_t *vect1, cs_real_t *vect2, cs_lnum_t nbpts)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  /* Compute the norm of the difference */
  cs_array_difference(3 * nbpts, vect1, vect2, cpl->tmp[0]);
  cs_real_t norm = cs_dot_xx(3 * nbpts, cpl->tmp[0]);

  /* Note that for vertices, vertices at shared parallel boundaries
     will appear multiple tiles, so have a higher "weight" than
     others, but the effect on the global norm should be minor,
     so we avoid a more complex test here */

  cs_real_t rescale = (cs_real_t)nbpts;

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    cs_real_t val[2] = { norm, rescale };
    cs_parall_sum(2, CS_DOUBLE, &val);

    norm = val[0], rescale = val[1];
  }
#endif

  norm = sqrt(norm / rescale);
  return norm;
}

/*----------------------------------------------------------------------------
 * Compute coefficient for aitken acceleration
 *----------------------------------------------------------------------------*/

static cs_real_t
_aitken(const cs_real_t *dp_k,
        const cs_real_t *dp_km,
        const cs_real_t *d_kp,
        const cs_real_t *d_k)
{
  cs_ast_coupling_t *cpl  = cs_glob_ast_coupling;
  const cs_lnum_t    size = 3 * cpl->n_vertices;

  cs_real_t *tmp1 = cpl->tmp[1], *tmp2 = cpl->tmp[2];

  /* Note that for vertices, vertices at shared parallel boundaries
   will appear multiple tiles, so have a higher "weight" than
   others, but the effect on the global direction should be minor,
   so we avoid a more complex test here */

  /* difference */
  cs_array_difference(size, dp_k, dp_km, cpl->tmp[0]);
  cs_array_difference(size, d_kp, dp_k, tmp1);
  cs_array_difference(size, d_k, dp_km, tmp2);

#pragma omp parallel for if (size > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < size; ii++)
    tmp1[ii] += tmp2[ii];

  cs_real_t xx, xy;
  cs_gdot_xx_xy(size, tmp1, cpl->tmp[0], &xx, &xy);

  return xy / xx;
}

/*----------------------------------------------------------------------------
 * Post process variables associated with code_aster couplings
 *
 * parameters:
 *   coupling        <--  Void pointer to code_aster coupling structure
 *   ts              <--  time step status structure, or nullptr
 *----------------------------------------------------------------------------*/

static void
_cs_ast_coupling_post_function(void *coupling, const cs_time_step_t *ts)
{
  auto cpl = static_cast<const cs_ast_coupling_t *>(coupling);

  if (cpl->post_mesh == nullptr)
    return;

  /* Note: since numbering in fvm_nodal_t structures (ordered by
     element type) may not align with the selection order, we need to project
     values on parent faces first */

  const cs_lnum_t *face_ids = cs_paramedmem_mesh_get_elt_list(cpl->mc_faces);
  const cs_lnum_t *vtx_ids =
    cs_paramedmem_mesh_get_vertex_list(cpl->mc_vertices);

  cs_real_t       *values;
  const cs_mesh_t *m      = cs_glob_mesh;
  cs_lnum_t        n_vals = CS_MAX(m->n_b_faces, m->n_vertices) * 3;
  BFT_MALLOC(values, n_vals, cs_real_t);
  cs_arrays_set_value<cs_real_t, 1>(n_vals, 0., values);

  _scatter_values_r3(cpl->n_vertices,
                     vtx_ids,
                     (const cs_real_3_t *)cpl->xast_curr[0],
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
                     (const cs_real_3_t *)cpl->vast_curr,
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
                     (const cs_real_3_t *)cpl->forc_curr,
                     (cs_real_3_t *)values);

  cs_post_write_var(cpl->post_mesh_id,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    "Stress",
                    3,
                    true, /* interlaced */
                    true, /* on parent */
                    CS_POST_TYPE_cs_real_t,
                    nullptr, /* cell values */
                    nullptr, /* interior face values */
                    values,
                    ts);

  BFT_FREE(values);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Query number of couplings with code_aster.
 *
 * Currently, a single coupling with code_aster is possible.
 */
/*----------------------------------------------------------------------------*/

int
cs_ast_coupling_n_couplings(void)
{
  int retval = 0;

  if (cs_glob_ast_coupling != nullptr)
    retval = 1;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initial exchange with code_aster.
 *
 * \param[in]  nalimx  maximum number of implicitation iterations of
 *                     the structure displacement
 * \param[in]  epalim  relative precision of implicitation of
 *                     the structure displacement
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_initialize(int nalimx, cs_real_t epalim)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  int       nbpdtm = ts->nt_max;
  cs_real_t ttinit = ts->t_prev;
  int       idtvar = cs_glob_time_step_options->idtvar;

  /* Allocate global coupling structure */

  cs_ast_coupling_t *cpl;

  BFT_MALLOC(cpl, 1, cs_ast_coupling_t);

  memset(&(cpl->aci), 0, sizeof(ple_coupling_mpi_set_info_t));
  cpl->aci.root_rank = -1;

  cpl->n_faces    = 0;
  cpl->n_vertices = 0;

  cpl->n_g_faces    = 0;
  cpl->n_g_vertices = 0;

#if defined(HAVE_PARAMEDMEM)

  cpl->mc_faces    = nullptr;
  cpl->mc_vertices = nullptr;

#endif

  cpl->verbosity     = _verbosity;
  cpl->visualization = _visualization;

  cpl->post_mesh = nullptr;

  cpl->iteration = 0; /* < 0 for disconnect */

  cpl->nbssit = nalimx; /* maximum number of sub-iterations */

  cpl->dt     = 0.;
  cpl->dtref  = ts->dt_ref; /* reference time step */
  cpl->epsilo = epalim;     /* scheme convergence threshold */

  cpl->aexxst = 1.0;
  cpl->bexxst = 0.5;
  cpl->rexxst = 0.5;
  cpl->cfopre = 2.0;

  cpl->icv1 = 0;
  cpl->icv2 = 0;
  cpl->lref = 0.;

  cpl->s_it_id = 0; /* Sub-iteration id */

  for (int i = 0; i < 2; i++) {
    cpl->xast_curr[i] = nullptr;
    cpl->xsat_pred[i] = nullptr;
  }
  cpl->vast_curr  = nullptr;
  cpl->vast_prev  = nullptr;
  cpl->vast_pprev = nullptr;

  cpl->forc_curr = nullptr;
  cpl->forc_prev = nullptr;
  cpl->forc_pred = nullptr;

  for (int i = 0; i < 3; i++) {
    cpl->tmp[i] = nullptr;
  }

  cs_glob_ast_coupling = cpl;

  cs_calcium_set_verbosity(cpl->verbosity);

  /* Find root rank of coupling */

#if defined(PLE_HAVE_MPI)
  const ple_coupling_mpi_set_t *mpi_apps = cs_coupling_get_mpi_apps();

  if (mpi_apps != nullptr) {
    int n_apps     = ple_coupling_mpi_set_n_apps(mpi_apps);
    int n_ast_apps = 0;

    /* First pass to count available code_aster couplings */

    for (int i = 0; i < n_apps; i++) {
      const ple_coupling_mpi_set_info_t ai =
        ple_coupling_mpi_set_get_info(mpi_apps, i);

      if (strncmp(ai.app_type, "code_aster", 10) == 0)
        n_ast_apps += 1;
    }

    /* In single-coupling mode, no identification necessary */

    if (n_ast_apps == 1) {
      for (int i = 0; i < n_apps; i++) {
        const ple_coupling_mpi_set_info_t ai =
          ple_coupling_mpi_set_get_info(mpi_apps, i);

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
      bft_error(__FILE__,
                __LINE__,
                0,
                "Detected %d code_aster instances; can handle exactly 1.",
                n_ast_apps);
  }
  else {
    bft_error(__FILE__, __LINE__, 0, "No PLE application detected.");
  }

#else

  bft_error(__FILE__,
            __LINE__,
            0,
            "code_aster coupling requires PLE with MPI support.");

#endif

  /* Calcium  (communication) initialization */

  if (cs_glob_rank_id <= 0) {
    int verbosity = _get_current_verbosity(cpl);

    if (verbosity > 0) {
      bft_printf("Send calculation parameters to code_aster\n");
    }

    /* Send data */

    cs_calcium_write_int(cpl->aci.root_rank, 0, "NBPDTM", 1, &nbpdtm);
    cs_calcium_write_int(cpl->aci.root_rank, 0, "NBSSIT", 1, &(cpl->nbssit));
    cs_calcium_write_int(cpl->aci.root_rank, 0, "TADAPT", 1, &idtvar);

    cs_calcium_write_double(cpl->aci.root_rank, 0, "EPSILO", 1, &(cpl->epsilo));
    cs_calcium_write_double(cpl->aci.root_rank, 0, "TTINIT", 1, &ttinit);
    cs_calcium_write_double(cpl->aci.root_rank, 0, "PDTREF", 1, &(cpl->dtref));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize coupling with code_aster.
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_finalize(void)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl == nullptr)
    return;

  for (int i = 0; i < 2; i++) {
    BFT_FREE(cpl->xast_curr[i]);
    BFT_FREE(cpl->xsat_pred[i]);
  }

  BFT_FREE(cpl->vast_curr);
  BFT_FREE(cpl->vast_prev);
  BFT_FREE(cpl->vast_pprev);

  BFT_FREE(cpl->forc_curr);
  BFT_FREE(cpl->forc_prev);
  BFT_FREE(cpl->forc_pred);

  for (int i = 0; i < 3; i++) {
    BFT_FREE(cpl->tmp[i]);
  }

  if (cpl->post_mesh != nullptr)
    cpl->post_mesh = fvm_nodal_destroy(cpl->post_mesh);

  cs_paramedmem_coupling_destroy(cpl->mc_vertices);
  cs_paramedmem_coupling_destroy(cpl->mc_faces);

  cpl->mc_vertices = nullptr;
  cpl->mc_faces    = nullptr;

  BFT_FREE(cpl);

  cs_glob_ast_coupling = cpl;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set coefficient for prediction
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_set_coefficients(cs_real_t aexxst,
                                 cs_real_t bexxst,
                                 cs_real_t cfopre)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl == nullptr)
    return;

  cpl->aexxst = aexxst;
  cpl->bexxst = bexxst;
  cpl->cfopre = cfopre;
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
cs_ast_coupling_geometry(cs_lnum_t        n_faces,
                         const cs_lnum_t *face_ids,
                         cs_real_t        almax)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl->aci.root_rank > -1) {
    cpl->mc_faces    = cs_paramedmem_coupling_create(nullptr,
                                                  cpl->aci.app_name,
                                                  "fsi_face_exchange");
    cpl->mc_vertices = cs_paramedmem_coupling_create(nullptr,
                                                     cpl->aci.app_name,
                                                     "fsi_vertices_exchange");
  }
  else {
    cpl->mc_faces =
      cs_paramedmem_coupling_create_uncoupled("fsi_face_exchange");
    cpl->mc_vertices =
      cs_paramedmem_coupling_create_uncoupled("fsi_vertices_exchange");
  }

  cs_paramedmem_add_mesh_from_ids(cpl->mc_faces, n_faces, face_ids, 2);

  cs_paramedmem_add_mesh_from_ids(cpl->mc_vertices, n_faces, face_ids, 2);

  cpl->n_faces    = n_faces;
  cpl->n_vertices = cs_paramedmem_mesh_get_n_vertices(cpl->mc_vertices);

  fvm_nodal_t *fsi_mesh =
    cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                   "FSI_mesh_1",
                                   true, /* include families */
                                   0,
                                   n_faces,
                                   nullptr,
                                   face_ids);

  cpl->n_g_faces = n_faces;
  cs_parall_counter(&(cpl->n_g_faces), 1);
  cpl->n_g_vertices = fvm_nodal_get_n_g_vertices(fsi_mesh);

  /* Creation of the information structure for code_saturne/code_aster
     coupling */

  static_assert(sizeof(cs_coord_t) == sizeof(cs_real_t),
                "Incorrect size of cs_coord_t");

  if (cpl->visualization > 0) {
    cpl->post_mesh = fsi_mesh;
  }
  else {
    fsi_mesh = fvm_nodal_destroy(fsi_mesh);
  }

  _allocate_arrays(cpl);

  if (almax <= 0) {
    bft_error(__FILE__,
              __LINE__,
              0,
              "%s: almax = %g, where a positive value is expected.",
              __func__,
              almax);
  }

  cpl->lref = almax;

  if (cs_glob_rank_id <= 0) {
    int verbosity = _get_current_verbosity(cpl);

    if (verbosity > 0) {
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
  }

  /* Define coupled fields */

  cs_paramedmem_def_coupled_field(cpl->mc_vertices,
                                  _name_m_d,
                                  3,
                                  CS_MEDCPL_FIELD_INT_MAXIMUM,
                                  CS_MEDCPL_ON_NODES,
                                  CS_MEDCPL_ONE_TIME);

  cs_paramedmem_def_coupled_field(cpl->mc_vertices,
                                  _name_m_v,
                                  3,
                                  CS_MEDCPL_FIELD_INT_MAXIMUM,
                                  CS_MEDCPL_ON_NODES,
                                  CS_MEDCPL_ONE_TIME);

  cs_paramedmem_def_coupled_field(cpl->mc_faces,
                                  _name_f_f,
                                  3,
                                  CS_MEDCPL_FIELD_INT_CONSERVATION,
                                  CS_MEDCPL_ON_CELLS,
                                  CS_MEDCPL_ONE_TIME);

  /* Post-processing */

  if (cpl->visualization > 0) {
    const int writer_ids[] = { CS_POST_WRITER_DEFAULT };
    cpl->post_mesh_id      = cs_post_get_free_mesh_id();

    cs_post_define_existing_mesh(cpl->post_mesh_id,
                                 cpl->post_mesh,
                                 0,     /* dim_shift */
                                 false, /* transfer */
                                 false, /* auto variables */
                                 1,     /* number of associated writers */
                                 writer_ids);

    cs_post_add_time_dep_output(_cs_ast_coupling_post_function, (void *)cpl);
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
cs_ast_coupling_exchange_time_step(cs_real_t c_dt[])
{
  if (DEBUG_LEVEL_CPL > 0) {
    bft_printf(_("cs_ast_coupling_exchange_time_step function\n"));
    bft_printf_flush();
  }

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  cs_real_t dttmp  = cpl->dtref;
  cs_real_t dt_ast = dttmp;

  if (cpl->iteration < 0)
    return;

  int err_code = 0;

  cpl->iteration += 1;

  if (cs_glob_rank_id <= 0) {
    cs_real_t dt_sat     = c_dt[0];
    int       n_val_read = 0;

    /* Receive time step sent by code_aster */

    err_code = cs_calcium_read_double(cpl->aci.root_rank,
                                      &(cpl->iteration),
                                      "DTAST",
                                      1,
                                      &n_val_read,
                                      &dt_ast);

    if (err_code >= 0) {
      assert(n_val_read == 1);

      /* Choose smallest time step */

      if (dt_ast < dttmp)
        dttmp = dt_ast;
      if (dt_sat < dttmp)
        dttmp = dt_sat;

      err_code = cs_calcium_write_double(cpl->aci.root_rank,
                                         cpl->iteration,
                                         "DTCALC",
                                         1,
                                         &dttmp);
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

  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  cs_arrays_set_value<cs_real_t, 1>(n_cells_ext, dttmp, c_dt);

  cpl->dt = dttmp;

  int verbosity = _get_current_verbosity(cpl);
  if (verbosity > 0)
    bft_printf("----------------------------------\n"
               "reference time step:     %4.2le\n"
               "code_saturne time step:  %4.2le\n"
               "code_aster time step:    %4.2le\n"
               "selected time step:      %4.2le \n"
               "----------------------------------\n\n",
               cpl->dtref,
               c_dt[0],
               dt_ast,
               cpl->dt);

  /* Reset sub-iteration count */
  cpl->s_it_id = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to array of fluid forces at faces coupled with
 *        code_aster.
 *
 * \return  array of forces from fluid at coupled faces
 */
/*----------------------------------------------------------------------------*/

cs_real_3_t *
cs_ast_coupling_get_fluid_forces_pointer(void)
{
  cs_real_3_t *f_forces = nullptr;

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl != nullptr)
    f_forces = (cs_real_3_t *)cpl->forc_curr;

  return f_forces;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send stresses acting on the fluid/structure interface.
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_send_fluid_forces(void)
{
  if (DEBUG_LEVEL_CPL > 0) {
    bft_printf(_("cs_ast_coupling_send_fluid_forces function\n"));
    bft_printf_flush();
  }

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  int verbosity = _get_current_verbosity(cpl);

  const cs_lnum_t n_faces = cpl->n_faces;

  /* Prediction ared defined in Fabien Huvelin PhD*/
  cs_real_t c1, c2;
  if (cpl->s_it_id == 0) {
    /* Explicit synchrone prediction */
    c1 = cpl->cfopre;
    c2 = 1.0 - cpl->cfopre;
    _pred(cpl->forc_pred,
          cpl->forc_curr,
          cpl->forc_prev,
          nullptr,
          c1,
          c2,
          0.0,
          n_faces);
  }
  else {
    /* Implicit prediction */
    c1 = 1.0;
    c2 = 0.0;
    cs_array_copy(3 * n_faces, cpl->forc_curr, cpl->forc_pred);
  }

  if (verbosity > 0)
    bft_printf("--------------------------------------\n"
               "Forces prediction coefficients\n"
               " C1: %4.2le\n"
               " C2: %4.2le\n"
               "--------------------------------------\n\n",
               c1,
               c2);

  /* Send forces */

  if (verbosity > 1) {
    bft_printf(_("code_aster: starting MEDCoupling send of values "
                 "at coupled faces..."));
    bft_printf_flush();
  }

  cs_paramedmem_send_field_vals_l(cpl->mc_faces, _name_f_f, cpl->forc_pred);

  if (verbosity > 1) {
    bft_printf(_("[ok]\n"));
    bft_printf_flush();
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Evaluate convergence of the coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_evaluate_cvg(void)
{
  if (DEBUG_LEVEL_CPL > 0) {
    bft_printf(_("cs_ast_coupling_evaluate_cvg function\n"));
    bft_printf_flush();
  }

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  int verbosity = _get_current_verbosity(cpl);

  int icv = 1;

  if (cpl->nbssit > 1) {
    /* implicit case: requires a convergence test */

    /* compute icv */

    cpl->rcv1 = _dinorm(cpl->xast_curr[0], cpl->xsat_pred[0], cpl->n_vertices) /
                cpl->lref;

    if (verbosity > 0)
      bft_printf("--------------------------------\n"
                 "convergence test:\n"
                 "residual = %4.2le\n",
                 cpl->rcv1);

    if (cpl->rcv1 <= cpl->epsilo) {
      icv = 1;

      if (verbosity > 0)
        bft_printf("icv = %d\n"
                   "convergence of sub iteration\n"
                   "----------------------------\n",
                   icv);
    }
    else {
      icv = 0;
      if (verbosity > 0)
        bft_printf("icv = %i\n"
                   "non convergence of sub iteration\n"
                   "--------------------------------\n",
                   icv);
    }
  }

  cpl->icv1 = icv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Receive displacements of the fluid/structure interface.
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_recv_displacement(void)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  int verbosity = _get_current_verbosity(cpl);

  if (verbosity > 1) {
    bft_printf(_("code_aster: starting MEDCouping receive of values "
                 "at coupled vertices..."));
    bft_printf_flush();
  }

  /* Update values */
  const cs_lnum_t nb_dyn = cpl->n_vertices;

  cs_array_copy(3 * nb_dyn, cpl->xast_curr[1], cpl->xast_curr[2]);
  cs_array_copy(3 * nb_dyn, cpl->xast_curr[0], cpl->xast_curr[1]);

  /* Received discplacement and velocity field */
  cs_paramedmem_recv_field_vals_l(cpl->mc_vertices,
                                  _name_m_d,
                                  cpl->xast_curr[0]);
  cs_paramedmem_recv_field_vals_l(cpl->mc_vertices, _name_m_v, cpl->vast_curr);

  if (verbosity > 1) {
    bft_printf(_("[ok]\n"));
    bft_printf_flush();
  }

  /* For dry run, reset values to zero to avoid uninitialized values */
  if (cpl->aci.root_rank < 0) {
    cs_arrays_set_value<cs_real_t, 1>(3 * nb_dyn,
                                      0.,
                                      cpl->xast_curr[0],
                                      cpl->vast_curr);
  }

  cpl->s_it_id += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Save previous values
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_save_values(void)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  const cs_lnum_t nb_dyn = cpl->n_vertices;
  const cs_lnum_t nb_for = cpl->n_faces;

  /* record efforts */
  cs_array_copy(3 * nb_for, cpl->forc_pred, cpl->forc_prev);

  /* record dynamic data */
  cs_array_copy(3 * nb_dyn, cpl->vast_prev, cpl->vast_pprev);
  cs_array_copy(3 * nb_dyn, cpl->vast_curr, cpl->vast_prev);

  cs_arrays_set_value<cs_real_t, 1>(3 * nb_dyn,
                                    0.,
                                    cpl->xast_curr[1],
                                    cpl->xsat_pred[0],
                                    cpl->xsat_pred[1]);
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
cs_ast_coupling_compute_displacement(cs_real_t disp[][3])
{
  if (DEBUG_LEVEL_CPL > 0) {
    bft_printf(_("cs_ast_coupling_compute_displacement function\n"));
    bft_printf_flush();
  }

  assert(disp != nullptr);

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  if (cpl->iteration < 0)
    return;

  const cs_lnum_t nb_dyn = cpl->n_vertices;

  /* Predict displacements */

  cs_real_t c1, c2, c3;

  /* Prediction ared defined in Fabien Huvelin PhD*/

  /* separate prediction for explicit/implicit cases */
  if (cpl->s_it_id == 0) {
    /* Adams-Bashforth scheme of order 2 if aexxst = 1, bexxst = 0.5 */
    /* Euler explicit scheme of order 1 if aexxst = 1, bexxst = 0 */
    const cs_real_t dt_curr = cs_glob_time_step->dt[0];
    const cs_real_t dt_prev = cs_glob_time_step->dt[1];

    c1 = 1.;
    c2 = dt_curr * (cpl->aexxst + cpl->bexxst * dt_curr / dt_prev);
    c3 = -cpl->bexxst * dt_curr * dt_curr / dt_prev;
    _pred(cpl->xsat_pred[0],
          cpl->xast_curr[0],
          cpl->vast_prev,
          cpl->vast_pprev,
          c1,
          c2,
          c3,
          nb_dyn);
  }
  else {
    /* rexxst could be defined differently to have a better convergence */
    cs_real_t rexxst;
    if (cpl->rexxst < 0.) {
      rexxst = _aitken(cpl->xsat_pred[0],
                       cpl->xsat_pred[1],
                       cpl->xast_curr[0],
                       cpl->xast_curr[1]);
    }
    else {
      rexxst = cpl->rexxst;
    }
    c1 = rexxst;
    c2 = 1. - rexxst;
    c3 = 0.;

    cs_array_copy(3 * nb_dyn, cpl->xsat_pred[0], cpl->xsat_pred[1]);

    _pred(cpl->xsat_pred[0],
          cpl->xast_curr[0],
          cpl->xsat_pred[1],
          nullptr,
          c1,
          c2,
          c3,
          nb_dyn);
  }

  int verbosity = _get_current_verbosity(cpl);

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
               c1,
               c2,
               c3);
  }

  /* Set in disp the values of prescribed displacements */

  const cs_lnum_t *vtx_ids =
    cs_paramedmem_mesh_get_vertex_list(cpl->mc_vertices);

  _scatter_values_r3(cpl->n_vertices,
                     vtx_ids,
                     (const cs_real_3_t *)cpl->xsat_pred[0],
                     disp);
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
cs_ast_coupling_get_current_cvg(void)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  return cpl->icv1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Receive convergence value of code_saturne/code_aster coupling
 *
 * \return  residual of the convergence
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ast_coupling_get_current_residual(void)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  return cpl->rcv1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Send global convergence value of FSI calculations
 *
 * \param[in]  icved  convergence indicator (1: converged, 0: not converged)
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_set_final_cvg(int icved)
{
  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;

  cpl->icv2 = icved;

  if (cs_glob_rank_id > 0)
    return;

  cs_calcium_write_int(cpl->aci.root_rank,
                       cpl->iteration,
                       "ICVAST",
                       1,
                       &cpl->icv2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get verbosity level for code_aster coupling.
 *
 * \return  verbosity level for code_aster coupling
 */
/*----------------------------------------------------------------------------*/

int
cs_ast_coupling_get_verbosity(void)
{
  return _verbosity;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set verbosity level for code_aster coupling.
 *
 * \param[in]  verbosity      verbosity level for code_aster coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_set_verbosity(int verbosity)
{
  _verbosity = verbosity;

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;
  if (cpl != nullptr) {
    cpl->verbosity = verbosity;
    cs_calcium_set_verbosity(verbosity);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get visualization level for code_aster coupling.
 *
 * \return  visualization level for code_aster coupling
 */
/*----------------------------------------------------------------------------*/

int
cs_ast_coupling_get_visualization(void)
{
  return _visualization;
}
/*----------------------------------------------------------------------------*/
/*!
 * \brief Set visualization level for code_aster coupling.
 *
 * \param[in]  visualization  visualization level for code_aster coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_ast_coupling_set_visualization(int visualization)
{
  _visualization = visualization;

  cs_ast_coupling_t *cpl = cs_glob_ast_coupling;
  if (cpl != nullptr) {
    cpl->visualization = visualization;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
