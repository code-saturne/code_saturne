/*============================================================================
 * Routines to handle common features for building algebraic system in CDO
 * schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_local.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_hho_scaleq.h"
#include "cs_log.h"
#include "cs_xdef_eval.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_common.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_EQUATION_COMMON_DBG  0

/* Main categories to consider for high-level matrix structures */
#define CS_EQ_COMMON_VERTEX   0
#define CS_EQ_COMMON_FACE     1
#define CS_EQ_N_COMMONS       2

/*============================================================================
 * Local private variables
 *============================================================================*/

/* Temporary buffers useful during the building of all algebraic systems */
static size_t  cs_equation_common_work_buffer_size = 0;
static cs_real_t  *cs_equation_common_work_buffer = NULL;

/* Store the matrix structure and its assembler structures for each family
   of space discretizations */
static cs_matrix_assembler_t  **cs_equation_common_ma = NULL;
static cs_matrix_structure_t  **cs_equation_common_ms = NULL;

/* Structure related to the index of a matrix for vertex-based schemes
   vertex --> vertices through cell connectivity */
static cs_adjacency_t  *cs_connect_v2v = NULL;

/* Structure related to the index of a matrix for face-based schemes
   face --> faces through cell connectivity */
static cs_adjacency_t  *cs_connect_f2f = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/* Monitoring */
static cs_timer_counter_t  tca; // assembling process
static cs_timer_counter_t  tcc; // connectivity building

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a vertex -> vertices connectivity index
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_get_v2v(const cs_cdo_connect_t     *connect)
{
  /* Build a (sorted) v2v connectivity index */
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_adjacency_t  *c2v = connect->c2v;

  cs_adjacency_t  *v2c = cs_adjacency_transpose(n_vertices, c2v);
  cs_adjacency_t  *v2v = cs_adjacency_compose(n_vertices, v2c, c2v);

  cs_adjacency_sort(v2v);

  /* Update index (v2v has a diagonal entry. We remove it since we have in
     mind an matrix structure stored using the MSR format (with diagonal terms
     counted outside the index) */
  cs_lnum_t  shift = 0;
  cs_lnum_t  prev_start = v2v->idx[0];
  cs_lnum_t  prev_end = v2v->idx[1];

  for (cs_lnum_t i = 0; i < n_vertices; i++) {

    for (cs_lnum_t j = prev_start; j < prev_end; j++)
      if (v2v->ids[j] != i)
        v2v->ids[shift++] = v2v->ids[j];

    if (i != n_vertices - 1) { // Update prev_start and prev_end
      prev_start = v2v->idx[i+1];
      prev_end = v2v->idx[i+2];
    }
    v2v->idx[i+1] = shift;

  } // Loop on vertices

  BFT_REALLOC(v2v->ids, v2v->idx[n_vertices], cs_lnum_t);

  /* Free temporary buffers */
  cs_adjacency_free(&v2c);

  return v2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a face -> faces connectivity index
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_adjacency_t *
_get_f2f(const cs_cdo_connect_t     *connect)
{
  cs_adjacency_t  *f2f = NULL;

  const cs_lnum_t  n_faces = connect->n_faces[0];

  /* Build a face -> face connectivity */
  f2f = cs_adjacency_compose(n_faces, connect->f2c, connect->c2f);
  cs_adjacency_sort(f2f);

  /* Update index (v2v has a diagonal entry. We remove it since we have in
     mind an index structure for a  matrix stored using the MSR format */
  cs_lnum_t  shift = 0;
  cs_lnum_t  prev_start = f2f->idx[0];
  cs_lnum_t  prev_end = f2f->idx[1];

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    for (cs_lnum_t j = prev_start; j < prev_end; j++)
      if (f2f->ids[j] != i)
        f2f->ids[shift++] = f2f->ids[j];

    if (i != n_faces - 1) { // Update prev_start and prev_end
      prev_start = f2f->idx[i+1];
      prev_end = f2f->idx[i+2];
    }
    f2f->idx[i+1] = shift;

  } // Loop on faces

  BFT_REALLOC(f2f->ids, f2f->idx[n_faces], cs_lnum_t);

  return f2f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_matrix_assembler_t structure
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      x2x      pointer to a cs_adjacency_t struct.
 * \param[in]      rs       pointer to a range set or NULL if sequential
 * \param[in, out] ma       pointer to the cs_matrix_assembler_t to update
 */
/*----------------------------------------------------------------------------*/

static void
_build_matrix_assembler(cs_lnum_t                n_elts,
                        const cs_adjacency_t    *x2x,
                        const cs_range_set_t    *rs,
                        cs_matrix_assembler_t   *ma)
{
  cs_gnum_t  *grows = NULL, *gcols = NULL;

  /* First loop to count max size of the buffer */
  cs_lnum_t  max_size = 0;
  for (cs_lnum_t id = 0; id < n_elts; id++)
    max_size = CS_MAX(max_size, x2x->idx[id+1] - x2x->idx[id]);
  BFT_MALLOC(grows, max_size + 1, cs_gnum_t); // +1 for the diagonal entry
  BFT_MALLOC(gcols, max_size + 1, cs_gnum_t);

  for (cs_lnum_t id = 0; id < n_elts; id++) {

    cs_lnum_t  start = x2x->idx[id];
    cs_lnum_t  end = x2x->idx[id+1];
    cs_gnum_t  grow_id = rs->g_id[id];

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    grows[0] = grow_id, gcols[0] = grow_id;
    for (cs_lnum_t j = start, i = 1; j < end; j++, i++) {
      grows[i] = grow_id;
      gcols[i] = rs->g_id[x2x->ids[j]];
    }

    cs_matrix_assembler_add_g_ids(ma, end-start+1, grows, gcols);

  } // Loop on entities

  /* Now compute structure */
  cs_matrix_assembler_compute(ma);

  /* Free temporary buffers */
  BFT_FREE(grows);
  BFT_FREE(gcols);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 * \param[in]  quant         pointer to additional mesh quantities struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_allocate_common_structures(const cs_cdo_connect_t     *connect,
                                       const cs_cdo_quantities_t  *quant,
                                       const cs_time_step_t       *time_step,
                                       cs_flag_t                   scheme_flag)
{
  assert(connect != NULL); // Sanity check

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(tca); // assembling system
  CS_TIMER_COUNTER_INIT(tcc); // connectivity

  /* Two types of mat. ass. are considered:
     - The one related to matrix based on vertices
     - The one related to matrix based on faces
  */
  BFT_MALLOC(cs_equation_common_ma, CS_EQ_N_COMMONS, cs_matrix_assembler_t *);
  for (int i = 0; i < CS_EQ_N_COMMONS; i++)
    cs_equation_common_ma[i] = NULL;

  BFT_MALLOC(cs_equation_common_ms, CS_EQ_N_COMMONS, cs_matrix_structure_t *);
  for (int i = 0; i < CS_EQ_N_COMMONS; i++)
    cs_equation_common_ms[i] = NULL;

  /* Allocate cell-wise and face-wise view of a mesh */
  cs_cdo_local_initialize(connect);

  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_vertices = connect->n_vertices;

  /* Allocate and initialize matrix assembler and matrix structures */
  if (scheme_flag & CS_SCHEME_FLAG_CDOVB ||
      scheme_flag & CS_SCHEME_FLAG_CDOVCB) {

    cs_timer_t t0 = cs_timer_time();

    /* Build the "v2v" connectivity index */
    cs_connect_v2v = _get_v2v(connect);

    /* Monitoring */
    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(&tcc, &t0, &t1);

    cs_matrix_assembler_t  *ma =
      cs_matrix_assembler_create(connect->v_rs->l_range, true); // sep_diag

    _build_matrix_assembler(n_vertices, cs_connect_v2v, connect->v_rs, ma);

    cs_matrix_structure_t  *ms =
      cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

    /* Monitoring */
    cs_timer_t t2 = cs_timer_time();
    cs_timer_counter_add_diff(&tca, &t1, &t2);

    cs_equation_common_ma[CS_EQ_COMMON_VERTEX] = ma;
    cs_equation_common_ms[CS_EQ_COMMON_VERTEX] = ms;

  } // Vertex-based schemes and related ones

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB || scheme_flag & CS_SCHEME_FLAG_HHO) {

    cs_timer_t t0 = cs_timer_time();

    /* Build the "f2f" connectivity index */
    cs_connect_f2f = _get_f2f(connect);

    /* Monitoring */
    cs_timer_t t1 = cs_timer_time();
    cs_timer_counter_add_diff(&tcc, &t0, &t1);

    cs_matrix_assembler_t  *ma =
      cs_matrix_assembler_create(connect->f_rs->l_range, true); // sep_diag

    _build_matrix_assembler(n_faces, cs_connect_f2f, connect->f_rs, ma);

    cs_matrix_structure_t  *ms =
      cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

    /* Monitoring */
    cs_timer_t t2 = cs_timer_time();
    cs_timer_counter_add_diff(&tca, &t1, &t2);

    cs_equation_common_ma[CS_EQ_COMMON_FACE] = ma;
    cs_equation_common_ms[CS_EQ_COMMON_FACE] = ms;

  } // Face-based schemes and related ones

  /* Allocate shared buffer and initialize shared structures */
  size_t  cwb_size = 2*n_cells; // initial cell-wise buffer size

  if (scheme_flag & CS_SCHEME_FLAG_CDOVB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_vertices);

    /* Initialize additional structures */
    cs_cdovb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovb_scaleq_initialize();

  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOVCB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)2*(n_vertices + n_cells));

    /* Initialize additional structures */
    cs_cdovcb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovcb_scaleq_initialize();
  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    /* Initialize additional structures */
    cs_cdofb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdofb_scaleq_initialize();

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_faces);

  }

  if (scheme_flag & CS_SCHEME_FLAG_HHO &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    /* Initialize additional structures */
    cs_hho_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_hho_scaleq_initialize();

    // TODO: Update this value accordingly to HHO needs
    cwb_size = CS_MAX(cwb_size, (size_t)3*n_faces);

  }

  /* Assign static const pointers: shared pointers with a cs_domain_t */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Common buffer for temporary usage */
  cs_equation_common_work_buffer_size = cwb_size;
  BFT_MALLOC(cs_equation_common_work_buffer, cwb_size, double);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_common_structures(cs_flag_t   scheme_flag)
{
  /* Free cell-wise and face-wise view of a mesh */
  cs_cdo_local_finalize();

  cs_timer_t t0 = cs_timer_time();

  if (scheme_flag & CS_SCHEME_FLAG_CDOVB || scheme_flag & CS_SCHEME_FLAG_CDOVCB)
    cs_adjacency_free(&(cs_connect_v2v));

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB || scheme_flag & CS_SCHEME_FLAG_HHO)
    cs_adjacency_free(&(cs_connect_f2f));

    /* Free common structures specific to a numerical scheme */
  if ((scheme_flag & CS_SCHEME_FLAG_CDOVB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOVCB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovcb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOFB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdofb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_HHO) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_hho_scaleq_finalize();

  BFT_FREE(cs_equation_common_work_buffer);

  /* Monitoring */
  cs_timer_t t1 = cs_timer_time();
  cs_timer_counter_add_diff(&tcc, &t0, &t1);

  /* matrix assemblers and structures */
  for (int i = 0; i < CS_EQ_N_COMMONS; i++) {
    cs_matrix_structure_destroy(&(cs_equation_common_ms[i]));
    cs_matrix_assembler_destroy(&(cs_equation_common_ma[i]));
  }
  BFT_FREE(cs_equation_common_ms);
  BFT_FREE(cs_equation_common_ma);

  /* Monitoring */
  cs_timer_t t2 = cs_timer_time();
  cs_timer_counter_add_diff(&tca, &t1, &t2);
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %10s %10s\n",
                " ", "Connectivity", "Assembly");
  cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f %9.3f seconds\n",
                "<CDO/CommonEq> Runtime",
                tcc.wall_nsec*1e-9, tca.wall_nsec*1e-9);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are scalar-valued
 *          and attached to vertices
 *
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      eqp          pointer to a cs_equation_param_t
 * \param[in]      dir          pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 *
 * \return a pointer to a new allocated array storing the dirichlet values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_compute_dirichlet_sv(const cs_mesh_t            *mesh,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_list_t     *dir,
                                 cs_cell_builder_t          *cb)
{
  cs_flag_t  *flag = NULL, *counter = NULL;
  cs_real_t  *dir_val = NULL;

  const cs_lnum_t  *face_vtx_idx = mesh->b_face_vtx_idx;
  const cs_lnum_t  *face_vtx_lst = mesh->b_face_vtx_lst;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Initialization */
  BFT_MALLOC(dir_val, quant->n_vertices, cs_real_t);
  BFT_MALLOC(counter, quant->n_vertices, cs_flag_t);
  BFT_MALLOC(flag, quant->n_vertices, cs_flag_t);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

    dir_val[v_id] = 0;
    flag[v_id] = 0;    // No flag by default
    counter[v_id] = 0; // Number of faces with a Dir. related to a vertex

  }

  /* Define array storing the Dirichlet values */
  for (cs_lnum_t i = 0; i < dir->n_nhmg_elts; i++) {

    const cs_lnum_t  f_id = dir->elt_ids[i];
    const cs_lnum_t  *f2v_idx = face_vtx_idx + f_id;
    const int  n_vf = f2v_idx[1] - f2v_idx[0];
    const cs_lnum_t  *f2v_lst = face_vtx_lst + f2v_idx[0];

    const short int  def_id = dir->def_ids[i];
    const cs_xdef_t  *def = eqp->bc_desc[def_id];

    switch(def->type) {

    case CS_XDEF_BY_VALUE:
      {
        const cs_real_t  *constant_val = (cs_real_t *)def->input;

        for (short int v = 0; v < n_vf; v++) {

          const cs_lnum_t  v_id = f2v_lst[v];
          dir_val[v_id] += constant_val[0];
          flag[v_id] |= CS_CDO_BC_DIRICHLET;
          counter[v_id] += 1;

        }

      }
      break;

    case CS_XDEF_BY_ARRAY:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_array(n_vf,
                                          f2v_lst,
                                          true, // compact ouput
                                          cs_glob_mesh,
                                          cs_shared_connect,
                                          cs_shared_quant,
                                          cs_shared_time_step,
                                          def->input,
                                          eval);

        for (short int v = 0; v < n_vf; v++) {

          const cs_lnum_t  v_id = f2v_lst[v];
          dir_val[v_id] += eval[v];
          flag[v_id] |= CS_CDO_BC_DIRICHLET;
          counter[v_id] += 1;

        }

      }
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        cs_real_t  *eval = cb->values;

        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_vertices_by_analytic(n_vf,
                                             f2v_lst,
                                             true, // compact output
                                             cs_glob_mesh,
                                             cs_shared_connect,
                                             cs_shared_quant,
                                             cs_shared_time_step,
                                             def->input,
                                             eval);

        for (short int v = 0; v < n_vf; v++) {

          const cs_lnum_t  v_id = f2v_lst[v];
          dir_val[v_id] += eval[v];
          flag[v_id] |= CS_CDO_BC_DIRICHLET;
          counter[v_id] += 1;

        }

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition.\n"
                  " Stop computing the Dirichlet value.\n"));

    } // switch def_type

  } // Loop on faces with a non-homogeneous Dirichlet BC

  /* Define array storing the Dirichlet values */
  for (cs_lnum_t i = dir->n_nhmg_elts; i < dir->n_elts; i++) {

    const cs_lnum_t  f_id = dir->elt_ids[i];
    const cs_lnum_t  *f2v_idx = face_vtx_idx + f_id;
    const int  n_vf = f2v_idx[1] - f2v_idx[0];
    const cs_lnum_t  *f2v_lst = face_vtx_lst + f2v_idx[0];

    for (short int v = 0; v < n_vf; v++)
      flag[f2v_lst[v]] |= CS_CDO_BC_HMG_DIRICHLET;

  } // Loop on faces with a non-homogeneous Dirichlet BC

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_interface_set_max(cs_shared_connect->v_rs->ifs,
                         quant->n_vertices,
                         1,            // stride
                         false,        // interlace (not useful here)
                         CS_FLAG_TYPE, // unsigned short int
                         flag);

    cs_interface_set_sum(cs_shared_connect->v_rs->ifs,
                         quant->n_vertices,
                         1,            // stride
                         false,        // interlace (not useful here)
                         CS_FLAG_TYPE, // unsigned short int
                         counter);

    cs_interface_set_sum(cs_shared_connect->v_rs->ifs,
                         quant->n_vertices,
                         1,            // stride
                         false,        // interlace (not useful here)
                         CS_REAL_TYPE,
                         dir_val);

  }

  /* Homogeneous Dirichlet are always enforced (even in case of multiple BCs).
     If multiple Dirichlet BCs are set, a weighted sum is used to set the
     Dirichlet value at each corresponding vertex */
# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

    if (flag[v_id] & CS_CDO_BC_HMG_DIRICHLET)
      dir_val[v_id] = 0.;
    else if (flag[v_id] & CS_CDO_BC_DIRICHLET) {
      assert(counter[v_id] > 0);
      if (counter[v_id] > 1)
        dir_val[v_id] /= counter[v_id];
    }

  } // Loop on vertices

  /* Free temporary buffers */
  BFT_FREE(counter);
  BFT_FREE(flag);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_COMMON_DBG > 1
  cs_dump_array_to_listing("DIRICHLET_VALUES", quant->n_vertices, dir_val, 8);
#endif

  return dir_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs when DoFs are scalar-valued
 *          and attached to faces
 *
 * \param[in]      mesh      pointer to a cs_mesh_t structure
 * \param[in]      eqp       pointer to a cs_equation_param_t
 * \param[in]      dir       pointer to a cs_cdo_bc_list_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 *
 * \return a pointer to a new allocated array storing the dirichlet values
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_compute_dirichlet_sf(const cs_mesh_t            *mesh,
                                 const cs_equation_param_t  *eqp,
                                 const cs_cdo_bc_list_t     *dir,
                                 cs_cell_builder_t          *cb)
{
  CS_UNUSED(mesh);
  CS_UNUSED(dir);
  CS_UNUSED(cb);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_val = NULL;

  /* Initialization */
  BFT_MALLOC(dir_val, quant->n_b_faces, cs_real_t);

# pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_b_faces; f_id++) dir_val[f_id] = 0;

  /* Define array storing the Dirichlet values */
  for (int def_id = 0; def_id < eqp->n_bc_desc; def_id++) {

    const cs_xdef_t  *def = eqp->bc_desc[def_id];
    if (def->meta & CS_CDO_BC_DIRICHLET) {

      assert(def->dim == 1); // scalar variable
      const cs_boundary_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);
      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

#         pragma omp parallel for if (bz->n_faces > CS_THR_MIN)
          for (cs_lnum_t i = 0; i < bz->n_faces; i++)
            dir_val[bz->face_ids[i]] = constant_val[0];
        }
        break;

      case CS_XDEF_BY_ARRAY:
        {
          cs_xdef_array_input_t  *array_input =
            (cs_xdef_array_input_t *)def->input;

          assert(eqp->n_bc_desc == 1); // Only one definition allowed
          assert(array_input->stride == 1); // other cases not managed up to now
          assert(cs_test_flag(array_input->loc, cs_cdo_primal_face));

          if (bz->n_faces > 0) // Not only interior
            memcpy(dir_val, array_input->values, sizeof(cs_real_t)*bz->n_faces);

        }
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        /* Evaluate the boundary condition at each boundary vertex */
        cs_xdef_eval_at_b_faces_by_analytic(bz->n_faces,
                                            bz->face_ids,
                                            false, // compact output
                                            cs_glob_mesh,
                                            cs_shared_connect,
                                            cs_shared_quant,
                                            cs_shared_time_step,
                                            def->input,
                                            dir_val);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid type of definition.\n"
                    " Stop computing the Dirichlet value.\n"));

      } // switch def_type

    } // Definition based on Dirichlet BC
  } // Loop on definitions


#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_COMMON_DBG > 1
  cs_dump_array_to_listing("DIRICHLET_VALUES", quant->n_b_faces, dir_val, 8);
#endif

  return dir_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Tag each face related to a Neumann BC with its definition id.
 *          Default tag is -1 (not a Neumann face)
 *
 * \param[in]      eqp        pointer to a cs_equation_param_t

 * \return an array with prescribed tags
 */
/*----------------------------------------------------------------------------*/

short int *
cs_equation_tag_neumann_face(const cs_equation_param_t  *eqp)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  short int  *face_tag = NULL;

  /* Initialization */
  BFT_MALLOC(face_tag, quant->n_b_faces, short int);

# pragma omp parallel for if (quant->n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_b_faces; f_id++)
    face_tag[f_id] = -1;

  /* Tag faces with Neumann BCs */
  for (int def_id = 0; def_id < eqp->n_bc_desc; def_id++) {

    const cs_xdef_t  *def = eqp->bc_desc[def_id];
    if (def->meta & CS_CDO_BC_NEUMANN) {

      const cs_boundary_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);

#     pragma omp parallel for if (bz->n_faces > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < bz->n_faces; i++)
        face_tag[bz->face_ids[i]] = def_id;

    }
  }

  return face_tag;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to vertices.
 *          Values in the cbc parameters are updated.
 *
 * \param[in]      def_id     id of the definition for setting the Neumann BC
 * \param[in]      f          local face number in the cs_cell_mesh_t
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cbc        pointer to a cs_cell_bc_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sv(short int                   def_id,
                               short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               cs_cell_bc_t               *cbc)
{
  assert(cbc != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  const cs_xdef_t  *def = eqp->bc_desc[def_id];

  assert(def->dim == 3);                 // flux is a vector in the scalar case
  assert(def->meta & CS_CDO_BC_NEUMANN); // Neuman BC

  /* Evaluate the boundary condition at each boundary vertex */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    cs_xdef_eval_cw_at_vtx_flux_by_val(cm, f, def->input, cbc->neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_eval_cw_at_vtx_flux_by_analytic(cm,
                                            f,
                                            cs_shared_time_step,
                                            def->input,
                                            def->qtype,
                                            cbc->neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input =
        (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_desc == 1); // Only one definition allowed
      assert(array_input->stride == 3);
      assert(cs_test_flag(array_input->loc, cs_cdo_primal_face));

      const cs_cdo_quantities_t  *quant = cs_shared_quant;

      cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
      assert(bf_id > -1);

      cs_real_t  *face_val = array_input->values + 3*bf_id;

      cs_xdef_eval_cw_at_vtx_flux_by_val(cm, f, face_val, cbc->neu_values);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } // switch def_type

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs when DoFs are scalar-valued
 *          and attached to faces.
 *          Values in the cbc parameters are set.
 *
 * \param[in]      def_id     id of the definition for setting the Neumann BC
 * \param[in]      f          local face number in the cs_cell_mesh_t
 * \param[in]      eqp        pointer to a cs_equation_param_t
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cbc        pointer to a cs_cell_bc_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_neumann_sf(short int                   def_id,
                               short int                   f,
                               const cs_equation_param_t  *eqp,
                               const cs_cell_mesh_t       *cm,
                               cs_cell_bc_t               *cbc)
{
  assert(cbc != NULL && cm != NULL && eqp != NULL);
  assert(def_id > -1);

  const cs_xdef_t  *def = eqp->bc_desc[def_id];

  assert(def->dim == 3);                 // flux is a vector in the scalar case
  assert(def->meta & CS_CDO_BC_NEUMANN); // Neuman BC

  /* Evaluate the boundary condition at each boundary vertex */
  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    cs_xdef_eval_cw_flux_by_val(cm, f, def->input, cbc->neu_values);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_eval_cw_flux_by_analytic(cm,
                                     f,
                                     cs_shared_time_step,
                                     def->input,
                                     def->qtype,
                                     cbc->neu_values);
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input =
        (cs_xdef_array_input_t *)def->input;

      assert(eqp->n_bc_desc == 1); // Only one definition allowed
      assert(array_input->stride == 3);
      assert(cs_test_flag(array_input->loc, cs_cdo_primal_face));

      const cs_cdo_quantities_t  *quant = cs_shared_quant;

      cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
      assert(bf_id > -1);

      cs_real_t  *face_val = array_input->values + 3*bf_id;

      cs_xdef_eval_cw_flux_by_val(cm, f, face_val, cbc->neu_values);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition.\n"
                " Stop computing the Neumann value.\n"));

  } // switch def_type

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system related to cell vertices into the global
 *         algebraic system
 *
 * \param[in]       csys      cellwise view of the algebraic system
 * \param[in]       rset      pointer to a cs_range_set_t structure on vertices
 * \param[in]       sys_flag  flag associated to the current system builder
 * \param[in, out]  rhs       array storing the right-hand side
 * \param[in, out]  sources   array storing the contribution of source terms
 * \param[in, out]  mav       pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_v(const cs_cell_sys_t            *csys,
                       const cs_range_set_t           *rset,
                       cs_flag_t                       sys_flag,
                       cs_real_t                      *rhs,
                       cs_real_t                      *sources,
                       cs_matrix_assembler_values_t   *mav)
{
  const short int  n_vc = csys->mat->n_ent;
  const cs_lnum_t  *v_ids = csys->mat->ids;

  cs_gnum_t  grows[CS_CDO_ASSEMBLE_BUF_SIZE], gcols[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_real_t  vals[CS_CDO_ASSEMBLE_BUF_SIZE];

  /* Assemble the matrix related to the advection/diffusion/reaction terms
     If advection is activated, the resulting system is not symmetric
     Otherwise, the system is symmetric with extra-diagonal terms. */
  /* TODO: Add a symmetric version for optimization */
  int  block_size = 0;
  for (short int i = 0; i < n_vc; i++) {

    const double  *mval_i = csys->mat->val + i*n_vc;
    const cs_lnum_t  vi_id = v_ids[i];

#   pragma omp atomic
    rhs[vi_id] += csys->rhs[i];

    if (sys_flag & CS_FLAG_SYS_SOURCETERM) {
#     pragma omp atomic
      sources[vi_id] += csys->source[i];
    }

    cs_gnum_t  grow_id = rset->g_id[vi_id];

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    for (short int j = 0; j < n_vc; j++) {

      grows[block_size] = grow_id;
      gcols[block_size] = rset->g_id[v_ids[j]];
      vals[block_size] = mval_i[j];
      block_size += 1;

      if (block_size == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
        cs_matrix_assembler_values_add_g(mav,
                                         CS_CDO_ASSEMBLE_BUF_SIZE,
                                         grows, gcols, vals);
        block_size = 0;
      }

    } /* Loop on cell vertices (local cols) */

  } /* Loop on cell vertices (local rows) */

  if (block_size > 0) {
#   pragma omp critical
    cs_matrix_assembler_values_add_g(mav, block_size, grows, gcols, vals);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system related to cell faces into the global
 *         algebraic system
 *
 * \param[in]       csys      cellwise view of the algebraic system
 * \param[in]       rset      pointer to a cs_range_set_t structure on vertices
 * \param[in]       sys_flag  flag associated to the current system builder
 * \param[in, out]  rhs       array storing the right-hand side
 * \param[in, out]  mav       pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_f(const cs_cell_sys_t            *csys,
                       const cs_range_set_t           *rset,
                       cs_flag_t                       sys_flag,
                       cs_real_t                      *rhs,
                       cs_matrix_assembler_values_t   *mav)
{
  CS_UNUSED(sys_flag);

  const short int  n_fc = csys->mat->n_ent;
  const cs_lnum_t  *f_ids = csys->mat->ids;

  cs_gnum_t  grows[CS_CDO_ASSEMBLE_BUF_SIZE], gcols[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_real_t  vals[CS_CDO_ASSEMBLE_BUF_SIZE];

  /* Assemble the matrix related to the advection/diffusion/reaction terms
     If advection is activated, the resulting system is not symmetric
     Otherwise, the system is symmetric with extra-diagonal terms. */
  /* TODO: Add a symmetric version for optimization */
  int  block_size = 0;
  for (short int i = 0; i < n_fc; i++) {

    const double  *val_rowi = csys->mat->val + i*n_fc;
    const cs_lnum_t  fi_id = f_ids[i];

#   pragma omp atomic
    rhs[fi_id] += csys->rhs[i];

  const cs_gnum_t  grow_id = rset->g_id[fi_id];

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    for (short int colj = 0; colj < n_fc; colj++) {

      grows[block_size] = grow_id;
      gcols[block_size] = rset->g_id[f_ids[colj]];
      vals[block_size] = val_rowi[colj];
      block_size += 1;

      if (block_size == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
        cs_matrix_assembler_values_add_g(mav,
                                         CS_CDO_ASSEMBLE_BUF_SIZE,
                                         grows, gcols, vals);
        block_size = 0;
      }

    } /* Loop on cell faces (local cols) */

  } /* Loop on cell faces (local rows) */

  if (block_size > 0) {
#   pragma omp critical
    cs_matrix_assembler_values_add_g(mav, block_size, grows, gcols, vals);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_structure_t according
 *         to the space scheme
 *
 * \param[in]  scheme       enum on the discretization scheme used
 *
 * \return  a pointer on a cs_matrix_structure_t *
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_structure_t *
cs_equation_get_matrix_structure(cs_space_scheme_t   scheme)
{
  cs_matrix_structure_t  *ms = NULL;

  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    ms = cs_equation_common_ms[CS_EQ_COMMON_VERTEX];
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    ms = cs_equation_common_ms[CS_EQ_COMMON_FACE];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Invalid space scheme.");
    break;

  } // Switch on space scheme

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_assembler_t according
 *         to the space scheme
 *
 * \param[in]  scheme       enum on the discretization scheme used
 *
 * \return  a pointer on a cs_matrix_assembler_t *
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_assembler_t *
cs_equation_get_matrix_assembler(cs_space_scheme_t   scheme)
{
  cs_matrix_assembler_t  *mav = NULL;

  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    mav = cs_equation_common_ma[CS_EQ_COMMON_VERTEX];
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    mav = cs_equation_common_ma[CS_EQ_COMMON_FACE];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Invalid space scheme.");
    break;

  } // Switch on space scheme

  return mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity vertex->vertices for the local rank
 *
 * \return  a pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t *
cs_equation_get_v2v_index(void)
{
  return cs_connect_v2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity face->faces for the local rank
 *
 * \return  a pointer to a cs_adjacency_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_adjacency_t *
cs_equation_get_f2f_index(void)
{
  return cs_connect_f2f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_equation_get_tmpbuf(void)
{
  return cs_equation_common_work_buffer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_equation_get_tmpbuf_size(void)
{
  return cs_equation_common_work_buffer_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a monitoring structure
 *
 * \return a cs_equation_monitor_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_monitor_t *
cs_equation_init_monitoring(void)
{
  cs_equation_monitor_t  *m = NULL;

  BFT_MALLOC(m, 1, cs_equation_monitor_t);

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(m->tcb); // build system

  CS_TIMER_COUNTER_INIT(m->tcd); // build diffusion terms
  CS_TIMER_COUNTER_INIT(m->tca); // build advection terms
  CS_TIMER_COUNTER_INIT(m->tcr); // build reaction terms
  CS_TIMER_COUNTER_INIT(m->tcs); // build source terms

  CS_TIMER_COUNTER_INIT(m->tce); // extra operations

  return m;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a message in the performance output file related to the
 *          monitoring of equation
 *
 * \param[in]  eqname    pointer to the name of the current equation
 * \param[in]  monitor   monitoring structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_write_monitoring(const char                    *eqname,
                             const cs_equation_monitor_t   *monitor)
{
  double t[6] = {monitor->tcb.wall_nsec, monitor->tcd.wall_nsec,
                 monitor->tca.wall_nsec, monitor->tcr.wall_nsec,
                 monitor->tcs.wall_nsec, monitor->tce.wall_nsec};
  for (int i = 0; i < 6; i++) t[i] *= 1e-9;

  if (eqname == NULL)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  " %-35s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f seconds\n",
                  "<CDO/Equation> Monitoring",
                  t[0], t[1], t[2], t[3], t[4], t[5]);
  else {
    char *msg = NULL;
    int len = 1 + strlen("<CDO/> Monitoring") + strlen(eqname);
    BFT_MALLOC(msg, len, char);
    sprintf(msg, "<CDO/%s> Monitoring", eqname);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  " %-35s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f seconds\n",
                  msg, t[0], t[1], t[2], t[3], t[4], t[5]);
    BFT_FREE(msg);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
