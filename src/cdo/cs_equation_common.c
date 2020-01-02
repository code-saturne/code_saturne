/*============================================================================
 * Routines to handle common features for building algebraic system in CDO
 * schemes
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_boundary_zone.h"
#include "cs_cdo_local.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_log.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_xdef_eval.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_common.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_EQUATION_COMMON_DBG               0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

/* Temporary buffers useful during the building of all algebraic systems */
static size_t  cs_equation_common_work_buffer_size = 0;
static cs_real_t  *cs_equation_common_work_buffer = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  quant        pointer to additional mesh quantities struct.
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  eb_flag      metadata for Edge-based schemes
 * \param[in]  fb_flag      metadata for Face-based schemes
 * \param[in]  vb_flag      metadata for Vertex-based schemes
 * \param[in]  vcb_flag     metadata for Vertex+Cell-basde schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_init(const cs_cdo_connect_t       *connect,
                        const cs_cdo_quantities_t    *quant,
                        const cs_time_step_t         *time_step,
                        cs_flag_t                     eb_flag,
                        cs_flag_t                     fb_flag,
                        cs_flag_t                     vb_flag,
                        cs_flag_t                     vcb_flag,
                        cs_flag_t                     hho_flag)
{
  assert(connect != NULL && quant != NULL); /* Sanity check */

  /* Allocate cell-wise and face-wise view of a mesh */
  cs_cdo_local_initialize(connect);

  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_lnum_t  n_edges = connect->n_edges;

  /* Allocate shared buffer and initialize shared structures */
  size_t  cwb_size = n_cells; /* initial cell-wise buffer size */

  /* Allocate and initialize matrix assembler and matrix structures */
  if (vb_flag > 0 || vcb_flag > 0) {

    if (vb_flag & CS_FLAG_SCHEME_SCALAR || vcb_flag & CS_FLAG_SCHEME_SCALAR) {

      if (vb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_vertices);

      if (vcb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)(n_vertices + n_cells));

    } /* scalar-valued equations */

    if (vb_flag & CS_FLAG_SCHEME_VECTOR || vcb_flag & CS_FLAG_SCHEME_VECTOR) {

      cwb_size = CS_MAX(cwb_size, (size_t)3*n_cells);
      if (vb_flag & CS_FLAG_SCHEME_VECTOR)
        cwb_size = CS_MAX(cwb_size, (size_t)3*n_vertices);

      if (vcb_flag & CS_FLAG_SCHEME_VECTOR)
        cwb_size = CS_MAX(cwb_size, (size_t)3*(n_vertices + n_cells));

    } /* vector-valued equations */

  } /* Vertex-based schemes and related ones */

  if (eb_flag > 0) {

    if (eb_flag & CS_FLAG_SCHEME_SCALAR) {

      /* This is a vector-valued equation but the DoF is scalar-valued since
       * it is a circulation associated to each edge */

      cwb_size = CS_MAX(cwb_size, (size_t)3*n_cells);
      cwb_size = CS_MAX(cwb_size, (size_t)n_edges);

    } /* vector-valued equations with scalar-valued DoFs */

  } /* Edge-based schemes */

  if (fb_flag > 0 || hho_flag > 0) {

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR)) {

      assert(n_faces > n_cells);
      if (fb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_faces);

      if (hho_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_faces);

    } /* Scalar-valued CDO-Fb or HHO-P0 */

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY1 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR)) {

      assert((CS_CDO_CONNECT_FACE_SP1 == CS_CDO_CONNECT_FACE_VP0) &&
             (CS_CDO_CONNECT_FACE_SP1 == CS_CDO_CONNECT_FACE_VHP0));

      cwb_size = CS_MAX(cwb_size, (size_t)CS_N_FACE_DOFS_1ST * n_faces);

    } /* Vector CDO-Fb or HHO-P1 or vector HHO-P0 */

    if (cs_flag_test(hho_flag,
                     CS_FLAG_SCHEME_POLY2 | CS_FLAG_SCHEME_SCALAR))
      cwb_size = CS_MAX(cwb_size, (size_t)CS_N_FACE_DOFS_2ND * n_faces);

    /* For vector equations and HHO */
    if (cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY1) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY2)) {

      if  (hho_flag & CS_FLAG_SCHEME_POLY1)
        cwb_size = CS_MAX(cwb_size, (size_t)3*CS_N_FACE_DOFS_1ST*n_faces);

      else if  (hho_flag & CS_FLAG_SCHEME_POLY2)
        cwb_size = CS_MAX(cwb_size, (size_t)3*CS_N_FACE_DOFS_2ND*n_faces);

    }

  } /* Face-based schemes (CDO or HHO) */

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
 * \brief  Free buffers shared among the equations solved with CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_common_finalize(void)
{
  /* Free cell-wise and face-wise view of a mesh */
  cs_cdo_local_finalize();

  /* Free common buffer */
  BFT_FREE(cs_equation_common_work_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a new structure to handle the building of algebraic system
 *         related to a cs_equation_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_init_builder(const cs_equation_param_t   *eqp,
                         const cs_mesh_t             *mesh)
{
  cs_equation_builder_t  *eqb = NULL;

  BFT_MALLOC(eqb, 1, cs_equation_builder_t);

  eqb->init_step = true;

  /* Initialize flags used to knows what kind of cell quantities to build */
  eqb->msh_flag = 0;
  eqb->bd_msh_flag = 0;
  eqb->st_msh_flag = 0;
  if (eqp->dim > 1)
    eqb->sys_flag = CS_FLAG_SYS_VECTOR;
  else
    eqb->sys_flag = 0;

  /* Handle properties */
  eqb->diff_pty_uniform = true;
  if (cs_equation_param_has_diffusion(eqp))
    eqb->diff_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);

  eqb->curlcurl_pty_uniform = true;
  if (cs_equation_param_has_curlcurl(eqp))
    eqb->curlcurl_pty_uniform = cs_property_is_uniform(eqp->curlcurl_property);

  eqb->graddiv_pty_uniform = true;
  if (cs_equation_param_has_graddiv(eqp))
    eqb->graddiv_pty_uniform = cs_property_is_uniform(eqp->graddiv_property);

  eqb->time_pty_uniform = true;
  if (cs_equation_param_has_time(eqp))
    eqb->time_pty_uniform = cs_property_is_uniform(eqp->time_property);

  if (eqp->n_reaction_terms > CS_CDO_N_MAX_REACTIONS)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Number of reaction terms for an equation is too high.\n"
              " Current value: %d (max: %d)\n"
              " Change the value of CS_CDO_N_MAX_REACTIONS in the code or\n"
              " modify your settings or contact the developpement team.",
              __func__, eqp->n_reaction_terms, CS_CDO_N_MAX_REACTIONS);

  for (int i = 0; i < eqp->n_reaction_terms; i++)
    eqb->reac_pty_uniform[i]
      = cs_property_is_uniform(eqp->reaction_properties[i]);

  /* Handle source terms */
  eqb->source_mask = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    /* Default initialization */
    eqb->st_msh_flag = cs_source_term_init(eqp->space_scheme,
                                           eqp->n_source_terms,
                       (cs_xdef_t *const *)eqp->source_terms,
                                           eqb->compute_source,
                                           &(eqb->sys_flag),
                                           &(eqb->source_mask));

  } /* There is at least one source term */

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */
  eqb->face_bc = cs_cdo_bc_face_define(eqp->default_bc,
                                       true, /* Steady BC up to now */
                                       eqp->dim,
                                       eqp->n_bc_defs,
                                       eqp->bc_defs,
                                       mesh->n_b_faces);

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(eqb->tcb); /* build system */
  CS_TIMER_COUNTER_INIT(eqb->tcs); /* solve system */
  CS_TIMER_COUNTER_INIT(eqb->tce); /* extra operations */

  return eqb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_builder_t structure
 *
 * \param[in, out]  p_builder  pointer of pointer to the cs_equation_builder_t
 *                             structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_builder(cs_equation_builder_t  **p_builder)
{
  if (p_builder == NULL)
    return;
  if (*p_builder == NULL)
    return;

  cs_equation_builder_t  *eqb = *p_builder;

  if (eqb->source_mask != NULL)
    BFT_FREE(eqb->source_mask);

  /* Free BC structure */
  eqb->face_bc = cs_cdo_bc_free(eqb->face_bc);

  BFT_FREE(eqb);

  *p_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the residual normalization at the cellwise level according
 *         to the type of renormalization
 *         Case of scalar-valued system.
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_c           cell volume
 * \param[in]      csys            pointer to a cs_cell_sys_t structure
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_cw_scal_res_normalization(cs_param_resnorm_type_t     type,
                                      cs_real_t                   vol_c,
                                      const cs_cell_sys_t        *csys,
                                      const cs_real_t             weight[],
                                      cs_real_t                  *normalization)
{
  if (type == CS_PARAM_RESNORM_WEIGHTED_RHS) {

    assert(weight != NULL);
    cs_real_t  _rhs_norm = 0;
    for (short int i = 0; i < csys->n_dofs; i++)
      _rhs_norm += weight[i] * csys->rhs[i]*csys->rhs[i];

    *normalization += vol_c * _rhs_norm;

  }
  else if (type == CS_PARAM_RESNORM_FILTERED_RHS) {

    cs_real_t  _rhs_norm = 0;
    if ((csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)   ||
        (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_VERTEX) ||
        csys->has_internal_enforcement) {

      for (short int i = 0; i < csys->n_dofs; i++) {
        if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
          continue;
        else if (csys->intern_forced_ids[i] > -1)
          continue;
        else
          _rhs_norm += csys->rhs[i]*csys->rhs[i];
      }

    }
    else { /* No need to apply a filter */

      for (short int i = 0; i < csys->n_dofs; i++)
        _rhs_norm += csys->rhs[i]*csys->rhs[i];

    }

    *normalization += _rhs_norm;

  } /* Type of residual normalization */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the residual normalization at the cellwise level according
 *         to the type of renormalization
 *         Case of vector-valued system
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_c           cell volume
 * \param[in]      csys            pointer to a cs_cell_sys_t structure
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_cw_vect_res_normalization(cs_param_resnorm_type_t     type,
                                      cs_real_t                   vol_c,
                                      const cs_cell_sys_t        *csys,
                                      const cs_real_t             weight[],
                                      cs_real_t                  *normalization)
{
  if (type == CS_PARAM_RESNORM_WEIGHTED_RHS) {

    assert(weight != NULL);
    cs_real_t  _rhs_norm = 0;
    short int  n_ent = csys->n_dofs/3;
    for (short int i = 0; i < n_ent; i++) {
      const cs_real_t  w = weight[i];
      const cs_real_t  *_rhs = csys->rhs + 3*i;
      _rhs_norm += w * (_rhs[0]*_rhs[0] + _rhs[1]*_rhs[1] + _rhs[2]*_rhs[2]);
    }

    *normalization += vol_c * _rhs_norm;

  }
  else if (type == CS_PARAM_RESNORM_FILTERED_RHS) {

    cs_real_t  _rhs_norm = 0;
    if ((csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)   ||
        (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_VERTEX) ||
        csys->has_internal_enforcement) {

      for (short int i = 0; i < csys->n_dofs; i++) {
        if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
          continue;
        else if (csys->intern_forced_ids[i] > -1)
          continue;
        else
          _rhs_norm += csys->rhs[i]*csys->rhs[i];
      }

    }
    else { /* No need to apply a filter */

      for (short int i = 0; i < csys->n_dofs; i++)
        _rhs_norm += csys->rhs[i]*csys->rhs[i];

    }

    *normalization += _rhs_norm;

  } /* Type of residual normalization */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage to compute of the renormalization coefficient for the
 *         the residual norm of the linear system
 *
 * \param[in]      type            type of renormalization
 * \param[in]      rhs_size        size of the rhs array
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the  residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_res_normalization(cs_param_resnorm_type_t    type,
                                   cs_lnum_t                  rhs_size,
                                   const cs_real_t            rhs[],
                                   cs_real_t                 *normalization)
{
  switch (type) {

  case CS_PARAM_RESNORM_NORM2_RHS:
  case CS_PARAM_RESNORM_FILTERED_RHS:
    *normalization = cs_dot_xx(rhs_size, rhs);
    cs_parall_sum(1, CS_REAL_TYPE, normalization);
    if (*normalization < 100*DBL_MIN)
      *normalization = 1.0;
    else
      *normalization = sqrt((*normalization));
    break;

  case CS_PARAM_RESNORM_WEIGHTED_RHS:
    cs_parall_sum(1, CS_REAL_TYPE, normalization);
    if (*normalization < 100*DBL_MIN)
      *normalization = 1.0;
    else
      *normalization = sqrt((*normalization)/cs_shared_quant->vol_tot);
    break;

  default:
    *normalization = 1.0;
    break;

  } /* Type of normalization */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare a linear system and synchronize buffers to handle parallelism.
 *        Transfer a mesh-based description of arrays x0 and rhs into an
 *        algebraic description for the linear system in x and b.
 *
 * \param[in]      stride     stride to apply to the range set operations
 * \param[in]      x_size     size of the vector unknowns (scatter view)
 * \param[in]      matrix     pointer to a cs_matrix_t structure
 * \param[in]      rset       pointer to a range set structure
 * \param[in]      rhs_redux  do or not a parallel sum reduction on the RHS
 * \param[in, out] x          array of unknowns (in: initial guess)
 * \param[in, out] b          right-hand side
 *
 * \returns the number of non-zeros in the matrix
 */
/*----------------------------------------------------------------------------*/

cs_gnum_t
cs_equation_prepare_system(int                     stride,
                           cs_lnum_t               x_size,
                           const cs_matrix_t      *matrix,
                           const cs_range_set_t   *rset,
                           bool                    rhs_redux,
                           cs_real_t              *x,
                           cs_real_t              *b)
{
  const cs_lnum_t  n_scatter_elts = x_size; /* size of x and rhs */
  const cs_lnum_t  n_gather_elts = cs_matrix_get_n_rows(matrix);

  /* Sanity checks */
  assert(n_gather_elts <= n_scatter_elts);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_COMMON_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                " n_gather_elts:    %d\n"
                " n_scatter_elts:   %d\n"
                " n_matrix_rows:    %d\n"
                " n_matrix_columns: %d\n",
                n_gather_elts, n_scatter_elts, cs_matrix_get_n_rows(matrix),
                cs_matrix_get_n_columns(matrix));
#endif

  if (cs_glob_n_ranks > 1) { /* Parallel mode */
                             /* ============= */

    /* x and b should be changed to have a "gathered" view through the range set
       operation.  Their size is equal to n_sles_gather_elts <=
       n_sles_scatter_elts */

    /* Compact numbering to fit the algebraic decomposition */
    cs_range_set_gather(rset,
                        CS_REAL_TYPE, /* type */
                        stride,       /* stride */
                        x,            /* in: size = n_sles_scatter_elts */
                        x);           /* out: size = n_sles_gather_elts */

    /* The right-hand side stems from a cellwise building on this rank.
       Other contributions from distant ranks may contribute to an element
       owned by the local rank */
    if (rhs_redux)
      cs_interface_set_sum(rset->ifs,
                           n_scatter_elts, stride, false, CS_REAL_TYPE,
                           b);

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        stride,      /* stride */
                        b,           /* in: size = n_sles_scatter_elts */
                        b);          /* out: size = n_sles_gather_elts */
  }

  /* Output information related to the linear system */
  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(matrix, &row_index, &col_id, &d_val, &x_val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_COMMON_DBG > 2
  cs_dbg_dump_linear_system("Dump linear system",
                            n_gather_elts, CS_EQUATION_COMMON_DBG,
                            x, b,
                            row_index, col_id, x_val, d_val);
#endif

  cs_gnum_t  nnz = row_index[n_gather_elts];
  cs_parall_counter(&nnz, 1);

  return nnz;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with scalar-valued
 *         degrees of freedom
 *
 * \param[in]  n_scatter_dofs local number of DoFs (may be != n_gather_elts)
 * \param[in]  eqp            pointer to a cs_equation_param_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  rs             pointer to a cs_range_set_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in]  rhs_redux      do or not a parallel sum reduction on the RHS
 * \param[in, out] sles       pointer to a cs_sles_t structure
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_solve_scalar_system(cs_lnum_t                     n_scatter_dofs,
                                const cs_equation_param_t    *eqp,
                                const cs_matrix_t            *matrix,
                                const cs_range_set_t         *rset,
                                cs_real_t                     normalization,
                                bool                          rhs_redux,
                                cs_sles_t                    *sles,
                                cs_real_t                    *x,
                                cs_real_t                    *b)
{
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_param_sles_t  slesp = eqp->sles_param;

  /* Set xsol */
  cs_real_t  *xsol = NULL;
  if (n_cols > n_scatter_dofs) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_dofs*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Retrieve the solving info structure stored in the cs_field_t structure */
  cs_field_t  *fld = cs_field_by_id(slesp.field_id);
  cs_solving_info_t  sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1, /* stride for scalar-valued */
                                              n_scatter_dofs,
                                              matrix,
                                              rset,
                                              rhs_redux,
                                              xsol,
                                              b);

  /* Solve the linear solver */
  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    slesp.eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (slesp.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d | n_iters %d"
                  " residual % -8.4e | normalization % -8.4e | nnz %lu\n",
                  eqp->name, code, sinfo.n_it, sinfo.res_norm, sinfo.rhs_norm,
                  nnz);

  if (cs_glob_n_ranks > 1) { /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, x);
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_EQUATION_COMMON_DBG > 1
  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        slesp.verbosity,
                        x, b, n_scatter_dofs);
#endif

  if (n_cols > n_scatter_dofs)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a message in the performance output file related to the
 *          monitoring of equation
 *
 * \param[in]  eqname    pointer to the name of the current equation
 * \param[in]  eqb       pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_write_monitoring(const char                    *eqname,
                             const cs_equation_builder_t   *eqb)
{
  double t[3] = {eqb->tcb.wall_nsec, eqb->tcs.wall_nsec, eqb->tce.wall_nsec};
  for (int i = 0; i < 3; i++) t[i] *= 1e-9;

  if (eqname == NULL)
    cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %10.4f %10.4f %10.4f (seconds)\n",
                  "<CDO/Equation> Monitoring", t[0], t[1], t[2]);
  else {

    char *msg = NULL;
    int len = 1 + strlen("<CDO/> Monitoring") + strlen(eqname);

    BFT_MALLOC(msg, len, char);
    sprintf(msg, "<CDO/%s> Monitoring", eqname);
    cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %10.4f %10.4f %10.4f (seconds)\n",
                  msg, t[0], t[1], t[2]);
    BFT_FREE(msg);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties for an algebraic system
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure (diffusion
 *                           property is stored inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_properties(const cs_equation_param_t     *eqp,
                            const cs_equation_builder_t   *eqb,
                            cs_real_t                      t_eval,
                            cs_cell_builder_t             *cb)
{
  /* Preparatory step for diffusion term */
  if (cs_equation_param_has_diffusion(eqp))
    if (eqb->diff_pty_uniform)
      /* One calls this function as if it's a boundary cell to scan all tests */
      cs_equation_set_diffusion_property(eqp,
                                         0, /* cell_id */
                                         t_eval,
                                         CS_FLAG_BOUNDARY_CELL_BY_FACE,
                                         cb);

  /* Preparatory step for curl-curl term */
  if (cs_equation_param_has_curlcurl(eqp))
    if (eqb->curlcurl_pty_uniform)
      /* One calls this function as if it's a boundary cell to scan all tests */
      cs_equation_set_curlcurl_property(eqp,
                                        0, /* cell_id */
                                        t_eval,
                                        CS_FLAG_BOUNDARY_CELL_BY_FACE,
                                        cb);

  /* Preparatory step for grad-div term */
  if (cs_equation_param_has_graddiv(eqp))
    if (eqb->graddiv_pty_uniform)
      cb->gpty_val = cs_property_get_cell_value(0, t_eval,
                                                eqp->graddiv_property);

  /* Preparatory step for unsteady term */
  if (cs_equation_param_has_time(eqp))
    if (eqb->time_pty_uniform)
      cb->tpty_val = cs_property_get_cell_value(0, t_eval, eqp->time_property);

  /* Preparatory step for reaction term */
  if (cs_equation_param_has_reaction(eqp)) {

    for (int i = 0; i < CS_CDO_N_MAX_REACTIONS; i++) cb->rpty_vals[i] = 1.0;

    for (int r = 0; r < eqp->n_reaction_terms; r++) {
      if (eqb->reac_pty_uniform[r]) {
        cb->rpty_vals[r] =
          cs_property_get_cell_value(0, t_eval, eqp->reaction_properties[r]);
      }
    } /* Loop on reaction properties */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties for a given cell when building the
 *         algebraic system. If the property is uniform, a first call has to
 *         be done before the loop on cells
 *         Call \ref cs_eqution_init_properties for instance
 *
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in]      eqb        pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in]      cell_flag  flag related to the current cell
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 *                            (diffusion property is stored inside)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_properties_cw(const cs_equation_param_t     *eqp,
                               const cs_equation_builder_t   *eqb,
                               const cs_real_t                t_eval,
                               const cs_flag_t                cell_flag,
                               const cs_cell_mesh_t          *cm,
                               cs_cell_builder_t             *cb)
{
  /* Set the diffusion property */
  if (cs_equation_param_has_diffusion(eqp))
    if (!(eqb->diff_pty_uniform))
      cs_equation_set_diffusion_property_cw(eqp, cm, t_eval, cell_flag, cb);

  /* Set the property related to the curl-curl operator */
  if (cs_equation_param_has_curlcurl(eqp)) {
    if (!(eqb->curlcurl_pty_uniform))
      cs_equation_set_curlcurl_property_cw(eqp, cm, t_eval, cell_flag, cb);
  }

  /* Set the (linear) reaction property */
  if (cs_equation_param_has_reaction(eqp)) {

    /* Define the local reaction property */
    cb->rpty_val = 0;
    for (int r = 0; r < eqp->n_reaction_terms; r++)
      if (eqb->reac_pty_uniform[r])
        cb->rpty_val += cb->rpty_vals[r];
      else
        cb->rpty_val += cs_property_value_in_cell(cm,
                                                  eqp->reaction_properties[r],
                                                  t_eval);

  }

  if (cs_equation_param_has_time(eqp))
    if (!(eqb->time_pty_uniform))
      cb->tpty_val = cs_property_value_in_cell(cm, eqp->time_property, t_eval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the list of degrees of freedom (DoFs) related to an internal
 *         enforcement. If set to NULL, the array dof_ids (storing the
 *         indirection) is allocated to n_x.
 *
 * \param[in]      n_x        number of entities where DoFs are defined
 * \param[in]      c2x        cell --> x connectivity
 * \param[in]      eqp        set of parameters related to an equation
 * \param[in, out] p_dof_ids  double pointer on DoF ids subject to enforcement
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_dof_enforcement(cs_lnum_t                     n_x,
                                  const cs_adjacency_t         *c2x,
                                  const cs_equation_param_t    *eqp,
                                  cs_lnum_t                    *p_dof_ids[])
{
  if (eqp->n_enforced_dofs == 0 && eqp->n_enforced_cells == 0)
    return;

  /* Initialize the indirection list (by default, no vertex selected) */
  cs_lnum_t  *dof_ids = *p_dof_ids;

  if (dof_ids == NULL)
    BFT_MALLOC(dof_ids, n_x, cs_lnum_t);

# pragma omp parallel for if (n_x > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_x; i++)
    dof_ids[i] = -1;     /* Not selected */

  if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_CELLS) {

    if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE) {

      for (cs_lnum_t i = 0; i < eqp->n_enforced_cells; i++) {
        const cs_lnum_t  c_id = eqp->enforced_cell_ids[i];
        for (cs_lnum_t j = c2x->idx[c_id]; j < c2x->idx[c_id+1]; j++)
          dof_ids[c2x->ids[j]] = i;
      }

    }
    else  /* This case can be tricky in parallel and can also impact the const
             specifier of cs_equation_param_t */
      bft_error(__FILE__, __LINE__, 0,
                "%s: Eq: %s\n"
                "Enforcement by a cell selection without a reference value"
                " is not implemented yet.", __func__, eqp->name);

  }
  else {

    assert(eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_DOFS);

    /* Build the indirection between vertices and enforced vertices */
    for (cs_lnum_t i = 0; i < eqp->n_enforced_dofs; i++)
      dof_ids[eqp->enforced_dof_ids[i]] = i;

  }

  /* Returns pointer */
  *p_dof_ids = dof_ids;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Apply an
 *          algebraic manipulation
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforced_internal_dofs(const cs_equation_param_t       *eqp,
                                   cs_cell_builder_t               *cb,
                                   cs_cell_sys_t                   *csys)
{
  /* Enforcement of the Dirichlet BCs */
  if (csys->has_internal_enforcement == false)
    return;  /* Nothing to do */

  double  *x_vals = cb->values;
  double  *ax = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  /* Build x_vals */
  if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE) {

    const cs_real_t  ref_val = eqp->enforcement_ref_value[0];
    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = ref_val;
    }

  }
  else if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_DOFS) {

    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = eqp->enforced_dof_values[csys->intern_forced_ids[i]];
    }

  }
  else {

    assert(eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_CELLS);
    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = eqp->enforced_cell_values[csys->intern_forced_ids[i]];
    }

  }

  /* Contribution of the DoFs which are enforced */
  cs_sdm_matvec(csys->mat, x_vals, ax);

  /* Second pass: Replace the block of enforced DoFs by a diagonal block */
  for (int i = 0; i < csys->n_dofs; i++) {

    if (csys->intern_forced_ids[i] > -1) {

      /* Reset row i */
      memset(csys->mat->val + csys->n_dofs*i, 0, csys->n_dofs*sizeof(double));
      /* Reset column i */
      for (int j = 0; j < csys->n_dofs; j++)
        csys->mat->val[i + csys->n_dofs*j] = 0;
      csys->mat->val[i*(1 + csys->n_dofs)] = 1;

      /* Set the RHS */
      csys->rhs[i] = x_vals[i];

    } /* DoF associated to a Dirichlet BC */
    else
      csys->rhs[i] -= ax[i];  /* Update RHS */

  } /* Loop on degrees of freedom */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Case of matrices
 *          defined by blocks. Apply an algebraic manipulation.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforced_internal_block_dofs(const cs_equation_param_t     *eqp,
                                         cs_cell_builder_t             *cb,
                                         cs_cell_sys_t                 *csys)
{
  /* Enforcement of the Dirichlet BCs */
  if (csys->has_internal_enforcement == false)
    return;  /* Nothing to do */

  double  *x_vals = cb->values;
  double  *ax = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE) {

    /* Build x_vals */
    const cs_real_t  *ref_val = eqp->enforcement_ref_value;
    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = ref_val[i%3];
    }

  }
  else if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_DOFS) {

    /* Build x_vals */
    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = eqp->enforced_dof_values[csys->intern_forced_ids[i]];
    }

  }
  else {

    /* Build x_vals */
    assert(eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_CELLS);
    for (int i = 0; i < csys->n_dofs; i++) {
      if (csys->intern_forced_ids[i] > -1)
        x_vals[i] = eqp->enforced_cell_values[csys->intern_forced_ids[i]];
    }

  }

  /* Contribution of the DoFs which are enforced */
  cs_sdm_block_matvec(csys->mat, x_vals, ax);

  /* Define the new right-hand side (rhs) */
  for (int i = 0; i < csys->n_dofs; i++) {
    if (csys->intern_forced_ids[i] > -1)
      csys->rhs[i] = x_vals[i];
    else
      csys->rhs[i] -= ax[i];  /* Update RHS */
  }

  const cs_sdm_block_t  *bd = csys->mat->block_desc;

  /* Second pass: Replace the block of enforced DoFs by a diagonal block */
  int s = 0;
  for (int ii = 0; ii < bd->n_row_blocks; ii++) {

    cs_sdm_t  *db = cs_sdm_get_block(csys->mat, ii, ii);
    const int  bsize = db->n_rows*db->n_cols;

    if (csys->intern_forced_ids[s] > -1) {

      /* Identity for the diagonal block */
      memset(db->val, 0, sizeof(cs_real_t)*bsize);
      for (int i = 0; i < db->n_rows; i++) {
        db->val[i*(1+db->n_rows)] = 1;
        assert(csys->intern_forced_ids[s+i] > -1);
      }

      /* Reset column and row block jj < ii */
      for (int jj = 0; jj < ii; jj++) {

        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, ii, jj);
        memset(bij->val, 0, sizeof(cs_real_t)*bsize);

        cs_sdm_t  *bji = cs_sdm_get_block(csys->mat, jj, ii);
        memset(bji->val, 0, sizeof(cs_real_t)*bsize);

      }

      /* Reset column and row block jj < ii */
      for (int jj = ii+1; jj < db->n_rows; jj++) {

        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, ii, jj);
        memset(bij->val, 0, sizeof(cs_real_t)*bsize);

        cs_sdm_t  *bji = cs_sdm_get_block(csys->mat, jj, ii);
        memset(bji->val, 0, sizeof(cs_real_t)*bsize);

      }

    } /* DoF associated to an enforcement of their values*/

    s += db->n_rows;

  } /* Loop on degrees of freedom */

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
 * \brief  Allocate a cs_equation_balance_t structure
 *
 * \param[in]  location   where the balance is performed
 * \param[in]  size       size of arrays in the structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_equation_balance_create(cs_flag_t    location,
                           cs_lnum_t    size)
{
  cs_equation_balance_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_equation_balance_t);

  b->size = size;
  b->location = location;
  if (cs_flag_test(location, cs_flag_primal_cell) == false &&
      cs_flag_test(location, cs_flag_primal_vtx) == false)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid location", __func__);

  BFT_MALLOC(b->balance, 7*size, cs_real_t);
  b->unsteady_term  = b->balance +   size;
  b->reaction_term  = b->balance + 2*size;
  b->diffusion_term = b->balance + 3*size;
  b->advection_term = b->balance + 4*size;
  b->source_term    = b->balance + 5*size;
  b->boundary_term  = b->balance + 6*size;

  /* Set to zero all members */
  cs_equation_balance_reset(b);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_equation_balance_t structure
 *
 * \param[in, out] b     pointer to a cs_equation_balance_t to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_reset(cs_equation_balance_t   *b)
{
  if (b == NULL)
    return;
  if (b->size < 1)
    return;

  if (b->balance == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: array is not allocated.", __func__);

  size_t  bufsize = b->size *sizeof(cs_real_t);

  memset(b->balance, 0, 7*bufsize);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize balance terms if this is a parallel computation
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in, out] b         pointer to a cs_equation_balance_t to rsync
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_sync(const cs_cdo_connect_t    *connect,
                         cs_equation_balance_t     *b)
{
  if (cs_glob_n_ranks < 2)
    return;
  if (b == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: structure not allocated", __func__);

  if (cs_flag_test(b->location, cs_flag_primal_vtx)) {

    assert(b->size == connect->n_vertices);
    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         b->size,
                         7,   /* stride: 1 for each kind of balance */
                         false,
                         CS_REAL_TYPE,
                         b->balance);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_balance_t structure
 *
 * \param[in, out]  p_balance  pointer to the pointer to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_balance_destroy(cs_equation_balance_t   **p_balance)
{
  cs_equation_balance_t *b = *p_balance;

  if (b == NULL)
    return;

  BFT_FREE(b->balance);

  BFT_FREE(b);
  *p_balance = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each vertex
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2v_idx   index array  to define
 * \param[in, out]  def2v_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_vol_def_at_vertices(const cs_cdo_connect_t  *connect,
                                     int                      n_defs,
                                     cs_xdef_t              **defs,
                                     cs_lnum_t                def2v_idx[],
                                     cs_lnum_t                def2v_ids[])
{
  if (n_defs == 0)
    return;

  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_adjacency_t  *c2v = connect->c2v;

  int  *v2def_ids = NULL;
  BFT_MALLOC(v2def_ids, n_vertices, int);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t v = 0; v < n_vertices; v++)
    v2def_ids[v] = -1;          /* default */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v = 0; v < n_vertices; v++)
        v2def_ids[v] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          v2def_ids[c2v->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL] != NULL);
    /* Last definition is used if there is a conflict between several
       definitions */
    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         n_vertices,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         v2def_ids);

  }

  /* 0. Initialization */
  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2v_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count number of vertices related to each definition */
  for (cs_lnum_t v = 0; v < n_vertices; v++)
    if (v2def_ids[v] > -1)
      def2v_idx[v2def_ids[v]+1] += 1;

  /* 2. Build index */
  for (int def_id = 0; def_id < n_defs; def_id++)
    def2v_idx[def_id+1] += def2v_idx[def_id];

  /* 3. Build list */
  for (cs_lnum_t v = 0; v < n_vertices; v++) {
    const int def_id = v2def_ids[v];
    if (def_id > -1) {
      def2v_ids[def2v_idx[def_id] + count[def_id]] = v;
      count[def_id] += 1;
    }
  }

  /* Free memory */
  BFT_FREE(v2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each edge
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2e_idx   index array  to define
 * \param[in, out]  def2e_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_vol_def_at_edges(const cs_cdo_connect_t  *connect,
                                  int                      n_defs,
                                  cs_xdef_t              **defs,
                                  cs_lnum_t                def2e_idx[],
                                  cs_lnum_t                def2e_ids[])
{
  if (n_defs == 0)
    return;

  const cs_lnum_t  n_edges = connect->n_edges;
  const cs_adjacency_t  *c2e = connect->c2e;

  int  *e2def_ids = NULL;
  BFT_MALLOC(e2def_ids, n_edges, int);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t e = 0; e < n_edges; e++)
    e2def_ids[e] = -1; /* default: not associated to a definition */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_edges; e++)
        e2def_ids[e] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++)
          e2def_ids[c2e->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_EDGE_SCAL] != NULL);
    /* Last definition is used if there is a conflict between several
       definitions */
    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_EDGE_SCAL],
                         n_edges,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         e2def_ids);

  }

  /* 0. Initialization */
  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2e_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count the number of edges related to each definition */
  for (cs_lnum_t e = 0; e < n_edges; e++)
    if (e2def_ids[e] > -1)
      def2e_idx[e2def_ids[e]+1] += 1;

  /* 2. Build the index */
  for (int def_id = 0; def_id < n_defs; def_id++)
    def2e_idx[def_id+1] += def2e_idx[def_id];

  /* 3. Build the list */
  for (cs_lnum_t e = 0; e < n_edges; e++) {
    const int def_id = e2def_ids[e];
    if (def_id > -1) {
      def2e_ids[def2e_idx[def_id] + count[def_id]] = e;
      count[def_id] += 1;
    }
  }

  /* Free memory */
  BFT_FREE(e2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each face
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2f_idx   index array  to define
 * \param[in, out]  def2f_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_vol_def_at_faces(const cs_cdo_connect_t    *connect,
                                  int                        n_defs,
                                  cs_xdef_t                **defs,
                                  cs_lnum_t                  def2f_idx[],
                                  cs_lnum_t                  def2f_ids[])
{
  if (n_defs == 0)
    return;

  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_adjacency_t  *c2f = connect->c2f;

  int  *f2def_ids = NULL;
  BFT_MALLOC(f2def_ids, n_faces, int);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_faces; f++)
    f2def_ids[f] = -1;          /* default */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */
    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_faces; f++)
        f2def_ids[f] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
          f2def_ids[c2f->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_FACE_SP0] != NULL);
    /* Last definition is used if there is a conflict between several
       definitions */
    cs_interface_set_max(connect->interfaces[CS_CDO_CONNECT_FACE_SP0],
                         n_faces,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         f2def_ids);

  }

  /* 0. Initialization */
  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2f_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count number of faces related to each definition */
  for (cs_lnum_t f = 0; f < n_faces; f++)
    if (f2def_ids[f] > -1)
      def2f_idx[f2def_ids[f]+1] += 1;

  /* 2. Build index */
  for (int def_id = 0; def_id < n_defs; def_id++)
    def2f_idx[def_id+1] += def2f_idx[def_id];

  /* 3. Build list */
  for (cs_lnum_t f = 0; f < n_faces; f++) {
    const int def_id = f2def_ids[f];
    if (def_id > -1) {
      def2f_ids[def2f_idx[def_id] + count[def_id]] = f;
      count[def_id] += 1;
    }
  }

  /* Free memory */
  BFT_FREE(f2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mean-value across ranks at each vertex
 *
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       dim         number of entries for each vertex
 * \param[in]       counter     number of occurences on this rank
 * \param[in, out]  values      array to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_sync_vertex_mean_values(const cs_cdo_connect_t     *connect,
                                    int                         dim,
                                    int                        *counter,
                                    cs_real_t                  *values)
{
  const cs_lnum_t  n_vertices = connect->n_vertices;

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    assert(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL] != NULL);

    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         n_vertices,
                         1,           /* stride */
                         false,       /* interlace (not useful here) */
                         CS_INT_TYPE, /* int */
                         counter);

    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         n_vertices,
                         dim,         /* stride */
                         true,        /* interlace */
                         CS_REAL_TYPE,
                         values);

  }

  if (dim == 1) {

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      if (counter[v_id] > 1)
        values[v_id] /= counter[v_id];

  }
  else { /* dim > 1 */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      if (counter[v_id] > 1) {
        const cs_real_t  inv_count = 1./counter[v_id];
        for (int k = 0; k < dim; k++)
          values[dim*v_id + k] *= inv_count;
      }
    }

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
