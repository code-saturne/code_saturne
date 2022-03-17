/*============================================================================
 * Set of helper functions to prepare or solve linear systems
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

#include <assert.h>
#include <stdlib.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_math.h"
#include "cs_parall.h"
#include "cs_parameters.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_solve.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_CDO_SOLVE_DBG        0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 1
static int cs_cdo_solve_dbg_counter = 0;  /* Id number for debug */
#endif

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the rhs norm used as a renormalization
 *         coefficient for the residual norm when solving the linear system. A
 *         pre-computation performed during the cellwise building of the
 *         algebraic has been done before this call according to the requested
 *         type of renormalization.
 *
 * \param[in]      type            type of renormalization
 * \param[in]      vol_tot         total volume of the computational domain
 * \param[in]      rhs_size        size of the rhs array
 * \param[in]      rhs             array related to the right-hand side
 * \param[in, out] normalization   value of the residual normalization
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_solve_sync_rhs_norm(cs_param_resnorm_type_t    type,
                           double                     vol_tot,
                           cs_lnum_t                  rhs_size,
                           const cs_real_t            rhs[],
                           double                    *normalization)
{
  switch (type) {

  case CS_PARAM_RESNORM_NORM2_RHS:
    *normalization = cs_dot_xx(rhs_size, rhs);
    cs_parall_sum(1, CS_REAL_TYPE, normalization);
    if (*normalization < 100*DBL_MIN)
      *normalization = 1.0;
    else
      *normalization = sqrt((*normalization));
    break;

  case CS_PARAM_RESNORM_FILTERED_RHS:
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
      *normalization = sqrt((*normalization)/vol_tot);
    break;

  default:
    *normalization = 1.0;
    break;

  } /* Type of normalization */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare a linear system and synchronize buffers in case of parallel
 *        or periodic computations. Transfer a mesh-based description of
 *        arrays x0 and rhs into an algebraic description for the linear system
 *        in x and b (scatter/gather views).
 *
 * \param[in]      stride     stride to apply to the range set operations
 * \param[in]      x_size     size of the vector unknowns (scatter view)
 * \param[in]      matrix     pointer to a cs_matrix_t structure
 * \param[in]      rset       pointer to a range set structure
 * \param[in]      rhs_redux  do or not a parallel sum reduction on the RHS
 * \param[in, out] x          array of unknowns (in: initial guess)
 * \param[in, out] b          right-hand side
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_solve_prepare_system(int                     stride,
                            cs_lnum_t               x_size,
                            const cs_matrix_t      *matrix,
                            const cs_range_set_t   *rset,
                            bool                    rhs_redux,
                            cs_real_t              *x,
                            cs_real_t              *b)
{
  const cs_lnum_t  n_scatter_elts = x_size; /* size of x and rhs */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 0
  const cs_lnum_t  n_gather_elts = cs_matrix_get_n_rows(matrix);
  assert(n_gather_elts <= n_scatter_elts);

  cs_log_printf(CS_LOG_DEFAULT,
                " n_gather_elts:    %d\n"
                " n_scatter_elts:   %d\n"
                " n_matrix_rows:    %d\n"
                " n_matrix_columns: %d\n",
                n_gather_elts, n_scatter_elts, cs_matrix_get_n_rows(matrix),
                cs_matrix_get_n_columns(matrix));
#else
  CS_UNUSED(matrix);
#endif

  if (rset != NULL) { /* Parallel or periodic mode
                         ========================= */

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

    /* TODO the system is presumed to have interlaced = true for vector
     * equations. No difference for scalar equations. Maybe its value
     * could be passed into by the caller of this function */

    if (rhs_redux && rset->ifs != NULL)
      cs_interface_set_sum(rset->ifs,
                           n_scatter_elts, stride, true, CS_REAL_TYPE,
                           b);

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        stride,      /* stride */
                        b,           /* in: size = n_sles_scatter_elts */
                        b);          /* out: size = n_sles_gather_elts */
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 2
  const cs_lnum_t  *row_index, *col_id;
  const cs_real_t  *d_val, *x_val;

  cs_matrix_get_msr_arrays(matrix, &row_index, &col_id, &d_val, &x_val);

  cs_dbg_dump_linear_system("Dump linear system",
                            n_gather_elts, CS_CDO_SOLVE_DBG,
                            x, b,
                            row_index, col_id, x_val, d_val);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising with scalar-valued cell-based DoFs*
 *         No rotation is taken into account when synchronizing halo.
 *
 * \param[in]  n_dofs         local number of DoFs
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
 * \param[in]  matrix         pointer to a cs_matrix_t structure
 * \param[in]  normalization  value used for the residual normalization
 * \param[in, out] sles       pointer to a cs_sles_t structure
 * \param[in, out] x          solution of the linear system (in: initial guess)
 * \param[in, out] b          right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdo_solve_scalar_cell_system(cs_lnum_t                n_dofs,
                                const cs_param_sles_t   *slesp,
                                const cs_matrix_t       *matrix,
                                cs_real_t                normalization,
                                cs_sles_t               *sles,
                                cs_real_t               *x,
                                cs_real_t               *b)
{
  /* Retrieve the solving info structure stored in the cs_field_t structure */

  cs_solving_info_t  sinfo;
  cs_field_t  *fld = NULL;
  if (slesp->field_id > -1) {
    fld = cs_field_by_id(slesp->field_id);
    cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);
  }

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

  const cs_halo_t  *halo = cs_matrix_get_halo(matrix);
  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(matrix);
  const cs_lnum_t  n_cols_ext = cs_matrix_get_n_columns(matrix);

  assert(n_dofs == n_rows);

  /* Handle parallelism */

  cs_real_t  *_x = x, *_b = b;
  if (n_cols_ext > n_rows) {

    BFT_MALLOC(_b, n_cols_ext, cs_real_t);
    BFT_MALLOC(_x, n_cols_ext, cs_real_t);

    memcpy(_x, x, n_dofs*sizeof(cs_real_t));
    memcpy(_b, b, n_dofs*sizeof(cs_real_t));

    cs_matrix_pre_vector_multiply_sync(matrix, _b);
    cs_halo_sync_var(halo, CS_HALO_STANDARD, _x);
  }

  /* Solve the linear solver */

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    slesp->eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    _b,
                                                    _x,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  if (n_cols_ext > n_rows) {
    BFT_FREE(_b);
    memcpy(x, _x, n_dofs*sizeof(cs_real_t));
    BFT_FREE(_x);
  }

  /* Output information about the convergence of the resolution */

  if (slesp->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d>"
                  " n_iter %3d | res.norm % -8.4e | rhs.norm % -8.4e\n",
                  slesp->name, code,
                  sinfo.n_it, sinfo.res_norm, sinfo.rhs_norm);

  if (slesp->field_id > -1)
    cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with scalar-valued
 *         degrees of freedom
 *
 * \param[in]  n_scatter_dofs local number of DoFs (may be != n_gather_elts)
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
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
cs_cdo_solve_scalar_system(cs_lnum_t                     n_scatter_dofs,
                           const cs_param_sles_t        *slesp,
                           const cs_matrix_t            *matrix,
                           const cs_range_set_t         *rset,
                           cs_real_t                     normalization,
                           bool                          rhs_redux,
                           cs_sles_t                    *sles,
                           cs_real_t                    *x,
                           cs_real_t                    *b)
{
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

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

  cs_field_t  *fld = cs_field_by_id(slesp->field_id);
  cs_solving_info_t  sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

  /* Prepare solving (handle parallelism)
   * stride = 1 for scalar-valued */

  cs_cdo_solve_prepare_system(1, n_scatter_dofs, matrix, rset, rhs_redux,
                              xsol, b);

  /* Solve the linear solver */

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    slesp->eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */

  if (slesp->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d>"
                  " n_iter %3d | res.norm % -8.4e | rhs.norm % -8.4e\n",
                  slesp->name, code,
                  sinfo.n_it, sinfo.res_norm, sinfo.rhs_norm);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       xsol, x);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 1, /* type and stride */
                       b, b);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 1
  cs_dbg_fprintf_system(slesp->name, cs_cdo_solve_dbg_counter++,
                        slesp->verbosity,
                        x, b, n_scatter_dofs);
#endif

  if (n_cols > n_scatter_dofs)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from CDO schemes with vector-valued
 *         degrees of freedom
 *
 * \param[in]  n_scatter_dofs local number of DoFs (may be != n_gather_elts)
 * \param[in]  slesp          pointer to a cs_param_sles_t structure
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
cs_cdo_solve_vector_system(cs_lnum_t                     n_scatter_dofs,
                           const cs_param_sles_t        *slesp,
                           const cs_matrix_t            *matrix,
                           const cs_range_set_t         *rset,
                           cs_real_t                     normalization,
                           bool                          rhs_redux,
                           cs_sles_t                    *sles,
                           cs_real_t                    *x,
                           cs_real_t                    *b)
{
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(matrix);

  /* Set xsol */

  cs_real_t  *xsol = NULL;
  if (n_cols > n_rows) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, 3*n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_dofs/3*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Retrieve the solving info structure stored in the cs_field_t structure */

  cs_field_t  *fld = cs_field_by_id(slesp->field_id);
  cs_solving_info_t  sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  sinfo.rhs_norm = normalization;

  /* Prepare solving (handle parallelism)
   * stride = 3 for vector-valued */

  cs_cdo_solve_prepare_system(3, n_scatter_dofs/3, matrix, rset, rhs_redux,
                              xsol, b);

  /* Solve the linear solver */

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    slesp->eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */

  if (slesp->verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%20s/sles_cvg_code=%-d>"
                  " n_iter %3d | res.norm % -8.4e | rhs.norm % -8.4e\n",
                  slesp->name, code,
                  sinfo.n_it, sinfo.res_norm, sinfo.rhs_norm);

  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 3, /* type and stride */
                       xsol, x);
  cs_range_set_scatter(rset,
                       CS_REAL_TYPE, 3, /* type and stride */
                       b, b);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_SOLVE_DBG > 1
  cs_dbg_fprintf_system(slesp->name, cs_cdo_solve_dbg_counter++,
                        slesp->verbosity,
                        x, b, n_scatter_dofs);
#endif

  if (n_cols > n_rows)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
